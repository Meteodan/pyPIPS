# pyPIPS_meteograms.py
#
# This script calculates DSD fits from PIPS data using the method of moments (netCDF version)
# Note: the script calc_DSD.py also does this, but this version allows one to select on the
# command line which ones to compute. The fits are saved as new variables in the PIPS data
# netCDF file.
import os, sys
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']


# Parse the command line options
description = "Calculates DSD fits using the method of moments from PIPS data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
# parser.add_argument('--ND-tag', dest='ND_tag', default=None,
#                     help='Tag for ND variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--QC-tags', dest='QC_tags', nargs='*', default=[''],
                    help='QC tag for DSD variables in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--calc-for-SATP', action='store_true', dest='calc_for_SATP',
                    help='calculate for the SATP-filtered dataset')
moment_combo_help_string = ('list of moment combos in the form XY or XYZ\n'
                            'where X,Y, and Z are each one of 2, 3, 4, or 6,\n'
                            'unique and in increasing order. In the case of XYZ, both regular\n'
                            'and truncated moment fits will be computed')
parser.add_argument('--moment-combos', dest='moment_combos', nargs='*', default=['246'],
                    help=moment_combo_help_string)

args = parser.parse_args()
# if not args.ND_tag:
#     ND_tag = ''
# else:
#     ND_tag = '_{}'.format(args.ND_tag)

QC_tags = args.QC_tags

# Should always do the derived calculations for the original non-QC'ed DSD
# This next if block makes sure of this
if (len(QC_tags) > 1) and ('' not in QC_tags):
    QC_tags.append('')

# Dynamically import the case configuration file
utils.log("Case config file is {}".format(args.case_config_path))
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")


# Extract needed lists and variables from PIPS_IO_dict configuration dictionary
dataset_name = config.PIPS_IO_dict.get('dataset_name', None)
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

if not args.calc_for_SATP:
    # Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
    parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in
                                  parsivel_combined_filenames]
else:
    # Just read in the single combined SATP dataset for each QC tag
    parsivel_combined_filelist = []
    for QC_tag in QC_tags:
        parsivel_combined_filename = 'ND_avg_{}{}_{:d}s.nc'.format(dataset_name, QC_tag,
                                                                   int(requested_interval))
        parsivel_combined_filelist.append(os.path.join(PIPS_dir, parsivel_combined_filename))

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)

    for QC_tag in QC_tags:
        if QC_tag != '':
            QC_tag = '_{}'.format(QC_tag)

        if not args.calc_for_SATP:
            DSD_interval = parsivel_combined_ds.DSD_interval
            PIPS_name = parsivel_combined_ds.probe_name
            deployment_name = parsivel_combined_ds.deployment_name
            ND = parsivel_combined_ds['ND{}'.format(QC_tag)]
            coord_to_combine = 'time'
        else:
            if 'D0_idx' not in parsivel_combined_ds.coords:
                parsivel_combined_ds = parsivel_combined_ds.rename({'D0_RR_level_0': 'D0_idx'})
            if 'RR_idx' not in parsivel_combined_ds.coords:
                parsivel_combined_ds = parsivel_combined_ds.rename({'D0_RR_level_1': 'RR_idx'})
            if 'D0_RR_level_0' in parsivel_combined_ds.coords:
                parsivel_combined_ds = parsivel_combined_ds.drop('D0_RR_level_0')
            if 'D0_RR_level_1' in parsivel_combined_ds.coords:
                parsivel_combined_ds = parsivel_combined_ds.drop('D0_RR_level_1')

            # ND = ND.rename({'D0_RR_level_0': 'D0_idx', 'D0_RR_level_1': 'RR_idx'})
            parsivel_combined_ds = pipsio.reconstruct_MultiIndex(parsivel_combined_ds,
                                                                ['D0_idx', 'RR_idx'], 'D0_RR')
            # NOTE 09/03/2021: for some reason I am suddenly getting errors further down relating
            # to the "diameter" coordinate not existing. I think it's related to me calling the
            # reconstruct_MultiIndex function on the parsivel_combined_ds instead of just the ND array
            # It wiped out all the original coordinates that weren't the same as the dimension
            # names and changed them back to variables. So, we need to fix that by reassigning
            # them here. Grrrrr....
            parsivel_combined_ds = parsivel_combined_ds.set_coords(['diameter', 'min_diameter',
                                                                    'max_diameter'])
            ND = parsivel_combined_ds['SATP_ND_{}{}'.format(dataset_name, QC_tag)]

            # NOTE: Below try-except block for backwards compatibility
            try:
                DSD_interval = ND.DSD_interval
            except AttributeError:
                DSD_interval = parsivel_combined_ds.DSD_interval

            coord_to_combine = 'D0_RR'
        ND = ND.where(ND > 0.)

        # Compute various fits using the MM and TMM

        M2, _ = dsd.calc_moment_bin(ND, moment=2)
        M3, _ = dsd.calc_moment_bin(ND, moment=3)
        M4, _ = dsd.calc_moment_bin(ND, moment=4)
        M6, _ = dsd.calc_moment_bin(ND, moment=6)
        moment_dict = {
            'M2': M2,
            'M3': M3,
            'M4': M4,
            'M6': M6
        }

        # Compute D_min and D_max
        # Gotcha, if at some point all the zeros in the ND bins have been turn to nan's, the
        # call below to calculate D_min and D_max won't work properly, so here we fill the nans
        # with zeros for the purposes of calculating it.
        D_min, D_max = dsd.get_max_min_diameters(ND.fillna(0.0), dim=coord_to_combine)

        fit_name_list = []
        fit_list = []
        for moment_combo in args.moment_combos:
            num_moments = len(moment_combo)
            moment_list = [moment_dict['M{}'.format(moment_combo[0])],
                           moment_dict['M{}'.format(moment_combo[1])]]
            if num_moments > 2:
                moment_list.append(moment_dict['M{}'.format(moment_combo[2])])

            # Compute regular MM fits
            DSD_MM_fit = dsd.fit_DSD_MMXYZ(moment_combo, moment_list)
            fit_name_list.append('DSD_MM{}'.format(moment_combo))
            fit_list.append(DSD_MM_fit)
            # Compute TMM fits (only for gamma fits)
            if num_moments > 2:
                DSD_TMM_fit = dsd.fit_DSD_TMMXYZ(moment_combo, moment_list, D_min, D_max)
                fit_name_list.append('DSD_TMM{}'.format(moment_combo))
                fit_list.append(DSD_TMM_fit)

        # Wrap fits into a DataSet and dump to netCDF file
        data_arrays = []
        for da_tuple in fit_list:
            da_concat = xr.concat(da_tuple, pd.Index(['N0', 'lamda', 'alpha'], name='parameter'))
            data_arrays.append(da_concat)

        fits_ds = xr.Dataset({'{}{}'.format(name, QC_tag): da for name, da in zip(fit_name_list,
                                                                                  data_arrays)})

        # Compute sigma and Dm and add to Dataset

        # D = ND['diameter']
        # dD = ND['max_diameter'] - ND['min_diameter']

        # Dm = dsd.calc_Dmpq_binned(4, 3, ND)
        # sigma = dsd.calc_sigma(D, dD, ND)

        # fits_ds = pipsio.combine_parsivel_data(fits_ds, Dm, name='Dm43{}'.format(QC_tag),
        #                                        coord=coord_to_combine)
        # fits_ds = pipsio.combine_parsivel_data(fits_ds, sigma, name='sigma{}'.format(QC_tag),
        #                                        coord=coord_to_combine)

        # Update parsivel_combined_ds with new fits_ds and save updated Dataset
        parsivel_combined_ds.update(fits_ds)

    if args.calc_for_SATP:
        # Save attrs (may not need to do this)
        attrs = parsivel_combined_ds.attrs
        parsivel_combined_ds = parsivel_combined_ds.reset_index('D0_RR')
        parsivel_combined_ds.attrs = attrs
        # parsivel_combined_ds.attrs = ND.attrs
    print("Dumping {}".format(parsivel_combined_file))
    parsivel_combined_ds.to_netcdf(parsivel_combined_file)
