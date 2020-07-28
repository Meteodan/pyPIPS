# pyPIPS_meteograms.py
#
# This script fits a CG mu-lambda relation to the PIPS data and plots it along with the mu-lambda
# pairs
import os
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
description = "Computes and plots CG relations from PIPS data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--ND-tag', dest='ND_tag', default=None,
                    help='Tag for ND variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--plot-SATP', action='store_true', dest='plot_SATP',
                    help='plot for the SATP-filtered dataset')
# parser.add_argument('--CG-name', dest='CG_name', default='',
#                     help='Name of CG fit we are plotting')
parser.add_argument('--filter_RR', dest='filter_RR', type=float, default=None,
                    help='filter rainrate < # (mm)')
parser.add_argument('--filter_counts', dest='filter_counts', type=int, default=None,
                    help='filter particle counts < #')
parser.add_argument('--use-parsivel-params', dest='use_parsivel_params', action='store_true',
                    default=False, help='Use parsivel RR and counts instead of computed')
parser.add_argument('--fig-fmt', dest='figfmt', default='png', help='format of saved figure')
parser.add_argument('--update-CG-coeff-attrs', dest='update_CG_coeff_attrs', action='store_true',
                    default=False,
                    help='Write CG coefficients as attributes back to PIPS nc files?')

args = parser.parse_args()
if not args.ND_tag:
    ND_tag = ''
else:
    ND_tag = '_{}'.format(args.ND_tag)

if args.plot_SATP:
    filter_RR = None
    filter_counts = None
else:
    filter_RR = args.filter_RR
    filter_counts = args.filter_counts

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

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

if not args.plot_SATP:

    # print(DSD_fits_filenames)

    # Open entire multi-file dataset
    print("Opening dataset {}".format(dataset_name))
    parsivel_combined_ds = xr.open_mfdataset(parsivel_combined_filelist, combine='nested',
                                             concat_dim='time', preprocess=pipsio.remove_unneeded)
    plot_tag = 'no_SATP'
    dim = 'time'
else:
    parsivel_combined_filename = 'ND_avg_{}{}_{:d}s.nc'.format(dataset_name, ND_tag,
                                                               int(requested_interval))
    parsivel_combined_filepath = os.path.join(PIPS_dir, parsivel_combined_filename)
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_filepath)
    plot_tag = 'SATP'
    dim = 'D0_RR'

DSD_interval = parsivel_combined_ds.DSD_interval

# Filter by RR and pcount if desired
filter_RR_tag = 'no_RR_filter'
filter_counts_tag = 'no_counts_filter'
if filter_RR:
    if args.use_parsivel_params:
        rainrate_key = 'precipintensity'
    else:
        rainrate_key = 'rainrate_derived{}'.format(ND_tag)
    parsivel_combined_ds = parsivel_combined_ds.where(
        parsivel_combined_ds[rainrate_key] >= filter_RR, drop=True)
    filter_RR_tag = 'RR_filter_{:d}mm'.format(int(filter_RR))
if filter_counts:
    if args.use_parsivel_params:
        counts_key = 'pcount'
    else:
        counts_key = 'pcount_derived{}'.format(ND_tag)
    parsivel_combined_ds = parsivel_combined_ds.where(
        parsivel_combined_ds[counts_key] >= filter_counts, drop=True)
    filter_counts_tag = 'counts_filter_{:d}'.format(filter_counts)


# Get the untruncated and truncated moment fits MM246.
# TODO: update this for greater flexibility later
DSD_MM246 = parsivel_combined_ds['DSD_MM246{}'.format(ND_tag)]
DSD_TMM246 = parsivel_combined_ds['DSD_TMM246{}'.format(ND_tag)]

# Drop points where mu > 30 or lambda > 20000
DSD_MM246 = DSD_MM246.where(DSD_MM246.sel(parameter='lamda') < 20000.)
DSD_MM246 = DSD_MM246.where(DSD_MM246.sel(parameter='alpha') < 30.)

DSD_MM246 = DSD_MM246.dropna(dim=dim, how='any')

DSD_TMM246 = DSD_TMM246.where(DSD_TMM246.sel(parameter='lamda') < 20000.)
DSD_TMM246 = DSD_TMM246.where(DSD_TMM246.sel(parameter='alpha') < 30.)
DSD_TMM246 = DSD_TMM246.dropna(dim=dim, how='any')

# Fit polynomials
lamda_MM246 = DSD_MM246.sel(parameter='lamda') / 1000.  # Get to mm^-1
mu_MM246 = DSD_MM246.sel(parameter='alpha')

MM246_poly_coeff, MM246_poly = dsd.calc_CG_polynomial(lamda_MM246, mu_MM246)
print("The MM246 polynomial coefficients are: ", MM246_poly_coeff)

# Plot the relationship on a scatterplot along with the Cao and Zhang relations
fig, ax = pm.plot_mu_lamda(MM246_poly_coeff, MM246_poly, lamda=lamda_MM246, mu=mu_MM246,
                           title='Untruncated')
plt.savefig(plot_dir +
            '/Untruncated_MM246_mu_lamda{}_{}_{}_{}.{}'.format(ND_tag, plot_tag, filter_RR_tag,
                                                               filter_counts_tag, args.figfmt),
            dpi=200, bbox_inches='tight')

# Fit polynomials
lamda_TMM246 = DSD_TMM246.sel(parameter='lamda') / 1000.  # Get to mm^-1
mu_TMM246 = DSD_TMM246.sel(parameter='alpha')

TMM246_poly_coeff, TMM246_poly = dsd.calc_CG_polynomial(lamda_TMM246, mu_TMM246)
print("The TMM246 polynomial coefficients are: ", TMM246_poly_coeff)

# Plot the relationship on a scatterplot along with the Cao and Zhang relations
pm.plot_mu_lamda(TMM246_poly_coeff, TMM246_poly, lamda=lamda_TMM246, mu=mu_TMM246,
                 title='Truncated')
plt.savefig(plot_dir +
            '/Truncated_MM246_mu_lamda{}_{}_{}_{}.{}'.format(ND_tag, plot_tag, filter_RR_tag,
                                                             filter_counts_tag, args.figfmt),
            dpi=200, bbox_inches='tight')

if args.update_CG_coeff_attrs:
    parsivel_combined_ds.close()
    print("Saving CG coefficients to nc files as global attributes.")
    filtering = False
    if filter_counts:
        filtering = (filter_counts > 100)
    if filter_RR:
        filtering = filtering | (filter_RR > 0.1)
    for parsivel_combined_file in parsivel_combined_filelist:
        print("Reading {}".format(parsivel_combined_file))
        parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
        if 'CG_coeff_Z01' not in parsivel_combined_ds.attrs:
            parsivel_combined_ds.attrs['CG_coeff_Z01'] = [-1.957, 1.213, -0.0160]
        if 'CG_coeff_C08' not in parsivel_combined_ds.attrs:
            parsivel_combined_ds.attrs['CG_coeff_C08'] = [-1.718, 0.902, -0.0201]
        # TODO: This is some ugly logic. Think about refactoring...
        if args.plot_SATP:
            parsivel_combined_ds.attrs['CG_coeff_SATP_TMM{}'.format(ND_tag)] = TMM246_poly_coeff
            parsivel_combined_ds.attrs['CG_coeff_SATP_MM{}'.format(ND_tag)] = MM246_poly_coeff
            # The following is for backwards compatibility
            # TODO: Still needed?
            parsivel_combined_ds.attrs['CG_coeff_SATP{}'.format(ND_tag)] = TMM246_poly_coeff
        else:
            if filtering:
                parsivel_combined_ds.attrs['CG_coeff_TMM_F{}'.format(ND_tag)] = TMM246_poly_coeff
                parsivel_combined_ds.attrs['CG_coeff_MM_F{}'.format(ND_tag)] = MM246_poly_coeff
                parsivel_combined_ds.attrs['Filter_for_MM_F_and_TMM_F_fits'] = \
                    [filter_RR_tag, filter_counts_tag]
            else:
                parsivel_combined_ds.attrs['CG_coeff_TMM{}'.format(ND_tag)] = TMM246_poly_coeff
                parsivel_combined_ds.attrs['CG_coeff_MM{}'.format(ND_tag)] = MM246_poly_coeff
        parsivel_combined_ds.close()
        # Save updated dataset back to file
        parsivel_combined_ds.to_netcdf(parsivel_combined_file)

