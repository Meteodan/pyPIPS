# pyPIPS_meteograms.py
#
# This script calculates radar retrievals from PIPS and collocated radar data (netCDF version)
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
description = "Calculates radar retreivals from PIPS and collocated radar data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--ND-tag', dest='ND_tag', default='qc',
                    help='tag for ND variable in file (either qc or RB15)')
parser.add_argument('--calc-for-SATP', action='store_true', dest='calc_for_SATP',
                    help='calculate for the SATP-filtered dataset')
parser.add_argument('--coefficients', nargs=3, metavar=('c1', 'c2', 'c3'), type=float,
                    dest='coefficients',
                    help='coefficients of mu-lambda polynomial (in decreasing order of exponent')


args = parser.parse_args()
ND_tag = args.ND_tag

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
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
comp_radar = config.radar_config_dict.get('comp_radar', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', 'NEXRAD')
radar_dir = config.radar_config_dict.get('radar_dir', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

scatt_file = os.path.join(scatt_dir, 'SCTT_RAIN_fw100.dat')

if not args.calc_for_SATP:
    # Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
    parsivel_combined_filenames = [
        'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name,
                                                  int(requested_interval))
        for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
    parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf)
                                  for pcf in parsivel_combined_filenames]
else:
    # Just read in the single combined SATP dataset
    parsivel_combined_filename = 'ND_avg_{}_{:d}s.nc'.format(dataset_name, int(requested_interval))
    parsivel_combined_filelist = [os.path.join(PIPS_dir, parsivel_combined_filename)]

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)

    if not args.calc_for_SATP:
        DSD_interval = parsivel_combined_ds.DSD_interval
        PIPS_name = parsivel_combined_ds.probe_name
        deployment_name = parsivel_combined_ds.deployment_name
        ND = parsivel_combined_ds['ND_{}'.format(ND_tag)]
        coord_to_combine = 'time'
    else:
        ND = parsivel_combined_ds['SATP_ND_{}'.format(dataset_name)]
        ND = ND.rename({'D0_RR_level_0': 'D0_idx', 'D0_RR_level_1': 'RR_idx'})
        DSD_interval = ND.DSD_interval
        ND = pipsio.reconstruct_MultiIndex(ND, ['D0_idx', 'RR_idx'], 'D0_RR')
        coord_to_combine = 'D0_RR'

    ND = ND.where(ND > 0.)
    D = parsivel_combined_ds.coords['diameter']
    dD = (parsivel_combined_ds.coords['max_diameter'] - parsivel_combined_ds.coords['min_diameter'])
    dualpol_dict = dp.calpolrain(wavelength, scatt_file, ND, dD)

    ZH = dualpol_dict['REF']
    ZDR = dualpol_dict['ZDR']
    fa2 = dualpol_dict['fa2']
    fb2 = dualpol_dict['fb2']

    mu_lamda_coeff = args.coefficients[::-1]
    retr_dict = dsd.retrieval_Cao_xr(ZH, ZDR, ND, D, dD, fa2, fb2, wavelength, mu_lamda_coeff)

    retr_ds = xr.Dataset(retr_dict)
    parsivel_combined_ds.update(retr_ds)
    if args.calc_for_SATP:
        parsivel_combined_ds = parsivel_combined_ds.reset_index('D0_RR')
        parsivel_combined_ds.attrs = ND.attrs
    parsivel_combined_ds.attrs['CG_coeff'] = mu_lamda_coeff
    parsivel_combined_ds.attrs['retrieval_wavelength'] = wavelength

    print("Dumping {}".format(parsivel_combined_file))
    parsivel_combined_ds.to_netcdf(parsivel_combined_file)
