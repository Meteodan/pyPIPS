# pyPIPS_meteograms.py
#
# This script combines DSDs from individual deployment netCDF files from the PIPS
import os
import sys
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
import logging

# TODO: Figure out how to turn off the netCDF logging messages! Below doesn't work....
# logging.getLogger("xarray").setLevel(logging.ERROR)

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Reads PIPS comma-delimited text data files and converts to netCDF"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
args = parser.parse_args()

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

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory

parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
print(parsivel_combined_filenames)

ND_list = []
rho_list = []
for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    if index == 0:
        DSD_interval = parsivel_combined_ds.DSD_interval
    else:
        if parsivel_combined_ds.DSD_interval != DSD_interval:
            print("Problem! One of the files in the list has a DSD interval that doesn't match!")
            print("Aborting!")
        sys.exit()
    ND = parsivel_combined_ds['ND_qc']
    ND.attrs = parsivel_combined_ds.attrs
    rho = parsivel_combined_ds['rho']
    ND = ND.where(ND > 0.)
    ND = ND.dropna(dim='time', how='all')
    rho = rho.reindex_like(ND)

    ND_list.append(ND)
    rho_list.append(rho)

# Concatenate everything together and copy over metadata
print("Combining ND data!")
ND_combined_da = xr.concat(ND_list, dim='time')
rho_combined_da = xr.concat(rho_list, dim='time')

ND_combined_ds = xr.Dataset({'ND_qc': ND_combined_da, 'rho': rho_combined_da})
ND_combined_ds.attrs = ND_combined_da.attrs

ND_combined_ncfile_name = 'ND_{}_combined_{}_{:d}s'.format(ND_tag, dataset_name, int(DSD_interval))
ND_combined_ncfile_path = os.path.join(PIPS_dir, ND_combined_ncfile_name)
print("Dumping {}".format(ND_combined_ncfile_path))
ND_combined_ds.to_netcdf(ND_combined_ncfile_path)

