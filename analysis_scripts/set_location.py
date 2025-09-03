# set_location.py
#
# This script sets the location metadata for the PIPS netCDF files
from __future__ import annotations

import argparse
import os

import xarray as xr

import pyPIPS.parsivel_params as pp
from pyPIPS import utils

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Removes select fields from nc files using xarray"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--PIPS-name', dest='PIPS_name', default=None,
                    help='PIPS name to set location for')
parser.add_argument('--location', nargs=3, dest='location',
                    help='location of the PIPS (latitude, longitude, altitude)')

args = parser.parse_args()

geo_loc = tuple(map(float, args.location))

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
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
PIPS_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)
conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', None)
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
PIPS_filelist = [os.path.join(PIPS_dir, pcf) for pcf in PIPS_filenames_nc]
conv_filelist = [os.path.join(PIPS_dir, cfn) for cfn in conv_filenames_nc]

PIPS_file = next(PIPS_filelist[i] for i, name in enumerate(PIPS_names) if name == args.PIPS_name)
conv_file = next(conv_filelist[i] for i, name in enumerate(PIPS_names) if name == args.PIPS_name)


print(f"Reading {PIPS_file}")  # noqa: T201
PIPS_ds = xr.load_dataset(PIPS_file)
print(f"Setting location attribute to {str(geo_loc)}")  # noqa: T201
PIPS_ds.attrs['location'] = str(geo_loc)
print(f"Dumping {PIPS_file}")  # noqa: T201
PIPS_ds.to_netcdf(PIPS_file)

print(f"Reading {conv_file}")  # noqa: T201
conv_ds = xr.load_dataset(conv_file)
print(f"Setting location attribute to {str(geo_loc)}")  # noqa: T201
conv_ds.attrs['location'] = str(geo_loc)
print(f"Dumping {conv_file}")  # noqa: T201
conv_ds.to_netcdf(conv_file)
