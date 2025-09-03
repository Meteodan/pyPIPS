# remove_fields_nc.py
#
# This script removes select fields from a nc file using xarray
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
parser.add_argument('--vars', nargs='*', dest='vars', default=None,
                    help='list of variables to remove')
parser.add_argument('--var-pattern', dest='var_pattern', default=None,
                    help='If provided, remove all variables containing this pattern in their name')
parser.add_argument('--output-tag', dest='output_tag', default='',
                    help='tag for output nc files to distinguish from original if desired')

args = parser.parse_args()

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
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', None)
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Get a list of the combined parsivel and conventional netCDF data files that are present in the
# PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames_nc]

for conv_file, parsivel_combined_file in zip(conv_filelist, parsivel_combined_filelist):
    print(f"Reading {parsivel_combined_file}")  # noqa: T201
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    print(f"Reading {conv_file}")  # noqa: T201
    conv_ds = xr.load_dataset(conv_file)

    varlist = [k for k, v in parsivel_combined_ds.items()]
    # vars_to_remove = [var for var in varlist if any(x in var for x in ['KGWX', 'retr', 'KHTX'])]
    # vars_to_remove = [var for var in varlist if any(x in var for x in ['_retr_'])]
    vars_to_remove = [var for var in varlist if var in args.vars] if args.vars else []
    if args.var_pattern:
        more_vars_to_remove = [var for var in varlist if args.var_pattern in var]
        vars_to_remove.extend(more_vars_to_remove)
        vars_to_remove = list(set(vars_to_remove))
    print("Variables to remove from parsivel combined file: ", vars_to_remove)  # noqa: T201
    parsivel_combined_ds = parsivel_combined_ds.drop_vars(vars_to_remove)

    varlist = [k for k, v in conv_ds.items()]
    # vars_to_remove = [var for var in varlist if any(x in var for x in ['KGWX', 'retr', 'KHTX'])]
    # vars_to_remove = [var for var in varlist if any(x in var for x in ['_retr_'])]
    vars_to_remove = [var for var in varlist if var in args.vars] if args.vars else []
    if args.var_pattern:
        more_vars_to_remove = [var for var in varlist if args.var_pattern in var]
        vars_to_remove.extend(more_vars_to_remove)
        vars_to_remove = list(set(vars_to_remove))
    print("Variables to remove from conventional file: ", vars_to_remove)  # noqa: T201
    conv_ds = conv_ds.drop_vars(vars_to_remove)

    conv_output_file = conv_file + args.output_tag
    parsivel_combined_output_file = parsivel_combined_file + args.output_tag
    print(f"Dumping {parsivel_combined_output_file}")  # noqa: T201
    parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
    print(f"Dumping {conv_output_file}")  # noqa: T201
    conv_ds.to_netcdf(conv_output_file)
