# radar_to_PIPS.py
#
# This script interpolates radar observations to the PIPS locations and times
from __future__ import annotations

import argparse
import os
from glob import glob

import xarray as xr

import pyPIPS.parsivel_params as pp
import pyPIPS.pips_io as pipsio
import pyPIPS.radarmodule as radar
from pyPIPS import utils

# NOTE: the following has *no effect* on this running script. Need to set the
# environment variable prior to running the script for some reason.
# os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

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
description = "interpolates radar observations to the PIPS locations and times"
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--average-gates', dest='average_gates', default=False, action='store_true',
                    help='Whether to average the nearest gates when interpolating to PIPS location')
parser.add_argument('--ngates2avg', dest='ngates2avg', default=1,
                    help='Width of the halo of gates to average around center gate')
help_msg = ('tag to determine filename variant for input nc files (V06 when produced by pyART, '
            'SUR when produced by RadxConvert)')
parser.add_argument('--fname-variant', dest='fname_variant', default='V06', help=help_msg)
parser.add_argument('--radar-input-tag', dest='radar_input_tag', default=None,
                    help='Input nametag to determine which radar files to read in')
parser.add_argument('--input-tag', dest='input_tag', default=None,
                    help='Input nametag to determine which PIPS files to read in')
parser.add_argument('--output-tag', dest='output_tag', default=None,
                    help='tag for output PIPS nc files to distinguish from original if desired')

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
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
input_parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
comp_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', 'NEXRAD')
radar_dir = config.radar_config_dict.get('radar_dir', None)
radar_fname_pattern = config.radar_config_dict.get('radar_fname_pattern', None)
# Add the input filename tag to the pattern if needed
if args.radar_input_tag:
    radar_fname_pattern = radar_fname_pattern.replace('.', f'_{args.radar_input_tag}.')
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Modify PIPS file names to account for an input tag (e.g., a variant version)
if args.input_tag is not None:
    parsivel_combined_filenames = []
    for parsivel_combined_filename_temp in input_parsivel_combined_filenames:
        parsivel_combined_filename = \
            parsivel_combined_filename_temp.replace(".nc", f"_{args.input_tag}.nc")
else:
    parsivel_combined_filenames = input_parsivel_combined_filenames

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# The following assumes that the same radar will be used for each PIPS in the deployment.
# TODO: make this more flexible
# Read radar sweeps
# First get a list of all potentially relevant radar files in the directory
if args.radar_input_tag is None:
    radar_paths = glob(radar_dir + f'/*{radar_name}*{args.fname_variant}.nc')
else:
    radar_paths = glob(radar_dir + f'/*{radar_name}*_{args.radar_input_tag}.nc')
# Then find only those between the requested times
radar_path_dict = radar.get_radar_paths_between_times(radar_paths, radar_start_timestamp,
                                                      radar_end_timestamp, radar_type=radar_type,
                                                      fname_format=radar_fname_pattern)
if radar_type == 'XTRRA':
    radar_path_dict = radar.get_radar_paths_single_elevation(radar_path_dict, el_req=el_req,
                                                             radar_type=radar_type)

# Outer file loop
# TODO: change this to interpolate the radar fields to each PIPS in one go, without having to
# reload
for parsivel_combined_path in parsivel_combined_filelist:
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_path)
    PIPS_name = parsivel_combined_ds.probe_name
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(float, geo_loc_str.strip('()').split(',')))
    radar_fields_at_PIPS_da, beam_height_da = \
        radar.interp_sweeps_to_one_PIPS_new(radar_name, radar_path_dict, PIPS_name, geo_loc,
                                            el_req=el_req, average_gates=args.average_gates)
    # Get rid of existing interpolated fields. NOTE: for some reason this doesn't always seem
    # to work. I've had to run this script *twice* in order for the removal of the previous
    # version to "stick". This seems like a bug and doesn't make a lot of sense.
    # NOTE: looks like the addition of the "drop_dims" call fixed this issue. And looks like
    # I don't even need the next 2 calls but will leave them in there just in case.
    if f'{radar_name}_at_PIPS' in parsivel_combined_ds:
        parsivel_combined_ds = parsivel_combined_ds.drop_vars(f'{radar_name}_at_PIPS',
                                                              errors='ignore')
        parsivel_combined_ds = parsivel_combined_ds.drop_vars(f'fields_{radar_name}',
                                                              errors='ignore')
    if f'{radar_name}_beam_height_at_PIPS' in parsivel_combined_ds:
        parsivel_combined_ds = \
            parsivel_combined_ds.drop_vars(f'{radar_name}_beam_height_at_PIPS', errors='ignore')
    if f'fields_{radar_name}' in parsivel_combined_ds.dims:
        parsivel_combined_ds = parsivel_combined_ds.drop_dims([f'fields_{radar_name}'],
                                                              errors='ignore')
    # Interpolate radar fields to the PIPS times
    radar_fields_at_PIPS_da = radar_fields_at_PIPS_da.interp_like(parsivel_combined_ds)
    radar_fields_at_PIPS_da.attrs['elevation_angle'] = el_req
    parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds,
                                                        radar_fields_at_PIPS_da,
                                                        name=f'{radar_name}_at_PIPS')
    beam_height_da = beam_height_da.interp_like(parsivel_combined_ds)
    beam_height_da.attrs['elevation_angle'] = el_req
    parsivel_combined_ds = \
        pipsio.combine_parsivel_data(parsivel_combined_ds, beam_height_da,
                                     name=f'{radar_name}_beam_height_at_PIPS')
    # parsivel_combined_ds.close()
    # Save updated dataset back to file
    # Add a new output name tag (variant) if requested
    if args.output_tag:
        parsivel_combined_filename = os.path.basename(parsivel_combined_path)
        parsivel_combined_filename = \
            parsivel_combined_filename.replace(".nc", f"_{args.output_tag}.nc")
        parsivel_combined_path_new = os.path.join(PIPS_dir, parsivel_combined_filename)
        parsivel_combined_ds.to_netcdf(parsivel_combined_path_new)
    else:
        parsivel_combined_ds.to_netcdf(parsivel_combined_path)
