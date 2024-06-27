# filter_radar_sweep.py
#
# This script filters radar data from CFRadial sweeps and adds the new filtered fields to the file
from __future__ import annotations

import argparse
import json
import os
from glob import glob
from pathlib import Path

import numpy as np
import pyart
from scipy import ndimage

import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.timemodule as tm
from pyPIPS import utils


def roundPartial(value, resolution, decimals=4):
    return np.around(np.round(value / resolution) * resolution, decimals=decimals)


min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Filters radar data for radar sweeps"
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--el-req', type=float, dest='el_req_cl', default=None,
                    help='Requested elevation angle (overrides value in config file).')
parser.add_argument('--dBZ-thresh', type=float, dest='dBZ_thresh', default=5.,
                    help='Threshold of reflectivity below which to exclude (dBZ)')
parser.add_argument('--RHV-thresh', type=float, dest='RHV_thresh', default=0.95,
                    help='Threshold of RHV below which to exclude')
parser.add_argument('--med-filter-footprint', type=json.loads, dest='med_filter_footprint',
                    default=[[1, 1, 1, 1, 1]],
                    help='Footprint of median filter (uses scipy.ndimage.median_filter)')
help_msg = ('tag to determine filename variant for input nc files (V06 when produced by pyART, '
            'SUR when produced by RadxConvert)')
parser.add_argument('--fname-variant', dest='fname_variant', default='V06', help=help_msg)
parser.add_argument('--input-tag', dest='input_tag', default=None,
                    help='Input nametag to determine which files to read in')
parser.add_argument('--output-tag', dest='output_tag', default='filt',
                    help='tag for output nc files to distinguish from original')

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
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', None)
radar_dir = config.radar_config_dict.get('radar_dir', None)
radar_fname_pattern = config.radar_config_dict.get('radar_fname_pattern', None)
# Add the input filename tag to the pattern if needed
if args.input_tag:
    radar_fname_pattern = radar_fname_pattern.replace('.', f'_{args.input_tag}.')
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = args.el_req_cl if args.el_req_cl else config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Read radar sweeps
# First get a list of all potentially relevant radar files in the directory
if args.input_tag is None:
    radar_paths = glob(radar_dir + f'/*{radar_name}*{args.fname_variant}.nc')
else:
    radar_paths = glob(radar_dir + f'/*{radar_name}*_{args.input_tag}.nc')
# Then find only those between the requested times
radar_dict = radar.get_radar_paths_between_times(radar_paths, radar_start_timestamp,
                                                 radar_end_timestamp, radar_type=radar_type,
                                                 fname_format=radar_fname_pattern)
# TODO: refactor this. Shouldn't need a separate function for XTRRA ideally
# And, anyway, this is broken since it yields a different format for the radar_dict dictionary
# that doesn't work anymore
if radar_type == 'XTRRA':
    radar_dict = radar.get_radar_paths_single_elevation(radar_dict, el_req=el_req,
                                                        radar_type=radar_type)

if el_req < 0.0:    # Only should be used for testing.
    radar_dict_out = radar.read_vols(radar_dict)
else:
    # TODO: really should try to avoid reading all the sweeps into memory at once. Change this
    # back to where it reads the radar object one at a time...
    radar_dict_out = radar.read_sweeps(radar_dict, el_req=el_req)
    # Now, loop through the radar sweeps and construct new file names for the filtered output
    # files using the file name format string in the config file, and tacking on
    # the output tag name as necessary.
    # First, treat the filename pattern as a path and extract the stem and extension from it
    if args.output_tag:
        radar_fname_pattern_path = Path(radar_fname_pattern)
        radar_fname_pattern_stem = radar_fname_pattern_path.stem
        radar_fname_pattern_ext = radar_fname_pattern_path.suffixes[-1]
        # Now modify the stem to tack on the output tag
        radar_fname_pattern_stem_new = radar_fname_pattern_stem + "_{output_tag}"
        # Now put the extension back on
        radar_fname_pattern_out = radar_fname_pattern_stem_new + radar_fname_pattern_ext
    else:
        radar_fname_pattern_out = radar_fname_pattern

    # Now, loop through the sweeps and construct new file names using the new format
    radar_sweep_list = radar_dict_out['rad_sweep_list']
    radar_output_fname_list = []
    radar_sweep_time_list = []
    for radar_sweep in radar_sweep_list:
        radar_sweep_time = pyart.graph.common.generate_radar_time_sweep(radar_sweep, 0)
        radar_sweep_time_list.append(radar_sweep_time)
        if args.output_tag:
            radar_sweep_fname = radar_fname_pattern_out.format(rad_name=radar_name,
                                                               year=radar_sweep_time.year,
                                                               month=radar_sweep_time.month,
                                                               day=radar_sweep_time.day,
                                                               hour=radar_sweep_time.hour,
                                                               min=radar_sweep_time.minute,
                                                               sec=radar_sweep_time.second,
                                                               output_tag=args.output_tag)
        else:
            radar_sweep_fname = radar_fname_pattern_out.format(rad_name=radar_name,
                                                               year=radar_sweep_time.year,
                                                               month=radar_sweep_time.month,
                                                               day=radar_sweep_time.day,
                                                               hour=radar_sweep_time.hour,
                                                               min=radar_sweep_time.minute,
                                                               sec=radar_sweep_time.second)
        radar_output_fname_list.append(radar_sweep_fname)
    radar_output_paths = [os.path.join(radar_dir, rof) for rof in radar_output_fname_list]

for radar_obj, radar_sweep_time, radar_output_path in zip(radar_sweep_list,
                                                          radar_sweep_time_list,
                                                          radar_output_paths):
    print(f"Working on time {radar_sweep_time.strftime(tm.timefmt2)}")  # noqa: T201
    print(f"Output filename will be {os.path.basename(radar_output_path)}")  # noqa: T201
    print("Getting fields")  # noqa: T201
    # Get polarimetric fields from the radar object
    ZH_rad_tuple = radar.get_field_to_plot(radar_obj, radar.REF_aliases)
    ZDR_rad_tuple = radar.get_field_to_plot(radar_obj, radar.ZDR_aliases)
    RHV_rad_tuple = radar.get_field_to_plot(radar_obj, radar.RHV_aliases)
    ZH_name = ZH_rad_tuple[0]
    ZDR_name = ZDR_rad_tuple[0]
    RHV_name = RHV_rad_tuple[0]

    # Make copies of polarimetric fields
    for field_name in [ZH_name, ZDR_name, RHV_name]:
        radar_obj.add_field_like(field_name, field_name + '_filtered',
                                 radar_obj.fields[field_name]['data'].copy(), replace_existing=True)

    print("Applying median filter with the following shape: ", args.med_filter_footprint)  # noqa: T201
    for field_name in [ZH_name, ZDR_name, RHV_name]:
        radar_obj.fields[field_name + '_filtered']['data'] = \
            ndimage.median_filter(radar_obj.fields[field_name + '_filtered']['data'].copy(),
                                  footprint=args.med_filter_footprint)
        # radar_obj.fields[field_name + '_filtered']['data'] = \
        #     medfilt2d(radar_obj.fields[field_name + '_filtered']['data'],
        #               kernel_size=args.med_filter_width)

    print("Creating dBZ and RHV gate filter")  # noqa: T201
    # Create a gate filter to mask out areas with dBZ and RHV below thresholds
    rhoHV_ref_filter = pyart.correct.moment_based_gate_filter(radar_obj,
                                                              rhv_field=RHV_name + '_filtered',
                                                              refl_field=ZH_name + '_filtered',
                                                              min_ncp=None,
                                                              min_rhv=args.RHV_thresh,
                                                              min_refl=args.dBZ_thresh,
                                                              max_refl=None)

    print("Applying gate filter")  # noqa: T201
    # Mask fields using the dBZ/RHV mask
    for field_name in [ZH_name, ZDR_name, RHV_name]:
        radar_obj.fields[field_name + '_filtered']['data'] = \
            np.ma.masked_where(rhoHV_ref_filter.gate_excluded,
                               radar_obj.fields[field_name + '_filtered']['data'])

    # Add metadata to the radar object regarding the filtering
    radar_obj.metadata['dBZ_threshold'] = args.dBZ_thresh
    radar_obj.metadata['RHV_threshold'] = args.RHV_thresh
    radar_obj.metadata['median_filter_footprint'] = repr(args.med_filter_footprint)

    # Now save the radar object to a new CFRadial file with the added filtered fields
    print(f"Saving {radar_output_path}")  # noqa: T201
    pyart.io.write_cfradial(radar_output_path, radar_obj)