# pyPIPS_meteograms.py
#
# This script calculates radar retrievals for entire radar sweeps. Reads in from previously created
# lookup tables
import os
from pathlib import Path
import argparse
from glob import glob
# from datetime import datetime, timedelta
import numpy as np
import pandas as pd
# import xarray as xr
# import matplotlib.ticker as ticker
# import matplotlib.dates as dates
# import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
# import pyPIPS.plotmodule as pm
# import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
# import pyPIPS.PIPS as pips
# import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
# import pyPIPS.DSDlib as dsd
# import pyPIPS.polarimetric as dp
import pyart


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
description = "Calculates radar retreivals for radar sweeps"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--el-req', type=float, dest='el_req_cl', default=None,
                    help='Requested elevation angle (overrides value in config file)')
parser.add_argument('--lookup-dirs', dest='lookup_dirs', nargs='*',
                    help='directory where lookup tables reside')
parser.add_argument('--use-filtered-fields', dest='use_filtered_fields', default=False,
                    action='store_true',
                    help='Whether to use previously filtered dBZ and ZDR fields for the retrieval')
parser.add_argument('--retrieval-tags', dest='retrieval_tags', nargs='*', default=[''],
                    help=('list of nametag for the names of the mu-lambda relation'
                          '(e.g. SATP, C08, Z01)'))
parser.add_argument('--input-tag', dest='input_tag', default=None,
                    help='Input nametag to determine which files to read in')
parser.add_argument('--output-tag', dest='output_tag', default=None,
                    help='tag for output nc files to distinguish from original')
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
    radar_fname_pattern = radar_fname_pattern.replace('.', '_{}.'.format(args.input_tag))
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
if args.el_req_cl:
    el_req = args.el_req_cl
else:
    el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Get a list of radar paths between the start and stop times and containing the elevation
# angle requested
if args.input_tag is None:
    radar_paths = glob(radar_dir + '/*{}*.nc'.format(radar_name))
else:
    radar_paths = glob(radar_dir + '/*{}*_{}.nc'.format(radar_name, args.input_tag))

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

if args.use_filtered_fields:
    tag = '_filtered'
else:
    tag = None

for radar_obj, radar_sweep_time, radar_output_path in zip(radar_sweep_list,
                                                          radar_sweep_time_list,
                                                          radar_output_paths):
    print("Working on time {}".format(radar_sweep_time.strftime(tm.timefmt2)))
    print("Output filename will be {}".format(os.path.basename(radar_output_path)))
    print("Getting ZH and ZDR fields")
    # Get the ZH and ZDR fields from the radar object
    ZH_rad_tuple = radar.get_field_to_plot(radar_obj, radar.REF_aliases, tag=tag)
    ZDR_rad_tuple = radar.get_field_to_plot(radar_obj, radar.ZDR_aliases, tag=tag)
    ZH_rad = ZH_rad_tuple[1]['data']
    ZDR_rad = ZDR_rad_tuple[1]['data']
    # Get the masks for both ZH_rad and ZDR_rad. These will be used later to mask the retrieved
    # values
    ZH_mask = ZH_rad.mask
    ZDR_mask = ZDR_rad.mask
    full_mask = np.ma.mask_or(ZH_mask, ZDR_mask)
    # Read in first lookup table to get the interval between reflectivity and ZDR
    lookup_path = os.path.join(args.lookup_dirs[0], 'D0.csv')
    retr_table = pd.read_csv(lookup_path, sep=',', header=0, index_col='dBZ')
    # Massage the index and column labels to get rid of extraneous zeros
    # Also convert column labels from strings to floats
    retr_table.index = retr_table.index.to_series().apply(np.around, decimals=4)
    retr_table.columns = [np.around(float(col), decimals=4) for col in retr_table.columns]

    dBZ_lookup_min = retr_table.index[0]
    dBZ_lookup_max = retr_table.index[-1]
    ZDR_lookup_min = retr_table.columns[0]
    ZDR_lookup_max = retr_table.columns[-1]
    dBZ_intv = retr_table.index[1] - retr_table.index[0]
    ZDR_intv = float(retr_table.columns[1]) - float(retr_table.columns[0])
    # Replace masked entries in ZH_rad and ZDR_rad with the minimum value of the lookup table
    ZH_rad = ZH_rad.filled(dBZ_lookup_min)
    ZDR_rad = ZDR_rad.filled(ZDR_lookup_min)
    # Now limit values to within lookup table limits
    ZH_rad = np.where(ZH_rad < dBZ_lookup_min, dBZ_lookup_min, ZH_rad)
    ZH_rad = np.where(ZH_rad > dBZ_lookup_max, dBZ_lookup_max, ZH_rad)
    ZDR_rad = np.where(ZDR_rad < ZDR_lookup_min, ZDR_lookup_min, ZDR_rad)
    ZDR_rad = np.where(ZDR_rad > ZDR_lookup_max, ZDR_lookup_max, ZDR_rad)
    # Round ZH and ZDR fields from the radar object to the nearest interval
    ZH_round = roundPartial(ZH_rad, dBZ_intv)
    ZDR_round = roundPartial(ZDR_rad, ZDR_intv)
    # Get the shape of the arrays for later, so we can reshape the flattened arrays of retrieved
    # values
    ZH_shape = ZH_round.shape
    ZH_flat = ZH_round.flatten()
    ZDR_flat = ZDR_round.flatten()

    for lookup_dir, retrieval_tag in zip(args.lookup_dirs, args.retrieval_tags):
        print("Retrieving using {}".format(retrieval_tag))
        for retr_varname in radar.retrieval_metadata.keys():
            print("Retrieving {} using lookup tables".format(retr_varname))
            lookup_path = os.path.join(lookup_dir, '{}.csv'.format(retr_varname))
            retr_table = pd.read_csv(lookup_path, sep=',', header=0, index_col='dBZ')
            # Round the indices and columns of the DataFrame (i.e. the dBZ values) to some sane
            # number of decimal places to facilitate using it as a lookup table. The floating point
            # precision gets in the way sometimes here. For example 56.4 is dumped out as
            # 56.4<some bunch of zeros>1
            retr_table.index = retr_table.index.to_series().apply(np.around, decimals=4)
            retr_table.columns = [np.around(float(col), decimals=4) for col in
                                  retr_table.columns]
            # Gah, for some reason DataFrame.lookup sometimes barfs on perfectly good floating point
            # values in columns, so convert them back to strings here. :rolleyes:
            retr_table.columns = [str(col) for col in retr_table.columns]
            # print(list(retr_table.index))
            # print(list(retr_table.columns))
            # Ok, now retrieve the desired retrieval variable
            # for each ZH/ZDR pair in the flattened radar sweep
            # for ZH, ZDR in zip(ZH_flat, ZDR_flat):
            #     print(ZH, ZDR)
            #     retr_val = retr_table.lookup([ZH], [str(ZDR)])
            retr_vals = retr_table.lookup(ZH_flat, ZDR_flat.astype('str'))
            # Reshape back to original shape
            retr_vals_data = retr_vals.reshape(ZH_shape)
            retr_vals_data = np.ma.masked_array(retr_vals_data, mask=full_mask)
            # Construct the dictionary. NOTE: Gotcha! Apparently have to make a copy here,
            # otherwise contents of 'data' don't get written when same dictionary is used later.
            retr_val_dict = radar.retrieval_metadata[retr_varname].copy()
            # If I don't make a copy on the previous line
            # The contents of retr_val_dict['data'] are
            # not updated after the first time. This is
            # very weird.
            retr_val_dict['data'] = retr_vals_data
            # Save the retrieved array as a new field in the radar sweep object
            retr_varname_out = '{}_{}'.format(retr_varname, retrieval_tag)
            radar_obj.add_field(retr_varname_out, retr_val_dict, replace_existing=True)
            # Now mask the retrieved values using the original combined mask of ZH and ZDR
            radar_obj.fields[retr_varname_out]['data'].mask = full_mask

    # Now save the radar object to a new CFRadial file with the added retrieved fields
    print("Saving {}".format(radar_output_path))
    pyart.io.write_cfradial(radar_output_path, radar_obj)
