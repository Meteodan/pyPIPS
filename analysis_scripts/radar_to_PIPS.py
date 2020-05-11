# radar_to_PIPS.py
#
# This script interpolates radar observations to the PIPS locations and times
import os
import sys
from glob import glob
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
description = "interpolates radar observations to the PIPS locations and times"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--average-gates', dest='average_gates', default=False, action='store_true',
                    help='Whether to average the nearest gates when interpolating to PIPS location')
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
load_radar_at_PIPS = config.radar_config_dict.get('load_radar_at_PIPS', False)
save_radar_at_PIPS = config.radar_config_dict.get('save_radar_at_PIPS', False)
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', 'NEXRAD')
radar_dir = config.radar_config_dict.get('radar_dir', None)
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory

parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# The following assumes that the same radar will be used for each PIPS in the deployment.
# TODO: make this more flexible
radar_paths = glob(radar_dir + '/*{}*.nc'.format(radar_name))
radar_dict = radar.read_sweeps(radar_paths, radar_start_timestamp,
                               radar_end_timestamp, field_names=field_names, el_req=el_req,
                               radar_type=radar_type)

# Outer file loop
for parsivel_combined_file in parsivel_combined_filelist:
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = parsivel_combined_ds.probe_name
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(np.float, geo_loc_str.strip('()').split(',')))
    rad_loc = radar.get_PIPS_loc_relative_to_radar(geo_loc, radar_dict['radarsweeplist'][0])
    radar_fields_at_PIPS_da, beam_height_da = \
        radar.interp_sweeps_to_one_PIPS(radar_name, radar_dict['radarsweeplist'], PIPS_name,
                                        rad_loc, sweeptime_list=radar_dict['sweeptimelist'],
                                        average_gates=args.average_gates)
    # Get rid of existing interpolated fields
    if '{}_at_PIPS'.format(radar_name) in parsivel_combined_ds:
        parsivel_combined_ds = parsivel_combined_ds.drop('{}_at_PIPS'.format(radar_name))
        parsivel_combined_ds = parsivel_combined_ds.drop('fields_{}'.format(radar_name))
    if '{}_beam_height_at_PIPS'.format(radar_name) in parsivel_combined_ds:
        parsivel_combined_ds = \
            parsivel_combined_ds.drop('{}_beam_height_at_PIPS'.format(radar_name))
    # Interpolate radar fields to the PIPS times
    radar_fields_at_PIPS_da = radar_fields_at_PIPS_da.interp_like(parsivel_combined_ds)
    radar_fields_at_PIPS_da.attrs['elevation_angle'] = el_req
    parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds,
                                                        radar_fields_at_PIPS_da,
                                                        name='{}_at_PIPS'.format(radar_name))
    beam_height_da = beam_height_da.interp_like(parsivel_combined_ds)
    beam_height_da.attrs['elevation_angle'] = el_req
    parsivel_combined_ds = \
        pipsio.combine_parsivel_data(parsivel_combined_ds, beam_height_da,
                                     name='{}_beam_height_at_PIPS'.format(radar_name))
    parsivel_combined_ds.close()
    # Save updated dataset back to file
    parsivel_combined_ds.to_netcdf(parsivel_combined_file)


    # if save_radar_at_PIPS:
    #     radar_at_PIPS_file_name = '{}_{}_{}_fields_at_PIPS.nc'.format(radar_name,
    #                                                                     radar_start_timestamp,
    #                                                                     radar_end_timestamp)
    #     radar_at_PIPS_path = os.path.join(PIPS_dir, radar_at_PIPS_file_name)
    #     radar.dump_radar_fields_at_PIPS_nc(radar_at_PIPS_path, radar_fields_at_PIPS_da)