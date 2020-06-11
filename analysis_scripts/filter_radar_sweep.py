# pyPIPS_meteograms.py
#
# This script filters radar data from CFRadial sweeps and adds the new filtered fields to the file
import os
import argparse
from glob import glob
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
import pyart
from scipy.signal import medfilt2d


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
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--dBZ-thresh', type=float, dest='dBZ_thresh', default=5.,
                    help='Threshold of reflectivity below which to exclude (dBZ)')
parser.add_argument('--RHV-thresh', type=float, dest='RHV_thresh', default=0.95,
                    help='Threshold of RHV below which to exclude')
parser.add_argument('--med-filter-width', type=float, dest='med_filter_width', default=3,
                    help='Width of median filter in gates')
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
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', None)
radar_dir = config.radar_config_dict.get('radar_dir', None)
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Read radar sweeps
if args.input_tag is None:
    radar_paths = glob(radar_dir + '/*{}*.nc'.format(radar_name))
else:
    radar_paths = glob(radar_dir + '/*{}*_{}.nc'.format(radar_name, args.input_tag))
radar_dict = radar.read_sweeps(radar_paths, radar_start_timestamp,
                               radar_end_timestamp, field_names=field_names, el_req=el_req,
                               radar_type=radar_type)

if args.output_tag:
    radar_output_paths = [radar_path.replace('.nc', '_{}.nc'.format(args.output_tag))
                          for radar_path in radar_dict['radarpathlist']]
else:
    radar_output_paths = radar_dict['radarpathlist']

for radar_obj, radar_output_path in zip(radar_dict['radarsweeplist'], radar_output_paths):
    print("Getting fields")
    # Get the ZH and ZDR fields from the radar object
    ZH_rad_tuple = radar.get_field_to_plot(radar_obj, radar.REF_aliases)
    ZDR_rad_tuple = radar.get_field_to_plot(radar_obj, radar.ZDR_aliases)
    RHV_rad_tuple = radar.get_field_to_plot(radar_obj, radar.RHV_aliases)
    ZH_name = ZH_rad_tuple[0]
    ZDR_name = ZDR_rad_tuple[0]
    RHV_name = RHV_rad_tuple[0]

    print("Creating dBZ and RHV gate filter")
    # Create a gate filter to mask out areas with dBZ and RHV below thresholds
    rhoHV_ref_filter = pyart.correct.moment_based_gate_filter(radar_obj,
                                                              rhv_field=RHV_name,
                                                              refl_field=ZH_name,
                                                              min_ncp=None,
                                                              min_rhv=args.RHV_thresh,
                                                              min_refl=args.dBZ_thresh,
                                                              max_refl=None)

    # Make copies of the ZH and ZDR fields
    radar_obj.add_field_like(ZH_name, ZH_name+'_filtered',
                             radar_obj.fields[ZH_name]['data'].copy(), replace_existing=True)
    radar_obj.add_field_like(ZDR_name, ZDR_name+'_filtered',
                             radar_obj.fields[ZDR_name]['data'].copy(), replace_existing=True)

    print("Applying gate filter")
    # Mask the ZH and ZDR fields using the dBZ/RHV mask
    radar_obj.fields[ZH_name+'_filtered']['data'] = \
        np.ma.masked_where(rhoHV_ref_filter.gate_excluded,
                           radar_obj.fields[ZH_name+'_filtered']['data'])
    radar_obj.fields[ZDR_name+'_filtered']['data'] = \
        np.ma.masked_where(rhoHV_ref_filter.gate_excluded,
                           radar_obj.fields[ZDR_name+'_filtered']['data'])

    print("Applying median filter")
    # Apply the median filter
    radar_obj.fields[ZH_name+'_filtered']['data'] = \
        medfilt2d(radar_obj.fields[ZH_name+'_filtered']['data'], kernel_size=args.med_filter_width)
    radar_obj.fields[ZH_name+'_filtered']['data'] = \
        np.ma.masked_array(data=radar_obj.fields[ZH_name+'_filtered']['data'],
                           mask=radar_obj.fields[ZH_name]['data'].mask)

    radar_obj.fields[ZDR_name+'_filtered']['data'] = \
        medfilt2d(radar_obj.fields[ZDR_name+'_filtered']['data'], kernel_size=args.med_filter_width)
    radar_obj.fields[ZDR_name+'_filtered']['data'] = \
        np.ma.masked_array(data=radar_obj.fields[ZDR_name+'_filtered']['data'],
                           mask=radar_obj.fields[ZDR_name]['data'].mask)

    # Now save the radar object to a new CFRadial file with the added filtered fields
    print("Saving {}".format(radar_output_path))
    pyart.io.write_cfradial(radar_output_path, radar_obj)