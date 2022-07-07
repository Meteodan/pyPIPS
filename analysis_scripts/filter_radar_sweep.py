# filter_radar_sweep.py
#
# This script filters radar data from CFRadial sweeps and adds the new filtered fields to the file
import argparse
from glob import glob
import numpy as np
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.utils as utils
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
parser.add_argument('--el-req', type=float, dest='el_req_cl', default=None,
                    help='Requested elevation angle (overrides value in config file)')
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
radar_fname_pattern = config.radar_config_dict('radar_config_dict', None)
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

# Read radar sweeps
# First get a list of all potentially relevant radar files in the directory
if args.input_tag is None:
    radar_paths = glob(radar_dir + '/*{}*SUR.nc'.format(radar_name))
else:
    radar_paths = glob(radar_dir + '/*{}*_{}.nc'.format(radar_name, args.input_tag))
# Then find only those between the requested times
radar_dict = radar.get_radar_paths_between_times(radar_paths, radar_start_timestamp,
                                                 radar_end_timestamp, radar_type=radar_type,
                                                 fname_format=radar_fname_pattern)
if radar_type == 'XTRRA':
    radar_dict = radar.get_radar_paths_single_elevation(radar_dict, el_req=el_req,
                                                        radar_type=radar_type)

# STOPPED HERE! Need to fix call to read_sweeps/update read_sweeps function in radar module to
# use the radar_dict as input.
radar_dict = radar.read_sweeps(radar_dict, el_req=el_req)

if args.output_tag:
    radar_output_paths = [radar_path.replace('.nc', '_{}.nc'.format(args.output_tag))
                          for radar_path in radar_dict['radarpathlist']]
else:
    radar_output_paths = radar_dict['radarpathlist']

for radar_obj, radar_output_path in zip(radar_dict['radarsweeplist'], radar_output_paths):
    print("Getting fields")
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

    print("Applying median filter")
    for field_name in [ZH_name, ZDR_name, RHV_name]:
        radar_obj.fields[field_name + '_filtered']['data'] = \
            medfilt2d(radar_obj.fields[field_name + '_filtered']['data'],
                      kernel_size=args.med_filter_width)

    print("Creating dBZ and RHV gate filter")
    # Create a gate filter to mask out areas with dBZ and RHV below thresholds
    rhoHV_ref_filter = pyart.correct.moment_based_gate_filter(radar_obj,
                                                              rhv_field=RHV_name + '_filtered',
                                                              refl_field=ZH_name + '_filtered',
                                                              min_ncp=None,
                                                              min_rhv=args.RHV_thresh,
                                                              min_refl=args.dBZ_thresh,
                                                              max_refl=None)

    print("Applying gate filter")
    # Mask fields using the dBZ/RHV mask
    for field_name in [ZH_name, ZDR_name, RHV_name]:
        radar_obj.fields[field_name + '_filtered']['data'] = \
            np.ma.masked_where(rhoHV_ref_filter.gate_excluded,
                               radar_obj.fields[field_name + '_filtered']['data'])

    # Now save the radar object to a new CFRadial file with the added filtered fields
    print("Saving {}".format(radar_output_path))
    pyart.io.write_cfradial(radar_output_path, radar_obj)