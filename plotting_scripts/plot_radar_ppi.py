# plot_radar_ppi.py
#
# This script plots radar PPI plots overlaid with the locations of the PIPS
import os
import argparse
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.signal import medfilt2d
import pyart
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.utils as utils
import pyPIPS.timemodule as tm

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
parser = argparse.ArgumentParser(description="Plots radar PPIs",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--el-req', type=float, dest='el_req_cl', default=None,
                    help='Requested elevation angle (overrides value in config file)')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--plot-filtered-fields', dest='plot_filtered', default=False,
                    action='store_true',
                    help=('Whether to also plot previously filtered dBZ and ZDR fields for the '
                          'retrieval'))
help_msg = ('tag to determine filename variant for input nc files (V06 when produced by pyART, '
            'SUR when produced by RadxConvert)')
parser.add_argument('--fname-variant', dest='fname_variant', default='V06', help=help_msg)
parser.add_argument('--input-tag', dest='input_tag', default=None,
                    help='Input nametag to determine which radar files to read in')
parser.add_argument('--filter-fix', dest='filter_fix', action='store_true',
                    help='Fix to properly filter the filtered fields that were affected by the bug')
parser.add_argument('--med-filter-width', type=float, dest='med_filter_width', default=3,
                    help='Width of median filter in gates')
parser.add_argument('--dBZ-thresh', type=float, dest='dBZ_thresh', default=5.,
                    help='Threshold of reflectivity below which to exclude (dBZ)')
parser.add_argument('--RHV-thresh', type=float, dest='RHV_thresh', default=0.95,
                    help='Threshold of RHV below which to exclude')
parser.add_argument('--image-fmt', dest='image_fmt', default='png',
                    help='Image file format (i.e. png, eps, pdf)')
parser.add_argument('--use-plot-ppi-map', dest='use_plot_ppi_map', default=False,
                    help='Whether to use pyart plot_ppi_map vs. simpler, faster solution.')
parser.add_argument('--dealias-vel', dest='dealias_vel', action='store_true',
                    help='whether to dealias velocity (uses pyART region-based method)')

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

# Dynamically import the plotting configuration file
utils.log("Plotting configuration file is {}".format(args.plot_config_path))
try:
    pc = utils.import_all_from(args.plot_config_path)
    utils.log("Successfully imported pyPIPS control parameters!")
except Exception:
    utils.warning(
        "Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
    import configs.plot_config_default as pc

# Extract needed lists and variables from PIPS_IO_dict configuration dictionary
dataset_name = config.PIPS_IO_dict.get('dataset_name', None)
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
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
if 'VEL' in field_names and args.dealias_vel and 'VEL_corrected' not in field_names:
    field_names.append("VEL_corrected")
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

# Create the directory for the PPI plots if it doesn't exist
radar_ppi_image_dir = os.path.join(plot_dir, 'radar_ppi')
if not os.path.exists(radar_ppi_image_dir):
    os.makedirs(radar_ppi_image_dir)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
# TODO: change this logic to actually use the files listed in the config file
# parsivel_combined_filenames = [
#     'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
#     for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
# parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# Read radar sweeps
# First get a list of all potentially relevant radar files in the directory
if args.input_tag is None:
    radar_paths = glob(radar_dir + '/*{}*{}.nc'.format(radar_name, args.fname_variant))
else:
    radar_paths = glob(radar_dir + '/*{}*_{}.nc'.format(radar_name, args.input_tag))
# Then find only those between the requested times
radar_path_dict = radar.get_radar_paths_between_times(radar_paths, radar_start_timestamp,
                                                      radar_end_timestamp, radar_type=radar_type,
                                                      fname_format=radar_fname_pattern)
if radar_type == 'XTRRA':
    radar_path_dict = radar.get_radar_paths_single_elevation(radar_path_dict, el_req=el_req,
                                                             radar_type=radar_type)

# Read first sweep in to get radar location
first_sweep_list = radar.readCFRadial_pyART(el_req, radar_path_dict['rad_path_list'][0],
                                            compute_kdp=False)
first_sweep = first_sweep_list[0]
rlat = first_sweep.latitude['data'][0]
rlon = first_sweep.longitude['data'][0]
ralt = first_sweep.altitude['data'][0]

geo_loc_dict = {}
rad_loc_dict = {}
PIPS_names = []
PIPS_ds_dict = {}
for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = parsivel_combined_ds.probe_name
    PIPS_ds_dict[PIPS_name] = parsivel_combined_ds
    PIPS_names.append(PIPS_name)
    deployment_name = parsivel_combined_ds.deployment_name
    image_dir = os.path.join(radar_ppi_image_dir, deployment_name)
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(float, geo_loc_str.strip('()').split(',')))
    geo_loc_dict[PIPS_name] = geo_loc
    rad_loc = radar.get_PIPS_loc_relative_to_radar(geo_loc, rlat, rlon, ralt)
    rad_loc_dict[PIPS_name] = rad_loc

# Find nice bounds for the radar PPI plots to center the PIPS deployments
PIPS_x = [rad_loc_dict[PIPS_name][0] for PIPS_name in PIPS_names]
PIPS_y = [rad_loc_dict[PIPS_name][1] for PIPS_name in PIPS_names]

# buffer zone in meters surrounding PIPS for radar plot
# TODO: make these command-line arguments
buffer_x = 20000.
buffer_y = 20000.

xmin = min(PIPS_x) - buffer_x
xmax = max(PIPS_x) + buffer_x
ymin = min(PIPS_y) - buffer_y
ymax = max(PIPS_y) + buffer_y

bounds = [xmin, xmax, ymin, ymax]

for radar_path in radar_path_dict['rad_path_list']:
    radar_obj_list = radar.readCFRadial_pyART(el_req, radar_path, compute_kdp=False)
    # Loop through the matching sweeps in the file (there may be more than one because of
    # split cuts/SAILS)
    for radar_obj in radar_obj_list:

        # TODO: temporarily including optional dealiasing of velocity. Later will split
        # this off into its own script
        if args.dealias_vel and 'VEL' in field_names:
            print("Dealiasing velocities using region-based method")
            # create a gate filter which specifies gates to exclude from dealiasing
            gatefilter = pyart.filters.GateFilter(radar_obj)
            gatefilter.exclude_transition()
            gatefilter.exclude_invalid("VEL")
            gatefilter.exclude_invalid("REF")
            gatefilter.exclude_outside("REF", 0, 80)

            # perform dealiasing
            dealias_data = pyart.correct.dealias_region_based(radar_obj, gatefilter=gatefilter,
                                                              vel_field='VEL')
            radar_obj.add_field("VEL_corrected", dealias_data, replace_existing=True)

        # TODO: temporary fix below for fixing bugged filtered fields.
        # Code borrowed from filter_radar_sweep.py
        # NOTE: Can we remove this now?
        if args.filter_fix:
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
                                         radar_obj.fields[field_name]['data'].copy(),
                                         replace_existing=True)

            print("Applying median filter")
            for field_name in [ZH_name, ZDR_name, RHV_name]:
                radar_obj.fields[field_name + '_filtered']['data'] = \
                    medfilt2d(radar_obj.fields[field_name + '_filtered']['data'],
                              kernel_size=args.med_filter_width)

            print("Creating dBZ and RHV gate filter")
            # Create a gate filter to mask out areas with dBZ and RHV below thresholds
            rhoHV_ref_filter = \
                pyart.correct.moment_based_gate_filter(radar_obj, rhv_field=RHV_name + '_filtered',
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
        # Extract time from start of sweep
        sweep_time = pyart.graph.common.generate_radar_time_sweep(radar_obj, 0)
        # Convert to np.datetime64
        sweep_time_dt64 = np.datetime64(sweep_time)
        sweep_time_string = sweep_time.strftime(tm.timefmt3)

        # Figure out which PIPS (if any) were deployed at this radar time and add them to the lists
        # to plot
        PIPS_names_toplot = []
        geo_locs_toplot = []
        rad_locs_toplot = []
        geo_locs_dict_toplot = {}
        rad_locs_dict_toplot = {}

        for PIPS_name, PIPS_ds in PIPS_ds_dict.items():
            PIPS_start_time = PIPS_ds.time[0].values
            PIPS_end_time = PIPS_ds.time[-1].values

            if sweep_time_dt64 >= PIPS_start_time and sweep_time_dt64 <= PIPS_end_time:
                PIPS_names_toplot.append(PIPS_name)
                geo_locs_toplot.append(geo_loc_dict[PIPS_name])
                rad_locs_toplot.append(rad_loc_dict[PIPS_name])
                geo_locs_dict_toplot[PIPS_name] = geo_loc_dict[PIPS_name]
                rad_locs_dict_toplot[PIPS_name] = rad_loc_dict[PIPS_name]

        if args.use_plot_ppi_map:
            figlist, axlist, fields_plotted = \
                radar.plotsweep_pyART(radar_obj, sweep_time, PIPS_names_toplot, geo_locs_toplot,
                                      rad_locs_toplot, field_names,
                                      plot_filtered=args.plot_filtered, bounds=bounds)
        else:
            figlist, axlist, fields_plotted = \
                radar.plotsweep_pcolor(radar_obj, field_names, sweep_time,
                                       PIPS_names=PIPS_names_toplot,
                                       PIPS_rad_loc_dict=rad_locs_dict_toplot,
                                       plot_filtered=args.plot_filtered, bounds=bounds)

        for fig, ax, field_name in zip(figlist, axlist, fields_plotted):
            PIPS_plot_name = '{}_{}_{}_{}_{}deg.{}'.format(field_name, deployment_name,
                                                           sweep_time_string, radar_name,
                                                           str(el_req),
                                                           args.image_fmt)
            PIPS_plot_path = os.path.join(image_dir, PIPS_plot_name)
            fig.savefig(PIPS_plot_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
