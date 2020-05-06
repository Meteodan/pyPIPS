# ARPS_xyslice_to_nc.py
#
# This script dumps out intermediate netCDF files of a 2D patch surrounding the PIPS from an ARPS
# ensemble
import os
import argparse
from datetime import datetime
import numpy as np
import xarray as xr
from mpl_toolkits.basemap import Basemap
import pyPIPS.parsivel_params as pp
import pyPIPS.utils as utils
import pyPIPS.simulator as sim
from pyCRMtools.modules import utils as CRMutils
from pyCRMtools.pycaps import arps_read

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']

# Parse the command line options
description = "Reads ARPS HDF files, extracts a 2D slice, and converts to netCDF"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--model-level', metavar='<int>', type=int, dest='model_level', default=1,
                    help='Model level to extract')
parser.add_argument('--output-dir', metavar='</path/to/output/dir>', dest='output_dir',
                    help='Output directory for nc files')
parser.add_argument('--parallel', action='store_true', dest='process_parallel', default=False,
                    help='Process in parallel')
parser.add_argument('--njobs', metavar='<int>', type=int, dest='n_jobs', default=5,
                    help='Number of parallel jobs')
args = parser.parse_args()

# Create output directory if it doesn't already exist
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# Dynamically import the case configuration file
utils.log("Case config file is {}".format(args.case_config_path))
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Extract needed parameters from the model_config_dict in the configuration file
# Timestamp of model initial time (i.e. zero seconds)
timestamp_model_init = config.model_config_dict['timestamp_model_init']
datetime_model_init = datetime.strptime(timestamp_model_init, '%Y%m%d%H%M%S')
# Start and stop time of desired time window
timestamp_start = config.model_config_dict['timestamp_model_start']
timestamp_stop = config.model_config_dict['timestamp_model_stop']
datetime_start = datetime.strptime(timestamp_start, '%Y%m%d%H%M%S')
datetime_stop = datetime.strptime(timestamp_stop, '%Y%m%d%H%M%S')
# Interval in seconds for model output
tintv = config.model_config_dict['model_dt']
# Interval in seconds for ensemble mean analysis
tintv_mean = config.model_config_dict['model_dt_mean']

datetime_range = CRMutils.get_datetime_range(datetime_start, datetime_stop, tintv)
trange_sec = CRMutils.modeltimes_from_datetimes(datetime_range, datetime_start=datetime_model_init)

datetime_range_mean = CRMutils.get_datetime_range(datetime_start, datetime_stop, tintv_mean)
trange_sec_mean = CRMutils.modeltimes_from_datetimes(datetime_range_mean,
                                                     datetime_start=datetime_model_init)

filetype = 'history'
fileformat = config.model_config_dict['fileformat']
expname = config.model_config_dict['runname']
basedir = config.model_config_dict['basedirname']
num_members = config.model_config_dict['nens']
nproc_x = config.model_config_dict['nproc_x']
nproc_y = config.model_config_dict['nproc_y']

# Extract needed lists and variables from PIPS_IO_dict configuration dictionary
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Load the ARPS grid
# Get file path for grdbas file (note that call to read_grid handles the reading of the individual
# patches)
# If the grdbas file doesn't exist, fall back to a history file
member = 1  # 0 is for ensemble mean
cycle = 'posterior'
member_dir, member_prefix = sim.get_ARPS_member_dir_and_prefix(member, cycle)
member_absdir = os.path.join(basedir, expname, member_dir)
trailer = ''
grdbas_path = arps_read.get_file_path(member_absdir, member_prefix, fileformat, filetype='grdbas')
patch_x = 1
patch_y = 1
grdbas_path_test = arps_read.add_patch_number(grdbas_path, patch_x, patch_y)
if not os.path.exists(grdbas_path_test):
    print("grdbas file doesn't exist, trying a history file!")
    grdbas_path = arps_read.get_file_path(member_absdir, member_prefix, fileformat,
                                          time=trange_sec[0], filetype='history')
    grdbas_path_test = arps_read.add_patch_number(grdbas_path, patch_x, patch_y)

    # print(grdbas_path_test)
    # print(os.path.exists(grdbas_path_test))

# Read in grid information
grid_dict = arps_read.readarpsgrid(grdbas_path, nproc_x=nproc_x, nproc_y=nproc_y)
# print(grid_dict.keys())

# Get map projection information and create a Basemap instance
# TODO: convert to use cartopy!

ctrlat, ctrlon, trulat1, trulat2, trulon = arps_read.readarpsmap(grdbas_path, nproc_x=nproc_x,
                                                                 nproc_y=nproc_y)
dx = grid_dict['dx']
dy = grid_dict['dy']
nx = grid_dict['nx']
ny = grid_dict['ny']
mapwidth = nx * dx
mapheight = ny * dy
bgmap = Basemap(projection='lcc', width=mapwidth, height=mapheight, lat_1=trulat1,
                lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h',
                area_thresh=10., suppress_ticks=False)
grid_dict['bgmap'] = bgmap

# Read in PIPS data (we just need the geographic locations)
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

geo_loc_list = []
for parsivel_combined_file in parsivel_combined_filelist:
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = parsivel_combined_ds.probe_name
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(np.float, geo_loc_str.strip('()').split(',')))
    geo_loc_list.append(geo_loc)

# A bit clunky, but extract the diameter bins as a DataArray from the last parsivel_combined_ds
# Later maybe can just promote avg_diameter to a DataArray
mid_diameters_da = parsivel_combined_ds['diameter']

# Find coordinates of PIPS stations in the model
modloc_list, coord_list = sim.get_dis_locs_arps_real_grid(grid_dict, geo_loc_list)
coord_array = np.array(coord_list)
dxlist = [i[0] for i in modloc_list]
dylist = [i[1] for i in modloc_list]
xc = grid_dict['xs']
yc = grid_dict['ys']
xe = grid_dict['x']
ye = grid_dict['y']
# Set model grid limits to center on the disdrometer locations
Dxmin = min(dxlist)
Dxmax = max(dxlist)
Dymin = min(dylist)
Dymax = max(dylist)
gridlims = [Dxmin - 25000., Dxmax + 25000., Dymin - 25000., Dymax + 25000.]
ibgn = np.searchsorted(xc, gridlims[0])
iend = np.searchsorted(xc, gridlims[1]) + 1
jbgn = np.searchsorted(yc, gridlims[2])
jend = np.searchsorted(yc, gridlims[3]) + 1

# print(gridlims)
# print(ibgn, iend, jbgn, jend)

# Read in ensemble and dump slice patches to netCDF files
# TODO: put varname list in config file
varnames = ['p', 'pt', 'qv', 'u', 'v', 'qr', 'nr', 'zr']
member_list = range(1, num_members)
klvls = [args.model_level]
xc_patch = xc[ibgn:iend+1]
yc_patch = yc[jbgn:jend+1]

# Pack the args and kwargs for the call to the parallel reader wrapper function
f_args = (basedir, expname, member, cycle, fileformat, trange_sec, varnames)
f_kwargs = {
    'filetype': filetype,
    'ibgn': ibgn,
    'iend': iend,
    'jbgn': jbgn,
    'jend': jend,
    'klvls': klvls,
    'nproc_x': nproc_x,
    'nproc_y': nproc_y,
    'ncdir': args.output_dir,
    'datetime_range': datetime_range,
    'x': xc_patch,
    'y': yc_patch,
    'mid_diameters': mid_diameters_da
    }

# Set up the parallel read of all the model members
runs_list = sim.read_ARPS_ensemble(sim.read_ARPS_member_data, member_list, f_args=f_args,
                                   f_kwargs=f_kwargs, iterate_over='member',
                                   process_parallel=args.process_parallel,
                                   n_jobs=args.n_jobs, verbose=10)
