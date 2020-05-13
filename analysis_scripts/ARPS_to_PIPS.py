# ARPS_to_PIPS.py
#
# This script reads in the intermediate nc files produced by ARPS_xyslice_to_nc.py and the PIPS data
# for the given case and interpolates the model DSD and other fields to the PIPS locations.
# Also interpolates the PIPS data to the model times. Then dumps out nc files containing the
# resulting timeseries at the PIPS locations, but at the model times.

import os
import sys
import argparse
from datetime import datetime
import numpy as np
import xarray as xr
import pandas as pd
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
description = "Interpolates ARPS model output to PIPS locations, and PIPS data to model times"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--input-dir', metavar='</path/to/input/dir>', dest='input_dir',
                    help='Input directory for intermediate model nc files')
parser.add_argument('--output-dir', metavar='</path/to/output/dir>', dest='output_dir',
                    help='Output directory for timeseries nc files')
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


cycle = 'posterior'
varnames = ['p', 'pt', 'qv', 'u', 'v', 'qr', 'nr', 'zr']
member_list = range(1, num_members+1)
# Read in the ensemble members
var_ds_list = []
for member in member_list:
    member_dir, member_prefix = sim.get_ARPS_member_dir_and_prefix(member, cycle)
    fileprefix = expname + '_' + member_prefix
    ncfilename = fileprefix + '_fields.nc'
    ncfilepath = os.path.join(args.input_dir, ncfilename)
    # Open the Dataset
    var_ds = xr.open_dataset(ncfilepath)
    # Get grid information from first file
    if member == member_list[0]:
        nx = var_ds.attrs['nx_full']
        ny = var_ds.attrs['ny_full']
        dx = var_ds.attrs['dx']
        dy = var_ds.attrs['dy']
        ctrlat = var_ds.attrs['ctrlat']
        ctrlon = var_ds.attrs['ctrlon']
        trulat1 = var_ds.attrs['trulat1']
        trulat2 = var_ds.attrs['trulat2']
        trulon = var_ds.attrs['trulon']
        xe = -dx + range(nx)*dx
        ye = -dy + range(ny)*dy
        xc = xe + 0.5*dx
        yc = ye + 0.5*dy
    var_ds_list.append(var_ds)

var_full_ds = xr.concat(var_ds_list, pd.Index(range(1, num_members+1), name='member'))

grid_dict = {
    'xs': xc,
    'x': xe,
    'ys': yc,
    'y': ye,
    'dx': dx,
    'dy': dy
}

mapwidth = nx * dx
mapheight = ny * dy
bgmap = Basemap(projection='lcc', width=mapwidth, height=mapheight, lat_1=trulat1,
                lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h',
                area_thresh=10., suppress_ticks=False)
grid_dict['bgmap'] = bgmap

# Read in PIPS data
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
geo_loc_list = []
parsivel_combined_ds_list = []
PIPS_names = []
for parsivel_combined_file in parsivel_combined_filelist:
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = parsivel_combined_ds.probe_name
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(np.float, geo_loc_str.strip('()').split(',')))

    PIPS_names.append(PIPS_name)
    geo_loc_list.append(geo_loc)
    parsivel_combined_ds_list.append(parsivel_combined_ds)
DSD_interval = parsivel_combined_ds.DSD_interval

# Find coordinates of PIPS stations in the model
modloc_list, coord_list = sim.get_dis_locs_arps_real_grid(grid_dict, geo_loc_list)

# Interpolate the model variables to the PIPS locations
x_coords = [a[0] for a in modloc_list]
y_coords = [a[1] for a in modloc_list]

x_coords_da = xr.DataArray(x_coords, coords=[PIPS_names], dims=['PIPS'])
y_coords_da = xr.DataArray(y_coords, coords=[PIPS_names], dims=['PIPS'])

var_ds_interp_list = []
for deployment_name, PIPS_name in zip(deployment_names, PIPS_names):
    var_ds_interp = var_full_ds.interp(xc=x_coords_da.sel({'PIPS': PIPS_name}),
                                       yc=y_coords_da.sel({'PIPS': PIPS_name}),
                                       xe=x_coords_da.sel({'PIPS': PIPS_name}),
                                       ye=y_coords_da.sel({'PIPS': PIPS_name}))
    print("Interpolating members to {}".format(PIPS_name))
    var_ds_interp_list.append(var_ds_interp)

# Interpolate PIPS variables to the model times
var_ds_interp_new_list = []
for parsivel_combined_ds, var_ds_interp in zip(parsivel_combined_ds_list, var_ds_interp_list):
    print("Interpolating {} to model times".format(parsivel_combined_ds.probe_name))
    parsivel_combined_ds_at_model_times = parsivel_combined_ds.interp_like(var_ds_interp)
    # Model fields and PIPS fields both have 'rho'. Rename 'rho' to 'density' in the PIPS timeseries
    # to distinguish them.
    # TODO: Probably should go through all the variable names and rename them to better distinguish
    # Model and PIPS fields...
    parsivel_combined_ds_at_model_times = \
        parsivel_combined_ds_at_model_times.rename({'rho': 'density'})
    var_ds_interp = xr.merge([var_ds_interp, parsivel_combined_ds_at_model_times])
    var_ds_interp.attrs = parsivel_combined_ds_at_model_times.attrs
    print(var_ds_interp)
    var_ds_interp_new_list.append(var_ds_interp)


# Dump interpolated timeseries to netCDF files
for deployment_name, PIPS_name, var_ds_interp in zip(deployment_names, PIPS_names,
                                                     var_ds_interp_new_list):
    output_file_name = 'ARPS_at_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name,
                                                       int(DSD_interval))
    output_file_path = os.path.join(args.output_dir, output_file_name)
    print("Writing {}".format(output_file_path))
    var_ds_interp.to_netcdf(output_file_path)