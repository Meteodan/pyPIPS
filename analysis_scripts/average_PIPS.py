# average_PIPS.py
#
# This script averages several PIPS together and saves the result in a new netCDF file
# It optionally will compute the PIPS timing relative to the arrival of a gust front
# The latter currently requires fine-tuning in the code but will be made more flexible in the future
from __future__ import annotations

import argparse
import contextlib
import os
from datetime import datetime

import numpy as np
import xarray as xr

import pyPIPS.parsivel_params as pp
from pyPIPS import utils


def adjust_time_coordinate(ds, gust_front_time):
    # Calculate the relative time and cast it to total seconds since the gust front
    relative_time = ds['time'] - gust_front_time
    relative_time = relative_time.dt.total_seconds()

    # Assign the new coordinate
    ds = ds.assign_coords(relative_time=relative_time)

    # Set the new time coordinate
    ds = ds.swap_dims({'time': 'relative_time'})

    # Drop the old time coordinate if desired
    # ds = ds.drop_vars('time')

    return ds  # noqa: RET504


def average_and_sum_for_PIPS_ds(ds, avg_vars, sum_vars, dim='PIPS_name', time_dim='time',
                                interp_na=True):

    # First handle the wind and compass variables
    result_vars = average_wind_and_compass(ds, dim=dim, time_dim=time_dim, interp_na=interp_na)

    # Now handle the rest of the variables. Note that some of these will need to/should be
    # overwritten by running the "calc_derived_params" script on the resulting dataset.
    for var in ds.data_vars:
        if interp_na:
            ds[var] = ds[var].interpolate_na(dim=time_dim)
        if var in avg_vars:
            result_vars[var] = ds[var].mean(dim=dim, skipna=True)
        elif var in sum_vars:
            result_vars[var] = ds[var].sum(dim=dim, skipna=True)

    return xr.Dataset(result_vars, coords={coord: ds[coord] for coord in ds.coords if coord != dim})


def average_wind_and_compass(ds, dim='PIPS_name', time_dim='time', interp_na=True):
    var_da_dict = {}

    if interp_na:
        ds['windspd'] = ds['windspd'].interpolate_na(dim=time_dim)
        with contextlib.suppress(KeyError):
            ds['windgust'] = ds['windgust'].interpolate_na(dim=time_dim)

    var_da_dict['windspd'] = ds['windspd'].mean(dim=dim, skipna=True)
    # We have to do the following because the conventional data doesn't have a windgust variable
    with contextlib.suppress(KeyError):
        var_da_dict['windgust'] = ds['windgust'].mean(dim=dim, skipna=True)

    # Check whether uavg and vavg already exist in the Dataset. If so, we are dealing with the
    # combined parsivel dataset (where the conventional data has been averaged to the parsivel
    # times)
    # EDIT: maybe we should just recompute the u and v components for the parsivel combined data
    # too, just to be safe?
    if "uavg" in ds.data_vars and "vavg" in ds.data_vars:
        u = ds['uavg']
        v = ds['vavg']
        if interp_na:
            u = u.interpolate_na(dim=time_dim)
            v = v.interpolate_na(dim=time_dim)
        var_da_dict['uavg'] = u.mean(dim=dim, skipna=True)
        var_da_dict['vavg'] = u.mean(dim=dim, skipna=True)
    else:
        # We are dealing with the conventional Dataset. Compute the u and v components of the wind
        # and their averages
        # First compute the u and v wind components
        u = ds['windspd'] * np.cos(np.deg2rad(-ds['winddirabs'] + 270.))
        v = ds['windspd'] * np.sin(np.deg2rad(-ds['winddirabs'] + 270.))
        if interp_na:
            u = u.interpolate_na(dim=time_dim)
            v = v.interpolate_na(dim=time_dim)
        # Compute the averages
        var_da_dict['uavg'] = u.mean(dim=dim, skipna=True)
        var_da_dict['vavg'] = v.mean(dim=dim, skipna=True)

    # Compute vector average wind speed
    var_da_dict['windspdavgvec'] = np.sqrt(var_da_dict['uavg']**2. + var_da_dict['vavg']**2.)
    # Compute vector average wind direction
    var_da_dict['winddirabs'] = (270.0 - (180. / np.pi) * np.arctan2(var_da_dict['vavg'],
                                                                     var_da_dict['uavg'])) % 360.
    var_da_dict['winddirabs'] = xr.where(var_da_dict['windspdavgvec'] == 0.,
                                         np.nan, var_da_dict['winddirabs'])

    # Compute unit average wind speed/direction
    unit_u = u / ds['windspd']
    unit_v = v / ds['windspd']
    unit_u_avg = unit_u.mean(dim=dim, skipna=True)
    unit_v_avg = unit_v.mean(dim=dim, skipna=True)

    wind_dir_unit_vec_avg = (270.0 - (180. / np.pi) * np.arctan2(unit_v_avg, unit_u_avg)) % 360.
    wind_dir_unit_vec_avg = xr.where((unit_u_avg == 0.) & (unit_v_avg == 0.), np.nan,
                                     wind_dir_unit_vec_avg)

    var_da_dict['unit_uavg'] = unit_u_avg
    var_da_dict['unit_vavg'] = unit_v_avg
    var_da_dict['winddirunitavgvec'] = wind_dir_unit_vec_avg

    # Handle the wind diagnostic variable. In this case, just take the maximum value
    var_da_dict['winddiag'] = ds['winddiag'].max(dim=dim, skipna=True)

    # Compute x- and y-components of compass direction
    x = np.cos(np.deg2rad(-ds['compass_dir'] + 270.))
    y = np.sin(np.deg2rad(-ds['compass_dir'] + 270.))

    if interp_na:
        x = x.interpolate_na(dim=time_dim)
        y = y.interpolate_na(dim=time_dim)

    # Compute averaged x- and y-components of compass direction
    x_avg = x.mean(dim=dim, skipna=True)
    y_avg = y.mean(dim=dim, skipna=True)

    # Calculate averaged compass direction
    var_da_dict['compass_dir'] = (270.0 - (180. / np.pi) * np.arctan2(y_avg, x_avg)) % 360.

    return var_da_dict


min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Averages several PIPS together and saves the result in a new netCDF file"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--gust-front-relative', dest='gust_front_relative', action='store_true',
                    help='Compute the PIPS timing relative to the arrival of a gust front')
parser.add_argument('--keep-GPS-for', dest='keep_GPS_for', type=str, default=None,
                    help=('Keep the GPS data for the specified PIPS probe (if any) for the final'
                          ' dataset'))
parser.add_argument('--window-size', dest='window_size', type=int, default=120,
                    help=('Window size for smoothing the fasttemp data to determine gust front '
                          'passage'))
parser.add_argument('--temp-drop-threshold', dest='temp_drop_threshold', type=float,
                    default=-0.0075,
                    help=('Temperature drop threshold for determining gust front passage'))
parser.add_argument('--reference-times', dest='reference_times', nargs='*', default=None,
                    help=('Reference times for the gust front calculations. If not provided, they '
                          'will be determined automatically'))

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

parsivel_ds_dict = {}
conv_ds_dict = {}

if args.reference_times is not None:
    reference_time_dict = {PIPS_name: args.reference_times[i]
                           for i, PIPS_name in enumerate(PIPS_names)}

for PIPS_name, parsivel_filepath, conv_filepath in zip(PIPS_names, parsivel_combined_filelist,
                                                       conv_filelist):
    parsivel_ds_dict[PIPS_name] = xr.load_dataset(parsivel_filepath)
    conv_ds_dict[PIPS_name] = xr.load_dataset(conv_filepath)

    # Check if we want to keep the GPS variables for one of the PIPS. This is useful when we are
    # averaging multiple collocated PIPS and there's no sense in averaging the GPS data. Otherwise
    # we drop these variables from the dataset
    GPS_varnames_parsivel = [var for var in parsivel_ds_dict[PIPS_name].data_vars if 'GPS' in var]
    GPS_varnames_conv = [var for var in conv_ds_dict[PIPS_name].data_vars if 'GPS' in var]
    GPS_vars_parsivel = []
    GPS_vars_conv = []

    if args.keep_GPS_for is not None and PIPS_name == args.keep_GPS_for:
        GPS_vars_parsivel = [parsivel_ds_dict[PIPS_name][var] for var in GPS_varnames_parsivel]
        GPS_vars_conv = [conv_ds_dict[PIPS_name][var] for var in GPS_varnames_conv]

if args.gust_front_relative and args.reference_times is None:
    # Compute the PIPS timing relative to the arrival of a gust front
    # This currently requires fine-tuning in the code but will be made more flexible in the
    # future

    # First smooth the fasttemp data to get rid of some of the noise
    window_size = args.window_size  # 2 minutes
    min_periods = 1
    temp_drop_threshold = args.temp_drop_threshold  # deg C. This is over one 1-s interval
    new_parsivel_ds_dict = {}
    new_conv_ds_dict = {}
    for PIPS_name in PIPS_names:
        fasttemp = conv_ds_dict[PIPS_name]['fasttemp']
        fasttemp_smoothed = fasttemp.rolling(time=window_size, center=True,
                                             min_periods=min_periods).mean()
        # parsivel_ds_dict[PIPS_name]['fasttemp_smoothed'] = fasttemp_smoothed

        # Take the 1st-order difference of the smoothed fasttemp data
        fasttemp_diff = fasttemp_smoothed.diff(dim='time')
        # parsivel_ds_dict[PIPS_name]['fasttemp_diff'] = fasttemp_diff

        # Now compute the gust front times based on a certain threshold of the temperature drop
        tindex = xr.where(fasttemp_diff < temp_drop_threshold,
                          True, False).argmax(dim='time').item()
        gust_front_time = conv_ds_dict[PIPS_name]['time'].isel(time=tindex).to_numpy()
        gust_front_time_parsivel = \
            parsivel_ds_dict[PIPS_name]['time'].sel(time=gust_front_time,
                                                    method='nearest').to_numpy()

        print(f"Gust front time for {PIPS_name}: {gust_front_time}")  # noqa: T201
        print(f"Gust front time for {PIPS_name} (parsivel): {gust_front_time_parsivel}")  # noqa: T201

        new_conv_ds_dict[PIPS_name] = adjust_time_coordinate(conv_ds_dict[PIPS_name],
                                                             gust_front_time)
        new_parsivel_ds_dict[PIPS_name] = adjust_time_coordinate(parsivel_ds_dict[PIPS_name],
                                                                 gust_front_time_parsivel)
elif args.reference_times is not None:
    new_parsivel_ds_dict = {}
    new_conv_ds_dict = {}
    for PIPS_name in PIPS_names:
        gust_front_time_dt = datetime.strptime(reference_time_dict[PIPS_name], "%Y%m%d%H%M%S")
        gust_front_time = conv_ds_dict[PIPS_name]['time'].sel(time=gust_front_time_dt)
        gust_front_time_parsivel = \
            parsivel_ds_dict[PIPS_name]['time'].sel(time=gust_front_time_dt, method='nearest')

        print(f"Gust front time for {PIPS_name}: {gust_front_time}")  # noqa: T201
        print(f"Gust front time for {PIPS_name} (parsivel): {gust_front_time_parsivel}")  # noqa: T201

        new_conv_ds_dict[PIPS_name] = adjust_time_coordinate(conv_ds_dict[PIPS_name],
                                                             gust_front_time)
        new_parsivel_ds_dict[PIPS_name] = adjust_time_coordinate(parsivel_ds_dict[PIPS_name],
                                                                 gust_front_time_parsivel)
else:
    new_parsivel_ds_dict = parsivel_ds_dict
    new_conv_ds_dict = conv_ds_dict

aligned_conv_ds_list = xr.align(*new_conv_ds_dict.values(), join='inner')
aggregated_conv_ds = xr.concat(aligned_conv_ds_list, dim='PIPS_name')
aggregated_conv_ds['PIPS_name'] = PIPS_names

# TODO: below may not work for straight-up averaging because the parsivel data don't necessarily
# line up in time.
aligned_parsivel_ds_list = xr.align(*new_parsivel_ds_dict.values(), join='inner')
aggregated_parsivel_ds = xr.concat(aligned_parsivel_ds_list, dim='PIPS_name')
aggregated_parsivel_ds['PIPS_name'] = PIPS_names

# Now, drop the GPS variables from the dataset. If we are keeping them from one of the PIPS, we
# will add them back in later

aggregated_conv_ds = aggregated_conv_ds.drop_vars(GPS_varnames_conv)
aggregated_parsivel_ds = aggregated_parsivel_ds.drop_vars(GPS_varnames_parsivel)

# Keep a copy of the attributes of the first DataSet. We will us a subset of these attributes
# for the average DataSet

parsivel_ds_attrs = new_parsivel_ds_dict[PIPS_names[0]].attrs
conv_ds_attrs = new_conv_ds_dict[PIPS_names[0]].attrs

# Compute the average of the aligned datasets. We actually want to compute the average of everything
# except the VD_matrix arrays. Instead, those should be summed, so we use a custom function
# for this.

# TODO: it's not as easy as this, actually. This won't properly compute the average wind speed
# and direction, for example. We need to do something similar to what we are doing when we are
# resampling to longer intervals in the PIPS_to_nc.py script. For now just do the raw average
# of the conventional dataset and all the resampled conventional data variables in the
# parsivel_combined dataset and we'll revisit this later.

wind_compass_vars = ['windspd', 'windspdavgvec', 'winddirabs', 'winddirunitavgvec', 'windgust',
                     'uavg', 'vavg', 'unit_uavg', 'unit_vavg', 'compass_dir', 'winddiag']
vars_to_average = []
vars_to_sum = []

for var in aggregated_parsivel_ds.data_vars:
    if 'VD_matrix' in var or 'pcount' in var:
        vars_to_sum.append(var)
    elif var not in wind_compass_vars:
        vars_to_average.append(var)

# Handle the wind and compass variables separately

time_dim = 'relative_time' if args.gust_front_relative else 'time'
conv_average_ds = average_and_sum_for_PIPS_ds(aggregated_conv_ds, vars_to_average, vars_to_sum,
                                              dim='PIPS_name', time_dim=time_dim)
parsivel_average_ds = average_and_sum_for_PIPS_ds(aggregated_parsivel_ds, vars_to_average,
                                                  vars_to_sum, dim='PIPS_name', time_dim=time_dim)

# conv_average_ds = aggregated_conv_ds.mean(dim='PIPS_name')
# parsivel_average_ds = aggregated_parsivel_ds.mean(dim='PIPS_name')

# Now add back the GPS variables if we are keeping them for one of the PIPS
if args.keep_GPS_for is not None:
    for var, da in zip(GPS_varnames_parsivel, GPS_vars_parsivel):
        parsivel_average_ds[var] = da
    for var, da in zip(GPS_varnames_conv, GPS_vars_conv):
        conv_average_ds[var] = da
    parsivel_average_ds.attrs['GPS_probe_name'] = args.keep_GPS_for
    conv_average_ds.attrs['GPS_probe_name'] = args.keep_GPS_for

# Add back appropriate attributes to the average DataSets
conv_average_ds.attrs['probe_name'] = 'PIPS_average'
conv_average_ds.attrs['deployment_name'] = conv_ds_attrs['deployment_name']

parsivel_average_ds.attrs['probe_name'] = 'PIPS_average'
for key in parsivel_ds_attrs:
    if key not in {'probe_name', 'parsivel_angle', 'location', 'starting_time', 'ending_time'}:
        parsivel_average_ds.attrs[key] = parsivel_ds_attrs[key]

# Add number of PIPS used in averaging as an attribute to both DataSets
conv_average_ds.attrs['num_PIPS'] = len(PIPS_names)
parsivel_average_ds.attrs['num_PIPS'] = len(PIPS_names)

# Now dump to files

conv_average_output_filename = f'conventional_raw_{dataset_name}_PIPS_avg.nc'
conv_average_output_path = os.path.join(PIPS_dir, conv_average_output_filename)
print(f"Dumping {conv_average_output_filename}")  # noqa: T201
conv_average_ds.to_netcdf(conv_average_output_path)

parsivel_average_output_filename = \
    f'parsivel_combined_{dataset_name}_PIPS_avg_{int(parsivel_average_ds.DSD_interval):d}s.nc'
parsivel_average_output_path = os.path.join(PIPS_dir, parsivel_average_output_filename)
print(f"Dumping {parsivel_average_output_filename}")  # noqa: T201
parsivel_average_ds.to_netcdf(parsivel_average_output_path)
