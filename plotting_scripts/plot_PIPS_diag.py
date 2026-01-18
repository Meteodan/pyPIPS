# plot_PIPS_diag.py
#
# This script plots meteograms for diagnostic variables from the Portable Integrated Precipitation
# Stations (PIPS)
from __future__ import annotations

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import pyPIPS.parsivel_params as pp
import pyPIPS.PIPS as pips
import pyPIPS.plotmodule as pm
import pyPIPS.timemodule as tm
from pyPIPS import utils

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
parser = argparse.ArgumentParser(description="Plots diagnostic meteograms from PIPS data")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--plot-dir', metavar='<path/to/plot/directory/>', dest='plot_dir',
                    default=None,
                    help='directory to store plots (overrides that in the config file')
parser.add_argument('--plot-start-time', metavar='YYYYmmDDHHMMSS', dest='plot_start_time',
                    default=None, help='start time for plotting (overrides those in config file)')
parser.add_argument('--plot-end-time', metavar='YYYYmmDDHHMMSS', dest='plot_end_time',
                    default=None, help='end time for plotting (overrides those in config file)')
parser.add_argument('--time-dim', dest='time_dim', default='time',
                    help='Name of the time dimension in the netCDF file')
parser.add_argument('--plot-conv', dest='plot_conv', action='store_true', default=False,
                    help='Whether to plot conventional diagnostic variables')
parser.add_argument('--plot-parsivel', dest='plot_parsivel', action='store_true', default=False,
                    help='Whether to plot parsivel diagnostic variables')

args = parser.parse_args()

use_plot_date = args.time_dim != 'relative_time'

# If neither flag is set, plot both
if not args.plot_conv and not args.plot_parsivel:
    args.plot_conv = True
    args.plot_parsivel = True

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Dynamically import the plotting configuration file
utils.log(f"Plotting configuration file is {args.plot_config_path}")
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
parsivel_combined_filenames = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)
conv_filenames = config.PIPS_IO_dict.get('conv_filenames_nc', parsivel_combined_filenames)
if args.plot_start_time:
    start_times = [args.plot_start_time] * len(PIPS_names)
else:
    start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
if args.plot_end_time:
    end_times = [args.plot_end_time] * len(PIPS_names)
else:
    end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Create the directory for the meteogram plots if it doesn't exist
if args.plot_dir:
    plot_dir = args.plot_dir
meteogram_image_dir = os.path.join(plot_dir, 'meteograms_diag')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Process conventional data files if requested
if args.plot_conv and conv_filenames is not None:
    conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames]

    for index, conv_file in enumerate(conv_filelist):
        print(f"Reading {conv_file}")  # noqa: T201
        conv_ds = xr.load_dataset(conv_file)
        PIPS_name = conv_ds.probe_name
        deployment_name = deployment_names[index]
        ptype = PIPS_types[index]
        image_dir = os.path.join(meteogram_image_dir, deployment_name, 'conv')
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)

        start_time = start_times[index]
        end_time = end_times[index]

        # Get times for PIPS meteogram plotting
        if args.time_dim == 'relative_time':
            conv_datetimes = conv_ds['relative_time'].to_numpy()
            conv_datetimes = conv_datetimes.astype('timedelta64[s]').astype('int')
            start_datetime = conv_datetimes[0]
            end_datetime = conv_datetimes[-1]
            start_time_string = str(start_datetime)
            end_time_string = str(end_datetime)
        else:
            conv_datetimes = pips.get_datetimes(conv_ds)
            try:
                start_datetime = datetime.strptime(start_time, tm.timefmt3)
            except (ValueError, TypeError):
                start_datetime = conv_datetimes[0]
            try:
                end_datetime = datetime.strptime(end_time, tm.timefmt3)
            except (ValueError, TypeError):
                end_datetime = conv_datetimes[-1]
            start_time_string = start_datetime.strftime(tm.timefmt3)
            end_time_string = end_datetime.strftime(tm.timefmt3)

        timelimits = [start_datetime, end_datetime]

        # Plot GPS variables
        fig, axes = pm.plot_GPS_variables(conv_datetimes, conv_ds, pc.PIPS_plotting_dict,
                                          xlimits=timelimits, ptype=ptype,
                                          use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_GPS.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot battery voltage
        fig, ax = pm.plot_voltage_meteogram(conv_datetimes, conv_ds, pc.PIPS_plotting_dict,
                                           xlimits=timelimits, ptype=ptype,
                                           use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_battery.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot compass direction
        fig, ax = pm.plot_compass_dir_meteogram(conv_datetimes, conv_ds, pc.PIPS_plotting_dict,
                                                xlimits=timelimits, ptype=ptype,
                                                use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_compass.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot wind diagnostic if available
        if 'winddiag' in conv_ds:
            fig, ax = pm.plot_winddiag_meteogram(conv_datetimes, conv_ds, pc.PIPS_plotting_dict,
                                                 xlimits=timelimits, ptype=ptype,
                                                 use_plot_date=use_plot_date)
            PIPS_plot_name = \
                f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_winddiag.png'
            plot_path = os.path.join(image_dir, PIPS_plot_name)
            fig.savefig(plot_path, dpi=300)
            plt.close(fig)

        print(f"Conventional diagnostic plots saved to {image_dir}")  # noqa: T201

# Process parsivel data files if requested
if args.plot_parsivel and parsivel_combined_filenames is not None:
    parsivel_filelist = [os.path.join(PIPS_dir, pf) for pf in parsivel_combined_filenames]

    for index, parsivel_file in enumerate(parsivel_filelist):
        print(f"Reading {parsivel_file}")  # noqa: T201
        parsivel_ds = xr.load_dataset(parsivel_file)
        PIPS_name = parsivel_ds.probe_name
        deployment_name = deployment_names[index]
        image_dir = os.path.join(meteogram_image_dir, deployment_name, 'parsivel')
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)

        start_time = start_times[index]
        end_time = end_times[index]

        # Get times for PIPS meteogram plotting
        if args.time_dim == 'relative_time':
            parsivel_datetimes = parsivel_ds['relative_time'].to_numpy()
            parsivel_datetimes = parsivel_datetimes.astype('timedelta64[s]').astype('int')
            start_datetime = parsivel_datetimes[0]
            end_datetime = parsivel_datetimes[-1]
            start_time_string = str(start_datetime)
            end_time_string = str(end_datetime)
        else:
            parsivel_datetimes = pips.get_PSD_datetimes(parsivel_ds['VD_matrix'],
                                                       dim_name=args.time_dim)
            try:
                start_datetime = datetime.strptime(start_time, tm.timefmt3)
            except (ValueError, TypeError):
                start_datetime = parsivel_datetimes[0]
            try:
                end_datetime = datetime.strptime(end_time, tm.timefmt3)
            except (ValueError, TypeError):
                end_datetime = parsivel_datetimes[-1]
            start_time_string = start_datetime.strftime(tm.timefmt3)
            end_time_string = end_datetime.strftime(tm.timefmt3)

        timelimits = [start_datetime, end_datetime]

        # Plot particle counts
        fig, ax = pm.plot_parsivel_counts_meteogram(parsivel_datetimes, parsivel_ds,
                                                    pc.PIPS_plotting_dict, xlimits=timelimits,
                                                    use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_pcount.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot reflectivity
        fig, ax = pm.plot_parsivel_reflectivity_meteogram(parsivel_datetimes, parsivel_ds,
                                                          pc.PIPS_plotting_dict,
                                                          xlimits=timelimits,
                                                          use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_reflectivity.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot rain rate
        fig, ax = pm.plot_parsivel_rainrate_meteogram(parsivel_datetimes, parsivel_ds,
                                                      pc.PIPS_plotting_dict, xlimits=timelimits,
                                                      use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_rainrate.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot sensor temperature and voltage
        fig, axes = pm.plot_parsivel_temp_voltage_meteogram(parsivel_datetimes, parsivel_ds,
                                                            pc.PIPS_plotting_dict,
                                                            xlimits=timelimits,
                                                            use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_temp_voltage.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot signal amplitude
        fig, ax = pm.plot_parsivel_signal_amplitude_meteogram(parsivel_datetimes, parsivel_ds,
                                                              pc.PIPS_plotting_dict,
                                                              xlimits=timelimits,
                                                              use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_amplitude.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot accumulated precipitation
        fig, ax = pm.plot_parsivel_accumulation_meteogram(parsivel_datetimes, parsivel_ds,
                                                          pc.PIPS_plotting_dict,
                                                          xlimits=timelimits,
                                                          use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_accumulation.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        # Plot sample interval
        fig, ax = pm.plot_parsivel_sample_interval_meteogram(parsivel_datetimes, parsivel_ds,
                                                             pc.PIPS_plotting_dict,
                                                             xlimits=timelimits,
                                                             use_plot_date=use_plot_date)
        PIPS_plot_name = \
            f'{PIPS_name}_{deployment_name}_{start_time_string}_{end_time_string}_sample_interval.png'
        plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

        print(f"Parsivel diagnostic plots saved to {image_dir}")  # noqa: T201

print("All diagnostic plots completed!")  # noqa: T201
