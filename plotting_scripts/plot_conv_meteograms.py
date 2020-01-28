# plot_conv_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp
import pyPIPS.timemodule as tm

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
parser = argparse.ArgumentParser(description="Plots conventional meteograms from PIPS data")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')

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
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Create the directory for the meteogram plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Read in the PIPS data for the deployment
conv_df_list = []
parsivel_df_list = []
vd_matrix_da_list = []

for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype in zip(range(
        0, len(PIPS_filenames)), PIPS_filenames, PIPS_names, start_times, end_times, geo_locs,
        PIPS_types):

    tripips = (ptype == 'TriPIPS')
    PIPS_data_file_path = os.path.join(PIPS_dir, PIPS_filename)
    conv_df, parsivel_df, vd_matrix_da = pipsio.read_PIPS(PIPS_data_file_path,
                                                          start_timestamp=start_time,
                                                          end_timestamp=end_time, tripips=tripips)
    # We need the disdrometer locations. If they aren't supplied in the input control file, find
    # them from the GPS data
    if not geo_loc:
        geo_locs[index] = pipsio.get_PIPS_loc(conv_df['GPS_status'], conv_df['GPS_lat'],
                                              conv_df['GPS_lon'], conv_df['GPS_alt'])

    print("Lat/Lon/alt of {}: {}".format(PIPS_name, str(geo_locs[index])))

    # Resample the parsivel data to a longer interval if desired
    if requested_interval > 10.:
        DSD_interval = pips.check_requested_resampling_interval(requested_interval, 10.)
        vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
        parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
    else:
        DSD_interval = 10.

    conv_df_list.append(conv_df)
    parsivel_df_list.append(parsivel_df)
    vd_matrix_da_list.append(vd_matrix_da)

# Outer disdrometer (and deployment) loop
for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype, conv_df, \
    parsivel_df, vd_matrix_da in zip(range(0, len(PIPS_filenames)), PIPS_filenames, PIPS_names,
                                     start_times, end_times, geo_locs, PIPS_types, conv_df_list,
                                     parsivel_df_list, vd_matrix_da_list):

    tripips = (ptype == 'TriPIPS')

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Get times for PIPS meteogram plotting
    conv_datetimes = pips.get_conv_datetimes(conv_df)
    conv_datetimes_nums = dates.date2num(conv_datetimes)

    try:
        start_datetime = datetime.strptime(start_time, tm.timefmt3)
    except (ValueError, TypeError):
        start_datetime = conv_datetimes[0]
    try:
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
    except (ValueError, TypeError):
        end_datetime = conv_datetimes[-1]
    timelimits = [dates.date2num(start_datetime), dates.date2num(end_datetime)]
    start_time_string = start_datetime.strftime(tm.timefmt3)
    end_time_string = end_datetime.strftime(tm.timefmt3)

    # Interval for plotting in seconds.  Setting to larger intervals is useful to avoid
    # plot clutter for long datasets, but information will be lost.
    # TODO: Changing the plot interval from 1 is currently not implemented!
    plotinterval = 1
    plotintervalstr = '{:d}S'.format(int(plotinterval))
    sec_offset = conv_datetimes[0].second

    # Plot wind meteogram
    windavgintv = 60
    windgustintv = 3

    # Resample wind diagnostics array to flag as bad any time in the interval given by
    # plotinterval. Note, have to use numpy max() as a lambda function because the
    # pandas resample.max() does not propagate NaN's!
    # if pc.plot_diagnostics and (ptype == 'PIPS' or ptype == 'TriPIPS'):
    #     winddiag_resampled = conv_df['winddiag'].resample(plotintervalstr, label='right',
    #                                                         closed='right', base=sec_offset,
    #                                                         how=lambda x: utils.trymax(x.values))
    #     conv_plot_df['winddiag'] = winddiag_resampled.loc[
    #         winddiag_resampled.index.intersection(plottimeindex)]

    # # Organize a bunch of stuff in a dictionary
    # convmeteodict = {
    #     'plotinterval': plotinterval,
    #     'plottimes': conv_datetimes_nums,
    #     'windavgintv': windavgintv,
    #     'windgustintv': windgustintv,
    #     'conv_plot_df': conv_df,
    #     'xaxislimits': timelimits,
    # }

    # pm.plotconvmeteograms(index, pc, ib, convmeteodict, plot_diagnostics=False, )

    # Plot wind meteogram
    fig, ax1, ax2 = pm.plot_wind_meteogram(conv_datetimes_nums, conv_df, pc.PIPS_plotting_dict,
                                           windavgintv=windavgintv, windgustintv=windgustintv,
                                           xlimits=timelimits, ptype=ptype)
    plot_path = os.path.join(meteogram_image_dir, '{}_{}_{}_wind.png'.format(PIPS_name,
                                                                             start_time_string,
                                                                             end_time_string))
    fig.savefig(plot_path, dpi=300)

    # Plot temperature and dewpoint meteogram
    fig, ax1 = pm.plot_temperature_dewpoint_meteogram(conv_datetimes_nums, conv_df,
                                                      pc.PIPS_plotting_dict, ptype=ptype)
    plot_path = os.path.join(meteogram_image_dir, '{}_{}_{}_T_Td.png'.format(PIPS_name,
                                                                             start_time_string,
                                                                             end_time_string))
    fig.savefig(plot_path, dpi=300)

    # Plot relative humidity meteogram
    fig, ax1 = pm.plot_RH_meteogram(conv_datetimes_nums, conv_df,
                                    pc.PIPS_plotting_dict, ptype=ptype)
    plot_path = os.path.join(meteogram_image_dir, '{}_{}_{}_RH.png'.format(PIPS_name,
                                                                           start_time_string,
                                                                           end_time_string))
    fig.savefig(plot_path, dpi=300)

    # Plot pressure meteogram
    fig, ax1 = pm.plot_pressure_meteogram(conv_datetimes_nums, conv_df,
                                          pc.PIPS_plotting_dict, ptype=ptype)
    plot_path = os.path.join(meteogram_image_dir, '{}_{}_{}_pressure.png'.format(PIPS_name,
                                                                                 start_time_string,
                                                                                 end_time_string))
    fig.savefig(plot_path, dpi=300)

    if pc.PIPS_plotting_dict['plot_diagnostics']:
        # Plot battery voltage
        fig, ax1 = pm.plot_voltage_meteogram(conv_datetimes_nums, conv_df,
                                             pc.PIPS_plotting_dict, ptype=ptype)
        plot_path = os.path.join(meteogram_image_dir,
                                 '{}_{}_{}_voltage.png'.format(PIPS_name, start_time_string,
                                                               end_time_string))
        fig.savefig(plot_path, dpi=300)

        # Plot GPS speed
        try:
            fig, ax1 = pm.plot_GPS_speed_meteogram(conv_datetimes_nums, conv_df,
                                                   pc.PIPS_plotting_dict, ptype=ptype)
            plot_path = os.path.join(meteogram_image_dir,
                                     '{}_{}_{}_GPS_speed.png'.format(PIPS_name, start_time_string,
                                                                     end_time_string))
            fig.savefig(plot_path, dpi=300)
        except TypeError:
            print("No GPS speed information in file!")
