# plot_conv_meteograms_nc.py
#
# This script plots meteograms for conventional data from the Portable Integrated Precipitation
# Stations (PIPS)
import os
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
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
parser = argparse.ArgumentParser(description="Plots DSD meteograms from PIPS data")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--plot-dir', metavar='<path/to/plot/directory/>', dest='plot_dir',
                    default=None,
                    help='directory to store plots (overrides that in the config file')

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
conv_filenames = config.PIPS_IO_dict.get('conv_filenames_nc', parsivel_combined_filenames)
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
load_radar_at_PIPS = config.radar_config_dict.get('load_radar_at_PIPS', False)
save_radar_at_PIPS = config.radar_config_dict.get('save_radar_at_PIPS', False)
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('clean_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
plot_retrieval = config.radar_config_dict.get('plot_retrieval', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_dir = config.radar_config_dict.get('radar_dir', None)
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Set plot averaging intervals, etc.
avgwind = pc.PIPS_plotting_dict.get('avgwind', False)
windavgintv = pc.PIPS_plotting_dict.get('windavgintv', 60)
windgustintv = pc.PIPS_plotting_dict.get('windgustintv', 3)

# Create the directory for the meteogram plots if it doesn't exist
if args.plot_dir:
    plot_dir = args.plot_dir
meteogram_image_dir = os.path.join(plot_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Get a list of the conventional netCDF data files that are present in the PIPS directory
conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames]

for index, conv_file in enumerate(conv_filelist):
    print("Reading {}".format(conv_file))
    conv_ds = xr.load_dataset(conv_file)
    PIPS_name = conv_ds.probe_name
    deployment_name = deployment_names[index]  # parsivel_combined_ds.deployment_name
    ptype = PIPS_types[index]
    image_dir = os.path.join(meteogram_image_dir, deployment_name)
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)

    start_time = start_times[index]
    end_time = end_times[index]

    # Get times for PIPS meteogram plotting
    conv_datetimes = pips.get_datetimes(conv_ds)
    try:
        start_datetime = datetime.strptime(start_time, tm.timefmt3)
        print('start_datetime', start_datetime)
    except (ValueError, TypeError):
        start_datetime = conv_datetimes[0]
    try:
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
        print('end_datetime', end_datetime)
    except (ValueError, TypeError):
        end_datetime = conv_datetimes[-1]
    timelimits = [start_datetime, end_datetime]
    start_time_string = start_datetime.strftime(tm.timefmt3)
    end_time_string = end_datetime.strftime(tm.timefmt3)

    # Interval for plotting in seconds.  Setting to larger intervals is useful to avoid
    # plot clutter for long datasets, but information will be lost.
    # TODO: Changing the plot interval from 1 is currently not implemented!
    plotinterval = 1
    plotintervalstr = '{:d}S'.format(int(plotinterval))
    sec_offset = conv_datetimes[0].second

    # Plot wind meteogram

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
    fig, ax1, ax2 = pm.plot_wind_meteogram(conv_datetimes, conv_ds, pc.PIPS_plotting_dict,
                                           avgwind=avgwind, windavgintv=windavgintv,
                                           windgustintv=windgustintv, xlimits=timelimits,
                                           ptype=ptype)
    PIPS_plot_name = '{}_{}_{}_{}_wind.png'.format(PIPS_name, deployment_name, start_time_string,
                                                   end_time_string)
    plot_path = os.path.join(image_dir, PIPS_plot_name)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    # Plot temperature and dewpoint meteogram
    fig, ax1 = pm.plot_temperature_dewpoint_meteogram(conv_datetimes, conv_ds,
                                                      pc.PIPS_plotting_dict, xlimits=timelimits,
                                                      ptype=ptype)
    PIPS_plot_name = '{}_{}_{}_{}_T_Td.png'.format(PIPS_name, deployment_name, start_time_string,
                                                   end_time_string)
    plot_path = os.path.join(image_dir, PIPS_plot_name)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    # Plot relative humidity meteogram
    fig, ax1 = pm.plot_RH_meteogram(conv_datetimes, conv_ds,
                                    pc.PIPS_plotting_dict, xlimits=timelimits, ptype=ptype)

    PIPS_plot_name = '{}_{}_{}_{}_RH.png'.format(PIPS_name, deployment_name, start_time_string,
                                                 end_time_string)
    plot_path = os.path.join(image_dir, PIPS_plot_name)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    # Plot pressure meteogram
    fig, ax1 = pm.plot_pressure_meteogram(conv_datetimes, conv_ds,
                                          pc.PIPS_plotting_dict, xlimits=timelimits, ptype=ptype)
    PIPS_plot_name = '{}_{}_{}_{}_pressure.png'.format(PIPS_name, deployment_name,
                                                       start_time_string, end_time_string)
    plot_path = os.path.join(image_dir, PIPS_plot_name)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    # Plot compass direction meteogram
    fig, ax1 = pm.plot_compass_dir_meteogram(conv_datetimes, conv_ds,
                                             pc.PIPS_plotting_dict, xlimits=timelimits, ptype=ptype)
    PIPS_plot_name = '{}_{}_{}_{}_compass.png'.format(PIPS_name, deployment_name,
                                                      start_time_string, end_time_string)
    plot_path = os.path.join(image_dir, PIPS_plot_name)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

