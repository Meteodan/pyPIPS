# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import argparse
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
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

# Create the directory for the meteogram plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
print(parsivel_combined_filenames)

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    fname_tag = ''
    if comp_radar:
        radar_fields_at_PIPS_da = parsivel_combined_ds['{}_at_PIPS'.format(radar_name)]
        dim_name = 'fields_{}'.format(radar_name)
        fname_tag = radar_name

    start_time = start_times[index]
    end_time = end_times[index]

    ND = parsivel_combined_ds['ND_qc']
    logND = np.log10(ND)

    # Get times for PIPS meteogram plotting
    PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    # PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_edgetimes = PSD_datetimes_dict['PSD_datetimes_edges']
    # PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])
    PSD_centertimes = PSD_datetimes_dict['PSD_datetimes_centers']

    # D0
    D0 = dsd.calc_D0_bin(ND) * 1000.
    mu_retr = parsivel_combined_ds['mu_retr'].load()
    lamda_retr = parsivel_combined_ds['lamda_retr'].load()
    D0_retr = (3.67 + mu_retr) / lamda_retr
    # D0_retr = parsivel_combined_ds['D0_retr'].load() # dsd.calc_D0_bin(ND_retr) * 1000.
    D0_rad = radar_fields_at_PIPS_da.loc[{dim_name: 'D0'}].load()
    RR_obs = parsivel_combined_ds['precipintensity'].load()
    D0_ds = xr.Dataset({'D0_obs': D0, 'D0_retr': D0_retr, 'RR_obs': RR_obs})
    D0_ds_rad = xr.Dataset({'D0_obs': D0, 'D0_rad': D0_rad, 'RR_obs': RR_obs})

    # Set up axis parameters
    try:
        start_datetime = datetime.strptime(start_time, tm.timefmt3)
        print('start_datetime', start_datetime)
    except (ValueError, TypeError):
        start_datetime = PSD_edgetimes[0]
    try:
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
        print('end_datetime', end_datetime)
    except (ValueError, TypeError):
        end_datetime = PSD_edgetimes[-1]
    timelimits = [start_datetime, end_datetime]
    start_time_string = start_datetime.strftime(tm.timefmt3)
    end_time_string = end_datetime.strftime(tm.timefmt3)

    locator = dates.MinuteLocator(byminute=[0, 15, 30, 45])
    minorlocator = dates.MinuteLocator(byminute=range(0, 60, 5))
    dateformat = '%H:%M'
    formatter = dates.DateFormatter(dateformat)

    axparamdict = {
        'majorxlocator': locator,
        'majorxformatter': formatter,
        'minorxlocator': minorlocator,
        'axeslimits': [timelimits, [0.0, 5.0]],
        'axeslabels': [None, r'D_0']
    }

    D0_dict = {
        'field': D0,
    }

    D0_retr_dict = {
        'field': D0_retr,
    }

    D0_rad_dict = {
        'field': D0_rad,
    }

    # Make the plot
    fig, ax = pm.plot_retr_timeseries(D0_dict, D0_retr_dict, D0_rad_dict, PSD_centertimes,
                                      axparamdict, name='D0')

    plot_name = '{}_{}_{}_{}_{}_{}_retr'.format(PIPS_name, deployment_name, start_time_string,
                                                end_time_string, 'D0', fname_tag)

    plot_path = os.path.join(meteogram_image_dir, plot_name)
    fig.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
