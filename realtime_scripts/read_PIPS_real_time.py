# This notebook is for testing the download of PIPS data using the web API.
# It uses urllib3 and BeautifulSoup4
import os, shutil
import time
import argparse
from datetime import datetime, timedelta
import bs4
import urllib3
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dates
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyPIPS.PIPS import avg_diameter, avg_fall_bins, max_diameter, \
    min_diameter, min_fall_bins, diameter_edges, fall_bins_edges
from pyPIPS.parsivel_params import parsivel_parameters
from pyPIPS import pips_realtime as prt
from pyPIPS import thermolib as thermo
from pyPIPS import timemodule as tm
from pyPIPS import pips_realtime_config as prc
from pyPIPS import PIPS as pips
import pyPIPS.utils as utils
import pyPIPS.plotmodule as pm


# From
# https://gist.github.com/blaylockbk/1677b446bc741ee2db3e943ab7e4cabd?permalink_comment_id=3775327
def to_datetime(date):
    """
    Converts a numpy datetime64 object to a python datetime object
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
    return datetime.utcfromtimestamp(timestamp)


eff_sensor_area = parsivel_parameters['eff_sensor_area_mm2'] * 1.e-6
interval_onesec = 1.
interval_tensec = 10.

# Parse the command line options
parser = argparse.ArgumentParser(description="Reads real-time PIPS data and dumps to netCDF files")
parser.add_argument('PIPS_name', default=None,
                    help='Name of the PIPS to query and plot data for')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--base-output-dir', dest='base_output_dir',
                    default='../web/perils_realtime/',
                    help='base output directory for web content')
args = parser.parse_args()

# Dynamically import the plotting configuration file
utils.log("Plotting configuration file is {}".format(args.plot_config_path))
try:
    pc = utils.import_all_from(args.plot_config_path)
    utils.log("Successfully imported pyPIPS control parameters!")
except Exception:
    utils.warning(
        "Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
    import configs.plot_config_default as pc

# Set up directories for output/plotting and URLs for grabbing the data
netcdf_output_dir = os.path.join(args.base_output_dir, 'netcdf')
csv_output_dir = os.path.join(args.base_output_dir, 'csv')
if not os.path.exists(netcdf_output_dir):
    os.makedirs(netcdf_output_dir)
if not os.path.exists(csv_output_dir):
    os.makedirs(csv_output_dir)

# Construct the URL for the PIPS
pips_ip = prc.PIPS_PERiLS_IPs[args.PIPS_name]
base_url = f"http://{pips_ip}/?command=TableDisplay&table="
url_onesec = f"{base_url}One_Hz"
print(url_onesec)
url_tensec = f"{base_url}Ten_Hz"

http = urllib3.PoolManager()

sleep_time_check = 10  # number of seconds to sleep between attempts to check if PIPS is online
data_interval = 60.0   # number of seconds of data per file
keep_data_for = data_interval + 10  # number of seconds to grab to ensure
                                    # overlap between loops and no missing data
keep_data_for_ts = pd.Timedelta(seconds=keep_data_for)
numrecords_onesec = int(keep_data_for / interval_onesec)
numrecords_tensec = int(keep_data_for / interval_tensec)

# Main loop

starttime = time.time()
online = False
while True:
    if not online:
        try:
            content = http.request('GET', f'http://{pips_ip}/default.html', retries=False,
                                   timeout=urllib3.Timeout(connect=5.0))
            online = True
            print(f"{args.PIPS_name} is online!")
        except:
            online = False
            print(f"{args.PIPS_name} is offline!")

    if online:
        # Grab and process the requested interval of data
        time1 = time.time()
        try:
            onesec_df = prt.scrape_onesec_data(url_onesec, numrecords=numrecords_onesec)
            onesec_df = pips.calc_thermo(onesec_df)
            print(onesec_df)
        except:
            online = False
            print(f"{args.PIPS_name} is offline or problem reading data!")
            continue
        # print(onesec_df.keys())
        # onesec_df['Dewpoint'] = thermo.calTdfromRH(onesec_df['Pressure'] * 100.,
        #                                            onesec_df['SlowTemp'] + 273.15,
        #                                            onesec_df['RH'] / 100.) - 273.15
        try:
            telegram_df, spectrum_da = prt.scrape_tensec_data(url_tensec,
                                                              numrecords=numrecords_tensec)
            print(telegram_df)
            ND_da = prt.calc_ND_da(spectrum_da)
        except:
            online = False
            print(f"{args.PIPS_name} is offline or problem reading data!")
            continue

        # TODO: dumping to netcdf might be overkill here. We might be able to
        # speed things up if needed by dumping to simple CSV files. Although
        # This does seem to be fast enough.

        # First dump a basic state vector to a CSV file for ingest by SASSI
        # Resample to 1-min interval and convert to xarray Dataset
         # Resample conventional data to the parsivel times
        sec_offset = 0
        onesec_rs_df = pips.resample_conv('PIPS', data_interval, sec_offset, onesec_df)
        onesec_rs_ds = onesec_rs_df.to_xarray()
        # Now just grab the last record
        onesec_rs_ds = onesec_rs_ds.isel(logger_datetime=-1)
        timestamp = to_datetime(onesec_rs_ds['logger_datetime'].values)
        print(timestamp)
        # Now construct the dictionary of the state variables

        PIPS_SASSI_state = {}
        PIPS_SASSI_state['ID'] = args.PIPS_name
        PIPS_SASSI_state['lon'] = onesec_rs_ds['GPS_lon'].values
        PIPS_SASSI_state['lat'] = onesec_rs_ds['GPS_lat'].values
        PIPS_SASSI_state['elev'] = onesec_rs_ds['GPS_alt'].values
        PIPS_SASSI_state['FastTemp'] = onesec_rs_ds['fasttemp'].values
        PIPS_SASSI_state['SlowTemp'] = onesec_rs_ds['slowtemp'].values
        PIPS_SASSI_state['RH'] = onesec_rs_ds['RH'].values
        PIPS_SASSI_state['p'] = onesec_rs_ds['pressure'].values
        PIPS_SASSI_state['dir'] = onesec_rs_ds['winddirabs'].values
        PIPS_SASSI_state['spd'] = onesec_rs_ds['windspd'].values


        # Output to file
        prt.dump_real_time_csv_for_SASSI(PIPS_SASSI_state, csv_output_dir, args.PIPS_name,
                                         timestamp)

        # Dump onesec dataframe to netCDF file (via xarray)
        prt.dump_real_time_df_netcdf(onesec_df, netcdf_output_dir, 'onesec',
                                     args.PIPS_name, copy_to_ts_file=True)
        # Dump raw spectrum to netcdf file
        prt.dump_real_time_da_netcdf(spectrum_da, netcdf_output_dir, 'spectrum',
                                     args.PIPS_name, copy_to_ts_file=True)
        # Dump ND_da to netcdf file
        prt.dump_real_time_da_netcdf(ND_da, netcdf_output_dir, 'ND',
                                     args.PIPS_name, copy_to_ts_file=True)
        # Dump telegram data to netcdf file
        prt.dump_real_time_df_netcdf(telegram_df, netcdf_output_dir, 'telegram',
                                     args.PIPS_name, copy_to_ts_file=True)

        time2 = time.time()
        duration = time2 - time1
        print(f"That took {duration} seconds.")
        time.sleep(data_interval - ((time.time() - starttime) % data_interval))
    else:
        time.sleep(sleep_time_check)

# # Set up dictionaries to control plotting parameters

# dateformat = '%D-%H:%M'
# datelabel = 'Time (date-H:M) UTC'

# # Temperature and dewpoint
# temp_dewp_ax_params = {
#     'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
#     'majorxformatter': dates.DateFormatter(dateformat),
#     'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
#                                          interval=1),
#     'axeslimits': [None, (-5., 35.)],
#     'axeslabels': [datelabel, r'Temperature ($^{\circ}$C)']
# }

# # Wind speed and direction
# windspd_ax_params = {
#     'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
#     'majorxformatter': dates.DateFormatter(dateformat),
#     'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
#                                          interval=1),
#     'axeslimits': [None, [0.0, 25.0]],
#     'axeslabels': [datelabel, r'wind speed (m s$^{-1}$)']
# }

# winddir_ax_params = {
#     'majorylocator': ticker.MultipleLocator(45.),
#     'axeslimits': [None, [0.0, 360.0]],
#     'axeslabels': [None, r'Wind direction ($^{\circ}$C)']
# }

# pressure_ax_params = {
#     'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
#     'majorxformatter': dates.DateFormatter(dateformat),
#     'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
#                                          interval=1),
#     'majorylocator': ticker.MultipleLocator(0.5),
#     'axeslimits': [None, [940., 980.]],
#     'axeslabels': [datelabel, r'Pressure (hPa)']
# }

# # Number concentration
# log_ND_params = {
#     'type': 'pcolor',
#     'vlimits': [-1.0, 3.0],
#     'clabel': r'log[N ($m^{-3} mm^{-1}$)]'
# }

# log_ND_ax_params = {
#     'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
#     'majorxformatter': dates.DateFormatter(dateformat),
#     'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
#                                          interval=1),
#     'axeslimits': [None, [0.0, 9.0]],
#     'majorylocator': ticker.MultipleLocator(base=1.0),
#     'axeslabels': [None, 'D (mm)'],
#     'axesautofmt': False,
# }

# interval_onesec = 1.
# interval_tensec = 10.

# plot_update_interval = 60  # seconds
# plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)
# # seconds; save slightly more than an hour so that we can be sure to have
# # enough overlap to save hourly nc files that start and end at the top of
# # the hour
# keep_data_for = 3900
# keep_data_for_ts = pd.Timedelta(seconds=keep_data_for)
# numrecords_onesec = int(keep_data_for / interval_onesec)
# numrecords_tensec = int(keep_data_for / interval_tensec)

# # Grab and process the previous hour's worth of data
# onesec_df = prt.scrape_onesec_data(url_onesec, numrecords=numrecords_onesec)
# print(onesec_df.keys())
# onesec_df['Dewpoint'] = thermo.calTdfromRH(onesec_df['Pressure'] * 100.,
#                                            onesec_df['SlowTemp'] + 273.15,
#                                            onesec_df['RH'] / 100.) - 273.15
# telegram_df, spectrum_da = scrape_tripips_tensec_data(url_tensec, numrecords=numrecords_tensec)
# ND_da = calc_ND_da(spectrum_da)

# # DTD 06/19/2021: getting some permission denied errors here for some reason. Need to find out
# # why but for now disable all dumping of netCDF files
# # EDIT: found the problem. The issue was that cron was not given full disk access. I followed
# # the instructions here to fix the problem: https://blog.bejarano.io/fixing-cron-jobs-in-mojave/

# # Dump onesec dataframe to netCDF file (via xarray)
# netcdf_filename = 'onesec_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_base_dir, netcdf_filename)
# onesec_ds = onesec_df.to_xarray()
# onesec_ds.to_netcdf(netcdf_path)
# onesec_ds.close()
# # Dump raw spectrum to netcdf file
# netcdf_filename = 'spectrum_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_base_dir, netcdf_filename)
# spectrum_ds = spectrum_da.to_dataset(name='spectrum')
# spectrum_ds.to_netcdf(netcdf_path)
# spectrum_ds.close()
# # Dump ND_da to netcdf file
# netcdf_filename = 'ND_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_base_dir, netcdf_filename)
# ND_ds = ND_da.to_dataset(name='ND')
# ND_ds.to_netcdf(netcdf_path)
# ND_ds.close()
# # Dump telegram data to netcdf file
# netcdf_filename = 'tensec_telegram_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_base_dir, netcdf_filename)
# telegram_ds = telegram_df.to_xarray()
# telegram_ds.to_netcdf(netcdf_path)
# telegram_ds.close()

# # Check if we are near the top of the hour and save top-of-the-hour netcdf files if we are
# last_timestamp_onesec = onesec_df.index[-1]
# last_timestamp_tensec = telegram_df.index[-1]
# last_timestamp_60sec = telegram_df.index[-6]
# oldest_timestamp_onesec = last_timestamp_onesec - keep_data_for_ts
# last_hour_timestamp = last_timestamp_onesec.replace(microsecond=0, second=0, minute=0)
# previous_hour_timestamp = last_hour_timestamp - timedelta(hours=1)
# diff = (last_timestamp_onesec - last_hour_timestamp).total_seconds()

# # If we are within 5 minutes past the top of the hour, see if we need to create the
# # hourly dump of data
# if diff > 0 and diff < 300:
#     last_hour_string = last_hour_timestamp.strftime('%Y%m%d%H%M%S')
#     previous_hour_string = previous_hour_timestamp.strftime('%Y%m%d%H%M%S')

#     # Save the files to a subdirectory structure of YYYY/MM/
#     year_string = previous_hour_timestamp.strftime('%Y')
#     month_string = previous_hour_timestamp.strftime('%m')
#     netcdf_output_dir = os.path.join(netcdf_output_base_dir, f'{year_string}/{month_string}')
#     if not os.path.exists(netcdf_output_dir):
#         os.makedirs(netcdf_output_dir)

#     onesec_hourly_filename = 'onesec_{}.nc'.format(previous_hour_string)
#     onesec_hourly_path = os.path.join(netcdf_output_dir, onesec_hourly_filename)
#     spectrum_hourly_filename = 'spectrum_{}.nc'.format(previous_hour_string)
#     spectrum_hourly_path = os.path.join(netcdf_output_dir, spectrum_hourly_filename)
#     ND_hourly_filename = 'ND_{}.nc'.format(previous_hour_string)
#     ND_hourly_path = os.path.join(netcdf_output_dir, ND_hourly_filename)
#     telegram_hourly_filename = 'tensec_{}.nc'.format(previous_hour_string)
#     telegram_hourly_path = os.path.join(netcdf_output_dir, telegram_hourly_filename)
#     # First see if the file already exists
#     if not os.path.exists(onesec_hourly_path):
#         onesec_df_dump = onesec_df.loc[previous_hour_timestamp:last_hour_timestamp]
#         spectrum_da_dump = spectrum_da.loc[previous_hour_timestamp:last_hour_timestamp]
#         ND_da_dump = ND_da.loc[previous_hour_timestamp:last_hour_timestamp]
#         telegram_df_dump = telegram_df.loc[previous_hour_timestamp:last_hour_timestamp]

#         onesec_ds_dump = onesec_df_dump.to_xarray()
#         onesec_ds_dump.to_netcdf(onesec_hourly_path)
#         onesec_ds_dump.close()
#         telegram_ds_dump = telegram_df_dump.to_xarray()
#         telegram_ds_dump.to_netcdf(telegram_hourly_path)
#         telegram_ds_dump.close()
#         spectrum_ds_dump = spectrum_da_dump.to_dataset(name='spectrum')
#         spectrum_ds_dump.to_netcdf(spectrum_hourly_path)
#         spectrum_ds_dump.close()
#         ND_ds_dump = ND_da_dump.to_dataset(name='spectrum')
#         ND_ds_dump.to_netcdf(ND_hourly_path)
#         ND_ds_dump.close()

# # Create the figures
# fig, ax = plt.subplots()
# fig_vd, ax_vd = plt.subplots()
# fig_dsd, ax_dsd = plt.subplots()
# fig_t_td, ax_t_td = plt.subplots()
# fig_wind, ax_windspd = plt.subplots()
# ax_winddir = ax_windspd.twinx()
# fig_pressure, ax_pressure = plt.subplots()

# plottimes_onesec = [onesec_df.index.to_pydatetime()]
# # Temperature and Dewpoint
# Tmin = np.nanmin(onesec_df['Dewpoint'].values)
# Tmax = np.nanmax(onesec_df['SlowTemp'].values)
# fields_to_plot_onesec = [onesec_df['SlowTemp'].values, onesec_df['Dewpoint'].values]
# temp_params = pm.temp_params.copy()
# dewpoint_params = pm.dewpoint_params.copy()
# temp_params['plotmin'] = Tmin - 5.0
# dewpoint_params['plotmin'] = Tmin - 5.0
# field_parameters_onesec = [temp_params, dewpoint_params]
# ax_t_td = pm.plotmeteogram(
#     ax_t_td,
#     plottimes_onesec,
#     fields_to_plot_onesec,
#     field_parameters_onesec)
# temp_dewp_ax_params['axeslimits'] = [[plottimes_onesec[0][0], plottimes_onesec[0][-1]],
#                                      [Tmin - 5.0, Tmax + 5.0]]
# ax_t_td, = pm.set_meteogram_axes([ax_t_td], [temp_dewp_ax_params])
# # Wind speed and direction
# ax_windspd = pm.plotmeteogram(
#     ax_windspd, plottimes_onesec, [
#         onesec_df['WS_ms'].values], [
#         pm.windspeed_params])
# ax_winddir = pm.plotmeteogram(
#     ax_winddir, plottimes_onesec, [
#         onesec_df['WindDir'].values], [
#         pm.winddir_params])
# windspd_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
# winddir_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
# ax_windspd, ax_winddir = pm.set_meteogram_axes(
#     [ax_windspd, ax_winddir], [windspd_ax_params, winddir_ax_params])
# # Pressure
# pmin = np.nanmin(onesec_df['Pressure'].values)
# pmax = np.nanmax(onesec_df['Pressure'].values)
# pressure_ax_params['axeslimits'] = [[plottimes_onesec[0][0], plottimes_onesec[0][-1]],
#                                     [pmin - 2.5, pmax + 2.5]]
# fields_to_plot_press = [onesec_df['Pressure'].values]
# field_parameters_press = [pm.pressure_params]
# ax_pressure = pm.plotmeteogram(
#     ax_pressure,
#     plottimes_onesec,
#     fields_to_plot_press,
#     field_parameters_press)
# ax_pressure, = pm.set_meteogram_axes([ax_pressure], [pressure_ax_params])

# # DSD plots
# plottimes_tmp = ND_da['time'].to_index().to_pydatetime()
# # Prepend an additional at the beginning of the array so that pcolor sees this as the
# # edges of the DSD intervals.
# plottimes = np.insert(plottimes_tmp, 0, plottimes_tmp[0] - timedelta(seconds=10))
# plottimes = [plottimes]
# # Do the same at the end of the "min_diameter" and "min_fall_bins" arrays
# # TODO: apparently this is a problem in the main DSD meteogram plotting script as well, so
# # go back and fix that soon! NOTE: this is now done in PIPS.py and stored in "diameter_edges"
# # and "fall_bins_edges"

# ND_arr = ND_da.values.T
# logND_arr = np.ma.log10(ND_arr)
# fields_to_plot = [logND_arr]
# field_parameters = [log_ND_params]
# ax = pm.plotmeteogram(ax, plottimes, fields_to_plot, field_parameters,
#                       yvals=[diameter_edges] * len(fields_to_plot))
# ax, = pm.set_meteogram_axes([ax], [log_ND_ax_params])
# ax_vd.set_title(
#     r'$v_T$ vs. $D$ for time {0} to {1}'.format(
#         (last_timestamp_60sec - timedelta(seconds=10)).strftime(tm.timefmt2),
#         last_timestamp_tensec.strftime(tm.timefmt2)))
# # spectrum_da.loc[last_timestamp_tensec]
# spectrum = spectrum_da.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).sum(dim='time')
# countsplot = np.ma.masked_where(spectrum.values <= 0, spectrum)
# C = ax_vd.pcolormesh(diameter_edges, fall_bins_edges, countsplot, vmin=1, vmax=50, edgecolors='w')
# divider = make_axes_locatable(ax_vd)
# cax = divider.append_axes("right", size="5%")
# cb = fig_vd.colorbar(C, cax=cax, orientation='vertical')
# cax.set_ylabel('# of drops')
# ax_vd.set_xlim(0.0, 10.0)
# ax_vd.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
# ax_vd.set_xlabel('diameter (mm)')
# ax_vd.set_ylim(0.0, 10.0)
# ax_vd.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
# ax_vd.set_ylabel('fall speed (m/s)')

# # Plot last 60-s DSD
# # Bugfix 03/12/22: changed from sum to mean in following line!
# ND_60s_da = ND_da.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).mean(dim='time')
# ND_60s = ND_60s_da.values
# ax_dsd.set_title('DSDs for time {0} to {1}'.format(
#     (last_timestamp_60sec - timedelta(seconds=10)).strftime(tm.timefmt2),
#     last_timestamp_tensec.strftime(tm.timefmt2)))
# ax_dsd.bar(min_diameter, ND_60s * 1000.0, max_diameter - min_diameter, 10.**2., align='edge',
#            log=True, color='tan', edgecolor='k')
# ax_dsd.set_xlim(0.0, 9.0)
# ax_dsd.set_ylim(10.**2., 10.**8.5)
# ax_dsd.set_xlabel('D (mm)')
# ax_dsd.set_ylabel(r'N(D) $(m^{-4})$')

# # fig.canvas.draw()
# # Try-except block to test if the mount is working and to reset if it isn't. This is a terrible
# # hack that I need to clean up later.
# try:
#     fig.savefig(os.path.join(image_output_dir, 'logND_current.png'), dpi=300)
# except:
#     print("that didn't work")
#     # Assume mount isn't working
#     try:
#         print("unmounting stormlab")
#         os.system("/Users/dawson29/scripts/stormlabunmount")
#     except:
#         print("mounting stormlab")
#         os.system("/Users/dawson29/scripts/stormlabmount")
#     finally:
#         print("mounting stormlab")
#         os.system("/Users/dawson29/scripts/stormlabmount")
# fig.savefig(os.path.join(image_output_dir, 'logND_current.png'), dpi=300)
# fig_vd.savefig(os.path.join(image_output_dir, 'VD_current.png'), dpi=300)
# fig_dsd.savefig(os.path.join(image_output_dir, 'DSD_current.png'), dpi=300)
# fig_t_td.savefig(os.path.join(image_output_dir, 'T_Td_current.png'), dpi=300)
# fig_wind.savefig(os.path.join(image_output_dir, 'wind_current.png'), dpi=300)
# fig_pressure.savefig(os.path.join(image_output_dir, 'pressure.png'), dpi=300)
