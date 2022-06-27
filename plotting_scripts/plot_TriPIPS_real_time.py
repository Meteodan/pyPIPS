# This notebook is for testing the download of PIPS data using the web API.
# It uses urllib3 and BeautifulSoup4
import os
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
from pyPIPS import thermolib as thermo
from pyPIPS import timemodule as tm
import pyPIPS.plotmodule as pm

eff_sensor_area = parsivel_parameters['eff_sensor_area_mm2'] * 1.e-6


# Function definitions. These will eventually go into their own module
def scrape_tripips_onesec_data(url, numrecords=3600):
    """Grabs one-second records from the TriPIPS http server. Uses urlib3 and beautifulsoup4"""
    http = urllib3.PoolManager()
    content = http.request('GET', url + '&records={:d}'.format(numrecords))
    soup = bs4.BeautifulSoup(content.data, "lxml")
    table = soup.find('table')
    rows = table.find_all('tr')
    headers = rows.pop(0).find_all('th')
    headers.pop(0)
    headers = [header.string for header in headers]
    data = []
    timestamps = []
    for row in rows:
        tokens = row.find_all('td')
        tokens = [token.string.strip(' "') for token in tokens]
        timestamp = tokens.pop(0)
        tokens = [np.float64(token) if token not in ('' or 'NAN') else np.nan for token in tokens]
        data.append(tokens)
        timestamps.append(timestamp)
    index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S')
    df = pd.DataFrame(data, columns=headers, index=index)
    return df


def scrape_tripips_tensec_data(url, numrecords=360):
    """Grabs ten-second records from the TriPIPS http server. Uses urlib3 and beautifulsoup4"""
    http = urllib3.PoolManager()
    content = http.request('GET', url + '&records={:d}'.format(numrecords))
    soup = bs4.BeautifulSoup(content.data, "lxml")
    table = soup.find('table')
    rows = table.find_all('tr')
    headers = rows.pop(0).find_all('th')
    headers.pop(0)
    headers = [header.string for header in headers]
    telegrams = []
    spectrum_list = []
    timestamps = []
    for row in rows:
        tokens = row.find_all('td')
        tokens = [token.string.strip(' "') for token in tokens]
        timestamp = tokens.pop(0)
        tokens.pop(0)
        parsivel_string = tokens.pop(0)
        try:
            telegram_dict = read_parsivel_telegram(parsivel_string)
            telegrams.append(telegram_dict)
            spectrum = read_parsivel_spectrum(parsivel_string)
            spectrum_list.append(spectrum)
            timestamps.append(timestamp)
        except ValueError:
            print("Bad Parsivel string detected. Ignoring this record!")

    index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S')
    df_telegram = pd.DataFrame(telegrams, index=index)
    # print(np.array(spectrum_list).shape)
    da_spectrum = xr.DataArray(spectrum_list, coords=[index, avg_fall_bins,
                                                      avg_diameter],
                               dims=['time', 'velocity', 'diameter'])
    return df_telegram, da_spectrum


def read_parsivel_telegram(parsivel_string):
    """
    Reads the parsivel telegram and returns a dictionary for each component,
    except for the spectrum, which is handled separately by read_spectrum.
    """
    # print(parsivel_string)
    parsivel_tokens = parsivel_string.strip(' "').split(';')
    # print(parsivel_tokens)
    parsivel_tokens = parsivel_tokens[:11]
    parsivel_telegram_dict = {
        'parsivel_id': int(parsivel_tokens[0]),
        'rain rate (mm per hr)': float(parsivel_tokens[1]),
        'rain accumulation (mm)': float(parsivel_tokens[2]),
        'radar reflectivity (dBZ)': float(parsivel_tokens[3]),
        'sample interval': int(parsivel_tokens[4]),
        'signal amplitude': int(parsivel_tokens[5]),
        'particle count': int(parsivel_tokens[6]),
        'sensor temp': int(parsivel_tokens[7]),
        'power supply voltage': float(parsivel_tokens[8]),
        'sensor time': parsivel_tokens[9],
        'sensor date': parsivel_tokens[10]
    }
    return parsivel_telegram_dict


def read_parsivel_spectrum(parsivel_string):
    """Given a raw string of Parsivel data, extract the DSD spectrum from it."""
    parsivel_tokens = parsivel_string.strip(' "').split(';')
    # parsivel_tokens = parsivel_tokens.strip('"')
    spectrum = [float(x) if ('' or '\n' or '\r') not in x else 0 for x in parsivel_tokens[11:]]
    try:
        spectrum = [float(x) if ('' or '\n' or '\r') not in x else 0 for x in parsivel_tokens[11:]]
    except ValueError:
        spectrum = [-999 for i in range(1025)]
    # Strip off bogus final value (why is it there?)
    if len(spectrum) == 1025:
        spectrum = spectrum[:-1]
    # Reshape the spectrum to a 32x32 matrix of integers
    spectrum = np.array(spectrum, dtype='int')
    # Assert that the size of the array is what we expect
    # Otherwise generate a spectrum of missing values (-999)
    if spectrum.size == 1024:
        spectrum = spectrum.reshape((32, 32))
    else:
        spectrum = -999 * np.ones((32, 32), dtype='int')
    return spectrum


def timestamp2datetime(timestamp):
    """Construct a datetime object from a timestamp of the form YYYY-MM-DD HH:MM:SS.XXX"""
    date, time = timestamp.strip().split()
    year = np.int(date[:4])
    month = np.int(date[5:7])
    day = np.int(date[8:])
    hour = np.int(time[:2])
    minute = np.int(time[3:5])
    sec = np.int(time[6:8])

    return datetime(year, month, day, hour, minute, sec)


def calc_ND_da(spectrum_da, interval=10, use_measured_fs=True):
    """Computes the number concentration from the 32x32 spectrum"""
    # index = spectrum_da.coords['time'].values
    if not use_measured_fs:
        raise NotImplementedError('Not yet implemented: use measured fall speed for now!')
    else:
        ND_da = spectrum_da.groupby('time').apply(calc_ND, interval=interval,
                                                  use_measured_fs=use_measured_fs)

    # ND_df = pd.DataFrame(ND_arr, columns=avg_diameter, index=index)
    return ND_da


def calc_ND_list(spectrum_list, interval=10, use_measured_fs=True):
    ND_list = [calc_ND(spectrum, interval=interval) for spectrum in spectrum_list]
    ND_arr = np.array(ND_list)
    return ND_arr


def calc_ND(spectrum, interval=10, use_measured_fs=True):
    _, vspectrum = np.meshgrid(avg_diameter, avg_fall_bins)
    dspectrum = spectrum
    if np.all(dspectrum == -999):
        # print('Here!')
        ND = -999. * np.ones_like(avg_diameter)
        ND = np.ma.masked_where(ND == -999., ND)
    else:
        ND = dspectrum / (vspectrum * interval * eff_sensor_area * (max_diameter - min_diameter))
        if use_measured_fs:
            # ND = ND.sum(axis=0)
            ND = ND.sum(dim='velocity')
    return ND


# Set up directories for output/plotting and URLs for grabbing the data
# base_output_dir = '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/webdata/'
base_output_dir = '/Users/dawson29/sshfs_mounts/stormlab_web/'
image_output_dir = os.path.join(base_output_dir, 'images')
# netcdf_output_dir = '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/netcdf_web/'
netcdf_output_dir = '/Volumes/scr_fast/Projects/TriPIPS/netcdf_web/'
if not os.path.exists(netcdf_output_dir):
    os.makedirs(netcdf_output_dir)
base_url = "http://10.163.29.26/?command=TableDisplay&table="
onesec_table = "One_Hz"
tensec_table = "Ten_Hz"
url_onesec = base_url + onesec_table
url_tensec = base_url + tensec_table

# Set up dictionaries to control plotting parameters

dateformat = '%D-%H:%M'
datelabel = 'Time (date-H:M) UTC'

# Temperature and dewpoint
temp_dewp_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'axeslimits': [None, (-5., 35.)],
    'axeslabels': [datelabel, r'Temperature ($^{\circ}$C)']
}

# Wind speed and direction
windspd_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'axeslimits': [None, [0.0, 25.0]],
    'axeslabels': [datelabel, r'wind speed (m s$^{-1}$)']
}

winddir_ax_params = {
    'majorylocator': ticker.MultipleLocator(45.),
    'axeslimits': [None, [0.0, 360.0]],
    'axeslabels': [None, r'Wind direction ($^{\circ}$C)']
}

pressure_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'majorylocator': ticker.MultipleLocator(0.5),
    'axeslimits': [None, [940., 980.]],
    'axeslabels': [datelabel, r'Pressure (hPa)']
}

# Number concentration
log_ND_params = {
    'type': 'pcolor',
    'vlimits': [-1.0, 3.0],
    'clabel': r'log[N ($m^{-3} mm^{-1}$)]'
}

log_ND_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'axeslimits': [None, [0.0, 9.0]],
    'majorylocator': ticker.MultipleLocator(base=1.0),
    'axeslabels': [None, 'D (mm)'],
    'axesautofmt': False,
}

interval_onesec = 1.
interval_tensec = 10.

plot_update_interval = 60  # seconds
plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)
# seconds; save slightly more than an hour so that we can be sure to have
# enough overlap to save hourly nc files that start and end at the top of
# the hour
keep_data_for = 3900
keep_data_for_ts = pd.Timedelta(seconds=keep_data_for)
numrecords_onesec = int(keep_data_for / interval_onesec)
numrecords_tensec = int(keep_data_for / interval_tensec)

# Grab and process the previous hour's worth of data
onesec_df = scrape_tripips_onesec_data(url_onesec, numrecords=numrecords_onesec)
print(onesec_df.keys())
onesec_df['Dewpoint'] = thermo.calTdfromRH(onesec_df['Pressure'] * 100.,
                                           onesec_df['SlowTemp'] + 273.15,
                                           onesec_df['RH'] / 100.) - 273.15
telegram_df, spectrum_da = scrape_tripips_tensec_data(url_tensec, numrecords=numrecords_tensec)
ND_da = calc_ND_da(spectrum_da)

# DTD 06/19/2021: getting some permission denied errors here for some reason. Need to find out
# why but for now disable all dumping of netCDF files
# EDIT: found the problem. The issue was that cron was not given full disk access. I followed
# the instructions here to fix the problem: https://blog.bejarano.io/fixing-cron-jobs-in-mojave/

# Dump onesec dataframe to netCDF file (via xarray)
netcdf_filename = 'onesec_current_hour.nc'
netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
onesec_ds = onesec_df.to_xarray()
onesec_ds.to_netcdf(netcdf_path)
onesec_ds.close()
# Dump raw spectrum to netcdf file
netcdf_filename = 'spectrum_current_hour.nc'
netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
spectrum_ds = spectrum_da.to_dataset(name='spectrum')
spectrum_ds.to_netcdf(netcdf_path)
spectrum_ds.close()
# Dump ND_da to netcdf file
netcdf_filename = 'ND_current_hour.nc'
netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
ND_ds = ND_da.to_dataset(name='ND')
ND_ds.to_netcdf(netcdf_path)
ND_ds.close()
# Dump telegram data to netcdf file
netcdf_filename = 'tensec_telegram_current_hour.nc'
netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
telegram_ds = telegram_df.to_xarray()
telegram_ds.to_netcdf(netcdf_path)
telegram_ds.close()

# Check if we are near the top of the hour and save top-of-the-hour netcdf files if we are
last_timestamp_onesec = onesec_df.index[-1]
last_timestamp_tensec = telegram_df.index[-1]
last_timestamp_60sec = telegram_df.index[-6]
oldest_timestamp_onesec = last_timestamp_onesec - keep_data_for_ts
last_hour_timestamp = last_timestamp_onesec.replace(microsecond=0, second=0, minute=0)
previous_hour_timestamp = last_hour_timestamp - timedelta(hours=1)
diff = (last_timestamp_onesec - last_hour_timestamp).total_seconds()

# If we are within 5 minutes past the top of the hour, see if we need to create the
# hourly dump of data
if diff > 0 and diff < 300:
    last_hour_string = last_hour_timestamp.strftime('%Y%m%d%H%M%S')
    previous_hour_string = previous_hour_timestamp.strftime('%Y%m%d%H%M%S')
    onesec_hourly_filename = 'onesec_{}.nc'.format(previous_hour_string)
    onesec_hourly_path = os.path.join(netcdf_output_dir, onesec_hourly_filename)
    spectrum_hourly_filename = 'spectrum_{}.nc'.format(previous_hour_string)
    spectrum_hourly_path = os.path.join(netcdf_output_dir, spectrum_hourly_filename)
    ND_hourly_filename = 'ND_{}.nc'.format(previous_hour_string)
    ND_hourly_path = os.path.join(netcdf_output_dir, ND_hourly_filename)
    telegram_hourly_filename = 'tensec_{}.nc'.format(previous_hour_string)
    telegram_hourly_path = os.path.join(netcdf_output_dir, telegram_hourly_filename)
    # First see if the file already exists
    if not os.path.exists(onesec_hourly_path):
        onesec_df_dump = onesec_df.loc[previous_hour_timestamp:last_hour_timestamp]
        spectrum_da_dump = spectrum_da.loc[previous_hour_timestamp:last_hour_timestamp]
        ND_da_dump = ND_da.loc[previous_hour_timestamp:last_hour_timestamp]
        telegram_df_dump = telegram_df.loc[previous_hour_timestamp:last_hour_timestamp]

        onesec_ds_dump = onesec_df_dump.to_xarray()
        onesec_ds_dump.to_netcdf(onesec_hourly_path)
        onesec_ds_dump.close()
        telegram_ds_dump = telegram_df_dump.to_xarray()
        telegram_ds_dump.to_netcdf(telegram_hourly_path)
        telegram_ds_dump.close()
        spectrum_ds_dump = spectrum_da_dump.to_dataset(name='spectrum')
        spectrum_ds_dump.to_netcdf(spectrum_hourly_path)
        spectrum_ds_dump.close()
        ND_ds_dump = ND_da_dump.to_dataset(name='spectrum')
        ND_ds_dump.to_netcdf(ND_hourly_path)
        ND_ds_dump.close()

# Create the figures
fig, ax = plt.subplots()
fig_vd, ax_vd = plt.subplots()
fig_dsd, ax_dsd = plt.subplots()
fig_t_td, ax_t_td = plt.subplots()
fig_wind, ax_windspd = plt.subplots()
ax_winddir = ax_windspd.twinx()
fig_pressure, ax_pressure = plt.subplots()

plottimes_onesec = [onesec_df.index.to_pydatetime()]
# Temperature and Dewpoint
Tmin = np.nanmin(onesec_df['Dewpoint'].values)
Tmax = np.nanmax(onesec_df['SlowTemp'].values)
fields_to_plot_onesec = [onesec_df['SlowTemp'].values, onesec_df['Dewpoint'].values]
temp_params = pm.temp_params.copy()
dewpoint_params = pm.dewpoint_params.copy()
temp_params['plotmin'] = Tmin - 5.0
dewpoint_params['plotmin'] = Tmin - 5.0
field_parameters_onesec = [temp_params, dewpoint_params]
ax_t_td = pm.plotmeteogram(
    ax_t_td,
    plottimes_onesec,
    fields_to_plot_onesec,
    field_parameters_onesec)
temp_dewp_ax_params['axeslimits'] = [[plottimes_onesec[0][0], plottimes_onesec[0][-1]],
                                     [Tmin - 5.0, Tmax + 5.0]]
ax_t_td, = pm.set_meteogram_axes([ax_t_td], [temp_dewp_ax_params])
# Wind speed and direction
ax_windspd = pm.plotmeteogram(
    ax_windspd, plottimes_onesec, [
        onesec_df['WS_ms'].values], [
        pm.windspeed_params])
ax_winddir = pm.plotmeteogram(
    ax_winddir, plottimes_onesec, [
        onesec_df['WindDir'].values], [
        pm.winddir_params])
windspd_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
winddir_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
ax_windspd, ax_winddir = pm.set_meteogram_axes(
    [ax_windspd, ax_winddir], [windspd_ax_params, winddir_ax_params])
# Pressure
pmin = np.nanmin(onesec_df['Pressure'].values)
pmax = np.nanmax(onesec_df['Pressure'].values)
pressure_ax_params['axeslimits'] = [[plottimes_onesec[0][0], plottimes_onesec[0][-1]],
                                    [pmin - 2.5, pmax + 2.5]]
fields_to_plot_press = [onesec_df['Pressure'].values]
field_parameters_press = [pm.pressure_params]
ax_pressure = pm.plotmeteogram(
    ax_pressure,
    plottimes_onesec,
    fields_to_plot_press,
    field_parameters_press)
ax_pressure, = pm.set_meteogram_axes([ax_pressure], [pressure_ax_params])

# DSD plots
plottimes_tmp = ND_da['time'].to_index().to_pydatetime()
# Prepend an additional at the beginning of the array so that pcolor sees this as the
# edges of the DSD intervals.
plottimes = np.insert(plottimes_tmp, 0, plottimes_tmp[0] - timedelta(seconds=10))
plottimes = [plottimes]
# Do the same at the end of the "min_diameter" and "min_fall_bins" arrays
# TODO: apparently this is a problem in the main DSD meteogram plotting script as well, so
# go back and fix that soon! NOTE: this is now done in PIPS.py and stored in "diameter_edges"
# and "fall_bins_edges"

ND_arr = ND_da.values.T
logND_arr = np.ma.log10(ND_arr)
fields_to_plot = [logND_arr]
field_parameters = [log_ND_params]
ax = pm.plotmeteogram(ax, plottimes, fields_to_plot, field_parameters,
                      yvals=[diameter_edges] * len(fields_to_plot))
ax, = pm.set_meteogram_axes([ax], [log_ND_ax_params])
ax_vd.set_title(
    r'$v_T$ vs. $D$ for time {0} to {1}'.format(
        (last_timestamp_60sec - timedelta(seconds=10)).strftime(tm.timefmt2),
        last_timestamp_tensec.strftime(tm.timefmt2)))
# spectrum_da.loc[last_timestamp_tensec]
spectrum = spectrum_da.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).sum(dim='time')
countsplot = np.ma.masked_where(spectrum.values <= 0, spectrum)
C = ax_vd.pcolormesh(diameter_edges, fall_bins_edges, countsplot, vmin=1, vmax=50, edgecolors='w')
divider = make_axes_locatable(ax_vd)
cax = divider.append_axes("right", size="5%")
cb = fig_vd.colorbar(C, cax=cax, orientation='vertical')
cax.set_ylabel('# of drops')
ax_vd.set_xlim(0.0, 10.0)
ax_vd.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
ax_vd.set_xlabel('diameter (mm)')
ax_vd.set_ylim(0.0, 10.0)
ax_vd.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax_vd.set_ylabel('fall speed (m/s)')

# Plot last 60-s DSD
# Bugfix 03/12/22: changed from sum to mean in following line!
ND_60s_da = ND_da.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).mean(dim='time')
ND_60s = ND_60s_da.values
ax_dsd.set_title('DSDs for time {0} to {1}'.format(
    (last_timestamp_60sec - timedelta(seconds=10)).strftime(tm.timefmt2),
    last_timestamp_tensec.strftime(tm.timefmt2)))
ax_dsd.bar(min_diameter, ND_60s * 1000.0, max_diameter - min_diameter, 10.**2., align='edge',
           log=True, color='tan', edgecolor='k')
ax_dsd.set_xlim(0.0, 9.0)
ax_dsd.set_ylim(10.**2., 10.**8.5)
ax_dsd.set_xlabel('D (mm)')
ax_dsd.set_ylabel(r'N(D) $(m^{-4})$')

# fig.canvas.draw()
# Try-except block to test if the mount is working and to reset if it isn't. This is a terrible
# hack that I need to clean up later.
try:
    fig.savefig(os.path.join(image_output_dir, 'logND_current.png'), dpi=300)
except:
    print("that didn't work")
    # Assume mount isn't working
    try:
        print("unmounting stormlab")
        os.system("/Users/dawson29/scripts/stormlabunmount")
    except:
        print("mounting stormlab")
        os.system("/Users/dawson29/scripts/stormlabmount")
    finally:
        print("mounting stormlab")
        os.system("/Users/dawson29/scripts/stormlabmount")
fig.savefig(os.path.join(image_output_dir, 'logND_current.png'), dpi=300)
fig_vd.savefig(os.path.join(image_output_dir, 'VD_current.png'), dpi=300)
fig_dsd.savefig(os.path.join(image_output_dir, 'DSD_current.png'), dpi=300)
fig_t_td.savefig(os.path.join(image_output_dir, 'T_Td_current.png'), dpi=300)
fig_wind.savefig(os.path.join(image_output_dir, 'wind_current.png'), dpi=300)
fig_pressure.savefig(os.path.join(image_output_dir, 'pressure.png'), dpi=300)
