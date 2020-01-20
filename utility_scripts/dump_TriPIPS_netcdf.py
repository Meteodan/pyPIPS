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
from pyPIPS.disdrometer_module import avg_diameter, fall_bins, eff_sensor_area, max_diameter, \
    min_diameter, min_fall_bins
from pyPIPS import thermolib as thermo
from pyPIPS import timemodule as tm
import pyPIPS.plotmodule as pm


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
        tokens = [np.float(token) if token not in ('' or 'NAN') else np.nan for token in tokens]
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
        telegram_dict = read_parsivel_telegram(parsivel_string)
        telegrams.append(telegram_dict)
        spectrum = read_parsivel_spectrum(parsivel_string)
        spectrum_list.append(spectrum)
        timestamps.append(timestamp)

    index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S')
    df_telegram = pd.DataFrame(telegrams, index=index)
    # print(np.array(spectrum_list).shape)
    da_spectrum = xr.DataArray(spectrum_list, coords=[index, fall_bins,
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
    _, vspectrum = np.meshgrid(avg_diameter, fall_bins)
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
netcdf_output_dir = '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/netcdf_web/'
base_url = "http://10.163.29.26/?command=TableDisplay&table="
onesec_table = "One_Hz"
tensec_table = "Ten_Hz"
url_onesec = base_url + onesec_table
url_tensec = base_url + tensec_table

# Set up dictionaries to control plotting parameters

dateformat = '%H:%M'

# Temperature and dewpoint
temp_dewp_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'axeslimits': [None, (-5., 35.)],
    'axeslabels': ['Time (H:M) UTC', r'Temperature ($^{\circ}$C)']
}

# Wind speed and direction
windspd_ax_params = {
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(byminute=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55],
                                         interval=1),
    'axeslimits': [None, [0.0, 25.0]],
    'axeslabels': ['Time (H:M) UTC', r'wind speed (m s$^{-1}$)']
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
    'axeslabels': ['Time (H:M) UTC', r'Pressure (hPa)']
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
    'axeslabels': [None, 'D (mm)']
}

interval_onesec = 1.
interval_tensec = 10.

plot_update_interval = 60  # seconds
plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)
# seconds; save slightly more than an hour so that we can be sure to have
# enough overlap to save hourly nc files that start and end at the top of
# the hour
keep_data_for = 5400
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

# # Dump onesec dataframe to netCDF file (via xarray)
# netcdf_filename = 'onesec_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
# onesec_df.to_xarray().to_netcdf(netcdf_path)
# # Dump raw spectrum to netcdf file
# netcdf_filename = 'spectrum_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
# spectrum_da.to_dataset(name='spectrum').to_netcdf(netcdf_path)
# # Dump ND_da to netcdf file
# netcdf_filename = 'ND_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
# ND_da.to_dataset(name='ND').to_netcdf(netcdf_path)
# # Dump telegram data to netcdf file
# netcdf_filename = 'tensec_telegram_current_hour.nc'
# netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
# print(telegram_df.to_xarray())
# telegram_df.to_xarray().to_netcdf(netcdf_path)

# Save top of the hour netCDF file
last_timestamp_onesec = onesec_df.index[-1]
last_timestamp_tensec = telegram_df.index[-1]
oldest_timestamp_onesec = last_timestamp_onesec-keep_data_for_ts
last_hour_timestamp = last_timestamp_onesec.replace(microsecond=0, second=0, minute=0)
previous_hour_timestamp = last_hour_timestamp - timedelta(hours=1)
diff = (last_timestamp_onesec - last_hour_timestamp).total_seconds()

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

    onesec_df_dump.to_xarray().to_netcdf(onesec_hourly_path)
    telegram_df_dump.to_xarray().to_netcdf(telegram_hourly_path)
    spectrum_da_dump.to_dataset(name='spectrum').to_netcdf(spectrum_hourly_path)
    ND_da_dump.to_dataset(name='spectrum').to_netcdf(ND_hourly_path)
