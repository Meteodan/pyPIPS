# This script grabs data from the TriPIPS in real-time and plots it, saving png files
# to a directory that is read by the real-time website (stormlab.eaps.purdue.edu/realtime)
import os
import bs4
import urllib3
import pandas as pd
import xarray as xr
import numpy as np
from pyPIPS.disdrometer_module import avg_diameter, fall_bins, eff_sensor_area, max_diameter,      \
                                      min_diameter
import time
import matplotlib.pyplot as plt
# from pyPIPS import thermolib as thermo
from datetime import datetime, timedelta
import matplotlib.dates as dates
import matplotlib.ticker as ticker
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
        'rain rate (mm/hr)': float(parsivel_tokens[1]),
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
    if(spectrum.size == 1024):
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
    min = np.int(time[3:5])
    sec = np.int(time[6:8])

    return datetime(year, month, day, hour, min, sec)


def calc_ND_da(spectrum_da, interval=10, use_measured_fs=True):
    """Computes the number concentration from the 32x32 spectrum"""
    # index = spectrum_da.coords['time'].values
    if not use_measured_fs:
        raise NotImplementedError('Not yet implemented: use measured fall speed for now!')
    else:
        ND_da = spectrum_da.groupby('time').apply(calc_ND)

    # ND_df = pd.DataFrame(ND_arr, columns=avg_diameter, index=index)
    return ND_da


# def calc_ND_list(spectrum_list, interval=10, use_measured_fs=True):
#     ND_list = [calc_ND(spectrum, interval=DSD_interval) for spectrum in spectrum_list]
#     ND_arr = np.array(ND_list)
#     return ND_arr


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
base_url = "http://10.163.29.26/?command=TableDisplay&table="
onesec_table = "One_Hz"
tensec_table = "Ten_Hz"
url_onesec = base_url + onesec_table
url_tensec = base_url + tensec_table

# Set up dictionaries to control plotting parameters

dateformat = '%H:%M'

# Temperature and dewpoint
temp_dewp_ax_params = {
    'majorxlocator': dates.MinuteLocator(interval=15),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(interval=5),
    'axeslimits': [None, (-10., 20.)],
    'axeslabels': ['Time (H:M) UTC', r'Temperature ($^{\circ}$C)']
}

# Number concentration
log_ND_params = {
    'type': 'pcolor',
    'vlimits': [-1.0, 3.0],
    'clabel': r'log[N ($m^{-3} mm^{-1}$)]'
}

log_ND_ax_params = {
    'majorxlocator': dates.MinuteLocator(interval=15),
    'majorxformatter': dates.DateFormatter(dateformat),
    'minorxlocator': dates.MinuteLocator(interval=5),
    'axeslimits': [None, [0.0, 9.0]],
    'majorylocator': ticker.MultipleLocator(base=1.0),
    'axeslabels': [None, 'D (mm)']
}

# # Testing real-time plotting loop
# interval_onesec = 1.
# interval_tensec = 10.

# plot_update_interval = 5 # seconds
# plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)
# keep_data_for = 3600 # seconds
# keep_data_for_ts = pd.Timedelta(seconds=keep_data_for)
# numrecords_onesec = int(keep_data_for / interval_onesec)
# numrecords_tensec = int(keep_data_for / interval_tensec)

# # Grab and process the initial period of data
# onesec_df = scrape_tripips_onesec_data(url_onesec, numrecords=numrecords_onesec)
# onesec_df['Dewpoint'] = thermo.calTdfromRH(onesec_df['Pressure'] * 100., onesec_df['SlowTemp'] + 273.15,
#                                      onesec_df['RH'] / 100.) - 273.15
# telegram_df, spectrum_da = scrape_tripips_tensec_data(url_tensec, numrecords=numrecords_tensec)
# ND_da = calc_ND_da(spectrum_da)

# # TODO: figure out an efficient way to dump data containing regular intervals (i.e. hourl) to netCDF files.
# # May eventually be able to replace current dumping of proprietary files to the flash drive on the CR6.

# # fig, ax = plt.subplots()
# # onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), color=['r', 'g'])

# fig, ax = plt.subplots()
# plt.ion()
# fig.show()
# fig.canvas.draw()

# numrecords_append=10
# numiter=10
# # onesec_new_df = scrape_tripips_onesec_data(url_onesec, numrecords=numrecords)
# # onesec_new_df['Dewpoint'] = thermo.calTdfromRH(onesec_new_df['Pressure'] * 100., onesec_new_df['SlowTemp'] + 273.15,
# #                                                onesec_new_df['RH'] / 100.) - 273.15
# # # Append new data onto onesec_df
# # onesec_df = onesec_df.append(onesec_new_df)
# # # Drop duplicate timestamps
# # onesec_df = onesec_df[~onesec_df.index.duplicated(keep='first')]
# # # Drop records older than specified time
# # last_timestamp = onesec_df.index[-1]
# # oldest_timestamp = last_timestamp-keep_data_for
# # onesec_df = onesec_df[oldest_timestamp:]
# # onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), colors=['r', 'g'])

# for i in range(numiter):
#     ax.clear()
#     onesec_new_df = scrape_tripips_onesec_data(url_onesec, numrecords=numrecords_append)
#     onesec_new_df['Dewpoint'] = thermo.calTdfromRH(onesec_new_df['Pressure'] * 100., onesec_new_df['SlowTemp'] + 273.15,
#                                      onesec_new_df['RH'] / 100.) - 273.15
#     # Append new data onto onesec_df
#     onesec_df = onesec_df.append(onesec_new_df)
#     # Drop duplicate timestamps
#     onesec_df = onesec_df[~onesec_df.index.duplicated(keep='first')]
#     # Drop records older than desired interval
#     last_timestamp = onesec_df.index[-1]
#     oldest_timestamp = last_timestamp-keep_data_for_ts
#     onesec_df = onesec_df[oldest_timestamp:]
#     plottimes = [onesec_df.index.to_pydatetime()]
#     fields_to_plot = [onesec_df['SlowTemp'].values, onesec_df['Dewpoint'].values]
#     field_parameters = [pm.temp_params, pm.dewpoint_params]
#     ax = pm.plotmeteogram(ax, plottimes, fields_to_plot, field_parameters)
#     ax, = pm.set_meteogram_axes([ax], [temp_dewp_ax_params])
#     #onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), color=['r', 'g'])
#     fig.canvas.draw()
#     # plt.pause(0.01)
#     # Sleep for the desired interval. This may not be perfectly accurate
#     time.sleep(plot_update_interval)

# Testing real-time plotting loop

interval_onesec = 1.
interval_tensec = 10.

plot_update_interval = 20  # seconds
plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)
keep_data_for = 3600  # seconds
keep_data_for_ts = pd.Timedelta(seconds=keep_data_for)
numrecords_onesec = int(keep_data_for / interval_onesec)
numrecords_tensec = int(keep_data_for / interval_tensec)

# Grab and process the initial period of data
telegram_df, spectrum_da = scrape_tripips_tensec_data(url_tensec, numrecords=numrecords_tensec)
ND_da = calc_ND_da(spectrum_da)

# TODO: figure out an efficient way to dump data containing regular intervals (i.e. hourl) to
# netCDF files. May eventually be able to replace current dumping of proprietary files to the flash
# drive on the CR6.

# fig, ax = plt.subplots()
# onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), color=['r', 'g'])

fig, ax = plt.subplots()
plt.ion()
fig.show()
fig.canvas.draw()
# display.display(fig)
# display.clear_output(wait=True)
numrecords_append = 10
numiter = 5
# onesec_new_df = scrape_tripips_onesec_data(url_onesec, numrecords=numrecords)
# onesec_new_df['Dewpoint'] = thermo.calTdfromRH(onesec_new_df['Pressure'] * 100., onesec_new_df['SlowTemp'] + 273.15,
#                                                onesec_new_df['RH'] / 100.) - 273.15
# # Append new data onto onesec_df
# onesec_df = onesec_df.append(onesec_new_df)
# # Drop duplicate timestamps
# onesec_df = onesec_df[~onesec_df.index.duplicated(keep='first')]
# # Drop records older than specified time
# last_timestamp = onesec_df.index[-1]
# oldest_timestamp = last_timestamp-keep_data_for
# onesec_df = onesec_df[oldest_timestamp:]
# onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), colors=['r', 'g'])

# for i in range(numiter):
while True:
    ax.clear()
    telegram_new_df, spectrum_new_da = scrape_tripips_tensec_data(url_tensec, numrecords=numrecords_append)
    ND_new_da = calc_ND_da(spectrum_new_da)
    # Append new data onto the data array
    ND_da = xr.concat([ND_da, ND_new_da], dim='time')
    # Drop duplicate timestamps
    ND_da = ND_da.groupby('time').first()
    # onesec_df = onesec_df[~onesec_df.index.duplicated(keep='first')]
    # Drop records older than desired interval
    last_timestamp = ND_da['time'].to_index()[-1]
    # print(last_timestamp)
    oldest_timestamp = last_timestamp-keep_data_for_ts
    ND_da = ND_da.loc[oldest_timestamp:]
    plottimes = [ND_da['time'].to_index().to_pydatetime()]
    ND_arr = ND_da.values.T
    logND_arr = np.ma.log10(ND_arr)
    fields_to_plot = [logND_arr]
    field_parameters = [log_ND_params]
    ax = pm.plotmeteogram(ax, plottimes, fields_to_plot, field_parameters,
                          yvals=[min_diameter] * len(fields_to_plot))
    ax, = pm.set_meteogram_axes([ax], [log_ND_ax_params])
    # onesec_df.plot(ax=ax, y=['SlowTemp', 'Dewpoint'], ylim=(-5.,5.), color=['r', 'g'])
    # display.display(fig)
    # display.clear_output(wait=True)
    fig.canvas.draw()
    fig.savefig(os.path.join(image_output_dir, 'logND_current.png'), dpi=300)
    # plt.pause(0.01)
    # Sleep for the desired interval. This may not be perfectly accurate
    time.sleep(plot_update_interval)
