"""Module for real-time monitoring functionality"""
from __future__ import annotations

import csv
import locale
import os
import shutil
from datetime import datetime

import bs4
import numpy as np
import pandas as pd
import urllib3
import xarray as xr

from pyPIPS import pips_io as pipsio
from pyPIPS.parsivel_params import parsivel_parameters
from pyPIPS.PIPS import (
    avg_diameter,
    avg_fall_bins,
    max_diameter,
    min_diameter,
)

eff_sensor_area = parsivel_parameters['eff_sensor_area_mm2'] * 1.e-6
# Function definitions


def get_data_table(http, url, numrecords=3600):
    """Grabs the data table from the PIPS http server given the requested http pool manager,
       requested URL, and desired number of records. Uses urlib3 and beautifulsoup4"""
    content = http.request('GET', url + f'&records={numrecords:d}', retries=False,
                           timeout=10.0)
    soup = bs4.BeautifulSoup(content.data, "lxml")
    table = soup.find('table')
    rows = table.find_all('tr')
    headers = rows.pop(0).find_all('th')
    return headers, rows


def scrape_onesec_data(url, http=None, numrecords=3600, tripips=False):
    """Grabs one-second records from the PIPS http server. Uses urlib3 and beautifulsoup4"""
    last_good_gps_dict = {}
    if http is None:
        http = urllib3.PoolManager()
    headers, rows = get_data_table(http, url, numrecords=numrecords)
    # headers.pop(0)
    headers = [header.string for header in headers]
    field_indices = pipsio.get_field_indices(headers, tripips=tripips)
    onesec_dict_list = []
    for row in rows:
        tokens = row.find_all('td')
        tokens = [token.string.strip(' "') for token in tokens]

        token_dict = pipsio.parse_PIPS_record(tokens, field_indices, tripips=tripips,
                                              include_parsivel_string=False)

        # Do some additional checks on the GPS data. Sometimes it gets garbled, so in case that it
        # does, we keep a running record of the last known "good" GPS info and insert it into
        # the bad fields when needed.

        if token_dict['GPS_status'] == "A":
            GPS_keys = [key for key in token_dict if "GPS" in key]
            for key in GPS_keys:
                last_good_gps_dict[key] = token_dict[key]
        elif token_dict['GPS_status'] != 'V' and last_good_gps_dict:
            for key, last_good_gps in last_good_gps_dict.items():
                token_dict[key] = last_good_gps

        # timestamp = token_dict['logger_datetime']
        # token_dict['logger_datetime'] = pd.to_datetime(timestamp, format='%Y-%m-%d %H:%M:%S')
        onesec_dict_list.append(token_dict)
        # timestamp = token_dict['TIMESTAMP']

        # timestamp = tokens.pop(0)
        # tokens = [np.float64(token) if token not in ('' or 'NAN') else np.nan for token in tokens]
        # data.append(tokens)
        # timestamps.append(timestamp)
    # index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S')
    # df = pd.DataFrame(data, columns=headers, index=index)
    onesec_df = pd.DataFrame(onesec_dict_list)
    onesec_df = onesec_df.set_index('logger_datetime')
    return onesec_df  # noqa: RET504


def scrape_tensec_data(url, http=None, numrecords=360):
    """Grabs ten-second records from the PIPS http server. Uses urlib3 and beautifulsoup4"""
    if http is None:
        http = urllib3.PoolManager()
    headers, rows = get_data_table(http, url, numrecords=numrecords)
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
            print("Bad Parsivel string detected. Ignoring this record!")  # noqa: T201
    try:
        index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S')
    except ValueError:  # Sometimes the timestring has microseconds
        index = pd.to_datetime(timestamps, format='%Y-%m-%d %H:%M:%S.%f')
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
    return parsivel_telegram_dict  # noqa: RET504


def read_parsivel_spectrum(parsivel_string):
    """Given a raw string of Parsivel data, extract the DSD spectrum from it."""
    parsivel_tokens = parsivel_string.strip(' "').split(';')
    # parsivel_tokens = parsivel_tokens.strip('"')
    spectrum = [float(x) if ('' or '\n' or '\r') not in x else 0 for x in parsivel_tokens[11:]]  # noqa: SIM222
    try:
        spectrum = [float(x) if ('' or '\n' or '\r') not in x else 0 for x in parsivel_tokens[11:]]  # noqa: SIM222
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
    year = int(date[:4])
    month = int(date[5:7])
    day = int(date[8:])
    hour = int(time[:2])
    minute = int(time[3:5])
    sec = int(time[6:8])

    return datetime(year, month, day, hour, minute, sec)


def calc_ND_da(spectrum_da, interval=10, use_measured_fs=True):
    """Computes the number concentration from the 32x32 spectrum"""
    # index = spectrum_da.coords['time'].values
    if not use_measured_fs:
        msg = 'Not yet implemented: use measured fall speed for now!'
        raise NotImplementedError(msg)
    ND_da = spectrum_da.groupby('time').apply(calc_ND, interval=interval,
                                              use_measured_fs=use_measured_fs)

    # ND_df = pd.DataFrame(ND_arr, columns=avg_diameter, index=index)
    return ND_da  # noqa: RET504


def calc_ND_list(spectrum_list, interval=10, use_measured_fs=True):  # noqa: ARG001
    ND_list = [calc_ND(spectrum, interval=interval) for spectrum in spectrum_list]
    ND_arr = np.array(ND_list)
    return ND_arr  # noqa: RET504


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


def dump_real_time_csv_for_SASSI(state_dict, csv_output_dir, PIPS_name, timestamp):
    """Dumps a one-line CSV file with the current state for a given PIPS for ingest by SASSI"""
    timestamp_str = timestamp.strftime('%Y%m%d%H%M%S')
    filename = f'{PIPS_name}_{timestamp_str}.txt'
    filepath = os.path.join(csv_output_dir, filename)
    field_names = ['ID', 'lon', 'lat', 'elev', 'FastTemp', 'SlowTemp', 'RH', 'p', 'dir', 'spd']

    with open(filepath, 'w', newline='', encoding=locale.getpreferredencoding(False)) as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=field_names)

        # writer.writeheader()
        writer.writerow(state_dict)


def dump_real_time_df_netcdf(df, netcdf_output_dir, ds_name, PIPS_name, copy_to_ts_file=True):
    """Dumps a netcdf file from real-time PIPS data to disk. Works on pandas dataframes
       and converts to xarray prior to saving to disk. Also optionally copies the
       dumped file to one with unique timestamps"""
    first_timestamp = df.index[0].strftime('%Y%m%d%H%M%S')
    last_timestamp = df.index[-1].strftime('%Y%m%d%H%M%S')
    netcdf_filename = f'{PIPS_name}_{ds_name}_current.nc'
    netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
    ds = df.to_xarray()
    ds.to_netcdf(netcdf_path)
    ds.close()
    if copy_to_ts_file:
        netcdf_ts_filename = \
            f'{PIPS_name}_{ds_name}_{first_timestamp}_{last_timestamp}.nc'
        netcdf_ts_path = os.path.join(netcdf_output_dir, netcdf_ts_filename)
        shutil.copy(netcdf_path, netcdf_ts_path)


def dump_real_time_da_netcdf(da, netcdf_output_dir, ds_name, PIPS_name, copy_to_ts_file=True):
    """Dumps a netcdf file from real-time PIPS data to disk. Works on pandas dataframes
       and converts to xarray prior to saving to disk. Also optionally copies the
       dumped file to one with unique timestamps"""
    # print(da)
    first_timestamp = pd.to_datetime(da.coords['time'].isel(time=0).item()).strftime('%Y%m%d%H%M%S')
    last_timestamp = pd.to_datetime(da.coords['time'].isel(time=-1).item()).strftime('%Y%m%d%H%M%S')
    netcdf_filename = f'{PIPS_name}_{ds_name}_current.nc'
    netcdf_path = os.path.join(netcdf_output_dir, netcdf_filename)
    ds = da.to_dataset(name=ds_name)
    ds.to_netcdf(netcdf_path)
    ds.close()
    if copy_to_ts_file:
        netcdf_ts_filename = \
            f'{PIPS_name}_{ds_name}_{first_timestamp}_{last_timestamp}.nc'
        netcdf_ts_path = os.path.join(netcdf_output_dir, netcdf_ts_filename)
        shutil.copy(netcdf_path, netcdf_ts_path)
