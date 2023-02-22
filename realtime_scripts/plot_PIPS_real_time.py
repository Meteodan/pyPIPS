# This notebook is for testing the download of PIPS data using the web API.
# It uses urllib3 and BeautifulSoup4
import os
import io
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
import pyPIPS.utils as utils
import pyPIPS.plotmodule as pm


def get_files(file_list, starttime, endtime, ftype='onesec'):
    ''' this seems like a janky way to do it, but is actually 3x faster than
        making a loop of files.'''

    len_prefix = 8 + len(ftype)

    day_str = [(starttime + timedelta(days=i)).strftime("%Y%m%d")
               for i in range((endtime - starttime).days + 1)]

    # find all PIPS files with days between starttime and endtime
    file_list = [file_name for file_name in file_list if any(day in file_name for day in day_str)]
    # file_list = [f for subf in file_list for f in subf]  # flatten list in case of multiple days

    if file_list:
        # sort files by date, then find nearest indices for all the dates, and loop over that
        sorted_files = sorted(file_list,
                              key=lambda f: datetime.strptime(f[len_prefix:len_prefix + 14],
                                                              '%Y%m%d%H%M%S'))
        starttimes = [datetime.strptime(f[len_prefix:len_prefix + 14], '%Y%m%d%H%M%S')
                      for f in sorted_files]
        endtimes = [datetime.strptime(f[len_prefix + 15:len_prefix + 29], '%Y%m%d%H%M%S')
                    for f in sorted_files]
        _, idx1 = min((abs(val - starttime), idx) for (idx, val) in enumerate(starttimes))
        _, idx2 = min((abs(val - endtime), idx) for (idx, val) in enumerate(endtimes))
        file_list = sorted_files[idx1:idx2 + 1]

    return file_list


eff_sensor_area = parsivel_parameters['eff_sensor_area_mm2'] * 1.e-6

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


# Parse the command line options
parser = argparse.ArgumentParser(description="Plots real-time PIPS data")
parser.add_argument('PIPS_name', default=None,
                    help='Name of the PIPS to query and plot data for')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--image-output-dir', dest='image_output_dir',
                    default='../web/perils_realtime/images',
                    help='output directory for web images')
parser.add_argument('--netcdf-dir-url', dest='netcdf_dir_url',
                    help='URL for directory containing netCDF files')
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
if not os.path.exists(args.image_output_dir):
    os.makedirs(args.image_output_dir)

# Construct the URL for the PIPS
pips_ip = prc.PIPS_PERiLS_IPs[args.PIPS_name]
base_url = f"http://{pips_ip}/?command=TableDisplay&table="
url_onesec = f"{base_url}One_Hz"
url_tensec = f"{base_url}Ten_Hz"

http = urllib3.PoolManager()

meteogram_duration = 900  # Length of meteogram plot in seconds
plot_update_interval = 60  # number of seconds between plotting activity
plot_update_interval_ts = pd.Timedelta(seconds=plot_update_interval)

starttime_loop = time.time()

# Main loop
while True:
    endtime = datetime.utcnow()
    starttime = endtime - timedelta(seconds=meteogram_duration)

    # Grab list of onesec files for the given PIPS
    content = http.request('GET', f'http://{args.netcdf_dir_url}')
    soup = bs4.BeautifulSoup(content.data, "lxml")
    table = soup.find('table')
    rows = table.find_all('tr')
    file_list = []
    for row in rows:
        try:
            file_link = row.find_all('a')[0]
            file_name = file_link.get('href')
        except IndexError:
            continue
        file_list.append(file_name)

    file_list_onePIPS = [file_name for file_name in file_list if args.PIPS_name in file_name]
    file_list_onesec = [file_name for file_name in file_list_onePIPS if 'onesec' in file_name]
    file_list_onesec = [file_name for file_name in file_list_onesec if 'current' not in file_name]
    file_list_onesec = get_files(file_list_onesec, starttime, endtime)

    # Now, load the list of files into memory and read them with xarray

    file_urls = [f'http://{args.netcdf_dir_url}/{file}' for file in file_list_onesec]
    print("Reading files from website")
    file_reqs = [http.request('GET', url) for url in file_urls]
    print("Extracting data from files into memory")
    fdata = [io.BytesIO(file_req.data) for file_req in file_reqs]
    print("Opening with xarray")
    onesec_ds = xr.open_mfdataset(fdata, combine='nested', concat_dim='logger_datetime')

    # Do the same for the ten-sec data files
    file_list_ND = [file_name for file_name in file_list_onePIPS if 'ND' in file_name]
    file_list_ND = [file_name for file_name in file_list_ND if 'current' not in file_name]
    file_list_ND = get_files(file_list_ND, starttime, endtime, ftype='ND')

    file_urls = [f'http://{args.netcdf_dir_url}/{file}' for file in file_list_ND]
    print("Reading ND files from website")
    file_reqs = [http.request('GET', url) for url in file_urls]
    print("Extracting data from files into memory")
    fdata = [io.BytesIO(file_req.data) for file_req in file_reqs]
    print("Opening with xarray")
    ND_ds = xr.open_mfdataset(fdata, combine='nested', concat_dim='time')

    file_list_spectrum = [file_name for file_name in file_list_onePIPS if 'spectrum' in file_name]
    file_list_spectrum = [file_name for file_name in file_list_spectrum if 'current' not in file_name]
    file_list_spectrum = get_files(file_list_spectrum, starttime, endtime, ftype='spectrum')

    file_urls = [f'http://{args.netcdf_dir_url}/{file}' for file in file_list_spectrum]
    print("Reading spectrum files from website")
    file_reqs = [http.request('GET', url) for url in file_urls]
    print("Extracting data from files into memory")
    fdata = [io.BytesIO(file_req.data) for file_req in file_reqs]
    print("Opening with xarray")
    spectrum_ds = xr.open_mfdataset(fdata, combine='nested', concat_dim='time')

    last_timestamp_onesec = onesec_ds['logger_datetime'].to_index().to_pydatetime()[-1]
    last_timestamp_tensec = spectrum_ds['time'].to_index().to_pydatetime()[-1]
    last_timestamp_60sec = spectrum_ds['time'].to_index().to_pydatetime()[-6]

    # Create the figures
    fig, ax = plt.subplots()
    fig_vd, ax_vd = plt.subplots()
    fig_dsd, ax_dsd = plt.subplots()
    fig_t_td, ax_t_td = plt.subplots()
    fig_wind, ax_windspd = plt.subplots()
    ax_winddir = ax_windspd.twinx()
    fig_pressure, ax_pressure = plt.subplots()

    plottimes_onesec = onesec_ds['logger_datetime'].to_index().to_pydatetime()
    plottimes_onesec = [plottimes_onesec]
    # Temperature and Dewpoint
    Tmin = np.nanmin(onesec_ds['dewpoint'].values)
    Tmax = np.nanmax(onesec_ds['fasttemp'].values)
    fields_to_plot_onesec = [onesec_ds['fasttemp'].values, onesec_ds['dewpoint'].values]
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
            onesec_ds['windspd'].values], [
            pm.windspeed_params])
    ax_winddir = pm.plotmeteogram(
        ax_winddir, plottimes_onesec, [
            onesec_ds['winddirabs'].values], [
            pm.winddir_params])
    windspd_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
    winddir_ax_params['axeslimits'][0] = (plottimes_onesec[0][0], plottimes_onesec[0][-1])
    ax_windspd, ax_winddir = pm.set_meteogram_axes(
        [ax_windspd, ax_winddir], [windspd_ax_params, winddir_ax_params])

    # Pressure
    pmin = np.nanmin(onesec_ds['pressure'].values)
    pmax = np.nanmax(onesec_ds['pressure'].values)
    pressure_ax_params['axeslimits'] = [[plottimes_onesec[0][0], plottimes_onesec[0][-1]],
                                        [pmin - 2.5, pmax + 2.5]]
    fields_to_plot_press = [onesec_ds['pressure'].values]
    field_parameters_press = [pm.pressure_params]
    ax_pressure = pm.plotmeteogram(
        ax_pressure,
        plottimes_onesec,
        fields_to_plot_press,
        field_parameters_press)
    ax_pressure, = pm.set_meteogram_axes([ax_pressure], [pressure_ax_params])

    # DSD plots
    plottimes_tmp = ND_ds['time'].to_index().to_pydatetime()
    # Prepend an additional at the beginning of the array so that pcolor sees this as the
    # edges of the DSD intervals.
    plottimes = np.insert(plottimes_tmp, 0, plottimes_tmp[0] - timedelta(seconds=10))
    plottimes = [plottimes]
    # Do the same at the end of the "min_diameter" and "min_fall_bins" arrays
    # TODO: apparently this is a problem in the main DSD meteogram plotting script as well, so
    # go back and fix that soon! NOTE: this is now done in PIPS.py and stored in "diameter_edges"
    # and "fall_bins_edges"

    ND_arr = ND_ds['ND'].values.T
    print(ND_arr.shape)
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
    spectrum = spectrum_ds.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).sum(dim='time')
    spectrum = spectrum['spectrum']
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
    ND_60s_da = ND_ds.sel(time=slice(last_timestamp_60sec, last_timestamp_tensec)).mean(dim='time')
    ND_60s = ND_60s_da['ND'].values
    ax_dsd.set_title('DSDs for time {0} to {1}'.format(
        (last_timestamp_60sec - timedelta(seconds=10)).strftime(tm.timefmt2),
        last_timestamp_tensec.strftime(tm.timefmt2)))
    ax_dsd.bar(min_diameter, ND_60s * 1000.0, max_diameter - min_diameter, 10.**2., align='edge',
            log=True, color='tan', edgecolor='k')
    ax_dsd.set_xlim(0.0, 9.0)
    ax_dsd.set_ylim(10.**2., 10.**8.5)
    ax_dsd.set_xlabel('D (mm)')
    ax_dsd.set_ylabel(r'N(D) $(m^{-4})$')

    # Main loop
    # while True:
    #     # Create the figures
    #     # fig, ax = plt.subplots()
    #     # fig_vd, ax_vd = plt.subplots()
    #     # fig_dsd, ax_dsd = plt.subplots()
    #     # fig_t_td, ax_t_td = plt.subplots()
    #     # fig_wind, ax_windspd = plt.subplots()
    #     # ax_winddir = ax_windspd.twinx()
    #     # fig_pressure, ax_pressure = plt.subplots()

    #     curtime = datetime.utcnow()
    #     starttime = curtime - timedelta(seconds=meteogram_duration)

    #     # Grab list of onesec files for the given PIPS
    #     content = http.request('GET', f'http://{args.netcdf_dir_url}')

    #     	PIPS1A_onesec_20230222170657_20230222170806.nc


    fig.savefig(os.path.join(args.image_output_dir, 'logND_current.png'), dpi=300)
    fig_vd.savefig(os.path.join(args.image_output_dir, 'VD_current.png'), dpi=300)
    fig_dsd.savefig(os.path.join(args.image_output_dir, 'DSD_current.png'), dpi=300)
    fig_t_td.savefig(os.path.join(args.image_output_dir, 'T_Td_current.png'), dpi=300)
    fig_wind.savefig(os.path.join(args.image_output_dir, 'wind_current.png'), dpi=300)
    fig_pressure.savefig(os.path.join(args.image_output_dir, 'pressure.png'), dpi=300)

    time.sleep(plot_update_interval - ((time.time() - starttime_loop) % plot_update_interval))
