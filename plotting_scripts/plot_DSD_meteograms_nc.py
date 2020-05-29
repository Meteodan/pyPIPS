# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import argparse
from datetime import datetime, timedelta
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
parser.add_argument('--plot-raw', action='store_true', dest='plot_raw',
                    help='plot raw ND')
parser.add_argument('--use-filtered-fields', dest='use_filtered_fields', default=False,
                    action='store_true',
                    help='Whether to use previously filtered dBZ and ZDR fields for the retrieval')
parser.add_argument('--plot-dir', metavar='<path/to/plot/directory/>', dest='plot_dir',
                    default=None,
                    help='directory to store plots (overrides that in the config file')
parser.add_argument('--filter_RR', dest='filter_RR', type=float, default=None,
                    help='filter rainrate < # (mm)')
parser.add_argument('--filter_counts', dest='filter_counts', type=int, default=None,
                    help='filter particle counts < #')
parser.add_argument('--use-parsivel-params', dest='use_parsivel_params', action='store_true',
                    default=False, help='Use parsivel RR and counts instead of computed')

args = parser.parse_args()

plot_raw = args.plot_raw

if plot_raw:
    tag = ''
else:
    tag = '_qc'

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
if args.plot_dir:
    plot_dir = args.plot_dir
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
    image_dir = os.path.join(meteogram_image_dir, deployment_name)
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    fname_tag = tag
    if comp_radar:
        radar_fields_at_PIPS_da = parsivel_combined_ds['{}_at_PIPS'.format(radar_name)]
        fname_tag = tag + '_{}'.format(radar_name)

    start_time = start_times[index]
    end_time = end_times[index]

    # Filter by RR and pcount if desired
    if args.filter_RR:
        if args.use_parsivel_params:
            rainrate_key = 'precipintensity'
        else:
            rainrate_key = 'rainrate_derived{}'.format(tag)
        parsivel_combined_ds = parsivel_combined_ds.where(
            parsivel_combined_ds[rainrate_key] >= args.filter_RR)
    if args.filter_counts:
        if args.use_parsivel_params:
            counts_key = 'pcount'
        else:
            counts_key = 'pcount_derived{}'.format(tag)
        parsivel_combined_ds = parsivel_combined_ds.where(
            parsivel_combined_ds[counts_key] >= args.filter_counts)

    ND = parsivel_combined_ds['ND{}'.format(tag)]
    logND = np.log10(ND)

    # Get times for PIPS meteogram plotting
    PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    # PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_edgetimes = PSD_datetimes_dict['PSD_datetimes_edges']
    # PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])
    PSD_centertimes = PSD_datetimes_dict['PSD_datetimes_centers']

    # Pack plotting variables into dictionary
    disvars = {
        'min_diameter': min_diameter,
        'PSDstarttimes': PSD_edgetimes,
        'PSDmidtimes': PSD_centertimes,
        'logND': logND.T
    }

    # Compute additional derived parameters
    disvars['D_0'] = dsd.calc_D0_bin(ND) * 1000.  # Get to mm
    if calc_dualpol:
        # Calculate polarimetric variables using the T-matrix technique
        # Note, may try to use pyDSD for this purpose.
        scattfile = os.path.join(scatt_dir, 'SCTT_RAIN_fw100.dat')
        # Observed DSD
        dualpol_dis = dp.calpolrain(wavelength, scattfile, ND, bin_width)
        for varname in ['REF', 'ZDR', 'KDP', 'RHO']:
            var = dualpol_dis.get(varname, np.empty((0)))
            if var.size:
                # If desired, perform centered running averages
                if pc.PIPS_plotting_dict['avgwindow'] and False:
                    window = int(pc.PIPS_plotting_dict['avgwindow'] / DSD_interval)
                    var = pd.Series(var).rolling(
                        window=window, center=True, win_type='triang',
                        min_periods=1).mean().values
                disvars[varname] = var
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
    try:
        diamlimits = pc.PIPS_plotting_dict['DSD_D_range']
        diamytick = pc.PIPS_plotting_dict['DSD_D_ytick']
    except KeyError:
        diamlimits = [0.0, 9.0]
        diamytick = 1.0

    DSDtype = 'observed'
    locator = dates.MinuteLocator(byminute=[0, 15, 30, 45])
    locator.MAXTICKS = 1500
    minorlocator = dates.MinuteLocator(byminute=range(0, 60, 5))
    minorlocator.MAXTICKS = 1500
    dateformat = '%H:%M'
    formatter = dates.DateFormatter(dateformat)

    axparams = {
        'majorxlocator': locator,
        'majorxformatter': formatter,
        'minorxlocator': minorlocator,
        'axeslimits': [timelimits, diamlimits],
        'majorylocator': ticker.MultipleLocator(base=diamytick),
        'axeslabels': [None, 'D (mm)']
    }
    radvars = {}
    # If we are comparing with radar, grab associated radar variables for plotting
    if comp_radar:
        # print(radar_fields_at_PIPS_da)
        # At least add the reflectivity
        # First get list of radar fields in DataArray
        dim_name = 'fields_{}'.format(radar_name)
        radar_fields = radar_fields_at_PIPS_da.coords[dim_name].values
        # Then find which one corresponds to reflectivity
        # ref_name = next((fname for fname in radar_fields if fname in radar.REF_aliases))
        ref_name = radar.find_radar_field_name(radar_fields, radar.REF_aliases)
        if args.use_filtered_fields:
            ref_name = ref_name + '_filtered'
        dBZ_D_plt = radar_fields_at_PIPS_da.loc[{dim_name: ref_name}]
        # indexrad = sb.outfieldnames.index('dBZ')
        # dBZ_D_plt = sb.fields_D_tarr[:, index, indexrad]
        radvars = {'radmidtimes': PSD_centertimes, 'REF': dBZ_D_plt}
        # Add other polarimetric fields
        if calc_dualpol:
            for radvar_name, radvar_aliases in zip(['ZDR', 'RHO', 'KDP'],
                                                   [radar.ZDR_aliases, radar.RHV_aliases,
                                                    radar.KDP_aliases]):
                radvar_name_in_file = radar.find_radar_field_name(radar_fields, radvar_aliases)
                if radvar_name == 'ZDR' and args.use_filtered_fields:
                    radvar_name_in_file = radvar_name_in_file + '_filtered'
                if radvar_name_in_file:
                    dualpol_rad_var = radar_fields_at_PIPS_da.loc[{dim_name: radvar_name_in_file}]
                    if dualpol_rad_var.size:
                        radvars[radvar_name] = dualpol_rad_var
            # TODO: this code below is obsolescent, as filtering is done elsewhere now.
            if clean_radar:
                # remove non-precipitation echoes from radar data
                gc_mask = np.where((radvars['RHO'] < 0.90), True, False)
                for radvarname in ['ZDR', 'REF', 'RHO']:
                    radvars[radvarname] = np.ma.masked_array(radvars[radvarname],
                                                             mask=gc_mask)
        if plot_retrieval:
            # Plot D_0 as overlay on other plots for now
            radvars['D_0_rad'] = radar_fields_at_PIPS_da.loc[{dim_name: 'D0'}]


    # Make the plot
    PIPS_plot_name = '{}_{}_{}_{}_{}{}'.format(PIPS_name, deployment_name, start_time_string,
                                               end_time_string, DSDtype, fname_tag)
    pm.plotDSDmeteograms(PIPS_plot_name, image_dir, axparams, disvars, radvars.copy())
