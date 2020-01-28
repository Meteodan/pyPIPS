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

# Extract needed lists and variables from the radar_dict configuration dictionary
load_radar_at_PIPS = config.radar_config_dict.get('load_radar_at_PIPS', False)
save_radar_at_PIPS = config.radar_config_dict.get('save_radar_at_PIPS', False)
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
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

# ------
# Grab radar data for comparison if desired

if comp_radar:
    # sb = radar.readsweeps2PIPS(fieldnames, pc, ib)
    if not load_radar_at_PIPS:
        radar_dict = radar.read_sweeps(radar_name, radar_dir, radar_start_timestamp,
                                       radar_end_timestamp, field_names=field_names, el_req=el_req)
        rad_locs = []
        for geo_loc in geo_locs:
            rad_loc = radar.get_PIPS_loc_relative_to_radar(geo_loc, radar_dict['radarsweeplist'][0])
            rad_locs.append(rad_loc)
        radar_fields_at_PIPS_da = radar.interp_sweeps_to_PIPS(radar_name,
                                                              radar_dict['radarsweeplist'],
                                                              PIPS_names, rad_locs)
        if save_radar_at_PIPS:
            radar_at_PIPS_file_name = '{}_{}_{}_fields_at_PIPS.nc'.format(radar_name,
                                                                          radar_start_timestamp,
                                                                          radar_end_timestamp)
            radar_at_PIPS_path = os.path.join(PIPS_dir, radar_at_PIPS_file_name)
            radar.dump_radar_fields_at_PIPS_nc(radar_at_PIPS_path, radar_fields_at_PIPS_da)
    else:
        radar_at_PIPS_file_name = '{}_{}_{}_fields_at_PIPS.nc'.format(radar_name,
                                                                      radar_start_timestamp,
                                                                      radar_end_timestamp)
        radar_at_PIPS_path = os.path.join(PIPS_dir, radar_at_PIPS_file_name)
        radar_fields_at_PIPS_da = xr.open_dataarray(radar_at_PIPS_path)

# Outer disdrometer (and deployment) loop
for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype, conv_df, \
    parsivel_df, vd_matrix_da in zip(range(0, len(PIPS_filenames)), PIPS_filenames, PIPS_names,
                                     start_times, end_times, geo_locs, PIPS_types, conv_df_list,
                                     parsivel_df_list, vd_matrix_da_list):

    tripips = (ptype == 'TriPIPS')

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Do some QC on the V-D matrix
    strongwindQC = pc.PIPS_qc_dict['strongwindQC']
    splashingQC = pc.PIPS_qc_dict['splashingQC']
    marginQC = pc.PIPS_qc_dict['marginQC']
    rainfallQC = pc.PIPS_qc_dict['rainfallQC']
    rainonlyQC = pc.PIPS_qc_dict['rainonlyQC']
    hailonlyQC = pc.PIPS_qc_dict['hailonlyQC']
    graupelonlyQC = pc.PIPS_qc_dict['graupelonlyQC']

    if pc.PIPS_qc_dict['basicQC']:
        strongwindQC = True
        splashingQC = True
        marginQC = True

    if strongwindQC:
        vd_matrix_da = pqc.strongwindQC(vd_matrix_da)
    if splashingQC:
        vd_matrix_da = pqc.splashingQC(vd_matrix_da)
    if marginQC:
        vd_matrix_da = pqc.marginQC(vd_matrix_da)
    if rainfallQC:
        fallspeedmask = pqc.get_fallspeed_mask(avg_diameter, avg_fall_bins)
        vd_matrix_da = pqc.rainfallspeedQC(vd_matrix_da, fallspeedmask)
    if rainonlyQC:
        vd_matrix_da = pqc.rainonlyQC(vd_matrix_da)
    if hailonlyQC:
        vd_matrix_da = pqc.hailonlyQC(vd_matrix_da)
    # if graupelonlyQC:
    #     vd_matrix_da = pqc.graupelonlyQC(vd_matrix_da)

    # Compute the number density ND from the V-D matrix
    # Use measured fallspeed by default
    empirical_fallspeed = pips.calc_empirical_fallspeed(avg_diameter)
    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    logND = np.log10(ND)

    # Get times for PIPS meteogram plotting
    PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])

    # If we are comparing with radar, get the radar timestamps as well
    if comp_radar:
        # plotx_rad = dates.date2num(sb.radtimes)
        plotx_rad = pd.to_datetime(radar_fields_at_PIPS_da['time'].values).to_pydatetime()
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
    except (ValueError, TypeError):
        start_datetime = PSD_edgetimes[0]
    try:
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
    except (ValueError, TypeError):
        end_datetime = PSD_edgetimes[-1]
    timelimits = [dates.date2num(start_datetime), dates.date2num(end_datetime)]
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
    minorlocator = dates.MinuteLocator(byminute=range(0, 60, 5))
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
    # If we are comparing with radar, grab associated radar variables for plotting
    if comp_radar:
        # print(radar_fields_at_PIPS_da)
        # At least add the reflectivity
        dBZ_D_plt = radar_fields_at_PIPS_da.sel(fields='REF', PIPS=PIPS_name)
        # indexrad = sb.outfieldnames.index('dBZ')
        # dBZ_D_plt = sb.fields_D_tarr[:, index, indexrad]
        radvars = {'radmidtimes': plotx_rad, 'REF': dBZ_D_plt}
        # Add other polarimetric fields
        if calc_dualpol:
            for radvarname in ['ZDR', 'KDP', 'RHO']:
                if radvarname in radar_fields_at_PIPS_da.fields:
                    # indexrad = sb.outfieldnames.index(radvarname)
                    # dualpol_rad_var = sb.fields_D_tarr[:, index, indexrad]
                    dualpol_rad_var = radar_fields_at_PIPS_da.sel(fields=radvarname, PIPS=PIPS_name)
                    if dualpol_rad_var.size:
                        radvars[radvarname] = dualpol_rad_var
            if clean_radar:
                # remove non-precipitation echoes from radar data
                gc_mask = np.where((radvars['RHO'] < 0.90), True, False)
                for radvarname in ['ZDR', 'REF', 'RHO']:
                    radvars[radvarname] = np.ma.masked_array(radvars[radvarname],
                                                             mask=gc_mask)

    # Make the plot
    PIPS_plot_name = '{}_{}_{}_{}'.format(PIPS_name, start_time_string, end_time_string, DSDtype)
    pm.plotDSDmeteograms(PIPS_plot_name, meteogram_image_dir, axparams, disvars, radvars.copy())
