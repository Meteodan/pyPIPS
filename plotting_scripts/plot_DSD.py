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
import matplotlib.pyplot as plt
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
parser = argparse.ArgumentParser(description="Plots DSD histograms from PIPS data")
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

# Create the directory for the meteogram plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Read in the PIPS data for the deployment
conv_df_list = []
conv_resampled_df_list = []
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
    print("Lat/Lon/alt of {}: {}".format(PIPS_name, str(geo_loc)))

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

    # Resample the parsivel data to a longer interval if desired
    if requested_interval > 10.:
        DSD_interval = pips.check_requested_resampling_interval(requested_interval, 10.)
        vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
        parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
    else:
        DSD_interval = 10.

    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    logND = np.log10(ND)

    # Read in the previously computed DSD fits
    # TODO: Make sure that the DSD interval used to create these fits
    # matches that of the requested interval!
    fits_filename = 'DSD_fits_{}_{:d}s_{}_{}.nc'.format(PIPS_name, int(DSD_interval), start_time,
                                                        end_time)
    fits_path = os.path.join(PIPS_dir, fits_filename)
    DSD_fits_ds = xr.open_dataset(fits_path)
    # TODO: allow for choosing of which fits to read and to plot
    DSD_MM246 = DSD_fits_ds['DSD_MM246']
    DSD_TMM246 = DSD_fits_ds['DSD_TMM246']
    ND_MM246 = dsd.calc_binned_DSD_from_params(DSD_MM246.loc['N0'], DSD_MM246.loc['lamda'],
                                               DSD_MM246.loc['alpha'], ND['diameter'])
    ND_TMM246 = dsd.calc_binned_DSD_from_params(DSD_TMM246.loc['N0'], DSD_TMM246.loc['lamda'],
                                                DSD_TMM246.loc['alpha'], ND['diameter'])

    # Get times for plotting
    PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])

    # Resample conventional data to the parsivel times
    sec_offset = PSD_datetimes[0].second
    conv_resampled_df = pips.resample_conv(ptype, DSD_interval, sec_offset, conv_df)
    conv_resampled_df_index = conv_resampled_df.index.intersection(parsivel_df.index)
    conv_resampled_df = conv_resampled_df.loc[conv_resampled_df_index]

    DSD_image_dir = os.path.join(plot_dir, 'DSDs/{}'.format(PIPS_name))
    if not os.path.exists(DSD_image_dir):
        os.makedirs(DSD_image_dir)

    axdict = {
        'xbin_left': min_diameter,
        'xbin_mid': avg_diameter,
        'xbin_right': max_diameter,
        'xlim': (0.0, 9.0),
        'ylim': (10.**2., 10.**8.5),
        'interval': int(DSD_interval),
        'PIPS_name': PIPS_name
    }

    for t, time in enumerate(vd_matrix_da['time'].to_index()):
        if parsivel_df['pcount'].loc[time] > 0 or True:
            print("Plotting for {} and time {}".format(PIPS_name, time.strftime(tm.timefmt3)))
            axdict['time'] = time
            PSDdict = {'ND': ND.loc[time]}
            PSDfitdict = {
                # 'Exponential_24': (ND_MM24.loc[time_to_plot], 'Exp fit (MM24)'),
                # 'Exponential_36': (ND_MM36.loc[time_to_plot], 'Exp fit (MM36)'),
                # 'Gamma_234': (ND_MM234.loc[time_to_plot], 'Gamma fit (MM234)'),
                'Gamma_246': (ND_MM246.loc[time], 'Gamma fit (MM246)'),
                # Gamma_346': (ND_MM346.loc[time_to_plot], 'Gamma fit (MM346)')
                'TruncGamma_246': (ND_TMM246.loc[time], 'Truncated Gamma fit (TMM246)'),
            }
            PSDparamdict = {}
            # FIXME
            # PSDparamdict = {'Shape_gam': (mu_gam[t], r'$\mu$ (gamma)'),
            #                 'Slope_gam': (lamda_gam[t], r'$\lambda$ (gamma)'),
            #                 'D0_gam': (D_med_gam[t], r'$D_0$ (gamma, mm)'),
            #                 'Slope_exp': (lamda_exp[t], r'$\lambda$ (exp)'),
            #                 'N0_exp': (N0_exp[t], r'$N_0$ (exp, m$^{-4}$)'),
            #                 'D0_dis': (D_med_disd[t], r'$D_0$ (obs, mm)'),
            #                 'dBZ_dis': (refl_disd[t], r'Z (obs, dBZ)'),
            #                 'RR_disd': (PSD_df['intensity'].values[t],
            #                             r'Rain rate (obs, mm hr$^{-1}$)'),
            #                 'Particle count': (PSD_df[pcountstr].values[t],
            #                                     r'Particle count (QC)')}

            # if pc.calc_dualpol:
            #     PSDparamdict['ZDR_dis'] = (dualpol_dis['ZDR'][t], r'$Z_{DR}$ (obs, dB)')
            fig, ax = pm.plot_DSD(axdict, PSDdict, PSDfitdict, PSDparamdict)
            image_name = PIPS_name + '_DSD_{:d}_{}_t{:04d}.png'.format(int(DSD_interval),
                                                                       time.strftime(tm.timefmt3),
                                                                       t)
            image_path = os.path.join(DSD_image_dir, image_name)
            fig.savefig(image_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
