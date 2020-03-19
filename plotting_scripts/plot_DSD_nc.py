# pyPIPS_meteograms.py
#
# This script plots DSD histograms from the PIPS (netCDF version)
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
parser.add_argument('--plot-raw', action='store_true', dest='plot_raw',
                    help='plot raw DSD')
parser.add_argument('--plot-qc', action='store_true', dest='plot_qc',
                    help='plot QC DSD')
parser.add_argument('--plot-full', action='store_true', dest='plot_full',
                    help='Plot full-deployment DSD')
parser.add_argument('--plot-series', action='store_true', dest='plot_series',
                    help='Plot time series of DSDs')

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

radar_name = config.radar_config_dict.get('radar_name', None)

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
    ND_list = []
    tag_list = []

    if args.plot_raw:
        ND_raw = parsivel_combined_ds['ND']
        ND_list.append(ND_raw)
        tag_list.append('raw')
    if args.plot_qc:
        ND_qc = parsivel_combined_ds['ND_qc']
        ND_list.append(ND_qc)
        tag_list.append('qc')

    for tag, ND in zip(tag_list, ND_list):

        DSD_MM24 = parsivel_combined_ds['DSD_MM24']
        DSD_MM246 = parsivel_combined_ds['DSD_MM246']
        DSD_TMM246 = parsivel_combined_ds['DSD_TMM246']
        DSD_retr_mu = parsivel_combined_ds['mu_retr']
        DSD_retr_N0 = parsivel_combined_ds['N0_retr'] * 1000**(1 + DSD_retr_mu)  # Get to m^-4
        DSD_retr_lamda = parsivel_combined_ds['lamda_retr'] * 1000.  # Get to m^-1
        rad_dim_name = 'fields_{}'.format(radar_name)
        rad_fields_key = '{}_at_PIPS'.format(radar_name)
        DSD_rad_mu = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'mu'}]
        DSD_rad_N0 = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'N0'}] * \
                     1000**(1 + DSD_rad_mu)
        DSD_rad_lamda = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'lamda'}] * 1000.

        ND_MM24 = dsd.calc_binned_DSD_from_params(DSD_MM24.loc['N0'], DSD_MM24.loc['lamda'], 0.,
                                                  ND['diameter'])
        ND_MM246 = dsd.calc_binned_DSD_from_params(DSD_MM246.loc['N0'], DSD_MM246.loc['lamda'],
                                                   DSD_MM246.loc['alpha'], ND['diameter'])
        ND_TMM246 = dsd.calc_binned_DSD_from_params(DSD_TMM246.loc['N0'], DSD_TMM246.loc['lamda'],
                                                    DSD_TMM246.loc['alpha'], ND['diameter'])
        ND_retr = dsd.calc_binned_DSD_from_params(DSD_retr_N0, DSD_retr_lamda, DSD_retr_mu,
                                                  ND['diameter'])
        ND_rad = dsd.calc_binned_DSD_from_params(DSD_rad_N0, DSD_rad_lamda, DSD_rad_mu,
                                                 ND['diameter'])
        ND_onedrop = pips.calc_ND_onedrop(DSD_interval, correct_rho=True,
                                          rho=parsivel_combined_ds['rho'])

        D0 = dsd.calc_D0_bin(ND) * 1000. # Get to mm
        D0_retr = parsivel_combined_ds['D0_retr']
        D0_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'D0'}]

        DSD_image_dir = os.path.join(plot_dir, 'DSDs/{}_{}_{:d}s'.format(deployment_name, PIPS_name,
                                                                         int(requested_interval)))
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

        if args.plot_series:
            for t, time in enumerate(parsivel_combined_ds['time'].to_index()):
                if parsivel_combined_ds['pcount'].loc[time] > 0:
                    print("Plotting for {} and time {}".format(PIPS_name,
                                                               time.strftime(tm.timefmt3)))
                    axdict['time'] = time

                    PSDdict = {
                        'ND': ND.loc[time],
                        'ND_onedrop': ND_onedrop.loc[time]
                    }
                    PSDfitdict = {
                        'Exponential_24': (ND_MM24.loc[time], 'Exp. fit (MM24)'),
                        # 'Exponential_36': (ND_MM36.loc[time_to_plot], 'Exp fit (MM36)'),
                        # 'Gamma_234': (ND_MM234.loc[time_to_plot], 'Gamma fit (MM234)'),
                        'Gamma_246': (ND_MM246.loc[time], 'Gamma fit (MM246)'),
                        # Gamma_346': (ND_MM346.loc[time_to_plot], 'Gamma fit (MM346)')
                        'TruncGamma_246': (ND_TMM246.loc[time], 'Truncated Gamma fit (TMM246)'),
                        'DSD_radar': (ND_rad.loc[time], 'Retrieved DSD (radar)'),
                        'DSD_retr': (ND_retr.loc[time], 'Retrieved DSD (disdrometer)')
                    }
                    PSDparamdict = {
                        'N0_gamma_TMM246': (DSD_TMM246.loc['N0'][t], r'$N_0$ (TMM, m$^{-4}$)'),
                        'lamda_gamma_TMM246': (DSD_TMM246.loc['lamda'][t],
                                               r'$\lambda$ (TMM, m$^{-1}$)'),
                        'mu_gamma_TMM246': (DSD_TMM246.loc['alpha'][t], r'$\mu$ (TMM)'),
                        'D0': (D0[t], r'$D_0$ (obs, mm)'),
                        'N0_rad': (DSD_rad_N0[t], r'$N_0$ (radar, m$^{-4}$)'),
                        'lamda_rad': (DSD_rad_lamda[t], r'$\lambda$ (radar, m$^{-1}$)'),
                        'mu_rad': (DSD_rad_mu[t], r'$\mu$ (radar)'),
                        'D0_rad': (D0_rad[t], r'$D_0$ (radar, mm)'),
                        'N0_retr': (DSD_retr_N0[t], r'$N_0$ (retr obs, m$^{-4}$)'),
                        'lamda_retr': (DSD_retr_lamda[t], r'$\lambda$ (retr obs, m$^{-1}$)'),
                        'mu_retr': (DSD_retr_mu[t], r'$\mu$ (retr obs)'),
                        'D0_retr': (D0_retr[t], r'$D_0$ (retr obs, mm)')
                    }
                    fig, ax = pm.plot_DSD(axdict, PSDdict, PSDfitdict, PSDparamdict)
                    image_name = \
                        '{}_{}_DSD_{}_{:d}_{}_t{:04d}.png'.format(PIPS_name, deployment_name, tag,
                                                                  int(DSD_interval),
                                                                  time.strftime(tm.timefmt3), t)
                    image_path = os.path.join(DSD_image_dir, image_name)
                    fig.savefig(image_path, dpi=200, bbox_inches='tight')
                    plt.close(fig)

        # if args.plot_full:
        #     vd_matrix_da_full = vd_matrix_da.sum(dim='time')
        #     print("Plotting full-deployment v-d matrix for {} and {}".format(deployment_name,
        #                                                                      PIPS_name))
        #     axdict['cblim'] = (1, 1500)
        #     PSDdict = {
        #         'vd_matrix_da': vd_matrix_da_full,
        #         'DSD_interval': len(vd_matrix_da['time']) * DSD_interval
        #         # FIXME
        #         # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
        #     }
        #     fig, ax = pm.plot_vel_D(axdict, PSDdict, parsivel_combined_ds['rho'].mean(dim='time'))
        #     image_name =  \
        #         '{}_{}_vel_D_{}_full.png'.format(PIPS_name, deployment_name, tag)
        #     image_path = os.path.join(vel_D_image_dir, image_name)
        #     fig.savefig(image_path, dpi=200, bbox_inches='tight')
        #     plt.close(fig)

