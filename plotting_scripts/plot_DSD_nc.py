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
parser.add_argument('--ND-tags', dest='ND_tags', nargs='*', default=None,
                    help='Tags for ND variable sin file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--plot-full', action='store_true', dest='plot_full',
                    help='Plot full-deployment DSD')
parser.add_argument('--plot-series', action='store_true', dest='plot_series',
                    help='Plot time series of DSDs')
parser.add_argument('--retr-tag', dest='retr_tag', default='SATP',
                    help='string to identify retrieval to plot')
parser.add_argument('--image-fmt', dest='image_fmt', default='png',
                    help='Image file format (i.e. png, eps, pdf)')

args = parser.parse_args()
if not args.ND_tags:
    ND_tags = ['']
else:
    ND_tags = ['_{}'.format(tag) for tag in args.ND_tags]

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
scatt_dir = config.radar_config_dict.get('scatt_dir', None)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
# print(parsivel_combined_filenames)

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    ND_list = [parsivel_combined_ds['ND{}'.format(ND_tag)] for ND_tag in ND_tags]

    for ND_tag, ND in zip(ND_tags, ND_list):
        DSD_MM24 = parsivel_combined_ds['DSD_MM24{}'.format(ND_tag)]
        DSD_MM246 = parsivel_combined_ds['DSD_MM246{}'.format(ND_tag)]
        DSD_TMM246 = parsivel_combined_ds['DSD_TMM246{}'.format(ND_tag)]
        DSD_retr_mu = parsivel_combined_ds['mu_retr_{}{}'.format(args.retr_tag, ND_tag)]
        DSD_retr_N0 = (parsivel_combined_ds['N0_retr_{}{}'.format(args.retr_tag, ND_tag)] *
                       1000**(1 + DSD_retr_mu))  # Get to m^-4
        DSD_retr_lamda = (parsivel_combined_ds['lamda_retr_{}{}'.format(args.retr_tag, ND_tag)] *
                          1000.)
        rad_dim_name = 'fields_{}'.format(radar_name)
        rad_fields_key = '{}_at_PIPS'.format(radar_name)
        # Annoying.. if the ND tag is 'qc', the radar fields don't have it as a suffix,
        # so remove it here. Also SATP_TMM is called just "SATP" here. Grumble.
        if ND_tag == '_qc':
            ND_rad_tag = ''
        else:
            ND_rad_tag = ND_tag
        if args.retr_tag == 'SATP_TMM':
            rad_retr_tag = 'SATP'
        else:
            rad_retr_tag = args.retr_tag
        DSD_rad_mu = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                               'mu_{}{}'.format(rad_retr_tag,
                                                                                ND_rad_tag)}]
        DSD_rad_N0 = (parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                                'N0_{}{}'.format(rad_retr_tag,
                                                                                 ND_rad_tag)}] *
                      1000.**(1. + DSD_rad_mu))
        DSD_rad_lamda = \
            (parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                       'lamda_{}{}'.format(rad_retr_tag,
                                                                           ND_rad_tag)}]
             * 1000.)

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

        # Commented out for now. Use Dm43 instead
        # D0 = dsd.calc_D0_bin(ND) * 1000.  # Get to mm
        # D0_retr = parsivel_combined_ds['D0_retr_{}'.format(args.retr_tag)]
        # D0_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'D0'}]

        Dm = parsivel_combined_ds['Dm43{}'.format(ND_tag)] * 1000.  # Get to mm
        Dm_retr = parsivel_combined_ds['Dm43_retr_{}{}'.format(args.retr_tag, ND_tag)]
        Dm_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                           'Dm_{}{}'.format(rad_retr_tag,
                                                                            ND_rad_tag)}]

        dualpol_dict = dp.calpolrain(10.7, os.path.join(scatt_dir, 'SCTT_RAIN_fw100.dat'), ND,
                                     bin_width)
        dBZ = dualpol_dict['REF']
        dBZ_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'REF_filtered'}]
        ZDR = dualpol_dict['ZDR']
        ZDR_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'ZDR_filtered'}]

        DSD_image_dir = os.path.join(plot_dir,
                                     'DSDs/{}_{}_{:d}s_{}'.format(deployment_name, PIPS_name,
                                                                  int(requested_interval),
                                                                  radar_name))
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
                if parsivel_combined_ds['pcount_derived_qc'].loc[time] > 0:
                    print("Plotting for {} and time {}".format(PIPS_name,
                                                               time.strftime(tm.timefmt3)))
                    axdict['time'] = time

                    print('N0, lambda, mu (radar): ', DSD_rad_N0[t], DSD_rad_lamda[t],
                          DSD_rad_mu[t])
                    print('ND (radar)', ND_rad.loc[time])

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
                        'DSD_radar': (ND_rad.loc[time],
                                      'Retrieved DSD (radar; {})'.format(args.retr_tag)),
                        'DSD_retr': (ND_retr.loc[time],
                                     'Retrieved DSD (disdrometer; {})'.format(args.retr_tag))
                    }
                    PSDparamdict = {
                        'dBZ': (dBZ[t], r'$Z_H$ (PIPS, dBZ)'),
                        'dBZ_rad': (dBZ_rad[t], r'$Z_H$ (radar, dBZ)'),
                        'ZDR': (ZDR[t], r'$Z_{DR}$ (PIPS, dB)'),
                        'ZDR_rad': (ZDR_rad[t], r'$Z_{DR}$ (radar, dB)'),
                        # 'N0_gamma_TMM246': (DSD_TMM246.loc['N0'][t], r'$N_0$ (TMM, m$^{-4}$)'),
                        # 'lamda_gamma_TMM246': (DSD_TMM246.loc['lamda'][t],
                        #                      r'$\lambda$ (TMM, m$^{-1}$)'),
                        # 'mu_gamma_TMM246': (DSD_TMM246.loc['alpha'][t], r'$\mu$ (TMM)'),
                        'Dm': (Dm[t], r'$D_m$ (PIPS, mm)'),
                        # 'N0_rad': (DSD_rad_N0[t], r'$N_0$ (radar, m$^{-4}$)'),
                        # 'lamda_rad': (DSD_rad_lamda[t], r'$\lambda$ (radar, m$^{-1}$)'),
                        # 'mu_rad': (DSD_rad_mu[t], r'$\mu$ (radar)'),
                        'Dm_rad': (Dm_rad[t], r'$D_m$ (radar, mm; {})'.format(args.retr_tag)),
                        # 'N0_retr': (DSD_retr_N0[t], r'$N_0$ (retr obs, m$^{-4}$)'),
                        # 'lamda_retr': (DSD_retr_lamda[t], r'$\lambda$ (retr obs, m$^{-1}$)'),
                        # 'mu_retr': (DSD_retr_mu[t], r'$\mu$ (retr obs)'),
                        'Dm_retr': (Dm_retr[t], r'$D_m$ (PIPS, mm; {})'.format(args.retr_tag))
                    }
                    fig, ax = pm.plot_DSD(axdict, PSDdict, PSDfitdict, PSDparamdict)
                    image_name = \
                        '{}_{}_DSD_{}_{:d}_{}_{}_t{:04d}.{}'.format(PIPS_name, deployment_name,
                                                                    ND_tag, int(DSD_interval),
                                                                    radar_name,
                                                                    time.strftime(tm.timefmt3), t,
                                                                    args.image_fmt)
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

