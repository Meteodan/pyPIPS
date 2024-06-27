# pyPIPS_meteograms.py
#
# This script plots DSD histograms from the PIPS (netCDF version)
from __future__ import annotations

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import xarray as xr

import pyPIPS.DSDlib as dsd
import pyPIPS.parsivel_params as pp
import pyPIPS.PIPS as pips
import pyPIPS.plotmodule as pm
import pyPIPS.polarimetric as dp
import pyPIPS.timemodule as tm
from pyPIPS import utils

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
                    help='Tags for ND variables in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--plot-raw', action='store_true', dest='plot_raw', default=False,
                    help='whether to plot the raw DSDs in addition to QC versions')
parser.add_argument('--plot-full', action='store_true', dest='plot_full',
                    help='Plot full-deployment DSD')
parser.add_argument('--plot-series', action='store_true', dest='plot_series',
                    help='Plot time series of DSDs')
parser.add_argument('--plot-MM-fits', action='store_true', dest='plot_MM_fits',
                    help='Plot MM fits and parameters?')
parser.add_argument('--plot-retr', action='store_true', dest='plot_retr',
                    help='Plot retrieval fits and parameters?')
parser.add_argument('--retr-tag', dest='retr_tag', default='SATP_TMM',
                    help='string to identify retrieval to plot')
parser.add_argument('--image-fmt', dest='image_fmt', default='png',
                    help='Image file format (i.e. png, eps, pdf)')
parser.add_argument('--time-dim', dest='time_dim', default='time',
                    help='Name of the time dimension in the netCDF file')

args = parser.parse_args()
plot_raw = args.plot_raw
if not args.ND_tags:
    ND_tags = ['']
else:
    ND_tags = [f'_{tag}' for tag in args.ND_tags]
    if plot_raw:
        ND_tags.insert(0, '')

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Dynamically import the plotting configuration file
utils.log(f"Plotting configuration file is {args.plot_config_path}")
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
PIPS_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

radar_name = config.radar_config_dict.get('radar_name', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
comp_radar = config.radar_config_dict.get('comp_radar', False)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
if PIPS_filenames_nc:
    parsivel_combined_filenames = PIPS_filenames_nc
else:
    parsivel_combined_filenames = [
        f'parsivel_combined_{deployment_name}_{PIPS_name}_{int(requested_interval):d}s.nc'
        for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
# print(parsivel_combined_filenames)

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print(f"Reading {parsivel_combined_file}")  # noqa: T201
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    # Set up desired start and end times if they are not "None" and we are not using relative times.
    # Otherwise just use all times in each file
    start_time = start_times[index]
    end_time = end_times[index]
    if start_time is not None and end_time is not None and args.time_dim == 'time':
        start_datetime = datetime.strptime(start_time, tm.timefmt3)
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
        start_time = start_datetime.strftime(tm.timefmt2)
        end_time = end_datetime.strftime(tm.timefmt2)
        print(f"Extracting subset of times between {start_time} and {end_time}")  # noqa: T201
        parsivel_combined_ds = parsivel_combined_ds.sel({
            args.time_dim: slice(start_time, end_time)})

    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    ND_list = [parsivel_combined_ds[f'ND{ND_tag}'] for ND_tag in ND_tags]
    try:
        num_PIPS = int(parsivel_combined_ds.num_PIPS)
    except AttributeError:
        num_PIPS = 1
    ND_onedrop = pips.calc_ND_onedrop(DSD_interval, correct_rho=True,
                                      rho=parsivel_combined_ds['rho'],
                                      num_PIPS=num_PIPS,
                                      time_dim=args.time_dim)

    for ND_tag, ND in zip(ND_tags, ND_list):
        rad_dim_name = f'fields_{radar_name}'
        rad_fields_key = f'{radar_name}_at_PIPS'

        if args.plot_MM_fits:
            DSD_MM24 = parsivel_combined_ds[f'DSD_MM24{ND_tag}']
            DSD_MM246 = parsivel_combined_ds[f'DSD_MM246{ND_tag}']
            DSD_TMM246 = parsivel_combined_ds[f'DSD_TMM246{ND_tag}']
        if args.plot_retr:
            DSD_retr_mu = parsivel_combined_ds[f'mu_retr_{args.retr_tag}{ND_tag}']
            DSD_retr_N0 = (parsivel_combined_ds[f'N0_retr_{args.retr_tag}{ND_tag}'] *
                           1000**(1 + DSD_retr_mu))  # Get to m^-4
            DSD_retr_lamda = (parsivel_combined_ds[f'lamda_retr_{args.retr_tag}{ND_tag}'] *
                              1000.)
            # Annoying.. if the ND tag is 'qc', the radar fields don't have it as a suffix,
            # so remove it here. Also SATP_TMM is called just "SATP" here. Grumble.

            rad_retr_tag = 'SATP' if args.retr_tag == 'SATP_TMM' else args.retr_tag
            try:
                ND_rad_tag = ND_tag
                DSD_rad_mu = \
                    parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                              f'mu_{rad_retr_tag}{ND_rad_tag}'}]
            except KeyError:
                ND_rad_tag = ''
                DSD_rad_mu = \
                    parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                              f'mu_{rad_retr_tag}{ND_rad_tag}'}]

            DSD_rad_N0 = \
                (parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                           f'N0_{rad_retr_tag}{ND_rad_tag}'}] *
                 1000.**(1. + DSD_rad_mu))
            DSD_rad_lamda = \
                (parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                        f'lamda_{rad_retr_tag}{ND_rad_tag}'}]
                * 1000.)
        if args.plot_MM_fits:
            ND_MM24 = dsd.calc_binned_DSD_from_params(DSD_MM24.loc['N0'], DSD_MM24.loc['lamda'], 0.,
                                                      ND['diameter'])
            # print('N0, lambda (MM24): ', DSD_MM24.loc['N0'], DSD_MM24.loc['lamda'])
            # print('ND (MM24)', ND_MM24)
            ND_MM246 = dsd.calc_binned_DSD_from_params(DSD_MM246.loc['N0'], DSD_MM246.loc['lamda'],
                                                       DSD_MM246.loc['alpha'], ND['diameter'])
            # print('N0, lambda (MM246): ', DSD_MM246.loc['N0'], DSD_MM24.loc['lamda'],  # noqa: T201
            #       DSD_MM246.loc['alpha'])
            # print('ND (MM246)', ND_MM246)  # noqa: T201
            ND_TMM246 = dsd.calc_binned_DSD_from_params(DSD_TMM246.loc['N0'],
                                                        DSD_TMM246.loc['lamda'],
                                                        DSD_TMM246.loc['alpha'], ND['diameter'])
        if args.plot_retr:
            ND_retr = dsd.calc_binned_DSD_from_params(DSD_retr_N0, DSD_retr_lamda, DSD_retr_mu,
                                                      ND['diameter'])
            ND_rad = dsd.calc_binned_DSD_from_params(DSD_rad_N0, DSD_rad_lamda, DSD_rad_mu,
                                                     ND['diameter'])

        # Commented out for now. Use Dm43 instead
        # D0 = dsd.calc_D0_bin(ND) * 1000.  # Get to mm
        # D0_retr = parsivel_combined_ds['D0_retr_{}'.format(args.retr_tag)]
        # D0_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'D0'}]

        Dm = parsivel_combined_ds[f'Dm43{ND_tag}'] * 1000.  # Get to mm
        if args.plot_retr:
            Dm_retr = parsivel_combined_ds[f'Dm43_retr_{args.retr_tag}{ND_tag}']
            Dm_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name:
                                                              f'Dm43_{rad_retr_tag}{ND_rad_tag}'}]

        dualpol_dict = dp.calpolrain(10.7, os.path.join(scatt_dir, 'SCTT_RAIN_fw100.dat'), ND,
                                     bin_width)
        dBZ = dualpol_dict['REF']
        ZDR = dualpol_dict['ZDR']
        if comp_radar:
            dBZ_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'REF_filtered'}]
            ZDR_rad = parsivel_combined_ds[rad_fields_key].loc[{rad_dim_name: 'ZDR_filtered'}]

        DSD_image_dir = os.path.join(plot_dir,
                                     f'DSDs/{deployment_name}_{PIPS_name}_{int(requested_interval):d}s_{radar_name}')
        if not os.path.exists(DSD_image_dir):
            os.makedirs(DSD_image_dir)

        axdict = {
            'xbin_left': min_diameter,
            'xbin_mid': avg_diameter,
            'xbin_right': max_diameter,
            'xlim': pc.PIPS_plotting_dict['diameter_range'],  # (0.0, 9.0),
            'ylim': (10.**2., 10.**8.5),
            'PIPS_name': PIPS_name
        }

        if args.plot_series:
            for t, ltime in enumerate(parsivel_combined_ds[args.time_dim].to_index()):
                time = int(ltime) if args.time_dim == 'relative_time' else ltime
                if parsivel_combined_ds[f'pcount_derived{ND_tag}'].loc[time] > 0:
                    if args.time_dim == 'time':
                        time_string = time.strftime(tm.timefmt3)
                    else:
                        time_string = f'{time:04d}'
                    print(f"Plotting for {PIPS_name} and time {time_string}")  # noqa: T201
                    axdict['interval'] = int(DSD_interval)
                    axdict['time'] = time
                    # if args.plot_retr:
                    #     print('N0, lambda, mu (radar): ', DSD_rad_N0[t], DSD_rad_lamda[t],  # noqa: T201, RUF100, E501
                    #           DSD_rad_mu[t])
                    #     print('ND (radar)', ND_rad.loc[time])  # noqa: T201, RUF100

                    PSDdict = {
                        'ND': ND.loc[time],
                        'ND_onedrop': ND_onedrop.loc[time]
                    }
                    if args.plot_MM_fits:
                        PSDfitdict = {
                            'Exponential_24': (ND_MM24.loc[time], 'Exp. fit (MM24)'),
                            # 'Exponential_36': (ND_MM36.loc[time_to_plot], 'Exp fit (MM36)'),
                            # 'Gamma_234': (ND_MM234.loc[time_to_plot], 'Gamma fit (MM234)'),
                            'Gamma_246': (ND_MM246.loc[time], 'Gamma fit (MM246)'),
                            # Gamma_346': (ND_MM346.loc[time_to_plot], 'Gamma fit (MM346)')
                            'TruncGamma_246': (ND_TMM246.loc[time], 'Truncated Gamma fit (TMM246)'),
                        }
                    else:
                        PSDfitdict = {}
                    if args.plot_retr:
                        PSDfitdict['DSD_radar'] = \
                            (ND_rad.loc[time], f'Retrieved DSD (radar; {args.retr_tag})')
                        PSDfitdict['DSD_retr'] = \
                            (ND_retr.loc[time],
                             f'Retrieved DSD (disdrometer; {args.retr_tag})')

                    PSDparamdict = {
                        'dBZ': (dBZ[t], r'$Z_H$ (PIPS, dBZ)'),
                        'ZDR': (ZDR[t], r'$Z_{DR}$ (PIPS, dB)'),
                        'Dm': (Dm[t], r'$D_m$ (PIPS, mm)'),
                    }
                    if args.plot_MM_fits:
                        PSDparamdict['dBZ_rad'] = (dBZ_rad[t], r'$Z_H$ (radar, dBZ)')
                        PSDparamdict['ZDR_rad'] = (ZDR_rad[t], r'$Z_{DR}$ (radar, dB)')
                    if args.plot_retr:
                        PSDparamdict['Dm_rad'] = (Dm_rad[t],
                                                  rf'$D_m$ (radar, mm; {args.retr_tag})')
                        PSDparamdict['Dm_retr'] = (Dm_retr[t],
                                                   rf'$D_m$ (PIPS, mm; {args.retr_tag})')
                        # PSDparamdict['N0_rad'] = (DSD_rad_N0[t], r'$N_0$ (radar, m$^{-4}$)')
                        # PSDparamdict['lamda_rad'] = \
                        #     (DSD_rad_lamda[t], r'$\lambda$ (radar, m$^{-1}$)')
                        # PSDparamdict['mu_rad'] = (DSD_rad_mu[t], r'$\mu$ (radar)')
                        # 'N0_gamma_TMM246': (DSD_TMM246.loc['N0'][t], r'$N_0$ (TMM, m$^{-4}$)'),
                        # 'lamda_gamma_TMM246': (DSD_TMM246.loc['lamda'][t],
                        #                      r'$\lambda$ (TMM, m$^{-1}$)'),
                        # 'mu_gamma_TMM246': (DSD_TMM246.loc['alpha'][t], r'$\mu$ (TMM)'),
                        # 'N0_rad': (DSD_rad_N0[t], r'$N_0$ (radar, m$^{-4}$)'),
                        # 'lamda_rad': (DSD_rad_lamda[t], r'$\lambda$ (radar, m$^{-1}$)'),
                        # 'mu_rad': (DSD_rad_mu[t], r'$\mu$ (radar)'),
                        # 'N0_retr': (DSD_retr_N0[t], r'$N_0$ (retr obs, m$^{-4}$)'),
                        # 'lamda_retr': (DSD_retr_lamda[t], r'$\lambda$ (retr obs, m$^{-1}$)'),
                        # 'mu_retr': (DSD_retr_mu[t], r'$\mu$ (retr obs)'),
                    fig, ax = pm.plot_DSD(axdict, PSDdict, PSDfitdict, PSDparamdict,
                                          time_dim=args.time_dim)
                    if args.time_dim == 'time':
                        time_string = time.strftime(tm.timefmt3)
                    else:
                        time_string = f'{t:04d}'
                    image_name = (f'{PIPS_name}_{deployment_name}_DSD_{ND_tag}_'
                                  f'{int(DSD_interval):d}_{radar_name}_'
                                  f'{time_string}_t{t:04d}.{args.image_fmt}')
                    image_path = os.path.join(DSD_image_dir, image_name)
                    fig.savefig(image_path, dpi=200, bbox_inches='tight')
                    plt.close(fig)
        if args.plot_full:
            ND_full = ND.mean(dim=args.time_dim)
            print(f"Plotting full-deployment DSD for {deployment_name} and {PIPS_name}")  # noqa: T201
            # FIXME: Getting some absurd values for the total duration for some of the plots.
            # Find out why
            full_DSD_interval = len(ND[args.time_dim]) * DSD_interval
            axdict['interval'] = int(full_DSD_interval)
            ND_onedrop_full = \
                pips.calc_ND_onedrop(full_DSD_interval, correct_rho=True,
                                     rho=parsivel_combined_ds['rho'].mean(dim=args.time_dim),
                                     time_dim=args.time_dim)
            PSDdict = {
                'ND': ND_full,
                'ND_onedrop': ND_onedrop_full
            }
            # TODO: re-compute parameters for full-deployment DSD and plot them
            fig, ax = pm.plot_DSD(axdict, PSDdict, {}, {}, time_dim=args.time_dim)
            image_name = (f'{PIPS_name}_{deployment_name}_DSD_{ND_tag}_{int(full_DSD_interval):d}_'
                          f'{radar_name}_full.{args.image_fmt}')
            image_path = os.path.join(DSD_image_dir, image_name)
            fig.savefig(image_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
