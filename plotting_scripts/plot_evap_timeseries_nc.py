# plot_evap_timeseries_nc.py
#
# Plots time series of surface rain evaporation and latent cooling rates computed by
# analysis_scripts/calc_evap_PIPS.py: binned (observed) vs binned + truncated-gamma extension vs
# bulk-scheme closures, closure-to-binned ratios, and the small-drop contribution to the
# evaporation rate. If time-integration output is present, also plots the integrated
# (rain-evaporation-only) temperature and water vapor changes against the observed changes. If
# rainshaft output is present, also plots time-height evaporation and
# cooling-rate sections with the PIPS surface series beneath.
#
# NOTE: the cooling rates are single-process (rain evaporation) tendencies, not a closed surface
# heat budget.
from __future__ import annotations

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import colors

import pyPIPS.plotmodule as pm
from pyPIPS import utils

mpl.rc('font', size=9)

description = "Plots surface evaporation/cooling-rate time series from calc_evap_PIPS.py output"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path', default='plot_config.py',
                    help='Location of the plot configuration file')
parser.add_argument('--QC-tags', dest='QC_tags', nargs='*', default=['qc'],
                    help="QC tag(s) to plot (e.g., qc, roqc). Use '' for raw ND.")
parser.add_argument('--input-tag', dest='input_tag', default='',
                    help='Output tag used when running calc_evap_PIPS.py')
parser.add_argument('--input-dir', dest='input_dir', default=None,
                    help='Directory of the evap_*.nc files (default: PIPS_dir)')
parser.add_argument('--plot-dir', dest='plot_dir', default=None,
                    help='Output directory for images (default: <plot_dir>/evaporation)')
parser.add_argument('--closures', dest='closures', nargs='*', default=None,
                    help='Closures to plot (default: all in the file)')
parser.add_argument('--image-format', dest='image_format', default='png')

args = parser.parse_args()

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal("Unable to import case configuration parameters! Aborting!")

# Dynamically import the plotting configuration file
try:
    pc = utils.import_all_from(args.plot_config_path)
    utils.log("Successfully imported plotting configuration parameters!")
except Exception:
    utils.warning("Unable to import user-defined plotting parameters! Reverting to defaults.")
    import configs.plot_config_default as pc

PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)
input_dir = args.input_dir or PIPS_dir
image_root = args.plot_dir or os.path.join(plot_dir, 'evaporation')

xformatter = pc.PIPS_plotting_dict.get('majorxformatter', None)
xlocator = pc.PIPS_plotting_dict.get('majorxlocator', None)
xminorlocator = pc.PIPS_plotting_dict.get('minorxlocator', None)
xlabel = pc.PIPS_plotting_dict.get('xlabel', 'Time (HH:MM)')


def axparams(ylabel, ylim=None):
    return {
        'majorxlocator': xlocator,
        'majorxformatter': xformatter,
        'minorxlocator': xminorlocator,
        'axeslimits': [None, ylim],
        'axeslabels': [None, ylabel],
        'axesautofmt': False,
    }


def plot_surface_timeseries(ds, tag, closures, title, outfile):
    """Four-panel surface evaporation time series for one QC tag."""
    times = ds['time'].to_numpy()
    rain = ds[f'E_binned{tag}'] > 0.
    fig, axes = plt.subplots(4, 1, figsize=(8, 10), sharex=True)

    # (a) Latent cooling rate
    ax = axes[0]
    ax.plot(times, ds[f'cooling_rate_binned_K_per_h{tag}'], 'k-', lw=1.5, label='binned')
    if f'cooling_rate_binned_ext_K_per_h{tag}' in ds:
        ax.plot(times, ds[f'cooling_rate_binned_ext_K_per_h{tag}'], 'k--', lw=1.,
                label='binned + gamma ext.')
    cool_closure = 3600. * ds[f'cooling_rate_closure{tag}']
    cmap = plt.get_cmap('tab10')
    for i, closure in enumerate(closures):
        ax.plot(times, cool_closure.sel(closure=closure).where(rain), lw=0.8, color=cmap(i),
                label=closure)
    ax.axhline(0., color='0.6', lw=0.5)
    ax.legend(fontsize=7, ncol=3, loc='lower left')
    ax.set_title(title)

    # (b) Closure-to-binned ratio
    ax = axes[1]
    ratio = ds[f'ratio_to_binned_closure{tag}']
    for i, closure in enumerate(closures):
        ax.plot(times, ratio.sel(closure=closure), lw=0.8, color=cmap(i))
    ax.axhline(1., color='k', lw=0.8)
    ax.set_yscale('log')

    # (c) Cumulative fraction of E from drops smaller than D
    ax = axes[2]
    cumfrac = ds[f'E_cumfrac{tag}'].where(rain)
    edges = np.append(ds['min_diameter'].to_numpy(), ds['max_diameter'].to_numpy()[-1])
    pcm = ax.pcolormesh(times, edges[1:], cumfrac.T.to_numpy(), cmap='viridis', vmin=0., vmax=1.,
                        shading='nearest')
    ax.contour(times, edges[1:], cumfrac.T.fillna(0.).to_numpy(), levels=[0.5], colors='w',
               linewidths=0.8)
    ax.set_ylim(0., 4.)
    # Colorbar in an inset outside the axes so all panels keep the same width
    fig.colorbar(pcm, cax=ax.inset_axes([1.01, 0., 0.015, 1.]),
                 label='fraction of E\nfrom drops < D')

    # (d) Context: rain mixing ratio and RH
    ax = axes[3]
    ax.plot(times, 1000. * ds[f'qr_obs{tag}'], 'b-', lw=1.)
    ax2 = ax.twinx()
    ax2.plot(times, 100. * ds['RH_used'], 'g-', lw=1.)
    ax2.set_ylabel(f"RH ({ds.attrs.get('RH_var', '')}; %)", color='g')
    if f'frac_E_sub{tag}' in ds:
        ax2.plot(times, 100. * ds[f'frac_E_sub{tag}'], 'm:', lw=1.,
                 label='% of E from D < D_cut (ext.)')
        ax2.legend(fontsize=7, loc='upper right')

    axparamdicts = [axparams(r'Cooling rate (K h$^{-1}$)'),
                    axparams('E(closure) / E(binned)'),
                    axparams('D (mm)'),
                    axparams(r'$q_r$ (g kg$^{-1}$)')]
    pm.set_meteogram_axes(list(axes), axparamdicts)
    axes[-1].set_xlabel(xlabel)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    fig.text(0.01, 0.005, 'Rain-evaporation tendency only; not a closed surface heat budget.',
             fontsize=7, color='0.4')
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_integrations(ds, tag, closures, title, outfile):
    """Integrated temperature, water vapor, and RH changes vs observed, for the static and
    dynamic integrations of the binned, extended, and closure estimates."""
    times = ds['time'].to_numpy()
    modes = [m for m in ('static', 'dynamic') if f'dT_{m}_binned{tag}' in ds]
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    styles = {'static': '--', 'dynamic': '-'}
    cmap = plt.get_cmap('tab10')
    axes[0].plot(times, ds[f'dT_obs_integration{tag}'], color='0.5', lw=2.,
                 label='observed (all processes)')
    axes[1].plot(times, 1000. * ds[f'dw_obs_integration{tag}'], color='0.5', lw=2.)
    axes[2].plot(times, 100. * ds['RH_used'].where(ds[f'dT_obs_integration{tag}'].notnull()),
                 color='0.5', lw=2.)
    for mode in modes:
        ls = styles[mode]
        for source, color in (('binned', 'k'), ('binned_ext', 'm')):
            if f'dT_{mode}_{source}{tag}' not in ds:
                continue
            n_miss = int(ds[f'n_missing_{mode}_{source}{tag}'])
            label = f'{mode} {source}' + (f' ({n_miss} missing)' if n_miss else '')
            axes[0].plot(times, ds[f'dT_{mode}_{source}{tag}'], ls=ls, color=color, lw=1.2,
                         label=label)
            axes[1].plot(times, 1000. * ds[f'dw_{mode}_{source}{tag}'], ls=ls, color=color,
                         lw=1.2)
            axes[2].plot(times, 100. * ds[f'RH_{mode}_{source}{tag}'], ls=ls, color=color,
                         lw=1.2)
        for i, closure in enumerate(closures):
            label = f'{mode} {closure}' if mode == modes[-1] else None
            axes[0].plot(times, ds[f'dT_{mode}_closure{tag}'].sel(closure=closure), ls=ls,
                         color=cmap(i), lw=0.7, label=label)
    axes[0].legend(fontsize=6, ncol=3, loc='lower left')
    axes[0].set_title(f"{title}\nwindow {ds.attrs.get('integration_start', '')[:19]} to "
                      f"{ds.attrs.get('integration_end', '')[:19]}; dashed: static, "
                      "solid: dynamic", fontsize=9)
    axparamdicts = [axparams(r'$\Delta T$ (K)'), axparams(r'$\Delta w$ (g kg$^{-1}$)'),
                    axparams('RH (%)')]
    pm.set_meteogram_axes(list(axes), axparamdicts)
    axes[-1].set_xlabel(xlabel)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    fig.text(0.01, 0.005, 'Integrated changes: rain evaporation only, same air assumed to stay at '
                          'the PIPS. Observed changes include all processes.', fontsize=7,
             color='0.4')
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


def plot_rainshaft_sections(ds, tag, title, outfile):
    """Time-height rainshaft E and cooling rate with the PIPS surface series beneath."""
    times = ds['time'].to_numpy()
    height = ds['height'].to_numpy()
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True,
                             gridspec_kw={'height_ratios': [2, 2, 1]})
    E = ds[f'rainshaft_E{tag}'] * 1.e6
    pcm = axes[0].pcolormesh(times, height, E.T.to_numpy(), cmap='YlOrRd', shading='nearest',
                             norm=colors.LogNorm(vmin=1.e-2, vmax=1.e1))
    fig.colorbar(pcm, cax=axes[0].inset_axes([1.01, 0., 0.015, 1.]),
                 label=r'E (mg m$^{-3}$ s$^{-1}$)')
    cool = 3600. * ds[f'rainshaft_cooling_rate{tag}']
    pcm = axes[1].pcolormesh(times, height, cool.T.to_numpy(), cmap='Blues_r', vmin=-10.,
                             vmax=0., shading='nearest')
    fig.colorbar(pcm, cax=axes[1].inset_axes([1.01, 0., 0.015, 1.]),
                 label=r'Cooling (K h$^{-1}$)')
    for ax in axes[:2]:
        ax.plot(times, ds[f'rainshaft_z_top{tag}'], 'k.', ms=1.5)
    axes[2].plot(times, 1.e6 * ds[f'E_binned{tag}'], 'k-', lw=1., label='PIPS binned (obs.)')
    axes[2].plot(times, 1.e6 * ds[f'rainshaft_E_surface{tag}'], 'r.', ms=3.,
                 label='rainshaft (surface)')
    axes[2].legend(fontsize=7)
    axparamdicts = [axparams('Height above PIPS (m)'), axparams('Height above PIPS (m)'),
                    axparams(r'E (mg m$^{-3}$ s$^{-1}$)')]
    pm.set_meteogram_axes(list(axes), axparamdicts)
    axes[0].set_title(title + '\n' + ds.attrs.get('rainshaft_thermo_profile_source', ''))
    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    plt.close(fig)


for deployment_name, PIPS_name in zip(deployment_names, PIPS_names):
    fname = (f'evap_{deployment_name}_{PIPS_name}_{int(requested_interval):d}s'
             f'{"_" + args.input_tag if args.input_tag else ""}.nc')
    evap_file = os.path.join(input_dir, fname)
    if not os.path.exists(evap_file):
        utils.log(f"{evap_file} not found; skipping")
        continue
    utils.log(f"Reading {evap_file}")
    ds = xr.load_dataset(evap_file)
    image_dir = os.path.join(image_root, deployment_name)
    os.makedirs(image_dir, exist_ok=True)

    for QC_tag_tmp in args.QC_tags:
        tag = f'_{QC_tag_tmp}' if QC_tag_tmp else ''
        if f'E_binned{tag}' not in ds:
            utils.log(f"  No E_binned{tag} in file; skipping")
            continue
        closures = args.closures or list(ds['closure'].to_numpy())
        title = f'{PIPS_name} {deployment_name} ND{tag}: surface rain evaporation'
        outfile = os.path.join(image_dir, f'evap_timeseries_{deployment_name}_{PIPS_name}{tag}.'
                                          f'{args.image_format}')
        plot_surface_timeseries(ds, tag, closures, title, outfile)
        utils.log(f"  Saved {outfile}")
        if f'dT_obs_integration{tag}' in ds:
            outfile = os.path.join(image_dir, f'evap_integration_{deployment_name}_{PIPS_name}'
                                              f'{tag}.{args.image_format}')
            plot_integrations(ds, tag, closures, f'{PIPS_name} {deployment_name} ND{tag}',
                              outfile)
            utils.log(f"  Saved {outfile}")
        if f'rainshaft_E{tag}' in ds:
            outfile = os.path.join(image_dir, f'evap_rainshaft_{deployment_name}_{PIPS_name}'
                                              f'{tag}.{args.image_format}')
            plot_rainshaft_sections(ds, tag, f'{PIPS_name} {deployment_name}', outfile)
            utils.log(f"  Saved {outfile}")
