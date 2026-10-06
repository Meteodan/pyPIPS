# calc_evap_PIPS.py
#
# Computes surface raindrop evaporation rates and latent cooling rates from PIPS DSDs and the
# collocated PIPS thermodynamic measurements (netCDF version):
#   - from the observed (binned) DSD, with and without a truncated-gamma extension below the
#     Parsivel detection limit, plus small-drop contribution diagnostics;
#   - from gamma DSDs built the way one-, two-, and three-moment bulk schemes would (and from the
#     MM/TMM gamma fits), with ratios to the binned rate;
#   - time integrations of the surface cooling and moistening over a chosen window (default:
#     the whole deployment): "static" (tendencies from the observed state accumulated) and
#     "dynamic" (evaporation recomputed from the observed DSD and the evolving T and water
#     vapor), for the binned, binned + extension, and closure estimates;
#   - optionally, a 1-D evaporation-only rainshaft from radar-retrieved DSDs above the PIPS to
#     the surface (requires radar retrievals from radar_to_PIPS.py and
#     calc_radar_retrieval_PIPS.py), with a well-mixed (ASSUMED) subcloud profile from the PIPS
#     surface state.
# Outputs go to separate netCDF files, one per PIPS: evap_<deployment>_<PIPS>_<N>s<tag>.nc.
#
# NOTE: the cooling rate is a single-process (rain evaporation) tendency at the surface. It is
# not a closed surface heat budget; advection and other terms are not included.
from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone

import numpy as np
import xarray as xr

import pyPIPS.DSDlib as dsd
import pyPIPS.evaporation as evap
import pyPIPS.rainshaft as rs
from pyPIPS import utils

# Parse the command line options
description = ("Calculates surface evaporation and latent cooling rates from PIPS DSDs, with "
               "bulk-closure comparisons and an optional rainshaft from radar (netCDF version)")
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--QC-tags', dest='QC_tags', nargs='*', default=['qc'],
                    help="QC tag(s) for the ND variables (e.g., qc, roqc). Use '' for raw ND.")
parser.add_argument('--input-file-tag', dest='input_file_tag', default='',
                    help='Tag inserted before .nc in the parsivel combined file names (e.g., '
                         'KLBB for files written by radar_to_PIPS.py with --output-tag KLBB)')
parser.add_argument('--output-tag', dest='output_tag', default='',
                    help='Tag appended to the output file names')
parser.add_argument('--output-dir', dest='output_dir', default=None,
                    help='Output directory (default: PIPS_dir from the case config)')
parser.add_argument('--T-var', dest='T_var', default='fasttemp',
                    help='Temperature variable (e.g., fasttemp or slowtemp)')
parser.add_argument('--p-var', dest='p_var', default='pressure', help='Pressure variable')
parser.add_argument('--RH-var', dest='RH_var', default='RH_derived',
                    help='Relative humidity variable (e.g., RH_derived or RH)')
parser.add_argument('--fits', dest='fits', nargs='*', default=['TMM246', 'MM246'],
                    help='Gamma fits to compare (names as in DSD_<fit><QC tag>, e.g., TMM246)')
parser.add_argument('--compute-missing-fits', dest='compute_missing_fits', action='store_true',
                    help='Compute requested fits that are not stored in the file (as in '
                         'calc_MM_fits.py); otherwise missing fits are skipped')
parser.add_argument('--extension-fit', dest='extension_fit', default='TMM246',
                    help='Fit used for the truncated-gamma extension below the detection limit '
                         "('none' to disable)")
parser.add_argument('--extension-fallback', dest='extension_fallback', action='store_true',
                    help='Where the extension is undefined (e.g., fit mu <= -2), fall back to the '
                         'measured-bins-only rate (flagged in extension_fallback_used)')
parser.add_argument('--D-detect-min', dest='D_detect_min', type=float, default=0.25,
                    help='Lower detection limit for the extension option (mm)')
parser.add_argument('--mu-fixed', dest='mu_fixed', type=float, nargs='*', default=[0., 2.5],
                    help='Fixed shape parameters for the two-moment closures')
parser.add_argument('--N0-1mom', dest='N0_1mom', type=float, default=8.e6,
                    help='Fixed intercept for the one-moment closure (m^-4)')
parser.add_argument('--thermo-form', dest='thermo_form', default='textbook',
                    choices=['textbook', 'milbrandt'])
parser.add_argument('--saturation-definition', dest='saturation_definition',
                    default='vapor_pressure', choices=['vapor_pressure', 'mixing_ratio'])
parser.add_argument('--ventilation', dest='ventilation', default=evap.DEFAULT_VENTILATION,
                    choices=list(evap.VENTILATION_COEFFS))
parser.add_argument('--closure-integration', dest='closure_integration', default='quadrature',
                    choices=['quadrature', 'closed_form'],
                    help='quadrature: same fall speed as the binned rate (isolates the closure); '
                         'closed_form: Ferrier fall speed, as a bulk scheme would compute')
parser.add_argument('--integration-modes', dest='integration_modes', nargs='*',
                    default=['static', 'dynamic'], choices=['static', 'dynamic'],
                    help='Time integrations to run (pass no values to skip)')
parser.add_argument('--integration-start', dest='integration_start', default=None,
                    help='Start of the integration window (e.g., 2025-06-06T23:50); default: '
                         'start of the record')
parser.add_argument('--integration-end', dest='integration_end', default=None,
                    help='End of the integration window; default: end of the record')
parser.add_argument('--missing-fill', dest='missing_fill', default='zero',
                    choices=['zero', 'last_good'],
                    help='Integration treatment of missing samples: zero tendency, or the last '
                         'known good sample')
parser.add_argument('--max-fill-samples', dest='max_fill_samples', type=int, default=None,
                    help="Maximum consecutive missing samples filled with 'last_good' "
                         "(default unlimited)")
parser.add_argument('--rainshaft', dest='rainshaft', action='store_true',
                    help='Run the rainshaft from radar-retrieved DSDs above each PIPS')
parser.add_argument('--radar-name', dest='radar_name', default=None,
                    help='Radar name (default: radar_config_dict["radar_name"])')
parser.add_argument('--retrieval-tag', dest='retrieval_tag', default='Z01',
                    help='Tag of the radar retrieval fields in <radar>_at_PIPS (e.g., Z01_qc)')
parser.add_argument('--rainshaft-stride', dest='rainshaft_stride', type=int, default=1,
                    help='Run the rainshaft every N valid times')
parser.add_argument('--rainshaft-dz', dest='rainshaft_dz', type=float, default=50.,
                    help='Rainshaft output level spacing (m)')

args = parser.parse_args()

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal("Unable to import case configuration parameters! Aborting!")

PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
radar_config_dict = getattr(config, 'radar_config_dict', {})
radar_name = args.radar_name or radar_config_dict.get('radar_name', None)


def tagged_filename(fname, tag):
    if not tag:
        return fname
    root, ext = os.path.splitext(fname)
    return f'{root}_{tag}{ext}'


def get_gamma_fit(ds, ND, fit, QC_tag):
    """Return (N0, lamda, mu) for a stored fit, computing it if requested and missing."""
    name = f'DSD_{fit}{QC_tag}'
    if name in ds:
        return evap.gamma_params_from_fit(ds[name])
    if not args.compute_missing_fits:
        utils.log(f"  Fit {name} not found; skipping (use --compute-missing-fits)")
        return None
    utils.log(f"  Computing fit {name} (as in calc_MM_fits.py)")
    combo = fit.replace('TMM', '').replace('MM', '')
    ND_pos = ND.where(ND > 0.)
    moments = [dsd.calc_moment_bin(ND_pos, moment=int(m))[0] for m in combo]
    if fit.startswith('TMM'):
        D_min, D_max = dsd.get_max_min_diameters(ND_pos.fillna(0.0), dim='time')
        N0, lamda, mu = dsd.fit_DSD_TMMXYZ(combo, moments, D_min, D_max)
    else:
        N0, lamda, mu = dsd.fit_DSD_MMXYZ(combo, moments)
    return N0, lamda, mu


def run_rainshaft(ds, thermo):
    """Rainshaft from radar-retrieved gamma DSDs at the beam height down to the PIPS."""
    dim = f'fields_{radar_name}'
    if f'{radar_name}_at_PIPS' not in ds or \
            f'N0_{args.retrieval_tag}' not in ds[f'{radar_name}_at_PIPS'][dim].values:
        utils.log(f"  No {radar_name} retrieval fields with tag '{args.retrieval_tag}'; "
                  "skipping rainshaft")
        return None
    start = rs.starting_conditions_from_radar_at_PIPS(ds, radar_name, args.retrieval_tag)
    valid = (start['N0'] > 0.) & np.isfinite(start['lamda']) & (start['z_top'] > 0.)
    itimes = np.nonzero(valid.to_numpy())[0][::args.rainshaft_stride]
    utils.log(f"  Running rainshaft for {itimes.size} times")
    z_max = float(start['z_top'].max())
    heights = np.arange(0., z_max + args.rainshaft_dz, args.rainshaft_dz)
    edges_start = rs.fine_diameter_edges_mm(10., 0.1)
    out_edges = rs.PARSIVEL_EDGES_MM
    nt, nz, nb = ds.sizes['time'], heights.size, out_edges.size - 1
    fields = ['E', 'cooling_rate', 'rainrate', 'reflectivity', 'Dm43', 'Dmv', 'Nt', 'qr']
    data = {f: np.full((nt, nz), np.nan) for f in fields}
    ND_sfc = np.full((nt, nb), np.nan)
    z_LCL = np.full(nt, np.nan)
    profile_md = {}
    for it in itimes:
        T_s, p_s, RH_s = (float(thermo[v].isel(time=it)) for v in ('T', 'p', 'RH'))
        if not np.isfinite([T_s, p_s, RH_s]).all():
            continue
        z_top = float(start['z_top'].isel(time=it))
        profile = rs.ThermoProfile.well_mixed_from_surface(T_s, p_s, min(RH_s, 1.), z_top + 100.)
        ND_top = rs.starting_dsd_from_gamma(float(start['N0'].isel(time=it)),
                                            float(start['lamda'].isel(time=it)),
                                            float(start['mu'].isel(time=it)), edges_start)
        z_out = np.sort(np.unique(np.append(heights[heights < z_top], z_top)))[::-1]
        col = rs.rainshaft(ND_top, edges_start, z_top, 0., profile, z_out=z_out, n_sub=4,
                           ventilation=args.ventilation, thermo_form=args.thermo_form,
                           saturation_definition=args.saturation_definition)
        for f in fields:
            data[f][it] = np.interp(heights, col['height'].values[::-1],
                                    col[f].values[::-1], right=np.nan)
        ND_sfc[it] = col['ND_surface'].values
        z_LCL[it] = profile.attrs['z_LCL']
        profile_md = {'thermo_profile_source': profile.source,
                      'thermo_profile_assumption': profile.attrs['assumption']}
    times = ds['time']
    out = xr.Dataset({f'rainshaft_{f}': (('time', 'height'), v) for f, v in data.items()},
                     coords={'time': times, 'height': heights})
    out['rainshaft_ND_surface'] = (('time', 'diameter_bin'), ND_sfc)
    out['rainshaft_E_surface'] = out['rainshaft_E'].isel(height=0, drop=True)
    out['rainshaft_z_top'] = start['z_top']
    out['rainshaft_z_LCL'] = (('time',), z_LCL)
    out['rainshaft_z_LCL'].attrs.update(units='m', description='LCL height of the assumed '
                                                               'well-mixed profile (RH capped at 1 '
                                                               'above it)')
    for f in fields:
        units, desc = rs._RAINSHAFT_ATTRS[f]
        out[f'rainshaft_{f}'].attrs.update(units=units, description=f'rainshaft {desc}')
    out['rainshaft_ND_surface'].attrs.update(units='m-3 mm-1',
                                             description='rainshaft-predicted surface DSD on '
                                                         'Parsivel bins')
    out['height'].attrs.update(units='m', description='height above the PIPS')
    out.attrs.update({f'rainshaft_{k}': v for k, v in profile_md.items()})
    out.attrs.update(rainshaft_radar_name=radar_name, rainshaft_retrieval_tag=args.retrieval_tag,
                     rainshaft_start='radar-retrieved gamma DSD at the nearest gate center above '
                                     'the PIPS (radar values interpolated to PIPS times)')
    return out


def run_integrations(ND, binned, closures, thermo, missing, sample_interval, extension):
    """Static and dynamic time integrations for the binned, extended, and closure estimates.

    Returns a dict {output variable name (without QC tag): DataArray} and the window attrs.
    """
    kw = dict(start=args.integration_start, end=args.integration_end,
              sample_interval=sample_interval, missing=missing, missing_fill=args.missing_fill,
              max_fill_samples=args.max_fill_samples)
    T, p, RH, rho = thermo['T'], thermo['p'], thermo['RH'], thermo['rho']
    rate_kw = dict(ventilation=args.ventilation, thermo_form=args.thermo_form,
                   saturation_definition=args.saturation_definition)
    sources = {'binned': binned['E']}
    if 'E_ext' in binned:
        sources['binned_ext'] = binned['E_ext']
    sources['closure'] = closures['E']
    results = {}
    attrs = {}
    for mode in args.integration_modes:
        for source, E in sources.items():
            if mode == 'static':
                res = evap.integrate_static(E, T, p, RH, rho=rho, **kw)
            else:
                if source == 'binned':
                    rate = evap.binned_rate_function(ND, **rate_kw)
                elif source == 'binned_ext':
                    rate = evap.binned_rate_function(
                        ND, extension=extension, extension_fallback=args.extension_fallback,
                        D_detect_min=args.D_detect_min * 1.e-3, **rate_kw)
                else:
                    rate = evap.gamma_rate_function(
                        closures['N0'], closures['lamda'], closures['mu'],
                        integration=args.closure_integration,
                        no_rain=(binned['E'] == 0.).to_numpy(), **rate_kw)
                res = evap.integrate_dynamic(rate, T, p, RH, **kw)
            utils.log(f"    {mode} integration ({source}): "
                      f"final dT = {np.atleast_1d(res['dT'].isel(time=-1).to_numpy())} K")
            for var in ('dT', 'dw', 'T', 'w', 'RH', 'status', 'filled', 'E', 'cooling_rate'):
                if var in res and not (mode == 'static' and var in ('E', 'cooling_rate')):
                    results[f'{var}_{mode}_{source}'] = res[var]
            for var in res.data_vars:
                if var.startswith(('n_', 'duration_')):
                    results[f'{var}_{mode}_{source}'] = res[var]
            results['dT_obs_integration'] = res['dT_obs']
            results['dw_obs_integration'] = res['dw_obs']
            attrs = {f'integration_{k}' if not k.startswith('integration_') else k: v
                     for k, v in res.attrs.items() if k != 'integration_mode'}
    attrs['integration_modes'] = ' '.join(args.integration_modes)
    return results, attrs


for parsivel_combined_filename in parsivel_combined_filenames:
    parsivel_combined_file = os.path.join(PIPS_dir, tagged_filename(parsivel_combined_filename,
                                                                    args.input_file_tag))
    utils.log(f"Reading {parsivel_combined_file}")
    ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = ds.probe_name
    deployment_name = ds.deployment_name
    DSD_interval = ds.DSD_interval

    thermo = evap.thermo_from_parsivel_ds(ds, T_var=args.T_var, p_var=args.p_var,
                                          RH_var=args.RH_var)
    out_ds = xr.Dataset(coords={'time': ds['time']})
    out_ds['T_used'] = thermo['T']
    out_ds['p_used'] = thermo['p']
    out_ds['RH_used'] = thermo['RH']
    out_ds['rho_used'] = thermo['rho']

    for QC_tag_tmp in args.QC_tags:
        QC_tag = f'_{QC_tag_tmp}' if QC_tag_tmp else ''
        ND_name = f'ND{QC_tag}'
        if ND_name not in ds:
            utils.log(f"  {ND_name} not found; skipping")
            continue
        utils.log(f"  Working on {ND_name}")
        ND = ds[ND_name]

        fits = {}
        for fit in args.fits:
            params = get_gamma_fit(ds, ND, fit, QC_tag)
            if params is not None:
                fits[fit] = params
        extension = None
        if args.extension_fit.lower() != 'none':
            extension = fits.get(args.extension_fit)
            if extension is None:
                extension = get_gamma_fit(ds, ND, args.extension_fit, QC_tag)

        binned = evap.evaporation_rate_binned(
            ND, thermo['T'], thermo['p'], thermo['RH'], rho=thermo['rho'],
            ventilation=args.ventilation, thermo_form=args.thermo_form,
            saturation_definition=args.saturation_definition,
            D_detect_min=args.D_detect_min * 1.e-3, extension=extension,
            extension_fallback=args.extension_fallback)
        closures = evap.closure_comparison(
            ND, thermo['T'], thermo['p'], thermo['RH'], rho=thermo['rho'], binned=binned,
            integration=args.closure_integration, ventilation=args.ventilation,
            thermo_form=args.thermo_form, saturation_definition=args.saturation_definition,
            N0_1mom=args.N0_1mom, mu_fixed=tuple(args.mu_fixed), fits=fits)

        rename = {
            'E': 'E_binned', 'E_ext': 'E_binned_ext', 'dqr_dt': 'dqr_dt_binned',
            'dqr_dt_ext': 'dqr_dt_binned_ext', 'cooling_rate': 'cooling_rate_binned',
            'cooling_rate_K_per_h': 'cooling_rate_binned_K_per_h',
            'cooling_rate_ext': 'cooling_rate_binned_ext',
            'cooling_rate_ext_K_per_h': 'cooling_rate_binned_ext_K_per_h',
        }
        for var in binned.data_vars:
            out_ds[f'{rename.get(var, var)}{QC_tag}'] = binned[var]
        for var in ('E', 'ratio_to_binned', 'ratio_to_binned_ext', 'N0', 'lamda', 'mu'):
            if var in closures:
                out_ds[f'{var}_closure{QC_tag}'] = closures[var]
        tend = evap.latent_tendencies(closures['E'], thermo['rho'])
        out_ds[f'cooling_rate_closure{QC_tag}'] = tend['cooling_rate']
        out_ds[f'cooling_rate_closure{QC_tag}'].attrs.update(
            units='K s-1', description='closure latent cooling rate (negative = cooling)')
        for var in ('qr', 'Nt', 'Z', 'Dmv'):
            out_ds[f'{var}_obs{QC_tag}'] = closures[var]
        if f'flagged_times{QC_tag}' in ds:
            out_ds[f'flagged_times{QC_tag}'] = ds[f'flagged_times{QC_tag}']
        out_ds.attrs[f'D_cut_m{QC_tag}'] = binned.attrs.get('D_cut_m', np.nan)

        if args.integration_modes:
            missing = ds['missing_times'] > 0 if 'missing_times' in ds else None
            results, int_attrs = run_integrations(ND, binned, closures, thermo, missing,
                                                  float(DSD_interval), extension)
            for name, da in results.items():
                out_ds[f'{name}{QC_tag}'] = da
            out_ds.attrs.update(int_attrs)

        if args.rainshaft and radar_name is not None:
            rs_ds = run_rainshaft(ds, thermo)
            if rs_ds is not None:
                out_ds = out_ds.merge(rs_ds.rename(
                    {v: f'{v}{QC_tag}' for v in rs_ds.data_vars}))
                out_ds.attrs.update(rs_ds.attrs)

    out_ds = out_ds.assign_coords({c: ds[c] for c in ('diameter', 'min_diameter', 'max_diameter')
                                   if c in ds.coords})
    out_ds.attrs.update({
        'probe_name': PIPS_name,
        'deployment_name': deployment_name,
        'DSD_interval': DSD_interval,
        'source_file': os.path.basename(parsivel_combined_file),
        'QC_tags': ' '.join(args.QC_tags),
        'T_var': args.T_var,
        'p_var': args.p_var,
        'RH_var': args.RH_var,
        'rho': 'recomputed from T_var, p_var, RH_var (stored rho not used)',
        'fall_speed_binned': binned.attrs['fall_speed'],
        'ventilation': binned.attrs['ventilation'],
        'thermo_form': args.thermo_form,
        'saturation_definition': args.saturation_definition,
        'D_detect_min_mm': args.D_detect_min,
        'extension_fit': args.extension_fit,
        'extension_fallback': int(args.extension_fallback),
        'fits': ' '.join(args.fits),
        'closure_integration': args.closure_integration,
        'closure_fall_speed': closures.attrs['fall_speed'],
        'mu_fixed': str(args.mu_fixed),
        'N0_1mom': args.N0_1mom,
        'sign_convention': 'E > 0: liquid lost to vapor; cooling_rate < 0: cooling',
        'note': 'Surface evaporation and latent cooling rates from rain evaporation only; not a '
                'closed surface heat budget (advection and other terms not included). Closures '
                'use moments of the measured (truncated) spectrum.',
        'history': f"{datetime.now(timezone.utc).isoformat()}: {' '.join(sys.argv)}",
    })
    output_file = os.path.join(
        args.output_dir or PIPS_dir, f'evap_{deployment_name}_{PIPS_name}_{int(DSD_interval):d}s'
                  f'{"_" + args.output_tag if args.output_tag else ""}.nc')
    utils.log(f"Writing {output_file}")
    out_ds.to_netcdf(output_file)
