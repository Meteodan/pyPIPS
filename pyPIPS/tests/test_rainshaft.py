"""Tests for pyPIPS.rainshaft."""
from __future__ import annotations

import numpy as np
import pytest
import xarray as xr
from scipy.special import gammaln

import pyPIPS.DSDlib as dsd
import pyPIPS.evaporation as ev
import pyPIPS.rainshaft as rs


def gamma_start(mu=2., lamda=4000., Nt=2000., D_max_mm=8., dD_mm=0.05):
    N0 = Nt * np.exp((1. + mu) * np.log(lamda) - gammaln(1. + mu))
    edges = rs.fine_diameter_edges_mm(D_max_mm, dD_mm)
    return rs.starting_dsd_from_gamma(N0, lamda, mu, edges), edges, N0


@pytest.fixture(scope='module')
def subsaturated_profile():
    return rs.ThermoProfile.well_mixed_from_surface(298.15, 97000., 0.5, 3000.)


@pytest.fixture(scope='module')
def saturated_profile():
    z = np.linspace(0., 3000., 31)
    T = 293.15 - 0.0065 * z
    p = 97000. * np.exp(-z / 8000.)
    return rs.ThermoProfile.from_arrays(z, T, p, np.ones_like(z), source='test_saturated')


@pytest.fixture(scope='module')
def subsat_column(subsaturated_profile):
    ND, edges, _ = gamma_start()
    # Start well below the LCL (~1410 m): near saturation, the slowing of drops in denser air
    # can slightly increase concentrations (and Z) downward before evaporation dominates.
    return rs.rainshaft(ND, edges, 1200., 0., subsaturated_profile, dz_out=100., n_sub=10)


def test_well_mixed_profile(subsaturated_profile):
    prof = subsaturated_profile
    import pyPIPS.thermolib as thermo
    theta = thermo.caltheta(prof.p, prof.T)
    np.testing.assert_allclose(theta, theta[0], rtol=1.e-10)
    below = prof.z < prof.attrs['z_LCL']
    np.testing.assert_allclose(prof.w[below], prof.w[0], rtol=1.e-8)
    assert np.all(np.diff(prof.p) < 0.)
    assert np.all(prof.RH <= 1.)
    assert prof.source == 'well_mixed_from_PIPS_surface'
    assert 'assumption' in prof.attrs


# ------------------------------------------------------------------------------------------------
# 6. Rainshaft under subsaturation (consistent with Xie et al., 2016)
# ------------------------------------------------------------------------------------------------

def test_small_drops_vanish_first(subsat_column):
    ND = subsat_column['ND']
    D = subsat_column['diameter']
    small = (D > 0.25) & (D < 0.5)
    large = (D > 2.) & (D < 4.)
    top, sfc = ND.isel(height=0), ND.isel(height=-1)
    ratio_small = float(sfc.where(small).sum() / top.where(small).sum())
    ratio_large = float(sfc.where(large).sum() / top.where(large).sum())
    # Small drops are depleted more than large drops. (The smallest bins are not empty at the
    # surface in steady state: larger drops continually shrink through them before vanishing.)
    assert ratio_small < ratio_large


def test_profiles_evolve_as_expected(subsat_column):
    ds = subsat_column
    assert np.all(np.diff(ds['reflectivity'].values) < 0.)  # Z decreases downward
    assert np.all(np.diff(ds['rainrate'].values) < 0.)
    assert np.all(np.diff(ds['Dm43'].values) > 0.)  # mass-weighted size increases downward
    assert np.all(ds['cooling_rate'].values[1:] < 0.)
    assert 0. < float(ds['frac_rain_reaching'].isel(height=-1)) < 1.


def test_Dmv_increases_once_small_drops_evaporate_completely(subsaturated_profile):
    """D_mv = (M3/M0)^(1/3) is not monotonic in general: it decreases while all drops shrink
    and increases once enough small drops vanish. For a small-drop-rich DSD in a deep, dry
    column, D_mv at the surface exceeds D_mv at the top."""
    ND, edges, _ = gamma_start(mu=2., lamda=6000., Nt=3000.)
    ds = rs.rainshaft(ND, edges, 1400., 0., subsaturated_profile, dz_out=100., n_sub=10)
    assert float(ds['Dmv'].isel(height=-1)) > float(ds['Dmv'].isel(height=0))
    assert float(ds['Nt'].isel(height=-1)) < float(ds['Nt'].isel(height=0))


def test_mass_budget_closes(subsat_column):
    ds = subsat_column
    z = ds['height'].values
    # Loss of liquid mass flux between top and bottom equals the column-integrated evaporation
    flux_loss = (ds['rainrate'].values[0] - ds['rainrate'].values[-1]) * ev.rhow / 3.6e6
    E_int = -np.trapezoid(ds['E'].values, z)
    np.testing.assert_allclose(flux_loss, E_int, rtol=0.02)


def test_number_flux_conserved_without_evaporation(saturated_profile):
    ND, edges, _ = gamma_start()
    fs = ev.AtlasFallSpeed()
    out_edges = rs.fine_diameter_edges_mm(8., 0.05)
    ds = rs.rainshaft(ND, edges, 2000., 0., saturated_profile, dz_out=500., n_sub=1,
                      fall_speed=fs, out_edges_mm=out_edges)
    assert np.all(np.abs(ds['E'].values) < 1.e-20)
    np.testing.assert_allclose(ds['rainrate'], ds['rainrate'].isel(height=0), rtol=1.e-10)
    D = ds['diameter'] * 1.e-3
    dD = ds['max_diameter'] - ds['min_diameter']
    flux = (ds['ND'] * fs(D, ds['rho']) * dD).sum('diameter_bin')
    np.testing.assert_allclose(flux, flux.isel(height=0), rtol=1.e-10)


def test_feedback_reduces_cooling(subsaturated_profile):
    ND, edges, _ = gamma_start()
    base = rs.rainshaft(ND, edges, 1400., 0., subsaturated_profile, dz_out=200., n_sub=4)
    fb = rs.rainshaft(ND, edges, 1400., 0., subsaturated_profile, dz_out=200., n_sub=4,
                      feedback=True, feedback_time=1800., n_feedback_steps=4)
    assert fb.attrs['feedback'] == 1 and fb.attrs['feedback_time_s'] == 1800.
    assert float(fb['RH'].isel(height=-1)) > float(base['RH'].isel(height=-1))
    assert abs(float(fb['cooling_rate'].isel(height=-1))) < \
        abs(float(base['cooling_rate'].isel(height=-1)))
    with pytest.raises(ValueError):
        rs.rainshaft(ND, edges, 1400., 0., subsaturated_profile, feedback=True)


def test_drift_hook_and_metadata(subsaturated_profile, subsat_column):
    ND, edges, _ = gamma_start()
    with pytest.raises(NotImplementedError):
        rs.rainshaft(ND, edges, 1400., 0., subsaturated_profile, drift_correction='pips_winds')
    assert subsat_column.attrs['thermo_profile_source'] == 'well_mixed_from_PIPS_surface'
    assert subsat_column['ND_surface'].dims == ('diameter_bin',)
    assert subsat_column['ND_surface'].sizes['diameter_bin'] == 32  # Parsivel bins


# ------------------------------------------------------------------------------------------------
# Bulk rainshaft
# ------------------------------------------------------------------------------------------------

def test_moment_weighted_fallspeed_matches_calc_VQR_ferrier():
    rho, qr, Nt, mu = 1.1, 1.e-3, 1000., 1.5
    lamda = dsd.calc_lamda_gamma(rho, qr, Nt, ev.cmr, mu)
    fs = ev.FerrierFallSpeed(rho0=1.204, exponent=0.5)
    V3 = rs.moment_weighted_fallspeed(3, lamda, mu, fs, rho)
    np.testing.assert_allclose(V3, dsd.calc_VQR_ferrier(rho, qr, Nt, ev.cmr, mu), rtol=1.e-10)


@pytest.mark.parametrize('n_moments,kwargs', [(1, {}), (2, {'mu_fixed': 2.}),
                                              (2, {'mu_closure': 'diag_alpha'}), (3, {})])
def test_bulk_rainshaft_no_evaporation_conserves_fluxes(saturated_profile, n_moments, kwargs):
    ND, edges, _ = gamma_start()
    ds = rs.bulk_rainshaft(2000., 0., saturated_profile, ND_top=ND, D_edges_mm=edges,
                           n_moments=n_moments, dz=50., **kwargs)
    assert ds.attrs['converged_to_steady_state'] == 1
    np.testing.assert_allclose(ds['rainrate'], ds['rainrate'].isel(height=0), rtol=1.e-4)


@pytest.mark.parametrize('n_moments', [1, 2, 3])
def test_bulk_rainshaft_subsaturated(subsaturated_profile, n_moments):
    ND, edges, _ = gamma_start()
    ds = rs.bulk_rainshaft(1400., 0., subsaturated_profile, ND_top=ND, D_edges_mm=edges,
                           n_moments=n_moments, mu_fixed=2., dz=50.)
    assert ds.attrs['converged_to_steady_state'] == 1
    assert float(ds['rainrate'].isel(height=-1)) < float(ds['rainrate'].isel(height=0))
    assert float(ds['E'].isel(height=-1)) > 0.
    assert ds['ND_surface'].sizes['diameter_bin'] == 32


def test_starting_conditions_from_radar_at_PIPS():
    times = np.array(['2024-01-01T00:00'], dtype='datetime64[ns]')
    tag = 'Z01_qc'
    fields = ['REF', 'ZDR', f'N0_{tag}', f'lamda_{tag}', f'mu_{tag}']
    mu, lamda_mm, N0_mm = 2., 3., 5000.
    da = xr.DataArray(np.array([[40., 1.5, N0_mm, lamda_mm, mu]]), dims=['time', 'fields_KXXX'],
                      coords={'time': times, 'fields_KXXX': fields})
    ds = xr.Dataset({'KXXX_at_PIPS': da,
                     'KXXX_beam_height_at_PIPS': ('time', [850.])}, coords={'time': times})
    sc = rs.starting_conditions_from_radar_at_PIPS(ds, 'KXXX', tag)
    np.testing.assert_allclose(sc['lamda'], lamda_mm * 1000.)
    np.testing.assert_allclose(sc['N0'], N0_mm * 1000.**(1. + mu))
    np.testing.assert_allclose(sc['z_top'], 850.)
    # Same N(D) in either unit system at D = 1 mm
    D_mm = 1.
    N_mm = N0_mm * D_mm**mu * np.exp(-lamda_mm * D_mm)  # m^-3 mm^-1
    N_SI = sc['N0'].item() * (D_mm * 1.e-3)**mu * np.exp(-sc['lamda'].item() * D_mm * 1.e-3)
    np.testing.assert_allclose(N_SI, 1000. * N_mm, rtol=1.e-12)
