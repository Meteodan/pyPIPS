"""Tests for the surface time-integration tools in pyPIPS.evaporation."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from scipy.special import gammaln

import pyPIPS.DSDlib as dsd
import pyPIPS.evaporation as ev
import pyPIPS.parsivel_params as pp
import pyPIPS.thermolib as thermo

DT = 10.  # sample interval (s)


def make_series(values, n=None, name=None):
    values = np.atleast_1d(np.asarray(values, dtype=float))
    if n is not None and values.size == 1:
        values = np.full(n, values[0])
    times = pd.date_range('2024-05-01T00:00:10', periods=values.size, freq='10s')
    return xr.DataArray(values, dims=['time'], coords={'time': times.to_numpy()}, name=name)


def gamma_ND_series(Nt, lamda, mu, n, edges_mm=None, zero_times=()):
    edges = pp.parsivel_parameters['diameter_bin_edges_mm'] if edges_mm is None else edges_mm
    D_mm = 0.5 * (edges[:-1] + edges[1:])
    N0 = Nt * np.exp((1. + mu) * np.log(lamda) - gammaln(1. + mu))
    row = 1.e-3 * dsd.calc_binned_DSD_from_params(N0, lamda, mu, D_mm)
    vals = np.tile(row, (n, 1))
    for i in zero_times:
        vals[i] = 0.
    times = make_series(np.zeros(n))['time']
    return xr.DataArray(vals, dims=['time', 'diameter_bin'], coords={
        'time': times,
        'diameter': ('diameter_bin', D_mm),
        'min_diameter': ('diameter_bin', edges[:-1]),
        'max_diameter': ('diameter_bin', edges[1:]),
    }), N0


@pytest.fixture
def thermo_obs():
    n = 30
    return make_series(293.15, n), make_series(97000., n), make_series(0.6, n)


# ------------------------------------------------------------------------------------------------
# Static integration
# ------------------------------------------------------------------------------------------------

def test_static_constant_tendency_integrates_linearly(thermo_obs):
    T, p, RH = thermo_obs
    E = make_series(2.e-6, T.size)
    rho = ev.air_density(T, p, RH)
    out = ev.integrate_static(E, T, p, RH, rho=rho, sample_interval=DT)
    k = np.arange(1, T.size + 1)
    dTdt = float(-thermo.Lv * 2.e-6 / (rho[0] * thermo.cp))
    np.testing.assert_allclose(out['dT'], dTdt * DT * k, rtol=1.e-12)
    w0 = float(ev.mixing_ratio_state(T[0], p[0], RH[0]))
    rho_d = float(rho[0]) * (1. - w0 / (1. + w0))
    np.testing.assert_allclose(out['dw'], 2.e-6 / rho_d * DT * k, rtol=1.e-12)
    # Energy per unit volume: rho c_p dT + rho_d L_v dw = 0
    np.testing.assert_allclose(float(rho[0]) * thermo.cp * out['dT'] +
                               rho_d * thermo.Lv * out['dw'], 0., atol=1.e-9)
    assert int(out['n_rain']) == T.size and int(out['n_missing']) == 0
    assert out.attrs['integration_mode'] == 'static'


def test_static_no_rain_and_missing_are_counted_not_integrated(thermo_obs):
    T, p, RH = thermo_obs
    vals = np.full(T.size, 1.e-6)
    vals[5:10] = 0.  # no rain
    vals[12:15] = np.nan  # missing
    E = make_series(vals)
    missing = make_series(np.zeros(T.size)).astype(bool)
    missing[20] = True  # flagged missing record
    out = ev.integrate_static(E, T, p, RH, sample_interval=DT, missing=missing)
    assert int(out['n_no_rain']) == 5
    assert int(out['n_missing']) == 4
    assert int(out['n_rain']) == T.size - 9
    assert float(out['duration_missing_s']) == 4 * DT
    dT = out['dT'].to_numpy()
    for i in list(range(5, 10)) + list(range(12, 15)) + [20]:
        assert dT[i] == dT[i - 1]  # zero tendency during no-rain and missing samples


def test_static_zero_tendency_gives_zero_change(thermo_obs):
    T, p, RH = thermo_obs
    out = ev.integrate_static(make_series(0., T.size), T, p, RH, sample_interval=DT)
    assert np.all(out['dT'] == 0.) and np.all(out['dw'] == 0.)
    np.testing.assert_allclose(out['T'], float(T[0]))
    np.testing.assert_allclose(out['RH'], 0.6, rtol=1.e-10)


def test_static_time_window_and_initial_state():
    n = 30
    T = make_series(np.linspace(290., 293., n))
    p = make_series(97000., n)
    RH = make_series(0.6, n)
    T[10] = np.nan  # first time in the window has missing T -> next valid time sets T0
    E = make_series(1.e-6, n)
    start, end = T['time'].values[10], T['time'].values[19]
    out = ev.integrate_static(E, T, p, RH, start=start, end=end, sample_interval=DT)
    assert out.sizes['time'] == 10
    np.testing.assert_allclose(out.attrs['T0_K'], float(T[11]))
    assert str(np.datetime64(out.attrs['integration_start'])) == str(start)


def test_static_with_closure_dimension(thermo_obs):
    T, p, RH = thermo_obs
    E = xr.concat([make_series(1.e-6, T.size), make_series(2.e-6, T.size)],
                  dim=xr.DataArray(['a', 'b'], dims=['closure']))
    out = ev.integrate_static(E, T, p, RH, sample_interval=DT)
    assert set(out['dT'].dims) == {'closure', 'time'}
    np.testing.assert_allclose(out['dT'].sel(closure='b'), 2. * out['dT'].sel(closure='a'),
                               rtol=1.e-12)


# ------------------------------------------------------------------------------------------------
# Dynamic integration
# ------------------------------------------------------------------------------------------------

def test_dynamic_matches_static_for_weak_evaporation(thermo_obs):
    T, p, RH = thermo_obs
    ND, _ = gamma_ND_series(50., 4000., 2., T.size)
    binned = ev.evaporation_rate_binned(ND, T, p, RH)
    static = ev.integrate_static(binned['E'], T, p, RH, sample_interval=DT)
    dynamic = ev.integrate_dynamic(ev.binned_rate_function(ND), T, p, RH, sample_interval=DT)
    assert float(dynamic['dT'][-1]) < 0.
    np.testing.assert_allclose(dynamic['dT'], static['dT'], rtol=1.e-2)
    np.testing.assert_allclose(dynamic['E'][0], binned['E'][0], rtol=1.e-10)


def test_dynamic_evaporation_is_self_limiting_and_approaches_wet_bulb():
    n = 1080  # 3 h of heavy rain at 10 s
    T, p, RH = make_series(298.15, n), make_series(97000., n), make_series(0.5, n)
    ND, _ = gamma_ND_series(3000., 2500., 2., n)
    binned = ev.evaporation_rate_binned(ND, T, p, RH)
    static = ev.integrate_static(binned['E'], T, p, RH, sample_interval=DT)
    dynamic = ev.integrate_dynamic(ev.binned_rate_function(ND), T, p, RH, sample_interval=DT)
    assert abs(float(dynamic['dT'][-1])) < abs(float(static['dT'][-1]))
    assert np.all(np.diff(dynamic['E'].values) <= 0.)  # evaporation slows as air moistens
    assert float(dynamic['RH'][-1]) > 0.99
    assert float(dynamic['RH'].max()) < 1.001
    # The end state is close to the wet-bulb temperature of the initial state
    mpcalc = pytest.importorskip('metpy.calc')
    from metpy.units import units
    Td = thermo.calTdfromRH(97000., 298.15, 0.5)
    Tw = mpcalc.wet_bulb_temperature(970. * units.hPa, 298.15 * units.K,
                                     Td * units.K).to('K').magnitude
    assert abs(float(dynamic['T'][-1]) - float(Tw)) < 0.5


def test_dynamic_no_rain_and_missing(thermo_obs):
    T, p, RH = thermo_obs
    ND, _ = gamma_ND_series(500., 3000., 1., T.size, zero_times=range(5, 10))
    ND[12:14] = np.nan  # missing DSD records
    missing = make_series(np.zeros(T.size)).astype(bool)
    missing[20] = True
    out = ev.integrate_dynamic(ev.binned_rate_function(ND), T, p, RH, sample_interval=DT,
                               missing=missing)
    assert int(out['n_no_rain']) == 5
    assert int(out['n_missing']) == 3
    dT = out['dT'].to_numpy()
    for i in list(range(5, 10)) + [12, 13, 20]:
        assert dT[i] == dT[i - 1]


def test_dynamic_extension_treats_no_drops_as_no_rain(thermo_obs):
    T, p, RH = thermo_obs
    ND, N0 = gamma_ND_series(500., 3000., 1., T.size, zero_times=range(5, 10))
    ND[:, :2] = 0.  # Parsivel does not report the two smallest bins
    nanfill = make_series(np.where(np.isin(np.arange(T.size), range(5, 10)), np.nan, 1.))
    ext = (N0 * nanfill, 3000. * nanfill, 1. * nanfill)  # fits undefined when no drops
    out = ev.integrate_dynamic(ev.binned_rate_function(ND, extension=ext), T, p, RH,
                               sample_interval=DT)
    assert int(out['n_no_rain']) == 5 and int(out['n_missing']) == 0
    base = ev.integrate_dynamic(ev.binned_rate_function(ND), T, p, RH, sample_interval=DT)
    assert float(out['dT'][-1]) < float(base['dT'][-1])  # extension adds small-drop cooling


def test_dynamic_gamma_members_and_consistency_with_binned(thermo_obs):
    T, p, RH = thermo_obs
    lamda, mu = 3000., 1.
    edges = np.arange(0., 12.0001, 0.02)
    ND, N0 = gamma_ND_series(800., lamda, mu, T.size, edges_mm=edges)
    closure = xr.DataArray(['truth', 'double'], dims=['closure'])
    N0s = xr.concat([N0 * xr.ones_like(T), 2. * N0 * xr.ones_like(T)], dim=closure)
    rate = ev.gamma_rate_function(N0s, lamda * xr.ones_like(N0s), mu * xr.ones_like(N0s))
    assert rate.members == ('closure', ['truth', 'double'])
    out = ev.integrate_dynamic(rate, T, p, RH, sample_interval=DT)
    assert set(out['dT'].dims) == {'closure', 'time'}
    binned = ev.integrate_dynamic(ev.binned_rate_function(ND), T, p, RH, sample_interval=DT)
    np.testing.assert_allclose(out['dT'].sel(closure='truth'), binned['dT'], rtol=2.e-3)
    assert float(out['dT'].sel(closure='double')[-1]) < float(out['dT'].sel(closure='truth')[-1])


def test_gamma_rate_function_no_rain_mask(thermo_obs):
    T, p, RH = thermo_obs
    nanfill = make_series(np.where(np.arange(T.size) < 5, np.nan, 1.))
    rate = ev.gamma_rate_function(1.e7 * nanfill, 2000. * nanfill, 0. * nanfill,
                                  no_rain=(np.arange(T.size) < 5))
    out = ev.integrate_dynamic(rate, T, p, RH, sample_interval=DT)
    assert int(out['n_no_rain']) == 5 and int(out['n_missing']) == 0


def test_nonfinite_rates_are_counted_as_missing(thermo_obs):
    T, p, RH = thermo_obs
    vals = np.full(T.size, 1.e-6)
    vals[3] = np.inf
    out = ev.integrate_static(make_series(vals), T, p, RH, sample_interval=DT)
    assert int(out['n_missing']) == 1 and np.all(np.isfinite(out['dT']))
    params = make_series(np.where(np.arange(T.size) == 3, np.inf, 1.e7))
    rate = ev.gamma_rate_function(params, 2000. * xr.ones_like(params), 0. * params)
    out = ev.integrate_dynamic(rate, T, p, RH, sample_interval=DT)
    assert int(out['n_missing']) == 1 and np.all(np.isfinite(out['dT']))


# ------------------------------------------------------------------------------------------------
# Extension fallback and last-good filling of missing samples
# ------------------------------------------------------------------------------------------------

def test_extension_fallback(thermo_obs):
    T, p, RH = thermo_obs
    ND, N0 = gamma_ND_series(500., 3000., 1., T.size)
    ND[:, :2] = 0.
    mu_vals = np.where(np.arange(T.size) == 4, -3., 1.)  # invalid extension at sample 4
    ext = (N0 * xr.ones_like(T), 3000. * xr.ones_like(T), make_series(mu_vals))
    off = ev.evaporation_rate_binned(ND, T, p, RH, extension=ext)
    on = ev.evaporation_rate_binned(ND, T, p, RH, extension=ext, extension_fallback=True)
    assert np.isnan(float(off['E_ext'][4])) and 'extension_fallback_used' not in off
    np.testing.assert_allclose(on['E_ext'][4], on['E'][4], rtol=1.e-12)
    assert int(on['extension_fallback_used'].sum()) == 1
    np.testing.assert_allclose(on['E_ext'].drop_isel(time=4), off['E_ext'].drop_isel(time=4))
    # The dynamic rate function applies the same fallback
    dyn = ev.integrate_dynamic(ev.binned_rate_function(ND, extension=ext, extension_fallback=True),
                               T, p, RH, sample_interval=DT)
    assert int(dyn['n_missing']) == 0
    np.testing.assert_allclose(dyn['E'][0], on['E_ext'][0], rtol=1.e-10)


def test_static_last_good_fill(thermo_obs):
    T, p, RH = thermo_obs
    vals = np.full(T.size, 1.e-6)
    gaps = [0, 1, 10, 11, 12]
    vals[gaps] = np.nan  # samples 0-1 have no earlier good sample
    E = make_series(vals)
    zero = ev.integrate_static(E, T, p, RH, sample_interval=DT)
    full = ev.integrate_static(make_series(1.e-6, T.size), T, p, RH, sample_interval=DT)
    fill = ev.integrate_static(E, T, p, RH, sample_interval=DT, missing_fill='last_good')
    assert int(fill['n_missing']) == 5 and int(fill['n_filled']) == 3
    assert int(zero['n_filled']) == 0
    # Constant E: filling restores everything except the two leading samples
    np.testing.assert_allclose(fill['dT'][-1], full['dT'][-1] * (T.size - 2) / T.size,
                               rtol=1.e-12)
    limited = ev.integrate_static(E, T, p, RH, sample_interval=DT, missing_fill='last_good',
                                  max_fill_samples=2)
    assert int(limited['n_filled']) == 2
    with pytest.raises(ValueError):
        ev.integrate_static(E, T, p, RH, sample_interval=DT, missing_fill='bogus')


def test_dynamic_last_good_fill_reproduces_gap_free_run_for_constant_dsd(thermo_obs):
    T, p, RH = thermo_obs
    ND_full, _ = gamma_ND_series(500., 3000., 1., T.size)
    ND_gaps = ND_full.copy()
    ND_gaps[[8, 9, 15]] = np.nan
    missing = make_series(np.zeros(T.size)).astype(bool)
    missing[20] = True
    full = ev.integrate_dynamic(ev.binned_rate_function(ND_full), T, p, RH, sample_interval=DT)
    fill = ev.integrate_dynamic(ev.binned_rate_function(ND_gaps), T, p, RH, sample_interval=DT,
                                missing=missing, missing_fill='last_good')
    assert int(fill['n_missing']) == 4 and int(fill['n_filled']) == 4
    np.testing.assert_allclose(fill['dT'], full['dT'], rtol=1.e-12)
    zero = ev.integrate_dynamic(ev.binned_rate_function(ND_gaps), T, p, RH, sample_interval=DT,
                                missing=missing)
    assert float(zero['dT'][-1]) > float(fill['dT'][-1])  # less cooling without filling
    limited = ev.integrate_dynamic(ev.binned_rate_function(ND_gaps), T, p, RH,
                                   sample_interval=DT, missing_fill='last_good',
                                   max_fill_samples=1)
    assert int(limited['n_filled']) == 2  # samples 8 and 15 (9 exceeds the limit)
