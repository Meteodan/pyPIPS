"""Tests for pyPIPS.evaporation."""
from __future__ import annotations

import warnings

import numpy as np
import pytest
import xarray as xr
from scipy.integrate import quad
from scipy.special import gammaln

import pyPIPS.DSDlib as dsd
import pyPIPS.evaporation as ev
import pyPIPS.parsivel_params as pp

FALL_SPEEDS = {
    'ferrier': ev.FerrierFallSpeed(),
    'powerlaw': ev.PowerLawFallSpeed(a=842., b=0.8),
}
THERMO_STATES = [(275., 90000., 0.4), (293., 97000., 0.7), (303., 100000., 0.95)]


def gamma_N0(Nt, lamda, mu):
    """Intercept (m^-(4+mu)) for total concentration Nt (m^-3)."""
    return Nt * np.exp((1. + mu) * np.log(lamda) - gammaln(1. + mu))


def make_ND(edges_mm, values, times=None):
    """pyPIPS-style ND DataArray (m^-3 mm^-1) on the given bin edges."""
    edges = np.asarray(edges_mm, dtype=float)
    coords = {
        'diameter': ('diameter_bin', 0.5 * (edges[:-1] + edges[1:])),
        'min_diameter': ('diameter_bin', edges[:-1]),
        'max_diameter': ('diameter_bin', edges[1:]),
    }
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return xr.DataArray(values, dims=['diameter_bin'], coords=coords)
    coords['time'] = times
    return xr.DataArray(values, dims=['time', 'diameter_bin'], coords=coords)


def gamma_ND(edges_mm, N0, lamda, mu):
    edges = np.asarray(edges_mm, dtype=float)
    D_mm = 0.5 * (edges[:-1] + edges[1:])
    return make_ND(edges, 1.e-3 * dsd.calc_binned_DSD_from_params(N0, lamda, mu, D_mm))


# ------------------------------------------------------------------------------------------------
# 1. Gamma closed form vs numerical quadrature
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('fs_name', list(FALL_SPEEDS))
@pytest.mark.parametrize('mu', [-0.5, 0., 2., 6., 12.])
@pytest.mark.parametrize('lamda', [1000., 3000., 9000.])
@pytest.mark.parametrize('state', THERMO_STATES)
def test_gamma_closed_form_matches_quadrature(fs_name, mu, lamda, state):
    fs = FALL_SPEEDS[fs_name]
    T, p, RH = state
    N0 = gamma_N0(1000., lamda, mu)
    env = ev.evaporation_environment(T, p, RH)
    E_closed = ev.evaporation_rate_gamma(N0, lamda, mu, env=env, fall_speed=fs)

    def integrand(D):
        return N0 * D**mu * np.exp(-lamda * D) * -ev.drop_mass_rate(D, fall_speed=fs, env=env)

    D_hi = (mu + 80.) / lamda
    E_quad = quad(integrand, 0., D_hi, limit=1000, epsabs=0., epsrel=1.e-11,
                  points=[0.1 / lamda, 1. / lamda, (mu + 1.) / lamda])[0]
    np.testing.assert_allclose(E_closed, E_quad, rtol=1.e-6)

    E_gl = ev.evaporation_rate_gamma_quadrature(N0, lamda, mu, env=env, fall_speed=fs,
                                                D_hi=D_hi)
    np.testing.assert_allclose(E_closed, E_gl, rtol=1.e-6)


@pytest.mark.parametrize('fs_name', list(FALL_SPEEDS))
def test_truncated_gamma_closed_form_matches_quadrature(fs_name):
    fs = FALL_SPEEDS[fs_name]
    mu, lamda = 1.5, 2500.
    N0 = gamma_N0(800., lamda, mu)
    env = ev.evaporation_environment(290., 95000., 0.6)
    D_lo, D_hi = 0.3e-3, 4.e-3
    E_closed = ev.evaporation_rate_gamma(N0, lamda, mu, env=env, fall_speed=fs, D_lo=D_lo,
                                         D_hi=D_hi)
    E_gl = ev.evaporation_rate_gamma_quadrature(N0, lamda, mu, env=env, fall_speed=fs, D_lo=D_lo,
                                                D_hi=D_hi)
    np.testing.assert_allclose(E_closed, E_gl, rtol=1.e-8)


def test_third_moment_tendency_consistent_with_evaporation_rate():
    mu, lamda = 2., 3000.
    N0 = gamma_N0(1000., lamda, mu)
    env = ev.evaporation_environment(290., 95000., 0.6)
    E = ev.evaporation_rate_gamma(N0, lamda, mu, env=env)
    dM3dt = ev.moment_tendency_gamma(3, N0, lamda, mu, env=env)
    np.testing.assert_allclose(-ev.cmr * dM3dt, E, rtol=1.e-12)
    with pytest.raises(ValueError):
        ev.moment_tendency_gamma(0, N0, lamda, mu, env=env)


def test_closed_form_rejects_fall_speed_without_closed_form():
    env = ev.evaporation_environment(290., 95000., 0.6)
    with pytest.raises(TypeError):
        ev.evaporation_rate_gamma(1.e7, 2000., 0., env=env, fall_speed=ev.AtlasFallSpeed())


# ------------------------------------------------------------------------------------------------
# 3. Binned sum converges to the gamma closed form
# ------------------------------------------------------------------------------------------------

def test_binned_converges_to_gamma_closed_form():
    fs = FALL_SPEEDS['powerlaw']
    mu, lamda = 2., 2500.
    N0 = gamma_N0(1000., lamda, mu)
    T, p, RH = 293., 97000., 0.6
    D_max_mm = 12.
    E_exact = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, RH, fall_speed=fs, D_lo=0.,
                                        D_hi=D_max_mm * 1.e-3)
    errors = []
    for width_mm in (0.4, 0.2, 0.1, 0.05):
        edges = np.arange(0., D_max_mm + 0.5 * width_mm, width_mm)
        ND = gamma_ND(edges, N0, lamda, mu)
        E_bin = ev.evaporation_rate_binned(ND, T, p, RH, fall_speed=fs)['E']
        errors.append(abs(float(E_bin) / E_exact - 1.))
    assert all(e1 < e0 for e0, e1 in zip(errors[:-1], errors[1:]))
    assert errors[-1] < 1.e-3


def test_binned_vs_gamma_on_parsivel_bins_reported():
    """On Parsivel bins the difference is reported, not asserted."""
    fs = FALL_SPEEDS['powerlaw']
    mu, lamda = 2., 2500.
    N0 = gamma_N0(1000., lamda, mu)
    T, p, RH = 293., 97000., 0.6
    edges = pp.parsivel_parameters['diameter_bin_edges_mm']
    ND = gamma_ND(edges, N0, lamda, mu)
    E_bin = float(ev.evaporation_rate_binned(ND, T, p, RH, fall_speed=fs)['E'])
    E_exact = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, RH, fall_speed=fs)
    print(f"Parsivel bins: binned/gamma - 1 = {E_bin / E_exact - 1.:.3e}")  # noqa: T201
    assert np.isfinite(E_bin)


# ------------------------------------------------------------------------------------------------
# 4. Saturation behavior and monotonicity in RH
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('saturation_definition', ['vapor_pressure', 'mixing_ratio'])
def test_zero_at_saturation_negative_when_supersaturated(saturation_definition):
    mu, lamda = 1., 2000.
    N0 = gamma_N0(1000., lamda, mu)
    T, p = 290., 95000.
    E_sat = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, 1.0,
                                      saturation_definition=saturation_definition)
    assert abs(E_sat) < 1.e-20
    E_super = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, 1.03,
                                        saturation_definition=saturation_definition)
    assert E_super < 0.

    ND = gamma_ND(pp.parsivel_parameters['diameter_bin_edges_mm'], N0, lamda, mu)
    out = ev.evaporation_rate_binned(ND, T, p, 1.0, saturation_definition=saturation_definition)
    assert abs(float(out['E'])) < 1.e-20
    out = ev.evaporation_rate_binned(ND, T, p, 1.03, saturation_definition=saturation_definition)
    assert float(out['E']) < 0.
    assert float(out['cooling_rate']) > 0.  # condensation warms


def test_evaporation_increases_monotonically_as_RH_decreases():
    mu, lamda = 1., 2000.
    N0 = gamma_N0(1000., lamda, mu)
    RH = np.linspace(0.99, 0.2, 25)
    E = ev.evaporation_rate_gamma(N0, lamda, mu, 293., 97000., RH)
    assert np.all(np.diff(E) > 0.)
    ND = gamma_ND(pp.parsivel_parameters['diameter_bin_edges_mm'], N0, lamda, mu)
    E_bin = [float(ev.evaporation_rate_binned(ND, 293., 97000., rh)['E']) for rh in RH]
    assert np.all(np.diff(E_bin) > 0.)


# ------------------------------------------------------------------------------------------------
# 5. At fixed q_r, more small drops evaporate faster
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('mu', [0., 2.5])
def test_more_small_drops_evaporate_faster_at_fixed_qr(mu):
    rho, qr = 1.15, 1.e-3
    T, p, RH = 293., 97000., 0.6
    E = []
    for Nt in (100., 1000., 10000.):
        lamda = dsd.calc_lamda_gamma(rho, qr, Nt, ev.cmr, mu)
        N0 = dsd.calc_N0_gamma(rho, qr, Nt, ev.cmr, mu)
        np.testing.assert_allclose(dsd.calc_qr_gamma(rho, N0, lamda, mu), qr, rtol=1.e-10)
        E.append(ev.evaporation_rate_gamma(N0, lamda, mu, T, p, RH, rho=rho))
    assert E[0] < E[1] < E[2]


# ------------------------------------------------------------------------------------------------
# Single-drop relations
# ------------------------------------------------------------------------------------------------

def test_single_drop_signs_and_consistency():
    env = ev.evaporation_environment(293., 97000., 0.6)
    D = np.array([0.1e-3, 0.5e-3, 1.e-3, 3.e-3])
    dmdt = ev.drop_mass_rate(D, env=env)
    dDdt = ev.drop_diameter_rate(D, env=env)
    dD2dt = ev.drop_diameter_squared_rate(D, env=env)
    dDdz = ev.drop_diameter_height_rate(D, env=env)
    assert np.all(dmdt < 0.) and np.all(dDdt < 0.)
    assert np.all(dDdz > 0.)  # drop shrinks as it descends (z positive upward)
    np.testing.assert_allclose(dmdt, ev.rhow * np.pi / 2. * D**2. * dDdt, rtol=1.e-12)
    np.testing.assert_allclose(dD2dt, 2. * D * dDdt, rtol=1.e-12)
    np.testing.assert_allclose(dDdz, -dDdt / ev.AtlasFallSpeed()(D, env['rho']), rtol=1.e-12)
    assert np.isfinite(ev.drop_diameter_squared_rate(np.array([0.]), env=env)).all()


def test_srivastava_coen_hook_not_implemented():
    env = ev.evaporation_environment(293., 97000., 0.6)
    with pytest.raises(NotImplementedError):
        ev.drop_mass_rate(1.e-3, env=env, method='srivastava_coen1992')


def test_atlas_fall_speed_matches_existing_function():
    import pyPIPS.PIPS as pips
    D_mm = pp.parsivel_parameters['avg_diameter_bins_mm']
    rho = 1.1
    v_existing = pips.calc_empirical_fallspeed_rain(D_mm, correct_rho=True, rho=rho)
    v_new = ev.AtlasFallSpeed()(D_mm * 1.e-3, rho)
    np.testing.assert_allclose(v_new, v_existing, rtol=1.e-12)


def test_thermo_forms_differ_by_documented_amount():
    T, p, RH = 293.15, 95000., 0.5
    mu, lamda = 1., 2000.
    N0 = gamma_N0(1000., lamda, mu)
    E_new = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, RH)
    E_old = ev.evaporation_rate_gamma(N0, lamda, mu, T, p, RH, thermo_form='milbrandt',
                                      saturation_definition='mixing_ratio')
    assert 1. < E_new / E_old < 1.03


# ------------------------------------------------------------------------------------------------
# Small-drop diagnostics and sub-detectable extension
# ------------------------------------------------------------------------------------------------

def test_cumulative_fraction_and_extension():
    mu, lamda = 1., 3000.
    N0 = gamma_N0(2000., lamda, mu)
    T, p, RH = 293., 97000., 0.6
    edges = pp.parsivel_parameters['diameter_bin_edges_mm']
    ND = gamma_ND(edges, N0, lamda, mu)
    ND[:2] = 0.  # Parsivel does not report the two smallest bins
    out = ev.evaporation_rate_binned(ND, T, p, RH, extension=(N0, lamda, mu))
    cum = out['E_cumfrac'].to_numpy()
    assert np.all(np.diff(cum) >= -1.e-12)
    np.testing.assert_allclose(cum[-1], 1., rtol=1.e-12)
    assert float(out['E_ext']) > float(out['E'])
    assert 0. < float(out['frac_E_sub']) < 1.
    np.testing.assert_allclose(out.attrs['D_cut_m'], 0.25e-3)
    E_true = ev.evaporation_rate_gamma_quadrature(N0, lamda, mu, T, p, RH)
    # The extension recovers the missing small-drop contribution (to bin-discretization error)
    assert abs(float(out['E_ext']) / E_true - 1.) < abs(float(out['E']) / E_true - 1.)


# ------------------------------------------------------------------------------------------------
# Closures and number tendencies
# ------------------------------------------------------------------------------------------------

def test_closure_comparison_recovers_gamma_with_three_moments():
    mu, lamda = 2., 2500.
    N0 = gamma_N0(1000., lamda, mu)
    times = np.array(['2024-01-01T00:00', '2024-01-01T00:01'], dtype='datetime64[ns]')
    edges = np.arange(0., 12.0001, 0.02)
    D_mm = 0.5 * (edges[:-1] + edges[1:])
    vals = 1.e-3 * dsd.calc_binned_DSD_from_params(N0, lamda, mu, D_mm)
    ND = make_ND(edges, np.vstack([vals, 2. * vals]), times=times)
    T = xr.DataArray([293., 293.], dims=['time'], coords={'time': times})
    p = 97000. * xr.ones_like(T)
    RH = 0.6 * xr.ones_like(T)
    out = ev.closure_comparison(ND, T, p, RH, mu_fixed=(0., 2.), fits={'truth': (N0, lamda, mu)})
    assert set(out.dims) == {'closure', 'time'}
    ratios = out['ratio_to_binned']
    # With the true shape parameter, the 2M closure reproduces the gamma from (q_r, N_t)
    np.testing.assert_allclose(out['mu'].sel(closure='MM03_mu2'), 2.)
    np.testing.assert_allclose(ratios.sel(closure='MM03_mu2'), 1., rtol=1.e-3)
    np.testing.assert_allclose(out['mu'].sel(closure='MM036'), 2., rtol=1.e-3)
    np.testing.assert_allclose(ratios.sel(closure='MM036'), 1., rtol=1.e-3)
    np.testing.assert_allclose(ratios.sel(closure='MM036_solve_alpha'), 1., rtol=0.05)
    np.testing.assert_allclose(ratios.sel(closure='truth').isel(time=0), 1., rtol=1.e-3)
    assert np.all(np.isfinite(ratios.sel(closure='MM3_N0fixed_mu0')))
    assert 'closure_description' in out.coords
    assert 'M0, M3, and M6' in str(out['closure_description'].sel(closure='MM036').item())


def test_number_tendency_methods():
    rho = 1.1
    q = np.array([1.e-3, 1.e-9])
    N = np.array([1.e3, 1.e3])
    dq_dt = np.array([-1.e-6, -1.e-12])
    prop = ev.number_tendency(q, N, dq_dt)  # default: 'proportional'
    np.testing.assert_allclose(prop, N * dq_dt / q)
    # Constant mean volume diameter: q/N is unchanged to first order over a small step
    dt = 1.
    np.testing.assert_allclose((q + dq_dt * dt) / (N + prop * dt), q / N, rtol=1.e-12)
    assert ev.number_tendency(np.array([0.]), np.array([10.]), np.array([0.]))[0] == 0.
    with pytest.raises(NotImplementedError):
        ev.number_tendency(q, N, dq_dt, rho, method='seifert2008')
    for bad in ('bogus', 'fixed_until_small'):
        with pytest.raises(ValueError):
            ev.number_tendency(q, N, dq_dt, rho, method=bad)


# ------------------------------------------------------------------------------------------------
# 7. xarray round trips and unit conversions
# ------------------------------------------------------------------------------------------------

def make_parsivel_like_ds():
    times = np.array(['2024-01-01T00:00', '2024-01-01T00:00:10', '2024-01-01T00:00:20'],
                     dtype='datetime64[ns]')
    edges = pp.parsivel_parameters['diameter_bin_edges_mm']
    mu, lamda = 1., 2500.
    D_mm = 0.5 * (edges[:-1] + edges[1:])
    vals = 1.e-3 * dsd.calc_binned_DSD_from_params(gamma_N0(1000., lamda, mu), lamda, mu, D_mm)
    ND = make_ND(edges, np.vstack([vals, 0. * vals, np.full_like(vals, np.nan)]), times=times)
    ds = xr.Dataset({'ND_qc': ND})
    ds['fasttemp'] = ('time', [20., 21., 22.], {'units': 'degrees Celsius'})
    ds['pressure'] = ('time', [970., 969., 968.], {'units': 'hectoPascals'})
    ds['RH_derived'] = ('time', [60., 65., 70.], {'units': 'percent'})
    return ds


def test_thermo_from_parsivel_ds_converts_units():
    ds = make_parsivel_like_ds()
    th = ev.thermo_from_parsivel_ds(ds)
    np.testing.assert_allclose(th['T'], ds['fasttemp'] + 273.15)
    np.testing.assert_allclose(th['p'], ds['pressure'] * 100.)
    np.testing.assert_allclose(th['RH'], ds['RH_derived'] / 100.)
    assert th.attrs['T_var'] == 'fasttemp'
    ds['fasttemp'].attrs = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        th2 = ev.thermo_from_parsivel_ds(ds)
    assert any('assuming degrees Celsius' in str(x.message) for x in w)
    np.testing.assert_allclose(th2['T'], th['T'])


def test_binned_xarray_roundtrip(tmp_path):
    ds = make_parsivel_like_ds()
    th = ev.thermo_from_parsivel_ds(ds)
    out = ev.evaporation_rate_binned(ds['ND_qc'], th['T'], th['p'], th['RH'], rho=th['rho'])
    assert out['E'].dims == ('time',)
    assert out['E_bin'].dims == ('time', 'diameter_bin')
    assert out['E'].isel(time=0) > 0.
    assert float(out['E'].isel(time=1)) == 0.  # no drops
    assert np.isnan(float(out['E'].isel(time=2)))  # missing DSD
    assert out['E'].attrs['units'] == 'kg m-3 s-1'
    np.testing.assert_allclose(out['cooling_rate_K_per_h'], 3600. * out['cooling_rate'])
    fname = tmp_path / 'evap.nc'
    out.to_netcdf(fname)
    back = xr.load_dataset(fname)
    xr.testing.assert_allclose(back['E'], out['E'])
    np.testing.assert_allclose(back['diameter'], ds['diameter'])


def test_gamma_params_from_fit():
    times = np.array(['2024-01-01T00:00'], dtype='datetime64[ns]')
    fit = xr.DataArray(np.array([[1.e7], [2000.], [1.5]]), dims=['parameter', 'time'],
                       coords={'parameter': ['N0', 'lamda', 'alpha'], 'time': times})
    N0, lamda, mu = ev.gamma_params_from_fit(fit)
    assert N0.item() == 1.e7 and lamda.item() == 2000. and mu.item() == 1.5
    assert 'parameter' not in N0.coords


class _AtlasNoClosedForm(ev.FallSpeed):
    """Atlas fall speed without the small-drop closed-form shortcut (pure quadrature)."""

    def __init__(self):
        self._atlas = ev.AtlasFallSpeed()
        self.rho0, self.exponent = self._atlas.rho0, self._atlas.exponent
        self.breakpoints = self._atlas.breakpoints

    def __call__(self, D, rho=None):
        return self._atlas(D, rho)


@pytest.mark.parametrize('mu', [0., 2., 6.])
def test_atlas_small_drop_closed_form_matches_pure_quadrature(mu):
    lamda = 2500.
    N0 = gamma_N0(1000., lamda, mu)
    env = ev.evaporation_environment(293., 97000., 0.6)
    E_hybrid = ev.evaporation_rate_gamma_quadrature(N0, lamda, mu, env=env)
    E_gl = ev.evaporation_rate_gamma_quadrature(N0, lamda, mu, env=env,
                                                fall_speed=_AtlasNoClosedForm(),
                                                max_width=0.05e-3)
    np.testing.assert_allclose(E_hybrid, E_gl, rtol=1.e-6)


def test_extension_divergent_for_mu_le_minus_two():
    lamda = 2500.
    env = ev.evaporation_environment(293., 97000., 0.6)
    assert np.isnan(ev.evaporation_rate_gamma_quadrature(1.e3, lamda, -2.5, env=env, D_hi=0.25e-3))
    E = ev.evaporation_rate_gamma_quadrature(1.e3, lamda, -1.8, env=env, D_hi=0.25e-3)
    assert np.isfinite(E) and E > 0.


def gamma_moment(k, N0, lamda, mu):
    return np.exp(np.log(N0) + gammaln(mu + k + 1.) - (mu + k + 1.) * np.log(lamda))


@pytest.mark.parametrize('mu', [-0.5, 0., 2., 10., 50.])
def test_fit_DSD_MM036_is_exact(mu):
    lamda, Nt = 3000., 500.
    N0 = gamma_N0(Nt, lamda, mu)
    M0, M3, M6 = (gamma_moment(k, N0, lamda, mu) for k in (0, 3, 6))
    N0_f, lamda_f, mu_f = dsd.fit_DSD_MM036(np.array([M0]), np.array([M3]), np.array([M6]))
    np.testing.assert_allclose(mu_f, mu, rtol=1.e-9, atol=1.e-9)
    np.testing.assert_allclose(lamda_f, lamda, rtol=1.e-9)
    np.testing.assert_allclose(N0_f, N0, rtol=1.e-8)


def test_closures_are_method_of_moments_fits():
    """Each closure reproduces exactly the observed moments it is constrained by (except the
    approximate solve_alpha variant, which matches M0 and M3 exactly and M6 approximately)."""
    rng = np.random.default_rng(1)
    n = 200
    rho = xr.DataArray(rng.uniform(1.0, 1.25, n), dims=['time'])
    mu_true = rng.uniform(0., 12., n)
    lamda = rng.uniform(1500., 8000., n)
    N0 = gamma_N0(rng.uniform(50., 3000., n), lamda, mu_true)
    # Perturb M6 so the observed moments are not exactly gamma-consistent
    M0, M3, M6 = (xr.DataArray(gamma_moment(k, N0, lamda, mu_true), dims=['time'])
                  for k in (0, 3, 6))
    M6 = M6 * rng.uniform(0.8, 1.25, n)
    closures = ev.build_closure_gammas(ev.cmr * M3 / rho, M0, M6, rho)
    for name, (N0_c, lamda_c, mu_c) in closures.items():
        np.testing.assert_allclose(gamma_moment(3, N0_c, lamda_c, mu_c), M3, rtol=1.e-10,
                                   err_msg=name)
        if not name.startswith('MM3_'):
            np.testing.assert_allclose(gamma_moment(0, N0_c, lamda_c, mu_c), M0, rtol=1.e-10,
                                       err_msg=name)
    N0_c, lamda_c, mu_c = closures['MM036']
    np.testing.assert_allclose(gamma_moment(6, N0_c, lamda_c, mu_c), M6, rtol=1.e-10)
    N0_c, lamda_c, mu_c = closures['MM036_solve_alpha']
    inside = (closures['MM036'][2] > 0.) & (closures['MM036'][2] < 40.)
    np.testing.assert_allclose(gamma_moment(6, N0_c, lamda_c, mu_c).where(inside),
                               M6.where(inside), rtol=0.1)
    # numpy inputs work too
    closures_np = ev.build_closure_gammas((ev.cmr * M3 / rho).to_numpy(), M0.to_numpy(),
                                          M6.to_numpy(), rho.to_numpy())
    np.testing.assert_allclose(closures_np['MM036'][1], closures['MM036'][1].to_numpy())
