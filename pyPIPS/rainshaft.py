"""
rainshaft.py: 1-D steady, evaporation-only rainshaft models for estimating evaporation aloft and
the evolution of the drop size distribution (DSD) between a starting level (e.g., a radar
sampling height) and the surface, for comparison with PIPS DSDs.

Two models are provided:

- :func:`rainshaft`: a Lagrangian bin model. Each initial size class conserves its steady-state
  number flux F_i = N_i V(D_i) dD_i while its diameter D_i(z) evolves by diffusional
  evaporation; the concentration at height z is F_i / V(D_i(z)). The approach follows the 1-D
  evaporation-only bin models of Kumjian and Ryzhkov (2010; https://doi.org/10.1175/
  2010JAMC2243.1), for light-to-moderate rain, and Xie et al. (2016;
  https://doi.org/10.1175/JTECH-D-15-0244.1), with per-diameter integration as in de Szoeke et
  al. (2024; https://doi.org/10.1029/2024JD041126). Seifert (2008;
  https://doi.org/10.1175/2008JAS2586.1) and Paukert et al. (2019;
  https://doi.org/10.1029/2018MS001512) used rainshaft models as bin benchmarks for bulk schemes.
  Ryzhkov et al. (2020; https://doi.org/10.3390/atmos11040362) and Xie et al. (2016) frame the
  polarimetric radar plus 1-D model approach for estimating evaporative cooling.
- :func:`bulk_rainshaft`: the same column with the rain represented by its moments (q_r, and
  optionally N_t and Z) and the gamma closures of one-, two-, and three-moment bulk schemes
  (see pyPIPS.evaporation.build_closure_gammas and number_tendency), integrated to a steady state
  with upwind sedimentation. Comparing it with :func:`rainshaft` from the same top DSD isolates
  the error a bulk scheme accumulates through the column.

Limitations (both models)
-------------------------
- No collisions, coalescence, or breakup. These shape the DSD and hence evaporation (Morrison et
  al., 2012, https://doi.org/10.1175/MWR-D-11-00283.1; Niquet et al., 2024,
  https://doi.org/10.1029/2024JD041269).
- Steady state: the rain flux entering the top is assumed constant in time. Lolli et al. (2017;
  https://doi.org/10.1175/JTECH-D-16-0146.1) attribute model-observation discrepancies partly to
  non-steady conditions and collisions.
- No horizontal drift. Wind drift of small drops matters when matching radar-aloft DSDs to the
  ground (Dai et al., 2019; https://doi.org/10.1029/2019WR025517). A hook for a future drift
  correction using the PIPS winds is provided (``drift_correction``), but not implemented.
- Starting DSDs retrieved from radar carry their own uncertainty (Huang et al., 2019;
  https://doi.org/10.1175/JTECH-D-18-0107.1), which propagates into the results.
- The PIPS measure thermodynamics only at the surface. The thermodynamic profile is either
  supplied by the user (sounding or model) or ASSUMED well mixed (constant theta and q_v) from the
  PIPS surface state; the choice is recorded in the output metadata.

Units: SI internally. Binned DSDs at the interfaces are in m^-3 mm^-1 with diameters in mm
(pyPIPS conventions); heights are in m (positive upward, typically above the PIPS).

References
----------
TODO(Dan): complete these APA entries (titles, full author lists).

- Dai, Q., et al. (2019). *Water Resources Research*, 55, 9211-9233.
  https://doi.org/10.1029/2019WR025517
- de Szoeke, S. P., et al. (2024). *Journal of Geophysical Research: Atmospheres*, 129.
  https://doi.org/10.1029/2024JD041126
- Huang, H., et al. (2019). *Journal of Atmospheric and Oceanic Technology*, 36, 585-605.
  https://doi.org/10.1175/JTECH-D-18-0107.1
- Kumjian, M. R., & Ryzhkov, A. V. (2010). *Journal of Applied Meteorology and Climatology*,
  49, 1247-1267. https://doi.org/10.1175/2010JAMC2243.1
- Lolli, S., et al. (2017). *Journal of Atmospheric and Oceanic Technology*, 34, 829-839.
  https://doi.org/10.1175/JTECH-D-16-0146.1
- Morrison, H., Tessendorf, S. A., Ikeda, K., & Thompson, G. (2012). *Monthly Weather Review*,
  140, 2437-2460. https://doi.org/10.1175/MWR-D-11-00283.1
- Niquet, L., et al. (2024). *Journal of Geophysical Research: Atmospheres*, 129.
  https://doi.org/10.1029/2024JD041269
- Paukert, M., et al. (2019). *Journal of Advances in Modeling Earth Systems*, 11, 257-277.
  https://doi.org/10.1029/2018MS001512
- Ryzhkov, A. V., et al. (2020). *Atmosphere*, 11, 362. https://doi.org/10.3390/atmos11040362
- Seifert, A. (2008). *Journal of the Atmospheric Sciences*, 65, 3608-3619.
  https://doi.org/10.1175/2008JAS2586.1
- Xie, X., et al. (2016). *Journal of Atmospheric and Oceanic Technology*, 33, 1779-1792.
  https://doi.org/10.1175/JTECH-D-15-0244.1
"""
from __future__ import annotations

import numpy as np
import xarray as xr
from scipy.special import gammaln

from . import DSDlib as dsd
from . import evaporation as evap
from . import parsivel_params as pp
from . import thermolib as thermo

g = 9.81  # Gravitational acceleration (m s^-2)
rhow = evap.rhow
cmr = evap.cmr

PARSIVEL_EDGES_MM = pp.parsivel_parameters['diameter_bin_edges_mm']


# =================================================================================================
# Thermodynamic profiles
# =================================================================================================

class ThermoProfile:
    """Vertical profile of temperature, pressure, and relative humidity.

    Parameters
    ----------
    z : array_like
        Height (m), e.g., above the PIPS
    T : array_like
        Temperature (K)
    p : array_like
        Pressure (Pa)
    RH : array_like
        Relative humidity (fraction)
    source : str
        Short label of where the profile came from (recorded in output metadata)
    attrs : dict, optional
        Additional metadata
    """

    def __init__(self, z, T, p, RH, source='user_supplied', attrs=None):
        z = np.asarray(z, dtype=float)
        order = np.argsort(z)
        self.z = z[order]
        self.T = np.asarray(T, dtype=float)[order]
        self.p = np.asarray(p, dtype=float)[order]
        self.RH = np.asarray(RH, dtype=float)[order]
        self.source = source
        self.attrs = dict(attrs or {})

    @classmethod
    def from_arrays(cls, z, T, p, RH, source='user_supplied', attrs=None):
        """Profile from a sounding or model output (T in K, p in Pa, RH as a fraction)."""
        return cls(z, T, p, RH, source=source, attrs=attrs)

    @classmethod
    def well_mixed_from_surface(cls, T_sfc, p_sfc, RH_sfc, z_top, z_sfc=0., dz=10.):
        """ASSUMED well-mixed subcloud layer built from the PIPS surface state.

        Potential temperature and water vapor mixing ratio are held constant with height and
        pressure is integrated hydrostatically. Above the lifting condensation level (LCL), RH is
        capped at 1 (no evaporation) while T remains dry adiabatic, which is not physical; the LCL
        height is recorded in ``attrs['z_LCL']``. This is an assumption, not an observation, and
        is labeled as such in the output metadata.

        Parameters
        ----------
        T_sfc, p_sfc, RH_sfc : float
            Surface temperature (K), pressure (Pa), and RH (fraction)
        z_top : float
            Top of the profile (m)
        z_sfc : float
            Surface height (m)
        dz : float
            Vertical resolution (m)
        """
        theta = thermo.caltheta(p_sfc, T_sfc)
        w = thermo.calw(RH_sfc, p_sfc, T_sfc)
        z = np.arange(z_sfc, z_top + dz, dz)
        p = np.empty_like(z)
        p[0] = p_sfc

        def Tv_of(pp_):
            T_ = thermo.calT(pp_, theta)
            return T_ * (1. + thermo.rvdrd * w) / (1. + w)

        for k in range(1, z.size):
            dzk = z[k] - z[k - 1]
            p_pred = p[k - 1] * np.exp(-g * dzk / (thermo.Rd * Tv_of(p[k - 1])))
            Tv_mid = 0.5 * (Tv_of(p[k - 1]) + Tv_of(p_pred))
            p[k] = p[k - 1] * np.exp(-g * dzk / (thermo.Rd * Tv_mid))
        T = thermo.calT(p, theta)
        e = p * w / (thermo.rddrv + w)
        RH_raw = e / thermo.cales(p, T)
        above_lcl = RH_raw >= 1.
        z_LCL = float(z[above_lcl][0]) if above_lcl.any() else np.nan
        attrs = {
            'assumption': 'well-mixed subcloud layer (constant theta and water vapor mixing ratio) '
                          'built from the PIPS surface state; not an observed profile',
            'theta_K': float(theta),
            'w_kg_per_kg': float(w),
            'z_LCL': z_LCL,
        }
        return cls(z, T, p, np.minimum(RH_raw, 1.), source='well_mixed_from_PIPS_surface',
                   attrs=attrs)

    @property
    def w(self):
        """Water vapor mixing ratio (kg/kg)."""
        return thermo.calw(self.RH, self.p, self.T)

    def interp(self, z):
        """Interpolate T, p (log-linear), and RH to heights z (m); held constant beyond the ends."""
        z = np.asarray(z, dtype=float)
        T = np.interp(z, self.z, self.T)
        p = np.exp(np.interp(z, self.z, np.log(self.p)))
        RH = np.interp(z, self.z, self.RH)
        return T, p, RH

    def copy(self):
        return ThermoProfile(self.z.copy(), self.T.copy(), self.p.copy(), self.RH.copy(),
                             source=self.source, attrs=dict(self.attrs))

    def metadata(self):
        """Attributes describing the profile, for output metadata."""
        md = {'thermo_profile_source': self.source}
        for key, val in self.attrs.items():
            md[f'thermo_profile_{key}'] = val
        return md


def _env_grid(profile, z_bottom, z_top, thermo_form, saturation_definition, dz_env=2.):
    """Evaporation environment precomputed on a fine height grid for fast interpolation."""
    nz = max(2, int(np.ceil((z_top - z_bottom) / dz_env)) + 1)
    z = np.linspace(z_bottom, z_top, nz)
    T, p, RH = profile.interp(z)
    env = evap.evaporation_environment(T, p, RH, thermo_form=thermo_form,
                                       saturation_definition=saturation_definition)
    keys = ('rho', 'S', 'AB', 'Dv', 'nu', 'Sc')
    return z, {k: np.asarray(env[k], dtype=float) for k in keys}


def _env_at(z_grid, env_grid, z):
    env = {k: np.interp(z, z_grid, v) for k, v in env_grid.items()}
    return env


# =================================================================================================
# Starting conditions
# =================================================================================================

def fine_diameter_edges_mm(D_max_mm=10., dD_mm=0.05):
    """Uniform diameter bin edges (mm) for gamma starting DSDs."""
    return np.arange(0., D_max_mm + 0.5 * dD_mm, dD_mm)


def starting_dsd_from_gamma(N0, lamda, mu, D_edges_mm):
    """Binned DSD (m^-3 mm^-1) at bin centers from SI gamma parameters (as in
    DSDlib.calc_binned_DSD_from_params, converted from m^-4 to m^-3 mm^-1)."""
    edges = np.asarray(D_edges_mm, dtype=float)
    D_mm = 0.5 * (edges[:-1] + edges[1:])
    return 1.e-3 * np.asarray(dsd.calc_binned_DSD_from_params(N0, lamda, mu, D_mm))


def starting_conditions_from_radar_at_PIPS(parsivel_ds, radar_name, retrieval_tag):
    """Starting gamma DSD parameters and height from radar data matched to a PIPS.

    Reads the radar-retrieved N0, lamda, and mu stored as ``N0_<tag>``, ``lamda_<tag>``, and
    ``mu_<tag>`` entries of ``<radar_name>_at_PIPS`` (written by
    analysis_scripts/calc_radar_retrieval_PIPS.py, from DSDlib.retrieval_Cao), and the beam
    height ``<radar_name>_beam_height_at_PIPS`` (written by analysis_scripts/radar_to_PIPS.py).

    Caveats of the current radar matching (to be addressed in Phase 2): the height is that of
    the nearest gate center above the PIPS (beam width and vertical extent not recorded); radar
    values have been interpolated to the PIPS times, so actual sweep times are not retained; the
    stored elevation angle is the requested one.

    Returns
    -------
    xarray.Dataset
        N0 (m^-(4+mu)), lamda (m^-1), mu, and z_top (m above the PIPS), on the PIPS time axis
    """
    da = parsivel_ds[f'{radar_name}_at_PIPS']
    dim = f'fields_{radar_name}'
    N0_mm = da.sel({dim: f'N0_{retrieval_tag}'}, drop=True)  # m^-3 mm^-(1+mu)
    lamda_mm = da.sel({dim: f'lamda_{retrieval_tag}'}, drop=True)  # mm^-1
    mu = da.sel({dim: f'mu_{retrieval_tag}'}, drop=True)
    out = xr.Dataset({
        'N0': N0_mm * 1000.**(1. + mu),
        'lamda': lamda_mm * 1000.,
        'mu': mu,
        'z_top': parsivel_ds[f'{radar_name}_beam_height_at_PIPS'],
    })
    out['N0'].attrs.update(units='m^-(4+mu)')
    out['lamda'].attrs.update(units='m-1')
    out['z_top'].attrs.update(units='m', description='nearest radar gate center height above '
                                                     'the PIPS')
    out.attrs.update(radar_name=radar_name, retrieval_tag=retrieval_tag,
                     caveats='gate-center height; radar values interpolated to PIPS times; '
                             'beam extent not accounted for')
    return out


def _remap_classes(n, D_lo_mm, D_hi_mm, out_edges_mm):
    """Conservatively remap class number concentrations n (m^-3), each spread uniformly over
    [D_lo, D_hi] (mm), onto output bins; returns number density (m^-3 mm^-1)."""
    out_edges = np.asarray(out_edges_mm, dtype=float)
    width = np.maximum(D_hi_mm - D_lo_mm, 1.e-12)
    lo = np.maximum(D_lo_mm[:, None], out_edges[None, :-1])
    hi = np.minimum(D_hi_mm[:, None], out_edges[None, 1:])
    overlap = np.clip(hi - lo, 0., None)
    counts = (n[:, None] * overlap / width[:, None]).sum(axis=0)
    return counts / np.diff(out_edges)


def _bin_coords(edges_mm):
    edges = np.asarray(edges_mm, dtype=float)
    return {
        'diameter': ('diameter_bin', 0.5 * (edges[:-1] + edges[1:])),
        'min_diameter': ('diameter_bin', edges[:-1]),
        'max_diameter': ('diameter_bin', edges[1:]),
    }


def _output_levels(z_top, z_bottom, dz_out, z_out):
    if z_out is not None:
        return np.sort(np.asarray(z_out, dtype=float))[::-1]
    levels = np.arange(z_top, z_bottom, -dz_out)
    return np.append(levels, z_bottom)


def _check_drift(drift_correction):
    if drift_correction is not None:
        # TODO(Dan): implement a drift correction using the PIPS (or sounding) winds; see Dai et
        # al. (2019; https://doi.org/10.1029/2019WR025517) on wind drift of small drops.
        raise NotImplementedError("drift_correction is not implemented yet")


# =================================================================================================
# Lagrangian bin rainshaft
# =================================================================================================

def _march_bin(ND_top, D_edges_mm, z_top, z_bottom, z_out, profile, fall_speed, ventilation,
               thermo_form, saturation_definition, n_sub, D_min, dz_step, out_edges_mm,
               max_steps):
    """Integrate the Lagrangian classes and diagnose outputs at levels z_out (descending)."""
    edges = np.asarray(D_edges_mm, dtype=float)
    ND_top = np.nan_to_num(np.asarray(ND_top, dtype=float))

    # Subdivide each bin into n_sub classes with the same number density (piecewise constant)
    sub_edges = np.concatenate([np.linspace(e0, e1, n_sub + 1)[:-1]
                                for e0, e1 in zip(edges[:-1], edges[1:])] + [edges[-1:]])
    N_sub = np.repeat(ND_top, n_sub)  # m^-3 mm^-1
    lo_mm, hi_mm = sub_edges[:-1], sub_edges[1:]
    keep = N_sub > 0.
    lo_mm, hi_mm, N_sub = lo_mm[keep], hi_mm[keep], N_sub[keep]
    c_mm = 0.5 * (lo_mm + hi_mm)
    ncls = c_mm.size

    z_grid, env_grid = _env_grid(profile, min(z_bottom, z_top), z_top, thermo_form,
                                 saturation_definition)
    a_v, b_v = evap._ventilation_coeffs(ventilation)

    def rates(z, x):
        D = np.sqrt(np.maximum(x, 0.))
        env = _env_at(z_grid, env_grid, z)
        V = fall_speed(D, env['rho'])
        fv = evap.ventilation_factor(D, V, env['nu'], env['Dv'], a_v, b_v)
        dxdt = 8. * fv * (env['S'] - 1.) / (rhow * env['AB'])
        return -V, dxdt

    # Number flux of each class at the top (m^-2 s^-1), conserved along the class trajectory
    env_top = _env_at(z_grid, env_grid, np.array([z_top]))
    n_top = 1000. * N_sub * (hi_mm - lo_mm) * 1.e-3  # m^-4 * m -> m^-3
    F = n_top * fall_speed(c_mm * 1.e-3, env_top['rho'])

    # Lagrangian points: lower edges, centers, upper edges of each class
    D0 = np.concatenate([lo_mm, c_mm, hi_mm]) * 1.e-3
    x = D0**2.
    z = np.full(x.shape, float(z_top))
    x_min = D_min**2.
    alive = x > x_min
    landed = np.zeros(x.shape, dtype=bool)
    Z_hist, X_hist = [z.copy()], [x.copy()]

    nsteps = 0
    while (alive & ~landed).any():
        active = alive & ~landed
        dzdt0, _ = rates(z[active], x[active])
        dt = dz_step / max(np.max(-dzdt0), 1.e-6)
        # Classic RK4 on (z, x) for active points
        za, xa = z[active], x[active]
        k1z, k1x = rates(za, xa)
        k2z, k2x = rates(za + 0.5 * dt * k1z, np.maximum(xa + 0.5 * dt * k1x, 0.))
        k3z, k3x = rates(za + 0.5 * dt * k2z, np.maximum(xa + 0.5 * dt * k2x, 0.))
        k4z, k4x = rates(za + dt * k3z, np.maximum(xa + dt * k3x, 0.))
        z[active] = za + dt / 6. * (k1z + 2. * k2z + 2. * k3z + k4z)
        x[active] = np.maximum(xa + dt / 6. * (k1x + 2. * k2x + 2. * k3x + k4x), 0.)
        died = active & (x <= x_min)
        x[died] = 0.
        alive &= ~died
        landed |= alive & (z <= z_bottom)
        Z_hist.append(z.copy())
        X_hist.append(x.copy())
        nsteps += 1
        if nsteps >= max_steps:
            raise RuntimeError(f"rainshaft did not finish within max_steps={max_steps}")

    Z_hist = np.array(Z_hist)
    X_hist = np.array(X_hist)

    # Interpolate each point's trajectory x(z) to the output levels
    npts = x.size
    X_out = np.zeros((z_out.size, npts))
    for i in range(npts):
        zs, xs = Z_hist[:, i], X_hist[:, i]
        # Trajectory is monotonic in z until the point died or landed (then frozen)
        moving = np.concatenate([[True], np.diff(zs) < 0.])
        last = np.nonzero(moving)[0][-1]
        zs, xs = zs[:last + 1], xs[:last + 1]
        xi = np.interp(-z_out, -zs, xs)
        xi[z_out < zs[-1]] = 0. if xs[-1] <= 0. else xs[-1]
        X_out[:, i] = xi

    D_lo = np.sqrt(X_out[:, :ncls])
    D_c = np.sqrt(X_out[:, ncls:2 * ncls])
    D_hi = np.sqrt(X_out[:, 2 * ncls:])
    cls_alive = D_c > D_min

    # Diagnostics at each output level
    T, p, RH = profile.interp(z_out)
    env = evap.evaporation_environment(T, p, RH, thermo_form=thermo_form,
                                       saturation_definition=saturation_definition)
    rho = np.asarray(env['rho'])
    nz = z_out.size
    E = np.zeros(nz)
    M = {k: np.zeros(nz) for k in (0, 3, 4, 6)}
    mass_flux = np.zeros(nz)
    D0_med = np.full(nz, np.nan)
    nbins_out = len(out_edges_mm) - 1
    ND_out = np.zeros((nz, nbins_out))
    for k in range(nz):
        a = cls_alive[k]
        if not a.any():
            continue
        Dk = D_c[k, a]
        env_k = {key: (val[k] if np.ndim(val) else val) for key, val in env.items()}
        Vk = fall_speed(Dk, env_k['rho'])
        n = F[a] / Vk  # m^-3
        dmdt = evap.drop_mass_rate(Dk, fall_speed=fall_speed, ventilation=ventilation, env=env_k)
        E[k] = np.sum(n * -dmdt)
        for mom in M:
            M[mom][k] = np.sum(n * Dk**mom)
        mass_flux[k] = np.sum(F[a] * cmr * Dk**3.)  # kg m^-2 s^-1
        order = np.argsort(Dk)
        cum = np.cumsum((n * Dk**3.)[order])
        if cum[-1] > 0.:
            D0_med[k] = np.interp(0.5 * cum[-1], cum, Dk[order])
        ND_out[k] = _remap_classes(n, D_lo[k, a] * 1.e3, D_hi[k, a] * 1.e3, out_edges_mm)

    with np.errstate(divide='ignore', invalid='ignore'):
        diag = {
            'E': E,
            'qr': cmr * M[3] / rho,
            'rainrate': mass_flux / rhow * 3.6e6,
            'Nt': M[0],
            'reflectivity': 10. * np.log10(1.e18 * M[6]),
            'Dm43': M[4] / M[3],
            'Dmv': (M[3] / M[0])**(1. / 3.),
            'D0': D0_med,
            'rho': rho,
            'T': T,
            'p': p,
            'RH': RH,
        }
    return diag, ND_out, nsteps


_RAINSHAFT_ATTRS = {
    'E': ('kg m-3 s-1', 'evaporation rate (positive = evaporation)'),
    'qr': ('kg kg-1', 'rain mixing ratio'),
    'rainrate': ('mm h-1', 'rain rate (liquid-equivalent mass flux)'),
    'Nt': ('m-3', 'total number concentration'),
    'reflectivity': ('dBZ', 'Rayleigh reflectivity factor'),
    'Dm43': ('m', 'mass-weighted mean diameter M4/M3'),
    'Dmv': ('m', 'mean volume diameter (M3/M0)^(1/3)'),
    'D0': ('m', 'median volume diameter'),
    'cooling_rate': ('K s-1', 'latent cooling rate (negative = cooling)'),
    'cooling_rate_K_per_h': ('K h-1', 'latent cooling rate (negative = cooling)'),
    'dqv_dt': ('kg kg-1 s-1', 'water vapor mixing-ratio tendency (moistening)'),
    'frac_rain_reaching': ('1', 'rain rate relative to the top level'),
    'rho': ('kg m-3', 'air density'),
    'T': ('K', 'temperature'),
    'p': ('Pa', 'pressure'),
    'RH': ('1', 'relative humidity'),
    'ND': ('m-3 mm-1', 'number density on the output diameter bins'),
    'N0': ('m^-(4+mu)', 'gamma intercept parameter'),
    'lamda': ('m-1', 'gamma slope parameter'),
    'mu': ('1', 'gamma shape parameter'),
}


def _to_dataset(diag, ND_out, z_out, out_edges_mm, height_dim='height'):
    tend = evap.latent_tendencies(diag['E'], diag['rho'])
    diag = dict(diag)
    diag['cooling_rate'] = tend['cooling_rate']
    diag['cooling_rate_K_per_h'] = tend['cooling_rate_K_per_h']
    diag['dqv_dt'] = tend['dqv_dt']
    with np.errstate(divide='ignore', invalid='ignore'):
        diag['frac_rain_reaching'] = diag['rainrate'] / diag['rainrate'][0]
    ds = xr.Dataset({k: ((height_dim,), np.asarray(v)) for k, v in diag.items()},
                    coords={height_dim: z_out})
    ds['ND'] = ((height_dim, 'diameter_bin'), ND_out)
    ds = ds.assign_coords(_bin_coords(out_edges_mm))
    ds['ND_surface'] = ds['ND'].isel({height_dim: -1}, drop=True)
    ds['ND_surface'].attrs.update(units='m-3 mm-1',
                                  description='predicted DSD at the bottom level on the output '
                                              'diameter bins')
    for var, (units, desc) in _RAINSHAFT_ATTRS.items():
        if var in ds:
            ds[var].attrs.update(units=units, description=desc)
    ds[height_dim].attrs.update(units='m', description='height (positive upward)')
    for c in ('diameter', 'min_diameter', 'max_diameter'):
        ds[c].attrs.update(units='mm')
    return ds


def rainshaft(ND_top, D_edges_mm, z_top, z_bottom, profile, fall_speed=None,
              ventilation=evap.DEFAULT_VENTILATION, thermo_form='textbook',
              saturation_definition='vapor_pressure', n_sub=10, dz_out=50., z_out=None,
              D_min=10.e-6, dz_step=2., out_edges_mm=None, feedback=False, feedback_time=None,
              n_feedback_steps=10, drift_correction=None, max_steps=500000):
    """1-D steady, evaporation-only Lagrangian bin rainshaft (aloft to surface).

    For each initial size class i, the steady-state number flux F_i = N_i V(D_i) dD_i is
    conserved while D_i(z) evolves downward by diffusional evaporation. The concentration at
    height z is F_i / V(D_i(z)). Classes whose diameter reaches D_min are removed (fully
    evaporated). Each class is followed as a Lagrangian drop in time, (dz/dt, d(D^2)/dt) =
    (-V, 8 f_v (S-1)/(rho_w (A+B))), with classic RK4. Integrating D^2 rather than D keeps the
    rate finite as D -> 0. Class edges are also followed so the DSD can be remapped
    conservatively onto the output bins (by default, the Parsivel bins).

    Parameters
    ----------
    ND_top : array_like
        Binned DSD at the starting height (m^-3 mm^-1)
    D_edges_mm : array_like
        Bin edges of ND_top (mm), length len(ND_top) + 1
    z_top, z_bottom : float
        Starting and bottom heights (m), z_top > z_bottom
    profile : ThermoProfile
        Thermodynamic profile covering [z_bottom, z_top]
    fall_speed : pyPIPS.evaporation.FallSpeed, optional
        Default AtlasFallSpeed() (as for the binned PIPS calculation)
    ventilation, thermo_form, saturation_definition :
        See pyPIPS.evaporation.evaporation_rate_gamma()
    n_sub : int
        Number of sub-classes per input bin (piecewise-constant number density within each bin)
    dz_out : float
        Output level spacing (m), used if z_out is None
    z_out : array_like, optional
        Output heights (m)
    D_min : float
        Diameter (m) at which a class is considered fully evaporated
    dz_step : float
        Integration step, as the vertical distance (m) traveled by the fastest active drop
    out_edges_mm : array_like, optional
        Output bin edges (mm); default Parsivel bin edges
    feedback : bool
        If True, apply the latent cooling and moistening to T and q_v (quasi-steady
        approximation; see below). Default False (instantaneous diagnostic for a fixed profile).
    feedback_time : float
        Exposure time tau (s) for feedback: how long the column has been under the rain.
    n_feedback_steps : int
        Number of quasi-steady steps over feedback_time
    drift_correction : optional
        Hook for a future horizontal drift correction; not implemented

    Feedback (optional, off by default)
    -----------------------------------
    A strictly steady column cannot keep cooling, so feedback requires an exposure time tau.
    The quasi-steady approximation splits tau into n_feedback_steps; each step solves the steady
    rainshaft for the current profile, then applies dT = cooling_rate * dt and dq_v = (E/rho) dt
    at the profile levels (pressure held fixed) before recomputing RH. This assumes the DSD adjusts
    to the changing air much faster than the air changes, which holds better for millimeter-size
    drops than for the smallest drops. tau and the approximation are recorded in the metadata.

    Notes on the diagnostics
    ------------------------
    Number flux is conserved per class, so the concentration F_i / V(D_i) of shrinking (slowing)
    drops can increase downward before they vanish; N_t need not decrease monotonically.
    Likewise, D_mv = (M3/M0)^(1/3) can first decrease near the top (all drops shrink, none has
    vanished yet) before increasing as small drops evaporate completely. Moment-weighted sizes
    such as Dm43 are better behaved. Discrete class removal adds small-scale noise to M0-based
    quantities; increase n_sub (or use finer input bins) to reduce it.

    Returns
    -------
    xarray.Dataset
        On dimension ``height``: E, qr, rainrate, Nt, reflectivity, Dm43, Dmv, D0, cooling_rate,
        cooling_rate_K_per_h, dqv_dt, frac_rain_reaching, T, p, RH, rho; ND (height x
        diameter_bin) and ND_surface (diameter_bin) on the output bins.
    """
    _check_drift(drift_correction)
    if z_top <= z_bottom:
        raise ValueError("z_top must be greater than z_bottom")
    if fall_speed is None:
        fall_speed = evap.AtlasFallSpeed()
    if out_edges_mm is None:
        out_edges_mm = PARSIVEL_EDGES_MM
    z_levels = _output_levels(z_top, z_bottom, dz_out, z_out)
    args = (fall_speed, ventilation, thermo_form, saturation_definition, n_sub, D_min, dz_step,
            out_edges_mm, max_steps)

    feedback_attrs = {'feedback': int(bool(feedback))}
    if feedback:
        if feedback_time is None or feedback_time <= 0.:
            raise ValueError("feedback=True requires a positive feedback_time (s)")
        profile = profile.copy()
        in_shaft = (profile.z >= z_bottom) & (profile.z <= z_top)
        z_prof = profile.z[in_shaft][::-1]
        dt_fb = feedback_time / n_feedback_steps
        for _ in range(n_feedback_steps):
            diag, _, _ = _march_bin(ND_top, D_edges_mm, z_top, z_bottom, z_prof, profile, *args)
            tend = evap.latent_tendencies(diag['E'], diag['rho'])
            w = profile.w
            T_new = profile.T.copy()
            T_new[in_shaft] += (tend['cooling_rate'] * dt_fb)[::-1]
            w[in_shaft] += (tend['dqv_dt'] * dt_fb)[::-1]
            e = profile.p * w / (thermo.rddrv + w)
            RH_new = e / thermo.cales(profile.p, T_new)
            profile = ThermoProfile(profile.z, T_new, profile.p, RH_new, source=profile.source,
                                    attrs=profile.attrs)
        feedback_attrs.update({
            'feedback_time_s': float(feedback_time),
            'n_feedback_steps': int(n_feedback_steps),
            'feedback_assumption': 'quasi-steady: steady DSD re-solved for the evolving '
                                   'profile; pressure held fixed; no advection or mixing',
        })

    diag, ND_out, nsteps = _march_bin(ND_top, D_edges_mm, z_top, z_bottom, z_levels, profile,
                                      *args)
    ds = _to_dataset(diag, ND_out, z_levels, out_edges_mm)
    ds.attrs.update({
        'model': 'Lagrangian bin rainshaft (steady, evaporation only)',
        'fall_speed': fall_speed.describe(),
        'ventilation': str(evap._ventilation_coeffs(ventilation)),
        'thermo_form': thermo_form,
        'saturation_definition': saturation_definition,
        'n_sub': n_sub,
        'D_min_m': D_min,
        'dz_step_m': dz_step,
        'n_integration_steps': nsteps,
        'limitations': 'no collisions/breakup; steady state; no horizontal drift',
        **feedback_attrs,
        **profile.metadata(),
    })
    return ds


# =================================================================================================
# Bulk-scheme rainshaft
# =================================================================================================

def moment_weighted_fallspeed(k, lamda, mu, fall_speed, rho):
    """M_k-weighted fall speed (m s^-1) of a gamma DSD for a power-law or Ferrier-type fall speed:

        V_k = a dens Gamma(1+mu+k+b) / Gamma(1+mu+k) * lamda^(1+mu+k) / (lamda+f)^(1+mu+k+b)

    with f = 0 for a power law (dens is the density correction). For k = 3 and the Ferrier
    constants this is the mass-weighted fall speed of DSDlib.calc_VQR_ferrier (which uses a
    (1.204/rho)^0.5 density correction).
    """
    if isinstance(fall_speed, evap.FerrierFallSpeed):
        f = fall_speed.f
    elif isinstance(fall_speed, evap.PowerLawFallSpeed):
        f = 0.
    else:
        raise TypeError("bulk_rainshaft requires a PowerLawFallSpeed or FerrierFallSpeed")
    a, b = fall_speed.a, fall_speed.b
    s = 1. + mu + k
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        log_v = (np.log(a) + gammaln(s + b) - gammaln(s) + s * np.log(lamda) -
                 (s + b) * np.log(lamda + f))
        return np.exp(log_v) * fall_speed.density_factor(rho)


def _diagnose_gamma(rho, Q, N, Zm, n_moments, N0_1mom, mu_1mom, mu_closure, mu_fixed,
                    q_min, N_min):
    """Gamma parameters (N0, lamda, mu) and validity mask from the prognostic moments."""
    q = Q / rho
    valid = q > q_min
    if n_moments >= 2:
        valid &= N > N_min
    qs = np.where(valid, q, 1.e-6)
    Ns = np.where(valid, N, 1.e3)
    if n_moments == 1:
        mu = np.full(q.shape, float(mu_1mom))
        N0 = np.full(q.shape, float(N0_1mom))
        Ns = np.asarray(dsd.calc_Nt_gamma(rho, qs, N0, cmr, mu))
        lamda = np.asarray(dsd.calc_lamda_gamma(rho, qs, Ns, cmr, mu))
    else:
        if n_moments == 3:
            Zs = np.where(valid, Zm, 1.e-9)
            mu = np.asarray(dsd.solve_alpha(rho, cmr, qs, Ns, Zs))
        elif callable(mu_closure):
            mu = np.asarray(mu_closure(np.asarray(dsd.calc_Dm_gamma(rho, qs, Ns, cmr))),
                            dtype=float) * np.ones(q.shape)
        elif mu_closure == 'diag_alpha':
            mu = np.asarray(dsd.diag_alpha('qr', np.asarray(dsd.calc_Dm_gamma(rho, qs, Ns, cmr))))
        elif mu_closure == 'fixed':
            mu = np.full(q.shape, float(mu_fixed))
        else:
            raise ValueError(f"Unknown mu_closure '{mu_closure}'")
        lamda = np.asarray(dsd.calc_lamda_gamma(rho, qs, Ns, cmr, mu))
        N0 = np.asarray(dsd.calc_N0_gamma(rho, qs, Ns, cmr, mu))
    return N0, lamda, mu, Ns, valid


def bulk_rainshaft(z_top, z_bottom, profile, qr_top=None, Nt_top=None, Z_top=None, ND_top=None,
                   D_edges_mm=None, n_moments=2, mu_closure='fixed', mu_fixed=0.,
                   N0_1mom=8.e6, mu_1mom=0., number_method='proportional',
                   fall_speed=None, ventilation=evap.DEFAULT_VENTILATION,
                   thermo_form='textbook', saturation_definition='vapor_pressure', dz=20.,
                   cfl=0.8, t_max=10800., check_interval=60., steady_tol=1.e-6,
                   out_edges_mm=None, q_min=1.e-10, N_min=1.e-3, drift_correction=None):
    """1-D evaporation-only rainshaft with rain represented as in a bulk microphysics scheme.

    The rain is represented by prognostic moments (Q = rho q_r; plus N_t for two- and
    three-moment; plus Z = M6 for three-moment) on a uniform height grid. At each level, the gamma
    DSD is diagnosed with the scheme's closure; evaporation of mass comes from the gamma closed
    form (pyPIPS.evaporation.evaporation_rate_gamma), the number tendency from
    pyPIPS.evaporation.number_tendency, and the Z tendency from
    pyPIPS.evaporation.moment_tendency_gamma(6) (computed consistently from the assumed gamma; a
    specific three-moment scheme may differ, TODO(Dan): compare with the scheme(s) of interest).
    Sedimentation uses moment-weighted fall speeds with first-order upwind fluxes. The top level
    is held fixed and the column is integrated in time until it reaches a steady state.

    Closures:
    - n_moments=1: fixed N0 (N0_1mom) and mu (mu_1mom); N_t is diagnosed.
    - n_moments=2: mu from mu_closure: 'fixed' (mu_fixed), 'diag_alpha' (DSDlib.diag_alpha), or
      a callable mu = f(D_mv) with D_mv in m.
    - n_moments=3: mu from DSDlib.solve_alpha(q_r, N_t, Z).

    The top values come from qr_top (kg/kg), Nt_top (m^-3), and Z_top (m^6 m^-3), or from the
    moments of a binned ND_top (m^-3 mm^-1) on D_edges_mm.

    Parameters
    ----------
    z_top, z_bottom : float
        Top and bottom heights (m)
    profile : ThermoProfile
        Thermodynamic profile
    fall_speed : PowerLawFallSpeed or FerrierFallSpeed
        Default FerrierFallSpeed(). Use the same relation in rainshaft() for a like-for-like
        bin-vs-bulk comparison.
    number_method : str
        See pyPIPS.evaporation.number_tendency (two- and three-moment only). The default,
        'proportional', holds the mean volume diameter constant during evaporation.
    dz : float
        Grid spacing (m)
    cfl : float
        Courant number for the time step
    t_max : float
        Maximum integration time (s)
    check_interval, steady_tol : float
        The column is considered steady when the maximum change of Q and N over check_interval
        (s), relative to their column maxima, is below steady_tol.

    Returns
    -------
    xarray.Dataset
        As rainshaft(), plus N0, lamda, mu on ``height``; attrs record the closure and whether
        a steady state was reached.
    """
    _check_drift(drift_correction)
    if z_top <= z_bottom:
        raise ValueError("z_top must be greater than z_bottom")
    if fall_speed is None:
        fall_speed = evap.FerrierFallSpeed()
    if out_edges_mm is None:
        out_edges_mm = PARSIVEL_EDGES_MM

    nlev = max(2, int(round((z_top - z_bottom) / dz)) + 1)
    z = np.linspace(z_top, z_bottom, nlev)
    dzg = (z_top - z_bottom) / (nlev - 1)
    T, p, RH = profile.interp(z)
    env = evap.evaporation_environment(T, p, RH, thermo_form=thermo_form,
                                       saturation_definition=saturation_definition)
    env = {k: (np.asarray(v, dtype=float) if not isinstance(v, str) else v)
           for k, v in env.items()}
    rho = env['rho']

    if ND_top is not None:
        edges = np.asarray(D_edges_mm, dtype=float)
        Dc = 0.5 * (edges[:-1] + edges[1:]) * 1.e-3
        dD = np.diff(edges) * 1.e-3
        n = 1000. * np.nan_to_num(np.asarray(ND_top, dtype=float)) * dD
        qr_top = cmr * np.sum(n * Dc**3.) / rho[0]
        Nt_top = np.sum(n)
        Z_top = np.sum(n * Dc**6.)
    if qr_top is None:
        raise ValueError("Provide qr_top (and Nt_top/Z_top as needed) or ND_top")

    Q = np.zeros(nlev)
    N = np.zeros(nlev)
    Zm = np.zeros(nlev)
    Q[0] = rho[0] * qr_top
    N[0] = Nt_top if Nt_top is not None else 0.
    Zm[0] = Z_top if Z_top is not None else 0.
    closure_args = (n_moments, N0_1mom, mu_1mom, mu_closure, mu_fixed, q_min, N_min)

    def diagnose(Q, N, Zm):
        N0, lamda, mu, Nd, valid = _diagnose_gamma(rho, Q, N, Zm, *closure_args)
        if n_moments == 1:
            N = np.where(valid, Nd, 0.)
        VQ = np.where(valid, moment_weighted_fallspeed(3, lamda, mu, fall_speed, rho), 0.)
        VN = np.where(valid, moment_weighted_fallspeed(0, lamda, mu, fall_speed, rho), 0.)
        VZ = np.where(valid, moment_weighted_fallspeed(6, lamda, mu, fall_speed, rho), 0.)
        return N0, lamda, mu, N, valid, VQ, VN, VZ

    N0, lamda, mu, Ndiag, valid, VQ, VN, VZ = diagnose(Q, N, Zm)
    if n_moments == 1:
        N[0] = Ndiag[0]
    if n_moments == 3 and Z_top is None:
        raise ValueError("n_moments=3 requires Z_top or ND_top")

    t = 0.
    t_next_check = check_interval
    Q_prev, N_prev = Q.copy(), N.copy()
    converged = False
    while t < t_max:
        N0, lamda, mu, Ndiag, valid, VQ, VN, VZ = diagnose(Q, N, Zm)
        if n_moments == 1:
            N = Ndiag
        vmax = max(np.max(VQ), np.max(VZ), 1.e-3)
        dt = cfl * dzg / vmax

        E = np.where(valid, evap.evaporation_rate_gamma(
            N0, lamda, mu, fall_speed=fall_speed, ventilation=ventilation, env=env), 0.)
        E = np.minimum(E, Q / dt)  # Cannot evaporate more than is present
        dQ_evap = -E
        if n_moments >= 2:
            q = Q / rho
            dN_evap = np.asarray(evap.number_tendency(q, N, dQ_evap / rho, rho,
                                                      method=number_method))
            dN_evap = np.maximum(dN_evap, -N / dt)
        if n_moments == 3:
            dZ_evap = np.where(valid, evap.moment_tendency_gamma(
                6, N0, lamda, mu, fall_speed=fall_speed, ventilation=ventilation, env=env), 0.)
            dZ_evap = np.maximum(dZ_evap, -Zm / dt)

        # Upwind sedimentation fluxes (downward positive); level 0 is the fixed top boundary
        FQ = VQ * Q
        Q[1:] += dt * ((FQ[:-1] - FQ[1:]) / dzg + dQ_evap[1:])
        Q = np.maximum(Q, 0.)
        if n_moments >= 2:
            FN = VN * N
            N[1:] += dt * ((FN[:-1] - FN[1:]) / dzg + dN_evap[1:])
            N = np.maximum(N, 0.)
        if n_moments == 3:
            FZ = VZ * Zm
            Zm[1:] += dt * ((FZ[:-1] - FZ[1:]) / dzg + dZ_evap[1:])
            Zm = np.maximum(Zm, 0.)
        t += dt

        if t >= t_next_check:
            t_next_check += check_interval
            dQ = np.max(np.abs(Q - Q_prev)) / max(np.max(Q), 1.e-30)
            dN = (np.max(np.abs(N - N_prev)) / max(np.max(N), 1.e-30)
                  if n_moments >= 2 else 0.)
            Q_prev, N_prev = Q.copy(), N.copy()
            if max(dQ, dN) < steady_tol and Q[-1] > 0.:
                converged = True
                break

    N0, lamda, mu, Ndiag, valid, VQ, VN, VZ = diagnose(Q, N, Zm)
    if n_moments == 1:
        N = Ndiag
    E = np.where(valid, evap.evaporation_rate_gamma(
        N0, lamda, mu, fall_speed=fall_speed, ventilation=ventilation, env=env), 0.)
    with np.errstate(divide='ignore', invalid='ignore'):
        M3 = Q / cmr
        Zout = np.where(valid, np.asarray(dsd.calc_Mp(6., lamda, N, mu)), 0.)
        diag = {
            'E': E,
            'qr': Q / rho,
            'rainrate': VQ * Q / rhow * 3.6e6,
            'Nt': np.where(valid, N, 0.),
            'reflectivity': 10. * np.log10(1.e18 * Zout),
            'Dm43': np.where(valid, (4. + mu) / lamda, np.nan),
            'Dmv': np.where(valid, (M3 / N)**(1. / 3.), np.nan),
            'D0': np.where(valid, (3.672 + mu) / lamda, np.nan),
            'rho': rho,
            'T': T,
            'p': p,
            'RH': RH,
            'N0': np.where(valid, N0, np.nan),
            'lamda': np.where(valid, lamda, np.nan),
            'mu': np.where(valid, mu, np.nan),
        }
    out_edges = np.asarray(out_edges_mm, dtype=float)
    D_mm = 0.5 * (out_edges[:-1] + out_edges[1:])
    ND_out = np.zeros((nlev, D_mm.size))
    for k in np.nonzero(valid)[0]:
        ND_out[k] = 1.e-3 * dsd.calc_binned_DSD_from_params(N0[k], lamda[k], mu[k], D_mm)
    ds = _to_dataset(diag, ND_out, z, out_edges_mm)
    ds['D0'].attrs['description'] = 'median volume diameter, (3.672 + mu)/lamda approximation'
    closure = {1: f'1M: N0={N0_1mom:g}, mu={mu_1mom:g}',
               2: f'2M: mu_closure={mu_closure if not callable(mu_closure) else "callable"}'
                  + (f', mu={mu_fixed:g}' if mu_closure == 'fixed' else '')
                  + f'; number_method={number_method}',
               3: f'3M: solve_alpha; number_method={number_method}'}[n_moments]
    ds.attrs.update({
        'model': 'bulk rainshaft (evaporation only, upwind sedimentation, run to steady state)',
        'closure': closure,
        'n_moments': n_moments,
        'fall_speed': fall_speed.describe(),
        'ventilation': str(evap._ventilation_coeffs(ventilation)),
        'thermo_form': thermo_form,
        'saturation_definition': saturation_definition,
        'dz_m': dzg,
        'integration_time_s': t,
        'converged_to_steady_state': int(converged),
        'limitations': 'no collisions/breakup; no horizontal drift; gamma DSD restriction',
        **profile.metadata(),
    })
    return ds
