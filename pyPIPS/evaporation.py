"""
evaporation.py: raindrop evaporation rates and latent cooling rates from binned (observed) and
gamma drop size distributions (DSDs), bulk-microphysics closure comparisons, and rain
number-concentration tendency closures.

Units
-----
All internal calculations are in SI units. At the interfaces:

- Binned DSDs follow pyPIPS conventions: ``ND`` in m^-3 mm^-1 with dimension ``diameter_bin``
  and coordinates ``diameter``, ``min_diameter``, and ``max_diameter`` in mm.
- Gamma DSDs, N(D) = N0 D^mu exp(-lamda D), use SI parameters: N0 in m^-(4+mu), lamda in m^-1,
  mu dimensionless (this matches the ``DSD_MM*``/``DSD_TMM*`` fits written by
  ``analysis_scripts/calc_MM_fits.py``, which call the shape parameter ``alpha``).
- Single-drop functions take the diameter D in m.
- Thermodynamic inputs are T in K, p in Pa, RH as a fraction (not percent), and air density rho
  in kg m^-3. Use :func:`thermo_from_parsivel_ds` to convert pyPIPS native variables
  (deg C, hPa, percent) and recompute rho consistently.

Sign conventions
----------------
- The evaporation rate E (kg m^-3 s^-1) is positive when liquid is lost to vapor.
- The rain mixing-ratio tendency is dqr/dt = -E / rho (kg kg^-1 s^-1).
- The latent cooling rate is dT/dt = -L_v E / (rho c_p) in K s^-1 (negative = cooling); it is
  also reported in K h^-1. It is a single-process tendency, not a closed heat budget: advection,
  mixing, surface fluxes, radiation, and other phase changes are not included.
- Single-drop dm/dt and dD/dt are negative for an evaporating drop.

These differ from the legacy :func:`pyPIPS.DSDlib.calc_evap`, which returns
``QREVP = dqr/dt`` (negative for evaporation) and ``COOL_RATE = (L_v/c_p) QREVP`` (negative for
cooling). ``calc_evap`` is kept as a backward-compatible wrapper around
:func:`evaporation_rate_gamma` with its original settings.

Scope
-----
A surface DSD gives the evaporation rate only at the surface. The vertical structure of
evaporation requires marching a DSD downward from aloft (see :mod:`pyPIPS.rainshaft`).

References
----------
Only the specific points noted where each is cited are attributed to these works. Titles are
omitted here. TODO(Dan): complete these APA entries (titles, full author lists).

- Brandes, E. A., Zhang, G., & Sun, J. (2006). *Journal of Applied Meteorology and
  Climatology*, 45, 259-268. https://doi.org/10.1175/JAM2335.1
- Dawson, D. T., Xue, M., Milbrandt, J. A., & Yau, M. K. (2010). *Monthly Weather Review*,
  138, 1152-1171. https://doi.org/10.1175/2009MWR2956.1
- Dearden, C., et al. (2010). *Atmospheric Chemistry and Physics Discussions*.
  https://doi.org/10.5194/acpd-10-23497-2010
- de Szoeke, S. P., et al. (2024). *Journal of Geophysical Research: Atmospheres*, 129.
  https://doi.org/10.1029/2024JD041126
- Kuba, N., et al. (2020). *Journal of Advances in Modeling Earth Systems*, 12.
  https://doi.org/10.1029/2019MS001891
- Kumjian, M. R., & Ryzhkov, A. V. (2010). *Journal of Applied Meteorology and Climatology*,
  49, 1247-1267. https://doi.org/10.1175/2010JAMC2243.1
- Liu, X. (2024). *Remote Sensing*, 16, 4297. https://doi.org/10.3390/rs16224297
- Morrison, H., & Milbrandt, J. A. (2011). *Monthly Weather Review*, 139, 1103-1130.
  https://doi.org/10.1175/2010MWR3433.1
- Paukert, M., et al. (2019). *Journal of Advances in Modeling Earth Systems*, 11, 257-277.
  https://doi.org/10.1029/2018MS001512
- Seifert, A. (2008). *Journal of the Atmospheric Sciences*, 65, 3608-3619.
  https://doi.org/10.1175/2008JAS2586.1
- Snook, N., & Xue, M. (2008). *Geophysical Research Letters*, 35.
  https://doi.org/10.1029/2008GL035866
- Srivastava, R. C., & Coen, J. L. (1992). *Journal of the Atmospheric Sciences*, 49,
  1643-1651. https://doi.org/10.1175/1520-0469(1992)049<1643:NEEFTA>2.0.CO;2
- Varble, A., et al. (2014). *Journal of Geophysical Research: Atmospheres*, 119.
  https://doi.org/10.1002/2013JD021372

Attributions already present in pyPIPS are retained: Ferrier (1994) ventilation coefficients
and rain fall-speed constants, the Atlas et al. (1973)-type empirical fall-speed relation, the
Foote and du Toit (1969) density correction, and Jason Milbrandt's microphysics code.
"""
from __future__ import annotations

import warnings

import numpy as np
import xarray as xr
from scipy.special import gammainc, gammaln

from . import DSDlib as dsd
from . import PIPS
from . import thermolib as thermo

rhow = thermo.rhow  # Density of liquid water (kg m^-3)
cmr = dsd.cmr  # pi * rhow / 6, mass-diameter coefficient for spherical drops (kg m^-3)

# =================================================================================================
# Shape-parameter conventions (comment block)
# =================================================================================================
# Bulk schemes and studies differ in how they define the DSD shape parameter:
# - Kuba et al. (2020; https://doi.org/10.1029/2019MS001891) summarize Seifert (2008) as
#   recommending nu_r = -1/3 for a generalized gamma distribution in mass space, while finding
#   nu_r = 1 adequate for shallow cumulus.
# - Varble et al. (2014; https://doi.org/10.1002/2013JD021372) found that mu = 0 causes excessive
#   size sorting in two-moment schemes and that mu = 2.5 improves one-moment biases.
# - Dearden et al. (2010; https://doi.org/10.5194/acpd-10-23497-2010) found diagnosed mu larger
#   below cloud base.
# The mass-space nu and the diameter-space mu used here are NOT interchangeable. For spherical
# drops (m proportional to D^3), a mass-space shape parameter maps to a diameter-space one only
# under specific functional forms. Do not convert one to the other without checking the paper.
# Seifert (2008; https://doi.org/10.1175/2008JAS2586.1) parameterized mu as a function of the
# mean volume diameter. TODO(Dan): add the mu(D_mv) relation and coefficients from Seifert (2008);
# they are deliberately not reproduced here (unverified). Pass it as ``mu_of_Dmv`` to
# closure_comparison() once added.

# =================================================================================================
# Why the DSD shape matters for evaporation (comment block)
# =================================================================================================
# The bulk evaporation rate depends on low-order moments of the DSD: roughly the 1st moment for
# the non-ventilated term and the (b+3)/2 moment (about 1.9 for b ~ 0.8) for the ventilated term.
# Rain mass is the 3rd moment. DSDs with the same rainwater content can therefore evaporate at
# very different rates, with DSDs having more small drops evaporating faster. This motivates the
# closure comparison below (Brandes et al., 2006, https://doi.org/10.1175/JAM2335.1; Snook and
# Xue, 2008, https://doi.org/10.1029/2008GL035866; Kumjian and Ryzhkov, 2010,
# https://doi.org/10.1175/2010JAMC2243.1; see also Dawson et al., 2010,
# https://doi.org/10.1175/2009MWR2956.1, and Morrison and Milbrandt, 2011,
# https://doi.org/10.1175/2010MWR3433.1). Because small drops matter disproportionately and
# Parsivel sampling of the smallest drops is the least reliable part of the spectrum, the binned
# calculation reports per-bin contributions and cumulative fractions.

# =================================================================================================
# Ventilation coefficients
# =================================================================================================
# f_v = a_v + b_v Sc^(1/3) Re^(1/2), with Sc = nu / D_v (Schmidt number) and Re = V D / nu.
VENTILATION_COEFFS = {
    # Ferrier (1994) values, as used in Milbrandt's code and DSDlib.calc_evap.
    'ferrier1994': (0.78, 0.30),  # TODO(Dan): confirm source
    # Commonly used values.
    'common': (0.78, 0.308),  # TODO(Dan): confirm source
}
DEFAULT_VENTILATION = 'ferrier1994'


def _ventilation_coeffs(ventilation):
    """Return (a_v, b_v) from a key of VENTILATION_COEFFS or an (a_v, b_v) pair."""
    if isinstance(ventilation, str):
        return VENTILATION_COEFFS[ventilation]
    a_v, b_v = ventilation
    return float(a_v), float(b_v)


# =================================================================================================
# Transport properties of air
# =================================================================================================
# These linear fits are those in Milbrandt's microphysics code as used in DSDlib.calc_evap.
# TODO(Dan): confirm the source of these fits.

def diffusivity_vapor(T, p):
    """Diffusivity of water vapor in air (m^2 s^-1) given T (K) and p (Pa)."""
    return (2.2157e-5 + 0.0155e-5 * (T - 273.15)) * 1.e5 / p


def thermal_conductivity_air(T):
    """Thermal conductivity of air (J m^-1 s^-1 K^-1) given T (K)."""
    return 2.3971e-2 + 0.0078e-2 * (T - 273.15)


def dynamic_viscosity_air(T):
    """Dynamic viscosity of air (kg m^-1 s^-1) given T (K)."""
    return 1.7153e-5 + 0.005e-5 * (T - 273.15)


def air_density(T, p, RH):
    """Moist air density (kg m^-3) from T (K), p (Pa), and RH (fraction).

    Uses the same thermolib functions (and inputs) as PIPS.calc_thermo.
    """
    qv = thermo.calqv(RH, p, T)
    pt = thermo.caltheta(p, T)
    return thermo.calrho(p, pt, qv)


# =================================================================================================
# Thermodynamic denominator (A + B) and saturation ratio
# =================================================================================================

def thermo_denominator(T, p, rho=None, form='textbook', Lv=thermo.Lv, Rv=thermo.Rv):
    """Heat-conduction (A) and vapor-diffusion (B) terms of the diffusional growth equation.

    The steady-state diffusional growth/evaporation rate of a drop is (e.g. see derivation in
    Rogers and Yau 1996, Chapter 7)

        dm/dt = 2 pi D f_v (S - 1) / (A + B)

    with A and B in m s kg^-1.

    form='textbook':
        A = (L_v / (R_v T) - 1) L_v / (K_a T)
        B = R_v T / (D_v e_s(T))
    form='milbrandt' (the ``ABw`` term of the legacy DSDlib.calc_evap):
        A = L_v^2 / (K_a R_v T^2)
        B = 1 / (rho q_vs D_v)

    Reconciliation of the two forms
    -------------------------------
    - A: the legacy form drops the -L_v/(K_a T) term, i.e., it approximates
      (L_v/(R_v T) - 1) by L_v/(R_v T). Since L_v/(R_v T) ~ 18-20 near the surface, the legacy A
      is ~5-6% larger.
    - B: with q_vs the saturation mixing ratio from thermolib.calqvs and rho_d the dry-air
      density, rho_d q_vs = rho_d eps e_s / (p - e_s) = (e_s / (R_v T)) (p - e) / (p - e_s),
      which equals the saturation vapor density e_s/(R_v T) only at saturation. If rho is the
      moist density (as in pyPIPS), an extra factor (1 + w) appears. The two B terms therefore
      agree exactly only when rho is the dry-air density and e = e_s; they differ by ~1-4% for
      typical subsaturated surface conditions.
    - Combined with the saturation-ratio definition (see saturation_ratio()), the textbook form
      gives E about 0-2.5% larger than the legacy form for the same DSD (5-30 deg C, 50-90% RH).
    Both forms use thermolib.cales/calqvs, so e_s is consistent between the two paths.

    Parameters
    ----------
    T : array_like
        Temperature (K)
    p : array_like
        Pressure (Pa)
    rho : array_like, optional
        Air density (kg m^-3); required for form='milbrandt'
    form : {'textbook', 'milbrandt'}
        Which form to use
    Lv : float
        Latent heat of vaporization (J kg^-1)
    Rv : float
        Gas constant for water vapor (J kg^-1 K^-1)

    Returns
    -------
    (array_like, array_like)
        A and B (m s kg^-1)
    """
    Ka = thermal_conductivity_air(T)
    Dv = diffusivity_vapor(T, p)
    if form == 'textbook':
        es = thermo.cales(p, T)
        A = (Lv / (Rv * T) - 1.) * Lv / (Ka * T)
        B = Rv * T / (Dv * es)
    elif form == 'milbrandt':
        if rho is None:
            raise ValueError("rho is required for form='milbrandt'")
        qvs = thermo.calqvs(p, T)
        A = Lv**2. / (Ka * Rv * T**2.)
        B = 1. / (rho * qvs * Dv)
    else:
        raise ValueError(f"Unknown thermo_denominator form '{form}'")
    return A, B


def saturation_ratio(T, p, RH, definition='vapor_pressure'):
    """Ambient saturation ratio S (fraction; S < 1 is subsaturated).

    definition='vapor_pressure': S = e / e_s = RH (the textbook definition).
    definition='mixing_ratio': S = w / w_s, as in the legacy DSDlib.calc_evap. Since
        w / w_s = (e / e_s) (p - e_s) / (p - e), this is slightly smaller than RH when
        subsaturated (by ~1-4% in 1 - S for typical surface conditions).

    Parameters
    ----------
    T : array_like
        Temperature (K)
    p : array_like
        Pressure (Pa)
    RH : array_like
        Relative humidity (fraction)
    definition : {'vapor_pressure', 'mixing_ratio'}
    """
    if definition == 'vapor_pressure':
        return RH * 1.
    if definition == 'mixing_ratio':
        return thermo.calw(RH, p, T) / thermo.calqvs(p, T)
    raise ValueError(f"Unknown saturation ratio definition '{definition}'")


def evaporation_environment(T, p, RH, rho=None, thermo_form='textbook',
                            saturation_definition='vapor_pressure', Lv=thermo.Lv, Rv=thermo.Rv):
    """Compute the environmental quantities needed by the evaporation rate calculations.

    Parameters
    ----------
    T, p, RH : array_like
        Temperature (K), pressure (Pa), relative humidity (fraction)
    rho : array_like, optional
        Air density (kg m^-3). Computed from T, p, and RH if not given.
    thermo_form : {'textbook', 'milbrandt'}
        See thermo_denominator()
    saturation_definition : {'vapor_pressure', 'mixing_ratio'}
        See saturation_ratio()

    Returns
    -------
    dict
        rho (kg m^-3), S (fraction), A and B (m s kg^-1), AB = A + B, Dv (m^2 s^-1),
        nu (kinematic viscosity, m^2 s^-1), and Sc (Schmidt number).
    """
    if rho is None:
        rho = air_density(T, p, np.minimum(RH, 1.))
    A, B = thermo_denominator(T, p, rho=rho, form=thermo_form, Lv=Lv, Rv=Rv)
    Dv = diffusivity_vapor(T, p)
    nu = dynamic_viscosity_air(T) / rho
    return {
        'rho': rho,
        'S': saturation_ratio(T, p, RH, definition=saturation_definition),
        'A': A,
        'B': B,
        'AB': A + B,
        'Dv': Dv,
        'nu': nu,
        'Sc': nu / Dv,
        'thermo_form': thermo_form,
        'saturation_definition': saturation_definition,
    }


# =================================================================================================
# Fall-speed relations
# =================================================================================================
# NOTE on air-density corrections: existing pyPIPS code mixes conventions for V = V0 (rho0/rho)^e:
#   - PIPS.calc_empirical_fallspeed_rain: rho0 = 1.204 kg m^-3, e = 0.4 (Foote and du Toit 1969)
#   - DSDlib.calc_evap: rho0 = 1.225 kg m^-3, e = 0.5
#   - DSDlib.calc_VQR_ferrier: rho0 = 1.204 kg m^-3, e = 0.5
#   - PIPS.calc_empirical_fallspeed_hail: rho0 = 1.225 kg m^-3, e = 0.5
# The classes below make rho0 and the exponent explicit parameters. Their defaults match the
# existing function each one corresponds to; no existing default has been changed.

class FallSpeed:
    """Base class for terminal fall-speed relations V(D, rho) in SI units.

    Subclasses implement ``__call__(D, rho=None)`` returning V (m s^-1) for diameter D (m) and
    air density rho (kg m^-3). If rho is None, no density correction is applied. D and rho must
    be broadcast-compatible (xarray DataArrays with named dimensions are recommended).
    """

    rho0 = 1.225
    exponent = 0.5
    breakpoints = ()  # Diameters (m) where V(D) is not smooth (used for quadrature)

    def density_factor(self, rho):
        """(rho0 / rho)^exponent, or 1 if rho is None."""
        if rho is None:
            return 1.
        return (self.rho0 / rho)**self.exponent

    def __call__(self, D, rho=None):
        raise NotImplementedError

    def describe(self):
        """Short description for output metadata."""
        return repr(self)


class AtlasFallSpeed(FallSpeed):
    """Empirical rain fall speed of PIPS.calc_empirical_fallspeed_rain (Atlas et al. 1973-type
    relation, via code from David Dowell and Terry Schuur), with an explicit density correction.

    V = 3.78 D^0.67 for D < 3 mm, 9.65 - 10.3 exp(-0.6 D) otherwise (D in mm, V in m s^-1),
    times (rho0 / rho)^exponent. Defaults rho0 = 1.204 kg m^-3 and exponent = 0.4 reproduce the
    Foote and du Toit (1969) correction used by calc_empirical_fallspeed_rain(correct_rho=True).
    """

    breakpoints = (3.0e-3,)

    def __init__(self, rho0=1.204, exponent=0.4):
        self.rho0 = rho0
        self.exponent = exponent

    def __call__(self, D, rho=None):
        v = PIPS.calc_empirical_fallspeed_rain(np.asarray(D) * 1000.)
        if isinstance(D, xr.DataArray):
            v = D.copy(data=v)
        return v * self.density_factor(rho)

    def small_drop_power_law(self):
        """The D < 3 mm branch as an exactly equivalent PowerLawFallSpeed in SI units:
        3.78 (1000 D)^0.67 = 3.78 * 1000^0.67 * D^0.67. Valid for D below self.breakpoints[0]."""
        return PowerLawFallSpeed(a=3.78 * 1000.**0.67, b=0.67, rho0=self.rho0,
                                 exponent=self.exponent)

    def __repr__(self):
        return f"AtlasFallSpeed(rho0={self.rho0}, exponent={self.exponent})"


class PowerLawFallSpeed(FallSpeed):
    """Power-law fall speed V = a D^b (rho0 / rho)^exponent (D in m, V in m s^-1).

    Parameters
    ----------
    a : float
        Coefficient (m^(1-b) s^-1)
    b : float
        Exponent
    rho0 : float
        Reference air density (kg m^-3), default 1.225
    exponent : float
        Density-correction exponent, default 0.5
    """

    def __init__(self, a, b, rho0=1.225, exponent=0.5):
        self.a = a
        self.b = b
        self.rho0 = rho0
        self.exponent = exponent

    def __call__(self, D, rho=None):
        return self.a * D**self.b * self.density_factor(rho)

    def __repr__(self):
        return (f"PowerLawFallSpeed(a={self.a}, b={self.b}, rho0={self.rho0}, "
                f"exponent={self.exponent})")


class FerrierFallSpeed(FallSpeed):
    """Ferrier-type fall speed V = a D^b exp(-f D) (rho0 / rho)^exponent (D in m, V in m s^-1).

    Defaults are the Ferrier (1994) rain constants and density correction used in
    DSDlib.calc_evap: a = 4854 m^(1-b) s^-1, b = 1, f = 195 m^-1, rho0 = 1.225 kg m^-3,
    exponent = 0.5.
    """

    def __init__(self, a=4854., b=1., f=195., rho0=1.225, exponent=0.5):
        self.a = a
        self.b = b
        self.f = f
        self.rho0 = rho0
        self.exponent = exponent

    def __call__(self, D, rho=None):
        return self.a * D**self.b * np.exp(-self.f * D) * self.density_factor(rho)

    def __repr__(self):
        return (f"FerrierFallSpeed(a={self.a}, b={self.b}, f={self.f}, rho0={self.rho0}, "
                f"exponent={self.exponent})")


def ventilation_factor(D, V, nu, Dv, a_v=0.78, b_v=0.30):
    """Ventilation factor f_v = a_v + b_v Sc^(1/3) Re^(1/2).

    Sc = nu / D_v is the Schmidt number (the legacy calc_evap comment erroneously called Sc^(1/3)
    the "Scorer parameter"); Re = V D / nu is the drop Reynolds number.

    Parameters
    ----------
    D : array_like
        Drop diameter (m)
    V : array_like
        Fall speed (m s^-1)
    nu : array_like
        Kinematic viscosity of air (m^2 s^-1)
    Dv : array_like
        Diffusivity of water vapor (m^2 s^-1)
    a_v, b_v : float
        Ventilation coefficients (see VENTILATION_COEFFS)
    """
    return a_v + b_v * (nu / Dv)**(1. / 3.) * np.sqrt(V * D / nu)


# =================================================================================================
# Single-drop evaporation
# =================================================================================================
# The steady-state diffusional form below is a standard textbook relation
# (e.g. Rogers and Yau 1996). Srivastava and Coen (1992; https://doi.org/10.1175/1520-0469(1992)049
# <1643:NEEFTA>2.0.CO;2) give more accurate explicit equations; the Thompson scheme's evaporation
# efficiency follows them (Liu, 2024; https://doi.org/10.3390/rs16224297). Their equations are
# not reproduced here; method='srivastava_coen1992' is a placeholder hook.
# Drops evaporate near the wet-bulb temperature (de Szoeke et al., 2024;
# https://doi.org/10.1029/2024JD041126); the steady-state form accounts for this through A.

_DROP_METHODS = ('steady_diffusion', 'srivastava_coen1992')


def _check_method(method):
    if method == 'srivastava_coen1992':
        # TODO(Dan): implement from Srivastava & Coen (1992)
        raise NotImplementedError("method='srivastava_coen1992' is not implemented yet")
    if method != 'steady_diffusion':
        raise ValueError(f"Unknown method '{method}'; expected one of {_DROP_METHODS}")


def _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition):
    if env is not None:
        return env
    return evaporation_environment(T, p, RH, rho=rho, thermo_form=thermo_form,
                                   saturation_definition=saturation_definition)


def drop_mass_rate(D, T=None, p=None, RH=None, rho=None, fall_speed=None,
                   ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                   saturation_definition='vapor_pressure', method='steady_diffusion', env=None):
    """Mass growth rate of a single drop, dm/dt = 2 pi D f_v (S - 1) / (A + B) (kg s^-1).

    Negative for an evaporating drop.

    Parameters
    ----------
    D : array_like
        Drop diameter (m)
    T, p, RH, rho : array_like
        Temperature (K), pressure (Pa), relative humidity (fraction), air density (kg m^-3).
        Not needed if ``env`` is given.
    fall_speed : FallSpeed, optional
        Fall-speed relation; default AtlasFallSpeed()
    ventilation : str or (float, float)
        Key of VENTILATION_COEFFS or (a_v, b_v)
    thermo_form, saturation_definition : str
        See thermo_denominator() and saturation_ratio()
    method : {'steady_diffusion', 'srivastava_coen1992'}
        Growth equation ('srivastava_coen1992' is not implemented)
    env : dict, optional
        Precomputed output of evaporation_environment()
    """
    _check_method(method)
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    a_v, b_v = _ventilation_coeffs(ventilation)
    V = fall_speed(D, env['rho'])
    fv = ventilation_factor(D, V, env['nu'], env['Dv'], a_v, b_v)
    return 2. * np.pi * D * fv * (env['S'] - 1.) / env['AB']


def drop_diameter_squared_rate(D, T=None, p=None, RH=None, rho=None, fall_speed=None,
                               ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                               saturation_definition='vapor_pressure', method='steady_diffusion',
                               env=None):
    """d(D^2)/dt = 8 f_v (S - 1) / (rho_w (A + B)) (m^2 s^-1).

    Unlike dD/dt, this stays finite as D -> 0, which makes it the preferred integration variable.
    Arguments are as in drop_mass_rate().
    """
    _check_method(method)
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    a_v, b_v = _ventilation_coeffs(ventilation)
    V = fall_speed(D, env['rho'])
    fv = ventilation_factor(D, V, env['nu'], env['Dv'], a_v, b_v)
    return 8. * fv * (env['S'] - 1.) / (rhow * env['AB'])


def drop_diameter_rate(D, *args, **kwargs):
    """dD/dt = 2 (dm/dt) / (rho_w pi D^2) (m s^-1), negative for evaporation.

    Arguments are as in drop_mass_rate(). Diverges as D -> 0 (see drop_diameter_squared_rate()).
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        return drop_diameter_squared_rate(D, *args, **kwargs) / (2. * D)


def drop_diameter_height_rate(D, T=None, p=None, RH=None, rho=None, fall_speed=None,
                              ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                              saturation_definition='vapor_pressure', method='steady_diffusion',
                              env=None):
    """dD/dz = (dD/dt) / (dz/dt) = -(dD/dt) / V(D) (dimensionless), with z positive upward.

    Lagrangian-to-height transformation as in de Szoeke et al. (2024;
    https://doi.org/10.1029/2024JD041126). For a falling, evaporating drop dD/dt < 0, so
    dD/dz > 0: the drop shrinks as it descends. Arguments are as in drop_mass_rate().
    """
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    dDdt = drop_diameter_rate(D, fall_speed=fall_speed, ventilation=ventilation, method=method,
                              env=env)
    return -dDdt / fall_speed(D, env['rho'])


def latent_tendencies(E, rho, Lv=thermo.Lv, cp=thermo.cp):
    """Rain and vapor mixing-ratio tendencies and the latent cooling rate implied by E.

    Parameters
    ----------
    E : array_like
        Evaporation rate (kg m^-3 s^-1, positive = evaporation)
    rho : array_like
        Air density (kg m^-3)

    Returns
    -------
    dict
        dqr_dt = -E/rho and dqv_dt = E/rho (kg kg^-1 s^-1); cooling_rate = -L_v E/(rho c_p)
        (K s^-1, negative = cooling); cooling_rate_K_per_h (K h^-1).
    """
    cooling_rate = -Lv * E / (rho * cp)
    return {
        'dqr_dt': -E / rho,
        'dqv_dt': E / rho,
        'cooling_rate': cooling_rate,
        'cooling_rate_K_per_h': 3600. * cooling_rate,
    }


# =================================================================================================
# Interface with pyPIPS datasets
# =================================================================================================

def _units_of(da):
    return str(da.attrs.get('units', '')).strip()


def thermo_from_parsivel_ds(ds, T_var='fasttemp', p_var='pressure', RH_var='RH_derived'):
    """Extract SI thermodynamic variables from a pyPIPS dataset and recompute air density.

    The caller chooses which temperature (e.g., ``fasttemp`` or ``slowtemp``) and RH (e.g.,
    measured ``RH`` or derived ``RH_derived``) to use. Units are read from the ``units``
    attributes; if missing, pyPIPS native units (deg C, hPa, percent) are assumed with a
    warning. Air density is recomputed from the chosen variables (as in PIPS.calc_thermo), rather
    than read from the stored ``rho``, which may predate thermodynamic QC.

    Returns
    -------
    xarray.Dataset
        T (K), p (Pa), RH (fraction), and rho (kg m^-3), with provenance attributes.
    """
    T = ds[T_var]
    p = ds[p_var]
    RH = ds[RH_var]

    T_units = _units_of(T).lower()
    if 'kelvin' in T_units or T_units == 'k':
        T_K = T * 1.
    else:
        if 'celsius' not in T_units and T_units not in ('degc', 'c', 'deg c'):
            warnings.warn(f"No recognized units for '{T_var}'; assuming degrees Celsius")
        T_K = T + 273.15

    p_units = _units_of(p).lower()
    if p_units in ('pa', 'pascal', 'pascals'):
        p_Pa = p * 1.
    else:
        if not any(s in p_units for s in ('hectopascal', 'hpa', 'mb')):
            warnings.warn(f"No recognized units for '{p_var}'; assuming hPa")
        p_Pa = p * 100.

    RH_units = _units_of(RH).lower()
    if RH_units in ('fraction', '1', 'unitless'):
        RH_frac = RH * 1.
    else:
        if 'percent' not in RH_units and RH_units != '%':
            warnings.warn(f"No recognized units for '{RH_var}'; assuming percent")
        RH_frac = RH / 100.

    rho = air_density(T_K, p_Pa, RH_frac)
    out = xr.Dataset({'T': T_K, 'p': p_Pa, 'RH': RH_frac, 'rho': rho})
    for var in out.data_vars:
        out[var].attrs = {}
    out['T'].attrs.update(units='K', source_variable=T_var)
    out['p'].attrs.update(units='Pa', source_variable=p_var)
    out['RH'].attrs.update(units='1', source_variable=RH_var)
    out['rho'].attrs.update(units='kg m-3',
                            description=f'moist air density recomputed from {T_var}, {p_var}, '
                                        f'and {RH_var} (as in PIPS.calc_thermo)')
    out.attrs.update(T_var=T_var, p_var=p_var, RH_var=RH_var)
    return out


def gamma_params_from_fit(fit_da):
    """Split a ``DSD_MM*``/``DSD_TMM*`` fit DataArray (dims parameter x time) into SI
    (N0, lamda, mu) DataArrays. The fits call the shape parameter ``alpha``."""
    N0 = fit_da.sel(parameter='N0', drop=True)
    lamda = fit_da.sel(parameter='lamda', drop=True)
    mu = fit_da.sel(parameter='alpha', drop=True)
    return N0, lamda, mu


# =================================================================================================
# Gamma-DSD evaporation rates
# =================================================================================================

def _log_or_neg_inf(x):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.log(x)


def _gamma_partial_integral(N0, lamda, s, D_lo=None, D_hi=None):
    """N0 * integral of D^(s-1) exp(-lamda D) dD over [D_lo, D_hi] (or [0, inf) by default),
    = N0 Gamma(s) / lamda^s [P(s, lamda D_hi) - P(s, lamda D_lo)], computed in log space."""
    with np.errstate(over='ignore', invalid='ignore'):
        full = np.exp(_log_or_neg_inf(N0) + gammaln(s) - s * np.log(lamda))
    if D_lo is None and D_hi is None:
        return full
    P_hi = 1. if D_hi is None else gammainc(s, lamda * D_hi)
    P_lo = 0. if D_lo is None else gammainc(s, lamda * D_lo)
    return full * (P_hi - P_lo)


def _closed_form_fallspeed_params(fall_speed, lamda):
    """Return (a, b, lamda_v) for the ventilated term, where lamda_v = lamda + f/2 for a
    Ferrier-type fall speed and lamda for a power law."""
    if isinstance(fall_speed, FerrierFallSpeed):
        return fall_speed.a, fall_speed.b, lamda + 0.5 * fall_speed.f
    if isinstance(fall_speed, PowerLawFallSpeed):
        return fall_speed.a, fall_speed.b, lamda
    raise TypeError("The gamma closed form requires a PowerLawFallSpeed or FerrierFallSpeed; "
                    "use evaporation_rate_gamma_quadrature() for other fall-speed relations.")


def evaporation_rate_gamma(N0, lamda, mu, T=None, p=None, RH=None, rho=None, fall_speed=None,
                           ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                           saturation_definition='vapor_pressure', D_lo=None, D_hi=None,
                           env=None):
    """Closed-form bulk evaporation rate for a gamma DSD N(D) = N0 D^mu exp(-lamda D).

    For a power-law fall speed V = a D^b (rho0/rho)^e,

        E = 2 pi (1 - S) / (A + B) * N0 * [ a_v Gamma(mu+2) / lamda^(mu+2)
              + b_v Sc^(1/3) (a (rho0/rho)^e / nu)^(1/2)
                Gamma(mu + (b+5)/2) / lamda^(mu + (b+5)/2) ]

    For the Ferrier-type fall speed V = a D^b exp(-f D) (rho0/rho)^e, lamda is replaced by
    (lamda + f/2) in the second (ventilated) term, as in the legacy DSDlib.calc_evap. Gamma
    functions are evaluated with gammaln in log space for stability. If D_lo and/or D_hi are
    given, the integral is truncated to [D_lo, D_hi] using regularized incomplete gamma functions.

    Parameters
    ----------
    N0 : array_like
        Intercept parameter (m^-(4+mu))
    lamda : array_like
        Slope parameter (m^-1)
    mu : array_like
        Shape parameter
    T, p, RH, rho : array_like
        Temperature (K), pressure (Pa), relative humidity (fraction), air density (kg m^-3)
    fall_speed : PowerLawFallSpeed or FerrierFallSpeed, optional
        Default FerrierFallSpeed()
    ventilation : str or (float, float)
        Key of VENTILATION_COEFFS or (a_v, b_v)
    thermo_form, saturation_definition : str
        See thermo_denominator() and saturation_ratio()
    D_lo, D_hi : float, optional
        Truncation limits (m)
    env : dict, optional
        Precomputed output of evaporation_environment()

    Returns
    -------
    array_like
        Evaporation rate E (kg m^-3 s^-1, positive = evaporation)
    """
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    if fall_speed is None:
        fall_speed = FerrierFallSpeed()
    a, b, lamda_v = _closed_form_fallspeed_params(fall_speed, lamda)
    a_v, b_v = _ventilation_coeffs(ventilation)
    dens = fall_speed.density_factor(env['rho'])

    term1 = a_v * _gamma_partial_integral(N0, lamda, mu + 2., D_lo, D_hi)
    term2 = (b_v * env['Sc']**(1. / 3.) * np.sqrt(a * dens / env['nu']) *
             _gamma_partial_integral(N0, lamda_v, mu + 0.5 * (b + 5.), D_lo, D_hi))
    return 2. * np.pi * (1. - env['S']) / env['AB'] * (term1 + term2)


def moment_tendency_gamma(k, N0, lamda, mu, T=None, p=None, RH=None, rho=None, fall_speed=None,
                          ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                          saturation_definition='vapor_pressure', env=None):
    """Diffusional tendency of the k-th moment of a gamma DSD (m^k m^-3 s^-1), k > 0.

    dM_k/dt = integral of k D^(k-1) (dD/dt) N(D) dD, with dD/dt = 4 f_v (S-1) / (rho_w (A+B) D):

        dM_k/dt = 4 k (S - 1) / (rho_w (A + B)) * N0 * [ a_v Gamma(mu+k-1) / lamda^(mu+k-1)
                    + b_v Sc^(1/3) (a dens / nu)^(1/2) Gamma(mu+k-1+(b+1)/2) / lamda_v^(...) ]

    For k = 3, -cmr dM3/dt equals evaporation_rate_gamma(). The number tendency (k = 0) is not
    defined by diffusional growth alone (it depends on how drops vanish at D -> 0) and requires
    a closure; see number_tendency(). Arguments are as in evaporation_rate_gamma().
    """
    if k <= 0:
        raise ValueError("moment_tendency_gamma requires k > 0; use number_tendency() for k = 0")
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    if fall_speed is None:
        fall_speed = FerrierFallSpeed()
    a, b, lamda_v = _closed_form_fallspeed_params(fall_speed, lamda)
    a_v, b_v = _ventilation_coeffs(ventilation)
    dens = fall_speed.density_factor(env['rho'])
    term1 = a_v * _gamma_partial_integral(N0, lamda, mu + k - 1.)
    term2 = (b_v * env['Sc']**(1. / 3.) * np.sqrt(a * dens / env['nu']) *
             _gamma_partial_integral(N0, lamda_v, mu + k - 1. + 0.5 * (b + 1.)))
    return 4. * k * (env['S'] - 1.) / (rhow * env['AB']) * (term1 + term2)


def _quadrature_nodes(D_lo, D_hi, breakpoints=(), max_width=0.5e-3, n_gl=16):
    """Composite Gauss-Legendre nodes and weights (m) on [D_lo, D_hi].

    The interval is split at the fall-speed breakpoints and into pieces no wider than
    max_width. On a piece starting at D = 0, the substitution D = D_lo + w u^2 removes the
    D^(mu+1) behavior of the integrand at the origin.
    """
    edges = [D_lo] + [bp for bp in sorted(breakpoints) if D_lo < bp < D_hi] + [D_hi]
    fine_edges = [edges[0]]
    for e0, e1 in zip(edges[:-1], edges[1:]):
        npieces = max(1, int(np.ceil((e1 - e0) / max_width)))
        fine_edges.extend(list(np.linspace(e0, e1, npieces + 1)[1:]))
    x, w = np.polynomial.legendre.leggauss(n_gl)
    nodes, weights = [], []
    for e0, e1 in zip(fine_edges[:-1], fine_edges[1:]):
        u = 0.5 * (x + 1.)  # (0, 1)
        wu = 0.5 * w
        if e0 == 0.:
            nodes.append(e1 * u**2.)
            weights.append(wu * 2. * u * e1)
        else:
            nodes.append(e0 + (e1 - e0) * u)
            weights.append(wu * (e1 - e0))
    return np.concatenate(nodes), np.concatenate(weights)


def evaporation_rate_gamma_quadrature(N0, lamda, mu, T=None, p=None, RH=None, rho=None,
                                      fall_speed=None, ventilation=DEFAULT_VENTILATION,
                                      thermo_form='textbook',
                                      saturation_definition='vapor_pressure', D_lo=0.,
                                      D_hi=26.e-3, max_width=0.5e-3, n_gl=16, env=None,
                                      node_dim='quad_node'):
    """Evaporation rate of a gamma DSD over [D_lo, D_hi] by numerical quadrature.

    Works with any fall-speed callable (e.g., AtlasFallSpeed, which has no gamma closed form).
    Uses composite Gauss-Legendre quadrature (see _quadrature_nodes). The default upper limit is
    the largest Parsivel bin edge (26 mm). Arguments are as in evaporation_rate_gamma(); D_lo and
    D_hi must be scalars (m).

    For AtlasFallSpeed, whose D < 3 mm branch is exactly a power law, the part of the interval
    below 3 mm is evaluated with the truncated closed form instead (exact, and robust for
    strongly negative mu, where the integrand is nearly singular at D = 0). When D_lo = 0 the
    integral diverges for mu <= -2; E is then returned as NaN.

    Returns
    -------
    array_like
        Evaporation rate E (kg m^-3 s^-1, positive = evaporation)
    """
    env = _env_or_compute(env, T, p, RH, rho, thermo_form, saturation_definition)
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    N0 = xr.DataArray(N0) if not isinstance(N0, xr.DataArray) else N0

    E = xr.zeros_like(N0 * lamda * mu, dtype=float)
    D_gl_lo = D_lo
    if hasattr(fall_speed, 'small_drop_power_law') and D_lo < fall_speed.breakpoints[0]:
        D_split = min(D_hi, fall_speed.breakpoints[0])
        E = E + evaporation_rate_gamma(N0, lamda, mu, fall_speed=fall_speed.small_drop_power_law(),
                                       ventilation=ventilation, env=env, D_lo=D_lo, D_hi=D_split)
        D_gl_lo = D_split

    if D_gl_lo < D_hi:
        nodes, weights = _quadrature_nodes(D_gl_lo, D_hi, getattr(fall_speed, 'breakpoints', ()),
                                           max_width=max_width, n_gl=n_gl)
        D = xr.DataArray(nodes, dims=[node_dim])
        W = xr.DataArray(weights, dims=[node_dim])
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            ND = np.exp(_log_or_neg_inf(N0) + mu * np.log(D) - lamda * D)
        dmdt = drop_mass_rate(D, fall_speed=fall_speed, ventilation=ventilation, env=env)
        E = E + (ND * (-dmdt) * W).sum(dim=node_dim, skipna=False)

    if D_lo == 0.:
        E = E.where(mu > -2.)  # Divergent integral at D = 0 for mu <= -2
    if E.ndim == 0:
        return float(E)
    return E


# =================================================================================================
# Binned (observed) DSD evaporation rates
# =================================================================================================

def _bin_geometry_m(ND, diameter_dim):
    """Bin centers, lower edges, and widths in m from pyPIPS ND coordinates (mm)."""
    D = ND['diameter'] / 1000.
    Dl = ND['min_diameter'] / 1000.
    Dr = ND['max_diameter'] / 1000.
    return D, Dl, Dr - Dl


def evaporation_rate_binned(ND, T, p, RH, rho=None, fall_speed=None,
                            ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                            saturation_definition='vapor_pressure', D_detect_min=0.25e-3,
                            extension=None, extension_fallback=False,
                            diameter_dim='diameter_bin', Lv=thermo.Lv, cp=thermo.cp,
                            quad_kwargs=None):
    """Evaporation rate from a binned (observed) DSD, E = sum_i N(D_i) [-dm/dt(D_i)] dD_i.

    The time-matched T, p, and RH must be supplied by the caller (see thermo_from_parsivel_ds()
    for selecting, e.g., fast vs slow temperature or measured vs derived RH). Per-drop rates are
    evaluated at bin centers (midpoint rule, as in DSDlib.calc_moment_bin).

    Sensitivity to small drops: evaporation is weighted toward low-order moments (roughly the
    1st to 1.9th), so small drops matter disproportionately, while Parsivel sampling of the
    smallest drops is the least reliable part of the spectrum. The outputs therefore include
    per-bin contributions (E_bin) and the cumulative fraction of E from drops smaller than each
    bin's upper edge (E_cumfrac). Two estimates are available:
    (a) measured bins only (E), and
    (b) if ``extension`` is given, measured bins with lower edge >= D_cut plus a gamma
        extrapolation over the sub-detectable range [0, D_cut) (E_ext), where D_cut is the
        smallest bin lower edge >= D_detect_min. The gamma is typically a truncated-moment fit
        from DSDlib.fit_DSD_TMM* (stored as DSD_TMM* by calc_MM_fits.py).

    The surface DSD gives the evaporation rate at the surface only; the vertical structure
    requires marching a DSD from aloft (see pyPIPS.rainshaft).

    Parameters
    ----------
    ND : xarray.DataArray
        Binned number density (m^-3 mm^-1), dims (time, diameter_bin) with coordinates
        diameter, min_diameter, and max_diameter (mm)
    T, p, RH : xarray.DataArray
        Temperature (K), pressure (Pa), relative humidity (fraction), on the ND time axis
    rho : xarray.DataArray, optional
        Air density (kg m^-3); computed from T, p, RH if not given
    fall_speed : FallSpeed, optional
        Default AtlasFallSpeed() (PIPS.calc_empirical_fallspeed_rain with its density correction)
    ventilation, thermo_form, saturation_definition : see evaporation_rate_gamma()
    D_detect_min : float
        Lower detection limit (m) for the extension option, default 0.25 mm (Parsivel does not
        report the two smallest bins)
    extension : tuple of array_like, optional
        (N0, lamda, mu) of a gamma DSD in SI units for the sub-detectable extension
    extension_fallback : bool
        If True, at times where the extension is undefined (gamma fit missing or invalid, e.g.,
        mu <= -2, for which the small-drop integral diverges), E_ext falls back to the
        measured-bins-only E, and ``extension_fallback_used`` flags those times. Default False
        (E_ext is NaN there).
    diameter_dim : str
        Name of the diameter dimension
    quad_kwargs : dict, optional
        Passed to evaporation_rate_gamma_quadrature() for the extension

    Returns
    -------
    xarray.Dataset
        E, dqr_dt, cooling_rate (K s^-1), cooling_rate_K_per_h, E_bin, E_cumfrac, S; and if
        extension is given, E_sub, E_ext, frac_E_sub, E_cumfrac_ext, dqr_dt_ext,
        cooling_rate_ext, cooling_rate_ext_K_per_h (and extension_fallback_used if
        extension_fallback is True).
    """
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    env = evaporation_environment(T, p, RH, rho=rho, thermo_form=thermo_form,
                                  saturation_definition=saturation_definition)
    rho = env['rho']
    D, Dl, dD = _bin_geometry_m(ND, diameter_dim)

    dmdt = drop_mass_rate(D, fall_speed=fall_speed, ventilation=ventilation, env=env)
    has_dsd = ND.notnull().any(dim=diameter_dim)
    E_bin = (1000. * ND.fillna(0.)) * (-dmdt) * dD  # m^-4 * kg s^-1 * m
    E = E_bin.sum(dim=diameter_dim, skipna=False).where(has_dsd)
    with np.errstate(divide='ignore', invalid='ignore'):
        E_cumfrac = E_bin.cumsum(dim=diameter_dim) / E.where(E != 0.)

    tend = latent_tendencies(E, rho, Lv=Lv, cp=cp)
    out = xr.Dataset({
        'E': E,
        'dqr_dt': tend['dqr_dt'],
        'cooling_rate': tend['cooling_rate'],
        'cooling_rate_K_per_h': tend['cooling_rate_K_per_h'],
        'E_bin': E_bin.where(has_dsd),
        'E_cumfrac': E_cumfrac,
        'S': env['S'] * xr.ones_like(E),
    })

    if extension is not None:
        N0, lamda, mu = extension
        edges_mm = ND['min_diameter'].to_numpy()
        candidates = edges_mm[edges_mm >= D_detect_min * 1000. - 1.e-9]
        D_cut = float(candidates.min()) / 1000.
        included = ND['min_diameter'] >= D_cut * 1000. - 1.e-9
        qk = {} if quad_kwargs is None else dict(quad_kwargs)
        qk.setdefault('max_width', 0.05e-3)
        E_sub = evaporation_rate_gamma_quadrature(N0, lamda, mu, fall_speed=fall_speed,
                                                  ventilation=ventilation, env=env, D_lo=0.,
                                                  D_hi=D_cut, **qk) * xr.ones_like(E)
        E_meas_incl = E_bin.where(included, 0.).sum(dim=diameter_dim, skipna=False)
        # No drops measured -> no evaporation, even where the gamma fit is undefined
        E_ext = xr.where(E == 0., 0., E_meas_incl + E_sub).where(has_dsd)
        with np.errstate(divide='ignore', invalid='ignore'):
            E_cumfrac_ext = ((E_bin.where(included, 0.).cumsum(dim=diameter_dim) + E_sub) /
                             E_ext.where(E_ext != 0.))
            frac_E_sub = E_sub / E_ext.where(E_ext != 0.)
        if extension_fallback:
            fallback = ~np.isfinite(E_ext) & np.isfinite(E)
            E_ext = xr.where(fallback, E, E_ext)
            E_cumfrac_ext = xr.where(fallback, E_cumfrac, E_cumfrac_ext)
            out['extension_fallback_used'] = fallback.astype(np.int8)
            out['extension_fallback_used'].attrs.update(
                units='1', description='1 where E_ext falls back to the measured-bins-only E '
                                       '(extension undefined)')
        tend_ext = latent_tendencies(E_ext, rho, Lv=Lv, cp=cp)
        out['E_sub'] = E_sub.where(has_dsd)
        out['E_ext'] = E_ext
        out['frac_E_sub'] = frac_E_sub
        out['E_cumfrac_ext'] = E_cumfrac_ext
        out['dqr_dt_ext'] = tend_ext['dqr_dt']
        out['cooling_rate_ext'] = tend_ext['cooling_rate']
        out['cooling_rate_ext_K_per_h'] = tend_ext['cooling_rate_K_per_h']
        out.attrs['D_cut_m'] = D_cut

    _add_evap_attrs(out)
    out.attrs.update({
        'fall_speed': fall_speed.describe(),
        'ventilation': str(_ventilation_coeffs(ventilation)),
        'thermo_form': thermo_form,
        'saturation_definition': saturation_definition,
        'D_detect_min_m': D_detect_min,
        'sign_convention': 'E > 0: liquid lost to vapor; cooling_rate < 0: cooling',
        'note': 'Surface evaporation rate only; cooling rate is a single-process tendency, '
                'not a closed surface heat budget (advection and other terms not included).',
    })
    return out


_EVAP_ATTRS = {
    'E': ('kg m-3 s-1', 'evaporation rate from measured bins (positive = evaporation)'),
    'E_ext': ('kg m-3 s-1', 'evaporation rate from measured bins >= D_cut plus gamma '
                            'extension below D_cut'),
    'E_sub': ('kg m-3 s-1', 'gamma-extension contribution from D < D_cut'),
    'E_bin': ('kg m-3 s-1', 'per-bin contribution to E'),
    'E_cumfrac': ('1', 'cumulative fraction of E from drops smaller than max_diameter'),
    'E_cumfrac_ext': ('1', 'cumulative fraction of E_ext from drops smaller than max_diameter'),
    'frac_E_sub': ('1', 'fraction of E_ext from the gamma extension below D_cut'),
    'dqr_dt': ('kg kg-1 s-1', 'rain mixing-ratio tendency (-E/rho)'),
    'dqr_dt_ext': ('kg kg-1 s-1', 'rain mixing-ratio tendency (-E_ext/rho)'),
    'cooling_rate': ('K s-1', 'latent cooling rate -Lv E/(rho cp) (negative = cooling)'),
    'cooling_rate_K_per_h': ('K h-1', 'latent cooling rate (negative = cooling)'),
    'cooling_rate_ext': ('K s-1', 'latent cooling rate from E_ext (negative = cooling)'),
    'cooling_rate_ext_K_per_h': ('K h-1', 'latent cooling rate from E_ext (negative = cooling)'),
    'S': ('1', 'ambient saturation ratio'),
}


def _add_evap_attrs(ds):
    for var, (units, desc) in _EVAP_ATTRS.items():
        if var in ds:
            ds[var].attrs.update(units=units, description=desc)


# =================================================================================================
# Bulk-scheme closures applied to observed moments
# =================================================================================================

def _gamma_from_q_N(rho, q, Nt, mu):
    """N0 (m^-(4+mu)) and lamda (m^-1) from q (kg/kg), Nt (m^-3), and mu, using DSDlib."""
    lamda = dsd.calc_lamda_gamma(rho, q, Nt, cmr, mu)
    N0 = dsd.calc_N0_gamma(rho, q, Nt, cmr, mu)
    return N0, lamda


def observed_moments(ND, rho, diameter_dim='diameter_bin'):
    """q_r (kg/kg), N_t (m^-3), Z (m^6 m^-3), and the mean volume diameter
    D_mv = (M3/M0)^(1/3) (m) from measured bins (DSDlib.calc_moment_bin).

    These are the moments of the measured (truncated) spectrum; closures built from them see
    only what the instrument sampled.
    """
    M0, _ = dsd.calc_moment_bin(ND.fillna(0.), moment=0, diameter_dim_name=diameter_dim)
    M3, _ = dsd.calc_moment_bin(ND.fillna(0.), moment=3, diameter_dim_name=diameter_dim)
    M6, _ = dsd.calc_moment_bin(ND.fillna(0.), moment=6, diameter_dim_name=diameter_dim)
    with np.errstate(divide='ignore', invalid='ignore'):
        Dmv = (M3 / M0.where(M0 > 0.))**(1. / 3.)
    return {'qr': cmr * M3 / rho, 'Nt': M0, 'Z': M6, 'Dmv': Dmv}


def build_closure_gammas(qr, Nt, Z, rho, N0_1mom=8.e6, mu_1mom=0., mu_fixed=(0., 2.5),
                         use_diag_alpha=True, mu_of_Dmv=None, use_MM036=True,
                         use_solve_alpha=True, fits=None, return_descriptions=False):
    """Gamma DSD parameters that one-, two-, and three-moment bulk schemes would assume, given
    the observed moments.

    Applied to observed moments, each closure is a method-of-moments (MM) fit to the moments a
    scheme predicts (M3 from q_r, M0 = N_t, M6 = Z), with the scheme's shape assumption, and is
    named after its moment set. The fitted gamma reproduces the constraining moments exactly,
    except for MM036_solve_alpha (see below).

    Closures (keys of the returned dict):
    - 'MM3_N0fixed_mu<value>': one-moment; M3 with fixed N0 and mu. Default N0 = 8e6 m^-4 and
      mu = 0, the Marshall-Palmer exponential (Marshall and Palmer 1948).
    - 'MM03_mu<value>': two-moment; M0 and M3 with fixed mu (one per value in mu_fixed).
    - 'MM03_diag_alpha': two-moment; M0 and M3 with mu diagnosed by DSDlib.diag_alpha('qr', D_m),
      where D_m = (M3/M0)^(1/3) is the mean-mass (= mean volume, for spheres) diameter from
      DSDlib.calc_Dm_gamma. Note diag_alpha caps mu at 80.
    - 'MM03_mu_of_Dmv': two-moment; M0 and M3 with mu = mu_of_Dmv(D_mv) (user-supplied, e.g., the
      Seifert (2008) relation once added; see the shape-parameter comment block).
    - 'MM036': three-moment; exact MM fit to M0, M3, and M6 (DSDlib.fit_DSD_MM036).
    - 'MM036_solve_alpha': three-moment; M0 and M3 matched exactly with mu from
      DSDlib.solve_alpha, a piecewise-polynomial approximation of the MM036 relation clipped to
      [0, 40] (as in Milbrandt's three-moment code). M6 is reproduced only approximately.
    - '<name>': (N0, lamda, mu) passed through ``fits``, e.g., the stored MM246 or TMM246 fits.

    Motivation: Brandes et al. (2006; https://doi.org/10.1175/JAM2335.1) found evaporation rates
    from Marshall-Palmer vs constrained-gamma DSDs differ by factors of ~2-3, with the sign of the
    difference depending on regime; see also Dawson et al. (2010;
    https://doi.org/10.1175/2009MWR2956.1), Paukert et al. (2019;
    https://doi.org/10.1029/2018MS001512), and Morrison and Milbrandt (2011;
    https://doi.org/10.1175/2010MWR3433.1).

    Parameters
    ----------
    qr, Nt, Z : array_like
        Rain mixing ratio (kg/kg), total number concentration (m^-3), reflectivity factor M6
        (m^6 m^-3)
    rho : array_like
        Air density (kg m^-3)
    fits : dict, optional
        {name: (N0, lamda, mu)} or {name: fit DataArray with a 'parameter' dim}
    return_descriptions : bool
        If True, also return {closure_name: description}

    Returns
    -------
    dict
        {closure_name: (N0, lamda, mu)} in SI units (and the descriptions if requested)
    """
    closures = {}
    descriptions = {}
    zeros = qr * 0.  # Broadcasts constants for numpy or xarray inputs
    M0 = Nt
    M3 = rho * qr / cmr

    name = f'MM3_N0fixed_mu{mu_1mom:g}'
    Nt_1m = dsd.calc_Nt_gamma(rho, qr, N0_1mom, cmr, mu_1mom)
    lamda_1m = dsd.calc_lamda_gamma(rho, qr, Nt_1m, cmr, mu_1mom)
    closures[name] = (N0_1mom + zeros, lamda_1m, mu_1mom + zeros)
    descriptions[name] = (f'one-moment: M3 with fixed N0 = {N0_1mom:g} m^-4 and mu = '
                          f'{mu_1mom:g}')

    for mu in mu_fixed:
        name = f'MM03_mu{mu:g}'
        N0, lamda = _gamma_from_q_N(rho, qr, Nt, mu)
        closures[name] = (N0, lamda, mu + zeros)
        descriptions[name] = f'two-moment: M0 and M3 with fixed mu = {mu:g}'

    if use_diag_alpha:
        Dm = dsd.calc_Dm_gamma(rho, qr, Nt, cmr)
        mu = dsd.diag_alpha('qr', Dm)
        N0, lamda = _gamma_from_q_N(rho, qr, Nt, mu)
        closures['MM03_diag_alpha'] = (N0, lamda, mu)
        descriptions['MM03_diag_alpha'] = ('two-moment: M0 and M3 with mu = '
                                           'DSDlib.diag_alpha((M3/M0)^(1/3))')

    if mu_of_Dmv is not None:
        Dmv = dsd.calc_Dm_gamma(rho, qr, Nt, cmr)
        mu = mu_of_Dmv(Dmv)
        N0, lamda = _gamma_from_q_N(rho, qr, Nt, mu)
        closures['MM03_mu_of_Dmv'] = (N0, lamda, mu)
        descriptions['MM03_mu_of_Dmv'] = 'two-moment: M0 and M3 with user-supplied mu(D_mv)'

    if use_MM036:
        closures['MM036'] = tuple(dsd.fit_DSD_MM036(M0, M3, Z))
        descriptions['MM036'] = 'three-moment: exact MM fit to M0, M3, and M6'

    if use_solve_alpha:
        mu = dsd.solve_alpha(rho, cmr, qr, Nt, Z)
        N0, lamda = _gamma_from_q_N(rho, qr, Nt, mu)
        closures['MM036_solve_alpha'] = (N0, lamda, mu)
        descriptions['MM036_solve_alpha'] = ('three-moment: M0 and M3 with mu from '
                                             'DSDlib.solve_alpha (approximate MM036, mu in '
                                             '[0, 40])')

    for name, fit in (fits or {}).items():
        if isinstance(fit, xr.DataArray) and 'parameter' in fit.dims:
            fit = gamma_params_from_fit(fit)
        closures[name] = tuple(fit)
        descriptions[name] = f'supplied fit {name}'

    if return_descriptions:
        return closures, descriptions
    return closures


def closure_comparison(ND, T, p, RH, rho=None, binned=None, fall_speed=None,
                       integration='quadrature', closed_form_fall_speed=None,
                       ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                       saturation_definition='vapor_pressure', diameter_dim='diameter_bin',
                       closure_dim='closure', quad_kwargs=None, **closure_kwargs):
    """Compare evaporation rates from bulk-scheme gamma closures with the binned rate.

    Builds the gamma DSDs a bulk scheme would assume from the PIPS-observed moments
    (build_closure_gammas()), computes E for each, and returns the ratio to the binned E.

    integration='quadrature' (default) integrates each gamma over [0, 26 mm] with the SAME
    fall-speed relation as the binned calculation, isolating the effect of the DSD closure.
    integration='closed_form' uses evaporation_rate_gamma() with ``closed_form_fall_speed``
    (default FerrierFallSpeed()), i.e., what a bulk scheme would compute, which also includes
    the difference in fall-speed relation.

    Parameters
    ----------
    ND, T, p, RH, rho : see evaporation_rate_binned()
    binned : xarray.Dataset, optional
        Output of evaporation_rate_binned(); computed if not given. If it contains E_ext, the
        ratio to E_ext is also returned.
    closure_kwargs :
        Passed to build_closure_gammas() (N0_1mom, mu_1mom, mu_fixed, use_diag_alpha,
        mu_of_Dmv, use_MM036, use_solve_alpha, fits)

    Returns
    -------
    xarray.Dataset
        N0, lamda, mu, E, and ratio_to_binned (and ratio_to_binned_ext) with a ``closure`` dim,
        plus the observed moments qr, Nt, Z, Dmv.
    """
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    env = evaporation_environment(T, p, RH, rho=rho, thermo_form=thermo_form,
                                  saturation_definition=saturation_definition)
    rho = env['rho']
    if binned is None:
        binned = evaporation_rate_binned(ND, T, p, RH, rho=rho, fall_speed=fall_speed,
                                         ventilation=ventilation, thermo_form=thermo_form,
                                         saturation_definition=saturation_definition,
                                         diameter_dim=diameter_dim)
    mom = observed_moments(ND, rho, diameter_dim=diameter_dim)
    closures, descriptions = build_closure_gammas(mom['qr'], mom['Nt'], mom['Z'], rho,
                                                  return_descriptions=True, **closure_kwargs)

    names, N0s, lamdas, mus, Es = [], [], [], [], []
    for name, (N0, lamda, mu) in closures.items():
        if integration == 'quadrature':
            E = evaporation_rate_gamma_quadrature(N0, lamda, mu, fall_speed=fall_speed,
                                                  ventilation=ventilation, env=env,
                                                  **(quad_kwargs or {}))
        elif integration == 'closed_form':
            E = evaporation_rate_gamma(N0, lamda, mu, fall_speed=closed_form_fall_speed,
                                       ventilation=ventilation, env=env)
        else:
            raise ValueError(f"Unknown integration '{integration}'")
        names.append(name)
        N0s.append(N0 * xr.ones_like(binned['E']))
        lamdas.append(lamda * xr.ones_like(binned['E']))
        mus.append(mu * xr.ones_like(binned['E']))
        Es.append(E * xr.ones_like(binned['E']))

    index = xr.DataArray(names, dims=[closure_dim], name=closure_dim)
    out = xr.Dataset({
        'N0': xr.concat(N0s, dim=index),
        'lamda': xr.concat(lamdas, dim=index),
        'mu': xr.concat(mus, dim=index),
        'E': xr.concat(Es, dim=index),
    })
    out = out.assign_coords(closure_description=(closure_dim, [descriptions[n] for n in names]))
    has_rain = binned['E'] > 0.
    # Non-finite rates come from invalid gamma parameters (e.g., N0 overflowing to inf for very
    # large mu in a moment fit) and are set to NaN
    out['E'] = out['E'].where(np.isfinite(out['E']))
    # No drops measured -> no evaporation, even where the closure gamma is undefined
    out['E'] = xr.where(binned['E'] == 0., 0., out['E']).where(binned['E'].notnull())
    out['ratio_to_binned'] = (out['E'] / binned['E']).where(has_rain)
    if 'E_ext' in binned:
        out['ratio_to_binned_ext'] = (out['E'] / binned['E_ext']).where(binned['E_ext'] > 0.)
    for key, val in mom.items():
        out[key] = val
    out['N0'].attrs.update(units='m^-(4+mu)')
    out['lamda'].attrs.update(units='m-1')
    out['mu'].attrs.update(units='1')
    out['E'].attrs.update(units='kg m-3 s-1', description='closure evaporation rate')
    out['qr'].attrs.update(units='kg kg-1', description='observed rain mixing ratio')
    out['Nt'].attrs.update(units='m-3', description='observed total number concentration')
    out['Z'].attrs.update(units='m6 m-3', description='observed Rayleigh reflectivity factor')
    out['Dmv'].attrs.update(units='m', description='observed mean volume diameter (M3/M0)^(1/3)')
    out.attrs.update({
        'integration': integration,
        'fall_speed': (fall_speed.describe() if integration == 'quadrature'
                       else (closed_form_fall_speed or FerrierFallSpeed()).describe()),
        'note': 'Closures are built from moments of the measured (truncated) spectrum.',
    })
    return out


# =================================================================================================
# Number-concentration tendency closures
# =================================================================================================
# Needed to compare rainshaft-evolved observed DSDs with what two- and three-moment bulk schemes
# would predict (Seifert, 2008, https://doi.org/10.1175/2008JAS2586.1; Paukert et al., 2019,
# https://doi.org/10.1029/2018MS001512).

NUMBER_TENDENCY_METHODS = ('proportional', 'seifert2008')


def number_tendency(q, N, dq_dt, rho=None, method='proportional'):
    """Rain number-concentration tendency due to evaporation, dN/dt (m^-3 s^-1).

    Methods
    -------
    'proportional' (default):
        dN/N = dq/q, i.e., the mean volume diameter D_mv = (rho q / (cmr N))^(1/3) is held
        constant during evaporation. This is the standard assumption in multi-moment bulk
        schemes.
    'seifert2008':
        Not implemented. TODO(Dan): implement the evaporation number-tendency parameterization
        of Seifert (2008; https://doi.org/10.1175/2008JAS2586.1).

    Parameters
    ----------
    q : array_like
        Rain mixing ratio (kg/kg)
    N : array_like
        Rain number concentration (m^-3)
    dq_dt : array_like
        Rain mixing-ratio tendency from evaporation (kg kg^-1 s^-1, negative for evaporation)
    rho : array_like, optional
        Air density (kg m^-3). Not used by 'proportional'; reserved for 'seifert2008'.
    method : str
        One of NUMBER_TENDENCY_METHODS
    """
    if method == 'seifert2008':
        # TODO(Dan): implement from Seifert (2008)
        raise NotImplementedError("method='seifert2008' is not implemented yet")
    if method not in NUMBER_TENDENCY_METHODS:
        raise ValueError(f"Unknown method '{method}'; expected one of {NUMBER_TENDENCY_METHODS}")
    with np.errstate(divide='ignore', invalid='ignore'):
        return xr.where(q > 0., N * dq_dt / q, 0.)


# =================================================================================================
# Time integration of surface evaporation tendencies
# =================================================================================================
# These tools integrate the rain-evaporation tendencies of temperature and water vapor at the PIPS
# over a chosen time window (default: the whole deployment). They answer the question "how much
# would rain evaporation alone have cooled and moistened the air at the PIPS?" under the
# assumption that the same air stays at the PIPS (no advection, mixing, surface fluxes, radiation,
# or other phase changes). The observed changes contain all of those processes, so differences
# between the observed and integrated changes are NOT errors in the evaporation calculation.
#
# Two variants:
# - integrate_static(): the tendencies computed from the observed thermodynamic state are
#   accumulated. No feedback of the accumulated cooling/moistening on the evaporation rate.
# - integrate_dynamic(): the evaporation rate is recomputed at each step from the observed DSD
#   and the evolving (integrated) T and water vapor, so evaporation slows as the air approaches
#   saturation. Pressure follows the observations.
#
# Each DSD sample is treated as representing the sample interval that ENDS at its (right-labeled)
# time stamp, as produced by the pyPIPS resampling (label='right', closed='right'). The state
# reported at time t_i is the state after integrating sample i. The initial state is the
# observed state at the first valid time in the window.
#
# Each time in the window is classified as rain (E != 0), no rain (E == 0), or missing (E not
# computable: missing DSD or thermodynamic data, invalid gamma parameters, or flagged by
# ``missing``). Missing samples are counted and reported separately from no-rain samples. By
# default (missing_fill='zero') they contribute zero tendency; with missing_fill='last_good' they
# use the last known good sample (its tendency in the static integration, its DSD in the dynamic
# integration), optionally limited to max_fill_samples consecutive samples. Filled samples remain
# classified as missing and are also flagged as filled.

INTEGRATION_STATUS = {'rain': 0, 'no_rain': 1, 'missing': 2}


def _select_window(da, start, end, time_dim='time'):
    return da.sel({time_dim: slice(start, end)})


def _infer_sample_interval(time):
    dt = np.diff(time.to_numpy()).astype('timedelta64[ns]').astype(float) * 1.e-9
    if dt.size == 0:
        raise ValueError("Need at least two times to infer the sample interval")
    return float(np.median(dt))


def mixing_ratio_state(T, p, RH):
    """Water vapor mixing ratio w (kg/kg) from T (K), p (Pa), and RH (fraction)."""
    return thermo.calw(RH, p, T)


def _state_from_T_p_w(T, p, w):
    """RH (fraction), moist and dry air density (kg m^-3) from T (K), p (Pa), w (kg/kg)."""
    e = p * w / (thermo.rddrv + w)
    RH = e / thermo.cales(p, T)
    qv = thermo.w2qv(w)
    rho = thermo.calrho(p, thermo.caltheta(p, T), qv)
    return RH, rho, rho * (1. - qv)


def _status_from_E(E, missing=None):
    status = xr.where(~np.isfinite(E), INTEGRATION_STATUS['missing'],
                      xr.where(E == 0., INTEGRATION_STATUS['no_rain'], INTEGRATION_STATUS['rain']))
    if missing is not None:
        status = xr.where(missing, INTEGRATION_STATUS['missing'], status)
    return status.astype(np.int8)


def _status_summary(status, sample_interval, time_dim='time'):
    """Counts and durations (s) of rain, no-rain, and missing samples (per member if present)."""
    out = {}
    for name, code in INTEGRATION_STATUS.items():
        n = (status == code).sum(dim=time_dim)
        out[f'n_{name}'] = n
        out[f'duration_{name}_s'] = n * sample_interval
    return out


def _ffill_last_good(da, valid, max_fill_samples=None, time_dim='time'):
    """Replace values where ``valid`` is False by the last valid value along time (NaN where
    there is none, or where the gap exceeds max_fill_samples). Returns (filled, was_filled)."""
    da, valid = xr.broadcast(da, valid)
    dims = [d for d in da.dims if d != time_dim] + [time_dim]
    arr = da.transpose(*dims).to_numpy()
    ok = valid.transpose(*dims).to_numpy()
    shape = arr.shape
    arr = arr.reshape(-1, shape[-1])
    ok = ok.reshape(-1, shape[-1])
    nt = shape[-1]
    idx = np.where(ok, np.arange(nt)[None, :], -1)
    idx = np.maximum.accumulate(idx, axis=1)
    has_prev = idx >= 0
    if max_fill_samples is not None:
        has_prev &= (np.arange(nt)[None, :] - idx) <= max_fill_samples
    rows = np.arange(arr.shape[0])[:, None]
    filled = np.where(ok, arr, np.where(has_prev, arr[rows, np.maximum(idx, 0)], np.nan))
    was_filled = ~ok & has_prev
    out = da.transpose(*dims).copy(data=filled.reshape(shape))
    flag = da.transpose(*dims).copy(data=was_filled.reshape(shape).astype(np.int8))
    return out.transpose(*da.dims), flag.transpose(*da.dims)


def _initial_state(T_obs, p_obs, RH_obs, time_dim='time'):
    valid = (T_obs.notnull() & p_obs.notnull() & RH_obs.notnull()).to_numpy()
    if not valid.any():
        raise ValueError("No valid initial thermodynamic state in the integration window")
    i0 = int(np.nonzero(valid)[0][0])
    T0 = float(T_obs.isel({time_dim: i0}))
    p0 = float(p_obs.isel({time_dim: i0}))
    RH0 = float(RH_obs.isel({time_dim: i0}))
    return i0, T0, mixing_ratio_state(T0, p0, RH0)


def _integration_attrs(mode, start, end, sample_interval, time):
    return {
        'integration_mode': mode,
        'integration_start': str(time.to_numpy()[0]) if start is None else str(start),
        'integration_end': str(time.to_numpy()[-1]) if end is None else str(end),
        'sample_interval_s': sample_interval,
        'integration_assumption': 'same air remains at the PIPS; rain evaporation is the only '
                                  'process (no advection, mixing, surface fluxes, radiation, or '
                                  'other phase changes); missing samples contribute zero tendency',
        'status_codes': '0 = rain, 1 = no rain, 2 = missing',
    }


def _check_missing_fill(missing_fill):
    if missing_fill not in ('zero', 'last_good'):
        raise ValueError(f"missing_fill must be 'zero' or 'last_good', not '{missing_fill}'")


def integrate_static(E, T_obs, p_obs, RH_obs, rho=None, start=None, end=None,
                     sample_interval=None, missing=None, missing_fill='zero',
                     max_fill_samples=None, Lv=thermo.Lv, cp=thermo.cp, time_dim='time'):
    """Accumulate the surface rain-evaporation tendencies computed from the observed state.

    dT/dt = -L_v E / (rho c_p) and dw/dt = E / rho_d are evaluated with the observed thermodynamic
    state at each sample and summed over the sample intervals in [start, end]. There is no
    feedback of the accumulated changes on E (see integrate_dynamic()).

    Parameters
    ----------
    E : xarray.DataArray
        Evaporation rate (kg m^-3 s^-1), dims (time,) or (..., time), e.g., (closure, time)
    T_obs, p_obs, RH_obs : xarray.DataArray
        Observed temperature (K), pressure (Pa), RH (fraction) on the same time axis
    rho : xarray.DataArray, optional
        Air density (kg m^-3) used to compute E; computed from the observations if not given
    start, end : str or datetime-like, optional
        Integration window (inclusive); default the whole record
    sample_interval : float, optional
        Duration (s) represented by each sample (e.g., the DSD_interval attribute); default
        the median time step
    missing : xarray.DataArray of bool, optional
        True for missing records (e.g., ``missing_times == 1``)
    missing_fill : {'zero', 'last_good'}
        Tendency used for missing samples: zero (default), or that of the last known good
        (non-missing) sample
    max_fill_samples : int, optional
        With 'last_good', the maximum number of consecutive samples filled (default unlimited)

    Returns
    -------
    xarray.Dataset
        dT and dw (cumulative changes, K and kg/kg), T and w (initial state plus cumulative
        change), RH (predicted, using the observed pressure), dT_obs and dw_obs (observed
        changes since the initial time, for reference), status, filled (1 where a missing
        sample was filled), and per-status counts/durations (plus n_filled).
    """
    _check_missing_fill(missing_fill)
    E = _select_window(E, start, end, time_dim)
    T_obs = _select_window(T_obs, start, end, time_dim)
    p_obs = _select_window(p_obs, start, end, time_dim)
    RH_obs = _select_window(RH_obs, start, end, time_dim)
    if missing is not None:
        missing = _select_window(missing, start, end, time_dim).astype(bool)
    if sample_interval is None:
        sample_interval = _infer_sample_interval(E[time_dim])
    w_obs = mixing_ratio_state(T_obs, p_obs, RH_obs)
    if rho is None:
        rho = air_density(T_obs, p_obs, RH_obs)
    else:
        rho = _select_window(rho, start, end, time_dim)
    rho_d = rho * (1. - thermo.w2qv(w_obs))

    status = _status_from_E(E, missing)
    dTdt = -Lv * E / (rho * cp)
    dwdt = E / rho_d
    valid = (status != INTEGRATION_STATUS['missing']) & np.isfinite(dTdt) & np.isfinite(dwdt)
    if missing_fill == 'last_good':
        dTdt, filled = _ffill_last_good(dTdt, valid, max_fill_samples, time_dim)
        dwdt, _ = _ffill_last_good(dwdt, valid, max_fill_samples, time_dim)
    else:
        dTdt, dwdt = dTdt.where(valid), dwdt.where(valid)
        filled = xr.zeros_like(status)
    dT = dTdt.fillna(0.).cumsum(dim=time_dim) * sample_interval
    dw = dwdt.fillna(0.).cumsum(dim=time_dim) * sample_interval

    _, T0, w0 = _initial_state(T_obs, p_obs, RH_obs, time_dim)
    T = T0 + dT
    w = w0 + dw
    RH, _, _ = _state_from_T_p_w(T, p_obs, w)
    out = xr.Dataset({
        'dT': dT, 'dw': dw, 'T': T, 'w': w, 'RH': RH,
        'dT_obs': T_obs - T0, 'dw_obs': w_obs - w0,
        'status': status,
        'filled': filled.astype(np.int8),
    })
    for key, val in _status_summary(status, sample_interval, time_dim).items():
        out[key] = val
    out['n_filled'] = out['filled'].sum(dim=time_dim)
    _add_integration_var_attrs(out)
    out.attrs.update(_integration_attrs('static', start, end, sample_interval, E[time_dim]))
    out.attrs.update(T0_K=T0, w0_kg_per_kg=w0, missing_fill=missing_fill,
                     max_fill_samples=-1 if max_fill_samples is None else max_fill_samples)
    return out


def _gamma_rate_numpy(N0, lamda, mu, env, fall_speed, ventilation, integration, D_lo, D_hi,
                      gl_nodes):
    """Gamma-DSD evaporation rate on plain numpy arrays (members along axis 0).

    Same calculation as evaporation_rate_gamma_quadrature() / evaporation_rate_gamma(), without
    xarray overhead, for use inside the time-integration loop. ``gl_nodes`` = (nodes, weights,
    D_gl_lo) from _quadrature_setup().
    """
    if integration == 'closed_form':
        return np.asarray(evaporation_rate_gamma(N0, lamda, mu, fall_speed=fall_speed,
                                                 ventilation=ventilation, env=env))
    E = np.zeros(np.shape(N0))
    nodes, weights, D_gl_lo = gl_nodes
    if D_gl_lo > D_lo:  # Exact closed form below the Atlas breakpoint
        E = E + evaporation_rate_gamma(N0, lamda, mu, fall_speed=fall_speed.small_drop_power_law(),
                                       ventilation=ventilation, env=env, D_lo=D_lo, D_hi=D_gl_lo)
    if nodes.size:
        D = nodes[None, :]
        env2 = {k: (v[:, None] if isinstance(v, np.ndarray) and v.ndim == 1 else v)
                for k, v in env.items()}
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            ND = np.exp(_log_or_neg_inf(N0)[:, None] + mu[:, None] * np.log(D) -
                        lamda[:, None] * D)
        dmdt = drop_mass_rate(D, fall_speed=fall_speed, ventilation=ventilation, env=env2)
        E = E + np.sum(ND * -dmdt * weights[None, :], axis=1)
    if D_lo == 0.:
        E = np.where(mu > -2., E, np.nan)  # Divergent integral at D = 0 for mu <= -2
    return E


def _quadrature_setup(fall_speed, D_lo, D_hi, max_width, n_gl):
    """Nodes/weights for the Gauss-Legendre part of the gamma quadrature (see
    evaporation_rate_gamma_quadrature()) and the lower limit of that part."""
    D_gl_lo = D_lo
    if hasattr(fall_speed, 'small_drop_power_law') and D_lo < fall_speed.breakpoints[0]:
        D_gl_lo = min(D_hi, fall_speed.breakpoints[0])
    if D_gl_lo < D_hi:
        nodes, weights = _quadrature_nodes(D_gl_lo, D_hi, getattr(fall_speed, 'breakpoints', ()),
                                           max_width=max_width, n_gl=n_gl)
    else:
        nodes, weights = np.array([]), np.array([])
    return nodes, weights, D_gl_lo


def _state_numpy(T, p, RH, rho):
    return tuple(np.atleast_1d(np.asarray(x, dtype=float)) for x in (T, p, RH, rho))


def binned_rate_function(ND, fall_speed=None, ventilation=DEFAULT_VENTILATION,
                         thermo_form='textbook', saturation_definition='vapor_pressure',
                         extension=None, extension_fallback=False, D_detect_min=0.25e-3,
                         diameter_dim='diameter_bin', time_dim='time'):
    """Build a rate function E(i, T, p, RH, rho) for integrate_dynamic() from a binned DSD.

    The rate function takes the sample index i and numpy arrays T (K), p (Pa), RH (fraction),
    and rho (kg m^-3) of shape (n_members,), and returns E (kg m^-3 s^-1) of the same shape.
    The DSD at sample i is held fixed while the thermodynamic state varies. The calculation is
    that of evaporation_rate_binned() (E, or E_ext with ``extension`` = (N0, lamda, mu) on the
    time axis, including ``extension_fallback``), evaluated with numpy for speed. The rate
    function returns NaN where the DSD at sample i is missing (all bins NaN) and 0 where no drops
    were measured (even if the extension's gamma fit is undefined).
    """
    if fall_speed is None:
        fall_speed = AtlasFallSpeed()
    ND = ND.transpose(time_dim, diameter_dim)
    has_dsd = ND.notnull().any(dim=diameter_dim).to_numpy()
    ND_np = ND.fillna(0.).to_numpy()
    no_drops = np.all(ND_np <= 0., axis=1)
    D, _, dD = _bin_geometry_m(ND, diameter_dim)
    D_np = D.to_numpy()
    weight = 1000. * dD.to_numpy()
    weight_all = weight.copy()
    if extension is not None:
        edges_mm = ND['min_diameter'].to_numpy()
        D_cut = float(edges_mm[edges_mm >= D_detect_min * 1000. - 1.e-9].min()) / 1000.
        weight = np.where(edges_mm >= D_cut * 1000. - 1.e-9, weight, 0.)
        N0_e, lamda_e, mu_e = (np.asarray(x, dtype=float) for x in extension)
        gl = _quadrature_setup(fall_speed, 0., D_cut, 0.05e-3, 16)

    def rate(i, T, p, RH, rho):
        Tn, pn, RHn, rhon = _state_numpy(T, p, RH, rho)
        if not has_dsd[i]:
            return np.full(Tn.shape, np.nan)
        if no_drops[i]:
            return np.zeros(Tn.shape)
        env = evaporation_environment(Tn[:, None], pn[:, None], RHn[:, None], rho=rhon[:, None],
                                      thermo_form=thermo_form,
                                      saturation_definition=saturation_definition)
        dmdt = drop_mass_rate(D_np[None, :], fall_speed=fall_speed, ventilation=ventilation,
                              env=env)
        E = np.sum(ND_np[i] * weight * -dmdt, axis=1)
        if extension is not None:
            env1 = {k: (v[:, 0] if isinstance(v, np.ndarray) else v) for k, v in env.items()}
            ones = np.ones_like(Tn)
            E_sub = _gamma_rate_numpy(N0_e[i] * ones, lamda_e[i] * ones, mu_e[i] * ones, env1,
                                      fall_speed, ventilation, 'quadrature', 0., D_cut, gl)
            if extension_fallback and not np.all(np.isfinite(E_sub)):
                # Extension undefined: fall back to the measured-bins-only rate
                return np.sum(ND_np[i] * weight_all * -dmdt, axis=1)
            E = E + E_sub
        return E

    return rate


def gamma_rate_function(N0, lamda, mu, integration='quadrature', fall_speed=None,
                        ventilation=DEFAULT_VENTILATION, thermo_form='textbook',
                        saturation_definition='vapor_pressure', D_hi=26.e-3, no_rain=None,
                        max_width=0.5e-3, n_gl=16, time_dim='time'):
    """Build a rate function E(i, T, p, RH, rho) for integrate_dynamic() from gamma parameters.

    N0, lamda, and mu (SI) are DataArrays on the time axis, optionally with one additional
    dimension (e.g., ``closure``) for an ensemble of DSDs. In that case the rate function carries
    ``members = (dim, labels)``, which integrate_dynamic() uses to label its output. The rate
    function interface is as in binned_rate_function().
    integration='quadrature' uses ``fall_speed`` (default AtlasFallSpeed, as for the binned
    rate) over [0, D_hi] (as in evaporation_rate_gamma_quadrature()); integration='closed_form'
    uses evaporation_rate_gamma() with ``fall_speed`` (default FerrierFallSpeed). Returns NaN where
    the parameters are invalid, except at times flagged by ``no_rain`` (boolean array on the
    time axis, e.g., measured E == 0), where E = 0 (gamma parameters are undefined when no drops
    were observed).
    """
    if integration not in ('quadrature', 'closed_form'):
        raise ValueError(f"Unknown integration '{integration}'")
    params = []
    for x in (N0, lamda, mu):
        x = xr.DataArray(x) if not isinstance(x, xr.DataArray) else x
        if len([d for d in x.dims if d != time_dim]) > 1:
            raise ValueError("Gamma parameters may have at most one dimension besides time")
        params.append(x)
    ntimes = params[0].sizes[time_dim]
    no_rain_np = (np.zeros(ntimes, dtype=bool) if no_rain is None
                  else np.asarray(no_rain, dtype=bool))
    if fall_speed is None:
        fall_speed = AtlasFallSpeed() if integration == 'quadrature' else FerrierFallSpeed()
    gl = _quadrature_setup(fall_speed, 0., D_hi, max_width, n_gl)

    # Parameters as numpy arrays (time, member)
    template = xr.broadcast(*params)[0]
    extra = [d for d in template.dims if d != time_dim]
    arrays = []
    for x in xr.broadcast(*params):
        if extra:
            arrays.append(x.transpose(time_dim, extra[0]).to_numpy().astype(float))
        else:
            arrays.append(x.to_numpy().astype(float)[:, None])
    N0_np, lamda_np, mu_np = arrays

    def rate(i, T, p, RH, rho):
        Tn, pn, RHn, rhon = _state_numpy(T, p, RH, rho)
        if no_rain_np[i]:
            return np.zeros(Tn.shape)
        shape = np.broadcast_shapes(Tn.shape, N0_np[i].shape)
        n0, lam, m = (np.broadcast_to(a[i], shape) for a in (N0_np, lamda_np, mu_np))
        env = evaporation_environment(Tn, pn, RHn, rho=rhon, thermo_form=thermo_form,
                                      saturation_definition=saturation_definition)
        valid = (n0 > 0.) & np.isfinite(n0) & (lam > 0.) & np.isfinite(lam) & np.isfinite(m)
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            E = _gamma_rate_numpy(np.where(valid, n0, 1.), np.where(valid, lam, 1.e3),
                                  np.where(valid, m, 0.), env, fall_speed, ventilation,
                                  integration, 0., D_hi, gl)
        return np.where(valid, E, np.nan)

    if extra:
        labels = (template[extra[0]].to_numpy() if extra[0] in template.coords
                  else np.arange(template.sizes[extra[0]]))
        rate.members = (extra[0], list(labels))
    return rate


def integrate_dynamic(rate_function, T_obs, p_obs, RH_obs, start=None, end=None,
                      sample_interval=None, missing=None, missing_fill='zero',
                      max_fill_samples=None, Lv=thermo.Lv, cp=thermo.cp, time_dim='time'):
    """Integrate surface rain-evaporation tendencies with feedback on the evaporation rate.

    Starting from the observed state at the first valid time in [start, end], T and the water
    vapor mixing ratio w evolve as dT/dt = -L_v E / (rho c_p) and dw/dt = E / rho_d, where
    E = rate_function(i, T, p, RH, rho) is recomputed from the observed DSD at sample i and the
    evolving (integrated) T, w, and the observed pressure. Each sample interval is integrated with
    classic RK4 with the DSD held fixed. Evaporation therefore slows as the air approaches
    saturation (and condensation occurs if it overshoots).

    Parameters
    ----------
    rate_function : callable
        From binned_rate_function() or gamma_rate_function(). If it has a ``members`` attribute
        (dim, labels), e.g., closures, the integration is done for each member.
    T_obs, p_obs, RH_obs : xarray.DataArray
        Observed temperature (K), pressure (Pa), RH (fraction) on the DSD time axis. Gaps in
        pressure are filled by linear interpolation in time.
    start, end : str or datetime-like, optional
        Integration window (inclusive); default the whole record
    sample_interval : float, optional
        Duration (s) of each sample; default the median time step
    missing : xarray.DataArray of bool, optional
        True for missing records (e.g., ``missing_times == 1``)
    missing_fill : {'zero', 'last_good'}
        For missing samples: zero tendency (default), or E recomputed (with the current state)
        from the DSD of the last known good (non-missing) sample
    max_fill_samples : int, optional
        With 'last_good', the maximum number of consecutive samples filled (default unlimited)

    Returns
    -------
    xarray.Dataset
        As integrate_static(), plus E, cooling_rate, and rho recomputed along the integration
        (E and cooling_rate at the start of each sample interval).
    """
    _check_missing_fill(missing_fill)
    times = T_obs[time_dim]
    i_window = np.nonzero(np.isin(times.to_numpy(),
                                  _select_window(times, start, end, time_dim).to_numpy()))[0]
    T_obs = T_obs.isel({time_dim: i_window})
    RH_obs = RH_obs.isel({time_dim: i_window})
    p_obs = p_obs.isel({time_dim: i_window})
    p_fill = p_obs.interpolate_na(dim=time_dim).ffill(time_dim).bfill(time_dim).to_numpy()
    if missing is not None:
        missing_np = missing.isel({time_dim: i_window}).astype(bool).to_numpy()
    else:
        missing_np = np.zeros(i_window.size, dtype=bool)
    if sample_interval is None:
        sample_interval = _infer_sample_interval(T_obs[time_dim])
    dt = sample_interval

    member_dim, members = getattr(rate_function, 'members', (None, None))
    nmem = 1 if members is None else len(members)

    _, T0, w0 = _initial_state(T_obs, p_obs, RH_obs, time_dim)
    T = np.full(nmem, T0)
    w = np.full(nmem, w0)
    nt = i_window.size
    out = {k: np.full((nmem, nt), np.nan) for k in ('T', 'w', 'RH', 'E', 'cooling_rate', 'rho')}
    status = np.full((nmem, nt), INTEGRATION_STATUS['missing'], dtype=np.int8)
    filled = np.zeros((nmem, nt), dtype=np.int8)
    last_good = np.full(nmem, -1)

    def tendencies(idx, T, w, p):
        """Tendencies with each member's DSD taken from window sample idx (-1: none)."""
        RH, rho, rho_d = _state_from_T_p_w(T, p, w)
        E = np.full(nmem, np.nan)
        for j in np.unique(idx[idx >= 0]):
            sel = idx == j
            Ej = np.broadcast_to(rate_function(i_window[j], T, np.full(nmem, p), RH, rho),
                                 (nmem,))
            E[sel] = Ej[sel]
        return -Lv * E / (rho * cp), E / rho_d, E

    for i in range(nt):
        p = float(p_fill[i])
        idx = np.full(nmem, i)
        dTdt1, dwdt1, E1 = tendencies(idx, T, w, p)
        is_missing = ~np.isfinite(E1) | missing_np[i]
        status[:, i] = np.where(is_missing, INTEGRATION_STATUS['missing'],
                                np.where(E1 == 0., INTEGRATION_STATUS['no_rain'],
                                         INTEGRATION_STATUS['rain']))
        idx = np.where(is_missing, -1, i)
        if missing_fill == 'last_good' and is_missing.any():
            can_fill = is_missing & (last_good >= 0)
            if max_fill_samples is not None:
                can_fill &= (i - last_good) <= max_fill_samples
            if can_fill.any():
                idx = np.where(can_fill, last_good, idx)
                filled[:, i] = can_fill
                dTdt_f, dwdt_f, E_f = tendencies(idx, T, w, p)
                dTdt1 = np.where(can_fill, dTdt_f, dTdt1)
                dwdt1 = np.where(can_fill, dwdt_f, dwdt1)
                E1 = np.where(can_fill, E_f, E1)
        active = idx >= 0
        out['E'][:, i] = np.where(active, E1, np.nan)
        out['cooling_rate'][:, i] = np.where(active, dTdt1, np.nan)
        if active.any():
            dTdt2, dwdt2, _ = tendencies(idx, T + 0.5 * dt * dTdt1, w + 0.5 * dt * dwdt1, p)
            dTdt3, dwdt3, _ = tendencies(idx, T + 0.5 * dt * dTdt2, w + 0.5 * dt * dwdt2, p)
            dTdt4, dwdt4, _ = tendencies(idx, T + dt * dTdt3, w + dt * dwdt3, p)
            dT = dt / 6. * (dTdt1 + 2. * dTdt2 + 2. * dTdt3 + dTdt4)
            dw = dt / 6. * (dwdt1 + 2. * dwdt2 + 2. * dwdt3 + dwdt4)
            T = np.where(active, T + np.nan_to_num(dT), T)
            w = np.where(active, w + np.nan_to_num(dw), w)
        last_good = np.where(is_missing, last_good, i)
        RH, rho, _ = _state_from_T_p_w(T, p, w)
        out['T'][:, i] = T
        out['w'][:, i] = w
        out['RH'][:, i] = RH
        out['rho'][:, i] = rho

    mdim = member_dim or 'member'
    coords = {mdim: members if members is not None else [0],
              time_dim: T_obs[time_dim].to_numpy()}
    dims = [mdim, time_dim]
    ds = xr.Dataset({k: (dims, v) for k, v in out.items()}, coords=coords)
    ds['status'] = (dims, status)
    ds['filled'] = (dims, filled)
    ds['dT'] = ds['T'] - T0
    ds['dw'] = ds['w'] - w0
    w_obs = mixing_ratio_state(T_obs, p_obs, RH_obs).reset_coords(drop=True)
    ds['dT_obs'] = (T_obs.reset_coords(drop=True) - T0)
    ds['dw_obs'] = w_obs - w0
    for key, val in _status_summary(ds['status'], sample_interval, time_dim).items():
        ds[key] = val
    ds['n_filled'] = ds['filled'].sum(dim=time_dim)
    if members is None:
        ds = ds.isel({mdim: 0}, drop=True)
    _add_integration_var_attrs(ds)
    ds.attrs.update(_integration_attrs('dynamic', start, end, sample_interval, ds[time_dim]))
    ds.attrs.update(T0_K=T0, w0_kg_per_kg=w0, pressure='observed (gaps linearly interpolated)',
                    missing_fill=missing_fill,
                    max_fill_samples=-1 if max_fill_samples is None else max_fill_samples)
    return ds


_INTEGRATION_ATTRS = {
    'dT': ('K', 'cumulative temperature change from rain evaporation'),
    'dw': ('kg kg-1', 'cumulative water vapor mixing-ratio change from rain evaporation'),
    'T': ('K', 'initial temperature plus the integrated change'),
    'w': ('kg kg-1', 'initial mixing ratio plus the integrated change'),
    'RH': ('1', 'relative humidity of the integrated state (observed pressure)'),
    'rho': ('kg m-3', 'air density of the integrated state'),
    'E': ('kg m-3 s-1', 'evaporation rate recomputed along the integration'),
    'cooling_rate': ('K s-1', 'latent cooling rate recomputed along the integration'),
    'dT_obs': ('K', 'observed temperature change since the initial time (all processes)'),
    'dw_obs': ('kg kg-1', 'observed mixing-ratio change since the initial time (all processes)'),
    'status': ('1', '0 = rain, 1 = no rain, 2 = missing'),
    'filled': ('1', '1 where a missing sample was filled from the last known good sample'),
}


def _add_integration_var_attrs(ds):
    for var, (units, desc) in _INTEGRATION_ATTRS.items():
        if var in ds:
            ds[var].attrs.update(units=units, description=desc)
    for var in ds.data_vars:
        if var.startswith('duration_'):
            ds[var].attrs.update(units='s')
        elif var.startswith('n_'):
            ds[var].attrs.update(units='1', description='number of samples')
