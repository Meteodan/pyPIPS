"""
DSDlib.py: This is a library of functions to calculate various DSD-related parameters
Some of these were originally written in fortran but re-written using python and numpy
"""

import numpy as np
from scipy.special import gamma as gamma_, gammainc as gammap, gammaln as gammln
from . import thermolib as thermo
from . import PIPS as PIPS
from . import radarmodule as radar
from .utils import first_nonzero, last_nonzero, enable_xarray_wrapper
import xarray as xr
from numba import jit

# Some physical constants
rhol = 1000.  # density of liquid water (kg / m^3)
cmr = np.pi * rhol / 6.  # Constant in mass-diameter relationship for spherical raindrops

# Save some common values of the gamma function for convenience
gamma3 = gamma_(3.)
gamma4 = gamma_(4.)
gamma5 = gamma_(5.)
gamma7 = gamma_(7.)


def calc_Mp(p, lamda, Nt, alpha):
    """Given lamda, Nt, and alpha, computes the pth moment of the gamma distribution"""
    return (Nt / lamda**p) * (gamma_(1. + alpha + p) / gamma_(1. + alpha))


def calc_Dmpq(p, q, lamda, Nt, alpha):
    """Given two moments, p and q, and lamda, Nt, and alpha of the gamma distribution, compute the
       mean diameter associated with the ratio of p to q

    Parameters
    ----------
    p : array_like
        The first moment
    q : array_like
        The second moment
    lamda : array_like
        Slope parameter of the gamma distribution
    Nt : array_like
        Total number concentration
    alpha : array_like
        Shape parameter of the gamma distribution

    Returns
    -------
    array_like
        The appropriate moment-weighted mean diameter
    """

    M_p = calc_Mp(p, lamda, Nt, alpha)
    M_q = calc_Mp(q, lamda, Nt, alpha)

    return (M_p / M_q)**(1. / (p - q))


def calc_qr_gamma(rhoa, N0, lamda, alpha):
    """Computes the rain mass mixing ratio given the standard gamma distribution parameters. All
    units are SI unless otherwise specified.

    Parameters
    ----------
    rhoa : array_like
        Air density
    N0 : array_like
        Intercept parameter
    lamda : array_like
        Slope parameter
    alpha : array_like
        Shape parameter

    Returns
    -------
    array_like
        Rain water mass mixing ratio
    """
    GR2 = gamma_(4. + alpha)
    return (cmr / rhoa) * N0 * GR2 / lamda**(alpha + 4.)


def calc_Zr_gamma(rhoa, q, Nt, alpha):
    """Computes the rain radar reflectivity factor given q, Nt, and alpha. All units are SI unless
    otherwise specified.

    Parameters
    ----------
    rhoa : array_like
        Air density
    q : array_like
        rain water mass mixing ratio
    Nt : array_like
        Total number concentration
    alpha : array_like
        Shape parameter

    Returns
    -------
    array_like
        The rain radar reflectivity factor
    """
    Gr = ((6. + alpha) * (5. + alpha) * (4. + alpha)) / \
        ((3. + alpha) * (2. + alpha) * (1. + alpha))
    Zr = ((1. / cmr)**2.) * Gr * ((rhoa * q)**2.) / Nt
    return 10.0 * np.log10(1.e18 * Zr)


def calc_Zr_lin_gamma(rhoa, q, Nt, alpha):
    """Computes the rain radar reflectivity factor given q, Nt, and alpha. All units are SI unless
    otherwise specified. Returns linear units

    Parameters
    ----------
    rhoa : array_like
        Air density
    q : array_like
        rain water mass mixing ratio
    Nt : array_like
        Total number concentration
    alpha : array_like
        Shape parameter

    Returns
    -------
    array_like
        The rain radar reflectivity factor
    """
    Gr = ((6. + alpha) * (5. + alpha) * (4. + alpha)) / \
        ((3. + alpha) * (2. + alpha) * (1. + alpha))
    Zr = ((1. / cmr)**2.) * Gr * ((rhoa * q)**2.) / Nt
    return Zr


def calc_Nt_gamma(rhoa, q, N0, cx, alpha):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates total number concentration for a gamma distribution
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !
    !  03/31/08 - converted intermediate calculations to double precision
    !             as well as a few of the input arguments.
    !
    !  04/24/09 - rewritten in Python/Numpy
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    gamma1alp = gamma_(1.0 + alpha)
    gamma4alp = gamma_(4.0 + alpha)

    Ntx = (N0 * gamma1alp)**(3.0 / (4.0 + alpha)) * \
        ((gamma1alp / gamma4alp) * rhoa * q / cx)**((1.0 + alpha) / (4.0 + alpha))

    return Ntx


@enable_xarray_wrapper
def calc_lamda_gamma(rhoa, q, Ntx, cx, alpha):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates slope parameter lamda for a gamma distribution
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Converted intermediate calculations and arrays alpha and lamda to
    !  double precision.
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !"""

    gamma1alp = gamma_(1.0 + alpha)
    gamma4alp = gamma_(4.0 + alpha)

    lamda = ((gamma4alp / gamma1alp) * cx * Ntx / (rhoa * q))**(1.0 / 3.0)
    lamda = np.where(rhoa * q > 0.0, lamda, 0.0)

    return lamda


@enable_xarray_wrapper
def calc_N0_gamma(rhoa, q, Ntx, cx, alpha):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates intercept parameter for a gamma distribution
    !
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !
    !  (03/26/2008)
    !  Recast N0 as a double precision variable, and used double precision for
    !  all intermediate calculations.  The calling subroutine should
    !  also define it as double precision.  For situations with large alpha,
    !  N0 can become very large, and loss of precision can result.
    !  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
    !  With Jason Milbrandt's calculation of N0 just before evaporation in
    !  the multi-moment code.
    !
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    gamma1alp = gamma_(1.0 + alpha)

    lamda = calc_lamda_gamma(rhoa, q, Ntx, cx, alpha)
    lamda = lamda.astype(np.float64)
    try:
        alpha = alpha.astype(np.float64)
    except Exception:
        pass

    N0 = Ntx * lamda**(0.50 * (1.0 + alpha)) * \
        (1.0 / gamma1alp) * lamda**(0.50 * (1.0 + alpha))
    N0 = np.where(lamda >= 0.0, N0, 0.0)
    return N0


@enable_xarray_wrapper
def calc_N0_norm_gamma(N0, alpha, lamda):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates normalized intercept parameter for a gamma distribution
    !            after Testud et al. (2001)
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !
    !  (03/26/2008)
    !  Recast N0 as a double precision variable, and used double precision for
    !  all intermediate calculations.  The calling subroutine should
    !  also define it as double precision.  For situations with large alpha,
    !  N0 can become very large, and loss of precision can result.
    !  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
    !  With Jason Milbrandt's calculation of N0 just before evaporation in
    !  the multi-moment code.
    !
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """
    gamma4alp = gamma_(4.0 + alpha)

    N0_norm = N0 * (((4.0 + alpha) / lamda) **
                    alpha) * gamma4alp * (128.0 / 3.0) / \
        ((4.0 + alpha)**(4.0 + alpha))
    N0_norm = np.where(lamda >= 0.0, N0_norm, 0.0)
    return N0_norm


@enable_xarray_wrapper
def calc_Dm_gamma(rhoa, q, Ntx, cx):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates mean-mass diameter
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    Dm = (rhoa * q / (cx * Ntx))**(1. / 3.)
    Dm = np.where(Ntx > 0.0, Dm, 0.0)

    return Dm


@enable_xarray_wrapper
def calc_D0_gamma(rhoa, q, Ntx, cx, alpha):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates median volume diameter (m) for a gamma distribution
    !            Note, assumes spherical particles.
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    # Compute lamda
    lamda = calc_lamda_gamma(rhoa, q, Ntx, cx, alpha)

    return (3.672 + alpha) / lamda


def calc_moment_bin(ND, moment=0):
    """Compute the desired moment from number density bins.

    Parameters
    ----------
    ND : array_like
        Number density as a function of diameter
    moment : int, optional
        The desired moment, by default 0

    Returns
    -------
    (array_like, array_like)
        The total moment and the binned moment
    """
    avg_diameter_m = ND['diameter'] / 1000.  # Convert mm to m
    bin_width_m = (ND['max_diameter'] - ND['min_diameter']) / 1000.  # Convert mm to m

    moment_binned = avg_diameter_m**moment * (1000. * ND) * bin_width_m  # binned moments
    moment = moment_binned.sum(dim='diameter_bin')  # Sum over bins to get final moment

    return moment, moment_binned


# TODO: this isn't right. Use old function for now. It works.
# def calc_D0_bin_new(ND):

#     Dm = ND['diameter'] / 1000.  # Convert mm to m
#     M3, M3_binned = calc_moment_bin(ND, moment=3)

#     M3_med = M3.quantile(0.5)
#     # Cumulative sum of M3 with increasing bin size
#     M3_cumsum = M3_binned.cumsum(dim='diameter_bin')
#     # Proportion of M3 in each bin
#     pro = M3_binned / M3
#     # Cumulative proportion of M3 in each bin
#     pro_cumsum = M3_cumsum / M3
#     # Multiply cumulative proportion by the midpoint diameter of that bin
#     pro_Dm = pro * Dm
#     print(pro_Dm.loc[dict(time='2016-03-31T22:30:00')])
#     D0 = pro_Dm.quantile(0.5, dim='diameter_bin')
#     print(D0.loc[dict(time='2016-03-31T22:30:00')])
#     return D0


def calc_D0_bin(ND):
    """Calculate D0 for a binned distribution"""
    Dl = ND['min_diameter'] / 1000.  # Convert mm to m
    Dm = ND['diameter'] / 1000.  # Convert mm to m
    Dr = ND['max_diameter'] / 1000.  # Convert mm to m
    M3, M3_binned = calc_moment_bin(ND, moment=3)

    # Cumulative sum of M3 with increasing bin size
    M3_cumsum = M3_binned.cumsum(dim='diameter_bin')
    # Proportion of M3 in each bin
    pro = M3_binned / M3
    # Cumulative proportion of M3 in each bin
    pro_cumsum = M3_cumsum / M3
    # Compute median volume diameter using a linear interpolation within the "mass-midpoint" bin.
    # Source: http://www.dropletmeasurement.com/PADS_Help/MVD_(um).htm
    # (originally from FAA Electronic Aircraft Icing Handbook)

    # Find indices of bin where cumulative sum exceeds 1/2 of total and the indices
    # of the bin just below that, thresholding on the smallest bin
    medindices = (pro_cumsum > 0.5).argmax(dim='diameter_bin')
    medindices_m1 = medindices - 1
    medindices_m1 = medindices_m1.where(medindices_m1 >= 0, other=0)
    # medindices_m1[medindices_m1 < 0] = 0
    b1 = Dl[medindices]  # Lower boundaries of mass-midpoint bin
    b2 = Dr[medindices]  # Upper boundaries of mass-midpoint bin
    # print(pro)
    # print(medindices)
    # pro_med = pro.loc[dict(diameter_bin=medindices)]
    # pro_cumsum_med_m1 = pro_cumsum.loc[dict(diameter_bin=medindices_m1)]
    pro_med = pro[medindices, :]
    pro_cumsum_med_m1 = pro_cumsum[medindices_m1, :]

    # Now we can calculate D0 by linearly interpolating diameters within
    # a given bin bounded by Dl and Dr which contains the half-mass point
    D0 = b1 + ((0.5 - pro_cumsum_med_m1) / pro_med) * (b2 - b1)
    # Don't let D0 be any smaller than the midpoint of the smallest bin
    D0 = D0.where(D0 >= Dm[0], other=Dm[0])
    # Set to NaN whereever there is no DSD
    D0 = D0.where(M3 > 0.)
    # Finally remove coordinates that we don't need (there's some issue with xarray where
    # their dimensions are reset to an incorrect one anyway. i.e. diameter(diameter_bin) becomes
    # diameter(time) for some reason)
    # TODO 05/12/2020. After above modifications to use where() instead of fancy indexing
    # this step doesn't appear to be needed anymore. Suspicious...
    # D0 = D0.reset_coords(names=['diameter', 'min_diameter', 'max_diameter'], drop=True)
    return D0


@enable_xarray_wrapper
def diag_alpha(varid_qscalar, Dm):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates shape parameter alpha
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Changed alpha array to double precision.
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    alphaMAX = 80.0

    c1 = np.array([19.0, 12.0, 4.5, 5.5, 3.7])
    c2 = np.array([0.6, 0.7, 0.5, 0.7, 0.3])
    c3 = np.array([1.8, 1.7, 5.0, 4.5, 9.0])
    c4 = np.array([17.0, 11.0, 5.5, 8.5, 6.5])

    if varid_qscalar == 'qr':
        nq = 0
    elif varid_qscalar == 'qi':
        nq = 1
    elif varid_qscalar == 'qs':
        nq = 2
    elif varid_qscalar == 'qg':
        nq = 3
    elif varid_qscalar == 'qh':
        nq = 4

    alpha = c1[nq] * np.tanh(c2[nq] * (1.e3 * Dm - c3[nq])) + c4[nq]
    if nq == 4:
        alpha = np.where(Dm > 0.008, 1.e3 * Dm - 2.6, alpha)

    alpha = np.minimum(alpha, alphaMAX)

    return alpha


@enable_xarray_wrapper
def solve_alpha(rhoa, cx, q, Ntx, Z):
    """!
    !-----------------------------------------------------------------------
    !  PURPOSE:  Calculates shape parameter alpha
    !-----------------------------------------------------------------------
    !
    !  AUTHOR: Dan Dawson
    !  (02/06/2008)
    !
    !  MODIFICATION HISTORY:
    !  (03/31/2008)
    !  Changed alpha array to double precision
    !-----------------------------------------------------------------------
    !  Variable Declarations:
    !-----------------------------------------------------------------------
    !     """

    epsQ = 1.e-14
    epsN = 1.e-3
    epsZ = 1.e-32

    alphaMax = 40.

    tmp1 = cx / (rhoa * q)
    g = tmp1 * Z * tmp1 * Ntx
    g = np.where((q > epsQ) & (Ntx > epsN) & (Z > epsZ), g, -99.0)

    a = np.empty_like(q)
    a = np.where(g == -99.0, 0.0, a)
    a = np.where(g >= 20.0, 0.0, a)
    g2 = g * g

    a = np.where((g < 20.0) & (g >= 13.31), 3.3638e-3 * g2 - 1.7152e-1 * g + 2.0857e+0, a)
    a = np.where((g < 13.31) & (g >= 7.123), 1.5900e-2 * g2 - 4.8202e-1 * g + 4.0108e+0, a)
    a = np.where((g < 7.123) & (g >= 4.200), 1.0730e-1 * g2 - 1.7481e+0 * g + 8.4246e+0, a)
    a = np.where((g < 4.200) & (g >= 2.946), 5.9070e-1 * g2 - 5.7918e+0 * g + 1.6919e+1, a)
    a = np.where((g < 2.946) & (g >= 1.793), 4.3966e+0 * g2 - 2.6659e+1 * g + 4.5477e+1, a)
    a = np.where((g < 1.793) & (g >= 1.405), 4.7552e+1 * g2 - 1.7958e+2 * g + 1.8126e+2, a)
    a = np.where((g < 1.405) & (g >= 1.230), 3.0889e+2 * g2 - 9.0854e+2 * g + 6.8995e+2, a)
    a = np.where(g < 1.230, alphaMax, a)

    alpha = np.maximum(0., np.minimum(a, alphaMax))

    return alpha


def calc_evap(rho, T, p, RH, N0, lamda, mu):
    """Computes bulk evaporation rate given the thermodynamic state of the air and gamma
       distribution parameters.  Much of this code was adapted from Jason Milbrandt's microphysics
       code.  Note, all thermodynamic variables are assumed to be in SI units. TODO: refactor this
       !"""

    # Thermodynamic constants

    Rv = thermo.Rv
    cp = thermo.cp
    Lv = thermo.Lv

    # Correction for fallspeeds at different pressures
    rho_ref = 1.225
    gamma_factor = np.sqrt(rho_ref / rho)

    # Calculate a bunch of quantities in the ventilation coefficient for the bulk evaporation rate
    # Ventilation coefficients (Ferrier 1994)
    Avx = 0.78
    Bvx = 0.30
    afr = 4854.0
    # constants in fallspeed relation for rain from Ferrier 1994
    bfr = 1.0
    ffr = 195.0

    Cdiff = (2.2157e-5 + 0.0155e-5 * (T - 273.15)) * 1e5 / p

    # Dynamic viscosity
    MUdyn = 1.72e-5 * (393. / (T + 120.)) * (T / 273.16)**1.5
    # Kinematic viscosity
    MUkin = MUdyn / rho
    # Scorer parameter to the one-third power
    ScTHRD = (MUkin / Cdiff)**(1. / 3.)

    GR16 = gamma_(2. + mu)
    GR17 = gamma_(2.5 + mu + 0.5 * bfr)

    cexr5 = 2. + mu
    cexr6 = 2.5 + mu + 0.5 * bfr

    # Now we can calculate the bulk ventilation coefficient for rain (Woohoo!)
    VENTr = Avx * GR16 / (lamda**cexr5) + Bvx * ScTHRD * \
        np.sqrt(gamma_factor * afr / MUkin) * GR17 / (lamda + ffr)**cexr6

    # Calculate the thermodynamic function in the denominator of the bulk evaporation rate
    # Thermal conductivity of air
    Ka = 2.3971e-2 + 0.0078e-2 * (T - 273.15)

    # Calculate saturation mixing ratio
    QSS = thermo.calqvs(p, T)

    ABw = (Lv**2.) / (Ka * Rv * T**2.) + 1. / (rho * QSS * Cdiff)

    # Calculate saturation deficit S
    RH_frac = RH / 100.0

    # Calculate water vapor mixing ratio
    qv = thermo.calqv(RH_frac, p, T)
    S = qv / QSS - 1.0

    # With the above quantities, we can calculate the bulk evaporation rate (Yay!)
    QREVP = 2. * np.pi * S * N0 * VENTr / ABw

    # With QREVP we can calculate the latent cooling rate
    COOL_RATE = (Lv / cp) * QREVP

    return QREVP, COOL_RATE


def calc_VQR_ferrier(rhoa, q, Ntx, cx, alpha):
    """Calculates the mass-weighted rain fall speed using the formula in Ferrier (1994),
       assuming a gamma-diameter distribution and given q, N, and alpha"""

    # Ferrier (1994)
    afr = 4854.0
    bfr = 1.0
    ffr = 195.0

    # First, compute slope parameter
    lamda = calc_lamda_gamma(rhoa, q, Ntx, cx, alpha)
    # Compute density correction factor
    gamfact = (1.204 / rhoa)**0.5

    ckQr1 = afr * gamma_(4.0 + alpha + bfr) / gamma_(4.0 + alpha)
    cexr1 = 4.0 + alpha + bfr
    cexr2 = 4.00 + alpha

    VQR = gamfact * ckQr1 * lamda**(0.5 * cexr2) / (lamda + ffr)**(0.5 * cexr1) * \
        lamda**(0.5 * cexr2) / (lamda + ffr)**(0.5 * cexr1)

    return VQR


def calc_VQR_gamvol(rhoa, q, Nt, nu):
    """Calculates the mass-weighted rain fall speed for the NSSL scheme gamma-volume rain DSD"""

    # Compute density correction factor
    gamfact = (1.204 / rhoa)**0.5
    # Compute mean volume
    vr = rhoa * q / (1000. * Nt)

    VQR = gamfact * (0.0911229 * (1 + nu)**1.3333333333333333 * gamma_(2. + nu) +
                     5430.313059683277 * (1 + nu) * vr**0.3333333333333333 *
                     gamma_(2.333333333333333 + nu) -
                     1.0732802065650471e6 * (1 + nu)**0.6666666666666666 * vr**0.6666666666666666 *
                     gamma_(2.6666666666666667 + nu) +
                     8.584110982429507e7 * (1 + nu)**0.3333333333333333 * vr * gamma_(3 + nu) -
                     2.3303765697228556e9 * vr**1.3333333333333333 *
                     gamma_(3.333333333333333 + nu)) /  \
                    ((1 + nu)**2.333333333333333 * gamma_(1 + nu))

    return VQR


def calc_VQR_gamdiam(rhoa, q, Nt, cx, alpha):
    """Calculates the mass-weighted rain fall speed for the NSSL scheme gamma-diameter rain DSD"""

    arx = 10.
    # raind fit parameters for arx*(1 - Exp(-fx*d)), where d is rain diameter in meters.
    frx = 516.575

    # First, compute slope parameter
    lamda = calc_lamda_gamma(rhoa, q, Nt, cx, alpha)
    Dn = 1. / lamda
    # Compute density correction factor
    gamfact = (1.204 / rhoa)**0.5

    VQR = gamfact * arx * (1.0 - (1.0 + frx * Dn)**(-alpha - 4.0))

    return VQR


def calc_VQG(rhoa, ax, bx, q, Ntx, cx, alpha):
    """Calculates mass-weighted graupel fall speed"""

    lamda = calc_lamda_gamma(rhoa, q, Ntx, cx, alpha)
    Dn = 1. / lamda
    gamfact = (1.204 / rhoa)**0.5
    VQG = gamfact * ax * (gamma_(4.0 + alpha + bx) / gamma_(4.0 + alpha)) * Dn**0.5

    return VQG


def calc_VNG(rhoa, ax, bx, q, Ntx, cx, alpha):
    """Calculates number-weighted graupel fall speed"""

    lamda = calc_lamda_gamma(rhoa, q, Ntx, cx, alpha)
    Dn = 1. / lamda
    gamfact = (1.204 / rhoa)**0.5
    VNG = gamfact * ax * (gamma_(1.0 + alpha + bx) / gamma_(1.0 + alpha)) * Dn**0.5

    return VNG


def power_mom(power, cx, t, q, moment):
    """!
    !-----------------------------------------------------------------------
    !
    ! PURPOSE:
    !
    ! Calculates moments of the PSD based on the Field et al. 2005 power law
    ! relations. Used for Thompson scheme.
    !
    ! AUTHOR:  Bryan Putnam, 4/16/2013
    !
    !-----------------------------------------------------------------------
    !"""

    T_c = t - 273.16

    second_moment = (q / cx)

    if power == 1:
        moment = second_moment
    else:
        log_a = 5.065339 - .062659 * T_c - 3.032362 * power + 0.029469 * T_c * power - \
            0.000285 * (T_c**2) + 0.312550 * (power**2) + 0.0000204 * (T_c**2) * power + \
            0.003199 * T_c * (power**2) + 0.000000 * (T_c**3) - 0.015952 * (power**3)

        a = 10.**log_a

        b = 0.476221 - 0.015896 * T_c + 0.165977 * power + 0.007468 * T_c * power -   \
            0.000141 * (T_c**2) + 0.060366 * (power**2) + 0.000079 * (T_c**2) * power + \
            0.000594 * T_c * (power**2) + 0.000000 * (T_c**3) - 0.003577 * (power**3)

    moment = a * (second_moment)**b

    return moment


def calc_gamma_DSD(rhoa, D, cx, q, Nt=None, N0=None, alpha=0):
    """Given cx, q, Nt or N0, and alpha, compute the gamma DSD for the sequence of diameters in D"""
    if N0 is None:
        try:
            N0, _ = calc_N0_gamma(rhoa, q / 1000., Nt, cx, alpha)
        except Exception:
            return None
    else:
        try:
            Nt = calc_Nt_gamma(rhoa, q / 1000., N0, cx, alpha)
        except Exception:
            return None

    lamda = calc_lamda_gamma(rhoa, q / 1000., Nt, cx, alpha)

    return N0 * D**alpha * np.exp(-lamda * D), Nt, lamda, alpha


def calc_NT_from_bins(ND):
    """Computes total number concentration from binned number density ND.

    Parameters
    ----------
    ND : array_like
        binned number density

    Returns
    -------
    array_like
        Total number concentration
    """
    Nt, _ = calc_moment_bin(ND, moment=0)
    return Nt


def calc_lwc_qr_from_bins(ND, rho):
    """Computes liquid water content and rain mass mixing ratio from binned number density

    Parameters
    ----------
    ND : array_like
        binned number density
    rho : array_like
        air density

    Returns
    -------
    (array_like, array_like)
        liquid water content and rain mass mixing ratio
    """
    M3, _ = calc_moment_bin(ND, moment=3)
    LWC = cmr * M3
    qr = LWC / rho

    return LWC, qr


def calc_dBZ_from_bins(ND):
    """Computes radar reflectivity (dBZ) from binned number density

    Parameters
    ----------
    ND : array_like
        binned number density

    Returns
    -------
    array_like
        logarithmic radar reflectivity (dBZ)
    """
    M6, _ = calc_moment_bin(ND, moment=6)
    return 10. * np.log10(1.e18 * M6)


# def calc_synthetic_bins(D_min, D_max, num_bins):
#     """TODO: deprecate this

#     Parameters
#     ----------
#     D_min : [type]
#         [description]
#     D_max : [type]
#         [description]
#     num_bins : [type]
#         [description]

#     Returns
#     -------
#     [type]
#         [description]
#     """
#     return np.linspace(D_min, D_max, num=num_bins)


def fit_DSD_MM24(M2, M4):
    """Uses the Method-of-Moments to fit an exponential distribution using M2 and M4

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M4 : array_like
        The fourth moment of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution. In this case
        the shape parameter is always zero.
    """
    # lamda = np.where(M4 == 0.0, 0.0, ((M2 * gamma5) / (M4 * gamma3))**(1./2.))
    lamda = ((M2 * gamma5) / (M4 * gamma3))**(1./2.)
    N0 = (M2 * lamda**3.) / gamma3
    mu = lamda.copy(data=np.zeros_like(lamda.values))
    return N0, lamda, mu


def fit_DSD_MM36(M3, M6):
    """Uses the Method-of-Moments to fit an exponential distribution using M3 and M6

    Parameters
    ----------
    M3 : array_like
        The third moment of the distribution
    M6 : array_like
        The sixth moment of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution. In this case
        the shape parameter is always zero.
    """
    # lamda = np.where(M6 == 0.0, 0.0, ((M3 * gamma7) / (M6 * gamma4))
    #                  ** (1. / 3.))
    lamda = ((M3 * gamma7) / (M6 * gamma4))** (1. / 3.)
    N0 = (M3 * lamda**4.) / gamma4
    mu = lamda.copy(data=np.zeros_like(lamda.values))
    return N0, lamda, mu


def fit_DSD_MM346(M3, M4, M6):
    """Uses the Method-of-Moments to fit a gamma distribution using M3, M4, and M6

    Parameters
    ----------
    M3 : array_like
        The third moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = (M4**3.) / ((M3**2.) * M6)
    # G = np.ma.masked_invalid(G)
    mu = (11. * G - 8. + (G * (G + 8.))**(1. / 2.)) / (2. * (1. - G))
    # mu = np.ma.masked_invalid(mu)
    lamda = (M3 * (mu + 4.)) / M4
    # lamda = np.ma.masked_invalid(lamda)
    N0 = (M3 * lamda**(mu + 4.)) / (gamma_(mu + 4.))

    return N0, lamda, mu


def fit_DSD_MM246_old(M2, M4, M6):
    """Uses the Method-of-Moments to fit a gamma distribution using M2, M4, and M6. Older method?

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = np.where((M2 == 0.0) | (M6 == 0.0), 0.0, (M4**2.) / (M2*M6))
    mu = np.where(G == 1.0, 0.0, ((7. - 11. * G) -
                                  ((7. - 11. * G)**2. -
                                   4. * (G - 1.) * (30. * G - 12.))**(1. / 2.)) /
                  (2. * (G - 1.)))
    mu = np.where(mu <= -4., -3.99, mu)
    mu = np.where(mu > 30., 30., mu)
    mu = np.ma.masked_where(M4 is np.ma.masked, mu)
    lamda = np.where(M4 == 0.0, 0.0, ((M2 * (mu + 3.) * (mu + 4.)) / (M4))**(1. / 2.))
    N0 = (M4*lamda**(mu + 5.))/(gamma_(mu + 5.))

    return N0, lamda, mu


def fit_DSD_MM246(M2, M4, M6, lamda_limit=20000., mu_limit=30.):
    """Uses the Method-of-Moments to fit a gamma distribution using M2, M4, and M6.

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution
    lamda_limit : array_like, optional
        Upper limit of lamda, by default 20000.
    mu_limit : array_like, optional
        Upper limit of mu, by default 30.

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = (M4**2.) / (M2 * M6)
    # G = np.ma.masked_invalid(G)
    mu = ((7. - 11. * G) - ((7. - 11. * G)**2. - 4. * (G - 1.)
                            * (30. * G - 12.))**(1. / 2.)) / (2. * (G - 1.))
    # TODO: should we just set mu and lambda greater than the limit to their limits instead of
    # masking them?
    # mu = np.ma.masked_where(mu > mu_limit, mu)
    # mu = np.ma.masked_invalid(mu)
    lamda = ((M2 * (mu + 3.) * (mu + 4.)) / (M4))**(1. / 2.)
    # mu = np.ma.masked_where((lamda > lamda_limit) | (mu > 30.), mu)
    # mu = np.ma.masked_invalid(mu)
    # lamda = np.ma.masked_where(lamda > lamda_limit, lamda)
    # lamda = np.ma.masked_invalid(lamda)
    N0 = (M4 * lamda**(mu + 5.)) / (gamma_(mu + 5.))

    return N0, lamda, mu


def fit_DSD_MM234(M2, M3, M4):
    """Uses the Method-of-Moments to fit a gamma distribution using M2, M3, and M4.

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M3 : array_like
        The third moment of the distribution
    M4 : array_like
        The fourth moment of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    mu = (3. * M2 * M4 - 4. * M3**2.) / (M3**2. - M2 * M4)
    # mu = np.ma.masked_invalid(mu)
    lamda = (M3 * (mu + 4.)) / M4
    # lamda = np.ma.masked_invalid(lamda)
    N0 = (M3 * lamda**(mu + 4.)) / (gamma_(mu + 4.))

    return N0, lamda, mu


def fit_DSD_MMXYZ(moment_combo, moment_list):
    """Given a string of the form 'XY(Z)' and a list [X, Y, (Z)], where Z is optional,
       use the Method-of-Moments to fit an exponential or gamma distribution using the combination
       of moments given by X, Y, and possibly Z.
    Parameters
    ----------
    moment_combo: str
        a string with the moment combination
    moment_list: list
        A list of the individual moments to fit
    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution. If the third
        moment (Z) is not provided, the shape parameter will be zero.
    """

    num_moments = len(moment_combo)

    if num_moments < 2 or num_moments > 3:
        print("Incorrect number of moments. Must be 2 or 3!")
        return
    try:
        X = moment_list[0]
        Y = moment_list[1]
    except IndexError:
        print("Not enough moments in list!")
        return

    if num_moments == 2:
        if moment_combo == '24':
            return fit_DSD_MM24(X, Y)
        elif moment_combo == '36':
            return fit_DSD_MM36(X, Y)
        else:
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return
    else:
        try:
            Z = moment_list[2]
        except IndexError:
            print("Not enough moments in list!")
            return
        if moment_combo == '234':
            return fit_DSD_MM234(X, Y, Z)
        elif moment_combo == '246':
            return fit_DSD_MM246(X, Y, Z)   # Todo allow for lambda limits in here through kwargs
        elif moment_combo == '346':
            return fit_DSD_MM346(X, Y, Z)
        else:
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return


def get_max_min_diameters(ND, dim='time'):
    """Gets the maximum and minimum diameters of an array of binned distributions

    Parameters
    ----------
    ND : array_like
        binned number density
    dim : str
        name of the dimension along which ND is changing, default 'time'

    Returns
    -------
    2-tuple of array_like
        The minimum and maximum diameters of an array of binned number density distributions
    """

    # This is surprisingly more of a PITA than I thought it would be.... For each time,
    # need to select the appropriate diameter, so we need to broadcast the diameter array
    # along the time dimension, and then index using two DataArrays, both dimensioned by time
    # The first is just the index of each time, the second is the index of the diameter we
    # want (for each time). It's clunky, but it works. Not sure there is a better way.
    ntimes = ND.sizes[dim]
    tindices = xr.DataArray(range(ntimes), dims=dim)
    D_min_indices = xr.DataArray(first_nonzero(ND, 1), dims=dim)
    D_max_indices = xr.DataArray(last_nonzero(ND, 1), dims=dim)
    D_min = ND['diameter'].expand_dims({dim: ntimes})[tindices, D_min_indices] / 1000.
    D_max = ND['diameter'].expand_dims({dim: ntimes})[tindices, D_max_indices] / 1000.

    return D_min, D_max


def fit_DSD_TMM246_xr(M2, M4, M6, D_min, D_max):
    """Fits gamma distributions using the Truncated Method of Moments (TMM) using M2, M4, and M6
    and D_min, and D_max. Wrapper function for fit_DSD_TMM which operates on numpy arrays and is
    sped up using numba.

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution
    D_min : array_like
        The minimum diameter of the distribution
    D_max : array_like
        The maximum diameter of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    M2_arr = M2.values
    M4_arr = M4.values
    M6_arr = M6.values
    D_min_arr = D_min.values
    D_max_arr = D_max.values

    N0, lamda, alpha = fit_DSD_TMM246(M2_arr, M4_arr, M6_arr, D_min_arr, D_max_arr)

    N0_da = M2.copy(data=N0)
    lamda_da = M2.copy(data=lamda)
    alpha_da = M2.copy(data=alpha)

    return N0_da, lamda_da, alpha_da


# @jit(parallel=True)
@jit
def fit_DSD_TMM246(M2, M4, M6, D_min, D_max):
    """Fits gamma distributions using the Truncated Method of Moments (TMM) using M2, M4, and M6
    and D_min, and D_max. Operates on numpy arrays and uses numba jit to speed things up
    considerably

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution
    D_min : array_like
        The minimum diameter of the distribution
    D_max : array_like
        The maximum diameter of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = (M4**2.) / (M2 * M6)  # moment ratio based on untruncated moments
    # G = np.ma.masked_invalid(G)

    # Get initial estimate of lamda and mu from the untruncated fit:
    _, lamda_init, mu_init = fit_DSD_MM246(M2, M4, M6)

    # TODO: make sure everything works with xarray and see if there's an efficient way to
    # vectorize this algorithm across time. Maybe wrap Guifu's original fortran code is the
    # best approach.
    numtimes = np.size(M2)
    # print(numtimes)
    # print(lamda_init.size)
    # print(D_max.size)

    mu_tmf = []
    lamda_tmf = []
    # N0_tmf = []

    for t in range(numtimes):
        # print("Working on time {:d}".format(t))
        if M2[t] > 0. and M4[t] > 0. and M6[t] > 0.:
            # print("D_max = ", D_max[t])
            print("Working on time {:d}/{:d}".format(t, numtimes))
            LDmx = lamda_init[t] * D_max[t]
            LDmn = lamda_init[t] * D_min[t]
            for _ in range(10):
                mu = mu_init[t]
                # truncated moment ratio below. Equation A8 from Thurai
                gm3 = (gammap(3. + mu, LDmx) - gammap(3. + mu, LDmn)) * np.exp(gammln(3. + mu))
                gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                gm7 = (gammap(7. + mu, LDmx) - gammap(7. + mu, LDmn)) * np.exp(gammln(7. + mu))
                z0 = G[t] - gm5**2. / (gm3 * gm7)
                z1 = z0

                while z1 / z0 > 0.0:
                    mu = mu - 0.01
                    gm3 = (gammap(3. + mu, LDmx) - gammap(3. + mu, LDmn)) * np.exp(gammln(3. + mu))
                    gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                    gm7 = (gammap(7. + mu, LDmx) - gammap(7. + mu, LDmn)) * np.exp(gammln(7. + mu))
                    z1 = G[t] - gm5**2. / (gm3 * gm7)

                lam_tmf = ((M2[t] * gm5) / (M4[t] * gm3))**0.5
                LDmx = lam_tmf * D_max[t]
                LDmn = lam_tmf * D_min[t]
        else:
            mu = np.nan
            lam_tmf = np.nan
        mu_tmf.append(mu)
        lamda_tmf.append(lam_tmf)

    mu_tmf = np.array(mu_tmf)
    lamda_tmf = np.array(lamda_tmf)
    LDmx = lamda_tmf * D_max
    LDmn = lamda_tmf * D_min
    N0_tmf = ((M4 * lamda_tmf**(5. + mu_tmf)) /
              ((gammap(5. + mu_tmf, LDmx) - gammap(5. + mu_tmf, LDmn)) *
               np.exp(gammln(5. + mu_tmf))))

    return N0_tmf, lamda_tmf, mu_tmf


@jit
def fit_DSD_TMM234(M2, M3, M4, D_min, D_max):
    """Fits gamma distributions using the Truncated Method of Moments (TMM) using M2, M3, and M4
    and D_min, and D_max. Operates on numpy arrays and uses numba jit to speed things up
    considerably

    Parameters
    ----------
    M2 : array_like
        The second moment of the distribution
    M3 : array_like
        The third moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    D_min : array_like
        The minimum diameter of the distribution
    D_max : array_like
        The maximum diameter of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = (M3**2.) / (M2 * M4)  # moment ratio based on untruncated moments
    # G = np.ma.masked_invalid(G)

    # Get initial estimate of lamda and mu from the untruncated fit:
    _, lamda_init, mu_init = fit_DSD_MM234(M2, M3, M4)

    # TODO: make sure everything works with xarray and see if there's an efficient way to
    # vectorize this algorithm across time. Maybe wrap Guifu's original fortran code is the
    # best approach.
    numtimes = np.size(M2)
    # print(numtimes)
    # print(lamda_init.size)
    # print(D_max.size)

    mu_tmf = []
    lamda_tmf = []
    # N0_tmf = []

    for t in range(numtimes):
        # print("Working on time {:d}".format(t))
        if M2[t] > 0. and M3[t] > 0. and M4[t] > 0.:
            # print("D_max = ", D_max[t])
            print("Working on time {:d}/{:d}".format(t, numtimes))
            LDmx = lamda_init[t] * D_max[t]
            LDmn = lamda_init[t] * D_min[t]
            for _ in range(10):
                mu = mu_init[t]
                # truncated moment ratio below.
                gm3 = (gammap(3. + mu, LDmx) - gammap(3. + mu, LDmn)) * np.exp(gammln(3. + mu))
                gm4 = (gammap(4. + mu, LDmx) - gammap(4. + mu, LDmn)) * np.exp(gammln(4. + mu))
                gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                z0 = G[t] - gm4**2. / (gm3 * gm5)
                z1 = z0

                while z1 / z0 > 0.0:
                    mu = mu - 0.01
                    gm3 = (gammap(3. + mu, LDmx) - gammap(3. + mu, LDmn)) * np.exp(gammln(3. + mu))
                    gm4 = (gammap(4. + mu, LDmx) - gammap(4. + mu, LDmn)) * np.exp(gammln(4. + mu))
                    gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                    z1 = G[t] - gm4**2. / (gm3 * gm5)

                lam_tmf = (M2[t] * gm4) / (M3[t] * gm3)
                LDmx = lam_tmf * D_max[t]
                LDmn = lam_tmf * D_min[t]
        else:
            mu = np.nan
            lam_tmf = np.nan
        mu_tmf.append(mu)
        lamda_tmf.append(lam_tmf)

    mu_tmf = np.array(mu_tmf)
    lamda_tmf = np.array(lamda_tmf)
    LDmx = lamda_tmf * D_max
    LDmn = lamda_tmf * D_min
    # TODO: should we use M3 for N0?
    N0_tmf = ((M4 * lamda_tmf**(5. + mu_tmf)) /
              ((gammap(5. + mu_tmf, LDmx) - gammap(5. + mu_tmf, LDmn)) *
               np.exp(gammln(5. + mu_tmf))))

    return N0_tmf, lamda_tmf, mu_tmf


def fit_DSD_TMM346(M3, M4, M6, D_min, D_max):
    """Fits gamma distributions using the Truncated Method of Moments (TMM) using M3, M4, and M6
    and D_min, and D_max. Operates on numpy arrays and uses numba jit to speed things up
    considerably

    Parameters
    ----------
    M3 : array_like
        The third moment of the distribution
    M4 : array_like
        The fourth moment of the distribution
    M6 : array_like
        The sixth moment of the distribution
    D_min : array_like
        The minimum diameter of the distribution
    D_max : array_like
        The maximum diameter of the distribution

    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution.
    """
    G = (M4**3.) / (M3**2. * M6)  # moment ratio based on untruncated moments
    # G = np.ma.masked_invalid(G)

    # Get initial estimate of lamda and mu from the untruncated fit:
    _, lamda_init, mu_init = fit_DSD_MM346(M3, M4, M6)

    # TODO: make sure everything works with xarray and see if there's an efficient way to
    # vectorize this algorithm across time. Maybe wrap Guifu's original fortran code is the
    # best approach.
    numtimes = np.size(M3)
    # print(numtimes)
    # print(lamda_init.size)
    # print(D_max.size)

    mu_tmf = []
    lamda_tmf = []
    # N0_tmf = []

    for t in range(numtimes):
        # print("Working on time {:d}".format(t))
        if M3[t] > 0. and M4[t] > 0. and M6[t] > 0.:
            # print("D_max = ", D_max[t])
            print("Working on time {:d}/{:d}".format(t, numtimes))
            LDmx = lamda_init[t] * D_max[t]
            LDmn = lamda_init[t] * D_min[t]
            for _ in range(10):
                mu = mu_init[t]
                # truncated moment ratio below.
                gm4 = (gammap(4. + mu, LDmx) - gammap(4. + mu, LDmn)) * np.exp(gammln(4. + mu))
                gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                gm7 = (gammap(7. + mu, LDmx) - gammap(7. + mu, LDmn)) * np.exp(gammln(7. + mu))
                z0 = G[t] - gm5**3. / (gm4**2. * gm7)
                z1 = z0

                while z1 / z0 > 0.0:
                    mu = mu - 0.01
                    gm4 = (gammap(4. + mu, LDmx) - gammap(4. + mu, LDmn)) * np.exp(gammln(4. + mu))
                    gm5 = (gammap(5. + mu, LDmx) - gammap(5. + mu, LDmn)) * np.exp(gammln(5. + mu))
                    gm7 = (gammap(7. + mu, LDmx) - gammap(7. + mu, LDmn)) * np.exp(gammln(7. + mu))
                    z1 = G[t] - gm5**3. / (gm4**2. * gm7)

                lam_tmf = (M3[t] * gm5) / (M4[t] * gm4)
                LDmx = lam_tmf * D_max[t]
                LDmn = lam_tmf * D_min[t]
        else:
            mu = np.nan
            lam_tmf = np.nan
        mu_tmf.append(mu)
        lamda_tmf.append(lam_tmf)

    mu_tmf = np.array(mu_tmf)
    lamda_tmf = np.array(lamda_tmf)
    LDmx = lamda_tmf * D_max
    LDmn = lamda_tmf * D_min
    # TODO: should we use a different moment for N0? Like M3?
    N0_tmf = ((M4 * lamda_tmf**(5. + mu_tmf)) /
              ((gammap(5. + mu_tmf, LDmx) - gammap(5. + mu_tmf, LDmn)) *
               np.exp(gammln(5. + mu_tmf))))

    return N0_tmf, lamda_tmf, mu_tmf


def fit_DSD_TMMXYZ(moment_combo, moment_list, D_min, D_max):
    """Given a string of the form 'XY(Z)' and a list [X, Y, (Z)], where Z is optional,
       use the Truncated Method-of-Moments to fit an exponential or gamma distribution using the
       combination of moments given by X, Y, and possibly Z.
    Parameters
    ----------
    moment_combo: str
        a string with the moment combination
    moment_list: list
        A list of the individual moments to fit
    Returns
    -------
    3-tuple of array_like
        The intercept, slope, and shape parameters of the fitted gamma distribution. If the third
        moment (Z) is not provided, the shape parameter will be zero.
    """

    # TODO: make default D_min = 0?
    D_min_arr = D_min.values
    D_max_arr = D_max.values

    num_moments = len(moment_combo)

    if num_moments < 2 or num_moments > 3:
        print("Incorrect number of moments. Must be 2 or 3!")
        return
    try:
        X = moment_list[0]
        Y = moment_list[1]
        X_arr = X.values
        Y_arr = Y.values
    except IndexError:
        print("Not enough moments in list!")
        return

    if num_moments == 2:
        if moment_combo == '24':
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return
            # return fit_DSD_TMM24(X, Y)
        elif moment_combo == '36':
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return
            # return fit_DSD_MM36(X, Y)
        else:
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return
    else:
        try:
            Z = moment_list[2]
            Z_arr = Z.values
        except IndexError:
            print("Not enough moments in list!")
            return

        if moment_combo == '234':
            N0, lamda, alpha = fit_DSD_TMM234(X_arr, Y_arr, Z_arr, D_min_arr, D_max_arr)
        elif moment_combo == '246':
            # TODO allow for lambda limits in here through kwargs
            N0, lamda, alpha = fit_DSD_TMM246(X_arr, Y_arr, Z_arr, D_min_arr, D_max_arr)
        elif moment_combo == '346':
            N0, lamda, alpha = fit_DSD_TMM346(X_arr, Y_arr, Z_arr, D_min_arr, D_max_arr)
        else:
            print("Sorry, {} not implemented yet!".format(moment_combo))
            return

        N0_da = X.copy(data=N0)
        lamda_da = X.copy(data=lamda)
        alpha_da = X.copy(data=alpha)

        return N0_da, lamda_da, alpha_da


def calc_binned_DSD_from_params(N0, lamda, alpha, D):
    """Discretizes the gamma distribution given N0, lamda, and alpha into the bins defined by
    the array of diameters D

    Parameters
    ----------
    N0 : array_like
        Intercept parameter
    lamda : array_like
        Slope parameter
    alpha : array_like
        Shape parameter
    D : array_like
        array of diameters

    Returns
    -------
    array_like
        The binned distribution.
    """

    # Get D to m
    D = D / 1000.
    # return 1.e-3 * N0 * D**alpha * np.exp(-lamda * D)  # Units are m^-3 mm^-1
    return N0 * D**alpha * np.exp(-lamda * D)  # Units are m^-4


def calc_rain_axis_ratio(D, fit_name='Brandes_2002'):
    """Computes the axis ratio of rain as a function of diameter using one of a number of
    approximations and fits.

    Parameters
    ----------
    D : array_like
        array of diameters
    fit_name : str, optional
        The type of fit, by default 'Brandes_2002'

    Returns
    -------
    array_like
        The array of axis ratios
    """
    if fit_name == 'Brandes_2002':
        ar = 0.9951 + 0.0251*D - 0.03644*D**2. + 0.005303*D**3. - 0.0002492*D**4.
    elif fit_name == 'Green_1975_exact':
        ar = 1.0148 - 0.020465*D - 0.020048*D**2. + 3.095e-3*D**3. - 1.453e-4*D**4.
    elif fit_name == 'experimental':
        ar = 1.0162 + 0.009714*D - 0.033*D**2. + 5.0424e-3*D**3. - 2.458e-4*D**4.
    elif fit_name == 'experimental_no_beard':
        ar = 1.0 + 0.01367*D - 0.03299*D**2. + 4.853e-3*D**3. - 2.26e-4*D**4.
    elif fit_name == 'Beard_Chuang_1987':
        ar = 1.0048 + 5.7e-4*D - 2.628e-2*D**2. + 3.682e-3*D**3. - 1.627e-4*D**4.

    return ar


def calc_Dmpq_binned(p, q, ND):
    """Computes the requested moment-weighted mean diameter

    Parameters
    ----------
    p : int
        Moment in numerator
    q : int
        Moment in denominator
    ND : array_like
        binned DSD

    Returns
    -------
    array_like
        The appropriate moment-weighted mean diameter
    """
    Mp, _ = calc_moment_bin(ND, moment=p)
    Mq, _ = calc_moment_bin(ND, moment=q)

    return (Mp / Mq)**(1. / (p - q))


# TODO: change all calculations everywhere to use SI units and only change them for plots!
def calc_sigma(D, dD, ND):
    """Computes the width of the mass distribution (from a binned number distribution)

    Parameters
    ----------
    D : array_like
        central values of diameter bins (mm)
    ND : array_like
        The binned number distribution (m^-3 mm^-1)

    Returns
    -------
    array_like
        The width of the mass distribution
    """
    D = D / 1000.  # Get to m
    dD = dD / 1000.  # Get to m
    M3, _ = calc_moment_bin(ND, moment=3)
    Dm = calc_Dmpq_binned(4, 3, ND)
    ND = ND * 1000.  # Get to m^-4
    sigma_numerator = np.sum((D - Dm)**2. * ND * D**3. * dD, axis=0)

    return np.sqrt(sigma_numerator / M3)


def calc_rainrate_from_bins(ND, correct_rho=False, rho=None):
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]
    correct_rho : bool, optional
        [description], by default False
    rho : [type], optional
        [description], by default None

    Returns
    -------
    [type]
        [description]
    """

    avg_diameter_mm = ND['diameter']
    bin_width_mm = (ND['max_diameter'] - ND['min_diameter'])
    fallspeed = PIPS.calc_empirical_fallspeed(avg_diameter_mm, correct_rho=correct_rho, rho=rho)
    rainrate_bin = (6. * 10.**-4.) * np.pi * fallspeed * avg_diameter_mm**3. * ND * bin_width_mm
    rainrate = rainrate_bin.sum(dim='diameter_bin')
    return rainrate


def calc_empirical_polyfit(var_x, var_y, order=3):
    """Compute a 3rd-order polynomial fit to two variables

    Parameters
    ----------
    var_x : array_like
        first variable
    var_y : array_like
        second variable

    Returns
    -------
    tuple (np.ndarray, np.polynomial.polynomial.Polynomial)
        The coefficients of the polynomial along with the Polynomial class object
    """

    poly_coeff = np.polynomial.polynomial.polyfit(var_x, var_y, order)
    poly = np.polynomial.polynomial.Polynomial(poly_coeff)
    return poly_coeff, poly


def calc_CG_polynomial(lamda, mu):
    """Computes a least-squares quadratic polynomial fit to the mu-lamda points

    Parameters
    ----------
    lamda : array_like
        Array of lamda values
    mu : array_like
        Corresponding array of mu values

    Returns
    -------
    tuple (np.ndarray, np.polynomial.polynomial.Polynomial)
        The coefficients of the polynomial along with the Polynomial class object
    """
    CG_poly_coeff = np.polynomial.polynomial.polyfit(lamda, mu, 2)
    CG_poly = np.polynomial.polynomial.Polynomial(CG_poly_coeff)
    return CG_poly_coeff, CG_poly


def calc_mu_lamda(lamda, coefficients):
    """Calculates mu from lamda given the C-G coefficients

    Parameters
    ----------
    lamda : array_like
        array of lamda values
    coefficients : 3-tuple of floats
        coefficients of the polynomial (as output of numpy.polynomial.polynomial.polyfit)

    Returns
    -------
    array_like
        array of mu values
    """
    polynomial = np.polynomial.polynomial.Polynomial(coefficients)
    return polynomial(lamda)

@jit
def calc_W_Cao_empirical(ZH_lin, ZDR):
    return ZH_lin * 10.**(-0.0493*ZDR**3. + 0.430*ZDR**2. - 1.524*ZDR - 3.019)

@jit
def calc_D0_Cao_empirical(ZDR):
    return 0.0436*ZDR**3. - 0.216*ZDR**2. + 1.076*ZDR + 0.659

@jit
def calc_Nt_Cao_empirical(ZH_lin, ZDR):
    return ZH_lin * 10.**(-0.0837*ZDR**3. + 0.702*ZDR**2. - 2.062*ZDR + 0.794)

@jit
def calc_RR_Cao_empirical(ZH_lin, ZDR):
    return ZH_lin * 10.**(-0.0363*ZDR**3. + 0.316*ZDR**2. - 1.178*ZDR - 1.1964)

@jit
def calc_Nt_Brandes2004_empirical(ZH_lin, ZDR):
    return 2.085 * ZH_lin * 10.**(0.728 * ZDR**2. - 2.066 * ZDR)

@jit
def calc_W_Brandes2004_empirical(ZH_lin, ZDR):
    return 5.589 * 10**-4. * ZH_lin * 10.**(0.223 * ZDR**2. - 1.124 * ZDR)

@jit
def calc_RR_Brandes2004_empirical(ZH_lin, ZDR):
    return 7.6 * 10**-3 * ZH_lin * 10.**(0.165 * ZDR**2. - 0.897 * ZDR)

@jit
def calc_D0_Brandes2004_empirical(ZDR):
    return 0.717 + 1.479 * ZDR - 0.725 * ZDR**2. + 0.171 * ZDR**3.

@jit
def calc_sigma_Brandes2004_empirical(ZDR):
    return 0.163 + 0.519 * ZDR - 0.0247 * ZDR**2.


def retrieval_Cao_xr(ZH, ZDR, ND, D, dD, fa2, fb2, wavelength, mu_lamda_coeff, ZDR_thresh=3.,
                     retrieval_tag=''):

    RR_key = 'RR_retr_' + retrieval_tag
    D0_key = 'D0_retr_' + retrieval_tag
    mu_key = 'mu_retr_' + retrieval_tag
    lamda_key = 'lamda_retr_' + retrieval_tag
    N0_key = 'N0_retr_' + retrieval_tag
    Nt_key = 'Nt_retr_' + retrieval_tag
    W_key = 'W_retr_' + retrieval_tag
    sigma_key = 'sigma_retr_' + retrieval_tag
    Dm_key = 'Dm43_retr_' + retrieval_tag
    ND_key = 'ND_retr_' + retrieval_tag

    ZH_arr = ZH.values
    ZDR_arr = ZDR.values
    full_len = len(D.values)
    trunc_len = len(fa2)
    D_arr = D.values[:trunc_len]
    dD_arr = dD.values[:trunc_len]

    ntimes = len(ZH)
    retr_dict_list = []
    for t in range(ntimes):
        # retr_dict = retrieval_Cao(ZH_arr[t].item(), ZDR_arr[t].item(), D_arr, dD_arr,
        #                           fa2, fb2, wavelength, mu_lamda_coeff, ZDR_thresh=ZDR_thresh,
        #                           retrieval_tag=retrieval_tag)
        retr_tuple = retrieval_Cao(ZH_arr[t].item(), ZDR_arr[t].item(), D_arr, dD_arr,
                                   fa2, fb2, wavelength, mu_lamda_coeff, ZDR_thresh=ZDR_thresh,
                                   retrieval_tag=retrieval_tag)
        retr_dict = {
            RR_key: retr_tuple[0],
            D0_key: retr_tuple[1],
            mu_key: retr_tuple[2],
            lamda_key: retr_tuple[3],
            N0_key: retr_tuple[4],
            Nt_key: retr_tuple[5],
            W_key: retr_tuple[6],
            sigma_key: retr_tuple[7],
            Dm_key: retr_tuple[8],
            ND_key: retr_tuple[9]
        }
        # DTD: test weird numba error by temporarily removing ND from dictionary, so set it to
        # empty here
        ND_retr = np.array(retr_dict['ND_retr_{}'.format(retrieval_tag)])
        # ND_retr = np.empty_like(D_arr)
        #print(ND_retr)
        ND_retr = np.append(ND_retr, [np.nan] * (full_len - trunc_len))

        retr_dict['ND_retr_{}'.format(retrieval_tag)] = ND_retr
        retr_dict_list.append(retr_dict)

    retr_dict_alltimes = {k: [dic[k] for dic in retr_dict_list] for k in retr_dict_list[0]}
    ND_retr = retr_dict_alltimes['ND_retr_{}'.format(retrieval_tag)]
    ND_retr_da = ND.copy(data=ND_retr)
    retr_dict_alltimes['ND_retr_{}'.format(retrieval_tag)] = ND_retr_da
    for name, values in retr_dict_alltimes.items():
        if name not in 'ND_retr_{}'.format(retrieval_tag):
            retr_dict_alltimes[name] = ZH.copy(data=values)
    return retr_dict_alltimes


# @jit(parallel=True)
@jit
def retrieval_Cao(ZH, ZDR, D, dD, fa2, fb2, wavelength, mu_lamda_coeff, ZDR_thresh=3.,
                  retrieval_tag=''):
    Dmx = 9.     # maximum diameter of 9 mm
    wavelength_mm = wavelength * 10.    # radar wavelength in mm
    # TODO: think about allowing for correction of fallspeed based on air density
    v_terminal = PIPS.calc_empirical_fallspeed(D)

    fp = -1.0   # # used to test loop condition for ZDR
    delta_lamda = 0.1   # step through lamda values by this amount (cm)
    min_lamda = 1.
    lamda_1 = min_lamda
    max_lamda = 30.     # Maximum lamda value
    lamda_range = np.arange(min_lamda, max_lamda + delta_lamda, delta_lamda)
    ND = np.empty_like(D)

    ZDR_lin = 10.**(ZDR / 10.)  # ZDR in linear units
    ZH_lin = 10.**(ZH / 10.)      # ZH in linear units

#    if np.all([ZH_lin >= 1., Dmx >= 0.4, ZDR_lin >= 1.]):
    if ZH_lin >= 1. and Dmx >= 0.4 and ZDR_lin >= 1.:
        if ZDR < 0.1:   # Use emperical relations when ZDR is small
            W = calc_W_Cao_empirical(ZH_lin, ZDR)
            D0 = calc_D0_Cao_empirical(ZDR)
            Nt = calc_Nt_Cao_empirical(ZH_lin, ZDR)
            RR = calc_RR_Cao_empirical(ZH_lin, ZDR)
            Dm = D0     # Why?
            mu = np.nan
            lamda = np.nan
            N0 = np.nan
            sigma = np.nan
        elif ZDR > ZDR_thresh:
            W = np.nan
            D0 = np.nan
            Dm = np.nan
            mu = np.nan
            lamda = np.nan
            N0 = np.nan
            sigma = np.nan
            RR = np.nan
            Nt = np.nan
        else:
            # Do the retrieval iteration
            for lamda in lamda_range:
                tflag = 0
                # Compute mu from mu_lamda relation
                mu = mu_lamda_coeff[0] + mu_lamda_coeff[1]*lamda + \
                     mu_lamda_coeff[2] * lamda**2.
                #print(lamda, mu)
                ND[:] = D**mu * np.exp(-lamda*D)
                for i in range(ND.shape[0]):
                    if D[i] > Dmx or D[i] < 0.3125:
                        ND[i] = 0.0
                # ND[D > Dmx] = 0.0   # Zero out bins greater than Dmx
                # Test removal of diameters less than lowest recorded parsivel bin:
                # ND[D < 0.3125] = 0.0
                ZH_lin_tmp = np.nansum(fa2*ND*dD)
                ZV_lin_tmp = np.nansum(fb2*ND*dD)
                if ZV_lin_tmp == 0.:
                    continue
                ZDR_lin_tmp = ZH_lin_tmp / ZV_lin_tmp
                ZDR_lin_diff = ZDR_lin - ZDR_lin_tmp
                #print(ZDR_lin, ZDR_lin_tmp)
                # if the observed linear zdr is greater than the calculated one, exit loop with
                # these values
                if ZDR_lin_diff * fp < 0:
                    #print("Here!")
                    lamda = lamda_1 - np.abs(fp / (ZDR_lin_diff - fp)) * delta_lamda
                    #print(lamda)
                    tflag = 1
                    break
                else:
                    fp = ZDR_lin_diff
                    lamda_1 = lamda

            # think this means a good lambda value was never found or that calculated zdr larger
            # than observed?
            if tflag == 0.:
                N0 = np.nan
                mu = np.nan
                lamda = np.nan
                D0 = np.nan
                Dm = np.nan
                sigma = np.nan

                # TODO: Zhang's code uses the relations from Brandes et al. 2004.
                # But since these are based on *their* CG relation, it seems we'll have to
                # update these for ours...
                Nt = calc_Nt_Brandes2004_empirical(ZH_lin, ZDR)
                W = calc_W_Brandes2004_empirical(ZH_lin, ZDR)
                RR = calc_RR_Brandes2004_empirical(ZH_lin, ZDR)
                if ZDR <= 0.35:
                    D0 = calc_D0_Brandes2004_empirical(ZDR)
                    sigma = calc_sigma_Brandes2004_empirical(ZDR)
            else:
                # Recompute mu from CG relation
                mu = mu_lamda_coeff[0] + mu_lamda_coeff[1]*lamda + \
                        mu_lamda_coeff[2] * lamda**2.

                # TODO: below is transcription from Cao and Zhang's original code. Will clean up
                # to make things clearer...
                ND[:] = D**mu * np.exp(-lamda*D)   # N0 is multiplied here later... confusing
                for i in range(ND.shape[0]):
                    if D[i] > Dmx or D[i] < 0.3125:
                        ND[i] = 0.0
                # ND[D > Dmx] = 0.0
                # Test removal of diameters less than lowest recorded parsivel bin:
                # ND[D < 0.3125] = 0.0
                ZH_0 = np.nansum(fa2*ND*dD)
                M3 = np.nansum(ND*D**3.*dD)
                M4 = np.nansum(ND*D**4.*dD)
                RR_0 = np.nansum(np.pi*6.e-4*v_terminal*D**3.*ND*dD)
                Dm_0 = np.nansum(ND*dD)

                Dm = M4 / M3
                sigma2 = np.nansum((D - Dm)**2. * ND * D**3. * dD)
                sigma = np.sqrt(sigma2 / M3)
                #print(ZH_lin, radar.K2, wavelength_mm, ZH_0)
                fhh2 = ZH_lin * np.pi**4. * radar.K2 / 4. / wavelength_mm**4.
                N0 = fhh2 / ZH_0
                Nt = Dm_0 * N0
                RR = RR_0 * N0
                W = 1.e-3 * np.pi / 6. * N0 * M3

                # Compute ND again...this is convoluted... this time *actual* ND
                ND[:] = N0 * D**mu * np.exp(-lamda * D)
                for i in range(ND.shape[0]):
                    if D[i] > Dmx or D[i] < 0.3125:
                        ND[i] = 0.0
                # ND[D > Dmx] = 0.0
                # Test removal of diameters less than lowest recorded parsivel bin:
                # ND[D < 0.3125] = 0.0
                M3_bin = ND*D**3.*dD
                mtp1 = np.nansum(M3_bin)
                mtp2 = 0.
                m = 0
                mas = 0.

                # Compute D0
                # TODO: make this consistent with calc_D0_bin
                if mtp1 > 0.:
                    while mas < 0.5:
                        masp = mtp2 / mtp1
                        mtp2 = mtp2 + M3_bin[m]
                        m = m + 1
                        mas = mtp2 / mtp1
                    D0 = D[m] - (mas - 0.5) / (mas - masp) * dD[m]
                else:
                    D0 = 0.
    else:
        RR = np.nan
        D0 = np.nan
        Dm = np.nan
        mu = np.nan
        lamda = np.nan
        N0 = np.nan
        Nt = np.nan
        W = np.nan
        sigma = np.nan

    ND[:] = N0 * D**mu * np.exp(-lamda * D)
    for i in range(ND.shape[0]):
        if D[i] > Dmx or D[i] < 0.3125:
            ND[i] = 0.0
    # ND[D > Dmx] = 0.0
    # ND[D < 0.3125] = 0.0

    # Have to do this because numba doesn't like str.format()...
    # RR_key = 'RR_retr_' + retrieval_tag
    # D0_key = 'D0_retr_' + retrieval_tag
    # mu_key = 'mu_retr_' + retrieval_tag
    # lamda_key = 'lamda_retr_' + retrieval_tag
    # N0_key = 'N0_retr_' + retrieval_tag
    # Nt_key = 'Nt_retr_' + retrieval_tag
    # W_key = 'W_retr_' + retrieval_tag
    # sigma_key = 'sigma_retr_' + retrieval_tag
    # Dm_key = 'Dm_retr_' + retrieval_tag
    # ND_key = 'ND_retr_' + retrieval_tag

    # retr_dict = {
    #     RR_key: RR,
    #     D0_key: D0,
    #     mu_key: mu,
    #     lamda_key: lamda,
    #     N0_key: N0,
    #     Nt_key: Nt,
    #     W_key: W,
    #     sigma_key: sigma,
    #     Dm_key: Dm,
    #     # ND_key: ND[:]
    # }
    # retr_dict = {
    #     'RR_retr_{}'.format(retrieval_tag): RR,
    #     'D0_retr_{}'.format(retrieval_tag): D0,
    #     'mu_retr_{}'.format(retrieval_tag): mu,
    #     'lamda_retr_{}'.format(retrieval_tag): lamda,
    #     'N0_retr_{}'.format(retrieval_tag): N0,
    #     'Nt_retr_{}'.format(retrieval_tag): Nt,
    #     'W_retr_{}'.format(retrieval_tag): W,
    #     'sigma_retr_{}'.format(retrieval_tag): sigma,
    #     'Dm_retr_{}'.format(retrieval_tag): Dm,
    #     'ND_retr_{}'.format(retrieval_tag): ND
    # }

    # return retr_dict

    return RR, D0, mu, lamda, N0, Nt, W, sigma, Dm, ND








