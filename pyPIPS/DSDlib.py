"""
DSDlib.py: This is a library of functions to calculate various DSD-related parameters
Some of these were originally written in fortran but re-written using python and numpy
"""

import numpy as np
from scipy.special import gamma as gamma_, gammainc as gammap, gammaln as gammln
from . import thermolib as thermo
from .utils import first_nonzero, last_nonzero, enable_xarray_wrapper

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
    """Given two moments, p and q, and lamda,Nt, and alpha, compute the
       mean diameter associated with the ratio of p to q"""

    M_p = calc_Mp(p, lamda, Nt, alpha)
    M_q = calc_Mp(q, lamda, Nt, alpha)

    return (M_p / M_q)**(1. / (p - q))


def calc_qr_gamma(rhoa, N0, lamda, alpha):
    """[summary]

    Parameters
    ----------
    rhoa : [type]
        [description]
    D : [type]
        [description]
    N0 : [type]
        [description]
    lamda : [type]
        [description]
    alpha : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    GR2 = gamma_(4. + alpha)
    return (cmr / rhoa) * N0 * GR2 / lamda**(alpha + 4.)


def calc_Zr_gamma(rhoa, q, Nt, alpha):
    """[summary]

    Parameters
    ----------
    rhoa : [type]
        [description]
    q : [type]
        [description]
    Nt : [type]
        [description]
    alpha : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    Gr = ((6. + alpha) * (5. + alpha) * (4. + alpha)) / \
        ((3. + alpha) * (2. + alpha) * (1. + alpha))
    Zr = ((1. / cmr)**2.) * Gr * ((rhoa * q)**2.) / Nt
    return 10.0 * np.log10(1.e18 * Zr)


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
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]
    moment : int, optional
        [description], by default 0

    Returns
    -------
    [type]
        [description]
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
#     print(pro_Dm.loc[dict(time_10s='2016-03-31T22:30:00')])
#     D0 = pro_Dm.quantile(0.5, dim='diameter_bin')
#     print(D0.loc[dict(time_10s='2016-03-31T22:30:00')])
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
    medindices_m1[medindices_m1 < 0] = 0
    b1 = Dl[medindices]  # Lower boundaries of mass-midpoint bin
    b2 = Dr[medindices]  # Upper boundaries of mass-midpoint bin
    pro_med = pro.loc[dict(diameter_bin=medindices)]
    pro_cumsum_med_m1 = pro_cumsum.loc[dict(diameter_bin=medindices_m1)]
    # Now we can calculate D0 by linearly interpolating diameters within
    # a given bin bounded by Dl and Dr which contains the half-mass point
    D0 = b1 + ((0.5 - pro_cumsum_med_m1) / pro_med) * (b2 - b1)
    # Don't let D0 be any smaller than the midpoint of the smallest bin
    D0[D0 < Dm[0]] = Dm[0]

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

    alphaMax = 40.0

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
       code.  Note, all thermodynamic variables are assumed to be in SI units."""

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
       assuming a gamma-diameter distribution and given q,N, and alpha"""

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
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return calc_moment_bin(ND, moment=0)


def calc_lwc_qr_from_bins(ND, rho):
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]
    rho : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    M3 = calc_moment_bin(ND, moment=3)
    LWC = cmr * M3
    qr = LWC / rho

    return LWC, qr


def calc_dBZ_from_bins(ND):
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    M6 = calc_moment_bin(ND, moment=6)
    return 10. * np.log10(1.e18 * M6)


def calc_synthetic_bins(D_min, D_max, num_bins):
    """[summary]

    Parameters
    ----------
    D_min : [type]
        [description]
    D_max : [type]
        [description]
    num_bins : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return np.linspace(D_min, D_max, num=num_bins)


def fit_DSD_MM24(M2, M4):
    """[summary]

    Parameters
    ----------
    M2 : [type]
        [description]
    M4 : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    lamda = np.where(M4 == 0.0, 0.0, ((M2 * gamma5) / (M4 * gamma3))**(1./2.))
    N0 = (M2 * lamda**3.) / gamma3
    mu = 0.
    return N0, lamda, mu


def fit_DSD_MM36(M3, M6):
    """[summary]

    Parameters
    ----------
    M3 : [type]
        [description]
    M6 : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    lamda = np.where(M6 == 0.0, 0.0, ((M3 * gamma7) / (M6 * gamma4))
                     ** (1. / 3.))
    N0 = (M3 * lamda**4.) / gamma4
    mu = 0.
    return N0, lamda, mu


def fit_DSD_MM346(M3, M4, M6):
    """[summary]

    Parameters
    ----------
    M3 : [type]
        [description]
    M4 : [type]
        [description]
    M6 : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    G = (M4**3.) / ((M3**2.) * M6)
    G = np.ma.masked_invalid(G)
    mu = (11. * G - 8. + (G * (G + 8.))**(1. / 2.))/(2. * (1. - G))
    mu = np.ma.masked_invalid(mu)
    lamda = (M3 * (mu + 4.)) / M4
    lamda = np.ma.masked_invalid(lamda)
    N0 = (M3 * lamda**(mu + 4.))/(gamma_(mu + 4.))

    return N0, lamda, mu


def fit_DSD_MM246_old(M2, M4, M6):
    """[summary]

    Parameters
    ----------
    M2 : [type]
        [description]
    M4 : [type]
        [description]
    M6 : [type]
        [description]

    Returns
    -------
    [type]
        [description]
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
    """[summary]

    Parameters
    ----------
    M2 : [type]
        [description]
    M4 : [type]
        [description]
    M6 : [type]
        [description]
    lamda_limit : [type], optional
        [description], by default 20000.
    mu_limit : [type], optional
        [description], by default 30.

    Returns
    -------
    [type]
        [description]
    """
    G = (M4**2.) / (M2 * M6)
    G = np.ma.masked_invalid(G)
    mu = ((7. - 11. * G) - ((7. - 11. * G)**2. - 4. * (G - 1.)
                            * (30. * G - 12.))**(1. / 2.)) / (2. * (G - 1.))
    mu = np.ma.masked_where(mu > mu_limit, mu)
    mu = np.ma.masked_invalid(mu)
    lamda = ((M2 * (mu + 3.) * (mu + 4.)) / (M4))**(1. / 2.)
    mu = np.ma.masked_where((lamda > lamda_limit) | (mu > 30.), mu)
    mu = np.ma.masked_invalid(mu)
    lamda = np.ma.masked_where(lamda > lamda_limit, lamda)
    lamda = np.ma.masked_invalid(lamda)
    N0 = (M4 * lamda**(mu + 5.)) / (gamma_(mu + 5.))

    return N0, lamda, mu


def fit_DSD_MM234(M2, M3, M4):
    """[summary]

    Parameters
    ----------
    M2 : [type]
        [description]
    M3 : [type]
        [description]
    M4 : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    mu = (3. * M2 * M4 - 4. * M3**2.) / (M3**2. - M2 * M4)
    mu = np.ma.masked_invalid(mu)
    lamda = (M3 * (mu + 4.)) / M4
    lamda = np.ma.masked_invalid(lamda)
    N0 = (M3 * lamda**(mu + 4.))/(gamma_(mu + 4.))

    return N0, lamda, mu


def get_max_min_diameters(ND):
    """[summary]

    Parameters
    ----------
    ND : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    D_min = first_nonzero(ND, 0)
    D_max = last_nonzero(ND, 0)

    return D_min, D_max


def fit_DSD_TMM(M2, M4, M6, D_min, D_max):
    """[summary]

    Parameters
    ----------
    M2 : [type]
        [description]
    M4 : [type]
        [description]
    M6 : [type]
        [description]
    D_min : [type]
        [description]
    D_max : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    G = (M4**2.) / (M2 * M6)  # moment ratio based on untruncated moments
    G = np.ma.masked_invalid(G)

    # Get initial estimate of lamda and mu from the untruncated fit:
    _, lamda_init, mu_init = fit_DSD_MM246(M2, M4, M6)

    # TODO: make sure everything works with xarray and see if there's an efficient way to
    # vectorize this algorithm across time. Maybe wrap Guifu's original fortran code is the
    # best approach.
    numtimes = np.size(M2)

    mu_tmf = []
    lamda_tmf = []
    # N0_tmf = []

    for t in range(numtimes):
        LDmx = lamda_init[t] * D_max[t]
        for x in range(10):
            mu = mu_init[t]
            # truncated moment ratio below. Equation A8 from Thurai
            gm3 = gammap(3. + mu, LDmx) * np.exp(gammln(3. + mu))
            gm5 = gammap(5. + mu, LDmx) * np.exp(gammln(5. + mu))
            gm7 = gammap(7. + mu, LDmx) * np.exp(gammln(7. + mu))
            z0 = G[t] - gm5**2. / gm3 / gm7
            z1 = G[t] - gm5**2. / gm3 / gm7

            while z1 / z0 > 0.0:
                mu = mu - 0.01
                gm3 = gammap(3. + mu, LDmx) * np.exp(gammln(3. + mu))
                gm5 = gammap(5. + mu, LDmx) * np.exp(gammln(5. + mu))
                gm7 = gammap(7. + mu, LDmx) * np.exp(gammln(7. + mu))
                z1 = G[t] - gm5**2. / gm3 / gm7

            lam_tmf = (M2[t] * gm5 / M4[t] / gm3)**0.5
            LDmx = lam_tmf * D_max[t]
        mu_tmf.append(mu)
        lamda_tmf.append(lam_tmf)

    mu_tmf = np.array(mu_tmf)
    lamda_tmf = np.array(lamda_tmf)
    LDmx = lamda_tmf * D_max
    N0_tmf = (M4 * lamda_tmf**(5. + mu_tmf)) / \
        (gammap(5. + mu_tmf, LDmx) * np.exp(gammln(5. + mu_tmf)))

    return N0_tmf, lamda_tmf, mu_tmf


def calc_binned_DSD_from_params(N0, lamda, alpha, D):
    """[summary]

    Parameters
    ----------
    N0 : [type]
        [description]
    lamda : [type]
        [description]
    alpha : [type]
        [description]
    D : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    # Get D to m
    D = D / 1000.

    return 1.e-3 * N0 * D**alpha * np.exp(-lamda * D)  # Units are m^-3 mm^-1


def calc_rain_axis_ratio(D, fit_name='Brandes_2002'):
    """[summary]

    Parameters
    ----------
    D : [type]
        [description]
    fit_name : str, optional
        [description], by default 'Brandes_2002'

    Returns
    -------
    [type]
        [description]
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
