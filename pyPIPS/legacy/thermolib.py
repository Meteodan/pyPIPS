# thermolib.py
# A collection of thermodynamic functions

import numpy as np

Rd = 287.0  # Gas constant for dry air (J/kg/K)
Rv = 461.51  # Gas constant for water vapor (J/kg/K)
cp = 1005.6  # Specific heat at constant pressure of dry air (J/kg/K)
rddcp = Rd / cp
rddrv = Rd / Rv
rvdrd = Rv / Rd
Lv = 2.501e6  # Latent heat of vaporization (J/kg)
Lf = 3.34e5  # Latent heat of freezing (J/kg)
Ls = Lv + Lf  # Latent heat of sublimation (J/kg)

p0 = 100000.0  # Reference pressure in Pa
rho_ref = 1.0  # Reference air density (kg/m^3)
rhow = 1000.0  # Density of liquid water (kg/m^3)

satfwa = 1.0007
satfwb = 3.46e-8
satewa = 611.21
satewb = 17.502
satewc = 32.18
satfia = 1.0003
satfib = 4.18e-8
sateia = 611.15
sateib = 22.452
sateic = 0.6


def caltheta(p, T):
    """ Calculate potential temperature from temperature and pressure"""
    return T * (p0 / p)**rddcp


def calthetav(pt, qv):
    """ Calculate virtual potential temperature from potential temperature and water vapor
        specific humidity"""
    # Compute water vapor mixing ratio from specific humidity

    wv = qv / (1.0 - qv)

    return pt * (1.0 + rvdrd * wv) / (1.0 + wv)


def calptr(pt, qv, qli):
    """ Calculate density potential temperature from potential temperature,
            water vapor specific humidity, and liquid+ice mixing ratio"""
    # Compute water vapor mixing ratio

    wv = qv / (1.0 - qv)

    return pt * (1.0 + rvdrd * wv) / (1.0 + wv + qli)


def calT(p, pt):
    """ Calculate temperature from potential temperature and pressure"""
    return pt * (p / p0)**rddcp


def calTv(p, pt, qv):
    """ Calculate virtual temperature from pressure, potential temperature,
        and water vapor specific humidity"""
    T = calT(p, pt)
    # Compute water vapor mixing ratio

    wv = qv / (1.0 - qv)

    return T * (1.0 + rvdrd * wv) / (1.0 + wv)


def calTr(p, pt, qv, qli):
    """ Calculate density temperature from pressure, potential temperature,
        water vapor specific humidity, and liquid+ice mixing ratio"""
    T = calT(p, pt)

    # Compute water vapor mixing ratio
    # Question: is qv in ARPS mixing ratio or specific humidity?

    wv = qv / (1.0 - qv)

    return T * (1.0 + rvdrd * wv) / (1.0 + wv + qli)


def calrho(p, pt, qv):
    """ Calculate moist density from pressure, potential temperature, and
        water vapor specific humidity"""
    Tv = calTv(p, pt, qv)
    return p / (Rd * Tv)


def calrhod(p, pt):
    """ Calculate dry air density from pressure and potential temperature"""
    T = calT(p, pt)
    return p / (Rd * T)


def calTd(p, qv):
    """ Calculate dewpoint temperature from pressure and water vapor specific humidity.
    This uses the forumulae from Buck (1981, JAM) as used in ARPS thermolib3d.f90. """
    f = satfwa + satfwb * p

    # Calculate vapor pressure from pressure and water vapor specific humidity
    e = p * qv / (rddrv + (1.0 - rddrv) * qv)

    A = np.log(e / (satewa * f))

    Td = (satewc * A - 273.15 * satewb) / (A - satewb)
    return Td


def calTdfromRH(p, T, RH):
    """ Calculate dewpoint temperature from pressure, temperature, and relative humidity.
    This uses the forumulae from Buck (1981, JAM) as used in ARPS thermolib3d.f90. """
    f = satfwa + satfwb * p
    es = cales(p, T)
    e = cale(RH, es)
    A = np.log(e / (satewa * f))
    Td = (satewc * A - 273.15 * satewb) / (A - satewb)
    return Td


def calpte(p, pt, qv):
    """ Calculate equivalent potential temperature based on Bolton's (1980) approximation
    from T, p, and qv.  Note, uses the calTd and calT functions defined previously."""
    # First, calculate temperature and dewpoint temperature

    T = calT(p, pt)
    Td = calTd(p, qv)

    # Calculate temperature at LCL

    Tl = 1.0 / ((1.0 / (Td - 56.)) + np.log(T / Td) / 800.0) + 56.0

    # Convert qv to water vapor mixing ratio
    # TODO: this isn't used right now. Should it be?
    # qvm = qv / (1.0 - qv)

    # Finally compute theta_e

    pte = (T * (p0 / p)**(0.2854 * (1.0 - (0.28e-3) * 1000.0 * qv))) * \
        np.exp((3.376 / Tl - 0.00254) * 1000.0 * qv * (1.0 + (0.81e-3) * 1000.0 * qv))

    return pte


def calptefromRHpT(RH, p, T):
    """
    Calculate equivalent potential temperature given pressure, temperature, and relative humidity.
    """
    # Compute qv

    qv = calqv(RH, p, T)

    # Compute dewpoint temperature

    Td = calTd(p, qv)

    # Calculate temperature at LCL

    Tl = 1.0 / ((1.0 / (Td - 56.)) + np.log(T / Td) / 800.0) + 56.0

    # Convert qv to water vapor mixing ratio
    # TODO: this isn't used right now. Should it be?
    # qvm = qv / (1.0 - qv)

    # Finally compute theta_e

    pte = (T * (p0 / p)**(0.2854 * (1.0 - (0.28e-3) * 1000.0 * qv))) * \
        np.exp((3.376 / Tl - 0.00254) * 1000.0 * qv * (1.0 + (0.81e-3) * 1000.0 * qv))

    return pte


def calqvs(p, T):
    """ Calculate saturation water vapor mixing ratio from temperature and pressure"""
    # First calculate saturation vapor pressure

    f = satfwa + satfwb * p
    es = f * satewa * np.exp(satewb * (T - 273.15) / (T - satewc))

    # Now calculate saturation water vapor mixing ratio

    ws = rddrv * (es / (p - es))

    return ws


def cales(p, T):
    """ Calculate saturation vapor pressure from pressure and temperature """
    # First calculate saturation vapor pressure

    f = satfwa + satfwb * p
    es = f * satewa * np.exp(satewb * (T - 273.15) / (T - satewc))

    return es


def cale(RH, es):
    """ Calculate vapor pressure from saturation vapor pressure and relative humidity """
    return RH * es


def calqv(RH, p, T):
    """ Calculate water vapor mixing ratio from relative humidity, pressure, and temperature"""
    es = cales(p, T)
    e = cale(RH, es)
    qv = (Rd / Rv) * (e / (p - e))

    return qv


def calRH(p, T, Td):
    """ Calculate RH from pressure, temperature, and dewpoint temperature"""
    es = cales(p, T)
    e = cales(p, Td)
    return e / es
