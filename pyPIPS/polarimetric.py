"""
polarimetric.py: a set of functions dealing with computing polarimetric radar variables for PSDs
"""
import numpy as np

# TODO: consider interfacing with pytmatrix and pyDSD


def readtmatrix(filename):
    """Reads a scattering amplitude lookup table created using the ARPS tmatrix program"""

    data = np.loadtxt(filename, skiprows=1)
    d = data[:, 0]
    far_b = data[:, 1] + 1j * data[:, 2]
    fbr_b = data[:, 3] + 1j * data[:, 4]
    far_f = data[:, 5] + 1j * data[:, 6]
    fbr_f = data[:, 7] + 1j * data[:, 8]

    return d, far_b, fbr_b, far_f, fbr_f


def calbackscatterrain(far_b, fbr_b, far_f, fbr_f):
    """Calculates rain backscattering amplitudes.
       Based on similar code in dualpara.f90"""

    fa2 = (np.abs(far_b))**2.
    fb2 = (np.abs(fbr_b))**2.
    fab = far_b * np.conjugate(fbr_b)
    fba = fbr_b * np.conjugate(far_b)
    far = np.real(far_f - fbr_f)

    return fa2, fb2, fab, fba, far


def calpolrain(wavelength, filename, Nd, intv):
    """Given backscattering amplitudes and a discrete distribution N(D) (m^-4), compute
       polarimetric variables for each bin."""

    d, far_b, fbr_b, far_f, fbr_f = readtmatrix(filename)
    fa2, fb2, fab, fba, far = calbackscatterrain(far_b, fbr_b, far_f, fbr_f)

    # There may be more bins in the given Nd than are read in from the file.
    # This is because the Nd contains bins above the maximum size of rain (~9 mm).
    # Truncate the first dimension of Nd to account for this.
    # Also transpose Nd to allow for numpy broadcasting
    Nd = Nd[:np.size(fa2), :].T
    intv = intv[:np.size(fa2)]

    lamda = wavelength * 10.  # Get radar wavelength in mm
    Kw2 = 0.93  # Dielectric factor for water
    sar_h = fa2 * Nd * intv
    sar_v = fb2 * Nd * intv
    sar_hv = fab * Nd * intv
    fsar = far * Nd * intv

    Zh = 4. * lamda**4. / (np.pi**4. * Kw2) * np.sum(sar_h, axis=1)
    Zv = 4. * lamda**4. / (np.pi**4. * Kw2) * np.sum(sar_v, axis=1)
    Zhv = 4. * lamda**4. / (np.pi**4. * Kw2) * np.abs(np.sum(sar_hv, axis=1))
    Kdp = 180. * lamda / np.pi * np.sum(fsar, axis=1) * 1.e-3
    dBZ = 10. * np.log10(Zh)
    temp = Zh / Zv
    ZDR = 10. * np.log10(np.maximum(1.0, temp))
    temp = Zh * Zv
    # Added by Jess (was temp > 0).  Find out why...
    rhv = np.where(Zh != Zv, Zhv / (np.sqrt(temp)), 0.0)
    # np.savetxt('temp.txt', temp)

    dualpol_dict = {'ZH': Zh, 'ZV': Zv, 'ZHV': Zhv, 'dBZ': dBZ, 'ZDR': ZDR, 'KDP': Kdp, 'RHV': rhv,
                    'intv': intv, 'd': d, 'fa2': fa2, 'fb2': fb2}

    return dualpol_dict
