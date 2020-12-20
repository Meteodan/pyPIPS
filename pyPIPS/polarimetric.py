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
    # Truncate the diameter dimension of Nd to account for this.
    # Also transpose Nd to allow for numpy broadcasting
    Nd = Nd[:, :np.size(fa2)] # Nd[:np.size(fa2), :].T
    intv = intv[:np.size(fa2)]

    lamda = wavelength * 10.  # Get radar wavelength in mm
    Kw2 = 0.93  # Dielectric factor for water
    sar_h = fa2 * Nd * intv
    sar_v = fb2 * Nd * intv
    sar_hv = fab * Nd * intv
    fsar = far * Nd * intv

    # TODO: return binned values (by diameter) in addition to total
    Zh_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_h
    Zv_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_v
    Zhv_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_hv
    Kdp_bin = 180. * lamda / np.pi * fsar * 1.e-3
    dBZ_bin = 10. * np.log10(Zh_bin)
    ZDR_lin_bin = Zh_bin / Zv_bin
    ZDR_bin = 10. * np.log10(np.maximum(1.0, ZDR_lin_bin))
    temp = Zh_bin * Zv_bin
    # Added by Jess (was temp > 0).  Find out why...
    rhv_bin = np.where(Zh_bin != Zv_bin, Zhv_bin / (np.sqrt(temp)), np.nan)
    Zh = np.sum(Zh_bin, axis=1)
    Zv = np.sum(Zv_bin, axis=1)
    Zhv = np.abs(np.sum(Zhv_bin, axis=1))
    Kdp = np.sum(Kdp_bin, axis=1)
    dBZ = 10. * np.log10(np.sum(Zh_bin, axis=1))
    temp = Zh / Zv
    ZDR = 10. * np.log10(np.maximum(1.0, temp))
    temp = Zh * Zv
    # Added by Jess (was temp > 0).  Find out why...
    rhv = np.where(Zh != Zv, Zhv / (np.sqrt(temp)), np.nan)
    # np.savetxt('temp.txt', temp)

    dualpol_dict = {'ZH_bin': Zh_bin, 'ZV_bin': Zv_bin, 'ZHV_bin': Zhv_bin, 'REF_bin': dBZ_bin,
                    'ZDR_lin_bin': ZDR_lin_bin, 'ZDR_bin': ZDR_bin, 'KDP_bin': Kdp_bin,
                    'RHO_bin': rhv_bin, 'ZH': Zh, 'ZV': Zv, 'ZHV': Zhv, 'REF': dBZ, 'ZDR': ZDR,
                    'KDP': Kdp, 'RHO': rhv, 'intv': intv, 'd': d, 'fa2': fa2, 'fb2': fb2}

    return dualpol_dict

def calpolrain_xr(wavelength, filename, Nd, intv):
    """Given backscattering amplitudes and a discrete distribution N(D) (m^-4), compute
       polarimetric variables for each bin."""

    d, far_b, fbr_b, far_f, fbr_f = readtmatrix(filename)
    fa2, fb2, fab, fba, far = calbackscatterrain(far_b, fbr_b, far_f, fbr_f)

    # There may be more bins in the given Nd than are read in from the file.
    # This is because the Nd contains bins above the maximum size of rain (~9 mm).
    # Truncate the diameter dimension of Nd to account for this.
    # Also transpose Nd to allow for numpy broadcasting
    # Nd = Nd[:, :np.size(fa2)] # Nd[:np.size(fa2), :].T
    Nd = Nd.isel(diameter_bin=slice(0, np.size(fa2)))
    intv = intv[:np.size(fa2)]

    lamda = wavelength * 10.  # Get radar wavelength in mm
    Kw2 = 0.93  # Dielectric factor for water
    sar_h = fa2 * Nd * intv
    sar_v = fb2 * Nd * intv
    sar_hv = fab * Nd * intv
    fsar = far * Nd * intv

    # TODO: return binned values (by diameter) in addition to total
    Zh_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_h
    Zv_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_v
    Zhv_bin = 4. * lamda**4. / (np.pi**4. * Kw2) * sar_hv
    Kdp_bin = 180. * lamda / np.pi * fsar * 1.e-3
    dBZ_bin = 10. * np.log10(Zh_bin)
    ZDR_lin_bin = Zh_bin / Zv_bin
    ZDR_bin = 10. * np.log10(np.maximum(1.0, ZDR_lin_bin))
    temp = Zh_bin * Zv_bin
    # Added by Jess (was temp > 0).  Find out why...
    rhv_bin = np.where(Zh_bin != Zv_bin, Zhv_bin / (np.sqrt(temp)), np.nan)
    Zh = Zh_bin.sum(dim='diameter_bin') # np.sum(Zh_bin, axis=1)
    Zv = Zv_bin.sum(dim='diameter_bin') # np.sum(Zv_bin, axis=1)
    Zhv = np.abs(Zhv_bin.sum(dim='diameter_bin'))
    Kdp = Kdp_bin.sum(dim='diameter_bin')
    dBZ = 10. * np.log10(Zh)
    temp = Zh / Zv
    ZDR = 10. * np.log10(np.maximum(1.0, temp))
    temp = Zh * Zv
    # Added by Jess (was temp > 0).  Find out why...
    rhv = np.where(Zh != Zv, Zhv / (np.sqrt(temp)), np.nan)
    # np.savetxt('temp.txt', temp)

    dualpol_dict = {'ZH_bin': Zh_bin, 'ZV_bin': Zv_bin, 'ZHV_bin': Zhv_bin, 'REF_bin': dBZ_bin,
                    'ZDR_lin_bin': ZDR_lin_bin, 'ZDR_bin': ZDR_bin, 'KDP_bin': Kdp_bin,
                    'RHO_bin': rhv_bin, 'ZH': Zh, 'ZV': Zv, 'ZHV': Zhv, 'REF': dBZ, 'ZDR': ZDR,
                    'KDP': Kdp, 'RHO': rhv, 'intv': intv, 'd': d, 'fa2': fa2, 'fb2': fb2}

    return dualpol_dict
