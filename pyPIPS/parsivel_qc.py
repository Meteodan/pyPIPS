"""pyPIPS.parsivel_qc: contains functions for QC of Parsivel data. Mainly based on IDL routines
   from Katja Friedrich (CU)
"""
from __future__ import annotations

import numpy as np
from numpy import ma

from . import PIPS, parsivel_params
from .utils import enable_xarray_wrapper

# This dictionary contains mappings of QC "groups" as keys with the values being
# dictionaries of the appropriate QC flags for that group. Replaces the old way of using
# QC dictionaries in the individual case config files since we often want to compute several
# different QC'ed versions of the parsivel data at once and save the associated variables to
# the netCDF files.
PIPS_qc_dict = {
    'qc': {
        'strongwindQC': True,
        'splashingQC': True,
        'marginQC': True,
        'rainfallQC': False,
        'rainonlyQC': False,
        'hailonlyQC': False,
        'graupelonlyQC': False,
    },
    'roqc': {
        'strongwindQC': True,
        'splashingQC': True,
        'marginQC': True,
        'rainfallQC': False,
        'rainonlyQC': True,
        'hailonlyQC': False,
        'graupelonlyQC': False,
    }
}


# Some rain QC parameters

rain_QC_params = {
    'rain_fall_tol': 0.6,     # Rain fallspeed tolerance (fractional)
    'mask_fall_high': True,   # Whether to mask rain fall speeds larger than tolerance threshold
    'mask_fall_low': True,    # Whether to mask rain fall speeds smaller than tolerance threshold
    'mask_diam_high': False,  # Whether to mask rain diameters above 'high_diam_thresh'
    'mask_diam_low': False,   # Whether to mask rain diameters below 'low_diam_thresh'
    'high_diam_thresh': 9.0,  # High diameter threshold (mm)
    'low_diam_thresh': 1.0,   # Low diameter threshold (mm)
}

# TODO: wrap these in functions
# Create mask for splashing drops
bottom = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = [0, 1, 4, 7, 9, 11, 12, 13, 14, 14, 15, 16, 16, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
splashingmask = [[bottom[j] <= i <= top[j] for i in range(32)] for j in range(32)]
splashingmask = np.array(splashingmask).T

# Create mask for margin falls
bottom = [0, 8, 14, 17, 20, 21, 22, 23, 24, 25, 26, 26, 27, 27,
          28, 28, 29, 29, 29, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = np.arange((32), dtype='int')
top[:] = 31
top[23:32] = 0
marginmask = [[bottom[j] <= i <= top[j] for i in range(32)] for j in range(32)]
marginmask = np.array(marginmask).T

# Create mask for non-raindrops
bottom = [0, 1, 4, 7, 9, 11, 12, 13, 14, 14, 15, 16, 16, 19, 19,
          20, 20, 21, 21, 21, 23, 24, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = [0, 8, 14, 17, 20, 21, 22, 23, 24, 25, 26, 26, 27, 27, 28,
       28, 29, 29, 29, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0]
rainonlymask = [[not bottom[j] < i < top[j] for i in range(32)] for j in range(32)]
rainonlymask = np.array(rainonlymask).T

# Create mask for non-hail
bottom = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 15, 15, 16, 16, 17, 17, 18, 19, 19, 20, 20, 20]
top = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 24, 25, 25, 31, 32, 32, 32, 32, 32, 32, 32, 32]
hailonlymask = [[not bottom[j] < i < top[j] for i in range(32)] for j in range(32)]
hailonlymask = np.array(hailonlymask).T

# Create mask for strong wind conditions
strongwindmask = np.zeros((32, 32), dtype=bool)
strongwindmask[0:11, 20:32] = True

avg_diameter = parsivel_params.parsivel_parameters['avg_diameter_bins_mm']


def interpnan1D(a):
    """Replaces NaN's in a 1D array by interpolating from good values on either side"""
    ind = np.where(~np.isnan(a))[0]  # indices of valid values
    # Use valid values to interpolate to invalid values
    return np.interp(list(range(len(a))), ind, a[ind])


def get_fallspeed_mask(diam_bins, fall_bins, masklow=True, maskhigh=True, falltol=0.6):
    """Computes the fallspeed mask based on +/- the fractional tolerance given by falltol

    Parameters
    ----------
    diam_bins : list
        list of the middle points of the diameter bins
    fall_bins : list
        list of the middle points of the fallspeed bins
    masklow : bool, optional
        mask below the lower fallspeed threshold, by default True
    maskhigh : bool, optional
        mask above the upper fallspeed threshold, by default True
    falltol : float, optional
        fractional tolerance (+/-) for fallspeed, by default 0.6

    Returns : np.array
        A 2D boolean array indicating fallspeed vs. diameter bins that are masked
    """

    # Create mask for all particles with fall speeds outside of fractional tolerance
    rainvd = PIPS.calc_empirical_fallspeed(diam_bins)
    X, Y = np.meshgrid(rainvd, fall_bins)

    if maskhigh and masklow:
        fallspeedmask = np.where(np.abs((X - Y) / X) < falltol, False, True)
    elif masklow:     # Mask out speeds lower than tolerance
        fallspeedmask = np.where((Y - X) / X < -falltol, True, False)  # noqa: SIM300
    elif maskhigh:               # Mask out speeds higher than tolerance
        fallspeedmask = np.where((Y - X) / X > falltol, True, False)  # noqa: SIM300
    else:
        fallspeedmask = None

    return fallspeedmask


@enable_xarray_wrapper
def truncatedspectrumQC(countsMatrix):
    """Masks out bad records where tokens have been set to -999 because the record
       was truncated"""
    countsMatrix = ma.masked_array(countsMatrix, mask=np.where(countsMatrix == -999, True, False))

    return countsMatrix  # noqa: RET504


# FIXME: I can't get this function to work with allowing either xarray or regular numpy arrays,
# because it outputs a second array of flagged times that doesn't have the same dimensions. I could
# probably fix the decorator to handle it but honestly for now it's better if I just convert
# this to use only xarray and keep the other one around for numpy arrays.
# @enable_xarray_wrapper(output_core_dims=((), ()))
def strongwindQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine.  Removes drops affected by strong winds.
       Specifically, removes the entire 1-min DSD whenever there are large drops (> 3 mm) that have
       a low fall velocity (< 1 m/s)."""

    numtimes = np.size(countsMatrix, axis=0)
    # flaggedtimes = []

    countsMatrix['flagged_times'] = ('time', range(numtimes))
    # Flag times that contain wind contamination
    for t in range(numtimes):
        baddrops = \
            countsMatrix.isel(time=t, fallspeed_bin=slice(0, 11),
                              diameter_bin=slice(20, 32)).sum()
        # baddrops = np.sum(countsMatrix[t, 0:11, 20:32])
        # bigdrops = np.sum(countsMatrix[t, :, 23:32])
        # totaldrops = np.sum(countsMatrix[t, :])
#         if(baddrops > 0.02 * totaldrops):  # Try relaxing criterion to allow up to 2% of drops
#                                          # to be in mask area
#             print "Severe Wind contamination, masking entire PSD!"
#             countsMatrix[t, :] = -999.
#             flaggedtimes.append(2)
#         elif(baddrops > 0):  # Let the PSD through QC, but mask the offending drops
#             print "Wind contamination!"
#             countsMatrix[t, 0:11, 20:32] = -999.
#             flaggedtimes.append(1)
#         else:
#             flaggedtimes.append(0)

        if baddrops > 0:
            print("Severe Wind contamination, masking entire PSD!")  # noqa: T201
            countsMatrix[{'time': t}] = np.nan  # [t, :] = np.nan
            countsMatrix['flagged_times'][{'time': t}] = 2
            # flaggedtimes.append(2)
        else:
            countsMatrix['flagged_times'][{'time': t}] = 0
            # flaggedtimes.append(0)

    # countsMatrix = ma.masked_array(countsMatrix, mask=np.where(countsMatrix == -999., True,
    #                                                            False))

    return countsMatrix  # , np.array(flaggedtimes)


@enable_xarray_wrapper
def splashingQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from splashing"""

    numtimes = np.size(countsMatrix, axis=0)

    # Remove drops that likely result from splashing (use mask to index the count array)

    for t in range(numtimes):
        countsMatrix[t, splashingmask] = 0.0
        # countsMatrix[t, :] = ma.masked_array(countsMatrix[t, :], mask=splashingmask)

    return countsMatrix


@enable_xarray_wrapper
def marginQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from margin falls"""

    numtimes = np.size(countsMatrix, axis=0)

    # Remove drops that likely result from margin falls

    for t in range(numtimes):
        countsMatrix[t, marginmask] = 0.0

    return countsMatrix


@enable_xarray_wrapper
def rainonlyQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes particles that are probably
       not raindrops"""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    # Remove particles that are probably not rain
    for t in range(numtimes):
        masktimes[t, :] = rainonlymask

    countsMatrix = ma.masked_array(countsMatrix, mask=masktimes)

    return countsMatrix  # noqa: RET504


@enable_xarray_wrapper
def hailonlyQC(countsMatrix, returnmasked=True):
    """Based on Katja Friedrich's IDL QC subroutine. Removes particles that are probably
       not hail. Also returns number of particles remaining"""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    # Remove particles that are probably not hail
    for t in range(numtimes):
        masktimes[t, :] = hailonlymask

    masked = ma.masked_array(countsMatrix, mask=masktimes)
    if returnmasked:
        countsMatrix = masked
    total = masked.sum(axis=2).sum(axis=1)
    return countsMatrix, total


@enable_xarray_wrapper
def rainfallspeedQC(countsMatrix, fallspeedmask):
    """Removes all drops fall speeds +/- tolerance relative to rain fall speed relation."""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    for t in range(numtimes):
        masktimes[t, :] = fallspeedmask

    countsMatrix = ma.masked_array(countsMatrix, mask=masktimes)

    return countsMatrix  # noqa: RET504


@enable_xarray_wrapper
def maskhighdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""

    numtimes = np.size(countsMatrix, axis=0)

    diamindex = np.where(avg_diameter > rain_QC_params['high_diam_thresh'])[0][0]
    mask = np.zeros((numtimes, 32, 32))

    mask[:, :, diamindex:] = 1
    countsMatrix = ma.masked_array(countsMatrix, mask=mask)

    return countsMatrix  # noqa: RET504


@enable_xarray_wrapper
def masklowdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""

    numtimes = np.size(countsMatrix, axis=0)

    diamindex = np.where(avg_diameter > rain_QC_params['low_diam_thresh'])[0][0]
    mask = np.zeros((numtimes, 32, 32))

    mask[:, :, :diamindex] = 1
    countsMatrix = ma.masked_array(countsMatrix, mask=mask)

    return countsMatrix  # noqa: RET504


def get_qr_mask(qr_thresh, qr_bin):
    """[summary]

    Parameters
    ----------
    qr_thresh : [type]
        [description]
    qr_bin : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return np.where(qr_bin > qr_thresh)


def get_lowcount_mask(pcount_thresh, rainrate_thresh, pcounts_bin, rainrate_bin):
    """[summary]

    Parameters
    ----------
    pcount_thresh : [type]
        [description]
    rainrate_thresh : [type]
        [description]
    pcounts_bin : [type]
        [description]
    rainrate_bin : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return np.where((pcounts_bin < pcount_thresh) | (rainrate_bin < rainrate_thresh))
