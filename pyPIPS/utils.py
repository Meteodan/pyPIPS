"""[summary]

Returns
-------
[type]
    [description]
"""
import sys
import functools
import imp
from datetime import datetime
import numpy as np
import matplotlib.dates as dates
import xarray as xr

logdest = sys.stdout


class Bunch(object):
    """[summary]

    Parameters
    ----------
    object : [type]
        [description]
    """
    def __init__(self, adict):
        self.__dict__.update(adict)


# Modified from https://stackoverflow.com/questions/52043669/
# decorators-to-pass-numpy-arrays-or-xarray-arrays-to-functions
# Also see https://realpython.com/primer-on-python-decorators/
# FIXME: this doesn't work. reverting to original form
# def enable_xarray_wrapper(_func=None, *, output_core_dims=(())):
#     def _enable_xarray_wrapper(func):
#         """Adds an xarray wrapper for a function without core dimensions."""
#         @functools.wraps(func)
#         def wrapper(*args, **kwargs):
#             return xr.apply_ufunc(func, output_core_dims=output_core_dims, *args, kwargs=kwargs)
#         return wrapper

#     if _func is None:
#         return _enable_xarray_wrapper
#     else:
#         return _enable_xarray_wrapper(_func)
def enable_xarray_wrapper(func):
    """Adds an xarray wrapper for a function without core dimensions."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return xr.apply_ufunc(func, *args, kwargs=kwargs)
    return wrapper


def log(msg):
    """[summary]

    Parameters
    ----------
    msg : [type]
        [description]
    """
    logdest.write("%s\n" % msg)
    return


def warning(msg):
    """[summary]

    Parameters
    ----------
    msg : [type]
        [description]
    """
    logdest.write("Warning: %s\n" % msg)
    return


def fatal(msg):
    """[summary]

    Parameters
    ----------
    msg : [type]
        [description]
    """
    logdest.write("Error: %s\n" % msg)
    sys.exit()


def interpnan1D(a):
    """Replaces NaN's in a 1D array by interpolating from good values on either side"""
    ind = np.where(~np.isnan(a))[0]  # indices of valid values
    # Use valid values to interpolate to invalid values
    return np.interp(list(range(len(a))), ind, a[ind])


def trymax(a, default=0):
    """Tries to take a maximum of an array and returns 0 if the array is empty"""
    try:
        return a.max()
    except BaseException:
        return default


def import_all_from(module_path):
    """Modified from
       http://grokbase.com/t/python/python-list/1172ahxp0s/from-module-import-using-import
       Loads python file at "module_path" as module and adds contents to global namespace."""
    mod = imp.load_source('mod', module_path)
    return mod


# Below from https://stackoverflow.com/questions/47269390/
# numpy-how-to-find-first-non-zero-value-in-every-column-of-a-numpy-array
def first_nonzero(arr, axis, invalid_val=-1):
    """[summary]

    Parameters
    ----------
    arr : [type]
        [description]
    axis : [type]
        [description]
    invalid_val : int, optional
        [description], by default -1

    Returns
    -------
    [type]
        [description]
    """
    mask = (arr != 0)
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)


def last_nonzero(arr, axis, invalid_val=-1):
    """[summary]

    Parameters
    ----------
    arr : [type]
        [description]
    axis : [type]
        [description]
    invalid_val : int, optional
        [description], by default -1

    Returns
    -------
    [type]
        [description]
    """
    mask = (arr != 0)
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)


def readpyPIPSinput(path):
    """Reads and parses a pyPIPS text input file
       TODO: Rewrite this to use configparser?"""
    inputdict = {}
    with open(path, 'r') as inputfile:
        # Read and discard header
        inputfile.readline()
        inputfile.readline()
        inputfile.readline()
        inputdict['dis_dir'] = inputfile.readline().strip()
        inputdict['image_dir'] = inputfile.readline().strip()
        inputfile.readline()
        # Read in number of disdrometers
        inputdict['numdis'] = np.int(inputfile.readline().strip().split()[0])

        # Read in disdrometer locations and names
        dis_list = []
        # dis_ftype_list = []
        # tprh_filenames = []
        # wind_filenames = []
        dlocs = []
        dis_name_list = []
        starttimes = []
        stoptimes = []
        centertimes = []
        types = []
        for l in range(inputdict['numdis']):
            line = inputfile.readline().strip().split(',')
            dname = line[0]  # Disdrometer name
            dfile = line[1]  # Disdrometer filename
            starttime = line[2]  # np.int(line[2])
            stoptime = line[3]  # np.int(line[3])
            centertime = line[4]  # np.int(line[4])
            lat = np.float(line[5])
            lon = np.float(line[6])
            alt = np.float(line[7])
            try:
                type = line[8]
            except Exception:
                type = 'PIPS'
            # Append to appropriate lists
            dis_list.append(dfile)
            dis_name_list.append(dname)
            starttimes.append(starttime)
            stoptimes.append(stoptime)
            centertimes.append(centertime)
            dlocs.append((lat, lon, alt))
            types.append(type)

        # In the future, need to think about making a Disdrometer class to encapsulate
        # this stuff
        inputdict['dis_list'] = dis_list
        inputdict['dis_name_list'] = dis_name_list
        inputdict['starttimes'] = starttimes
        inputdict['stoptimes'] = stoptimes
        inputdict['centertimes'] = centertimes
        inputdict['dlocs'] = dlocs
        inputdict['type'] = types

        inputfile.readline()
        inputfile.readline()
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        inputdict['platform'] = line[0]
        inputdict['wavelength'] = np.float(line[1])
        inputfile.readline()
        inputdict['radar_dir'] = inputfile.readline().strip()
        inputfile.readline()
        inputdict['scattdir'] = inputfile.readline().strip()

        # Read in start and end times for the radar data analysis
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_int = list(map(int, line))
        inputdict['starttimerad'] = datetime(line_int[0], line_int[1], line_int[2], line_int[3],
                                             line_int[4], line_int[5])
        line = inputfile.readline().strip().split(',')
        line_int = list(map(int, line))
        inputdict['stoptimerad'] = datetime(line_int[0], line_int[1], line_int[2], line_int[3],
                                            line_int[4], line_int[5])

        # Read in plot window bounds
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_float = list(map(float, line))
        plotxmin = line_float[0]
        plotxmax = line_float[1]
        plotymin = line_float[2]
        plotymax = line_float[3]

        # Read in radar name,lat,lon
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        inputdict['radar_name'] = line[0]
        try:
            rlat = np.float(line[1])
        except BaseException:
            rlat = None
        inputdict['rlat'] = rlat
        try:
            rlon = np.float(line[2])
        except BaseException:
            rlon = None
        inputdict['rlon'] = rlon
        try:
            ralt = np.float(line[3])
        except BaseException:
            ralt = None
        inputdict['ralt'] = ralt
        try:
            el_req = np.float(line[4])
            print("requested elevation angle", el_req)
        except BaseException:
            el_req = 0.5    # Default to 0.5 degrees
        inputdict['el_req'] = el_req
        try:
            heading = np.float(line[5])
            print("Radar heading: ", heading)
        except BaseException:
            heading = None
        inputdict['heading'] = heading

        # Read in min range,max range, min azimuth, and max azimuth for radar
        # plotting (may deprecate this)
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_float = list(map(float, line))
        minrange = line_float[0]
        maxrange = line_float[1]
        minazim = line_float[2]
        maxazim = line_float[3]

        inputdict['radlims'] = [minrange, maxrange, minazim, maxazim]
        inputdict['plotlims'] = [plotxmin, plotxmax, plotymin, plotymax]

    return Bunch(inputdict)


def getTimeWindow(starttime, stoptime, timestampsnums):
    """Finds the start and end time indices of the data that best fit the requested times.
       This functionality is probably redundant now that we are converting things
       over to pandas and will probably be deprecated in the future."""

    if(np.int(starttime) == -1):
        startindex = 0
    else:
        starttime = dates.date2num(datetime(np.int(starttime[:4]), np.int(starttime[4:6]), np.int(
            starttime[6:8]), np.int(starttime[8:10]), np.int(starttime[10:12])))
        try:
            startindex = next(i for i, t in enumerate(
                timestampsnums) if t >= starttime)
        except Exception:
            startindex = 0

    if(np.int(stoptime) == -1):
        stopindex = np.size(timestampsnums) - 1
    else:
        stoptime = dates.date2num(datetime(np.int(stoptime[:4]), np.int(stoptime[4:6]), np.int(
            stoptime[6:8]), np.int(stoptime[8:10]), np.int(stoptime[10:12])))
        try:
            stopindex = next(i for i, t in enumerate(
                timestampsnums) if t >= stoptime)
        except Exception:
            stopindex = np.size(timestampsnums) - 1

    return startindex, stopindex


def mtokm(val, pos):
    """Convert m to km for formatting axes tick labels"""
    val = val / 1000.0
    return '%i' % val


def mtokmr1(val, pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val = val / 1000.0
    # return '{:2.{prec}f}'.format(val,prec=prec)
    return '%2.1f' % val


def mtokmr2(val, pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val = val / 1000.0
    # return '{:2.{prec}f}'.format(val,prec=prec)
    return '%2.2f' % val


def DDMtoDD(DDM, hem):
    """Converts from 'Degrees + Decimal Minutes' format to Decimal Degrees"""

    degrees = np.floor(DDM)
    DM = (DDM - degrees) * 100.

    if hem == 'N' or hem == 'E':
        sign = 1.
    else:
        sign = -1.

    return sign * (degrees + DM / 60.)


def interp_along_1D(ds, da1D, dim_to_interp, dim_fixed):
    """Interpolates 2D+ arrays in a Dataset ds to values given by a 1D DataArray da1D.
    The DataArray is dimensioned by "dim_fixed" and contains values to which to interpolate along
    "dim_to_interp" for each label of "dim_fixed".
    The result for each DataArray in ds (that is dimensioned by at least dim_fixed and
    dim_to_interp) is to reduce the dimension by one, yielding the interpolated values of that
    DataArray to the desired points in "dim_to_interp" for each label of "dim_fixed".

    Parameters
    ----------
    ds : xr.DataSet
        The DataSet containing the arrays to interpolate
    da1D : xr.DataArray
        1D DataArray containing the points to interpolate to along dimension "dim_to_interp".
        Dimensioned by dim_fixed.
    dim_to_interp : str
        Name of dimension along which to interpolate
    dim_fixed : str
        Name of dimension to hold fixed. The interpolation will be done along dim_to_interp for
        each value in da1D that corresponds to each label of dim_fixed.

    Returns
    -------
    xr.DataSet
        The DataSet containing the interpolated fields.
    """
    # First, rename the fixed dimension to 'temp' in da1D
    da1D = da1D.rename({dim_fixed: 'temp'})
    # Extract the associated coordinate
    temp_dim = da1D.coords['temp']

    # Interpolate ds along the fixed dimension using the values in da1D
    ds_interp = ds.interp({dim_to_interp: da1D, dim_fixed: temp_dim})
    # Drop the old "fixed_dim" coords, and then rename the "temp" coord to the original
    # name of the "fixed_dim" coords
    ds_interp = ds_interp.reset_coords(dim_fixed, drop=True)
    ds_interp = ds_interp.rename({'temp': dim_fixed})

    return ds_interp

