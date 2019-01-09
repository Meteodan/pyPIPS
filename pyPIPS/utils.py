import sys
import numpy as N
import imp
from datetime import datetime
import matplotlib.dates as dates

logdest = sys.stdout


class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)


def log(msg):
    logdest.write("%s\n" % msg)
    return


def warning(msg):
    logdest.write("Warning: %s\n" % msg)
    return


def fatal(msg):
    logdest.write("Error: %s\n" % msg)
    sys.exit()


def interpnan1D(a):
    """Replaces NaN's in a 1D array by interpolating from good values on either side"""
    ind = N.where(~N.isnan(a))[0]  # indices of valid values
    # Use valid values to interpolate to invalid values
    return N.interp(range(len(a)), ind, a[ind])


def trymax(a, default=0):
    """Tries to take a maximum of an array and returns 0 if the array is empty"""
    try:
        return a.max()
    except BaseException:
        return default


def import_all_from(module_path):
    """Modified from http://grokbase.com/t/python/python-list/1172ahxp0s/from-module-import-using-import
       Loads python file at "module_path" as module and adds contents to global namespace."""
    mod = imp.load_source('mod', module_path)
    return mod


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
        inputdict['numdis'] = N.int(inputfile.readline().strip().split()[0])

        # Read in disdrometer locations and names
        dis_list = []
        dis_ftype_list = []
        tprh_filenames = []
        wind_filenames = []
        dlocs = []
        dis_name_list = []
        starttimes = []
        stoptimes = []
        centertimes = []
        types = []
        for l in xrange(inputdict['numdis']):
            line = inputfile.readline().strip().split(',')
            dname = line[0]  # Disdrometer name
            dfile = line[1]  # Disdrometer filename
            starttime = line[2]  # N.int(line[2])
            stoptime = line[3]  # N.int(line[3])
            centertime = line[4]  # N.int(line[4])
            lat = N.float(line[5])
            lon = N.float(line[6])
            alt = N.float(line[7])
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
        inputdict['wavelength'] = N.float(line[1])
        inputfile.readline()
        inputdict['radar_dir'] = inputfile.readline().strip()
        inputfile.readline()
        inputdict['scattdir'] = inputfile.readline().strip()

        # Read in start and end times for the radar data analysis
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_int = map(int, line)
        inputdict['starttimerad'] = datetime(line_int[0], line_int[1], line_int[2], line_int[3],
                                             line_int[4], line_int[5])
        line = inputfile.readline().strip().split(',')
        line_int = map(int, line)
        inputdict['stoptimerad'] = datetime(line_int[0], line_int[1], line_int[2], line_int[3],
                                            line_int[4], line_int[5])

        # Read in plot window bounds
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_float = map(float, line)
        plotxmin = line_float[0]
        plotxmax = line_float[1]
        plotymin = line_float[2]
        plotymax = line_float[3]

        # Read in radar name,lat,lon
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        inputdict['radar_name'] = line[0]
        try:
            rlat = N.float(line[1])
        except BaseException:
            rlat = None
        inputdict['rlat'] = rlat
        try:
            rlon = N.float(line[2])
        except BaseException:
            rlon = None
        inputdict['rlon'] = rlon
        try:
            ralt = N.float(line[3])
        except BaseException:
            ralt = None
        inputdict['ralt'] = ralt
        try:
            el_req = N.float(line[4])
            print "requested elevation angle", el_req
        except BaseException:
            el_req = 0.5    # Default to 0.5 degrees
        inputdict['el_req'] = el_req
        try:
            heading = N.float(line[5])
            print "Radar heading: ", heading
        except BaseException:
            heading = None
        inputdict['heading'] = heading

        # Read in min range,max range, min azimuth, and max azimuth for radar
        # plotting (may deprecate this)
        inputfile.readline()
        line = inputfile.readline().strip().split(',')
        line_float = map(float, line)
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

    if(N.int(starttime) == -1):
        startindex = 0
    else:
        starttime = dates.date2num(datetime(N.int(starttime[:4]), N.int(starttime[4:6]), N.int(
            starttime[6:8]), N.int(starttime[8:10]), N.int(starttime[10:12])))
        try:
            startindex = next(i for i, t in enumerate(
                timestampsnums) if t >= starttime)
        except Exception:
            startindex = 0

    if(N.int(stoptime) == -1):
        stopindex = N.size(timestampsnums) - 1
    else:
        stoptime = dates.date2num(datetime(N.int(stoptime[:4]), N.int(stoptime[4:6]), N.int(
            stoptime[6:8]), N.int(stoptime[8:10]), N.int(stoptime[10:12])))
        try:
            stopindex = next(i for i, t in enumerate(
                timestampsnums) if t >= stoptime)
        except Exception:
            stopindex = N.size(timestampsnums) - 1

    return startindex, stopindex
