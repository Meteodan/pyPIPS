# disdrometer_module.py
# A collection of functions for reading VORTEX 2 disdrometer data

import netCDF4 as netcdf
import Nio
import numpy as np
from numpy import ma as ma
import matplotlib.dates as dates
from datetime import datetime, timedelta
import os
import shlex
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from . import obanmodule as oban
from . import radarmodule as radar
from . import thermolib as thermo
from . import utils
from scipy import special
from scipy.special import gammainc as gammap
from scipy.special import gammaln as gammln
import pandas as pd
import xarray as xr

deg2rad = np.pi / 180.

# ID #'s for PIPS 1A,1B,2A,2B
parsivel_ids = ['304545', '295153', '295166', '304543', '390654']

parsivel_names = {'304545': 'PIPS1A', '295153': 'PIPS1B', '295166': 'PIPS2A', '304543': 'PIPS2B',
                  '390654': 'TriPIPS'}

fieldnames = ['TIMESTAMP', 'RECORD', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms', 'WSDiag',
              'FastTemp', 'SlowTemp', 'RH', 'Pressure', 'FluxDirection', 'GPSTime',
              'GPSStatus', 'GPSLat', 'GPSLatHem', 'GPSLon', 'GPSLonHem', 'GPSSpd', 'GPSDir',
              'GPSDate', 'GPSMagVar', 'GPSAlt', 'WindDirAbs', 'Dewpoint', 'RHDer',
              'ParsivelStr']

fieldnames_TriPIPS = fieldnames[:]
fieldnames_TriPIPS.remove('FastTemp')
fieldnames_TriPIPS.remove('RECORD')

# Min diameter of bins (mm)
min_diameter_bins = [0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000, 1.125, 1.250,
                     1.500, 1.750, 2.000, 2.250, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000,
                     7.000, 8.000, 9.000, 10.000, 12.000, 14.000, 16.000, 18.000, 20.000, 23.000]
max_diameter_bins = [0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000, 1.125, 1.250, 1.500,
                     1.750, 2.000, 2.250, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 6.000, 7.000,
                     8.000, 9.000, 10.000, 12.000, 14.000, 16.000, 18.000, 20.000, 23.000, 26.000]

# Average diameter of bins (mm)
avg_diameter_bins = [0.5 * (x + y) for x, y in zip(min_diameter_bins, max_diameter_bins)]


fall_bins = [0.050, 0.150, 0.250, 0.350, 0.450, 0.550, 0.650, 0.750, 0.850, 0.950, 1.100, 1.300,
             1.500, 1.700, 1.900, 2.200, 2.600, 3.000, 3.400, 3.800, 4.400, 5.200, 6.000, 6.800,
             7.600, 8.800, 10.400, 12.000, 13.600, 15.200, 17.600, 20.800]

min_fall_bins = [0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900, 1.000, 1.200,
                 1.400, 1.600, 1.800, 2.000, 2.400, 2.800, 3.200, 3.600, 4.000, 4.800, 5.600, 6.400,
                 7.200, 8.000, 9.600, 11.200, 12.800, 14.400, 16.000, 19.200]

# Parsivel sampling area and sampling period
sensor_area = 0.0054    # (m^2)
sensor_width = 0.03     # (m)
sensor_length = 0.18    # (m)
sampling_period = 10.0  # (s)

min_diameter = np.array(min_diameter_bins)
avg_diameter = np.array(avg_diameter_bins)
max_diameter = np.array(max_diameter_bins)
bin_width = max_diameter - min_diameter
fall_bins = np.array(fall_bins)
min_fall_bins = np.array(min_fall_bins)

# Jaffrain and Berne (2011), Tokay et al. (2011)
eff_sensor_area = (180. * (30. - avg_diameter / 2.)) * 1.e-6

use_DD_QC = False   # Use my QC methods (see below)

# If True, use the raw measured fall speeds to compute number
# concentration
# If False, use the fall speed curve of Terry Schuur
use_measured_fs = True

use_strongwindQC = True      # Remove time records that are contaminated by strong wind?
use_splashingQC = False      # Remove drops that result from splashing?
use_marginQC = False         # Remove drops that result from margin falls?
use_rainonlyQC = False       # Remove all particles that are probably not rain?
use_hailonlyQC = False      # Remove all particles that are probably not hail?
use_graupelonlyQC = False   # Remove all particles that are probably not graupel?

# Remove all fall speeds +/- a percentage of the rain fall speed relation
use_rainfallspeedQC = False
# Percentage is set by falltol
falltol = 0.6
maskhigh = True
masklow = True

# True to mask out high diameter particles (given by threshold below), False to keep them
maskhighdiam = False
# True to mask out low diameter particles (given by threshold below), False to keep them
masklowdiam = False

highdiamthresh = 9.0
lowdiamthresh = 1.0

plot_QC = True        # True to plot fallspeed vs. diameter plot with QC information for each DSD
plot_splashingQC = False
plot_marginQC = False
plot_strongwindQC = False
plot_rainonlyQC = False
plot_rainfallspeedQC = False


def interpnan1D(a):
    """Replaces NaN's in a 1D array by interpolating from good values on either side"""
    ind = np.where(~np.isnan(a))[0]  # indices of valid values
    # Use valid values to interpolate to invalid values
    return np.interp(list(range(len(a))), ind, a[ind])


def DDMtoDD(DDM, hem):
    """Converts from 'Degrees + Decimal Minutes' format to Decimal Degrees"""

    degrees = np.floor(DDM)
    DM = (DDM - degrees) * 100.

    if(hem == 'N' or hem == 'E'):
        sign = 1.
    else:
        sign = -1.

    return sign * (degrees + DM / 60.)


def assignfallspeed(d, rhocorrect=False, rho=None):
    """Assigns a fall speed for a range of diameters based on code
       from David Dowell (originally from Terry Schuur).  It appears that
       the formulas originate from Atlas et al. (1973), but this took a bit of sleuthing!"""

    # Note, this appears to be valid at sea level.  For higher altitudes, a fall speed correction
    # should probably be applied based on Foote and duToit (1969): v = v0*(rho0/rho)^(0.4)
    # where rho0 = 1.204 kg/m^3 -- that corresponding to a T of 20 C and pressure of 1013 mb.

    v = np.where(d < 3.0, 3.78 * d**0.67, 9.65 - 10.3 * np.exp(-0.6 * d))
    v = np.where(d < 0.0, 0.0, v)

    # Correct fall speed based on air density
    if(rhocorrect and rho is not None):
        v = v[:, None] * (1.204 / rho)**(0.4)
        v = v.squeeze()
        v = np.atleast_1d(v)

    return v


# Create mask for splashing drops
bottom = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = [0, 1, 4, 7, 9, 11, 12, 13, 14, 14, 15, 16, 16, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
splashingmask = [[True if bottom[j] <= i <= top[j]
                  else False for i in range(32)] for j in range(32)]
splashingmask = np.array(splashingmask).T

# Create mask for margin falls
bottom = [0, 8, 14, 17, 20, 21, 22, 23, 24, 25, 26, 26, 27, 27,
          28, 28, 29, 29, 29, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = np.arange((32), dtype='int')
top[:] = 31
top[23:32] = 0
marginmask = [[True if bottom[j] <= i <= top[j] else False for i in range(32)] for j in range(32)]
marginmask = np.array(marginmask).T

# Create mask for non-raindrops
bottom = [0, 1, 4, 7, 9, 11, 12, 13, 14, 14, 15, 16, 16, 19, 19,
          20, 20, 21, 21, 21, 23, 24, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0]
top = [0, 8, 14, 17, 20, 21, 22, 23, 24, 25, 26, 26, 27, 27, 28,
       28, 29, 29, 29, 30, 30, 30, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0]
rainonlymask = [[False if bottom[j] < i < top[j] else True for i in range(32)] for j in range(32)]
rainonlymask = np.array(rainonlymask).T

# Create mask for non-hail
bottom = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0, 0, 0, 15, 15, 16, 16, 17, 17, 18, 19, 19, 20, 20, 20]
top = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 24, 25, 25, 31, 32, 32, 32, 32, 32, 32, 32, 32]
hailonlymask = [[False if bottom[j] < i < top[j] else True for i in range(32)] for j in range(32)]
hailonlymask = np.array(hailonlymask).T

# Create mask for all particles with fall speeds outside of fractional tolerance
rainvd = assignfallspeed(avg_diameter)
X, Y = np.meshgrid(rainvd, fall_bins)

if(maskhigh and masklow):
    fallspeedmask = np.where(np.abs((X - Y) / X) < falltol, False, True)
elif(masklow):     # Mask out speeds lower than tolerance
    fallspeedmask = np.where((Y - X) / X < -falltol, True, False)
elif(maskhigh):               # Mask out speeds higher than tolerance
    fallspeedmask = np.where((Y - X) / X > falltol, True, False)
else:
    fallspeedmask = None

# Create mask for strong wind conditions
strongwindmask = np.zeros((32, 32), dtype=bool)
strongwindmask[0:11, 20:32] = True


def truncatedspectrumQC(countsMatrix):
    """Masks out bad records where tokens have been set to -999 because the record
       was truncated"""
    countsMatrix = ma.masked_array(countsMatrix, mask=np.where(countsMatrix == -999, True, False))

    return countsMatrix


def strongwindQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine.  Removes drops affected by strong winds.
       Specifically, removes the entire 1-min DSD whenever there are large drops (> 3 mm) that have
       a low fall velocity (< 1 m/s)."""

    numtimes = np.size(countsMatrix, axis=0)
    flaggedtimes = []

    # Flag times that contain wind contamination
    for t in range(numtimes):
        baddrops = np.sum(countsMatrix[t, 0:11, 20:32])
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

        if(baddrops > 0):
            print("Severe Wind contamination, masking entire PSD!")
            countsMatrix[t, :] = -999.
            flaggedtimes.append(2)
        else:
            flaggedtimes.append(0)

    countsMatrix = ma.masked_array(countsMatrix, mask=np.where(countsMatrix == -999., True, False))

    return countsMatrix, flaggedtimes


def splashingQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from splashing"""

    numtimes = np.size(countsMatrix, axis=0)

    # Remove drops that likely result from splashing (use mask to index the count array)

    for t in range(numtimes):
        # countsMatrix[t,mask] = 0.0
        countsMatrix[t, :] = ma.masked_array(countsMatrix[t, :], mask=splashingmask)

    return countsMatrix


def marginQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from margin falls"""

    numtimes = np.size(countsMatrix, axis=0)

    # Remove drops that likely result from margin falls

    for t in range(numtimes):
        countsMatrix[t, marginmask] = 0.0

    return countsMatrix


def rainonlyQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes particles that are probably
       not raindrops"""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    # Remove particles that are probably not rain
    for t in range(numtimes):
        masktimes[t, :] = rainonlymask

    countsMatrix = ma.masked_array(countsMatrix, mask=masktimes)

    return countsMatrix


def hailonlyQC(countsMatrix, returnmasked=True, count=True):
    """Based on Katja Friedrich's IDL QC subroutine. Removes particles that are probably
       not hail. Also returns number of particles remaining"""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    # Remove particles that are probably not hail
    for t in range(numtimes):
        masktimes[t, :] = hailonlymask

    masked = ma.masked_array(countsMatrix, mask=masktimes)
    if(returnmasked):
        countsMatrix = masked
    total = masked.sum(axis=2).sum(axis=1)
    return countsMatrix, total


def rainfallspeedQC(countsMatrix, rainvd, falltol, maskhigh, masklow):
    """Removes all drops fall speeds +/- tolerance relative to rain fall speed relation."""

    numtimes = np.size(countsMatrix, axis=0)
    masktimes = np.zeros((numtimes, 32, 32), dtype=bool)

    for t in range(numtimes):
        masktimes[t, :] = fallspeedmask

    countsMatrix = ma.masked_array(countsMatrix, mask=masktimes)

    return countsMatrix


def maskhighdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""

    numtimes = np.size(countsMatrix, axis=0)

    diamindex = np.where(avg_diameter > highdiamthresh)[0][0]
    mask = np.zeros((numtimes, 32, 32))

    mask[:, :, diamindex:] = 1
    countsMatrix = ma.masked_array(countsMatrix, mask=mask)

    return countsMatrix


def masklowdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""

    numtimes = np.size(countsMatrix, axis=0)

    diamindex = np.where(avg_diameter > lowdiamthresh)[0][0]
    mask = np.zeros((numtimes, 32, 32))

    mask[:, :, :diamindex] = 1
    countsMatrix = ma.masked_array(countsMatrix, mask=mask)

    return countsMatrix


def correctPIPS(serialnum, infile, outfile):
    """Corrects offset Parsivel strings in a PIPS data file"""
    disfile = open(infile, 'r')

    lines_out = []
    parsivel_string_out_list = []
    parsivel_string_out = ''
    truncated = False
    first = True
    # ll = 0
    for line in disfile:
        line = line.rstrip('\r\n')
        tokens = line.strip().split(',')
        header_info = ",".join(tokens[:26])
        parsivel_string = tokens[26]
        parsivel_tokens = parsivel_string.strip().split(';')
        if(len(parsivel_tokens) > 1):
            if(truncated):  # Previous record was truncated
                # Look for the parsivel serial number somewhere in the next record
                # and find the index if it exists
                # print parsivel_tokens
                try:
                    line_out = header_info
                    sindex = [i for i, x in enumerate(parsivel_tokens) if serialnum in x][0]
                    # sindex = parsivel_tokens.index(serialnum)
                    # print sindex,parsivel_tokens[sindex]
                    # Concatenate portion of string before serial number onto end of previous
                    # record
                    parsivel_string_out = parsivel_string_out + ";".join(parsivel_tokens[:sindex])
                    # print parsivel_string_out
                    lines_out.append(line_out + ',' + parsivel_string_out)
                    parsivel_string_out_list.append(parsivel_string_out)
                    parsivel_string_out = ";".join(parsivel_tokens[sindex:]).lstrip()
                    truncated = False
                except BaseException:
                    print("Something went wrong!")
            elif(not first):
                # Should be able to just rip through the rest of the records and piece
                # them together
                sindex = [i for i, x in enumerate(parsivel_tokens) if serialnum in x][0]
                parsivel_string_out = parsivel_string_out + ";".join(parsivel_tokens[:sindex])
                parsivel_string_out_list.append(parsivel_string_out)
                lines_out.append(line_out + ',' + parsivel_string_out)
                line_out = header_info
                parsivel_string_out = ";".join(parsivel_tokens[sindex:]).lstrip()
            elif(first):
                if(len(parsivel_tokens) < 1036):    # Likely a truncated record
                    truncated = True
                    first = False
                    parsivel_string_out = parsivel_string
                    line_out = header_info
                else:
                    lines_out.append(line)
        else:
            lines_out.append(line)

    # Sort the output lines by record #
    sorted_lines_out = sorted(lines_out, key=lambda record: int(record.strip().split(',')[1]))
    # sorted_lines_out = lines_out

    outdisfile = open(outfile, 'w')
    for line in sorted_lines_out:
        outdisfile.write(line + '\n')


def readPIPS(filename, fixGPS=True, basicqc=False, rainfallqc=False, rainonlyqc=False,
             hailonlyqc=False, strongwindqc=False, detecthail=True, requested_interval=10.0,
             starttime=None, stoptime=None, tripips=False):
    """Reads data from Purdue-OU PIPS"""

    # Figure out which version we are reading in

    if tripips:
        curfieldnames = fieldnames_TriPIPS
    else:
        curfieldnames = fieldnames

    pdatetimes = []
    intensities = []
    preciptots = []
    reflectivities = []
    sampleintervals = []
    pcounts = []
    pcounts2 = []
    sensortemps = []
    amplitudes = []
    voltages = []
    pvoltages = []

    # Air density at the PSD times (needed for rain fall speed curve correction)
    rhoatpdatetimes = []

    # Will contain the number of drops in the diameter vs. fall speed matrix for each time
    # i.e. a (numtimes,32,32) array
    countsMatrix = []
    # Will contain computed number concentrations for each size bin for each time
    # i.e. a (numtimes,32) array
    concentrations = []

    onedrop_concentrations = []

    # dates = []
    # times = []
    datetimes = []
    windspds = []
    winddirrels = []
    winddirabss = []
    winddiags = []
    fasttemps = []
    slowtemps = []
    dewpoints = []
    RHs_derived = []
    RHs = []
    pressures = []
    pts = []
    qvs = []
    rhos = []
    compass_dirs = []
    GPS_lats = []
    GPS_lons = []
    GPS_stats = []
    GPS_alts = []

    disfile = open(filename, 'r')

    firstgoodGPS = False

    for line in disfile:
        tokens = line.strip().split(',')
        # Check for header line (older versions don't have it)
        if(tokens[0] == 'TIMESTAMP'):
            continue

        timestamp = tokens[curfieldnames.index('TIMESTAMP')]
        timestring = timestamp.strip().split()
        date = timestring[0]  # .strip('-')
        time = timestring[1]  # .strip(':')

        # Construct datetime object
        year = np.int(date[:4])
        month = np.int(date[5:7])
        day = np.int(date[8:])
        hour = np.int(time[:2])
        min = np.int(time[3:5])
        sec = np.int(time[6:8])

        datetimelogger = datetime(year, month, day, hour, min, sec)

        # recordnum = np.int(tokens[1])
        voltage = np.float(tokens[curfieldnames.index('BattV')])
        # paneltemp = np.float(tokens[3])
        winddirrel = np.float(tokens[curfieldnames.index('WindDir')])
        windspd = np.float(tokens[curfieldnames.index('WS_ms')])
        winddiag = np.float(tokens[curfieldnames.index('WSDiag')])
        if not tripips:
            fasttemp = np.float(tokens[curfieldnames.index('FastTemp')])
        slowtemp = np.float(tokens[curfieldnames.index('SlowTemp')])
        if tripips:
            fasttemp = slowtemp
        RH = np.float(tokens[curfieldnames.index('RH')])
        pressure = np.float(tokens[curfieldnames.index('Pressure')])
        compass_dir = np.float(tokens[curfieldnames.index('FluxDirection')])
        GPS_time = tokens[curfieldnames.index('GPSTime')]
        GPS_status = tokens[curfieldnames.index('GPSStatus')]
        GPS_lat = np.float(tokens[curfieldnames.index('GPSLat')])
        GPS_lat_hem = tokens[curfieldnames.index('GPSLatHem')]
        GPS_lat = DDMtoDD(GPS_lat, GPS_lat_hem)
        GPS_lon = np.float(tokens[curfieldnames.index('GPSLon')])
        GPS_lon_hem = tokens[curfieldnames.index('GPSLonHem')]
        GPS_lon = DDMtoDD(GPS_lon, GPS_lon_hem)
        # GPS_spd = np.float(tokens[18])
        # GPS_dir = np.float(tokens[19])
        GPS_date = tokens[curfieldnames.index('GPSDate')]
        # try:
        #     GPS_magvar = np.float(tokens[21])
        # except BaseException:
        #     GPS_magvar = np.nan
        try:
            GPS_alt = np.float(tokens[curfieldnames.index('GPSAlt')])
        except ValueError:
            GPS_alt = np.nan

        try:
            winddirabs = np.float(tokens[curfieldnames.index('WindDirAbs')])
            if np.isnan(winddirabs):
                winddirabs = winddirrel
        except ValueError:
            winddirabs = np.nan
        try:
            dewpoint = np.float(tokens[curfieldnames.index('Dewpoint')])
            if np.isnan(dewpoint):
                dewpoint = (thermo.calTdfromRH(pressure * 100., fasttemp + 273.15, RH / 100.) -
                            273.15)
        except ValueError:
            dewpoint = thermo.calTdfromRH(pressure * 100., fasttemp + 273.15, RH / 100.) - 273.15
        try:
            RH_derived = np.float(tokens[curfieldnames.index('RHDer')])
            if np.isnan(RH_derived):
                RH_derived = RH
        except ValueError:
            RH_derived = RH

        # Find the first good GPS time and date and use that
        # to construct the time offset for the data logger
        if not np.isnan(GPS_alt) and not firstgoodGPS and GPS_status == 'A':
            firstgoodGPS = True
            print(GPS_date, GPS_time)
            print(date, time)

            print(year, month, day, hour, min, sec)

            # Construct datetime object
            gyear = np.int('20' + GPS_date[4:])
            gmonth = np.int(GPS_date[2:4])
            gday = np.int(GPS_date[:2])
            ghour = np.int(GPS_time[:2])
            gmin = np.int(GPS_time[2:4])
            gsec = np.int(GPS_time[4:])

            datetimeGPS = datetime(gyear, gmonth, gday, ghour, gmin, gsec)
            GPS_offset = datetimeGPS - datetimelogger
            print(datetimeGPS, datetimelogger)
            print("GPS Offset", GPS_offset)

        datetimes.append(datetimelogger)
        windspds.append(windspd)
        winddirrels.append(winddirrel)
        winddirabss.append(winddirabs)
        winddiags.append(winddiag)
        fasttemps.append(fasttemp)
        slowtemps.append(slowtemp)
        dewpoints.append(dewpoint)
        RHs_derived.append(RH_derived)
        RHs.append(RH)
        pressures.append(pressure)
        compass_dirs.append(compass_dir)
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)
        voltages.append(voltage)

        pt = thermo.caltheta(pressure * 100., fasttemp + 273.15)
        qv = thermo.calqv(RH_derived / 100., pressure * 100., fasttemp + 273.15)
        rho = thermo.calrho(pressure * 100., pt, qv)

        pts.append(pt)
        qvs.append(qv)
        rhos.append(rho)

        parsivel_string = tokens[curfieldnames.index('ParsivelStr')]
        parsivel_tokens = parsivel_string.strip().split(';')
        serialnum = parsivel_tokens[0]
        if(serialnum in parsivel_ids and len(parsivel_tokens) >= 11):
            # parsivel_name = parsivel_names[serialnum]
            # print timestring
            precipintensity = np.float(parsivel_tokens[1])
            precipaccum = np.float(parsivel_tokens[2])
            parsivel_dBZ = np.float(parsivel_tokens[3])
            sample_interval = np.float(parsivel_tokens[4])
            signal_amplitude = np.float(parsivel_tokens[5])
            pcount = np.int(parsivel_tokens[6])
            sensor_temp = np.float(parsivel_tokens[7])
            pvoltage = np.float(parsivel_tokens[8])
            # sensor_time = parsivel_tokens[9]
            # sensor_date = parsivel_tokens[10]
            try:
                spectrum = [float(x) if x != '' else 0 for x in parsivel_tokens[11:]]
            except BaseException:
                spectrum = [-999 for i in range(1025)]
            # print "spectrum length = ",len(spectrum)
            if(len(spectrum) < 1024):
                print("Problem with Parsivel spectrum.  Flagging as bad!")
                print("Time: ", timestring)
#                 spectrum = [np.nan for i in xrange(1025)]
            else:
                if(len(spectrum) == 1025):
                    spectrum = spectrum[:-1]  # Strip off bogus last value

#             if(GPS_time):
#                 print GPS_date,GPS_time
#
            pdatetimes.append(datetimelogger)
            preciptots.append(precipaccum)
            intensities.append(precipintensity)
            reflectivities.append(parsivel_dBZ)
            sampleintervals.append(sample_interval)
            amplitudes.append(signal_amplitude)
            pcounts.append(pcount)
            sensortemps.append(sensor_temp)
            pvoltages.append(pvoltage)

#             pt = thermo.caltheta(pressure*100., fasttemp+273.15)
#             qv = thermo.calqv(RH_derived/100., pressure*100.,
#                               fasttemp+273.15)
#             rhoatpdatetimes.append(thermo.calrho(pressure*100., pt, qv))

            # Now create an array out of the spectrum and reshape it to 32x32
            spectrum = np.array(spectrum, dtype='int')
            if(spectrum.size == 1024):
                spectrum = spectrum.reshape((32, 32))
            else:
                spectrum = -999 * np.ones((32, 32), dtype='int')
                print(spectrum.size)
            # Append spectrum (corrected or not) to spectrum time list

            countsMatrix.append(spectrum)

    # Find indices in conventional time list that match those in the Parsivel time list
    pindices = np.searchsorted(datetimes, pdatetimes)
    rhoatpdatetimes = np.array(rhos)[pindices]

    # Recast countsMatrix as numpy array

    countsMatrix = np.dstack(countsMatrix)
    countsMatrix = np.rollaxis(countsMatrix, 2, 0)

    # Perform Katja's QC routines if desired (should not be used in combination with my methods
    # above, most are redundant anyway).

    X, Y = np.meshgrid(avg_diameter, fall_bins)
    flaggedtimes = np.zeros(len(pdatetimes), dtype=bool)
    # splashmask = np.zeros_like(countsMatrix)
    # marginmask = np.zeros_like(countsMatrix)

    countsMatrix = truncatedspectrumQC(countsMatrix)

    if(use_strongwindQC or basicqc or strongwindqc):
        countsMatrix, flaggedtimes = strongwindQC(countsMatrix)
        # if(detecthail and False):  # This doesn't seem to be working correctly, so disable for now
        #     # Unflag hail detection for those times where wind contamination was also detected
        #     hailflag = np.where(flaggedtimes > 0, False, hailflag)

    if(use_splashingQC or basicqc):
        countsMatrix = splashingQC(countsMatrix)

    if(use_marginQC or basicqc):
        countsMatrix = marginQC(countsMatrix)

    if(detecthail):
        countsMatrix, hailcounts = hailonlyQC(countsMatrix, returnmasked=False)
        hailflag = np.where(hailcounts > 0, True, False)

    # FIXME
    # if(use_rainfallspeedQC or rainfallqc):
    #     countsMatrix = rainfallspeedQC(countsMatrix, rainvd, falltol, maskhigh, masklow)

    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)

    if(use_rainonlyQC or rainonlyqc):
        countsMatrix = rainonlyQC(countsMatrix)

    if(use_hailonlyQC or hailonlyqc):
        countsMatrix, __ = hailonlyQC(countsMatrix)

    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)

    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)

    # Find total number of non-masked particles

    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)

    counts_1drop = np.ones_like(avg_diameter)

    # print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t, time in enumerate(pdatetimes):
        rainvd = assignfallspeed(avg_diameter, rhocorrect=True, rho=rhoatpdatetimes[t])
        spectrum = countsMatrix[t, :]
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy, vspectrum = np.meshgrid(avg_diameter, fall_bins)
            dspectrum = spectrum
        else:
            vspectrum = rainvd
            dspectrum = spectrum.sum(axis=0)   # Sum up particles for each diameter bin

            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]

        # Now compute the number concentration using the assumed fall speeds, sensor area, and
        # sampling interval
        # Units are #/m^3/mm
#         if time == '171320':
#             print spectrum.size,spectrum,countsMatrix.data[t,:],spectrum.data
        if(spectrum.size == 1024 and flaggedtimes[t] < 2):
            concentration = dspectrum / \
                (vspectrum * sampleintervals[t] * eff_sensor_area * (max_diameter - min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
#                 onedropvel = np.ma.average(vspectrum,axis=0,weights=dspectrum)
#             else:
#                 onedropvel = rainvd
            onedropvel = rainvd
            onedrop_concentration = counts_1drop / \
                (onedropvel * sampleintervals[t] * eff_sensor_area * (max_diameter - min_diameter))
            # print "concentration.shape"
            # print concentration.shape
        elif(flaggedtimes[t] < 2):
            concentration = np.zeros_like(avg_diameter)
        else:
            concentration = -999. * np.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999., concentration)

        # Throw out particles above and below a certain diameter if desired

#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)

        # print "Number of particles counted vs. summed number: "
        # print np.int(line[7]),np.sum(dspectrum)

        # print concentration

        concentrations.append(concentration)
        onedrop_concentrations.append(onedrop_concentration)

    concentrations = ma.array(concentrations)
    onedrop_concentrations = ma.array(onedrop_concentrations)
    # print "concentrations: ",concentrations
    pcounts = ma.array(pcounts)
    # pcounts2 = np.array(pcounts2)
    amplitudes = np.array(amplitudes)

    # Correct the logger time and date using the GPS time

    # Sometimes the GPS has major issues for the whole deployment.
    # In this case, set the offset to 0
    if(not firstgoodGPS):
        GPS_offset = timedelta(seconds=0)

    datetimes_corrected = []
    pdatetimes_corrected = []
    for datetimelogger in datetimes:
        datetimes_corrected.append(datetimelogger + GPS_offset)
    for pdatetimelogger in pdatetimes:
        pdatetimes_corrected.append(pdatetimelogger + GPS_offset)

    # Create a dataframe out of the 1-s data
    # First create a dictionary to organize the 1-s data
    convdict = {'windspd': windspds, 'winddirrel': winddirrels, 'winddirabs': winddirabss,
                'winddiag': winddiags, 'fasttemp': fasttemps, 'slowtemp': slowtemps,
                'RH': RHs, 'RH_derived': RHs_derived, 'pressure': pressures,
                'compassdir': compass_dirs, 'GPS_lat': GPS_lats, 'GPS_lon': GPS_lons,
                'GPS_stat': GPS_stats, 'GPS_alt': GPS_alts, 'voltage': voltages,
                'dewpoint': dewpoints, 'pt': pts, 'qv': qvs, 'rho': rhos}

    conv_df = pd.DataFrame(convdict, index=datetimes_corrected)

    # Create one dataframe for N(D), one for the "one-drop" N(D), and one for all the
    # other 10-s data

    ND_df = pd.DataFrame(data=concentrations, index=pdatetimes_corrected, columns=avg_diameter)
    ND_onedrop_df = pd.DataFrame(data=onedrop_concentrations, index=pdatetimes_corrected,
                                 columns=avg_diameter)

    # Create a dictionary to organize the 10-s data
    PSDdict = {'intensity': intensities, 'preciptot': preciptots, 'reflectivity': reflectivities,
               'pcount': pcounts, 'pcount2': pcounts2, 'amplitude': amplitudes,
               'flaggedtimes': flaggedtimes, 'hailflag': hailflag}

    PSD_df = pd.DataFrame(PSDdict, index=pdatetimes_corrected)

    # For the countsMatrix, need to use xarray
    countsMatrix_da = xr.DataArray(countsMatrix, coords=[pdatetimes_corrected, avg_diameter,
                                                         fall_bins],
                                   dims=['time', 'diameter', 'velocity'])

    # Average and thin the DSD data with the desired interval

    DSD_interval, intervalstr, ND_df, ND_onedrop_df, PSD_df, countsMatrix_da = \
        resamplePSD(
            requested_interval,
            sampleintervals[0],
            ND_df,
            ND_onedrop_df,
            PSD_df,
            countsMatrix_da)

    # Pandas apparently (gotcha!) converts missing values to NaN when extracting the numpy array
    # representation using .values. Since I want the original masked array functionality for now for
    # further computations, I need to remask the array here.
    # Otherwise, the NaN's propagate in further computations...
    # In the future, another solution that uses Pandas more natively should be
    # pursued, but this will work for now
    ND = ND_df.values
    mask = ND_df.isnull()
    ND = ma.array(ND, mask=mask)
    ND_onedrop = ND_onedrop_df.values
    mask = ND_onedrop_df.isnull()
    ND_onedrop = ma.array(ND_onedrop, mask=mask)
    # Argh, have to convert back to datetime objects.  This one from
    # http://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
    datetimes_corrected = conv_df.index.to_pydatetime()
    pdatetimes_corrected = ND_df.index.to_pydatetime()
    DSD_index = ND_df.index
    countsMatrix = countsMatrix_da.values
    mask = countsMatrix_da.isnull()
    countsMatrix = ma.array(countsMatrix, mask=mask)

    # Aggregate everything into a dictionary to return

    PIPS_dict = {'convtimestamps': datetimes_corrected, 'PSDtimestamps': pdatetimes_corrected,
                 'DSD_index': DSD_index, 'DSD_interval': DSD_interval, 'intervalstr': intervalstr,
                 'conv_df': conv_df, 'PSD_df': PSD_df, 'ND_df': ND_df,
                 'ND_onedrop_df': ND_onedrop_df, 'countsMatrix': countsMatrix, 'ND': ND,
                 'ND_onedrop': ND_onedrop}

    return PIPS_dict


def readCU(filename_conv, filename_dsd, fixGPS=True, basicqc=False, rainfallqc=False,
           rainonlyqc=True, hailonlyqc=False, strongwindqc=False, detecthail=True,
           requested_interval=60.0, starttime=None, stoptime=None, matchtimes=True):
    """Reads data from the CU VORTEX-2 disdrometer probes. Currently only works for 2009 data"""

    if starttime is not None:
        starttimestamp = datetime(np.int(starttime[:4]), np.int(starttime[4:6]), np.int(
            starttime[6:8]), np.int(starttime[8:10]), np.int(starttime[10:12]))
    else:
        starttimestamp = None
    if stoptime is not None:
        stoptimestamp = datetime(np.int(stoptime[:4]), np.int(stoptime[4:6]), np.int(
            stoptime[6:8]), np.int(stoptime[8:10]), np.int(stoptime[10:12]))
    else:
        stoptimestamp = None

    pdatetimes = []
    intensities = []
    preciptots = []
    reflectivities = []
    # sampleintervals = []
    pcounts = []
    pcounts2 = []
    sensortemps = []
    # amplitudes = []
    voltages = []
    # pvoltages = []
    weathercodes = []
    visibilities = []

    # Air density at the PSD times (needed for rain fall speed curve correction)
    rhoatpdatetimes = []

    # Will contain the number of drops in the diameter vs. fall speed matrix for each time
    # i.e. a (numtimes,32,32) array
    countsMatrix = []
    # Will contain computed number concentrations for each size bin for each time
    # i.e. a (numtimes,32) array
    concentrations = []

    onedrop_concentrations = []

    # dates = []
    # times = []
    datetimes = []
    bwindspds = []
    bwinddirabss = []
    swindspds = []
    swinddirabss = []
    swinddiags = []
    # fasttemps = []
    slowtemps = []
    dewpoints = []
    pts = []
    qvs = []
    rhos = []
    # RHs_derived = []
    RHs = []
    pressures = []
    # compass_dirs = []
    GPS_lats = []
    GPS_lons = []
    GPS_stats = []
    GPS_alts = []
    GPS_speeds = []

    # First, read in the disdrometer data

    disfile = open(filename_dsd, 'r')

    for line in disfile:
        tokens = line.strip().split(';')
        # print line
        date = tokens[0]
        time = tokens[1]

        # Construct datetime object
        day = int(date[0:2])
        month = int(date[3:5])
        year = int(date[6:10])
        hour = int(time[0:2])
        min = int(time[3:5])
        sec = int(time[6:8])

        pdatetimelogger = datetime(year, month, day, hour, min, sec)

        if starttimestamp is not None:
            if pdatetimelogger < starttimestamp:
                continue
        if stoptimestamp is not None:
            if pdatetimelogger > stoptimestamp:
                continue

        pdatetimes.append(pdatetimelogger)
        intensities.append(np.float(tokens[2]))
        preciptots.append(np.float(tokens[3]))
        weathercodes.append(np.int(tokens[4]))
        reflectivities.append(np.float(tokens[5]))
        visibilities.append(np.float(tokens[6]))
        pcounts.append(np.int(tokens[7]))
        sensortemps.append(np.float(line[8]))

        # spectra start with token position 10

        spectrum = [int(x) if x != '' else 0 for x in tokens[10:-1]]
        # Looks like the first value might be left off: prepend a 0 to the
        # beginning if the spectrum isn't empty
        if(len(spectrum) > 0):
            spectrum.insert(0, 0)

        # Now create an array out of the spectrum and reshape it to 32x32
        spectrum = np.array(spectrum)
        if(spectrum.size > 0):
            spectrum = spectrum.reshape((32, 32))
        else:
            spectrum = np.zeros((32, 32), dtype='int')

        # Append spectrum (corrected or not) to spectrum time list

        countsMatrix.append(spectrum)

    # Calculate sampling interval from successive timestamps in the Parsivel data file
    # TODO: need to come up with another solution: this could be wrong if there are missing records
    # at the beginning....
    sample_interval = np.float((pdatetimes[1] - pdatetimes[0]).seconds)

    disfile.close()
    # Next, read in the conventional data

    convfile = open(filename_conv, 'r')

    # Read in the header information

    dummy = convfile.readline()
    dummy = convfile.readline()
    dummy = convfile.readline()
    dummy = convfile.readline()

    for line in convfile:
        mysplit = shlex.shlex(line)
        mysplit.whitespace += ','
        mysplit.whitespace_split = True
        # line = l.strip().split(',')
        tokens = list(mysplit)
        datetimestring = tokens[0]

        # Construct datetime object
        year = int(datetimestring[1:5])
        month = int(datetimestring[6:8])
        day = int(datetimestring[9:11])
        hour = int(datetimestring[12:14])
        min = int(datetimestring[15:17])
        sec = int(datetimestring[18:20])

        datetimelogger = datetime(year, month, day, hour, min, sec)
        # Skip this record if it lies before or after the Parsivel file time period
        if matchtimes and (datetimelogger < pdatetimes[0] or datetimelogger > pdatetimes[-1]):
            continue
        # Skip this record if it lies before or after the desired period
        if starttimestamp is not None:
            if datetimelogger < starttimestamp:
                continue
        if stoptimestamp is not None:
            if datetimelogger > stoptimestamp:
                continue

        # recordnum = np.int(tokens[1])
        # Not sure what temperature sensor the CU probes are using. Assuming it's "slow response"
        # for now
        slowtemp = np.float(tokens[2])
        RH = np.float(tokens[3])
        bwindspd = np.float(tokens[4])
        bwinddirabs = np.float(tokens[5])
        try:
            swindspd = np.float(tokens[6])
        except Exception:
            swindspd = np.nan
        try:
            swinddirabs = np.float(tokens[7])
        except Exception:
            swinddirabs = np.nan
        try:
            swinddiag = np.float(tokens[8])
        except BaseException:
            swinddiag = np.nan
        voltage = np.float(tokens[9])
        try:
            GPS_status = np.int(tokens[15].strip('\"'))
        except Exception:
            GPS_status = 0
        try:
            GPS_alt = np.float(tokens[16].strip('\"'))
        except Exception:
            GPS_alt = np.nan
        try:
            GPS_lat = np.float(tokens[11].strip('\"'))
            GPS_lat = DDMtoDD(GPS_lat / 100., 'N')
        except Exception:
            GPS_lat = np.nan
            GPS_status = 0
        try:
            GPS_lon = np.float(tokens[12].strip('\"'))
            GPS_lon = DDMtoDD(GPS_lon / 100., 'W')
        except Exception:
            GPS_lon = np.nan
            GPS_status = 0
        try:
            GPS_speed = tokens[13].strip('\"')
            GPS_speed = np.float(GPS_speed)
        except Exception:
            GPS_speed = np.nan
        try:
            GPS_alt = np.float(tokens[16].strip('\"'))
        except Exception:
            GPS_alt = np.nan

        try:
            pressure = np.float(tokens[17].strip('\"'))
        except Exception:
            pressure = np.nan

        # Compute dewpoint
        try:
            dewpoint = thermo.calTdfromRH(pressure * 100., slowtemp + 273.15, RH / 100.) - 273.15
        except BaseException:
            # Can actually calculate dewpoint without pressure but we'll not mess with that
            # for now
            dewpoint = np.nan

        pt = thermo.caltheta(pressure * 100., slowtemp + 273.15)
        qv = thermo.calqv(RH / 100., pressure * 100., slowtemp + 273.15)

        pts.append(pt)
        qvs.append(qv)
        rhos.append(thermo.calrho(pressure * 100., pt, qv))

        # Append values to their respective lists
        datetimes.append(datetimelogger)
        bwindspds.append(bwindspd)
        bwinddirabss.append(bwinddirabs)
        swindspds.append(swindspd)
        swinddirabss.append(swinddirabs)
        swinddiags.append(swinddiag)
        slowtemps.append(slowtemp)
        dewpoints.append(dewpoint)
        RHs.append(RH)
        pressures.append(pressure)
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_speeds.append(GPS_speed)
        GPS_alts.append(GPS_alt)
        GPS_stats.append(GPS_status)
        voltages.append(voltage)

    convfile.close()

    # Find indices in conventional time list that match those in the Parsivel time list
    pindices = np.searchsorted(datetimes, pdatetimes)
    rhoatpdatetimes = np.array(rhos)[pindices]

    # Recast countsMatrix as numpy array

    countsMatrix = np.dstack(countsMatrix)
    countsMatrix = np.rollaxis(countsMatrix, 2, 0)

    # Perform Katja's QC routines if desired (should not be used in combination with my methods
    # above, most are redundant anyway).

    X, Y = np.meshgrid(avg_diameter, fall_bins)
    flaggedtimes = np.zeros(len(pdatetimes), dtype=bool)
    # splashmask = np.zeros_like(countsMatrix)
    # marginmask = np.zeros_like(countsMatrix)

    countsMatrix = truncatedspectrumQC(countsMatrix)

    if(use_strongwindQC or basicqc or strongwindqc):
        countsMatrix, flaggedtimes = strongwindQC(countsMatrix)
        # if(detecthail and False):  # This doesn't seem to be working correctly, so disable for now
        #     # Unflag hail detection for those times where wind contamination was also detected
        #     hailflag = np.where(flaggedtimes > 0, False, hailflag)

    if(use_splashingQC or basicqc):
        countsMatrix = splashingQC(countsMatrix)

    if(use_marginQC or basicqc):
        countsMatrix = marginQC(countsMatrix)

    if(detecthail):
        countsMatrix, hailcounts = hailonlyQC(countsMatrix, returnmasked=False)
        hailflag = np.where(hailcounts > 0, True, False)

    # FIXME
    # if(use_rainfallspeedQC or rainfallqc):
    #     countsMatrix = rainfallspeedQC(countsMatrix, rainvd, falltol, maskhigh, masklow)

    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)

    if(use_rainonlyQC or rainonlyqc):
        countsMatrix = rainonlyQC(countsMatrix)

    if(use_hailonlyQC or hailonlyqc):
        countsMatrix, __ = hailonlyQC(countsMatrix)

    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)

    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)

    # Find total number of non-masked particles

    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)

    counts_1drop = np.ones_like(avg_diameter)

    # print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t, time in enumerate(pdatetimes):
        rainvd = assignfallspeed(avg_diameter, rhocorrect=True, rho=rhoatpdatetimes[t])
        spectrum = countsMatrix[t, :]
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy, vspectrum = np.meshgrid(avg_diameter, fall_bins)
            dspectrum = spectrum
        else:
            vspectrum = rainvd
            dspectrum = spectrum.sum(axis=0)   # Sum up particles for each diameter bin

            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]

        # Now compute the number concentration using the assumed fall speeds, sensor area, and
        # sampling interval
        # Units are #/m^3/mm
#         if time == '171320':
#             print spectrum.size,spectrum,countsMatrix.data[t,:],spectrum.data
        onedropvel = rainvd
        onedrop_concentration = counts_1drop / \
            (onedropvel * sample_interval * eff_sensor_area * (max_diameter - min_diameter))

        if(spectrum.size == 1024 and flaggedtimes[t] < 2):
            concentration = dspectrum / \
                (vspectrum * sample_interval * eff_sensor_area * (max_diameter - min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
#                 onedropvel = np.ma.average(vspectrum,axis=0,weights=dspectrum)
#             else:
#                 onedropvel = rainvd
            # print "concentration.shape"
            # print concentration.shape
        elif(flaggedtimes[t] < 2):
            concentration = np.zeros_like(avg_diameter)
        else:
            concentration = -999. * np.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999., concentration)

        # Throw out particles above and below a certain diameter if desired

#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)

        # print "Number of particles counted vs. summed number: "
        # print np.int(line[7]),np.sum(dspectrum)

        # print concentration

        concentrations.append(concentration)
        onedrop_concentrations.append(onedrop_concentration)

    concentrations = ma.array(concentrations)
    onedrop_concentrations = ma.array(onedrop_concentrations)
    # print "concentrations: ",concentrations
    pcounts = ma.array(pcounts)
    # pcounts2 = np.array(pcounts2)

    # Correct the logger time and date using the GPS time (commented out for now since
    # probably not needed for CU data)

    # Sometimes the GPS has major issues for the whole deployment.
    # In this case, set the offset to 0
#     if(not firstgoodGPS):
#         GPS_offset = timedelta(seconds=0)
#
#     datetimes_corrected = []
#     pdatetimes_corrected = []
#     for datetimelogger in datetimes:
#         datetimes_corrected.append(datetimelogger + GPS_offset)
#     for pdatetimelogger in pdatetimes:
#         pdatetimes_corrected.append(pdatetimelogger + GPS_offset)

    datetimes_corrected = pd.Index(datetimes)
    pdatetimes_corrected = pd.Index(pdatetimes)

    # Create a dataframe out of the 1-s data
    # First create a dictionary to organize the 1-s data
    convdict = {'bwindspd': bwindspds, 'bwinddirabs': bwinddirabss, 'swindspd': swindspds,
                'swinddirabs': swinddirabss, 'swinddiag': swinddiags, 'slowtemp': slowtemps,
                'RH': RHs, 'pressure': pressures, 'GPS_lat': GPS_lats, 'GPS_lon': GPS_lons,
                'GPS_speed': GPS_speeds, 'GPS_alt': GPS_alts, 'voltage': voltages,
                'dewpoint': dewpoints, 'pt': pts, 'qv': qvs, 'rho': rhos}

    conv_df = pd.DataFrame(convdict, index=datetimes_corrected)

    # Create one dataframe for N(D), one for the "one-drop" N(D), and one for all the
    # other 10-s data

    ND_df = pd.DataFrame(data=concentrations, index=pdatetimes_corrected, columns=avg_diameter)
    ND_onedrop_df = pd.DataFrame(data=onedrop_concentrations, index=pdatetimes_corrected,
                                 columns=avg_diameter)

    # Create a dictionary to organize the 10-s data
    PSDdict = {'intensity': intensities, 'preciptot': preciptots, 'reflectivity': reflectivities,
               'pcount': pcounts, 'pcount2': pcounts2, 'flaggedtimes': flaggedtimes,
               'hailflag': hailflag}

    PSD_df = pd.DataFrame(PSDdict, index=pdatetimes_corrected)

    # For the countsMatrix, need to use xarray
    countsMatrix_da = xr.DataArray(countsMatrix, coords=[pdatetimes_corrected, avg_diameter,
                                                         fall_bins],
                                   dims=['time', 'diameter', 'velocity'])

    # Some of the files have duplicated times (though they seem to contain different data: need to
    # check with Katja about this). Remove them here.
    ND_df = ND_df[~ND_df.index.duplicated()]
    ND_onedrop_df = ND_onedrop_df[~ND_onedrop_df.index.duplicated()]
    PSD_df = PSD_df[~PSD_df.index.duplicated()]
    countsMatrix_da = countsMatrix_da.sel(time=~pdatetimes_corrected.duplicated())
    pdatetimes_corrected = pdatetimes_corrected.drop_duplicates()
    # If desired, clean up the Parsivel data so that the start and end times are within the
    # datetimes index. This option is recommended!
    if matchtimes:
        pdatetimes_corrected = pdatetimes_corrected.intersection(datetimes_corrected)
        ND_df = ND_df.reindex(pdatetimes_corrected)
        ND_onedrop_df = ND_onedrop_df.reindex(pdatetimes_corrected)
        PSD_df = PSD_df.reindex(pdatetimes_corrected)
        countsMatrix_da.reindex(time=pdatetimes_corrected)

    # Average and thin the DSD data with the desired interval

    DSD_interval, intervalstr, ND_df, ND_onedrop_df, PSD_df, countsMatrix_da = \
        resamplePSD(
            requested_interval,
            sample_interval,
            ND_df,
            ND_onedrop_df,
            PSD_df,
            countsMatrix_da)

    # Pandas apparently (gotcha!) converts missing values to NaN when extracting the numpy array
    # representation using .values. Since I want the original masked array functionality for now for
    # further computations, I need to remask the array here.
    # Otherwise, the NaN's propagate in further computations...
    # In the future, another solution that uses Pandas more natively should be
    # pursued, but this will work for now
    ND = ND_df.values
    mask = ND_df.isnull()
    ND = ma.array(ND, mask=mask)
    ND_onedrop = ND_onedrop_df.values
    mask = ND_onedrop_df.isnull()
    ND_onedrop = ma.array(ND_onedrop, mask=mask)
    # Argh, have to convert back to datetime objects.  This one from
    # http://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
    datetimes_corrected = datetimes_corrected.to_pydatetime()
    pdatetimes_corrected = ND_df.index.to_pydatetime()
    DSD_index = ND_df.index
    countsMatrix = countsMatrix_da.values
    mask = countsMatrix_da.isnull()
    countsMatrix = ma.array(countsMatrix, mask=mask)

    # Aggregate everything into a dictionary to return

    PIPS_dict = {'convtimestamps': datetimes_corrected, 'PSDtimestamps': pdatetimes_corrected,
                 'DSD_index': DSD_index, 'DSD_interval': DSD_interval, 'intervalstr': intervalstr,
                 'conv_df': conv_df, 'PSD_df': PSD_df, 'ND_df': ND_df,
                 'ND_onedrop_df': ND_onedrop_df, 'countsMatrix': countsMatrix, 'ND': ND,
                 'ND_onedrop': ND_onedrop}

    return PIPS_dict


def resampleconv(type, interval, sec_offset, conv_df, gusts=False, gustintvstr='3S', center=False):
    """Resamples the conventional data to a longer interval"""

    if type == 'PIPS':
        winddirkey = 'winddirabs'
        windspdkey = 'windspd'
        other_data = ['fasttemp', 'slowtemp', 'RH', 'RH_derived', 'pressure',
                      'GPS_lat', 'GPS_lon', 'GPS_alt', 'voltage', 'dewpoint', 'rho']
    elif type == 'CU':
        winddirkey = 'bwinddirabs'
        windspdkey = 'bwindspd'
        other_data = ['slowtemp', 'RH', 'pressure',
                      'GPS_lat', 'GPS_lon', 'GPS_alt', 'voltage', 'dewpoint', 'rho']
    elif type == 'NV2':
        winddirkey = 'swinddirabs'
        windspdkey = 'swindspd'
        other_data = ['slowtemp', 'RH', 'pressure', 'dewpoint', 'rho']

    intervalstr = '{:d}S'.format(int(interval))

    # First, resample the winds

    conv_resampled_dict = resamplewind(conv_df.index, sec_offset, conv_df[winddirkey],
                                       conv_df[windspdkey], intervalstr, gusts=gusts,
                                       gustintvstr=gustintvstr, center=center)
    conv_resampled_df = pd.DataFrame(conv_resampled_dict)

    # Special treatment for wind diagnostic flags
    # Note, have to use numpy max() as a lambda function because the
    # pandas resample.max() does not propagate NaN's!
    if type == 'PIPS':
        winddiags_rs = pd.Series(data=conv_df['winddiag'], index=conv_df.index).resample(
            intervalstr, label='right', closed='right',
            base=sec_offset, how=lambda x: utils.trymax(x.values))

        conv_resampled_df['winddiag'] = winddiags_rs

    # Now, resample the other one-sec data
    other_resampled_df = \
        conv_df[other_data].resample(intervalstr, label='right', closed='right',
                                     base=sec_offset).mean()

    conv_resampled_df = conv_resampled_df.join(other_resampled_df)
    return conv_resampled_df


def resamplePSD(requested_interval, actual_interval, ND_df, ND_onedrop_df, PSD_df,
                countsMatrix_da=None):
    """Resamples Parsivel data to a longer interval"""

    DSD_index_interval = int(requested_interval / actual_interval)
    print("Requested DSD interval: {:.1f}. Actual DSD interval: {:.1f}".format(
        requested_interval, DSD_index_interval * actual_interval))
    DSD_interval = DSD_index_interval * actual_interval
    intervalstr = '{:d}S'.format(int(DSD_interval))

    if(DSD_interval == actual_interval):
        return DSD_interval, intervalstr, ND_df, ND_onedrop_df, PSD_df, countsMatrix_da
    else:
        # We need to find the offset corresponding to the starting second and then
        # generate the frequency string. Seems like there should be an easier way...
        sec_offset = ND_df.index.to_pydatetime()[0].second

        # When resampling, fill in missing values with zeros: TEMPORARY FIX. later will deal
        # with missing values differently

        ND_df = ND_df.resample(
            intervalstr, label='right', closed='right', base=sec_offset).mean().fillna(0)
        ND_onedrop_df = ND_onedrop_df.resample(
            intervalstr, label='right', closed='right', base=sec_offset).mean().fillna(0)
        # Concentration in each bin assuming only one drop over the new interval
        ND_onedrop_df = ND_onedrop_df.applymap(lambda x: 10. * x / DSD_interval)

        PSD_df_rs = PSD_df.resample(intervalstr, label='right', closed='right',
                                    base=sec_offset)

        # Each column of the dataframe needs to downsample differently. For example, the
        # precipitation totals need to be summed, while the reflectivity should be averaged, etc.
        # We can do this with the Pandas "agg" function on the resampler. Neat, huh?
        PSD_df = PSD_df_rs.agg({'intensity': np.mean, 'preciptot': np.sum,
                                'reflectivity': np.mean, 'pcount': np.sum, 'pcount2': np.sum,
                                'amplitude': np.mean, 'flaggedtimes': np.any,
                                'hailflag': np.any}).fillna(0)

        # Use xarray for the counts matrix
        if countsMatrix_da is not None:
            countsMatrix_da = countsMatrix_da.resample(time=intervalstr, label='right',
                                                       closed='right',
                                                       base=sec_offset).sum(dim='time').fillna(0)
        # Keep the following code on standby until the above is tested!
#
#         # I'm not sure what's going on here. Sometimes I get a valueerror suggesting the data are
#         # empty, so the except clause is a temporary
#         # bandaid.
#         try:
#             flaggedtimes_df = flaggedtimes_df.resample(
#                 intervalstr,
#                 label='right',
#                 closed='right',
#                 base=sec_offset,
#                 how=lambda x: x.values.max())
#         except BaseException:
#             flaggedtimes_df = pd.Series(
#                 data=np.zeros_like(
#                     amplitudes_df.values),
#                 index=amplitudes_df.index)
#         try:
#             hailflag_df = hailflag_df.resample(
#                 intervalstr,
#                 label='right',
#                 closed='right',
#                 base=sec_offset,
#                 how=lambda x: x.values.max())
#         except BaseException:
#             hailflag_df = pd.Series(
#                 data=np.zeros_like(
#                     amplitudes_df.values),
#                 index=amplitudes_df.index)

        return DSD_interval, intervalstr, ND_df, ND_onedrop_df, PSD_df, countsMatrix_da


def readPIPSdf(filename, fixGPS=True, basicqc=False, rainfallqc=False,
               rainonlyqc=False, strongwindqc=False, DSD_interval=10.0):
    """Reads data from Purdue-OU PIPS. This version returns pandas dataframes"""

    pdatetimes = []
    intensities = []
    preciptots = []
    reflectivities = []
    sampleintervals = []
    pcounts = []
    pcounts2 = []
    sensortemps = []
    amplitudes = []
    voltages = []
    pvoltages = []

    # Will contain the number of drops in the diameter vs. fall speed matrix for each time
    # i.e. a (numtimes,32,32) array
    countsMatrix = []
    # Will contain computed number concentrations for each size bin for each time
    # i.e. a (numtimes,32) array
    concentrations = []

    onedrop_concentrations = []

    # dates = []
    # times = []
    datetimes = []
    windspds = []
    winddirrels = []
    winddirabss = []
    winddiags = []
    fasttemps = []
    slowtemps = []
    dewpoints = []
    RHs_derived = []
    RHs = []
    pressures = []
    compass_dirs = []
    GPS_lats = []
    GPS_lons = []
    GPS_stats = []
    GPS_alts = []

    disfile = open(filename, 'r')

    firstgoodGPS = False

    for line in disfile:
        tokens = line.strip().split(',')
        timestamp = tokens[0]
        timestring = timestamp.strip().split()
        date = timestring[0]  # .strip('-')
        time = timestring[1]  # .strip(':')

        # 2016-03-31 22:19:02

        # Construct datetime object
        year = np.int(date[:4])
        month = np.int(date[5:7])
        day = np.int(date[8:])
        hour = np.int(time[:2])
        min = np.int(time[3:5])
        sec = np.int(time[6:8])

        datetimelogger = datetime(year, month, day, hour, min, sec)

        # recordnum = np.int(tokens[1])
        voltage = np.float(tokens[2])
        # paneltemp = np.float(tokens[3])
        winddirrel = np.float(tokens[4])
        windspd = np.float(tokens[5])
        winddiag = np.float(tokens[6])
        fasttemp = np.float(tokens[7])
        slowtemp = np.float(tokens[8])
        RH = np.float(tokens[9])
        pressure = np.float(tokens[10])
        compass_dir = np.float(tokens[11])
        GPS_time = tokens[12]
        GPS_status = tokens[13]
        GPS_lat = np.float(tokens[14])
        GPS_lat_hem = tokens[15]
        GPS_lat = DDMtoDD(GPS_lat, GPS_lat_hem)
        GPS_lon = np.float(tokens[16])
        GPS_lon_hem = tokens[17]
        GPS_lon = DDMtoDD(GPS_lon, GPS_lon_hem)
        # GPS_spd = np.float(tokens[18])
        # GPS_dir = np.float(tokens[19])
        GPS_date = tokens[20]
        # GPS_magvar = np.float(tokens[21])
        GPS_alt = np.float(tokens[22])

        winddirabs = np.float(tokens[23])
        try:
            dewpoint = np.float(tokens[24])
        except BaseException:
            dewpoint = np.nan
        try:
            RH_derived = np.float(tokens[25])
        except BaseException:
            RH_derived = np.nan

        # Find the first good GPS time and date and use that
        # to construct the time offset for the data logger
        if(not np.isnan(GPS_alt) and not firstgoodGPS and GPS_status == 'A'):
            firstgoodGPS = True
            print(GPS_date, GPS_time)
            print(date, time)

            print(year, month, day, hour, min, sec)

            # Construct datetime object
            gyear = np.int('20' + GPS_date[4:])
            gmonth = np.int(GPS_date[2:4])
            gday = np.int(GPS_date[:2])
            ghour = np.int(GPS_time[:2])
            gmin = np.int(GPS_time[2:4])
            gsec = np.int(GPS_time[4:])

            datetimeGPS = datetime(gyear, gmonth, gday, ghour, gmin, gsec)
            GPS_offset = datetimeGPS - datetimelogger
            print(datetimeGPS, datetimelogger)
            print("GPS Offset", GPS_offset)

        datetimes.append(datetimelogger)
        windspds.append(windspd)
        winddirrels.append(winddirrel)
        winddirabss.append(winddirabs)
        winddiags.append(winddiag)
        fasttemps.append(fasttemp)
        slowtemps.append(slowtemp)
        dewpoints.append(dewpoint)
        RHs_derived.append(RH_derived)
        RHs.append(RH)
        pressures.append(pressure)
        compass_dirs.append(compass_dir)
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)
        voltages.append(voltage)

        parsivel_string = tokens[26]
        parsivel_tokens = parsivel_string.strip().split(';')
        serialnum = parsivel_tokens[0]
        if(serialnum in parsivel_ids and len(parsivel_tokens) >= 11):
            # print timestring
            precipintensity = np.float(parsivel_tokens[1])
            precipaccum = np.float(parsivel_tokens[2])
            parsivel_dBZ = np.float(parsivel_tokens[3])
            sample_interval = np.float(parsivel_tokens[4])
            signal_amplitude = np.float(parsivel_tokens[5])
            pcount = np.int(parsivel_tokens[6])
            sensor_temp = np.float(parsivel_tokens[7])
            pvoltage = np.float(parsivel_tokens[8])
            # sensor_time = parsivel_tokens[9]
            # sensor_date = parsivel_tokens[10]
            try:
                spectrum = [float(x) if x != '' else 0 for x in parsivel_tokens[11:]]
            except BaseException:
                spectrum = [-999 for i in range(1025)]
            # print "spectrum length = ",len(spectrum)
            if(len(spectrum) < 1024):
                print("Problem with Parsivel spectrum.  Flagging as bad!")
                print("Time: ", timestring)
#                 spectrum = [np.nan for i in xrange(1025)]
            else:
                if(len(spectrum) == 1025):
                    spectrum = spectrum[:-1]  # Strip off bogus last value

#             if(GPS_time):
#                 print GPS_date,GPS_time
#
            pdatetimes.append(datetimelogger)
            preciptots.append(precipaccum)
            intensities.append(precipintensity)
            reflectivities.append(parsivel_dBZ)
            sampleintervals.append(sample_interval)
            amplitudes.append(signal_amplitude)
            pcounts.append(pcount)
            sensortemps.append(sensor_temp)
            pvoltages.append(pvoltage)

            # Now create an array out of the spectrum and reshape it to 32x32
            spectrum = np.array(spectrum, dtype='int')
            if(spectrum.size == 1024):
                spectrum = spectrum.reshape((32, 32))
            else:
                spectrum = -999 * np.ones((32, 32), dtype='int')
                print(spectrum.size)
            # Append spectrum (corrected or not) to spectrum time list

            countsMatrix.append(spectrum)

    # Recast countsMatrix as numpy array

    countsMatrix = np.dstack(countsMatrix)
    countsMatrix = np.rollaxis(countsMatrix, 2, 0)

    # Perform Katja's QC routines if desired (should not be used in combination with my methods
    # above, most are redundant anyway).

    X, Y = np.meshgrid(avg_diameter, fall_bins)
    flaggedtimes = np.zeros(len(pdatetimes), dtype=bool)
    # splashmask = np.zeros_like(countsMatrix)
    # marginmask = np.zeros_like(countsMatrix)

    countsMatrix = truncatedspectrumQC(countsMatrix)

    if(use_strongwindQC or basicqc or strongwindqc):
        countsMatrix, flaggedtimes = strongwindQC(countsMatrix)

    if(use_splashingQC or basicqc):
        countsMatrix = splashingQC(countsMatrix)

    if(use_marginQC or basicqc):
        countsMatrix = marginQC(countsMatrix)

    if(use_rainfallspeedQC or rainfallqc):
        countsMatrix = rainfallspeedQC(countsMatrix, rainvd, falltol, maskhigh, masklow)

    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)

    if(use_rainonlyQC or rainonlyqc):
        countsMatrix = rainonlyQC(countsMatrix)

    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)

    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)
    # Find total number of non-masked particles

    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)

    counts_1drop = np.ones_like(avg_diameter)

    # print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t, time in enumerate(pdatetimes):
        spectrum = countsMatrix[t, :]
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy, vspectrum = np.meshgrid(avg_diameter, fall_bins)
            dspectrum = spectrum
        else:
            vspectrum = rainvd
            dspectrum = spectrum.sum(axis=0)   # Sum up particles for each diameter bin

            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]

        # Now compute the number concentration using the assumed fall speeds, sensor area, and
        # sampling interval
        # Units are #/m^3/mm
#         if time == '171320':
#             print spectrum.size,spectrum,countsMatrix.data[t,:],spectrum.data
        if(spectrum.size == 1024 and not flaggedtimes[t]):
            concentration = dspectrum / \
                (vspectrum * sampleintervals[t] * eff_sensor_area * (max_diameter - min_diameter))
            onedrop_concentration = counts_1drop / \
                (rainvd * sampleintervals[t] * eff_sensor_area * (max_diameter - min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
            # print "concentration.shape"
            # print concentration.shape
        elif(not flaggedtimes[t]):
            concentration = np.zeros_like(avg_diameter)
        else:
            concentration = -999. * np.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999., concentration)

        # Throw out particles above and below a certain diameter if desired

#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)

        # print "Number of particles counted vs. summed number: "
        # print np.int(line[7]),np.sum(dspectrum)

        # print concentration

        concentrations.append(concentration)
        onedrop_concentrations.append(onedrop_concentration)

    concentrations = ma.array(concentrations)
    onedrop_concentrations = ma.array(onedrop_concentrations)
    # print "concentrations: ",concentrations
    pcounts = np.array(pcounts)
    # pcounts2 = ma.array(pcounts2)

    # Correct the logger time and date using the GPS time

    datetimes_corrected = []
    pdatetimes_corrected = []
    for datetimelogger in datetimes:
        datetimes_corrected.append(datetimelogger + GPS_offset)
    for pdatetimelogger in pdatetimes:
        pdatetimes_corrected.append(pdatetimelogger + GPS_offset)

    # Create a DataFrame for the conventional data

    conv_data_dict = {'wind speed': windspds, 'relative wind direction': winddirrels,
                      'absolute wind direction': winddirabss, 'wind diagnostic': winddiags,
                      'fast temperature': fasttemps, 'slow temperature': slowtemps,
                      'dewpoint': dewpoints, 'derived RH': RHs_derived, 'measured RH': RHs,
                      'pressure': pressures, 'compass direction': compass_dirs, 'GPS lat': GPS_lats,
                      'GPS lon': GPS_lons, 'GPS status': GPS_stats, 'GPS altitude': GPS_alts,
                      'voltage': voltages}

    conv_data_df = pd.DataFrame(conv_data_dict, index=pd.Index(datetimes_corrected, name='Time'))

    # Average and thin the DSD data with the desired interval
    DSD_index_interval = int(DSD_interval / 10.0)
    print(
        "Requested DSD interval: {:.1f}. Actual DSD interval: {:.1f}".format(
            DSD_interval,
            DSD_index_interval *
            10.0))
    DSD_interval = DSD_index_interval * 10.0
    # We need to find the offset corresponding to the starting second and then
    # generate the frequency string
    sec_offset = pdatetimes_corrected[0].second
    intervalstr = '{:d}S'.format(int(DSD_interval))

    # Create dataframes for the Parsivel data
    # One dataframe for the 3D concentration array, one for the 3D "one-drop" concentration array,
    # and one to contain all the other derived data

    derived_DSD_dict = {'Precipitation intensity': intensities, 'Precipitation totals': preciptots,
                        'Radar reflectivity': reflectivities, 'Particle counts': pcounts,
                        'Derived particle counts': pcounts2, 'Bad wind flagged times': flaggedtimes}

    concentrations_df = pd.DataFrame(
        data=concentrations,
        index=pd.Index(
            pdatetimes_corrected,
            name='Time'),
        columns=avg_diameter)
    onedrop_concentrations_df = pd.DataFrame(
        data=onedrop_concentrations, index=pd.Index(
            pdatetimes_corrected, name='Time'), columns=avg_diameter)
    derived_DSD_df = pd.DataFrame(
        derived_DSD_dict, index=pd.Index(
            pdatetimes_corrected, name='Time'))

    if(DSD_interval > 10.0):
        # Resample at the new interval filling in missing values with zeros

        concentrations_df = concentrations_df.resample(
            intervalstr, label='right', closed='right', base=sec_offset).mean().fillna(0)
        onedrop_concentrations_df = onedrop_concentrations_df.resample(
            intervalstr, label='right', closed='right', base=sec_offset).mean().fillna(0)
        derived_DSD_df = derived_DSD_df.resample(
            intervalstr,
            label='right',
            closed='right',
            base=sec_offset).mean().fillna(0)

    return (conv_data_df, concentrations_df, onedrop_concentrations_df, derived_DSD_df,
            DSD_interval, intervalstr)


def readPIPSloc(filename, starttime=None, stoptime=None):
    """Reads the location of the PIPS from the data file.  Averages all valid GPS lat/lons"""

    if starttime is not None:
        starttimestamp = datetime(np.int(starttime[:4]), np.int(starttime[4:6]), np.int(
            starttime[6:8]), np.int(starttime[8:10]), np.int(starttime[10:12]))
    else:
        starttimestamp = None
    if stoptime is not None:
        stoptimestamp = datetime(np.int(stoptime[:4]), np.int(stoptime[4:6]), np.int(
            stoptime[6:8]), np.int(stoptime[8:10]), np.int(stoptime[10:12]))
    else:
        stoptimestamp = None

    GPS_lats = []
    GPS_lons = []
    GPS_stats = []
    GPS_alts = []

    disfile = open(filename, 'r')

    for line in disfile:
        tokens = line.strip().split(',')
        # Check for header line (older versions don't have it)
        if(tokens[0] == 'TIMESTAMP'):
            continue

        timestamp = tokens[0]
        timestring = timestamp.strip().split()
        date = timestring[0]  # .strip('-')
        time = timestring[1]  # .strip(':')

        # Construct datetime object
        year = np.int(date[:4])
        month = np.int(date[5:7])
        day = np.int(date[8:])
        hour = np.int(time[:2])
        min = np.int(time[3:5])
        sec = np.int(time[6:8])

        datetimelogger = datetime(year, month, day, hour, min, sec)
        # Skip this record if it lies before or after the desired period
        if starttimestamp is not None:
            if datetimelogger < starttimestamp:
                continue
        if stoptimestamp is not None:
            if datetimelogger > stoptimestamp:
                continue

        GPS_status = tokens[13]
        GPS_lat = np.float(tokens[14])
        GPS_lat_hem = tokens[15]
        GPS_lat = DDMtoDD(GPS_lat, GPS_lat_hem)
        GPS_lon = np.float(tokens[16])
        GPS_lon_hem = tokens[17]
        GPS_lon = DDMtoDD(GPS_lon, GPS_lon_hem)
        GPS_alt = np.float(tokens[22])

        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)

    # Find the disdrometer location by averaging the valid GPS lats,lons, and alts
    GPSnonvalid = (np.array(GPS_stats) != 'A')
    GPS_lats_masked = ma.masked_where(GPSnonvalid, np.array(GPS_lats))
    GPS_lons_masked = ma.masked_where(GPSnonvalid, np.array(GPS_lons))
    GPS_alts_masked = ma.masked_where(GPSnonvalid, np.array(GPS_alts))

    GPS_lats_masked = ma.masked_invalid(GPS_lats_masked)
    GPS_lons_masked = ma.masked_invalid(GPS_lons_masked)
    GPS_alts_masked = ma.masked_invalid(GPS_alts_masked)

    lat = GPS_lats_masked.mean()
    lon = GPS_lons_masked.mean()
    alt = GPS_alts_masked.mean()

    # There's a bug in numpy.ma that sometimes causes operations such as mean() to return
    # a 0d array instead of a scalar.  To remedy this, explicitly cast them as scalars
    # here.

    lat = lat.item()
    lon = lon.item()
    alt = alt.item()

    dloc = (lat, lon, alt)

    return GPS_lats, GPS_lons, GPS_stats, GPS_alts, dloc


def readCUloc(filename, starttime=None, stoptime=None):
    """Reads the location of a CU probe from the data file. Averages all valid GPS lats/lons"""

    if starttime is not None:
        starttimestamp = datetime(np.int(starttime[:4]), np.int(starttime[4:6]), np.int(
            starttime[6:8]), np.int(starttime[8:10]), np.int(starttime[10:12]))
    else:
        starttimestamp = None
    if stoptime is not None:
        stoptimestamp = datetime(np.int(stoptime[:4]), np.int(stoptime[4:6]), np.int(
            stoptime[6:8]), np.int(stoptime[8:10]), np.int(stoptime[10:12]))
    else:
        stoptimestamp = None

    GPS_lats = []
    GPS_lons = []
    GPS_alts = []
    GPS_stats = []

    convfile = open(filename, 'r')

    # Read in the header information

    convfile.readline()
    convfile.readline()
    convfile.readline()
    convfile.readline()

    for line in convfile:
        mysplit = shlex.shlex(line)
        mysplit.whitespace += ','
        mysplit.whitespace_split = True
        # line = l.strip().split(',')
        tokens = list(mysplit)
        datetimestring = tokens[0]
        # Construct datetime object
        year = int(datetimestring[1:5])
        month = int(datetimestring[6:8])
        day = int(datetimestring[9:11])
        hour = int(datetimestring[12:14])
        min = int(datetimestring[15:17])
        sec = int(datetimestring[18:20])

        datetimelogger = datetime(year, month, day, hour, min, sec)
        # Skip this record if it lies before or after the desired period
        if starttimestamp is not None:
            if datetimelogger < starttimestamp:
                continue
        if stoptimestamp is not None:
            if datetimelogger > stoptimestamp:
                continue
        try:
            GPS_status = np.int(tokens[15].strip('\"'))
        except Exception:
            GPS_status = 0
        try:
            GPS_alt = np.float(tokens[16].strip('\"'))
        except Exception:
            GPS_alt = np.nan
        try:
            GPS_lat = np.float(tokens[11].strip('\"'))
            GPS_lat = DDMtoDD(GPS_lat / 100., 'N')
        except Exception:
            GPS_lat = np.nan
            GPS_status = 0
        try:
            GPS_lon = np.float(tokens[12].strip('\"'))
            GPS_lon = DDMtoDD(GPS_lon / 100., 'W')
        except Exception:
            GPS_lon = np.nan
            GPS_status = 0

#         print datetimelogger, GPS_lat, GPS_lon, GPS_alt, GPS_status
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)

    # Find the disdrometer location by averaging the valid GPS lats,lons, and alts
    GPSnonvalid = (np.array(GPS_stats) == 0)
    GPS_lats_masked = ma.masked_where(GPSnonvalid, np.array(GPS_lats))
    GPS_lons_masked = ma.masked_where(GPSnonvalid, np.array(GPS_lons))
    GPS_alts_masked = ma.masked_where(GPSnonvalid, np.array(GPS_alts))

    GPS_lats_masked = ma.masked_invalid(GPS_lats_masked)
    GPS_lons_masked = ma.masked_invalid(GPS_lons_masked)
    GPS_alts_masked = ma.masked_invalid(GPS_alts_masked)

    lat = GPS_lats_masked.mean()
    lon = GPS_lons_masked.mean()
    alt = GPS_alts_masked.mean()

    # There's a bug in numpy.ma that sometimes causes operations such as mean() to return
    # a 0d array instead of a scalar.  To remedy this, explicitly cast them as scalars
    # here.

    lat = lat.item()
    lon = lon.item()
    alt = alt.item()

    dloc = (lat, lon, alt)

    return GPS_lats, GPS_lons, GPS_stats, GPS_alts, dloc


def readNV2loc(filename):
    """Reads the location information from one of the V2 NSSL probe netCDF data files."""
    disfile = netcdf.Dataset(filename, 'r')
    # Read in latitude and longitude of the probe (too bad altitude isn't here, too)
    lat = disfile.variables['latitude'][...]
    print("lat = ", lat)
    lon = disfile.variables['longitude'][...]
    disfile.close()

    return (lat, lon, None)


def readPIPStimerange(filename):
    """Reads in a PIPS data file and returns the range of times within the file"""

    # Below is adapted from https://stackoverflow.com/questions/3346430/
    # what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file
    with open(filename, "rb") as f:
        first = f.readline()        # Read the first line.
        second = f.readline()       # Read second line in case we need it (first may be header)
        f.seek(-2, os.SEEK_END)     # Jump to the second last byte.
        while f.read(1) != b"\n":   # Until EOL is found...
            f.seek(-2, os.SEEK_CUR)  # ...jump back the read byte plus one more.
        last = f.readline()         # Read last line.

    # First time
    tokens = first.strip().split(b',')
    # Check for header line (older versions don't have it)
    if(tokens[0] == b'TIMESTAMP'):
        tokens = second.strip().split(b',')

    timestamp = tokens[0]
    timestring = timestamp.strip().split()
    firstdate = timestring[0]  # .strip('-')
    firsttime = timestring[1]  # .strip(':')

    # 2016-03-31 22:19:02

    # Construct datetime object
    year = np.int(firstdate[:4])
    month = np.int(firstdate[5:7])
    day = np.int(firstdate[8:])
    hour = np.int(firsttime[:2])
    min = np.int(firsttime[3:5])
    sec = np.int(firsttime[6:8])

    datetimefirst = datetime(year, month, day, hour, min, sec)

    # Last time
    tokens = last.strip().split(b',')
    # Check for header line (older versions don't have it)

    timestamp = tokens[0]
    timestring = timestamp.strip().split()
    lastdate = timestring[0]  # .strip('-')
    lasttime = timestring[1]  # .strip(':')

    # 2016-03-31 22:19:02

    # Construct datetime object
    year = np.int(lastdate[:4])
    month = np.int(lastdate[5:7])
    day = np.int(lastdate[8:])
    hour = np.int(lasttime[:2])
    min = np.int(lasttime[3:5])
    sec = np.int(lasttime[6:8])

    datetimelast = datetime(year, month, day, hour, min, sec)

    return (firstdate, firsttime, datetimefirst), (lastdate, lasttime, datetimelast)


def readPIPSstation(filename, fixGPS=True, starttimestamp=None, stoptimestamp=None):
    """Reads in just the data from a PIPS file that is needed for a station plot"""

    if starttimestamp is not None:
        starttime = datetime.strptime(starttimestamp, '%Y%m%d%H%M%S')
    else:
        starttime = None
    if stoptimestamp is not None:
        stoptime = datetime.strptime(stoptimestamp, '%Y%m%d%H%M%S')
    else:
        stoptime = None

    # dates = []
    # times = []
    datetimes = []
    windspds = []
    winddirabss = []
    winddiags = []
    fasttemps = []
    dewpoints = []
    pressures = []
    GPS_lats = []
    GPS_lons = []
    GPS_stats = []
    GPS_alts = []

    disfile = open(filename, 'r')

    firstgoodGPS = False

    for line in disfile:
        tokens = line.strip().split(',')
        timestamp = tokens[0]
        timestring = timestamp.strip().split()
        date = timestring[0]  # .strip('-')
        time = timestring[1]  # .strip(':')

        # 2016-03-31 22:19:02

        # Construct datetime object
        year = np.int(date[:4])
        month = np.int(date[5:7])
        day = np.int(date[8:])
        hour = np.int(time[:2])
        min = np.int(time[3:5])
        sec = np.int(time[6:8])

        datetimelogger = datetime(year, month, day, hour, min, sec)

        # Skip this record if it lies before or after the desired period
        if starttimestamp is not None:
            if datetimelogger < starttime:
                continue
        if stoptimestamp is not None:
            if datetimelogger > stoptime:
                continue

        windspd = np.float(tokens[5])
        winddiag = np.float(tokens[6])
        fasttemp = np.float(tokens[7])
        pressure = np.float(tokens[10])
        GPS_time = tokens[12]
        GPS_status = tokens[13]
        GPS_lat = np.float(tokens[14])
        GPS_lat_hem = tokens[15]
        GPS_lat = DDMtoDD(GPS_lat, GPS_lat_hem)
        GPS_lon = np.float(tokens[16])
        GPS_lon_hem = tokens[17]
        GPS_lon = DDMtoDD(GPS_lon, GPS_lon_hem)
        # GPS_spd = np.float(tokens[18])
        # GPS_dir = np.float(tokens[19])
        GPS_date = tokens[20]
        GPS_alt = np.float(tokens[22])

        winddirabs = np.float(tokens[23])
        try:
            dewpoint = np.float(tokens[24])
        except BaseException:
            dewpoint = np.nan

        # Find the first good GPS time and date and use that
        # to construct the time offset for the data logger
        if(not np.isnan(GPS_alt) and not firstgoodGPS and GPS_status == 'A'):
            firstgoodGPS = True
            print(GPS_date, GPS_time)
            print(date, time)

            print(year, month, day, hour, min, sec)

            # Construct datetime object
            gyear = np.int('20' + GPS_date[4:])
            gmonth = np.int(GPS_date[2:4])
            gday = np.int(GPS_date[:2])
            ghour = np.int(GPS_time[:2])
            gmin = np.int(GPS_time[2:4])
            gsec = np.int(GPS_time[4:])

            datetimeGPS = datetime(gyear, gmonth, gday, ghour, gmin, gsec)
            GPS_offset = datetimeGPS - datetimelogger
            print(datetimeGPS, datetimelogger)
            print("GPS Offset", GPS_offset)

        datetimes.append(datetimelogger)
        windspds.append(windspd)
        winddirabss.append(winddirabs)
        winddiags.append(winddiag)
        fasttemps.append(fasttemp)
        dewpoints.append(dewpoint)
        pressures.append(pressure)
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)

    # Correct the logger time and date using the GPS time

    # Sometimes the GPS has major issues for the whole deployment.
    # In this case, set the offset to 0
    if(not firstgoodGPS):
        GPS_offset = timedelta(seconds=0)

    datetimes_corrected = []
    for datetimelogger in datetimes:
        datetimes_corrected.append(datetimelogger + GPS_offset)

    return datetimes_corrected, windspds, winddirabss, fasttemps, dewpoints, pressures,    \
        GPS_lats, GPS_lons, GPS_stats, GPS_alts


def PIPS2TTU(filename, fixGPS=True):
    """Reads conventional data from a PIPS text file and outputs it in the format used by TTU's
       StickNets"""

    (datetimes_corrected, pdatetimes_corrected, flaggedtimes, intensities, preciptots,
     reflectivities, pcounts, pcounts2, sensortemps, concentrations, gamma_concentrations,
     countsMatrix, windspds, winddirrels, winddirabss, winddiags, fasttemps, slowtemps, dewpoints,
     RHs_derived, RHs, pressures, compass_dirs, GPS_lats, GPS_lons, GPS_stats, GPS_alts,
     voltages) = readPIPS(filename)

    GPS_lats, GPS_lons, GPS_stats, GPS_alts, dloc = readPIPSloc(filename)

    splitpath = os.path.split(filename)
    inputfilename = splitpath[1]
    # outputfiledir = splitpath[0]
    outputfilename = inputfilename[:-4] + '_TTU.txt'
    print(outputfilename)
    # outputfilepath = outputfiledir + outputfilename

    # outputfile = open(outputfilepath, 'w')

    probeID = os.path.split(filename)[1][5:7]
    date = str(datetimes_corrected[0].year) + \
        str(datetimes_corrected[0].month) + str(datetimes_corrected[0].day)
    time = str(datetimes_corrected[0].hour) + ':' + \
        str(datetimes_corrected[0].minute) + ':' + str(datetimes_corrected[0].second)
    lat = dloc[0]
    lon = dloc[1]
    alt = dloc[2]
    heading = compass_dirs[0]
    dummy = 3

    header = [probeID, date, time, lat, lon, alt, heading, dummy]
    header = ','.join(map(str, header))
    print(header)

    return date


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
    # Added by Jessie (was temp > 0).  Find out why...
    rhv = np.where(Zh != Zv, Zhv / (np.sqrt(temp)), 0.0)
    # np.savetxt('temp.txt', temp)

    dualpol_dict = {'ZH': Zh, 'ZV': Zv, 'ZHV': Zhv, 'dBZ': dBZ, 'ZDR': ZDR, 'KDP': Kdp, 'RHV': rhv,
                    'intv': intv, 'd': d, 'fa2': fa2, 'fb2': fb2}

    return dualpol_dict


def calc_D0_bin(D, Dl, Dr, ND, bin_width):
    """Calculate D0 for a binned distribution"""
    temp_M3 = (D**3.) * ND * bin_width
    temp_M3_cumsum = np.cumsum(temp_M3)  # Cumulative sum of M3 with increasing bin size
    temp_M3_sum = np.sum(temp_M3)
    pro = temp_M3 / np.sum(temp_M3_sum)  # Proportion of M3 in each bin
    pro_cumsum = temp_M3_cumsum / np.sum(temp_M3_sum)  # Cumulative proportion of M3 in each bin
    # Compute median volume diameter using a linear interpolation within the "mass-midpoint" bin.
    # Source: http://www.dropletmeasurement.com/PADS_Help/MVD_(um).htm
    # (originally from FAA Electronic Aircraft Icing Handbook)
    try:
        # Find index of bin where cumulative sum exceeds 1/2 of total
        medindex = np.where(pro_cumsum > 0.5)[0][0]
        b1 = Dl[medindex]  # Lower boundary of mass-midpoint bin
        b2 = Dr[medindex]  # Upper boundary of mass-midpoint bin
        if(medindex == 0):
            pro_cumsum_medindexm1 = 0.0
        else:
            pro_cumsum_medindexm1 = pro_cumsum[medindex - 1]
        # Linearly-interpolated (within mass-midpoint bin) D0
        D0 = b1 + ((0.5 - pro_cumsum_medindexm1) / pro[medindex]) * (b2 - b1)
    except Exception:
        medindex = 0
        D0 = np.nan
    return D0


def calc_DSD(Nc_bin, rho, qrQC=False, qr_thresh=None,
             pcounts=None, intensities=None, pc=None):
    """
    Fits exponential and gamma DSDs to disdrometer data and returns several DSD related quantities.
    """

    # First calculate the required moment estimators
    M0 = []   # Not used in moment estimator, but used for mean diameter calculation
    M1 = []   # Ditto
    M2 = []
    M3 = []
    M4 = []
    M6 = []
    M7 = []   # Ditto
    D_med_disd = []    # Median volume diameter
    Dmax = []  # largest observed drop diameter bin
    Dmin = []

    try:
        numtimes = np.size(Nc_bin, axis=1)
    except BaseException:
        numtimes = 1
        Nc_bin = Nc_bin[:, np.newaxis]

    # FIXME: Change this to use functioned defined above with air density correction!
#     v = -0.1021 + 4.932 * avg_size - 0.9551 * avg_size**2. + \
#         0.07934 * avg_size**3. - 0.002362 * avg_size**4.

    v = assignfallspeed(avg_diameter, rhocorrect=True, rho=rho)
    v = v.T
    rainrate = []

    for t in range(numtimes):

        temp_M0 = (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M1 = ((avg_diameter[:] / 1000.)) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M2 = ((avg_diameter[:] / 1000.)**2.) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M3 = ((avg_diameter[:] / 1000.)**3.) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M4 = ((avg_diameter[:] / 1000.)**4.) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M6 = ((avg_diameter[:] / 1000.)**6.) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.
        temp_M7 = ((avg_diameter[:] / 1000.)**7.) * (1000. * Nc_bin[:, t]) * bin_width[:] / 1000.

        # Note, this is *observed* rainrate calculated directly from the disdrometer bins
        temp_rainrate = np.sum((6. * 10.**-4.) * np.pi * v *
                               avg_diameter[:]**3. * Nc_bin[:, t] * bin_width[:])
        rainrate.append(temp_rainrate)

        # Before summing up the 3rd moment, use the bin information to find the
        # (approximate) median volume diameter

        # print "temp_M3",temp_M3,temp_M3.shape
        # print "np.array(temp_M3)",np.array(temp_M3)
        temp_M3_cumsum = np.cumsum(temp_M3)  # Cumulative sum of M3 with increasing bin size
        temp_M3_sum = np.sum(temp_M3)
        pro = temp_M3 / np.sum(temp_M3_sum)  # Proportion of M3 in each bin
        # Cumulative proportion of M3 in each bin
        pro_cumsum = temp_M3_cumsum / np.sum(temp_M3_sum)

        # print "temp_M3_cumsum",temp_M3_cumsum

        # Compute median volume diameter using a linear interpolation within the "mass-midpoint"
        # bin.
        # Source: http://www.dropletmeasurement.com/PADS_Help/MVD_(um).htm
        # (originally from FAA Electronic Aircraft Icing Handbook)
        try:
            # Find index of bin where cumulative sum exceeds 1/2 of total
            medindex = np.where(pro_cumsum > 0.5)[0][0]
            b1 = min_diameter[medindex]  # Lower boundary of mass-midpoint bin
            b2 = max_diameter[medindex]  # Upper boundary of mass-midpoint bin
            if(medindex == 0):
                pro_cumsum_medindexm1 = 0.0
            else:
                pro_cumsum_medindexm1 = pro_cumsum[medindex - 1]
            # Linearly-interpolated (within mass-midpoint bin) D0
            temp_D_med = b1 + ((0.5 - pro_cumsum_medindexm1) / pro[medindex]) * (b2 - b1)
        except BaseException:
            medindex = 0
            temp_D_med = np.nan

        temp_Dmax = 0.  # initialize in case there are no drops ?
        temp_Dmin = 0.

        for index, value in reversed(list(enumerate(Nc_bin[:, t]))):
            if (value > 0.):
                temp_Dmax = avg_diameter[index]
                break
        for index2, value2 in list(enumerate(Nc_bin[:, t])):
            if (value2 > 0.):
                temp_Dmin = avg_diameter[index2]
                break

        # temp_M0_2 = ma.sum(temp_M0)
        temp_M0 = np.sum(temp_M0)
        temp_M1 = np.sum(temp_M1)
        temp_M2 = np.sum(temp_M2)
        temp_M4 = np.sum(temp_M4)
        temp_M3 = np.sum(temp_M3)
        temp_M6 = np.sum(temp_M6)
        temp_M7 = np.sum(temp_M7)

        M0.append(temp_M0)
        M1.append(temp_M1)
        M2.append(temp_M2)
        M3.append(temp_M3)
        M4.append(temp_M4)
        M6.append(temp_M6)
        M7.append(temp_M7)
        D_med_disd.append(temp_D_med)
        Dmax.append(temp_Dmax)
        Dmin.append(temp_Dmin)

    M0 = ma.array(M0, dtype=np.float64)
    M1 = ma.array(M1, dtype=np.float64)
    M2 = ma.array(M2, dtype=np.float64)
    M3 = ma.array(M3, dtype=np.float64)
    M4 = ma.array(M4, dtype=np.float64)
    M6 = ma.array(M6, dtype=np.float64)
    M7 = ma.array(M7, dtype=np.float64)
    D_med_disd = np.array(D_med_disd)
    rainrate = np.array(rainrate)
    Dmax = np.array(Dmax) / 1000.
    Dmin = np.array(Dmin) / 1000.

    # --- Compute various mean diameters directly from measured discrete distribution ---
    cmr = (np.pi / 6.) * \
        1000.                                     # Constant in mass-diameter relation for rain
    LWC_disd = cmr * 1000.0 * M3
    # print LWC_disd.shape
    # print rho.shape
    QR_disd = LWC_disd / rho
    if(qrQC):
        qrmask1D = np.where(QR_disd > qr_thresh, True, False)
        qr2D = QR_disd.reshape(1, numtimes).repeat(32, 0)
        qrmask2D = np.where(qr2D > qr_thresh, True, False)

        # print qrmask2D.shape

        # Mask out all the needed arrays, then all derived arrays below will be masked in the same
        # way

        D_med_disd = ma.masked_array(D_med_disd, mask=qrmask1D)
        Nc_bin = ma.masked_array(Nc_bin, mask=qrmask2D)
        # logNc_bin = ma.masked_array(logNc_bin, mask=qrmask2D)
        M0 = ma.masked_array(M0, mask=qrmask1D)
        M1 = ma.masked_array(M1, mask=qrmask1D)
        M2 = ma.masked_array(M2, mask=qrmask1D)
        M3 = ma.masked_array(M3, mask=qrmask1D)
        M4 = ma.masked_array(M4, mask=qrmask1D)
        M6 = ma.masked_array(M6, mask=qrmask1D)
        M7 = ma.masked_array(M7, mask=qrmask1D)
        rho = ma.masked_array(rho, mask=qrmask1D)
        QR_disd = ma.masked_array(QR_disd, mask=qrmask1D)
        LWC_disd = ma.masked_array(LWC_disd, mask=qrmask1D)
#
    # TODO: Make the masking by particle counts adjustable in the pyPIPScontrol file
    # There may be times when we don't want such stringent masking, and it is dependent on the
    # integration interval we choose, anyway.
    try:
        masklowcounts = pc.masklowcounts
    except BaseException:
        masklowcounts = False

    if(masklowcounts):
        nummask1D = np.where((pcounts < 50.) | (rainrate < 1.), True, False)
        num2D = pcounts.reshape(1, numtimes).repeat(32, 0)
        rain2D = rainrate.reshape(1, numtimes).repeat(32, 0)
        nummask2D = np.where((num2D < 50.) | (rain2D < 1.), True, False)

        D_med_disd = ma.masked_array(D_med_disd, mask=nummask1D)
        Nc_bin = ma.masked_array(Nc_bin, mask=nummask2D)
        # logNc_bin = ma.masked_array(logNc_bin, mask=nummask2D)
        M0 = ma.masked_array(M0, mask=nummask1D)
        M1 = ma.masked_array(M1, mask=nummask1D)
        M2 = ma.masked_array(M2, mask=nummask1D)
        M3 = ma.masked_array(M3, mask=nummask1D)
        M4 = ma.masked_array(M4, mask=nummask1D)
        M6 = ma.masked_array(M6, mask=nummask1D)
        M7 = ma.masked_array(M7, mask=nummask1D)
        rho = ma.masked_array(rho, mask=nummask1D)
        QR_disd = ma.masked_array(QR_disd, mask=nummask1D)
        LWC_disd = ma.masked_array(LWC_disd, mask=nummask1D)
        rainrate = ma.masked_array(rainrate, mask=nummask1D)

    # Compute reflectivity from the 6th moment
    refl_disd = 10.0 * np.log10(1e18 * M6)
    # print 'reflectivity (disdrometer) = ',refl_disd

    # Compute the mass-weighted mean diameter (M(4)/M(3) )

    D_m_disd = (M4 / M3) * 1000.0

    # Compute the mean-volume (or mean-mass if rho=constant) diameter ((M(3)/M(0))^1./3.)

    D_mv_disd = ((M3 / M0)**(1. / 3.)) * 1000.0

    # Compute reflectivity-weighted mean diameter (M(7)/M(6))

    D_ref_disd = (M7 / M6) * 1000.0

    # --- Exponential Distribution Fit by Method of Moments (Zhang et al. 2008) ---

    # Now that we have the moment estimates, we can calculate N0 and lamda
    # (eqns. 3 and 4 in Zhang et al. 2008)

    # gamma3 = special.gamma(3.)
    gamma4 = special.gamma(4.)
    # gamma5 = special.gamma(5.)
    gamma7 = special.gamma(7.)

    # Uncomment if you want estimates based on M2 and M4
    # lamda = np.where(M4 == 0.0, 0.0, ((M2*gamma5)/(M4*gamma3))**(1./2.))
    # N0 = (M2*lamda**3.)/gamma3

    # Uncomment if you want estimates based on M3 and M6
    lamda_exp = np.where(M6 == 0.0, 0.0, ((M3 * gamma7) / (M6 * gamma4))
                         ** (1. / 3.))  # Should be masked
    N0_exp = (M3 * lamda_exp**4.) / gamma4                                      # Should be masked
    mu_exp = 0.0

    # print 'lamda (exp) = ',lamda_exp
    # print 'N0 (exp) = ',N0_exp

    # Now create synthetic bins to plot the derived exponential DSD
    # Each index holds the value of the midpoint diameter of the bin
    # in mm

    synthbins = np.linspace(0.1, 16.0, num=160)

    N_expDSD = []
    # refl_expDSD=[]
    # Construct the derived DSD
    for t in range(numtimes):
        # temp_N_expDSD = N0_exp[t] * np.exp(-lamda_exp[t] * synthbins / 1000.0)
        temp_N_expDSD = N0_exp[t] * np.exp(-lamda_exp[t] * avg_diameter / 1000.0)
        N_expDSD.append(temp_N_expDSD)

    N_expDSD = ma.array(N_expDSD)

    # --- Gamma Distribution Fit by Method of Moments (Tokay and Short 1996) ---

    # Now do the same as above for a gamma distribution fit

    # Calculate G, mu (shape parameter), lambda, and N0

    # Uncomment if you want estimates based on M3,M4,M6
#     G = (M4**3.)/((M3**2.)*M6)
#     G = ma.masked_invalid(G)
#     mu_gam = (11.*G-8.+(G*(G+8.))**(1./2.))/(2.*(1.-G))
#     mu_gam = ma.masked_invalid(mu_gam)
#     lamda_gam = (M3*(mu_gam+4.))/M4
#     lamda_gam = ma.masked_invalid(lamda_gam)
#     N0_gam = (M3*lamda_gam**(mu_gam+4.))/(special.gamma(mu_gam+4.))

    # Uncomment if you want estimates based on M2,M4,M6 (old/alternate masking method)
#     G = np.where((M2 == 0.0) | (M6 == 0.0), 0.0, (M4**2.)/(M2*M6))
#     mu_gam = np.where(G == 1.0, 0.0, ((7.-11.*G) -
#                                       ((7.-11.*G)**2. - 4.*(G-1.)*(30.*G-12.))
#                                       **(1./2.))/(2.*(G-1.)))
#     mu_gam = np.where(mu_gam <= -4., -3.99, mu_gam)
#     mu_gam = np.where(mu_gam > 30.,30.,mu_gam)
#     mu_gam = ma.masked_where(M4 is ma.masked,mu_gam)
#     lamda_gam = np.where(M4 == 0.0,0.0, ((M2*(mu_gam+3.)*(mu_gam+4.))/(M4))**(1./2.))
#     N0_gam = (M4*lamda_gam**(mu_gam+5.))/(special.gamma(mu_gam+5.))

    # Uncomment if you want estimates based on M2,M3,M4
#     mu_gam = (3.*M2*M4 - 4.*M3**2.)/(M3**2. - M2*M4)
#     mu_gam = ma.masked_invalid(mu_gam)
#     lamda_gam = (M3*(mu_gam+4.))/(M4)
#     lamda_gam = ma.masked_invalid(lamda_gam)
#     N0_gam = (M3*lamda_gam**(mu_gam+4.))/(special.gamma(mu_gam+4.))

# Uncomment if you want estimates based on M2,M4,M6 #######
    G = (M4**2.) / (M2 * M6)  # moment ratio based on untruncated moments
    G = ma.masked_invalid(G)
    mu_gam = ((7. - 11. * G) - ((7. - 11. * G)**2. - 4. * (G - 1.)
                                * (30. * G - 12.))**(1. / 2.)) / (2. * (G - 1.))
    mu_gam = ma.masked_where(mu_gam > 30., mu_gam)
    mu_gam = ma.masked_invalid(mu_gam)
    lamda_gam = ((M2 * (mu_gam + 3.) * (mu_gam + 4.)) / (M4))**(1. / 2.)
    mu_gam = ma.masked_where((lamda_gam > 20000.) | (mu_gam > 30.), mu_gam)
    mu_gam = ma.masked_invalid(mu_gam)
    lamda_gam = ma.masked_where(lamda_gam > 20000., lamda_gam)
    lamda_gam = ma.masked_invalid(lamda_gam)
    N0_gam = (M4 * lamda_gam**(mu_gam + 5.)) / (special.gamma(mu_gam + 5.))


# Truncated moments (TMF)
    lam_gam = lamda_gam
    mu_tmf = []
    lamda_tmf = []
    N0_tmf = []
    # LDmxlist = []

    for t in range(numtimes):
        LDmx = lam_gam[t] * Dmax[t]
        for x in range(10):
            mu = mu_gam[t]
            # truncated moment ratio below. Equation A8 from Thurai
            gm3 = gammap(3. + mu, LDmx) * np.exp(gammln(3. + mu))
            gm5 = gammap(5. + mu, LDmx) * np.exp(gammln(5. + mu))
            gm7 = gammap(7. + mu, LDmx) * np.exp(gammln(7. + mu))
            z0 = G[t] - gm5**2. / gm3 / gm7
            z1 = G[t] - gm5**2. / gm3 / gm7

            while(z1 / z0 > 0.0):
                mu = mu - 0.01
                gm3 = gammap(3. + mu, LDmx) * np.exp(gammln(3. + mu))
                gm5 = gammap(5. + mu, LDmx) * np.exp(gammln(5. + mu))
                gm7 = gammap(7. + mu, LDmx) * np.exp(gammln(7. + mu))
                z1 = G[t] - gm5**2. / gm3 / gm7

            lam_tmf = (M2[t] * gm5 / M4[t] / gm3)**0.5
            LDmx = lam_tmf * Dmax[t]
        mu_tmf.append(mu)
        lamda_tmf.append(lam_tmf)

    mu_tmf = np.array(mu_tmf)
    lamda_tmf = np.array(lamda_tmf)
    LDmx = lamda_tmf * Dmax
    N0_tmf = (M4 * lamda_tmf**(5. + mu_tmf)) / \
        (gammap(5. + mu_tmf, LDmx) * np.exp(gammln(5. + mu_tmf)))

##################################

    N_gamDSD = []
    N_tmfDSD = []
    # refl_gamDSD=[]

    for t in range(numtimes):
        # temp_N_gamDSD = N0_gam[t]*((synthbins/1000.0)**mu_gam[t])*
        #                 np.exp(-lamda_gam[t]*synthbins/1000.0)
        temp_N_gamDSD = N0_gam[t] * ((avg_diameter / 1000.0)**mu_gam[t]) * \
            np.exp(-lamda_gam[t] * avg_diameter / 1000.)    # Use disdrometer bins
        temp_N_tmfDSD = N0_tmf[t] * ((avg_diameter / 1000.0)**mu_tmf[t]) * \
            np.exp(-lamda_tmf[t] * avg_diameter / 1000.)
        # temp_Z_gamDSD = ((synthbins/1000.)**6.)*(1000.*temp_N_gamDSD)*0.1/1000.
        N_gamDSD.append(temp_N_gamDSD)
        N_tmfDSD.append(temp_N_tmfDSD)
        # refl_gamDSD.append(10.0*np.log10(np.sum(temp_Z_gamDSD)))
    N_gamDSD = ma.array(N_gamDSD, dtype=np.float64)
    N_tmfDSD = ma.array(N_tmfDSD, dtype=np.float64)

    # Quantities based on exponential distribution

    GR1 = special.gamma(1. + mu_exp)
    GR2 = special.gamma(4. + mu_exp)

    qr_exp = (cmr / rho) * N0_exp * GR2 / lamda_exp**(mu_exp + 4.)
    Ntr_exp = N0_exp * GR1 / lamda_exp**(mu_exp + 1.)
    # Compute reflectivity for DSD fit
    Gr_exp = ((6. + mu_exp) * (5. + mu_exp) * (4. + mu_exp)) / \
        ((3. + mu_exp) * (2. + mu_exp) * (1. + mu_exp))
    Zr_exp = ((1. / cmr)**2.) * Gr_exp * ((rho * qr_exp)**2.) / Ntr_exp
    refl_DSD_exp = 10.0 * np.log10(1.e18 * Zr_exp)
    # print 'reflectivity (exp DSD) = ',refl_DSD
    # print 'reflectivity (gam DSD) = ',refl_gamDSD

    # Median volume diameter for exponential distribution
    D_med_exp = np.where(lamda_exp == 0., np.nan, (3.67 / lamda_exp) * 1000.0)
    # Mass-weighted mean diameter for exp dist.
    D_m_exp = np.where(lamda_exp == 0., np.nan, (4. / lamda_exp) * 1000.0)

    # Quantities based on gamma distribution

    GR1 = special.gamma(1. + mu_gam)
    GR2 = special.gamma(4. + mu_gam)
    GR3 = special.gamma(4.67 + mu_gam)

    qr_gam = (cmr / rho) * N0_gam * GR2 / lamda_gam**(mu_gam + 4.)
    Ntr_gam = N0_gam * GR1 / lamda_gam**(mu_gam + 1.)
    LWC_gam = qr_gam * rho * 1000.  # g/m^3

    # Compute reflectivity for DSD fit
    Gr_gam = ((6. + mu_gam) * (5. + mu_gam) * (4. + mu_gam)) / \
        ((3. + mu_gam) * (2. + mu_gam) * (1. + mu_gam))
    Zr_gam = ((1. / cmr)**2.) * Gr_gam * ((rho * qr_gam)**2.) / Ntr_gam
    refl_DSD_gam = 10.0 * np.log10(1.e18 * Zr_gam)

    # print 'reflectivity (exp DSD) = ',refl_DSD
    # print 'reflectivity (gam DSD) = ',refl_gamDSD

    D_med_gam = np.where(lamda_gam == 0., np.nan, ((3.67 + mu_gam) / lamda_gam) *
                         1000.0)    # Median volume diameter for gamma distribution
    # Mass-weighted mean diameter for gam. dist.
    D_m_gam = np.where(lamda_gam == 0., np.nan, ((4. + mu_gam) / lamda_gam) * 1000.0)

    # Rain rate for gamma distribution
    # Note, the original coefficient for RR from the Zhang papers is 7.125 x 10^-3
    # The units of this coefficient are mm hr^-1 m^3 mm^-3.67. We need the coefficient
    # in units of mm hr^-1 m^3 m^-3.67, because our lamda is in units of m^-1.
    # Thus we have 7.125e-3*1000.^3.67 = 7.29096257 x 10^8
    rainrate_gam = 7.29096257e8 * N0_gam * GR3 / lamda_gam**(4.67 + mu_gam)

    # Quantities based on truncated gamma distribution

    GR1 = special.gamma(1. + mu_tmf)
    GR2 = special.gamma(4. + mu_tmf)
    GR3 = special.gamma(4.67 + mu_tmf)
    GR4 = special.gamma(7. + mu_tmf)
    IGR1 = gammap(1. + mu_tmf, LDmx) * GR1
    IGR2 = gammap(4. + mu_tmf, LDmx) * GR2
    IGR3 = gammap(4.67 + mu_tmf, LDmx) * GR3
    IGR4 = gammap(7. + mu_tmf, LDmx) * GR4

    Ntr_tmf = N0_tmf * IGR1 / lamda_tmf**(mu_tmf + 1.)
    TM3 = N0_tmf * lamda_tmf**-(mu_tmf + 4) * IGR2
    LWC_tmf = cmr * 1000. * TM3  # g/m^3
    qr_tmf = LWC_tmf / rho
    # Can't use Gr in form of untruncated gamma case. Compute Ztr directly using incomplete gamma
    # function instead
    Ztr_tmf = ((rho * qr_tmf)**2.) / (cmr**2. * Ntr_tmf) * (IGR4 * IGR1) / IGR2
    refl_DSD_tmf = 10.0 * np.log10(1.e18 * Ztr_tmf)
    rainrate_tmf = 7.29096257e8 * N0_tmf * IGR3 / lamda_tmf**(4.67 + mu_tmf)

    # FIXME: need to compute D0 for truncated gamma distribution a bit differently, since the above
    # gamma formula is for an untruncated distribution. One possibility is to compute it the same
    # way we do for the observed DSD, but using the discretized truncated gamma DSD.
    # D_med_tmf = np.where(lamda_tmf == 0., np.nan, ((3.67 + mu_tmf) / lamda_tmf) *
    #                     1000.0)    # Median volume diameter for gamma distribution

    # Create several tuples to pack the data, and then return them
    # NOTE: Consider updating these to namedtuples
    # FIXME: This is ridiculous. Use a dictionary!

    exp_DSD = (N_expDSD, N0_exp, lamda_exp, mu_exp, qr_exp, Ntr_exp, refl_DSD_exp, D_med_exp,
               D_m_exp)
    gam_DSD = (N_gamDSD, N0_gam, lamda_gam, mu_gam, qr_gam, Ntr_gam, refl_DSD_gam, D_med_gam,
               D_m_gam, LWC_gam, rainrate_gam)
    tmf_DSD = (N_tmfDSD, N0_tmf, lamda_tmf, mu_tmf, qr_tmf, Ntr_tmf, refl_DSD_tmf,
               LWC_tmf, rainrate_tmf)
    dis_DSD = (Nc_bin, D_med_disd, D_m_disd, D_mv_disd, D_ref_disd, QR_disd, refl_disd,
               LWC_disd, M0, rainrate)

    return synthbins, exp_DSD, gam_DSD, tmf_DSD, dis_DSD


def rad2DD(filename, dlocs):
    """Given an 88D CFRadial file and a list of disdrometer locations (tuples of lat,lon),
       compute the reflectivity at the disdrometer location from the lowest sweep."""

    average_gates = False        # True if reflectivity should be averaged in the closest 9 gates
    # False if just picking the value from the closest gate

    Cressman = False             # Perform a Cressman analysis on nearby radar gates to determine
    # reflectivity value

    roi = 1750.  # 750.                   # Radius of influence of Cressman analysis in m

    sweepfile_netcdf = Nio.open_file(filename)

    # Grab some needed variables from the netCDF file

    # numgates = sweepfile_netcdf.dimensions['range']     # Number of gates
    # numtimes = sweepfile_netcdf.dimensions['time']     # Number of times
    # numsweeps = sweepfile_netcdf.dimensions['sweep']     # Number of sweeps
    # sweep_start_ray_index = sweepfile_netcdf.variables['sweep_start_ray_index'][:]
    sweep_end_ray_index = sweepfile_netcdf.variables['sweep_end_ray_index'][:]
    try:
        ray_n_gates = sweepfile_netcdf.variables['ray_n_gates'][:]
        twoDarray = False
    except BaseException:
        print("No ray_n_gates in file, assuming 2D arrays (azimuth,range)")
        twoDarray = True

    rlat = sweepfile_netcdf.variables['latitude'].get_value()
    rlon = sweepfile_netcdf.variables['longitude'].get_value()
    ralt = sweepfile_netcdf.variables['altitude'].get_value()

    # print "Number of gates: ",numgates
    # print "Radar lat,lon,alt",rlat,rlon,ralt

    rlat = rlat * deg2rad
    rlon = rlon * deg2rad

    elevs = sweepfile_netcdf.variables['elevation'][:]

    el = elevs[0]  # Just pick first element for now

    # print "Elevation angle ",el

    el = el * deg2rad

    temp = sweepfile_netcdf.variables['range']    # range to center of each gate
    gatewidth = temp.meters_between_gates[0]  # gate spacing for each gate
    range = temp[:]
    range_start = temp[:] - gatewidth / 2.  # We also want range to start of each gate

    # print "Gatewidth ",gatewidth

    # Get the Azimuth info

    azimuth_rad = sweepfile_netcdf.variables['azimuth'][:]
    azimuth_rad = azimuth_rad[:sweep_end_ray_index[0]]  # Just grab first sweep for now
    # print "number of azimuths in sweep ",np.size(azimuth_rad)

    # Roll the azimuth dimension around so that it starts at 0 and ends at 360
    try:
        shift = np.where(azimuth_rad < azimuth_rad[0])[0][0]
    except BaseException:
        shift = 0
    # print "shift is ",shift
    azimuth_rad = np.roll(azimuth_rad, shift=-shift)

    beamwidth = azimuth_rad[1:] - azimuth_rad[0:-1]

    azimuth_start_rad = np.zeros(np.size(azimuth_rad) + 1)
    # Find azimuth of "start" of each gate (in azimuth) -- approximate
    azimuth_start_rad[1:-1] = azimuth_rad[1:] - 0.5 * beamwidth[:]
    azimuth_start_rad[0] = azimuth_rad[0] - 0.5 * beamwidth[0]
    azimuth_start_rad[-1] = azimuth_rad[-1] + 0.5 * beamwidth[-1]

    azimuth_rad = azimuth_rad * deg2rad
    azimuth_start_rad = azimuth_start_rad * deg2rad

    # Get the Reflectivity field
    tempdBZ = sweepfile_netcdf.variables['REF']
    dBZ = sweepfile_netcdf.variables['REF'][:]

    # Grab first sweep out of dBZ array
    # Need to check if the array is dimensioned by total number of gates or just number of gates for
    # each sweep

    # print ray_n_gates[0]
    # print sweep_end_ray_index[0]
    if(not twoDarray):  # Note, this needs to be tested!
        numpointsinsweep = ray_n_gates[0] * sweep_end_ray_index[0]
        dBZ = dBZ[:numpointsinsweep]
        dBZ = dBZ.reshape((sweep_end_ray_index[0], -1))
    else:
        dBZ = dBZ[:sweep_end_ray_index[0], :]

    # print "scale_factor,add_offset",tempdBZ.scale_factor[0],tempdBZ.add_offset[0]

    # Unpack values
    dBZ = dBZ * tempdBZ.scale_factor[0] + tempdBZ.add_offset[0]

    # Adjust azimuth axis

    dBZ = np.roll(dBZ, shift=-shift, axis=0)

    # print azimuth_rad.shape
    # print range.shape
    # print dBZ.shape

    theta, rad = np.meshgrid(azimuth_start_rad, range_start)
    theta_c, rad_c = np.meshgrid(azimuth_rad, range)
    xplt_c, yplt_c = oban.xyloc(rad_c, theta_c, el, rlat, rlon, ralt, 1)
    dBZ = dBZ.swapaxes(0, 1)
    dBZ = ma.masked_invalid(dBZ)
    # Find x,y location of disdrometer

    dBZ_D_list = []

    for dloc in dlocs:
        Dx, Dy = oban.ll_to_xy(dloc[0] * deg2rad, dloc[1] * deg2rad, rlat, rlon, 1)

        # print "Dx,Dy",Dx,Dy

        # Find closest gate to disdrometer location
        # First find the closest azimuth

#         theta_dis = -np.arctan2(Dy,Dx)+np.pi/2.       # I *think* this is correct
#
#         #print "azimuth of disdrometer: ",theta_dis/deg2rad
#         #Find closest index of azimuth
#
#         theta_diff = theta_dis-theta_c[0,:]
#         theta_index = np.argmin(np.abs(theta_diff))
#
#         # Now find index along range; first compute slant range
#
#         srange_dis = oban.computeslantrange([Dx],[Dy],el)
#         rad_diff = srange_dis-rad_c[:,0]
#         srange_index = np.argmin(np.abs(rad_diff))
#         #print "range of disdrometer: ",srange_dis

        # Try another way to get the indices (the above seems to have problems).
        # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why
        # the above doesn't work in all situations, but it sometimes picks a gate adjacent to the
        # one we want...

        if(Cressman):
            # print xplt_c.shape,yplt_c.shape
            dBZ_D = oban.Cresmn(xplt_c, yplt_c, dBZ, Dx, Dy, roi)
            # print dBZ_D,dBZ_D.shape
        # First, compute the euclidian distance of the x,y location of the disdrometer to each of
        # the radar gate centers
        else:
            distance = np.sqrt(((Dx - xplt_c)**2. + (Dy - yplt_c)**2.))

            # Now, find the index of the closest radar gate
            srange_index, theta_index = np.unravel_index(distance.argmin(), distance.shape)
            # print "srange_index,theta_index",srange_index,theta_index

            # Finally, grab reflectivity at closest gate to disdrometer
            # Average the reflectivity in the closest gate and surrounding 8 gates if desired
            if(not average_gates):
                dBZ_D = dBZ[srange_index, theta_index]
            else:
                dBZ_D = (1. / 9.) * (dBZ[srange_index - 1, theta_index - 1]
                                     + dBZ[srange_index, theta_index - 1]
                                     + dBZ[srange_index - 1, theta_index]
                                     + dBZ[srange_index, theta_index]
                                     + dBZ[srange_index + 1, theta_index]
                                     + dBZ[srange_index, theta_index + 1]
                                     + dBZ[srange_index + 1, theta_index + 1]
                                     + dBZ[srange_index + 1, theta_index - 1]
                                     + dBZ[srange_index - 1, theta_index + 1])

        # print "Reflectivity at disdrometer = ","%.2f"%dBZ_D

        dBZ_D_list.append(dBZ_D)

    # dBZ_D_arr = np.array(dBZ_D_list)

    sweepfile_netcdf.close()

    return dBZ_D_list


# def rad2loc_pyart(fieldlist, radarsweep, )


def rad2DD2(fieldlist, range_start, rrange, azimuth_start_rad, azimuth_rad, rlat, rlon, ralt, el,
            dlocs, average_gates=True, Cressman=False, roi=750., map_proj=1):
    """
    Another version of rad2DD: assumes radar sweep has been read in and computes values of fields in
    fieldlist at the disdrometer location(s).  Eventually will allow for optional
    advection/sedimentation correction of radar reflectivity to account for height above surface of
    radar scan.  Returns list of lists of field values at each disdrometer location
    """

    # First find x,y locations of radar gates

    xrad, yrad, xrad_c, yrad_c = radar.sweep2xy(
        azimuth_start_rad, azimuth_rad, range_start, rrange, el, rlat, rlon, ralt, map_proj)

    field_D_list = []
    dxy_list = []

    for dloc in dlocs:
        Dx, Dy = oban.ll_to_xy(dloc[0] * deg2rad, dloc[1] * deg2rad, rlat, rlon, map_proj)
        dxy_list.append((Dx, Dy))
        # print "Dx,Dy",Dx,Dy

        # Find closest gate to disdrometer location
        # First find the closest azimuth

#         theta_dis = -np.arctan2(Dy,Dx)+np.pi/2.       # I *think* this is correct
#
#         #print "azimuth of disdrometer: ",theta_dis/deg2rad
#         #Find closest index of azimuth
#
#         theta_diff = theta_dis-theta_c[0,:]
#         theta_index = np.argmin(np.abs(theta_diff))
#
#         # Now find index along range; first compute slant range
#
#         srange_dis = oban.computeslantrange([Dx],[Dy],el)
#         rad_diff = srange_dis-rad_c[:,0]
#         srange_index = np.argmin(np.abs(rad_diff))
#         #print "range of disdrometer: ",srange_dis

        # Try another way to get the indices (the above seems to have problems).
        # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why
        # the above doesn't work in all situations, but it sometimes picks a gate adjacent to the
        # onewe want...

        field_list = []

        for field in fieldlist:
            field = field.swapaxes(0, 1)
            field = np.ma.masked_invalid(field)

            if(Cressman):
                # print xplt_c.shape,yplt_c.shape
                field_D = oban.Cresmn(xrad_c, yrad_c, field, Dx, Dy, roi)
                # print dBZ_D,dBZ_D.shape
                # First, compute the euclidian distance of the x,y location of the disdrometer to
                # each of the radar gate centers
            else:
                distance = np.sqrt(((Dx - xrad_c)**2. + (Dy - yrad_c)**2.))

                # Now, find the index of the closest radar gate
                srange_index, theta_index = np.unravel_index(distance.argmin(), distance.shape)
                # print "srange_index,theta_index",srange_index,theta_index

                print("Distance to closest gate: ", distance[srange_index, theta_index])

                # Finally, grab field at closest gate to disdrometer
                # Average the field in the closest gate and surrounding 8 gates if desired
                # Disdrometer is probably not covered by sweep, set value to np.nan
                if(distance[srange_index, theta_index] > 3000.):
                    field_D = np.nan
                else:
                    if(not average_gates):
                        field_D = field[srange_index, theta_index]
                    else:
                        # field_D = (1./9.)*(field[srange_index-1,theta_index-1]+
                        #                    field[srange_index,theta_index-1]+
                        #                    field[srange_index-1,theta_index]+
                        #                    field[srange_index,theta_index]+
                        #                    field[srange_index+1,theta_index]+
                        #                    field[srange_index,theta_index+1]+
                        #                    field[srange_index+1,theta_index+1]+
                        #                    field[srange_index+1,theta_index-1]+
                        #                    field[srange_index-1,theta_index+1])
                        field_D = np.nanmean(field[srange_index - 1:srange_index + 2,
                                             theta_index - 1:theta_index + 2])

            # print "Value of field at disdrometer = ","%.2f"%field_D

            field_list.append(field_D)

        field_D_list.append(field_list)

    field_D_arr = np.array(field_D_list)

    return dxy_list, field_D_arr


def avgwind(winddirs, windspds, avgintv, gusts=True, gustintv=3, center=True):
    """Given a timeseries of wind directions and speeds, and an interval for averaging,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""

    windspdsavg = pd.Series(windspds).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    if(gusts):
        windgusts = pd.Series(windspds).rolling(
            window=gustintv,
            center=center,
            min_periods=1).mean().values
        windgustsavg = pd.Series(windgusts).rolling(
            window=avgintv,
            center=center,
            min_periods=1).max().values
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds * np.cos(np.deg2rad(-winddirs + 270.))
    vs = windspds * np.sin(np.deg2rad(-winddirs + 270.))

    # Linearly interpolate for bad values of us,vs
    us = interpnan1D(us)
    vs = interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.Series(us).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    vsavg = pd.Series(vs).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    windspdsavgvec = np.sqrt(usavg**2. + vsavg**2.)
    winddirsavgvec = (270.0 - (180. / np.pi) * np.arctan2(vsavg, usavg)
                      ) % 360.  # Need to use %360 to keep wind dir between 0 and 360 degrees

    return windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg


def resamplewind(datetimes, offset, winddirs, windspds, intervalstr, gusts=True, gustintvstr='3S',
                 center=False):
    """Given a timeseries of wind directions and speeds, and an interval for resampling,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""

    windspdsavg = pd.Series(data=windspds, index=datetimes).resample(intervalstr, label='right',
                                                                     closed='right',
                                                                     base=offset).mean()
    if(gusts):
        windgusts = pd.Series(data=windspds, index=datetimes).resample(gustintvstr, label='right',
                                                                       closed='right',
                                                                       base=offset).mean()
        windgustsavg = windgusts.resample(intervalstr, label='right', closed='right',
                                          base=offset).max()
        windgusts = windgusts
        windgustsavg = windgustsavg
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds * np.cos(np.deg2rad(-winddirs + 270.))
    vs = windspds * np.sin(np.deg2rad(-winddirs + 270.))

    # Linearly interpolate for bad values of us,vs
    # us = interpnan1D(us)
    # vs = interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.Series(
        data=us,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    vsavg = pd.Series(
        data=vs,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    windspdsavgvec = np.sqrt(usavg**2. + vsavg**2.)
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    winddirsavgvec = (270.0 - (180. / np.pi) * np.arctan2(vsavg, usavg)) % 360.

    # unit average wind direction
    unit_us = us / windspds
    unit_vs = vs / windspds
    unit_usavg = pd.Series(
        data=unit_us,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    unit_vsavg = pd.Series(
        data=unit_vs,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    winddirsunitavgvec = (270.0 - (180. / np.pi) * np.arctan2(unit_vsavg, unit_usavg)) % 360.

    # Pack everything into a dictionary for return

    wind_dict = {'windspdavg': windspdsavg, 'windspdavgvec': windspdsavgvec,
                 'winddiravgvec': winddirsavgvec, 'winddirunitavgvec': winddirsunitavgvec,
                 'windgust': windgusts, 'windgustavg': windgustsavg, 'uavg': usavg, 'vavg': vsavg,
                 'unit_uavg': unit_usavg, 'unit_vavg': unit_vsavg}

    return wind_dict

# -----------------------------------------------------------------------------------------
#
# Below are several routines for reading and analyzing VORTEX-2 disdrometer probe data.
# Eventually will merge this with the current pyPIPS code package.  There are kept here
# for reference.
#
# -----------------------------------------------------------------------------------------


def readdataCU_2010(filename):
    """This function reads in disdrometer data from the CU V2 2010 mission and computes number
       concentration from the particle counts and fall velocities for each bin"""

    disfile = open(filename, 'r')

    # Create V-D relationship for rain based on Terry Schuur's relationship
    rainvd = assignfallspeed(avg_diameter)

    dates = []
    times = []
    intensities = []
    preciptots = []
    weathercodes1 = []
    weathercodes2 = []
    weathercodes3 = []
    weathercodes4 = []
    reflectivities = []
    visibilities = []
    sampleintervals = []
    pcounts = []
    pcounts2 = []
    sensortemps = []
    # Will contain the number of drops in the diameter vs. fall speed matrix for each time
    # i.e. a (numtimes,32,32) array
    countsMatrix = []
    # Will contain computed number concentrations for each size bin for each time
    # i.e. a (numtimes,32) array
    concentrations = []

    ntimes = 0
    sampleinterval = -1  # Have to give it some initial value
    for l in disfile:

        # Read in ATM4 number
        try:
            number = np.int(l[:2])
        except BaseException:
            number = -1  # Blank line
        try:
            line = l[3:]
        except BaseException:
            line = None  # Blank line

        line = line.strip()

        if(line == ''):
            line = 0

        if(number == 1):
            intensities.append(np.float(line))
        elif(number == 2):
            preciptots.append(np.float(line))
        elif(number == 3):
            weathercodes1.append(np.int(line))
        elif(number == 4):
            weathercodes2.append(np.int(line))
        elif(number == 5):
            weathercodes3.append(line)
        elif(number == 6):
            weathercodes4.append(line)
        elif(number == 7):
            reflectivities.append(np.float(line))
        elif(number == 8):
            visibilities.append(np.float(line))
        elif(number == 9):
            sampleinterval = np.int(line)
            # Don't read records that have a sample interval that isn't 10 seconds
            if(sampleinterval != 10):
                print("sample interval = ", sampleinterval)
                # Pop the values from the end of the previous lists
                intensities.pop()
                preciptots.pop()
                weathercodes1.pop()
                weathercodes2.pop()
                weathercodes3.pop()
                weathercodes4.pop()
                reflectivities.pop()
                visibilities.pop()
            else:
                sampleintervals.append(sampleinterval)
        if(sampleinterval == 10):  # Read the rest of the data record if sampleinterval = 10 s
            if(number == 11):
                pcount = np.float(line)
                try:
                    pcount = np.int(pcount)
                except BaseException:
                    pcount = 0
                pcounts.append(pcount)
            elif(number == 12):
                sensortemps.append(np.float(line))
            elif(number == 20):
                times.append(line)
                ntimes = ntimes + 1
                # print "time,ntimes = ",times[ntimes-1],ntimes
            elif(number == 21):
                dates.append(line)
            elif(number == 93):
                line = line.strip().split(';')
                spectrum = [float(x) if x != '' else 0 for x in line]
                spectrum = spectrum[:-1]  # Strip off bogus last value

                # Now create an array out of the spectrum and reshape it to 32x32
                spectrum = np.array(spectrum, dtype='int')
                if(spectrum.size > 0):
                    spectrum = spectrum.reshape((32, 32))
                else:
                    spectrum = np.zeros((32, 32), dtype='int')

                # Append spectrum (corrected or not) to spectrum time list

                countsMatrix.append(spectrum)

    # Recast countsMatrix as numpy array

    countsMatrix = np.dstack(countsMatrix)
    countsMatrix = np.rollaxis(countsMatrix, 2, 0)

    # Perform Katja's QC routines if desired (should not be used in combination with my methods
    # above, most are redundant anyway).

    X, Y = np.meshgrid(avg_diameter, fall_bins)
    flaggedtimes = np.zeros(len(times), dtype=bool)

    if(use_strongwindQC):
        countsMatrix, flaggedtimes = strongwindQC(countsMatrix)

    if(use_splashingQC):
        countsMatrix = splashingQC(countsMatrix)

    if(use_marginQC):
        countsMatrix = marginQC(countsMatrix)

    if(use_rainfallspeedQC):
        countsMatrix = rainfallspeedQC(countsMatrix, rainvd, falltol, maskhigh, masklow)

    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)

    if(use_rainonlyQC):
        countsMatrix = rainonlyQC(countsMatrix)

    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)

    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)

    # print flaggedtimes

    if(plot_QC):
        for t in range(np.size(countsMatrix, axis=0)):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plt.title('Fall speed vs. diameter for time ' + times[t])

            countsplot = ma.masked_where(countsMatrix[t, :] <= 0, countsMatrix[t, :])

            C = ax1.pcolor(
                min_diameter,
                min_fall_bins,
                countsplot,
                vmin=1,
                vmax=50,
                edgecolors='w')
            ax1.plot(avg_diameter, rainvd, c='r')
            # ax1.scatter(X[0:10,20:31],Y[0:10,20:31],c='r',marker='x')
            fig.colorbar(C)

            if(len(flaggedtimes) > 0 and flaggedtimes[t]):
                ax1.text(0.5, 0.5, 'Flagged for strong wind contamination!',
                         horizontalalignment='center',
                         verticalalignment='center', color='y',
                         transform=ax1.transAxes)
            if(plot_strongwindQC):
                ax1.scatter(X[strongwindmask], Y[strongwindmask], c='r', marker='x', alpha=1.0)
            # FIXME: figure out what happened to the splashmask
            # if(plot_splashingQC):
            #     ax1.scatter(X[splashmask], Y[splashmask], c='w', marker='o', alpha=0.75)
            #     # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(splashmask,
            #                  mask=-splashmask),cmap=cm.Reds,alpha=0.1)
            if(plot_marginQC):
                ax1.scatter(X[marginmask], Y[marginmask], c='g', marker='x', alpha=0.1)
                # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(marginmask,
                #            mask=-marginmask),cmap=cm.Reds,alpha=0.1)
            if(plot_rainfallspeedQC):
                ax1.scatter(X[fallspeedmask], Y[fallspeedmask], c='k', marker='x', alpha=0.5)
                # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(fallspeedmask,
                #            mask=-fallspeedmask),cmap=cm.gray,alpha=0.1)
            if(plot_rainonlyQC):
                ax1.scatter(X[rainonlymask], Y[rainonlymask], c='g', marker='x', alpha=0.5)

            ax1.set_xlim(0.0, 9.0)
            ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
            ax1.set_xlabel('diameter (mm)')
            ax1.set_ylim(0.0, 15.0)
            ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
            ax1.set_ylabel('fall speed (m/s)')

            plt.savefig('/Users/ddawson/Dropbox/temp6/temp_' + times[t] + '.png')

    # Find total number of non-masked particles

    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)

    # print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t, time in enumerate(times):
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy, vspectrum = np.meshgrid(avg_diameter, fall_bins)
            dspectrum = countsMatrix[t, :]
        else:
            vspectrum = rainvd
            dspectrum = countsMatrix[t, :].sum(axis=0)   # Sum up particles for each diameter bin

            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]

        # Now compute the number concentration using the assumed fall speeds, sensor area, and
        # sampling interval
        # Units are #/m^3/mm
        if(spectrum.size > 0 and not flaggedtimes[t]):
            concentration = dspectrum / \
                (vspectrum * sampleintervals[t] * eff_sensor_area * (max_diameter - min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
            # print "concentration.shape"
            # print concentration.shape
        elif(not flaggedtimes[t]):
            concentration = np.zeros_like(avg_diameter)
        else:
            concentration = -999. * np.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999., concentration)

        # Throw out particles above and below a certain diameter if desired

#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)

        # print "Number of particles counted vs. summed number: "
        # print np.int(line[7]),np.sum(dspectrum)

        # print concentration

        concentrations.append(concentration)

    concentrations = ma.array(concentrations)
    # print "concentrations: ",concentrations
    pcounts = np.array(pcounts)
    # pcounts2 = ma.array(pcounts2)

    return (dates, times, intensities, preciptots, weathercodes1, reflectivities, visibilities,
            pcounts, pcounts2, sensortemps, min_diameter, max_diameter, avg_diameter,
            concentrations)


def readNV2netCDF(filename_conv, filename_dsd, fixGPS=True, basicqc=False, rainfallqc=False,
                  rainonlyqc=False, hailonlyqc=False, strongwindqc=False, detecthail=True,
                  requested_interval=60.0, starttime=None, stoptime=None, matchtimes=True,
                  UTC=True):
    """Reads NSSL V2 disdrometer data in netCDF format."""

    sampling_interval = 60.  # Always the case for Glen's data

    # Create V-D relationship for rain based on Terry Schuur's relationship
    rainvd = assignfallspeed(avg_diameter)

    pdatetimes = []
    intensities = []
    # preciptots = []
    # reflectivities = []
    # sampleintervals = []
    # pcounts = []
    # pcounts2 = []
    # sensortemps = []
    # amplitudes = []
    # voltages = []
    # pvoltages = []
    # weathercodes = []
    # visibilities = []

    # Air density at the PSD times (needed for rain fall speed curve correction)
    rhoatpdatetimes = []
    onedrop_concentrations = []

    # dates = []
    # times = []
    datetimes = []
    # bwindspds = []
    # bwinddirabss = []
    swindspds = []
    swinddirabss = []
    # swinddiags = []
    # fasttemps = []
    slowtemps = []
    dewpoints = []
    pts = []
    qvs = []
    rhos = []
    # RHs_derived = []
    RHs = []
    pressures = []
    # compass_dirs = []
    # GPS_lats = []
    # GPS_lons = []
    # GPS_stats = []
    # GPS_alts = []
    # GPS_speeds = []
    ptvs = []
    ptes = []

    # First, read in the disdrometer data

    disfile = netcdf.Dataset(filename_dsd, 'r')

    # Read in the date and time information
    year_start = disfile.variables['year1'][:]
    year_end = disfile.variables['year2'][:]
    # month_start = disfile.variables['month1'][:]
    month_end = disfile.variables['month2'][:]
    # day_start = disfile.variables['day1'][:]
    day_end = disfile.variables['day2'][:]
    # hour_start = disfile.variables['hour1'][:]
    hour_end = disfile.variables['hour2'][:]
    # min_start = disfile.variables['minute1'][:]
    min_end = disfile.variables['minute2'][:]
    # sec_start = disfile.variables['second1'][:]
    sec_end = disfile.variables['second2'][:]

    # Construct the datetime objects (turns out starts and ends are the same time -- this is a bug
    # that will be corrected later -- the end times are correct and represent the
    # end of each 1-min DSD collection)

    # disdate_start = []
    # disdate_end = []

    # dis_startindex = 0
    # dis_endindex = len(year_end)

    for i in range(np.size(year_start)):
        pdatetime = datetime(year_end[i], month_end[i], day_end[i], hour_end[i], min_end[i],
                             sec_end[i])
        # Times in the NSSL V2 probe files are in CST. Convert to UTC here
        if UTC:
            pdatetime = pdatetime + timedelta(hours=6)
        pdatetimes.append(pdatetime)

    #         # Find index corresponding to chosen start and endtimes
    #         if(date2num(pdatetime) == starttime):
    #             dis_startindex = i
    #         if(date2num(pdatetime) == endtime):
    #             dis_endindex = i

    #     disdates_1min = disdate_end[dis_startindex:dis_endindex + 1]

    # Read in the bin size information
    # NOTE: commented out for now since this information is the same for all Parsivels
    # TODO: in the future, maybe read this in and check that it matches the bins defined
    # in the top level of the module.

    #     min_size = disfile.variables['min_size'][:]
    #     max_size = disfile.variables['max_size'][:]
    #     avg_size = (min_size + max_size) / 2.
    bin_width = max_diameter - min_diameter

    # Read in the concentration information
    # It appears that the data is in log (base 10) of number concentration.
    # Need to confirm with Glen!  Status: confirmed

    logND = disfile.variables['concentration'][:]
    ND = np.where(logND == 0.0, 0.0, 10.**(logND))
    logND = np.ma.masked_where(ND <= 0.0, logND)
#     ND = ND.swapaxes(0, 1)
#     logND = logND.swapaxes(0, 1)

    # Restrict time dimension to start and end times
    #     Nc_bin = Nc_bin[:, dis_startindex:dis_endindex + 1]
    #     logNc_bin = logNc_bin[:, dis_startindex:dis_endindex + 1]

    # Read in the total particle count for each 1-min interval
    pcount = disfile.variables['pcount'][:]
    #     pcount = pcount[dis_startindex:dis_endindex + 1]

    # Invert the particle counts from the # conc. in each bin
    # TODO: Update this to take into account surface air density. Not high priority since we
    # aren't really using this information right now.
    counts_bin = np.zeros_like(ND)
    for t in range(np.size(counts_bin, axis=0)):
        counts_bin[t, :] = ND[t, :] * rainvd[:] * \
            sampling_interval * eff_sensor_area * bin_width[np.newaxis, :]

    # Zero out bins outside desired size range if desired
    # TODO: are we using this?
#     if(maskhighdiam):
#         diamindex = np.where(avg_diameter > highdiamthresh)[0][0]
#         Nc_bin[diamindex:, :] = 0.0
#         logNc_bin = np.ma.masked_where(Nc_bin <= 0.0, logNc_bin)
#
#     if(masklowdiam):
#         diamindex = np.where(avg_diameter > lowdiamthresh)[0][0]
#         Nc_bin[:diamindex, :] = 0.0
#         logNc_bin = np.ma.masked_where(Nc_bin <= 0.0, logNc_bin)

    intensities = disfile.variables['rain_rate'][:]

    # Read in latitude and longitude of the probe (too bad altitude isn't here, too)
    latitude = disfile.variables['latitude'][...]
    longitude = disfile.variables['longitude'][...]

    disfile.close()

    # Now read in the conventional data

    convfile = netcdf.Dataset(filename_conv, 'r')

    # First read in the start date:
    datestring = convfile.start_date
    year = int(datestring[:4])
    month = int(datestring[5:7])
    day = int(datestring[8:10])
    hour = 0
    min = 0
    sec = 0
    datetimestart = datetime(year, month, day, hour, min, sec)

    # Read in the time variable, which contains the time in seconds since the start time for
    # each record and then convert to datetime objects
    # TODO: see if there's a way to easily vectorize this operation

    timestampssec = convfile.variables['time'][:].tolist()

    datetimes = []
    for i, second in enumerate(timestampssec):
        # Append values to their respective lists
        tempdatetime = datetimestart + timedelta(seconds=second)
        if UTC:
            tempdatetime = tempdatetime + timedelta(hours=6)
        datetimes.append(tempdatetime)

    slowtemps = convfile.variables['tmpc'][:]
    RHs = convfile.variables['relh'][:]
    swindspds = convfile.variables['wspd'][:]
    swinddirabss = np.mod(convfile.variables['wdir'][:] - convfile.sta_heading, 360.)
    swindgusts = convfile.variables['gust'][:]
    dewpoints = convfile.variables['dewc'][:]
    pressures = convfile.variables['pres'][:]
    pts = convfile.variables['thet'][:]
    ptvs = convfile.variables['thtv'][:]
    ptes = convfile.variables['thte'][:]
    qvs = convfile.variables['qvap'][:]
    rhos = thermo.calrho(pressures * 100., pts, qvs)

    convfile.close()

    # Find indices in conventional time list that match those in the Parsivel time list
    pindices = np.searchsorted(datetimes, pdatetimes)
    rhoatpdatetimes = rhos[pindices]

    # Note that Glen's data includes the internally computed number concentrations in each
    # diameter bin instead of the full velocity-diameter matrix, so we won't be able to use
    # any of the QC routines. We still compute the "One-drop" concentrations here.

    counts_1drop = np.ones_like(avg_diameter)
    for t, time in enumerate(pdatetimes):
        onedropvel = assignfallspeed(avg_diameter, rhocorrect=True, rho=rhoatpdatetimes[t])
        onedrop_concentration = counts_1drop / \
            (onedropvel * sampling_interval * eff_sensor_area * (max_diameter - min_diameter))
        onedrop_concentrations.append(onedrop_concentration)

    # Not really "corrected" times but these are used for symmetry with readCU and readPIPS
    # TODO: Make a unified reader routine to consolidate code?

    datetimes_corrected = pd.Index(datetimes)
    pdatetimes_corrected = pd.Index(pdatetimes)

    # Create a dataframe out of the conventional data (For Glen's data it's in 10-s intervals)
    # First create a dictionary to organize the conventional data
    convdict = {'swindspd': swindspds, 'swinddirabs': swinddirabss, 'swindgust': swindgusts,
                'slowtemp': slowtemps, 'RH': RHs, 'pressure': pressures, 'latitude': latitude,
                'longitude': longitude, 'dewpoint': dewpoints, 'pt': pts, 'qv': qvs, 'rho': rhos,
                'ptv': ptvs, 'pte': ptes}
    conv_df = pd.DataFrame(convdict, index=datetimes_corrected)

    # Create one dataframe for N(D), one for the "one-drop" N(D), and one for all the
    # other 10-s data

    ND_df = pd.DataFrame(data=ND, index=pdatetimes_corrected, columns=avg_diameter)
    ND_onedrop_df = pd.DataFrame(data=onedrop_concentrations, index=pdatetimes_corrected,
                                 columns=avg_diameter)

    # Create a dictionary to organize the 10-s data
    PSDdict = {'intensity': intensities, 'pcount': pcount}
    PSD_df = pd.DataFrame(PSDdict, index=pdatetimes_corrected)

    # If desired, clean up the Parsivel data so that the start and end times are within the
    # datetimes index. This option is recommended!
    if matchtimes:
        pdatetimes_corrected = pdatetimes_corrected.intersection(datetimes_corrected)
        ND_df = ND_df.reindex(pdatetimes_corrected)
        ND_onedrop_df = ND_onedrop_df.reindex(pdatetimes_corrected)
        PSD_df = PSD_df.reindex(pdatetimes_corrected)

    # Average and thin the DSD data with the desired interval

    DSD_interval, intervalstr, ND_df, ND_onedrop_df, PSD_df, __ = \
        resamplePSD(requested_interval, sampling_interval, ND_df, ND_onedrop_df, PSD_df)

    # Pandas apparently (gotcha!) converts missing values to NaN when extracting the numpy array
    # representation using .values. Since I want the original masked array functionality for now for
    # further computations, I need to remask the array here.
    # Otherwise, the NaN's propagate in further computations...
    # In the future, another solution that uses Pandas more natively should be
    # pursued, but this will work for now
    ND = ND_df.values
    mask = ND_df.isnull()
    ND = ma.array(ND, mask=mask)
    ND_onedrop = ND_onedrop_df.values
    mask = ND_onedrop_df.isnull()
    ND_onedrop = ma.array(ND_onedrop, mask=mask)
    # Argh, have to convert back to datetime objects.  This one from
    # http://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
    datetimes_corrected = datetimes_corrected.to_pydatetime()
    pdatetimes_corrected = ND_df.index.to_pydatetime()
    DSD_index = ND_df.index

    # Aggregate everything into a dictionary to return

    PIPS_dict = {'convtimestamps': datetimes_corrected, 'PSDtimestamps': pdatetimes_corrected,
                 'DSD_index': DSD_index, 'DSD_interval': DSD_interval, 'intervalstr': intervalstr,
                 'conv_df': conv_df, 'PSD_df': PSD_df, 'ND_df': ND_df,
                 'ND_onedrop_df': ND_onedrop_df, 'ND': ND, 'ND_onedrop': ND_onedrop}

    return PIPS_dict


def readdataGlennetCDF_old(filename, starttime=None, endtime=None):
    """Reads disdrometer data in netCDF format for Glens data."""

    # Create V-D relationship for rain based on Terry Schuur's relationship
    rainvd = assignfallspeed(avg_diameter)

    dis_file = netcdf.Dataset(filename, 'r')

    # Read in the date and time information
    year_start = dis_file.variables['year1'][:]
    year_end = dis_file.variables['year2'][:]
    month_start = dis_file.variables['month1'][:]
    month_end = dis_file.variables['month2'][:]
    day_start = dis_file.variables['day1'][:]
    day_end = dis_file.variables['day2'][:]
    hour_start = dis_file.variables['hour1'][:]
    hour_end = dis_file.variables['hour2'][:]
    min_start = dis_file.variables['minute1'][:]
    min_end = dis_file.variables['minute2'][:]
    sec_start = dis_file.variables['second1'][:]
    sec_end = dis_file.variables['second2'][:]

    # Construct the datetime objects (turns out starts and ends are the same time -- this is a bug
    # that will be corrected later -- the end times are correct and represent the
    # end of each 1-min DSD collection)

    disdate_start = []
    disdate_end = []

    dis_startindex = 0
    dis_endindex = len(year_end)

    for i in range(np.size(year_start)):
        distimestart = datetime(
            year_start[i],
            month_start[i],
            day_start[i],
            hour_start[i],
            min_start[i],
            sec_start[i])
        distimeend = datetime(
            year_end[i],
            month_end[i],
            day_end[i],
            hour_end[i],
            min_end[i],
            sec_end[i])
        # Find index corresponding to chosen start and endtimes
        if(dates.date2num(distimeend) == starttime):
            dis_startindex = i
        if(dates.date2num(distimeend) == endtime):
            dis_endindex = i
        disdate_start.append(distimestart)
        disdate_end.append(distimeend)

    disdates_1min = disdate_end[dis_startindex:dis_endindex + 1]

    # Read in the bin size information

    min_size = dis_file.variables['min_size'][:]
    max_size = dis_file.variables['max_size'][:]
    avg_size = (min_size + max_size) / 2.
    bin_width = max_size - min_size

    # Read in the concentration information
    # It appears that the data is in log (base 10) of number concentration.
    # Need to confirm with Glen!  Status: confirmed

    logNc_bin = dis_file.variables['concentration'][:]
    Nc_bin = np.where(logNc_bin == 0.0, 0.0, 10.**(logNc_bin))
    logNc_bin = np.ma.masked_where(Nc_bin <= 0.0, logNc_bin)
    Nc_bin = Nc_bin.swapaxes(0, 1)
    logNc_bin = logNc_bin.swapaxes(0, 1)

    # Restrict time dimension to start and end times
    Nc_bin = Nc_bin[:, dis_startindex:dis_endindex + 1]
    logNc_bin = logNc_bin[:, dis_startindex:dis_endindex + 1]

    # Read in the total particle count for each 1-min interval
    pcount = dis_file.variables['pcount'][:]
    pcount = pcount[dis_startindex:dis_endindex + 1]

    # Invert the particle counts from the # conc. in each bin
    counts_bin = np.zeros_like(Nc_bin)
    for t in range(np.size(counts_bin, axis=1)):
        counts_bin[:, t] = Nc_bin[:, t] * rainvd[:] * \
            sampling_period * eff_sensor_area * bin_width[:]

    # Zero out bins outside desired size range if desired
    if(maskhighdiam):
        diamindex = np.where(avg_diameter > highdiamthresh)[0][0]
        Nc_bin[diamindex:, :] = 0.0
        logNc_bin = np.ma.masked_where(Nc_bin <= 0.0, logNc_bin)

    if(masklowdiam):
        diamindex = np.where(avg_diameter > lowdiamthresh)[0][0]
        Nc_bin[:diamindex, :] = 0.0
        logNc_bin = np.ma.masked_where(Nc_bin <= 0.0, logNc_bin)

    return (pcount, counts_bin, Nc_bin, logNc_bin, min_size, avg_size, max_size, bin_width,
            disdates_1min)


def readthermoNV2(filename):
    """Reads thermodynamic data from a NSSL V2 disdrometer probe"""

    tprh_file = open(filename, 'r')

    # Throw out the first 4 lines, which are just header information

    for line in range(4):
        tprh_file.readline()
    # print dummyline

    # Now, read in the data and put them into numpy arrays

    dates = []
    records = []
    temps = []
    rhs = []
    pressures = []

    for line in tprh_file:
        tokens = line.strip().split(',')
        dates.append(tokens[0])
        records.append(np.int(tokens[1]))
        temps.append(np.float(tokens[2]))
        rhs.append(np.float(tokens[3]))
        pressures.append(np.float(tokens[4]))

    temps = np.array(temps)
    rhs = np.array(rhs)
    pressures = np.array(pressures)

    # Dates are in the following format: YYYY-MM-DD HH:MM:SS
    # where the time is CST (-6 GMT)
    # We can use the datetime module to convert the strings to a list of datetime objects
    # For some reason I can't get the pytz module to specify the time zone as CST: it always
    # converts the time zone to CDT.  Sigh...

    thermodates = []

    for date in dates:
        year = int(date[1:5])
        month = int(date[6:8])
        day = int(date[9:11])
        hour = int(date[12:14])
        minute = int(date[15:17])
        second = int(date[18:20])
        # Now that we have all the date and time elements, we can create a datetime object
        thermotime = datetime(year, month, day, hour, minute, second)
        # newtime1=Central.localize(newtime)
        # print newtime.strftime(fmt)
        thermodates.append(thermotime)

    return thermodates, temps, rhs, pressures


def readwindNV2(filename):
    """Reads wind data from one of the NSSL V2 disdrometer probes"""

    wind_file = open(filename, 'r')

    # Throw out the first 4 lines, which are just header information

    for line in range(4):
        wind_file.readline()
    # print dummyline

    # Now, read in the data and put them into numpy arrays

    dates = []
    records = []
    winds = []
    winddirs = []
    windgusts = []

    for line in wind_file:
        tokens = line.strip().split(',')
        dates.append(tokens[0])
        records.append(np.int(tokens[1]))
        winds.append(np.float(tokens[2]))
        winddirs.append(np.float(tokens[3]))
        windgusts.append(np.float(tokens[5]))

    winds = np.array(winds)
    winddirs = np.array(winddirs)
    windgusts = np.array(windgusts)

    # Dates are in the following format: YYYY-MM-DD HH:MM:SS
    # where the time is CST (-6 GMT)
    # We can use the datetime module to convert the strings to a list of datetime objects
    # For some reason I can't get the pytz module to specify the time zone as CST: it always
    # converts the time zone to CDT.  Sigh...

    winddates = []

    for date in dates:
        year = int(date[1:5])
        month = int(date[6:8])
        day = int(date[9:11])
        hour = int(date[12:14])
        minute = int(date[15:17])
        second = int(date[18:20])
        # Now that we have all the date and time elements, we can create a datetime object
        windtime = datetime(year, month, day, hour, minute, second)
        # newtime1=Central.localize(newtime)
        # print newtime.strftime(fmt)
        winddates.append(windtime)

    return winddates, winds, winddirs, windgusts
