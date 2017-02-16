# disdrometer_module.py
# A collection of functions for reading VORTEX 2 disdrometer data

import netCDF4 as netcdf
import Nio
import numpy as N
from numpy import ma as ma
from matplotlib.dates import *
import matplotlib.cm as cm
from datetime import datetime,timedelta
import pytz as pytz
import os
import shlex
import matplotlib.pyplot as plt
from matplotlib.ticker import *
from . import obanmodule as oban
from . import radarmodule as radar
from scipy import special
import pdb
import pandas as pd

deg2rad = N.pi/180.

parsivel_ids = ['304545','295153','295166','304543'] # ID #'s for PIPS 1A,1B,2A,2B

# Min diameter of bins (mm)
min_diameter_bins = [0.000,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.000,1.125,1.250,1.500,1.750,
                     2.000,2.250,2.500,3.000,3.500,4.000,4.500,5.000,6.000,7.000,8.000,9.000,10.000,
                     12.000,14.000,16.000,18.000,20.000,23.000]
max_diameter_bins = [0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.000,1.125,1.250,1.500,1.750,2.000,
                     2.250,2.500,3.000,3.500,4.000,4.500,5.000,6.000,7.000,8.000,9.000,10.000,12.000,
                     14.000,16.000,18.000,20.000,23.000,26.000]
                     
# Average diameter of bins (mm)
avg_diameter_bins = [0.5*(x+y) for x,y in zip(min_diameter_bins,max_diameter_bins)]


fall_bins = [0.050,0.150,0.250,0.350,0.450,0.550,0.650,0.750,0.850,0.950,1.100,1.300,1.500,1.700,
             1.900,2.200,2.600,3.000,3.400,3.800,4.400,5.200,6.000,6.800,7.600,8.800,10.400,12.000,
             13.600,15.200,17.600,20.800]

min_fall_bins = [0.000,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.000,1.200,1.400,
                 1.600,1.800,2.000,2.400,2.800,3.200,3.600,4.000,4.800,5.600,6.400,7.200,8.000,
                 9.600,11.200,12.800,14.400,16.000,19.200]

# Parsivel sampling area and sampling period
sensor_area = 0.0054    # (m^2)
sampling_period = 10.0  # (s)

min_diameter = N.array(min_diameter_bins)
avg_diameter = N.array(avg_diameter_bins)
max_diameter = N.array(max_diameter_bins)
fall_bins = N.array(fall_bins)
min_fall_bins = N.array(min_fall_bins)

eff_sensor_area = (180.*(30.-avg_diameter/2.))*1.e-6 # Jaffrain and Berne (2011), Tokay et al. (2011)


use_DD_QC = False   # Use my QC methods (see below)

use_measured_fs = True     # If True, use the raw measured fall speeds to compute number concentration
                            # If False, use the fall speed curve of Terry Schuur

use_strongwindQC = False      # Remove time records that are contaminated by strong wind?
use_splashingQC = False      # Remove drops that result from splashing?
use_marginQC = False         # Remove drops that result from margin falls?
use_rainonlyQC = False       # Remove all particles that are probably not rain?
use_hailonlyQC = False      # Remove all particles that are probably not hail?
use_graupelonlyQC = False   # Remove all particles that are probably not graupel?

use_rainfallspeedQC = True      # Remove all fall speeds +/- a percentage of the rain fall speed relation
                                # Percentage is set by falltol
falltol = 0.6
maskhigh = True
masklow = True

maskhighdiam = False     # True to mask out high diameter particles (given by threshold below), False to keep them
masklowdiam  = False    # True to mask out low diameter particles (given by threshold below), False to keep them

highdiamthresh = 9.0
lowdiamthresh = 1.0

plot_QC = False        # True to plot fallspeed vs. diameter plot with QC information for each DSD
plot_splashingQC = False
plot_marginQC = False
plot_strongwindQC = False
plot_rainonlyQC = False
plot_rainfallspeedQC = False

def DDMtoDD(DDM,hem):
    """Converts from 'Degrees + Decimal Minutes' format to Decimal Degrees"""
    
    degrees = N.floor(DDM)
    DM = (DDM-degrees)*100.
    
    if(hem == 'N' or hem == 'E'):
        sign = 1.
    else:
        sign = -1.
    
    return sign*(degrees+DM/60.)

def assignfallspeed(d):
    """Assigns a fall speed for a range of diameters based on code
       from David Dowell (originally from Terry Schuur).  It appears that
       the formulas originate from Atlas et al. (1973), but this took a bit of sleuthing!"""
    
    # Note, this appears to be valid at sea level.  For higher altitudes, a fall speed correction
    # should probably be applied based on Foote and duToit (1969): v = v0*(rho0/rho)^(0.4)
    # where rho0 = 1.204 kg/m^3 -- that corresponding to a T of 20 C and pressure of 1013 mb.
    
    v = N.where(d < 3.0, 3.78*d**0.67,9.65-10.3*N.exp(-0.6*d))
    v = N.where(d < 0.0, 0.0, v)
        
    return v

# Create mask for splashing drops
bottom = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
top=[0,1,4,7,9,11,12,13,14,14,15,16,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
splashingmask = [[True if bottom[j] <= i <= top[j] else False for i in range(32)] for j in range(32)]
splashingmask = N.array(splashingmask).T

# Create mask for margin falls
bottom = [0,7,13,16,19,20,21,22,23,24,25,25,26,26,27,27,28,28,28,29,29,29,29,0,0,0,0,0,0,0,0,0]
top = N.zeros((32),dtype='int')
top[:] = 31
top[23:32] = 0
marginmask = [[True if bottom[j] <= i <= top[j] else False for i in range(32)] for j in range(32)]
marginmask = N.array(marginmask).T

# Create mask for non-raindrops
bottom=[0,1,4,7,9,11,12,13,14,14,15,16,16,19,19,20,20,21,21,21,23,24,24,0,0,0,0,0,0,0,0,0]
top=[0,8,14,17,20,21,22,23,24,25,26,26,27,27,28,28,29,29,29,30,30,30,30,0,0,0,0,0,0,0,0,0]
rainonlymask = [[False if bottom[j] <= i <= top[j] else True for i in range(32)] for j in range(32)]
rainonlymask = N.array(rainonlymask).T

# Create mask for all particles with fall speeds outside of fractional tolerance
rainvd = assignfallspeed(avg_diameter)
X,Y = N.meshgrid(rainvd,fall_bins)

if(maskhigh and masklow):
    fallspeedmask = N.where(N.abs((X-Y)/X) < falltol,False,True)
elif(masklow):     # Mask out speeds lower than tolerance
    fallspeedmask = N.where((Y-X)/X < -falltol,True,False)
elif(maskhigh):               # Mask out speeds higher than tolerance
    fallspeedmask = N.where((Y-X)/X > falltol,True,False)
else:
    fallspeedmask = None

# Create mask for strong wind conditions
strongwindmask = N.zeros((32,32),dtype=bool)
strongwindmask[0:11,20:32] = True

def truncatedspectrumQC(countsMatrix):
    """Masks out bad records where tokens have been set to -999 because the record
       was truncated"""
    countsMatrix = ma.masked_array(countsMatrix,mask=N.where(countsMatrix == -999,True,False))
        
    return countsMatrix

def strongwindQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine.  Removes drops affected by strong winds.
       Specifically, removes the entire 1-min DSD whenever there are large drops (> 3 mm) that have
       a low fall velocity (< 1 m/s)."""
    
    numtimes = N.size(countsMatrix,axis=0)
    flaggedtimes = []
    
    # Flag times that contain wind contamination
    for t in range(numtimes):
        baddrops = N.sum(countsMatrix[t,0:11,20:32])
        bigdrops = N.sum(countsMatrix[t,:,23:32])
        totaldrops = N.sum(countsMatrix[t,:])
        if(baddrops > 0.02*totaldrops):  # Try relaxing criterion to allow up to 2% of drops
                                         # to be in mask area
            print "Severe Wind contamination, masking entire PSD!"
            countsMatrix[t,:] = -999.
            flaggedtimes.append(True)
        else: # Let the PSD through QC, but mask the offending drops
            print "Wind contamination!"
            countsMatrix[t,0:11,20:32] = -999.
            flaggedtimes.append(False)
    
    countsMatrix = ma.masked_array(countsMatrix,mask=N.where(countsMatrix == -999.,True,False))
    
    return countsMatrix,flaggedtimes

def splashingQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from splashing"""
    
    numtimes = N.size(countsMatrix,axis=0)
    
    #Remove drops that likely result from splashing (use mask to index the count array)
    
    for t in range(numtimes):
        #countsMatrix[t,mask] = 0.0
        countsMatrix[t,:] = ma.masked_array(countsMatrix[t,:],mask=splashingmask)
        
    return countsMatrix

def marginQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes drops from margin falls"""
    
    numtimes = N.size(countsMatrix,axis=0)
    
    #Remove drops that likely result from margin falls
    
    for t in range(numtimes):
        countsMatrix[t,marginmask] = 0.0
    
    return countsMatrix

def rainonlyQC(countsMatrix):
    """Based on Katja Friedrich's IDL QC subroutine. Removes particles that are probably
       not raindrops"""
    
    numtimes = N.size(countsMatrix,axis=0)
    masktimes = N.zeros((numtimes,32,32),dtype=bool)

    # Remove particles that are probably not rain
    for t in range(numtimes):
        masktimes[t,:] = rainonlymask
    
    countsMatrix = ma.masked_array(countsMatrix,mask=masktimes)
    
    return countsMatrix
    
def rainfallspeedQC(countsMatrix,rainvd,falltol,maskhigh,masklow):
    """Removes all drops fall speeds +/- tolerance relative to rain fall speed relation."""
    
    numtimes = N.size(countsMatrix,axis=0)
    masktimes = N.zeros((numtimes,32,32),dtype=bool)
    
    for t in range(numtimes):
        masktimes[t,:] = fallspeedmask

    countsMatrix = ma.masked_array(countsMatrix,mask=masktimes)
        
    return countsMatrix

def maskhighdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""
    
    numtimes = N.size(countsMatrix,axis=0)
    
    diamindex = N.where(avg_diameter > highdiamthresh)[0][0]
    mask = N.zeros((numtimes,32,32))
    
    mask[:,:,diamindex:] = 1
    countsMatrix = ma.masked_array(countsMatrix,mask=mask)
    
    return countsMatrix
    
def masklowdiamQC(countsMatrix):
    """Mask out particles above a certain diameter"""
    
    numtimes = N.size(countsMatrix,axis=0)
    
    diamindex = N.where(avg_diameter > lowdiamthresh)[0][0]
    mask = N.zeros((numtimes,32,32))
    
    mask[:,:,:diamindex] = 1
    countsMatrix = ma.masked_array(countsMatrix,mask=mask)
    
    return countsMatrix
    
def correctPIPS(serialnum,infile,outfile):
    """Corrects offset Parsivel strings in a PIPS data file"""
    disfile = open(infile,'r')
    
    lines_out = []
    parsivel_string_out_list = []
    parsivel_string_out = ''
    truncated = False
    first = True
    l = 0
    for line in disfile:
        line = line.rstrip('\r\n')
        tokens = line.strip().split(',')
        header_info = ",".join(tokens[:26])
        parsivel_string = tokens[26]
        parsivel_tokens = parsivel_string.strip().split(';')
        if(len(parsivel_tokens) > 1): 
            if(truncated): # Previous record was truncated
                # Look for the parsivel serial number somewhere in the next record
                # and find the index if it exists
                #print parsivel_tokens
                try:
                    sindex = [i for i,x in enumerate(parsivel_tokens) if serialnum in x][0]
                    #sindex = parsivel_tokens.index(serialnum)
                    #print sindex,parsivel_tokens[sindex]
                    # Concatenate portion of string before serial number onto end of previous record
                    parsivel_string_out = parsivel_string_out+";".join(parsivel_tokens[:sindex])
                    #print parsivel_string_out
                    lines_out.append(line_out+','+parsivel_string_out)
                    parsivel_string_out_list.append(parsivel_string_out)
                    parsivel_string_out = ";".join(parsivel_tokens[sindex:]).lstrip()
                    line_out = header_info
                    truncated = False
                except:
                    print "Something went wrong!"
            elif(not first):
                # Should be able to just rip through the rest of the records and piece them together
                sindex = [i for i,x in enumerate(parsivel_tokens) if serialnum in x][0]
                parsivel_string_out = parsivel_string_out+";".join(parsivel_tokens[:sindex])
                parsivel_string_out_list.append(parsivel_string_out)
                lines_out.append(line_out+','+parsivel_string_out)
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
    sorted_lines_out = sorted(lines_out,key=lambda record: int(record.strip().split(',')[1]))
    #sorted_lines_out = lines_out
    
    outdisfile = open(outfile,'w')
    for line in sorted_lines_out:
        outdisfile.write(line+'\n')

def readPIPS(filename,fixGPS=True,basicqc=False,rainfallqc=False,rainonlyqc=False,strongwindqc=False,DSD_interval=10.0):
    """Reads data from Purdue-OU PIPS"""
    
    pdatetimes=[]
    intensities=[]
    preciptots=[]
    reflectivities=[]
    sampleintervals=[]
    pcounts=[]
    pcounts2=[]
    sensortemps=[]
    amplitudes=[]
    voltages=[]
    pvoltages=[]
    
    countsMatrix=[]     # Will contain the number of drops in the diameter vs. fall speed matrix for each time
                        # i.e. a (numtimes,32,32) array
                
    concentrations=[]   # Will contain computed number concentrations for each size bin for each time
                        # i.e. a (numtimes,32) array
                        
    onedrop_concentrations=[]
                  
    dates=[]
    times=[]
    datetimes=[]             
    windspds=[]
    winddirrels=[]
    winddirabss=[]
    winddiags=[]
    fasttemps=[]
    slowtemps=[]
    dewpoints=[]
    RHs_derived=[]
    RHs=[]
    pressures=[]
    compass_dirs=[]
    GPS_lats=[]
    GPS_lons=[]
    GPS_stats=[]
    GPS_alts=[]
    
    disfile = open(filename,'r')
    
    firstgoodGPS = False
    
    for line in disfile:
        tokens = line.strip().split(',')
        timestamp = tokens[0]
        timestring = timestamp.strip().split()
        date = timestring[0] # .strip('-')
        time = timestring[1] # .strip(':')
        
        #2016-03-31 22:19:02
        
        #Construct datetime object
        year = N.int(date[:4])
        month = N.int(date[5:7])
        day = N.int(date[8:])
        hour = N.int(time[:2])
        min = N.int(time[3:5])
        sec = N.int(time[6:])

        datetimelogger = datetime(year,month,day,hour,min,sec)
        
        recordnum = N.int(tokens[1])
        voltage = N.float(tokens[2])
        paneltemp = N.float(tokens[3])
        winddirrel = N.float(tokens[4])
        windspd = N.float(tokens[5])
        winddiag = N.float(tokens[6])
        fasttemp = N.float(tokens[7])
        slowtemp = N.float(tokens[8])
        RH = N.float(tokens[9])
        pressure = N.float(tokens[10])
        compass_dir = N.float(tokens[11])
        GPS_time = tokens[12]
        GPS_status = tokens[13]
        GPS_lat = N.float(tokens[14])
        GPS_lat_hem = tokens[15]
        GPS_lat = DDMtoDD(GPS_lat,GPS_lat_hem)
        GPS_lon = N.float(tokens[16])
        GPS_lon_hem = tokens[17]
        GPS_lon = DDMtoDD(GPS_lon,GPS_lon_hem)
        GPS_spd = N.float(tokens[18])
        GPS_dir = N.float(tokens[19])
        GPS_date = tokens[20]
        GPS_magvar = N.float(tokens[21])
        GPS_alt = N.float(tokens[22])
        
        winddirabs = N.float(tokens[23])
        try:
            dewpoint = N.float(tokens[24])
        except:
            dewpoint = N.nan
        try:
            RH_derived = N.float(tokens[25])
        except:
            RH_derived = N.nan
                 
        # Find the first good GPS time and date and use that
        # to construct the time offset for the data logger
        if(not N.isnan(GPS_alt) and not firstgoodGPS and GPS_status == 'A'):
            firstgoodGPS = True
            print GPS_date,GPS_time
            print date,time
            
            print year,month,day,hour,min,sec
            
            #Construct datetime object
            gyear = N.int('20'+GPS_date[4:])
            gmonth = N.int(GPS_date[2:4])
            gday = N.int(GPS_date[:2])
            ghour = N.int(GPS_time[:2])
            gmin = N.int(GPS_time[2:4])
            gsec = N.int(GPS_time[4:])

            datetimeGPS = datetime(gyear,gmonth,gday,ghour,gmin,gsec)
            GPS_offset = datetimeGPS-datetimelogger
            print datetimeGPS,datetimelogger
            print "GPS Offset",GPS_offset
    
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
            #print timestring
            precipintensity = N.float(parsivel_tokens[1])
            precipaccum = N.float(parsivel_tokens[2])
            parsivel_dBZ = N.float(parsivel_tokens[3])
            sample_interval = N.float(parsivel_tokens[4])
            signal_amplitude = N.float(parsivel_tokens[5])
            pcount = N.int(parsivel_tokens[6])
            sensor_temp = N.float(parsivel_tokens[7])
            pvoltage = N.float(parsivel_tokens[8])
            sensor_time = parsivel_tokens[9]
            sensor_date = parsivel_tokens[10]
            try:
                spectrum = [float(x) if x!= '' else 0 for x in parsivel_tokens[11:]]
            except:
                spectrum = [-999 for i in xrange(1025)]
            #print "spectrum length = ",len(spectrum)
            if(len(spectrum) < 1024):
                print "Problem with Parsivel spectrum.  Flagging as bad!"
                print "Time: ",timestring
#                 spectrum = [N.nan for i in xrange(1025)]
            else:
                if(len(spectrum) == 1025):
                    spectrum = spectrum[:-1] # Strip off bogus last value
        
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
            spectrum = N.array(spectrum,dtype='int')
            if(spectrum.size == 1024):
                spectrum = spectrum.reshape((32,32))
            else:
                spectrum = -999*N.ones((32,32),dtype='int')
                print spectrum.size
            # Append spectrum (corrected or not) to spectrum time list

            countsMatrix.append(spectrum)
        
    # Recast countsMatrix as numpy array
    
    countsMatrix = N.dstack(countsMatrix)
    countsMatrix = N.rollaxis(countsMatrix,2,0)
    
    # Perform Katja's QC routines if desired (should not be used in combination with my methods above,
    # most are redundant anyway).

    X,Y = N.meshgrid(avg_diameter,fall_bins)
    flaggedtimes=N.zeros(len(pdatetimes),dtype=bool)
    splashmask=N.zeros_like(countsMatrix)
    marginmask=N.zeros_like(countsMatrix)
    
    countsMatrix = truncatedspectrumQC(countsMatrix)
    
    if(use_strongwindQC or basicqc or strongwindqc):
        countsMatrix,flaggedtimes = strongwindQC(countsMatrix)
    
    if(use_splashingQC or basicqc):
        countsMatrix = splashingQC(countsMatrix)
    
    if(use_marginQC or basicqc):
        countsMatrix = marginQC(countsMatrix)

    if(use_rainfallspeedQC or rainfallqc):
        countsMatrix = rainfallspeedQC(countsMatrix,rainvd,falltol,maskhigh,masklow)

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
    
    counts_1drop = N.ones_like(avg_diameter)
    
    #print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t,time in enumerate(pdatetimes):
        spectrum = countsMatrix[t,:]
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy,vspectrum = N.meshgrid(avg_diameter,fall_bins)
            dspectrum = spectrum
        else:
            vspectrum = rainvd
            dspectrum = spectrum.sum(axis=0)   # Sum up particles for each diameter bin
            
            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]
    
        # Now compute the number concentration using the assumed fall speeds, sensor area, and sampling interval
        # Units are #/m^3/mm
#         if time == '171320':
#             print spectrum.size,spectrum,countsMatrix.data[t,:],spectrum.data
        if(spectrum.size == 1024 and not flaggedtimes[t]):
            concentration = dspectrum/(vspectrum*sampleintervals[t]*eff_sensor_area*(max_diameter-min_diameter))
            onedrop_concentration = counts_1drop/(rainvd*sampleintervals[t]*eff_sensor_area*(max_diameter-min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
            #print "concentration.shape"
            #print concentration.shape
        elif(not flaggedtimes[t]):
            concentration = N.zeros_like(avg_diameter)
        else:
            concentration = -999.*N.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999.,concentration)
        
        # Throw out particles above and below a certain diameter if desired
        
#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)
        
        #print "Number of particles counted vs. summed number: "
        #print N.int(line[7]),N.sum(dspectrum)
        
        #print concentration
        
        concentrations.append(concentration)
        onedrop_concentrations.append(onedrop_concentration)
        
    concentrations = ma.array(concentrations)
    onedrop_concentrations = ma.array(onedrop_concentrations)
    #print "concentrations: ",concentrations
    pcounts = N.array(pcounts)
    #pcounts2 = ma.array(pcounts2)
    
    # Correct the logger time and date using the GPS time
    
    datetimes_corrected = []
    pdatetimes_corrected = []
    for datetimelogger in datetimes:
        datetimes_corrected.append(datetimelogger+GPS_offset)
    for pdatetimelogger in pdatetimes:
        pdatetimes_corrected.append(pdatetimelogger+GPS_offset)

    # Average and thin the DSD data with the desired interval
    DSD_index_interval = int(DSD_interval/10.0)    
    print "Requested DSD interval: {:.1f}. Actual DSD interval: {:.1f}".format(DSD_interval,DSD_index_interval*10.0)
    DSD_interval = DSD_index_interval*10.0
    # We need to find the offset corresponding to the starting second and then generate the frequency string
    sec_offset = pdatetimes_corrected[0].second
    intervalstr = '{:d}S'.format(int(DSD_interval))
    
    concentrations_df = pd.DataFrame(data=concentrations,index=pdatetimes_corrected,columns=avg_diameter)
    onedrop_concentrations_df = pd.DataFrame(data=onedrop_concentrations,index=pdatetimes_corrected,columns=avg_diameter)
    intensities_df = pd.Series(intensities,index=pdatetimes_corrected)
    preciptots_df = pd.Series(preciptots,index=pdatetimes_corrected)
    reflectivities_df = pd.Series(reflectivities,index=pdatetimes_corrected)
    pcounts_df = pd.Series(pcounts,index=pdatetimes_corrected)
    pcounts2_df = pd.Series(pcounts2,index=pdatetimes_corrected) # STOPPED HERE!

    
    if(DSD_interval > 10.0):
        # Create a dataframe of the concentrations in each diameter bin and then resample at the new interval,
        # filling in missing values with zeros
    
        concentrations_df = concentrations_df.resample(intervalstr,label='right',closed='right',base=sec_offset).mean().fillna(0)
        onedrop_concentrations_df = onedrop_concentrations_df.resample(intervalstr,label='right',closed='right',base=sec_offset).mean().fillna(0)
        intensities_df = intensities_df.resample(intervalstr,label='right',closed='right',base=sec_offset).mean().fillna(0)
        preciptots_df = preciptots_df.resample(intervalstr,label='right',closed='right',base=sec_offset).sum()
        reflectivities_df = reflectivities_df.resample(intervalstr,label='right',closed='right',base=sec_offset).mean().fillna(0)
        pcounts_df = pcounts_df.resample(intervalstr,label='right',closed='right',base=sec_offset).sum()
        pcounts2_df = pcounts2_df.resample(intervalstr,label='right',closed='right',base=sec_offset).sum()

    # Pandas apparently (gotcha!) converts missing values to NaN when extracting the numpy array representation using .values
    # Since I want the original masked array functionality for now for further computations, I need to remask the array here.
    # Otherwise, the NaN's propagate in further computations...
    # In the future, another solution that uses Pandas more natively should be pursued, but this will work for now
    concentrations = concentrations_df.values
    mask = concentrations_df.isnull()
    concentrations = ma.array(concentrations,mask=mask)
    onedrop_concentrations = onedrop_concentrations_df.values
    mask = onedrop_concentrations_df.isnull()
    onedrop_concentrations = ma.array(onedrop_concentrations,mask=mask)
    #print len(pdatetimes_corrected),pdatetimes_corrected
    # Argh, have to convert back to datetime objects.  This one from http://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64
    pdatetimes_corrected = concentrations_df.index.to_pydatetime()
    DSD_index = concentrations_df.index
    #print len(pdatetimes_corrected),pdatetimes_corrected
    onedrop_concentrations = 10.0*onedrop_concentrations/DSD_interval # Concentration in each bin assuming only one drop over the new interval
    intensities = intensities_df.values
    preciptots = preciptots_df.values
    reflectivities = reflectivities_df.values
    pcounts = pcounts_df.values
    pcounts2 = pcounts2_df.values
        
    return datetimes_corrected,pdatetimes_corrected,flaggedtimes,intensities,preciptots,reflectivities,pcounts,pcounts2, \
            sensortemps,concentrations,onedrop_concentrations,countsMatrix,windspds,winddirrels,winddirabss, \
            winddiags,fasttemps,slowtemps,dewpoints,RHs_derived,RHs,pressures,compass_dirs,    \
            GPS_lats,GPS_lons,GPS_stats,GPS_alts,voltages,DSD_interval,intervalstr,DSD_index

def readPIPSloc(filename):
    """Reads the location of the PIPS from the data file.  Averages all valid GPS lat/lons"""
    
    GPS_lats=[]
    GPS_lons=[]
    GPS_stats=[]
    GPS_alts=[]
    
    disfile = open(filename,'r')
        
    for line in disfile:
        tokens = line.strip().split(',')
        
        GPS_status = tokens[13]
        GPS_lat = N.float(tokens[14])
        GPS_lat_hem = tokens[15]
        GPS_lat = DDMtoDD(GPS_lat,GPS_lat_hem)
        GPS_lon = N.float(tokens[16])
        GPS_lon_hem = tokens[17]
        GPS_lon = DDMtoDD(GPS_lon,GPS_lon_hem)
        GPS_alt = N.float(tokens[22])
                     
        GPS_lats.append(GPS_lat)
        GPS_lons.append(GPS_lon)
        GPS_stats.append(GPS_status)
        GPS_alts.append(GPS_alt)
    
    # Find the disdrometer location by averaging the valid GPS lats,lons, and alts
    GPSnonvalid = (N.array(GPS_stats) != 'A')
    GPS_lats_masked = ma.masked_where(GPSnonvalid,N.array(GPS_lats))
    GPS_lons_masked = ma.masked_where(GPSnonvalid,N.array(GPS_lons))
    GPS_alts_masked = ma.masked_where(GPSnonvalid,N.array(GPS_alts))
        
    lat = GPS_lats_masked.mean()
    lon = GPS_lons_masked.mean()
    alt = GPS_alts_masked.mean()
    
    # There's a bug in numpy.ma that sometimes causes operations such as mean() to return
    # a 0d array instead of a scalar.  To adjust for this, explicitly cast them as scalars
    # here.
    
    lat = N.asscalar(lat)
    lon = N.asscalar(lon)
    alt = N.asscalar(alt)
    
    dloc = (lat,lon,alt)
        
    return GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc

def PIPS2TTU(filename,fixGPS=True):
    """Reads conventional data from a PIPS text file and outputs it in the format used by TTU's
       StickNets"""
    
    datetimes_corrected,pdatetimes_corrected,flaggedtimes,intensities,preciptots,reflectivities,pcounts,pcounts2, \
            sensortemps,concentrations,gamma_concentrations,countsMatrix,windspds,winddirrels,winddirabss, \
            winddiags,fasttemps,slowtemps,dewpoints,RHs_derived,RHs,pressures,compass_dirs,    \
            GPS_lats,GPS_lons,GPS_stats,GPS_alts,voltages = readPIPS(filename)
    
    GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc = readPIPSloc(filename)
    
    splitpath = os.path.split(filename)
    inputfilename = splitpath[1]
    outputfiledir = splitpath[0]
    outputfilename = inputfilename[:-4]+'_TTU.txt'
    print outputfilename
    outputfilepath = outputfiledir+outputfilename
    
    outputfile = open(outputfilepath,'w')
    
    probeID = os.path.split(filename)[1][5:7]
    date = str(datetimes_corrected[0].year)+str(datetimes_corrected[0].month)+str(datetimes_corrected[0].day)
    time = str(datetimes_corrected[0].hour)+':'+str(datetimes_corrected[0].minute)+':'+str(datetimes_corrected[0].second)
    lat = dloc[0]
    lon = dloc[1]
    alt = dloc[2]
    heading = compass_dirs[0]
    dummy = 3
    
    header=[probeID,date,time,lat,lon,alt,heading,dummy]
    header = ','.join(map(str, header))
    print header
    
    return date

def readtmatrix(filename):
    """Reads a scattering amplitude lookup table created using the ARPS tmatrix program"""
    
    data = N.loadtxt(filename,skiprows=1)
    d = data[:,0]
    far_b = data[:,1]+1j*data[:,2]
    fbr_b = data[:,3]+1j*data[:,4]
    far_f = data[:,5]+1j*data[:,6]
    fbr_f = data[:,7]+1j*data[:,8]
        
    return d,far_b,fbr_b,far_f,fbr_f
    
def calbackscatterrain(far_b,fbr_b,far_f,fbr_f):
    """Calculates rain backscattering amplitudes.
       Based on similar code in dualpara.f90"""
    
    fa2 = (N.abs(far_b))**2.
    fb2 = (N.abs(fbr_b))**2.
    fab = far_b*N.conjugate(fbr_b)
    fba = fbr_b*N.conjugate(far_b)
    far = N.real(far_f-fbr_f)
    
    return fa2,fb2,fab,fba,far
    
def calpolrain(wavelength,filename,Nd,intv):
    """Given backscattering amplitudes and a discrete distribution N(D) (m^-4), compute
       polarimetric variables for each bin."""
    
    d,far_b,fbr_b,far_f,fbr_f = readtmatrix(filename)
    fa2,fb2,fab,fba,far = calbackscatterrain(far_b,fbr_b,far_f,fbr_f)
    
    # There may be more bins in the given Nd than are read in from the file.
    # This is because the Nd contains bins above the maximum size of rain (~9 mm).  
    # Truncate the first dimension of Nd to account for this.
    # Also transpose Nd to allow for numpy broadcasting
    Nd = Nd[:N.size(fa2),:].T
    intv = intv[:N.size(fa2)]
        
    lamda = wavelength*10. # Get radar wavelength in mm
    Kw2 = 0.93 # Dielectric factor for water
    
    sar_h = fa2*Nd*intv
    sar_v = fb2*Nd*intv
    sar_hv = fab*Nd*intv
    fsar = far*Nd*intv 

    Zh = 4.*lamda**4./(N.pi**4.*Kw2)*N.sum(sar_h,axis=1)
    Zv = 4.*lamda**4./(N.pi**4.*Kw2)*N.sum(sar_v,axis=1)
    Zhv= 4.*lamda**4./(N.pi**4.*Kw2)*N.abs(N.sum(sar_hv,axis=1))
    Kdp= 180.*lamda/N.pi*N.sum(fsar,axis=1)*1.e-3
    dBZ = 10.*N.log10(Zh)
    temp = Zh/Zv
    ZDR = 10.*N.log10(N.maximum(1.0,temp))
    temp = Zh*Zv
    rhv = N.where(Zh != Zv ,Zhv/(N.sqrt(temp)),0.0) # Added by Jessie (was temp > 0).  Find out why...
    #N.savetxt('temp.txt', temp)


    return Zh,Zv,Zhv,dBZ,ZDR,Kdp,rhv
    

def calc_DSD(min_size,avg_size,max_size,bin_width,Nc_bin,logNc_bin,rho,qrQC,qr_thresh):
    """Fits exponential and gamma DSDs to disdrometer data and returns several DSD related quantities."""
    
    
    min_size = N.array(min_size)
    avg_size = N.array(avg_size)
    max_size = N.array(max_size)
    bin_width = N.array(bin_width)
    
    # First calculate the required moment estimators
    M0=[]   # Not used in moment estimator, but used for mean diameter calculation
    M1=[]   # Ditto
    M2=[]
    M3=[]
    M4=[]
    M6=[]
    M7=[]   # Ditto
    D_med_disd=[]    # Median volume diameter

    try:
        numtimes = N.size(Nc_bin,axis=1)
    except:
        numtimes = 1
        Nc_bin = Nc_bin[:,N.newaxis]
    
    
    
    for t in range(numtimes):
    
        temp_M0 = (1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M1 = ((avg_size[:]/1000.))*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M2 = ((avg_size[:]/1000.)**2.)*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M3 = ((avg_size[:]/1000.)**3.)*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M4 = ((avg_size[:]/1000.)**4.)*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M6 = ((avg_size[:]/1000.)**6.)*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
        temp_M7 = ((avg_size[:]/1000.)**7.)*(1000.*Nc_bin[:,t])*bin_width[:]/1000.
    
        # Before summing up the 3rd moment, use the bin information to find the (approximate) median volume diameter
        
        #print "temp_M3",temp_M3,temp_M3.shape
        #print "N.array(temp_M3)",N.array(temp_M3)
        temp_M3_cumsum = N.cumsum(temp_M3) # Cumulative sum of M3 with increasing bin size
        temp_M3_sum = N.sum(temp_M3)
        pro = temp_M3/N.sum(temp_M3_sum) # Proportion of M3 in each bin
        pro_cumsum = temp_M3_cumsum/N.sum(temp_M3_sum) # Cumulative proportion of M3 in each bin
    
        #print "temp_M3_cumsum",temp_M3_cumsum
        
        # Compute median volume diameter using a linear interpolation within the "mass-midpoint" bin.
        # Source: http://www.dropletmeasurement.com/PADS_Help/MVD_(um).htm (originally from FAA Electronic Aircraft Icing Handbook)
        try:
            medindex = N.where(pro_cumsum > 0.5)[0][0] # Find index of bin where cumulative sum exceeds 1/2 of total
            b1 = min_size[medindex] # Lower boundary of mass-midpoint bin
            b2 = max_size[medindex] # Upper boundary of mass-midpoint bin
            if(medindex == 0):
                pro_cumsum_medindexm1 = 0.0
            else:
                pro_cumsum_medindexm1 = pro_cumsum[medindex-1]
            temp_D_med = b1 + ((0.5 - pro_cumsum_medindexm1)/pro[medindex])*(b2-b1) # Linearly-interpolated (within mass-midpoint bin) D0
        except:
            medindex = 0
            temp_D_med = N.nan
        temp_M0_2 = ma.sum(temp_M0)
        temp_M0 = N.sum(temp_M0)
        temp_M1 = N.sum(temp_M1)
        temp_M2 = N.sum(temp_M2)
        temp_M4 = N.sum(temp_M4)
        temp_M3 = N.sum(temp_M3)
        temp_M6 = N.sum(temp_M6)
        temp_M7 = N.sum(temp_M7)
    
        M0.append(temp_M0)
        M1.append(temp_M1) 
        M2.append(temp_M2)
        M3.append(temp_M3)
        M4.append(temp_M4)
        M6.append(temp_M6)
        M7.append(temp_M7)
        D_med_disd.append(temp_D_med)
            
    M0 = ma.array(M0,dtype=N.float64)
    M1 = ma.array(M1,dtype=N.float64)
    M2 = ma.array(M2,dtype=N.float64)
    M3 = ma.array(M3,dtype=N.float64)
    M4 = ma.array(M4,dtype=N.float64)
    M6 = ma.array(M6,dtype=N.float64)
    M7 = ma.array(M7,dtype=N.float64)
    D_med_disd = N.array(D_med_disd)
        
    # --- Compute various mean diameters directly from measured discrete distribution ---
    cmr = (N.pi/6.)*1000.                                     # Constant in mass-diameter relation for rain
    LWC_disd = cmr*1000.0*M3
    #print LWC_disd.shape
    #print rho.shape
    QR_disd = LWC_disd/rho
    if(qrQC):
        qrmask1D = N.where(QR_disd > qr_thresh,True,False)
        qr2D = QR_disd.reshape(1,numtimes).repeat(32,0)
        qrmask2D = N.where(qr2D > qr_thresh,True,False)
        
        #print qrmask2D.shape
        
        # Mask out all the needed arrays, then all derived arrays below will be masked in the same
        # way
        
        D_med_disd = ma.masked_array(D_med_disd,mask=qrmask1D)
        Nc_bin = ma.masked_array(Nc_bin,mask=qrmask2D)
        logNc_bin = ma.masked_array(logNc_bin,mask=qrmask2D)
        M0 = ma.masked_array(M0,mask=qrmask1D)
        M1 = ma.masked_array(M1,mask=qrmask1D)
        M2 = ma.masked_array(M2,mask=qrmask1D)
        M3 = ma.masked_array(M3,mask=qrmask1D)
        M4 = ma.masked_array(M4,mask=qrmask1D)
        M6 = ma.masked_array(M6,mask=qrmask1D)
        M7 = ma.masked_array(M7,mask=qrmask1D)
        rho = ma.masked_array(rho,mask=qrmask1D)
        QR_disd = ma.masked_array(QR_disd,mask=qrmask1D)
    
    # Compute reflectivity from the 6th moment
    refl_disd=10.0*ma.log10(1e18*M6)
    #print 'reflectivity (disdrometer) = ',refl_disd
    
    # Compute the mass-weighted mean diameter (M(4)/M(3) )
    
    D_m_disd = (M4/M3)*1000.0
    
    # Compute the mean-volume (or mean-mass if rho=constant) diameter ((M(3)/M(0))^1./3.)
    
    D_mv_disd = ((M3/M0)**(1./3.))*1000.0
    
    # Compute reflectivity-weighted mean diameter (M(7)/M(6))
    
    D_ref_disd = (M7/M6)*1000.0
        
    # --- Exponential Distribution Fit by Method of Moments (Zhang et al. 2008) ---
    
    # Now that we have the moment estimates, we can calculate N0 and lamda (eqns. 3 and 4 in Zhang et al. 2008)
    
    gamma3 = special.gamma(3.)
    gamma4 = special.gamma(4.)
    gamma5 = special.gamma(5.)
    gamma7 = special.gamma(7.)
    
    # Uncomment if you want estimates based on M2 and M4
    #lamda = N.where(M4 == 0.0, 0.0, ((M2*gamma5)/(M4*gamma3))**(1./2.))
    #N0 = (M2*lamda**3.)/gamma3
        
    # Uncomment if you want estimates based on M3 and M6
    lamda_exp = N.where(M6 == 0.0, 0.0, ((M3*gamma7)/(M6*gamma4))**(1./3.)) # Should be masked
    N0_exp = (M3*lamda_exp**4.)/gamma4                                      # Should be masked
    mu_exp = 0.0
    
    #print 'lamda (exp) = ',lamda_exp
    #print 'N0 (exp) = ',N0_exp
    
    # Now create synthetic bins to plot the derived exponential DSD
    # Each index holds the value of the midpoint diameter of the bin
    # in mm
    
    synthbins = N.linspace(0.1,16.0,num=160)
    
    N_expDSD=[]
    #refl_expDSD=[]
    # Construct the derived DSD
    for t in range(numtimes):
        temp_N_expDSD = N0_exp[t]*N.exp(-lamda_exp[t]*synthbins/1000.0)
        N_expDSD.append(temp_N_expDSD)
    
    N_expDSD = ma.array(N_expDSD)
    
    # --- Gamma Distribution Fit by Method of Moments (Tokay and Short 1996) ---
    
    # Now do the same as above for a gamma distribution fit
    
    # Calculate G
    
    G = N.where((M3 == 0.0) | (M6 == 0.0), 0.0, (M4**3.)/((M3**2.)*M6))
    
    # Calculate mu (shape parameter), lamda, and N0
    
    mu_gam = N.where(G == 1.0, 0.0, (11.*G-8.+(G*(G+8.))**(1./2.))/(2.*(1.-G)))
    mu_gam = N.where(mu_gam <= -4.,-3.99,mu_gam)
    mu_gam = N.where(mu_gam > 40.,40.,mu_gam)
    mu_gam = ma.masked_where(M3 is ma.masked,mu_gam)     # Mask where the moment is masked
    lamda_gam = N.where(M4 == 0.0, 0.0, (M3*(mu_gam+4.))/M4)
    N0_gam = (M3*lamda_gam**(mu_gam+4.))/(special.gamma(mu_gam+4.))
        
    #print 'mu (gam) = ',mu_gam
    #print 'lamda (gam) = ',lamda_gam
    #print 'N0 (gam) = ',N0_gam
    
    #print N0_gam.dtype
    
    N_gamDSD=[]
    #refl_gamDSD=[]
    
    for t in range(numtimes):
        #temp_N_gamDSD = N0_gam[t]*((synthbins/1000.0)**mu_gam[t])*N.exp(-lamda_gam[t]*synthbins/1000.0)
        temp_N_gamDSD = N0_gam[t]*((avg_size/1000.0)**mu_gam[t])*N.exp(-lamda_gam[t]*avg_size/1000.)    # Use disdrometer bins
        #temp_Z_gamDSD = ((synthbins/1000.)**6.)*(1000.*temp_N_gamDSD)*0.1/1000.
        N_gamDSD.append(temp_N_gamDSD)
        #refl_gamDSD.append(10.0*N.log10(N.sum(temp_Z_gamDSD)))
    
    N_gamDSD = ma.array(N_gamDSD,dtype=N.float64)
    
    # Quantities based on exponential distribution    
    
    GR1 = special.gamma(1.+mu_exp)
    GR2 = special.gamma(4.+mu_exp)
    
    qr_exp = (cmr/rho)*N0_exp*GR2/lamda_exp**(mu_exp+4.)
    Ntr_exp = N0_exp*GR1/lamda_exp**(mu_exp+1.)
    # Compute reflectivity for DSD fit
    Gr_exp = ((6.+mu_exp)*(5.+mu_exp)*(4.+mu_exp))/((3.+mu_exp)*(2.+mu_exp)*(1.+mu_exp))
    Zr_exp = ((1./cmr)**2.)*Gr_exp*((rho*qr_exp)**2.)/Ntr_exp
    refl_DSD_exp = 10.0*ma.log10(1.e18*Zr_exp)
    #print 'reflectivity (exp DSD) = ',refl_DSD
    #print 'reflectivity (gam DSD) = ',refl_gamDSD
    
    D_med_exp = N.where(lamda_exp == 0., N.nan, (3.67/lamda_exp)*1000.0)     # Median volume diameter for exponential distribution
    D_m_exp = N.where(lamda_exp == 0., N.nan, (4./lamda_exp)*1000.0) # Mass-weighted mean diameter for exp dist.
    
    # Quantities based on gamma distribution    
    
    GR1 = special.gamma(1.+mu_gam)
    GR2 = special.gamma(4.+mu_gam)
    
    qr_gam = (cmr/rho)*N0_gam*GR2/lamda_gam**(mu_gam+4.)
    Ntr_gam = N0_gam*GR1/lamda_gam**(mu_gam+1.)
    # Compute reflectivity for DSD fit
    Gr_gam = ((6.+mu_gam)*(5.+mu_gam)*(4.+mu_gam))/((3.+mu_gam)*(2.+mu_gam)*(1.+mu_gam))
    Zr_gam = ((1./cmr)**2.)*Gr_gam*((rho*qr_gam)**2.)/Ntr_gam
    refl_DSD_gam = 10.0*ma.log10(1.e18*Zr_gam)
    #print 'reflectivity (exp DSD) = ',refl_DSD
    #print 'reflectivity (gam DSD) = ',refl_gamDSD
    
    D_med_gam = N.where(lamda_gam == 0., N.nan, ((3.67+mu_gam)/lamda_gam)*1000.0)    # Median volume diameter for gamma distribution
    D_m_gam = N.where(lamda_gam == 0., N.nan, ((4.+mu_gam)/lamda_gam)*1000.0) # Mass-weighted mean diameter for gam. dist.
                
    # Create several tuples to pack the data, and then return them
    
    exp_DSD = (N_expDSD,N0_exp,lamda_exp,mu_exp,qr_exp,Ntr_exp,refl_DSD_exp,D_med_exp,D_m_exp)
    gam_DSD = (N_gamDSD,N0_gam,lamda_gam,mu_gam,qr_gam,Ntr_gam,refl_DSD_gam,D_med_gam,D_m_gam)
    dis_DSD = (Nc_bin,logNc_bin,D_med_disd,D_m_disd,D_mv_disd,D_ref_disd,QR_disd,refl_disd)
    
    return synthbins,exp_DSD,gam_DSD,dis_DSD

    
def rad2DD(filename,dlocs):
    """Given an 88D CFRadial file and a list of disdrometer locations (tuples of lat,lon),
       compute the reflectivity at the disdrometer location from the lowest sweep."""
    
    average_gates = False        # True if reflectivity should be averaged in the closest 9 gates
                                # False if just picking the value from the closest gate
    
    Cressman = False             # Perform a Cressman analysis on nearby radar gates to determine
                                # reflectivity value
    
    roi = 1750. # 750.                   # Radius of influence of Cressman analysis in m
    
    
    sweepfile_netcdf = Nio.open_file(filename)
    
    # Grab some needed variables from the netCDF file
    
    numgates = sweepfile_netcdf.dimensions['range']     # Number of gates
    numtimes = sweepfile_netcdf.dimensions['time']     # Number of times
    numsweeps = sweepfile_netcdf.dimensions['sweep']     # Number of sweeps
    sweep_start_ray_index = sweepfile_netcdf.variables['sweep_start_ray_index'][:]
    sweep_end_ray_index = sweepfile_netcdf.variables['sweep_end_ray_index'][:]
    try:
        ray_n_gates = sweepfile_netcdf.variables['ray_n_gates'][:]
        twoDarray = False
    except:
        print "No ray_n_gates in file, assuming 2D arrays (azimuth,range)"
        twoDarray = True
            
    rlat=sweepfile_netcdf.variables['latitude'].get_value()
    rlon=sweepfile_netcdf.variables['longitude'].get_value()
    ralt=sweepfile_netcdf.variables['altitude'].get_value()
    
    #print "Number of gates: ",numgates
    #print "Radar lat,lon,alt",rlat,rlon,ralt
    
    rlat = rlat*deg2rad
    rlon = rlon*deg2rad
    
    elevs=sweepfile_netcdf.variables['elevation'][:]
    
    el = elevs[0]  # Just pick first element for now
    
    #print "Elevation angle ",el
    
    el = el*deg2rad
    
    temp = sweepfile_netcdf.variables['range']    # range to center of each gate
    gatewidth = temp.meters_between_gates[0] # gate spacing for each gate
    range = temp[:]
    range_start = temp[:]-gatewidth/2. # We also want range to start of each gate
    
    #print "Gatewidth ",gatewidth
    
    # Get the Azimuth info
    
    azimuth_rad = sweepfile_netcdf.variables['azimuth'][:]
    azimuth_rad = azimuth_rad[:sweep_end_ray_index[0]] # Just grab first sweep for now
    #print "number of azimuths in sweep ",N.size(azimuth_rad)
    
    # Roll the azimuth dimension around so that it starts at 0 and ends at 360
    try:
        shift = N.where(azimuth_rad < azimuth_rad[0])[0][0]
    except:
        shift = 0
    #print "shift is ",shift
    azimuth_rad = N.roll(azimuth_rad,shift=-shift)
    
    beamwidth = azimuth_rad[1:]-azimuth_rad[0:-1]
    
    azimuth_start_rad = N.zeros(N.size(azimuth_rad)+1)
    # Find azimuth of "start" of each gate (in azimuth) -- approximate
    azimuth_start_rad[1:-1] = azimuth_rad[1:]-0.5*beamwidth[:]
    azimuth_start_rad[0] = azimuth_rad[0]-0.5*beamwidth[0]
    azimuth_start_rad[-1] = azimuth_rad[-1]+0.5*beamwidth[-1]
    
    azimuth_rad = azimuth_rad*deg2rad
    azimuth_start_rad = azimuth_start_rad*deg2rad
    
    # Get the Reflectivity field
    tempdBZ = sweepfile_netcdf.variables['REF']
    dBZ    = sweepfile_netcdf.variables['REF'][:]
    
    # Grab first sweep out of dBZ array
    # Need to check if the array is dimensioned by total number of gates or just number of gates for
    # each sweep
    
    #print ray_n_gates[0]
    #print sweep_end_ray_index[0]
    if(not twoDarray):  # Note, this needs to be tested!
        numpointsinsweep = ray_n_gates[0]*sweep_end_ray_index[0]
        dBZ = dBZ[:numpointsinsweep]
        dBZ = dBZ.reshape((sweep_end_ray_index[0],-1))
    else:
        dBZ = dBZ[:sweep_end_ray_index[0],:]
    
    #print "scale_factor,add_offset",tempdBZ.scale_factor[0],tempdBZ.add_offset[0]
    
    #Unpack values
    dBZ = dBZ*tempdBZ.scale_factor[0]+tempdBZ.add_offset[0] 
    
    # Adjust azimuth axis
    
    dBZ = N.roll(dBZ,shift=-shift,axis=0)
        
    #print azimuth_rad.shape
    #print range.shape
    #print dBZ.shape
        
    theta,rad = N.meshgrid(azimuth_start_rad,range_start)
    theta_c,rad_c = N.meshgrid(azimuth_rad,range)
    xplt_c,yplt_c = oban.xyloc(rad_c,theta_c,el,rlat,rlon,ralt,1)
    dBZ = dBZ.swapaxes(0,1)
    dBZ = ma.masked_invalid(dBZ)
    # Find x,y location of disdrometer
    
    dBZ_D_list = []
    
    for dloc in dlocs:
        Dx,Dy = oban.ll_to_xy(dloc[0]*deg2rad,dloc[1]*deg2rad,rlat,rlon,1)
    
        #print "Dx,Dy",Dx,Dy
        
        # Find closest gate to disdrometer location
        # First find the closest azimuth
        
#         theta_dis = -N.arctan2(Dy,Dx)+N.pi/2.       # I *think* this is correct
#         
#         #print "azimuth of disdrometer: ",theta_dis/deg2rad
#         #Find closest index of azimuth
#         
#         theta_diff = theta_dis-theta_c[0,:]
#         theta_index = N.argmin(N.abs(theta_diff))
#         
#         # Now find index along range; first compute slant range
#         
#         srange_dis = oban.computeslantrange([Dx],[Dy],el)
#         rad_diff = srange_dis-rad_c[:,0]
#         srange_index = N.argmin(N.abs(rad_diff))
#         #print "range of disdrometer: ",srange_dis
        
        # Try another way to get the indices (the above seems to have problems).
        # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why the 
        # above doesn't work in all situations, but it sometimes picks a gate adjacent to the one
        # we want...
        
        if(Cressman):
            #print xplt_c.shape,yplt_c.shape
            dBZ_D = oban.Cresmn(xplt_c,yplt_c,dBZ,Dx,Dy,roi)
            #print dBZ_D,dBZ_D.shape
        else: # First, compute the euclidian distance of the x,y location of the disdrometer to each of the radar gate centers
            distance = N.sqrt(((Dx-xplt_c)**2.+(Dy-yplt_c)**2.))
            
            # Now, find the index of the closest radar gate
            srange_index,theta_index = N.unravel_index(distance.argmin(),distance.shape)
            #print "srange_index,theta_index",srange_index,theta_index
            
            # Finally, grab reflectivity at closest gate to disdrometer
            # Average the reflectivity in the closest gate and surrounding 8 gates if desired
            if(not average_gates):
                dBZ_D = dBZ[srange_index,theta_index]
            else:
                dBZ_D = (1./9.)*(dBZ[srange_index-1,theta_index-1]+dBZ[srange_index,theta_index-1]+
                              dBZ[srange_index-1,theta_index]+dBZ[srange_index,theta_index]+
                              dBZ[srange_index+1,theta_index]+dBZ[srange_index,theta_index+1]+
                              dBZ[srange_index+1,theta_index+1]+dBZ[srange_index+1,theta_index-1]+
                              dBZ[srange_index-1,theta_index+1])
        #print "Reflectivity at disdrometer = ","%.2f"%dBZ_D
        
        dBZ_D_list.append(dBZ_D)
    
    #dBZ_D_arr = N.array(dBZ_D_list)
    
    sweepfile_netcdf.close()
    
    return dBZ_D_list

def rad2DD2(fieldlist,range_start,range,azimuth_start_rad,azimuth_rad,rlat,rlon,ralt,el,dlocs,average_gates=True,Cressman=False,roi=750.,map_proj=1):
    """Another version of rad2DD: assumes radar sweep has been read in and computes values of fields in fieldlist
       at the disdrometer location(s).  Eventually will allow for optional advection/sedimentation correction of radar reflectivity
       to account for height above surface of radar scan.  Returns list of lists of field values at each disdrometer location"""
    
    # First find x,y locations of radar gates
    
    xrad,yrad,xrad_c,yrad_c = radar.sweep2xy(azimuth_start_rad,azimuth_rad,range_start,range,el,rlat,rlon,ralt,map_proj)
    
    field_D_list = []
    dxy_list = []
    
    for dloc in dlocs:
        Dx,Dy = oban.ll_to_xy(dloc[0]*deg2rad,dloc[1]*deg2rad,rlat,rlon,map_proj)
        dxy_list.append((Dx,Dy))
        #print "Dx,Dy",Dx,Dy
        
        # Find closest gate to disdrometer location
        # First find the closest azimuth
        
#         theta_dis = -N.arctan2(Dy,Dx)+N.pi/2.       # I *think* this is correct
#         
#         #print "azimuth of disdrometer: ",theta_dis/deg2rad
#         #Find closest index of azimuth
#         
#         theta_diff = theta_dis-theta_c[0,:]
#         theta_index = N.argmin(N.abs(theta_diff))
#         
#         # Now find index along range; first compute slant range
#         
#         srange_dis = oban.computeslantrange([Dx],[Dy],el)
#         rad_diff = srange_dis-rad_c[:,0]
#         srange_index = N.argmin(N.abs(rad_diff))
#         #print "range of disdrometer: ",srange_dis
        
        # Try another way to get the indices (the above seems to have problems).
        # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why the 
        # above doesn't work in all situations, but it sometimes picks a gate adjacent to the one
        # we want...
        
        field_list=[]
        
        for field in fieldlist:
            field = field.swapaxes(0,1)
            field = N.ma.masked_invalid(field)
        
            if(Cressman):
                #print xplt_c.shape,yplt_c.shape
                field_D = oban.Cresmn(xrad_c,yrad_c,field,Dx,Dy,roi)
                #print dBZ_D,dBZ_D.shape
            else: # First, compute the euclidian distance of the x,y location of the disdrometer to each of the radar gate centers
                distance = N.sqrt(((Dx-xrad_c)**2.+(Dy-yrad_c)**2.))
                
                # Now, find the index of the closest radar gate
                srange_index,theta_index = N.unravel_index(distance.argmin(),distance.shape)
                #print "srange_index,theta_index",srange_index,theta_index
                
                #print "Distance to closest gate: ",distance[srange_index,theta_index]
                
                # Finally, grab field at closest gate to disdrometer
                # Average the field in the closest gate and surrounding 8 gates if desired
                if(distance[srange_index,theta_index] > 3000.): # Disdrometer is probably not covered by sweep, set value to N.nan
                    field_D = N.nan
                else:
                    if(not average_gates):
                        field_D = field[srange_index,theta_index]
                    else:
    #                     field_D = (1./9.)*(field[srange_index-1,theta_index-1]+field[srange_index,theta_index-1]+
    #                                   field[srange_index-1,theta_index]+field[srange_index,theta_index]+
    #                                   field[srange_index+1,theta_index]+field[srange_index,theta_index+1]+
    #                                   field[srange_index+1,theta_index+1]+field[srange_index+1,theta_index-1]+
    #                                   field[srange_index-1,theta_index+1])
                        field_D = N.nanmean(field[srange_index-1:srange_index+2,theta_index-1:theta_index+2])
            
            #print "Value of field at disdrometer = ","%.2f"%field_D
            
            field_list.append(field_D)
        
        field_D_list.append(field_list)
    
    field_D_arr = N.array(field_D_list)
    
    return dxy_list,field_D_arr
    
def avgwind(winddirs,windspds,avgintv,gusts=True,gustintv=3,center=True):
    """Given a timeseries of wind directions and speeds, and an interval for averaging,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""
    
    windspdsavg = pd.rolling_mean(windspds,avgintv,center=center)
    if(gusts):
        windgusts = pd.rolling_mean(windspds,gustintv,center=center)
        windgustsavg = pd.rolling_max(windgusts,avgintv,center=center)
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds*N.cos(N.deg2rad(-winddirs+270.))
    vs = windspds*N.sin(N.deg2rad(-winddirs+270.))

    # Linearly interpolate for bad values of us,vs
    #us = interpnan1D(us)
    #vs = interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.rolling_mean(us,avgintv,center=center)
    vsavg = pd.rolling_mean(vs,avgintv,center=center)
    windspdsavgvec = N.sqrt(usavg**2.+vsavg**2.)
    winddirsavgvec = (270.0-(180./N.pi)*N.arctan2(vsavg,usavg))%360. # Need to use %360 to keep wind dir between 0 and 360 degrees

    return windspdsavg,windspdsavgvec,winddirsavgvec,windgusts,windgustsavg
    
def resamplewind(datetimes,offset,winddirs,windspds,DSD_intervalstr,gusts=True,gustintvstr='3S',center=False):
    """Given a timeseries of wind directions and speeds, and an interval for resampling,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""
    
    windspdsavg = pd.Series(data=windspds,index=datetimes).resample(DSD_intervalstr,label='right',closed='right',base=offset).mean()
    if(gusts):
        windgusts = pd.Series(data=windspds,index=datetimes).resample(gustintvstr,label='right',closed='right',base=offset).mean()
        windgustsavg = windgusts.resample(DSD_intervalstr,label='right',closed='right',base=offset).max()
        windgusts = windgusts
        windgustsavg = windgustsavg
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds*N.cos(N.deg2rad(-winddirs+270.))
    vs = windspds*N.sin(N.deg2rad(-winddirs+270.))

    # Linearly interpolate for bad values of us,vs
    #us = interpnan1D(us)
    #vs = interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.Series(data=us,index=datetimes).resample(DSD_intervalstr,label='right',closed='right',base=offset).mean()
    vsavg = pd.Series(data=vs,index=datetimes).resample(DSD_intervalstr,label='right',closed='right',base=offset).mean()
    windspdsavgvec = N.sqrt(usavg**2.+vsavg**2.)
    winddirsavgvec = (270.0-(180./N.pi)*N.arctan2(vsavg,usavg))%360. # Need to use %360 to keep wind dir between 0 and 360 degrees
    
    #unit average wind direction
    unit_us = us/windspds
    unit_vs = vs/windspds
    unit_usavg = pd.Series(data=unit_us,index=datetimes).resample(DSD_intervalstr,label='right',closed='right',base=offset).mean()
    unit_vsavg = pd.Series(data=unit_vs,index=datetimes).resample(DSD_intervalstr,label='right',closed='right',base=offset).mean()
    winddirsunitavgvec = (270.0-(180./N.pi)*N.arctan2(unit_vsavg,unit_usavg))%360. # Need to use %360 to keep wind dir between 0 and 360 degrees

    return windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg,usavg,vsavg,unit_usavg,unit_vsavg
    
#-----------------------------------------------------------------------------------------
#
# Below are several routines for reading and analyzing VORTEX-2 disdrometer probe data.
# Eventually will merge this with the current pyPIPS code package.  There are kept here
# for reference.
#
#-----------------------------------------------------------------------------------------
    
def readthermoCU(filename):
    """Reads in thermodynamic data from a CU probe"""
    
    thermofile = open(filename,'r')
    
    # Read in the header information
    
    dummy = thermofile.readline()
    dummy = thermofile.readline()
    dummy = thermofile.readline()
    dummy = thermofile.readline()
    
    datetimes=[]
    temps=[]
    rhs=[]
    pressures=[]
    bws=[]
    bwd=[]
    sws=[]
    swd=[]
    
    for l in thermofile:
        mysplit = shlex.shlex(l)
        mysplit.whitespace += ','
        mysplit.whitespace_split = True
        #line = l.strip().split(',')
        line = list(mysplit)
        datetimestring = line[0]
        #print datetimestring
        year = int(datetimestring[1:5])
        month = int(datetimestring[6:8])
        day = int(datetimestring[9:11])
        hour = int(datetimestring[12:14])
        min = int(datetimestring[15:17])
        sec = int(datetimestring[18:20])
        newdatetime = datetime(year,month,day,hour,min,sec)
        datetimes.append(newdatetime)
        temps.append(N.float(line[2]))
        rhs.append(N.float(line[3]))
        bws.append(N.float(line[4]))
        bwd.append(N.float(line[5]))
        try:
            temp = N.float(line[6])
        except:
            temp = N.nan
        sws.append(temp)
        try:
            temp = N.float(line[7])
        except:
            temp = N.nan
        swd.append(temp)
        try:
            pressures.append(N.float(line[17][1:-1]))
        except:
            # Temporary hack since pressures aren't in some of the files. Set it to something here.
            pressures.append(962.00)
        
    temps = N.array(temps)
    rhs = N.array(rhs)
    pressures = N.array(pressures)
    bws = N.array(bws)
    bwd = N.array(bwd)
    sws = N.array(sws)
    swd = N.array(swd)
    
    
    #N.set_printoptions(threshold='nan')
    #print bws.shape
    
    #print "temps = "
    #print temps
    return datetimes,temps,rhs,pressures,bws,bwd,sws,swd
    
def readdataCU(filename):
    """This function reads in disdrometer data from the CU V2 2009 mission and computes number
       concentration from the particle counts and fall velocities for each bin"""
       
    disfile = open(filename,'r')
    
    rainvd = assignfallspeed(avg_diameter)  # Create V-D relationship for rain based on Terry Schuur's relationship
    
    dates=[]
    times=[]
    intensities=[]
    preciptots=[]
    weathercodes=[]
    reflectivities=[]
    visibilities=[]
    pcounts=[]
    pcounts2=[]
    sensortemps=[]
    
    countsMatrix=[]     # Will contain the number of drops in the diameter vs. fall speed matrix for each time
                        # i.e. a (numtimes,32,32) array
                
    concentrations=[]   # Will contain computed number concentrations for each size bin for each time
                        # i.e. a (numtimes,32) array
                        
    # First, count number of records in file to determine size of time dimension
    
    for l in disfile:
        line = l.strip().split(';')
        #print line
        dates.append(line[0])
        time = line[1]
        times.append(time)
        intensities.append(N.float(line[2]))
        preciptots.append(N.float(line[3]))
        weathercodes.append(N.int(line[4]))
        reflectivities.append(N.float(line[5]))
        visibilities.append(N.float(line[6]))
        pcounts.append(N.int(line[7]))
        sensortemps.append(N.float(line[8]))
        
        # spectra start with token position 10
        
        spectrum = [int(x) if x != '' else 0 for x in line[10:-1]]
        # Looks like the first value might be left off: prepend a 0 to the beginning if the spectrum isn't empty
        if(len(spectrum) > 0):
            spectrum.insert(0,0)
        
        # Now create an array out of the spectrum and reshape it to 32x32
        spectrum = N.array(spectrum)
        if(spectrum.size > 0):
            spectrum = spectrum.reshape((32,32))
        else:
            spectrum = N.zeros((32,32),dtype='int')
        
        # Append spectrum (corrected or not) to spectrum time list
        
        countsMatrix.append(spectrum)
    
    # Recast countsMatrix as numpy array
    
    countsMatrix = N.dstack(countsMatrix)
    countsMatrix = N.rollaxis(countsMatrix,2,0)
        
    # Perform Katja's QC routines if desired (should not be used in combination with my methods above,
    # most are redundant anyway).
    
    X,Y = N.meshgrid(avg_diameter,fall_bins)
    flaggedtimes=N.zeros(len(times),dtype=bool)
        
    if(use_strongwindQC):
        countsMatrix,flaggedtimes = strongwindQC(countsMatrix)
        
    if(use_splashingQC):
        countsMatrix = splashingQC(countsMatrix)
        
    if(use_marginQC):
        countsMatrix = marginQC(countsMatrix)
    
    if(use_rainfallspeedQC):
        countsMatrix = rainfallspeedQC(countsMatrix,rainvd,falltol,maskhigh,masklow)
    
    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)
    
    if(use_rainonlyQC):
        countsMatrix = rainonlyQC(countsMatrix)
    
    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)
    
    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)
    
    #print flaggedtimes
    
    if(plot_QC):
        for t in range(N.size(countsMatrix,axis=0)):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plt.title('Fall speed vs. diameter for time '+times[t])
            
            countsplot = ma.masked_where(countsMatrix[t,:] <= 0, countsMatrix[t,:])
            
            C = ax1.pcolor(min_diameter,min_fall_bins,countsplot,vmin=1,vmax=50,edgecolors='w')
            ax1.plot(avg_diameter,rainvd,c='r')
            #ax1.scatter(X[0:10,20:31],Y[0:10,20:31],c='r',marker='x')
            fig.colorbar(C)
    
            if(len(flaggedtimes) > 0 and flaggedtimes[t]):
                ax1.text(0.5,0.5,'Flagged for strong wind contamination!',
                         horizontalalignment='center',
                         verticalalignment='center',color='y',
                         transform = ax1.transAxes)
            if(plot_strongwindQC):
                ax1.scatter(X[strongwindmask],Y[strongwindmask],c='r',marker='x',alpha=1.0)
            if(plot_splashingQC):
                ax1.scatter(X[splashmask],Y[splashmask],c='w',marker='o',alpha=0.75)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(splashmask,mask=-splashmask),cmap=cm.Reds,alpha=0.1)
            if(plot_marginQC):
                ax1.scatter(X[marginmask],Y[marginmask],c='g',marker='x',alpha=0.1)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(marginmask,mask=-marginmask),cmap=cm.Reds,alpha=0.1)
            if(plot_rainfallspeedQC):
                ax1.scatter(X[fallspeedmask],Y[fallspeedmask],c='k',marker='x',alpha=0.5)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(fallspeedmask,mask=-fallspeedmask),cmap=cm.gray,alpha=0.1)
            if(plot_rainonlyQC):
                ax1.scatter(X[rainonlymask],Y[rainonlymask],c='g',marker='x',alpha=0.5)
    
            ax1.set_xlim(0.0,9.0)
            ax1.xaxis.set_major_locator(MultipleLocator(1.0))
            ax1.set_xlabel('diameter (mm)')
            ax1.set_ylim(0.0,15.0)
            ax1.yaxis.set_major_locator(MultipleLocator(1.0))
            ax1.set_ylabel('fall speed (m/s)')
            
            #plt.savefig('/Users/ddawson/Dropbox/SLS_2012/vd_images/temp_'+times[t]+'.png')
            
    # Find total number of non-masked particles
    
    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)
    
    #print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t,time in enumerate(times):
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy,vspectrum = N.meshgrid(avg_diameter,fall_bins)
            dspectrum = countsMatrix[t,:]
        else:
            vspectrum = rainvd
            dspectrum = countsMatrix[t,:].sum(axis=0)   # Sum up particles for each diameter bin
            
            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]
    
        # Now compute the number concentration using the assumed fall speeds, sensor area, and sampling interval
        # Units are #/m^3/mm
        if(spectrum.size > 0 and not flaggedtimes[t]):
            concentration = dspectrum/(vspectrum*sampling_period*eff_sensor_area*(max_diameter-min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
            #print "concentration.shape"
            #print concentration.shape
        elif(not flaggedtimes[t]):
            concentration = N.zeros_like(avg_diameter)
        else:
            concentration = -999.*N.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999.,concentration)
        
        # Throw out particles above and below a certain diameter if desired
        
#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)
        
        #print "Number of particles counted vs. summed number: "
        #print N.int(line[7]),N.sum(dspectrum)
        
        #print concentration
        
        concentrations.append(concentration)
    
    concentrations = ma.array(concentrations)
    #print "concentrations: ",concentrations
    pcounts = N.array(pcounts)
    #pcounts2 = ma.array(pcounts2)
    
    return dates,times,intensities,preciptots,weathercodes,reflectivities,visibilities,pcounts,pcounts2, \
               sensortemps,min_diameter,max_diameter,avg_diameter,concentrations

def readdataCU_2010(filename):
    """This function reads in disdrometer data from the CU V2 2010 mission and computes number
       concentration from the particle counts and fall velocities for each bin"""
       
    disfile = open(filename,'r')
    
    rainvd = assignfallspeed(avg_diameter)  # Create V-D relationship for rain based on Terry Schuur's relationship
    
    dates=[]
    times=[]
    intensities=[]
    preciptots=[]
    weathercodes1=[]
    weathercodes2=[]
    weathercodes3=[]
    weathercodes4=[]
    reflectivities=[]
    visibilities=[]
    sampleintervals=[]
    pcounts=[]
    pcounts2=[]
    sensortemps=[]
    
    countsMatrix=[]     # Will contain the number of drops in the diameter vs. fall speed matrix for each time
                        # i.e. a (numtimes,32,32) array
                
    concentrations=[]   # Will contain computed number concentrations for each size bin for each time
                        # i.e. a (numtimes,32) array
    
    ntimes = 0
    sampleinterval = -1 # Have to give it some initial value
    for l in disfile:
        
        # Read in ATM4 number
        try:
            number = N.int(l[:2])
        except:
            number = -1 # Blank line
        try:
            line = l[3:]
        except:
            line = None # Blank line
        
        line = line.strip()
        
        if(line == ''):
            line = 0
        
        if(number == 1):
            intensities.append(N.float(line))
        elif(number == 2):
            preciptots.append(N.float(line))
        elif(number == 3):
            weathercodes1.append(N.int(line))
        elif(number == 4):
            weathercodes2.append(N.int(line))
        elif(number == 5):
            weathercodes3.append(line)
        elif(number == 6):
            weathercodes4.append(line)
        elif(number == 7):
            reflectivities.append(N.float(line))
        elif(number == 8):
            visibilities.append(N.float(line))
        elif(number == 9):
            sampleinterval = N.int(line)
            if(sampleinterval != 10):   # Don't read records that have a sample interval that isn't 10 seconds 
                print "sample interval = ",sampleinterval
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
        if(sampleinterval == 10): # Read the rest of the data record if sampleinterval = 10 s
            if(number == 11):
                pcount = N.float(line)
                try:
                    pcount = N.int(pcount)
                except:
                    pcount = 0
                pcounts.append(pcount)
            elif(number == 12):
                sensortemps.append(N.float(line))
            elif(number == 20):
                times.append(line)
                ntimes = ntimes + 1
                #print "time,ntimes = ",times[ntimes-1],ntimes
            elif(number == 21):
                dates.append(line)
            elif(number == 93):
                line = line.strip().split(';')
                spectrum = [float(x) if x!= '' else 0 for x in line]
                spectrum = spectrum[:-1] # Strip off bogus last value
                 
                # Now create an array out of the spectrum and reshape it to 32x32
                spectrum = N.array(spectrum,dtype='int')
                if(spectrum.size > 0):
                    spectrum = spectrum.reshape((32,32))
                else:
                    spectrum = N.zeros((32,32),dtype='int')
    
                # Append spectrum (corrected or not) to spectrum time list
    
                countsMatrix.append(spectrum)
    
    # Recast countsMatrix as numpy array
    
    countsMatrix = N.dstack(countsMatrix)
    countsMatrix = N.rollaxis(countsMatrix,2,0)
        
    # Perform Katja's QC routines if desired (should not be used in combination with my methods above,
    # most are redundant anyway).
    
    X,Y = N.meshgrid(avg_diameter,fall_bins)
    flaggedtimes=N.zeros(len(times),dtype=bool)
        
    if(use_strongwindQC):
        countsMatrix,flaggedtimes = strongwindQC(countsMatrix)
        
    if(use_splashingQC):
        countsMatrix = splashingQC(countsMatrix)
        
    if(use_marginQC):
        countsMatrix = marginQC(countsMatrix)
    
    if(use_rainfallspeedQC):
        countsMatrix = rainfallspeedQC(countsMatrix,rainvd,falltol,maskhigh,masklow)
    
    # if(use_strongwindQC):
#         countsMatrix,flaggedtimes = strongwindQC(countsMatrix)
    
    if(use_rainonlyQC):
        countsMatrix = rainonlyQC(countsMatrix)
    
    if(maskhighdiam):
        countsMatrix = maskhighdiamQC(countsMatrix)
    
    if(masklowdiam):
        countsMatrix = masklowdiamQC(countsMatrix)
    
    #print flaggedtimes
    
    if(plot_QC):
        for t in range(N.size(countsMatrix,axis=0)):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            plt.title('Fall speed vs. diameter for time '+times[t])
            
            countsplot = ma.masked_where(countsMatrix[t,:] <= 0, countsMatrix[t,:])
            
            C = ax1.pcolor(min_diameter,min_fall_bins,countsplot,vmin=1,vmax=50,edgecolors='w')
            ax1.plot(avg_diameter,rainvd,c='r')
            #ax1.scatter(X[0:10,20:31],Y[0:10,20:31],c='r',marker='x')
            fig.colorbar(C)
    
            if(len(flaggedtimes) > 0 and flaggedtimes[t]):
                ax1.text(0.5,0.5,'Flagged for strong wind contamination!',
                         horizontalalignment='center',
                         verticalalignment='center',color='y',
                         transform = ax1.transAxes)
            if(plot_strongwindQC):
                ax1.scatter(X[strongwindmask],Y[strongwindmask],c='r',marker='x',alpha=1.0)
            if(plot_splashingQC):
                ax1.scatter(X[splashmask],Y[splashmask],c='w',marker='o',alpha=0.75)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(splashmask,mask=-splashmask),cmap=cm.Reds,alpha=0.1)
            if(plot_marginQC):
                ax1.scatter(X[marginmask],Y[marginmask],c='g',marker='x',alpha=0.1)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(marginmask,mask=-marginmask),cmap=cm.Reds,alpha=0.1)
            if(plot_rainfallspeedQC):
                ax1.scatter(X[fallspeedmask],Y[fallspeedmask],c='k',marker='x',alpha=0.5)
                #ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(fallspeedmask,mask=-fallspeedmask),cmap=cm.gray,alpha=0.1)
            if(plot_rainonlyQC):
                ax1.scatter(X[rainonlymask],Y[rainonlymask],c='g',marker='x',alpha=0.5)
    
            ax1.set_xlim(0.0,9.0)
            ax1.xaxis.set_major_locator(MultipleLocator(1.0))
            ax1.set_xlabel('diameter (mm)')
            ax1.set_ylim(0.0,15.0)
            ax1.yaxis.set_major_locator(MultipleLocator(1.0))
            ax1.set_ylabel('fall speed (m/s)')
            
            plt.savefig('/Users/ddawson/Dropbox/temp6/temp_'+times[t]+'.png')
            
    # Find total number of non-masked particles
    
    pcount2 = countsMatrix.sum(axis=2)
    pcounts2 = pcount2.sum(axis=1)
    
    #print pcount2,pcount2.shape
    # Now, after QC, compute number concentrations for each size bin for each time
    for t,time in enumerate(times):
        # Determine whether to use actual measured fall speeds or assumed fallspeeds
        if(use_measured_fs):
            dummy,vspectrum = N.meshgrid(avg_diameter,fall_bins)
            dspectrum = countsMatrix[t,:]
        else:
            vspectrum = rainvd
            dspectrum = countsMatrix[t,:].sum(axis=0)   # Sum up particles for each diameter bin
            
            # print "time = ",time
#             for i,diam in enumerate(min_diameter):
#                 print diam,dspectrum[i]
    
        # Now compute the number concentration using the assumed fall speeds, sensor area, and sampling interval
        # Units are #/m^3/mm
        if(spectrum.size > 0 and not flaggedtimes[t]):
            concentration = dspectrum/(vspectrum*sampleintervals[t]*eff_sensor_area*(max_diameter-min_diameter))
            if(use_measured_fs):    # Sum up # concentration for each velocity bin
                concentration = concentration.sum(axis=0)
            #print "concentration.shape"
            #print concentration.shape
        elif(not flaggedtimes[t]):
            concentration = N.zeros_like(avg_diameter)
        else:
            concentration = -999.*N.ones_like(avg_diameter)
            concentration = ma.masked_where(concentration == -999.,concentration)
        
        # Throw out particles above and below a certain diameter if desired
        
#         if(maskhighdiam):
#             concentration = ma.masked_where(avg_diameter > highdiamthresh, concentration)
#             #print "avg_diameter,concentration",avg_diameter,concentration
#         if(masklowdiam):
#             concentration = ma.masked_where(avg_diameter < lowdiamthresh, concentration)
        
        #print "Number of particles counted vs. summed number: "
        #print N.int(line[7]),N.sum(dspectrum)
        
        #print concentration
        
        concentrations.append(concentration)
    
    concentrations = ma.array(concentrations)
    #print "concentrations: ",concentrations
    pcounts = N.array(pcounts)
    #pcounts2 = ma.array(pcounts2)
    
    return dates,times,intensities,preciptots,weathercodes1,reflectivities,visibilities,pcounts,pcounts2, \
               sensortemps,min_diameter,max_diameter,avg_diameter,concentrations
        
def readdataGlennetCDF(starttime,endtime,filename):
    """Reads disdrometer data in netCDF format for Glen's data."""
    
    rainvd = assignfallspeed(avg_diameter)  # Create V-D relationship for rain based on Terry Schuur's relationship
    
    dis_file = netcdf.Dataset(filename,'r')
    
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
    
    # Construct the datetime objects (turns out starts and ends are the same time -- this is a bug that
    # will be corrected later -- the end times are correct and represent the end of each 1-min DSD collection)
    
    disdate_start=[]
    disdate_end=[]
    
    dis_startindex = 0
    dis_endindex = len(year_end)
    
    for i in range(N.size(year_start)):
        distimestart = datetime(year_start[i],month_start[i],day_start[i],hour_start[i],min_start[i],sec_start[i])
        distimeend = datetime(year_end[i],month_end[i],day_end[i],hour_end[i],min_end[i],sec_end[i])
        # Find index corresponding to chosen start and endtimes
        if(date2num(distimeend) == starttime):
            dis_startindex = i
        if(date2num(distimeend) == endtime):
            dis_endindex = i
        disdate_start.append(distimestart)
        disdate_end.append(distimeend)
        
    disdates_1min = disdate_end[dis_startindex:dis_endindex+1]
            
    # Read in the bin size information
    
    min_size = dis_file.variables['min_size'][:]
    max_size = dis_file.variables['max_size'][:]
    avg_size = (min_size + max_size)/2.
    bin_width = max_size-min_size
            
    # Read in the concentration information
    # It appears that the data is in log (base 10) of number concentration.
    # Need to confirm with Glen!  Status: confirmed
    
    logNc_bin = dis_file.variables['concentration'][:]
    Nc_bin = N.where(logNc_bin == 0.0, 0.0, 10.**(logNc_bin))
    logNc_bin = N.ma.masked_where(Nc_bin <= 0.0,logNc_bin)
    Nc_bin = Nc_bin.swapaxes(0,1)
    logNc_bin = logNc_bin.swapaxes(0,1)
    
    #Restrict time dimension to start and end times
    Nc_bin = Nc_bin[:,dis_startindex:dis_endindex+1]
    logNc_bin = logNc_bin[:,dis_startindex:dis_endindex+1]
    
    # Read in the total particle count for each 1-min interval
    pcount = dis_file.variables['pcount'][:]
    pcount = pcount[dis_startindex:dis_endindex+1]
    
    #Invert the particle counts from the # conc. in each bin
    counts_bin = N.zeros_like(Nc_bin)
    for t in range(N.size(counts_bin,axis=1)):
        counts_bin[:,t] = Nc_bin[:,t]*rainvd[:]*sampling_period*eff_sensor_area*bin_width[:]
    
    # Zero out bins outside desired size range if desired
    if(maskhighdiam):
        diamindex = N.where(avg_diameter > highdiamthresh)[0][0]
        Nc_bin[diamindex:,:] = 0.0
        logNc_bin = N.ma.masked_where(Nc_bin <= 0.0,logNc_bin)
        
    if(masklowdiam):
        diamindex = N.where(avg_diameter > lowdiamthresh)[0][0]
        Nc_bin[:diamindex,:] = 0.0
        logNc_bin = N.ma.masked_where(Nc_bin <= 0.0,logNc_bin)
    
    return pcount,counts_bin,Nc_bin,logNc_bin,min_size,avg_size,max_size,bin_width,disdates_1min

def readthermoGlennetcdf(filename):
    """Reads thermodynamic and wind data from one of Glen's disdrometer probes"""

def readthermoGlen(filename):
    """Reads thermodynamic data from one of Glen's disdrometer probes"""
    
    tprh_file = open(filename,'r')
    
    # Throw out the first 4 lines, which are just header information
    
    for line in range(4):
        dummyline = tprh_file.readline()
    #print dummyline
    
    # Now, read in the data and put them into numpy arrays
    
    dates=[]
    records=[]
    temps=[]
    rhs=[]
    pressures=[]
    
    for line in tprh_file:
        tokens = line.strip().split(',')
        dates.append(tokens[0])
        records.append(N.int(tokens[1]))
        temps.append(N.float(tokens[2]))
        rhs.append(N.float(tokens[3]))
        pressures.append(N.float(tokens[4]))
    
    temps = N.array(temps)
    rhs = N.array(rhs)
    pressures = N.array(pressures)
    
    # Dates are in the following format: YYYY-MM-DD HH:MM:SS
    # where the time is CST (-6 GMT)
    # We can use the datetime module to convert the strings to a list of datetime objects
    # For some reason I can't get the pytz module to specify the time zone as CST: it always converts
    # the time zone to CDT.  Sigh...
    
    thermodates=[]
    
    for date in dates:
        year = int(date[1:5])
        month = int(date[6:8])
        day = int(date[9:11])
        hour = int(date[12:14])
        minute = int(date[15:17])
        second = int(date[18:20])
        # Now that we have all the date and time elements, we can create a datetime object
        thermotime=datetime(year,month,day,hour,minute,second)
        #newtime1=Central.localize(newtime)
        #print newtime.strftime(fmt)
        thermodates.append(thermotime)
    
    return thermodates,temps,rhs,pressures

def readwindGlen(filename):
    """Reads wind data from one of Glen's disdrometer probes"""
    
    wind_file = open(filename,'r')
    
    # Throw out the first 4 lines, which are just header information
    
    for line in range(4):
        dummyline = wind_file.readline()
    #print dummyline
    
    # Now, read in the data and put them into numpy arrays
    
    dates=[]
    records=[]
    winds=[]
    winddirs=[]
    windgusts=[]
    
    for line in wind_file:
        tokens = line.strip().split(',')
        dates.append(tokens[0])
        records.append(N.int(tokens[1]))
        winds.append(N.float(tokens[2]))
        winddirs.append(N.float(tokens[3]))
        windgusts.append(N.float(tokens[5]))
    
    winds = N.array(winds)
    winddirs = N.array(winddirs)
    windgusts = N.array(windgusts)
    
    # Dates are in the following format: YYYY-MM-DD HH:MM:SS
    # where the time is CST (-6 GMT)
    # We can use the datetime module to convert the strings to a list of datetime objects
    # For some reason I can't get the pytz module to specify the time zone as CST: it always converts
    # the time zone to CDT.  Sigh...
    
    winddates=[]
    
    for date in dates:
        year = int(date[1:5])
        month = int(date[6:8])
        day = int(date[9:11])
        hour = int(date[12:14])
        minute = int(date[15:17])
        second = int(date[18:20])
        # Now that we have all the date and time elements, we can create a datetime object
        windtime=datetime(year,month,day,hour,minute,second)
        #newtime1=Central.localize(newtime)
        #print newtime.strftime(fmt)
        winddates.append(windtime)
    
    return winddates,winds,winddirs,windgusts
