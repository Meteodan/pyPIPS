# Plot_Disdrometer
#
# This script plots several quantities based on disdrometer data from the Parsivel laser disdrometer

import netCDF4 as netcdf
import Nio
import matplotlib
#matplotlib.use('TkAgg')
import numpy as N
import numpy.ma as MA
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import *
from matplotlib.ticker import *
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
import pytz as pytz
from scipy import special
from scipy import signal
from scipy import ndimage
import modules.thermolib as thermo
import math as M
import imp
import glob
import os
import modules.obanmodule as oban
from mpl_toolkits.axes_grid1 import make_axes_locatable,host_subplot
from matplotlib import colorbar
import modules.disdrometer_module as dis
from matplotlib.font_manager import FontProperties
import modules.raymond_lowpass as raymond
import modules.DSDlib as dsd
import pickle
#import cressman as cress
import pdb
import sys
import modules.ctablesfrompyesviewer as ctables
import modules.radarmodule as radar
import pandas as pd
import modules.plotmodule as pm
from modules.utils import log,warning,fatal
import modules.DSDretrieval as DR
import modules.empirical_module as em

clevels_ref = N.arange(5.0,85.0,5.0)          # Contour levels for reflectivity (dBZ)
clevels_zdr = N.arange(0.0,6.25,0.25)         # Contour levels for Zdr (dB)
clevels_vr  = N.arange(-40.0,41.0,1.0)        # Contour levels for Vr (m/s)
cmapdBZ = ctables.__getattribute__('REF_default')
cmapzdr = cm.Reds
cmapvr = cm.RdBu_r

fieldnames = ['dBZ','ZDR','KDP','RHV','Vr']

fmt = '%Y-%m-%d %H:%M:%S %Z%z'
fmt2 = '%Y-%m-%d %H:%M:%S'
fmt3 = '%Y%m%d%H%M%S'

# Function definitions

def interpnan1D(a):
    """Replaces NaN's in a 1D array by interpolating from good values on either side"""
    ind = N.where(~N.isnan(a))[0] # indices of valid values
    return N.interp(list(range(len(a))),ind,a[ind]) # Use valid values to interpolate to invalid values

def trymax(a,default=0):
    """Tries to take a maximum of an array and returns 0 if the array is empty"""
    try:
        return a.max()
    except:
        return default


deg2rad = N.pi/180.

# Set global font size for axes and colorbar labels, etc.

font = {'size':10}
matplotlib.rc('font',**font)

fontP = FontProperties()
fontP.set_size('small')

# Thermodynamic constants

Rd=287.0		#Gas constant for dry air (J/kg/K)
Rv=461.51               #Gas constant for water vapor (J/kg/K)
cp=1005.6		#Specific heat at constant pressure of dry air (J/kg/K)
Lv=2.501e6		#Latent heat of vaporization (J/kg)
Lf=3.34e5		#Latent heat of freezing (J/kg)
Ls=Lv+Lf		#Latent heat of sublimation (J/kg)

rho_ref = 1.0           #Reference air density (kg/m^3)
PI = M.pi
rhow = 1000.0           #Density of liquid water (kg/m^3)

min_diameter = dis.min_diameter
max_diameter = dis.max_diameter
bin_width = max_diameter-min_diameter
avg_diameter = dis.avg_diameter
fall_bins = dis.fall_bins
min_fall_bins = dis.min_fall_bins

rainvd = dis.assignfallspeed(avg_diameter)  # Create V-D relationship for rain based on Terry Schuur's relationship

#Time zone stuff
Central = pytz.timezone('US/Central')
fmt = '%Y-%m-%d %H:%M:%S %Z%z'
fmt2 = '%Y-%m-%d %H:%M:%S'

timedelta5s = timedelta(seconds=5)
timedelta10s = timedelta(seconds=10)
timedelta30s = timedelta(seconds=30)
timedelta1min = timedelta(minutes=1)


#-----------------------------------------------------------------------
#
#   Dynamically import pyPIPScontrol.py or user version of it.
#
#-----------------------------------------------------------------------

def import_all_from(module_path):
    """Modified from http://grokbase.com/t/python/python-list/1172ahxp0s/from-module-import-using-import
       Loads python file at "module_path" as module and adds contents to global namespace."""
    #mod = __import__(module_name)
    mod = imp.load_source('mod',module_path)
    return mod

if len(sys.argv) > 1:   # Try to import user-defined plotcontrol module
    controlmodpath = sys.argv[1]
    log("Input file is "+controlmodpath)
    pc = import_all_from(controlmodpath)
    try:
        pc = import_all_from(controlmodpath)
        log("Successfully imported pyPIPS control parameters!")
    except:
        warning("Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
        import pyPIPScontrol as pc
else:   # Read in default plotcontrol.py
    import pyPIPScontrol as pc

# Parse command line argument and read in disdrometer/radar information from text file
# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
# version first, and if not, reading the textfile and pickling the read-in variables for future read-ins.  Think about this.

if(len(sys.argv) == 1):
    argindex = 1
elif(len(sys.argv) > 1):
    argindex = 2
else:
    sys.exit("No text input file defined! Quitting!")

with open(sys.argv[argindex],'r') as inputfile:
    # Read and discard header
    inputfile.readline()
    inputfile.readline()
    inputfile.readline()
    dis_dir = inputfile.readline().strip()
    image_dir = inputfile.readline().strip()
    inputfile.readline()
    numdis = N.int(inputfile.readline().strip().split()[0]) # Read in number of disdrometers
    
    # Read in disdrometer locations and names
    dis_list = []
    dis_ftype_list = []
    tprh_filenames = []
    wind_filenames = []
    dlocs=[]
    dis_name_list = []
    starttimes = []
    stoptimes = []
    centertimes = []
    for l in range(numdis):
        line = inputfile.readline().strip().split(',')
        dname = line[0] # Disdrometer name
        dfile = line[1] # Disdrometer filename
        starttime = line[2] # N.int(line[2])
        stoptime = line[3] # N.int(line[3])
        centertime = line[4] # N.int(line[4])
        lat = N.float(line[5])
        lon = N.float(line[6])
        alt = N.float(line[7])
        # Append to appropriate lists
        dis_list.append(dfile)
        dis_name_list.append(dname)
        starttimes.append(starttime)
        stoptimes.append(stoptime)
        centertimes.append(centertime)
        dlocs.append((lat,lon,alt))
        
    inputfile.readline()
    inputfile.readline()
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    platform = line[0]
    wavelength = N.float(line[1])
    inputfile.readline()
    radar_dir = inputfile.readline().strip()
    inputfile.readline()
    scattdir = inputfile.readline().strip()
    
     # Read in start and end times for the radar data analysis
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_int = list(map(int,line))
    starttimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    line = inputfile.readline().strip().split(',')
    line_int = list(map(int,line))
    stoptimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    
    # Read in plot window bounds
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_float = list(map(float,line))
    plotxmin = line_float[0]
    plotxmax = line_float[1]
    plotymin = line_float[2]
    plotymax = line_float[3]
    
    # Read in radar name,lat,lon
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    radar_name = line[0]
    try:
        rlat = N.float(line[1])
    except:
        rlat = None
    try:
        rlon = N.float(line[2])
    except:
        rlon = None
    try:
        ralt = N.float(line[3])
    except:
        ralt = None
    
    #print line[4]
    #print N.float(line[4])
    try:
        el_req = N.float(line[4])
        print("requested elevation angle",el_req)
    except:
        el_req = 0.5    # Default to 0.5 degrees
    
    try:
        heading = N.float(line[5])
        print("Radar heading: ",heading)
    except:
        heading = None
    
    # Read in min range,max range, min azimuth, and max azimuth for radar plotting (may deprecate this)
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_float = list(map(float,line))
    minrange = line_float[0]
    maxrange = line_float[1]
    minazim = line_float[2]
    maxazim = line_float[3]
    
    radlims = [minrange,maxrange,minazim,maxazim]
    plotlims = [plotxmin,plotxmax,plotymin,plotymax]
    
# We need the disdrometer locations. If they aren't supplied in the input control file, find them 
# from the GPS data

for index,dis_name,dis_filename,dloc in \
    zip(range(0,len(dis_list)),dis_name_list,dis_list,dlocs):
    
    if(N.int(dloc[0]) == -1):
        filepath = dis_dir+dis_filename
        GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc = dis.readPIPSloc(filepath)
        dlocs[index] = dloc
    
    print("Lat/Lon/alt of "+dis_name+": "+"{:7.4f}".format(dloc[0])+" {:7.4f}".format(dloc[1])+ \
                                          " {:3.0f}".format(dloc[2]))
                                          
    # Return time range of disdrometer deployment
    
    first,last = dis.readPIPStimerange(filepath)
    #print first[0],first[1],last[0],last[1]
    print("Start/End time of deployment for "+dis_name+": "+first[0]+" "+first[1]+" to "+ \
                                             last[0]+" "+last[1])