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
    return N.interp(range(len(a)),ind,a[ind]) # Use valid values to interpolate to invalid values

deg2rad = N.pi/180.

# Set global font size for axes and colorbar labels, etc.

font = {'size':10}
matplotlib.rc('font',**font)

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
    for l in xrange(numdis):
        line = inputfile.readline().strip().split(',')
        dname = line[0] # Disdrometer name
        dfile = line[1] # Disdrometer filename
        starttime = N.int(line[2])
        stoptime = N.int(line[3])
        centertime = N.int(line[4])
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
    line_int = map(int,line)
    starttimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    line = inputfile.readline().strip().split(',')
    line_int = map(int,line)
    stoptimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    
    # Read in plot window bounds
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_float = map(float,line)
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
        print "requested elevation angle",el_req
    except:
        el_req = 0.5    # Default to 0.5 degrees
        
    # Read in min range,max range, min azimuth, and max azimuth for radar plotting (may deprecate this)
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_float = map(float,line)
    minrange = line_float[0]
    maxrange = line_float[1]
    minazim = line_float[2]
    maxazim = line_float[3]
    
    radlims = [minrange,maxrange,minazim,maxazim]
    plotlims = [plotxmin,plotxmax,plotymin,plotymax]
    
# We need the disdrometer locations. If they aren't supplied in the input control file, find them 
# from the GPS data

for index,dis_name,dis_filename,dloc in \
    zip(xrange(0,len(dis_list)),dis_name_list,dis_list,dlocs):
    
    if(N.int(dloc[0]) == -1):
        filepath = dis_dir+dis_filename
        GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc = dis.readPIPSloc(filepath)
        dlocs[index] = dloc
    
    print "Lat/Lon/alt of "+dis_name+": "+str(dloc)

# ------
# Grab radar data for comparison if desired

if(pc.comp_radar):
    if(not pc.loadradopt):

        radar_filelist,radtimes = radar.getradarfilelist(radar_dir,starttime=starttimerad,
                                  stoptime=stoptimerad,platform=platform,el_req=el_req)        
        # Read in the radar data, sweep by sweep.
        
        fields_tlist = []
        range_start_tlist = []
        range_tlist = []
        azimuth_start_tlist = []
        azimuth_tlist = []
        rlat_tlist = []
        rlon_tlist = []
        ralt_tlist = []
        el_tlist = []
        fields_D_tlist = []
        dxy_tlist = []
        
        for index,path,sweeptime in zip(xrange(len(radar_filelist)),radar_filelist,radtimes):
            print "Processing file: "+path
            
            if(platform == 'NEXRAD'):
                outfieldnames,fields,range_start,radar_range,azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad = \
                radar.readCFRadial(True,el_req,rlat,rlon,ralt,path,sweeptime,fieldnames)
            else:
                outfieldnames,fields,range_start,radar_range,azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad = \
                radar.readUMXPnc(path,sweeptime,fieldnames)

        
            fields_tlist.append(fields)
            range_start_tlist.append(range_start)
            range_tlist.append(radar_range)
            azimuth_start_tlist.append(azimuth_start_rad)
            azimuth_tlist.append(azimuth_rad)
            rlat_tlist.append(rlat_rad)
            rlon_tlist.append(rlon_rad)
            ralt_tlist.append(ralt)
            el_tlist.append(el_rad)
            
            dxy_list,fields_D = dis.rad2DD2(fields,range_start,radar_range,
            azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad,dlocs)
            
            fields_D_tlist.append(fields_D)
            dxy_tlist.append(dxy_list)
        
        fields_tarr = N.array(fields_tlist)
        range_start_tarr = N.array(range_start_tlist)
        range_tarr = N.array(range_tlist)
        azimuth_start_tarr = N.array(azimuth_start_tlist)
        azimuth_tarr = N.array(azimuth_tlist)
        rlat_tarr = N.array(rlat_tlist)
        rlon_tarr = N.array(rlon_tlist)
        ralt_tarr = N.array(ralt_tlist)
        el_tarr = N.array(el_tlist)
        fields_D_tarr = N.array(fields_D_tlist)
        dxy_tarr = N.array(dxy_tlist)
        
        if(pc.saveradopt):
            raddate_file=open(radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
                              stoptimerad.strftime(fmt3).strip()+'.txt','w')
            pickle.dump(radtimes,raddate_file)
            pickle.dump(radar_filelist,raddate_file)
            pickle.dump(outfieldnames,raddate_file)
            raddate_file.close()
            
            radnpz_filename = radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
                              stoptimerad.strftime(fmt3).strip()+str(el_req)+'.npz'           
            savevars={}
            savevars['fields_tarr'] = fields_tarr
            savevars['range_start_tarr'] = range_start_tarr
            savevars['range_tarr'] = range_tarr
            savevars['azimuth_start_tarr'] = azimuth_start_tarr
            savevars['azimuth_tarr'] = azimuth_tarr
            savevars['rlat_tarr'] = rlat_tarr
            savevars['rlon_tarr'] = rlon_tarr
            savevars['ralt_tarr'] = ralt_tarr
            savevars['el_tarr'] = el_tarr
            savevars['fields_D_tarr'] = fields_D_tarr
            savevars['dxy_tarr'] = dxy_tarr
            
            N.savez(radnpz_filename,**savevars)
        
    if(pc.loadradopt):
        raddate_file=open(radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+        \
                              stoptimerad.strftime(fmt3).strip()+'.txt','r')
        radtimes = pickle.load(raddate_file)
        radar_filelist = pickle.load(raddate_file)
        outfieldnames = pickle.load(raddate_file)
        raddate_file.close()
        
        radnpz_filename = radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
                              stoptimerad.strftime(fmt3).strip()+str(el_req)+'.npz'
        radnpz_file = N.load(radnpz_filename)
        
        fields_tarr = radnpz_file['fields_tarr']
        range_start_tarr = radnpz_file['range_start_tarr']
        range_tarr = radnpz_file['range_tarr']
        azimuth_start_tarr = radnpz_file['azimuth_start_tarr']
        azimuth_tarr = radnpz_file['azimuth_tarr']
        rlat_tarr = radnpz_file['rlat_tarr']
        rlon_tarr = radnpz_file['rlon_tarr']
        ralt_tarr = radnpz_file['ralt_tarr']
        el_tarr = radnpz_file['el_tarr']
        fields_D_tarr = radnpz_file['fields_D_tarr']
        dxy_tarr = radnpz_file['dxy_tarr']        
        
# Plot radar scans with overlaid disdrometer locations if desired
    if(pc.plot_radar):
        if (not os.path.exists(image_dir+'radar_PPI/')):
            os.makedirs(image_dir+'radar_PPI/')

        print "Plotting radar scans with overlaid disdrometer locations and data."
        for index,path,sweeptime in zip(xrange(len(radar_filelist)),radar_filelist,radtimes):
        
            fields_arr = fields_tarr[index]
            fields_D_arr = fields_D_tarr[index]
            rlat_rad = rlat_tarr[index]
            rlon_rad = rlon_tarr[index]
            ralt = ralt_tarr[index]
            el_rad = el_tarr[index]
            range_start = range_start_tarr[index]
            radar_range = range_tarr[index]
            azimuth_start_rad = azimuth_start_tarr[index]
            azimuth_rad = azimuth_tarr[index]
            dxy_arr = dxy_tarr[index]
            
            fields_list = list(fields_arr)  # list(array) where array is 2D or higher yields a list of arrays!
            dxy_list = dxy_arr.tolist()
            fields_D_list = fields_D_arr.T.tolist()   # array.tolist() yields a nested list for a 2D+ array!
                                                      # Note, need to transpose here because plotsweep expects the first
                                                      # axis to be the field (i.e. dBZ, ZDR, etc.) and the second axis to be 
                                                      # the disdrometer.  What a tangled web I weave!
            
            # In most of our cases the radar location isn't going to change with time, but in the more
            # general case, this may not be true (i.e. if we are dealing with a mobile radar).        
            #print "sweeptime,rlat_rad,rlon_rad,ralt,el_rad",sweeptime.strftime(fmt),rlat_rad,rlon_rad,ralt,el_rad
            print "Radar sweep time: ",sweeptime.strftime(fmt)
            # Prepare masks for fields by reflectivity < some threshold
    
            for field,fieldname in zip(fields_list,outfieldnames):   # Probably should do this using a dictionary
                if(fieldname == 'dBZ'):
                    mask = N.where(field > 5.0,False,True)
            
            masklist = [mask,mask,mask]
        
            # Compute height and range of radar beam above disdrometer (assumes lambert conformal like plotsweep for now)
            for dloc,dname,dxy in zip(dlocs,dis_name_list,dxy_list):
                Dx,Dy = dxy
                print Dx,Dy
                h,r = oban.computeheightrangesingle(Dx,Dy,el_req*deg2rad)
                # In most of our cases the radar location isn't going to change with time, but in the more
                # general case, this may not be true (i.e. if we are dealing with a mobile radar). 
                #print "Disdrometer name,x,y,radar elevation angle,slant range, approximate beam height:"
                #print dname,Dx,Dy,el_rad/deg2rad,r,h
                if(dloc == dlocs[0] and plotxmin == -1):
                    plotlims = [Dx-25000.,Dx+25000.,Dy-30000.,Dy+20000.]
        

            figlist,gridlist = radar.plotsweep(radlims,plotlims,outfieldnames,fields_list,masklist,
                        range_start,radar_range,azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad,False,
                        pc.plot_radar,dis_name_list,dxy_list,fields_D_list)

                # Save figures
            for fieldname,fig in zip(outfieldnames,figlist):
                plt.figure(fig.number) 
            	plt.title(sweeptime.strftime(fmt2).strip())
                plt.savefig(image_dir+'radar_PPI/'+fieldname+sweeptime.strftime(fmt3).strip()+'el'+str(el_req)+'.png',dpi=200,bbox_inches='tight')
                plt.close(fig)

# Outer disdrometer (and deployment) loop

for index,dis_filename,dis_name,starttime,stoptime,centertime,dloc in \
    zip(xrange(0,len(dis_list)),dis_list,dis_name_list,starttimes,stoptimes,centertimes,dlocs):

    if(pc.timetospace):
        centertime = centertimes[index]

    if (not os.path.exists(image_dir)):
        os.makedirs(image_dir)
    
    if(pc.plot_conv_meteo or pc.plot_DSD_meteo):
        # Create the directory for the plots if it doesn't exist
        meteogram_image_dir = image_dir+'meteograms/'
        if (not os.path.exists(meteogram_image_dir)):
            os.makedirs(meteogram_image_dir)

    
    # Read in the disdrometer data file
    filepath = dis_dir+dis_filename
    datetimesUTC,pdatetimesUTC,flaggedtimes,intensities,preciptots,reflectivities,pcounts,pcounts2, \
               sensortemps,concentrations,onedrop_concentrations,countsMatrix,windspds,winddirrels,winddirabss, \
            winddiags,fasttemps,slowtemps,dewpoints,RHs_derived,RHs,pressures,compass_dirs,    \
            GPS_lats,GPS_lons,GPS_stats,GPS_alts,voltages,DSD_interval,DSD_intervalstr,DSD_index = dis.readPIPS(filepath,rainonlyqc=True,DSD_interval=pc.DSD_interval)

    DSD_interval_td = timedelta(seconds=DSD_interval)
    DSD_halfinterval_td = timedelta(seconds=DSD_interval/2.)

    # Convert to numpy arrays
    windspds = N.array(windspds)
    winddirrels = N.array(winddirrels)
    winddirabss = N.array(winddirabss)
    fasttemps = N.array(fasttemps)
    slowtemps = N.array(slowtemps)
    dewpoints = N.array(dewpoints)
    RHs_derived = N.array(RHs_derived)
    RHs = N.array(RHs)
    pressures = N.array(pressures)
    pcounts = N.array(pcounts)
    pcounts2 = N.array(pcounts2)
    intensities = N.array(intensities)

    Nc_bin = concentrations.T
    dropperbin = onedrop_concentrations.T
        
    logNc_bin = N.ma.log10(Nc_bin)
    logNc_bin = N.ma.masked_where(Nc_bin < dropperbin, logNc_bin)

    # Compute potential temperature, water vapor mixing ratio, and density
    pt = thermo.caltheta(pressures*100.,slowtemps+273.15)
    qv = thermo.calqv(RHs/100.,pressures*100.,slowtemps+273.15)
    rho=thermo.calrho(pressures*100.,pt,qv)
        
    # Determine start and end times/indices for analysis
    
    datetimesUTCnums = date2num(datetimesUTC)
    pdatetimesUTCnums = date2num(pdatetimesUTC)
    
    if(starttime == -1):
        startindex = 0
        pstartindex = 0
        starttime = datetimesUTCnums[startindex]
        pstarttime = pdatetimesUTCnums[startindex]
    else:
        try:
            starttime = date2num(datetime(N.int(starttime[:4]),N.int(starttime[4:6]),
                        N.int(starttime[6:8]),N.int(starttime[8:10]),N.int(starttime[10:12])))
            pstarttime = starttime
            startindex = next(i for i,t in enumerate(datetimesUTCnums) if t >= starttime)
            pstartindex = next(i for i,t in enumerate(pdatetimesUTCnums) if t >= starttime)
        except:
            startindex = 0
            pstartindex = 0
            starttime = datetimesUTCnums[startindex]
            pstarttime = pdatetimesUTCnums[startindex]
    
    if(stoptime == -1):
        stopindex = N.size(datetimesUTCnums)-1
        pstopindex = N.size(pdatetimesUTCnums)-1
        stoptime = datetimesUTCnums[stopindex]
        pstoptime = pdatetimesUTCnums[pstopindex]
    else:
        try:
            stoptime = date2num(datetime(N.int(stoptime[:4]),N.int(stoptime[4:6]),
                        N.int(stoptime[6:8]),N.int(stoptime[8:10]),N.int(stoptime[10:12])))
            pstoptime = stoptime
            stopindex = next(i for i,t in enumerate(datetimesUTCnums) if t >= stoptime)
            pstopindex = next(i for i,t in enumerate(pdatetimesUTCnums) if t >= stoptime)
        except:
            stopindex = N.size(datetimesUTCnums)-1
            pstopindex = N.size(pdatetimesUTCnums)-1
            stoptime = datetimesUTCnums[stopindex]
            pstoptime = pdatetimesUTCnums[pstopindex]
    
    disdates = pdatetimesUTC
    disdates_start = [x-DSD_interval_td for x in pdatetimesUTC]
    disdates_start.append(disdates_start[-1]+DSD_interval_td)    # Add an extra 10 sec for the last time bin boundary
    disdates_avg = [x-DSD_halfinterval_td for x in pdatetimesUTC]
    
    if(pc.timetospace and pc.comp_radar):
        delta_t_rad = [(x-centertime).total_seconds() for x in radtimes]
        delta_t_dis = [(x-centertime).total_seconds() for x in disdates]
        delta_t_dis_start = [(x-centertime).total_seconds() for x in disdates_start]
        delta_t_dis_avg = [(x-centertime).total_seconds() for x in disdates_avg]
#         if(plot_thermo):
#             delta_t_thermo = [(x-centertime).total_seconds() for x in thermodates]
#             delta_t_thermo_1min = [(x-centertime).total_seconds() for x in thermodates_1min]
#             delta_t_wind = [(x-centertime).total_seconds() for x in winddates]
#             delta_t_wind_1min = [(x-centertime).total_seconds() for x in winddates_1min]
        
        xloc_rad = [(x*pc.stormmotion)/1000.0 for x in delta_t_rad]
        xloc_dis = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis]
        xloc_dis_start = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis_start]
        xloc_dis_avg = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis_avg]
#         if(plot_thermo):
#             xloc_thermo = [(x*pc.stormmotion)/1000.0 for x in delta_t_thermo]
#             xloc_thermo_1min = [(x*pc.stormmotion)/1000.0 for x in delta_t_thermo_1min]
#             xloc_wind = [(x*pc.stormmotion)/1000.0 for x in delta_t_wind]
#             xloc_wind_1min = [(x*pc.stormmotion)/1000.0 for x in delta_t_wind_1min]
        
        plotx_rad = N.array(xloc_rad)
        plotx_dis = N.array(xloc_dis)
        plotx_dis_start = N.array(xloc_dis_start)
        plotx_dis_avg = N.array(xloc_dis_avg)
#         if(plot_thermo):
#             plotx_thermo = N.array(xloc_thermo)
#             plotx_thermo_1min = N.array(xloc_thermo_1min)
#             plotx_wind = N.array(xloc_wind)
#             plotx_wind_1min = N.array(xloc_wind_1min)
    else:
        if(pc.comp_radar):
            plotx_rad = date2num(radtimes)
        plotx_dis = date2num(disdates)
        plotx_dis_start = date2num(disdates_start)
        plotx_dis_avg = date2num(disdates_avg)
#         if(plot_thermo):
#             plotx_thermo = date2num(thermodates)
#             plotx_thermo_1min = date2num(thermodates_1min)
#             plotx_wind = date2num(winddates)
#             plotx_wind_1min = date2num(winddates_1min)

    
    # Resample the 1-s data corresponding to each integrated DSD (for the given DSD interval)
#     pindex1s = [i for i,x in enumerate(datetimesUTCnums) if x in pdatetimesUTCnums]
#     pstartindex1s = pindex1s[0]
#     pstopindex1s = len(datetimesUTC) # pindex1s[-1]
    sec_offset = pdatetimesUTC[0].second
    #print "startindex1,datetimesUTC[startindex1],startindex2,pdatetimesUTC[startindex2]",pstartindex1s,datetimesUTC[pstartindex1s],0,pdatetimesUTC[0]
#     dummy,windspds_tDSD,winddirs_tDSD,dummy,dummy = dis.resamplewind(datetimesUTC[pstartindex1s:pstopindex1s+1],winddirabss[pstartindex1s:pstopindex1s+1],
#                     windspds[pstartindex1s:pstopindex1s+1],DSD_intervalstr,gusts=False,gustintvstr='3S',center=False)
#     
#     fasttemps_tDSD = pd.Series(data=fasttemps[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     slowtemps_tDSD = pd.Series(data=slowtemps[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     dewpoints_tDSD = pd.Series(data=dewpoints[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     RHs_derived_tDSD = pd.Series(data=RHs_derived[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     RHs_tDSD = pd.Series(data=RHs[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     pressures_tDSD = pd.Series(data=pressures[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().values
#     rho_tDSD = pd.Series(data=rho[pstartindex1s:pstopindex1s+1],index=datetimesUTC[pstartindex1s:pstopindex1s+1]).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean()
#     tempdatetime = rho_tDSD.index.to_pydatetime()
#     rho_tDSD = rho_tDSD.loc[DSD_index].values
#     print len(rho_tDSD)
    #print len(tempdatetime),tempdatetime

    dummy,windspds_tDSD,winddirs_tDSD,dummy,dummy,dummy,dummy,dummy,dummy,dummy = dis.resamplewind(datetimesUTC,sec_offset,winddirabss,
                    windspds,DSD_intervalstr,gusts=False,gustintvstr='3S',center=False)
    
    windspds_tDSD = windspds_tDSD.loc[DSD_index].values
    winddirs_tDSD = winddirs_tDSD.loc[DSD_index].values
    
    fasttemps_tDSD = pd.Series(data=fasttemps,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    slowtemps_tDSD = pd.Series(data=slowtemps,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    dewpoints_tDSD = pd.Series(data=dewpoints,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    RHs_derived_tDSD = pd.Series(data=RHs_derived,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    RHs_tDSD = pd.Series(data=RHs,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    pressures_tDSD = pd.Series(data=pressures,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    rho_tDSD = pd.Series(data=rho,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
    
    # Conventional data meteogram plotting
    if(pc.plot_conv_meteo):    
        # Interval for plotting in seconds.  Setting to larger intervals is useful to avoid
        # plot clutter for long datasets, but information will be lost.
        plotinterval = 10
    
        # Plot wind meteogram
        windavgintv = 60
        windgustintv = 3
        
        times = datetimesUTCnums[startindex:stopindex+1:plotinterval]
        
        # Compute wind speed and direction, and wind gusts
        windspdsavg,windspdsavgvec,winddirsavgvec,windgusts,windgustsavg = dis.avgwind(winddirabss,
                windspds,windavgintv,gusts=True,gustintv=windgustintv,center=False)

        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        #plt.title('Wind speed (5-min mean) and gust (5-min max of 3-s mean)')

        fields = [windspdsavg[startindex:stopindex+1:plotinterval],
                  windgustsavg[startindex:stopindex+1:plotinterval]]
        fieldparamdict1 = {'type':'fill_between','linestyle':'-','color':'b','alpha':0.5,'plotmin':0}
        fieldparamdict2 = {'type':'line','linestyle':'none','marker':'o','color':'r','ms':2}
        fieldparamdicts = [fieldparamdict1,fieldparamdict2]
        xvals = [times] * len(fields) # Should really do this inside call to plotmeteogram
        ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
        
        fields = [winddirsavgvec[startindex:stopindex+1:plotinterval]]
        fieldparamdict1 = {'type':'line','linestyle':'none','marker':'o','color':'k','ms':2}
        fieldparamdicts = [fieldparamdict1]
        xvals = [times] * len(fields)
        ax2 = pm.plotmeteogram(ax2,xvals,fields,fieldparamdicts)
        
        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
                        'minorxlocator':pc.minorlocator,'axeslimits':[[starttime,stoptime],[0.0,15.0]],
                        'axeslabels':[pc.timelabel,r'wind speed (m s$^{-1}$)']}
        axparamdict2 = {'majorylocator':MultipleLocator(45.),'axeslimits':[None,[0.0,360.0]],
                        'axeslabels':[None,r'Wind direction ($^{\circ}$C)']}
        axparamdicts = [axparamdict1,axparamdict2]
        ax1,ax2 = pm.set_meteogram_axes([ax1,ax2],axparamdicts)
        
        plt.savefig(image_dir+'meteograms/'+dis_name+'_wind.png',dpi=300)
        
        # Plot temperature and dewpoint
        tavgintv = 10 # Currently not used

        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_subplot(111)
        
        fields = [fasttemps[startindex:stopindex+1:plotinterval],
                  dewpoints[startindex:stopindex+1:plotinterval]]
        fieldparamdict1 = {'type':'fill_between','linestyle':'-','color':'r','alpha':0.5,'plotmin':-10}
        fieldparamdict2 = {'type':'fill_between','linestyle':'-','color':'g','alpha':0.5,'plotmin':-10}
        fieldparamdicts = [fieldparamdict1,fieldparamdict2]
        xvals = [times] * len(fields)
        ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
        
        ax1.axhline(0.0,ls=':',color='k')
        
        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
                        'minorxlocator':pc.minorlocator,'axeslimits':[[starttime,stoptime],[-10.,40.0]],
                        'axeslabels':[pc.timelabel,r'Temperature ($^{\circ}$C)']}
        axparamdicts = [axparamdict1]
        ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
        
        plt.savefig(image_dir+'meteograms/'+dis_name+'_T_Td.png',dpi=300)
        
        # Plot relative humidity
        avgintv = 10 # Currently not used

        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_subplot(111)
        
        fields = [RHs[startindex:stopindex+1:plotinterval]]
        fieldparamdict1 = {'type':'fill_between','linestyle':'-','color':'g','alpha':0.5,'plotmin':0}
        fieldparamdicts = [fieldparamdict1]
        xvals = [times] * len(fields)
        ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
        
        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
                        'minorxlocator':pc.minorlocator,'axeslimits':[[starttime,stoptime],[0.,100.]],
                        'axeslabels':[pc.timelabel,'Relative Humidity (%)']}

        axparamdicts = [axparamdict1]
        ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
        
        plt.savefig(image_dir+'meteograms/'+dis_name+'_RH.png',dpi=300)
        
        # Plot station pressure
        avgintv = 1 # Currently not used

        fig = plt.figure(figsize=(5,3))
        ax1 = fig.add_subplot(111)
        
        fields = [pressures[startindex:stopindex+1:plotinterval]]
        fieldparamdict1 = {'type':'fill_between','linestyle':'-','color':'purple','alpha':0.5,'plotmin':0}
        fieldparamdicts = [fieldparamdict1]
        xvals = [times] * len(fields)
        ax1 = pm.plotmeteogram(ax1,xvals,fields,fieldparamdicts)
        
        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
                        'minorxlocator':pc.minorlocator,'axeslimits':[[starttime,stoptime],[0.,100.]],
                        'axeslabels':[pc.timelabel,'Relative Humidity (%)']}
        axparamdicts = [axparamdict1]
        ax1, = pm.set_meteogram_axes([ax1],axparamdicts)
        ax1.xaxis.set_major_locator(pc.locator)
        ax1.xaxis.set_major_formatter(pc.formatter)
        ax1.xaxis.set_minor_locator(pc.minorlocator)
        ax1.set_xlim(starttime,stoptime)
        ax1.set_xlabel(pc.timelabel)
        ax1.set_ylim(950.,1000.)
        ax1.set_ylabel('Station pressure (hPa)')
        fig.autofmt_xdate()
        plt.savefig(image_dir+'meteograms/'+dis_name+'_pres.png',dpi=300)
        

    #N.savetxt('Nc_bin1minavg.txt', Nc_bin1minavg) # Commented out for now
    
    if(pc.calc_DSD):
        # Now for the fun part.  Calculate exponential and gamma distribution size distribution parameters
        # using the method of moments, after Zhang et al. 2008 and Tokay and Short 1996
        
#         print bin_width.shape
#         print Nc_bin1minavg.shape
#         print logNc_bin1minavg.shape
#         print rho_tDSD.shape
    
        synthbins,exp_DSD,gam_DSD,dis_DSD = dis.calc_DSD(min_diameter,avg_diameter,
                        max_diameter,bin_width,Nc_bin,logNc_bin,rho_tDSD,pc.qrQC,pc.qr_thresh)
    
        # Unpack needed values from returned tuples
    
        N_expDSD,N0_exp,lamda_exp,mu_exp,qr_exp,Ntr_exp,refl_DSD_exp,D_med_exp,D_m_exp = exp_DSD
        N_gamDSD,N0_gam,lamda_gam,mu_gam,qr_gam,Ntr_gam,refl_DSD_gam,D_med_gam,D_m_gam = gam_DSD
        Nc_bin,logNc_bin,D_med_disd,D_m_disd,D_mv_disd,D_ref_disd,QR_disd,refl_disd = dis_DSD
        N_gamDSD = N_gamDSD.T
        logN_gamDSD = N.ma.log10(N_gamDSD/1000.) # Get to log(m^-3 mm^-1)
    	logN_gamDSD = N.ma.masked_where(N_gamDSD < dropperbin, logN_gamDSD)

        if(pc.calc_dualpol):
            # Calculate polarimetric variables using the T-matrix technique
            scattfile = scattdir+'SCTT_RAIN_fw100.dat'

            dualpol_dis = dis.calpolrain(wavelength,scattfile,Nc_bin,bin_width) # use for raw distribution
#            dualpol_dis = dis.calpolrain(wavelength,scattfile,(N_gamDSD/1000.),bin_width) # use for gamma fit 

            # Unpack needed values from returned tuple
            Zh,Zv,Zhv,dBZ,ZDR,KDP,RHV = dualpol_dis
        
    if(pc.plot_DSD_meteo):
 
        # Prepare arrays for plotting
        DSDstarttimes = date2num(disdates_start[pstartindex:pstopindex+1])
        DSDmidtimes = date2num(disdates_avg[pstartindex:pstopindex+1])

        logNc_bin_plot = logNc_bin[:,pstartindex:pstopindex+1] # use for raw distribution
#        logNc_bin_plot = logN_gamDSD[:,pstartindex:pstopindex+1] # use for gamma fit

        disvars = {'min_diameter':min_diameter,'DSDstarttimes':DSDstarttimes,
                   'DSDmidtimes':DSDmidtimes,'logNc_bin':logNc_bin_plot}
        
        # Plot one meteogram for each dualpol variable, otherwise just plot a meteogram for reflectivity
        # based on 6th moment of disdrometer DSD (Rayleigh approximation).
        # Important! Currently the dualpol calculation for the disdrometer data assumes rain only
        # and otherwise ignores data in all bins larger than 9 mm.  For this reason it is recommended
        # to first set the "rain_only_QC" flag to True in disdrometer_module.py
        
        if(pc.calc_DSD):
            # 3-min average of median volume diameter and radar reflectivity (Rayleigh approx.) from disdrometer

            window = int(180./DSD_interval)
            D_0_dis_avg = pd.Series(D_med_disd).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for raw distribution
#            D_0_dis_avg = pd.rolling_mean(D_med_gam,18,center=True,win_type='triang',min_periods=1) # use for gamma fit 
            dBZ_ray_dis_avg = pd.Series(refl_disd).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for raw distribution
#            dBZ_ray_dis_avg = pd.rolling_mean(refl_DSD_gam,18,center=True,win_type='triang',min_periods=1) # use for gamma fit
            disvars['D_0_dis'] = D_0_dis_avg[pstartindex:pstopindex+1]
            disvars['dBZ_ray'] = dBZ_ray_dis_avg[pstartindex:pstopindex+1]
            
            if(pc.calc_dualpol):
                # Computed centered running averages of dualpol variables
                for varname,var in zip(['dBZ','ZDR','KDP','RHV'],[dBZ,ZDR,KDP,RHV]):
                    if(var.size):
                        var_avg = pd.Series(var).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values        
                        disvars[varname] = var_avg[pstartindex:pstopindex+1]


        # Get radar variables together for comparison, if desired
        if(pc.comp_radar):
            # At least add the reflectivity
            indexrad = outfieldnames.index('dBZ')
            dBZ_D_plt = fields_D_tarr[:,index,indexrad]
            radvars = {'radmidtimes':plotx_rad,'dBZ':dBZ_D_plt}
            # Add other polarimetric fields
            if(pc.calc_dualpol):
                for radvarname in ['ZDR','KDP','RHV']:
                    if radvarname in outfieldnames:
                        indexrad = outfieldnames.index(radvarname)
                        dualpol_rad_var = fields_D_tarr[:,index,indexrad]
                        if(dualpol_rad_var.size):
                            radvars[radvarname] = dualpol_rad_var


        # Prepare axis parameters
        timelimits = [starttime,stoptime]
        diamlimits = [0.0,9.0]
        
        axparams = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
                    'minorxlocator':pc.minorlocator,'axeslimits':[timelimits,diamlimits],
                    'majorylocator':MultipleLocator(base=1.0),'axeslabels':[None,'D (mm)']}
                    
        # Ok, now we should have everything ready to go to plot the meteograms.  Let'er rip!
                
        pm.plotDSDmeteograms(dis_name,meteogram_image_dir,axparams,disvars,radvars.copy())

	if(pc.plot_DSDs):
		if (not os.path.exists(image_dir+'DSDs/'+dis_name)):
			os.makedirs(image_dir+'DSDs/'+dis_name)
		for t in range(N.size(Nc_bin,axis=1)):
			sum = N.ma.sum(Nc_bin[:,t])
			if(sum != 0.0):
				fig1=plt.figure(figsize=(8,6))
				ax1=fig1.add_subplot(111)
				plt.title('Disdrometer 10-sec DSD Fits for Time '+disdates[t].strftime(fmt2)+' EST')
				ax1.bar(min_diameter,Nc_bin[:,t]*1000.0,max_diameter-min_diameter,log=True,color='tan')
#				ax1.bar(min_diameter,ones_1minavg[:,t]*1000.0,max_diameter-min_diameter,log=True,color='None',edgecolor='k')
				ax1.plot(synthbins,N_expDSD[t],lw=2)
				ax1.plot(avg_diameter,N_gamDSD[:,t],lw=2)
				ax1.set_yscale('log')
				ax1.set_ylim(10.**2.0,10.**8.5)
				ax1.set_xlim(0.0,9.0)
				ax1.set_xlabel('Drop Diameter (mm)')
				ax1.set_ylabel(r'N(D) $(m^{-4})$')
				ax1.text(0.50,0.8,'Shape parameter (gamma) = %2.2f'%mu_gam[t],transform=ax1.transAxes)
				ax1.text(0.50,0.75,'Median Volume Diameter (gamma) =%2.2f'%D_med_gam[t],transform=ax1.transAxes)
				ax1.text(0.50,0.7,'Reflectivity =%2.2f'%refl_disd[t]+'dBZ',transform=ax1.transAxes)
				ax1.text(0.50,0.65,'ZDR = %2.2f'%ZDR[t]+'dB',transform=ax1.transAxes)
				ax1.text(0.50,0.6,'Intensity =%2.2f'%intensities[t]+'mm/hr',transform=ax1.transAxes)
				ax1.text(0.50,0.55,'Particle count (QC) = '+str(pcounts2[t]),transform=ax1.transAxes)
				plt.savefig(image_dir+'DSDs/'+dis_name+'/'+dis_name+'_t'+str(t)+'DSD_plot.png',dpi=200,bbox_inches='tight')
				plt.close(fig1)

	
	Zh_Cao = N.arange(20,61,1)
	Zdr_Cao = 10**((-2.6857*10**-4*Zh_Cao**2)+0.04892*Zh_Cao-1.4287)
	
	if(pc.plot_scat):
		if (not os.path.exists(image_dir+'scattergrams/')):
			os.makedirs(image_dir+'scattergrams/')
		fig1=plt.figure(figsize=(8,8))
		ax1=fig1.add_subplot(111)
		plt.title('Z vs. Zdr')
		ax1.scatter(dBZ, ZDR, color='k')
		ax1.scatter(radvars['dBZ'], radvars['ZDR'], color='r')
		ax1.plot(Zh_Cao,Zdr_Cao,lw=2)
		ax1.set_ylim(-2.0,6.0)
		ax1.set_xlim(20.0,60.0)
		ax1.set_xlabel('Zh in dBZ')
		ax1.set_ylabel('ZDR in dB')
		plt.savefig(image_dir+'scattergrams/'+dis_name+'scattergrams.png',dpi=200,bbox_inches='tight')			
#plt.show()

