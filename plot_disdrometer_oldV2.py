# Plot_Disdrometer
#
# This script plots several quantities based on disdrometer data from the Parsivel laser disdrometer

import netCDF4 as netcdf
import Nio
import matplotlib
matplotlib.use('TkAgg')
import numpy as N
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
import thermolib as thermo
import math as M
import glob
import os
import obanmodule as oban
from mpl_toolkits.axes_grid1 import make_axes_locatable,host_subplot
from matplotlib import colorbar
import disdrometer_module as dis
from matplotlib.font_manager import FontProperties
import raymond_lowpass as raymond
import DSDlib as dsd
import pickle
#import cressman as cress
import pdb
import sys

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

enhancesmalldrops = False        # For testing, artificially boost concentration of small drops?
enhancefactor = 2.0            
enhancethresh = 0.562           # Drop size at and below which to enhance the small drops

saveradopt = True      # Option to save reflectivity at DSD locations to file to avoid unneccessary
                        # recalculation
loadradopt = False       # Option to load the above

#Time zone stuff
Central = pytz.timezone('US/Central')
fmt = '%Y-%m-%d %H:%M:%S %Z%z'
fmt2 = '%Y-%m-%d %H:%M:%S'

timedelta5s = timedelta(seconds=5)
timedelta10s = timedelta(seconds=10)
timedelta30s = timedelta(seconds=30)
timedelta1min = timedelta(minutes=1)


# Plot/analysis control flags

plot_opt = True        # Plot or not?
calc_DSD = True         # Calculate DSD fits (exp and gamma)?
windQC = False          # Perform additional QC on disdrometer DSDs based on wind speed?
wind_thresh = 5.0       # wind speed threshold (1-min average) above which to throw out DSDs
qrQC = True             # Perform additional QC on disdrometer DSDs based on unrealistic qr?
qr_thresh = 15.0        # Threshold of qr above which to throw out DSD (g/kg)
plot_DSDs = False      # Plot individual 1-min (or 10-sec) DSDs and fits?
plot_thermo = False     # Plot thermodynamic fields?
calc_evap = False       # Calculate evaporation rates?
calc_radar = False       # Calculate polarimetric radar variables (currently only Z and Zdr)?
comp_radar = True       # Compare with reflectivity from 88D?
reverse_times = True     # Reverse time axis?
timetospace = True      # Perform a simple time-to-space conversion on the data?
stormmotion = 12.55 # 10.0, (5 June) 12.55 (9 June), 13.0 (7 June)     # Approximate storm motion in m/s
# Approximate center (zero) times for each deployment to compute time-to-space coordinates
centertimes = [datetime(2009,6,9,17,29,45),datetime(2009,6,9,17,31,15),datetime(2009,6,9,17,30,45),datetime(2009,6,9,17,31,15)]
#centertimes = [datetime(2009,6,5,16,11,28),datetime(2009,6,5,16,06,52)] # 5 June 2009
#centertimes = [datetime(2009,6,7,17,43,54),datetime(2009,6,7,17,39,21)] # 7 June 2009
#centertimes = [datetime(2009,6,7,17,39,21)]
xstart = -2.0
xstop = 13.0
xtickintv = 1.0
avg_data = True        # average disdrometer data in time and calculate new DSD fits?
                        # For now just sum over each deployment using start and end times.

# ------ 2010/05/19 Putnam, OK supercell ------
# 
# Directories containing disdrometer data files
dis_dir_glen = '/Users/ddawson/NSF_V2_study/Disdrometer_study/Glen_newdata/2010_05_19/'
image_dir_glen = dis_dir_glen+'images/'
dis_dir_CU = '/Users/ddawson/NSF_V2_study/Disdrometer_study/CU_data/disdrometer_data/'
image_dir_CU = dis_dir_CU+'images/'
# 
# # Glen Deployment 1
# dis_dir_list = [dis_dir_glen]
# image_dir_list = [image_dir_glen]
# 
# dis_list = ['2010_06_13_D1_DIS_P1.nc']
# dis_ftype_list = ['Romine_netcdf']
# 
# dlocs = [(36.45612,-100.6039)]
# dis_name_list = ['2010_06_13_D1_P1'] # i.e. deployment #, probe #
# 
# # Comma-delimited thermodynamic and wind data files
# 
# tprh_filenames = [None]
# wind_filenames = [None]

# CU Deployment 2
dis_dir_list = [dis_dir_CU,dis_dir_CU,dis_dir_CU,dis_dir_CU,dis_dir_CU,dis_dir_CU]
image_dir_list = [image_dir_CU,image_dir_CU,image_dir_CU,image_dir_CU,image_dir_CU,image_dir_CU]

dis_list = ['20100519_UF06_PSD.dat','20100519_UF04_PSD.dat','20100519_UF07_PSD.dat','20100519_UF05_PSD.dat','20100519_UF01_PSD.dat','20100519_CU01_PSD.dat']
dis_ftype_list = ['CU_txt_2010','CU_txt_2010','CU_txt_2010','CU_txt_2010','CU_txt_2010','CU_txt_2010']

dlocs = [(35.86570,-98.96898),(35.91366,-98.97351),(35.95699,-98.97334),
         (36.00051,-98.96441),(36.05960,-98.97359),(36.02948,-98.96420)]
dis_name_list = ['UF6','UF4','UF7','UF5','SD1','CU1']

# Comma-delimited thermodynamic and wind data files 
tprh_filenames = [None,None,None,None,None,None]
wind_filenames = [None,None,None,None,None,None]

# Directory containing NEXRAD CFRadial files
nexrad_dir = '/Users/ddawson/NSF_V2_study/Disdrometer_study/NEXRAD_CFradial/KTLX/20100519/'

# Start and endtimes of timeseries plots (CST) NOTE: These are valid at the start of each 1-min DSD interval
# The times recorded in the files are valid at the *end* of the intervals
# IMPORTANT NOTE: For now, these times *MUST* correspond to times in both the thermodynamic and disdrometer 
# data files.  That is, they *CANNOT* lie outside the window of times.  In the future, I'll make the
# reading and time-matching routines more robust.
starttimes = [date2num(datetime(2010,5,19,14,29,00)),date2num(datetime(2010,5,19,14,33,00)),
              date2num(datetime(2010,5,19,14,38,00)),date2num(datetime(2010,5,19,14,42,00)),
              date2num(datetime(2010,5,19,14,49,00)),date2num(datetime(2010,5,19,14,53,00))]
endtimes =   [date2num(datetime(2010,5,19,15,13,00)),date2num(datetime(2010,5,19,15,22,00)),
              date2num(datetime(2010,5,19,15,31,00)),date2num(datetime(2010,5,19,15,37,00)),
              date2num(datetime(2010,5,19,15,45,00)),date2num(datetime(2010,5,19,15,40,00))]

# Parse command line argument and read in disdrometer/radar information from text file
# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
# version first, and if not, reading the textfile and pickling the read-in variables for future read-ins.  Think about this.

with open(sys.argv[1],'r') as inputfile:
    # Read and discard header
    inputfile.readline()
    inputfile.readline()
    inputfile.readline()
    dis_dir_glen = inputfile.readline().strip()
    image_dir_glen = inputfile.readline().strip()
    inputfile.readline()
    dis_dir_CU = inputfile.readline().strip()
    image_dir_CU = inputfile.readline().strip()
    inputfile.readline()
    numdis = N.int(inputfile.readline().strip().split()[0]) # Read in number of disdrometers
    
    # Read in disdrometer locations and names
    dis_list = []
    dis_dir_list = []
    image_dir_list = []
    dis_ftype_list = []
    tprh_filenames = []
    wind_filenames = []
    dlocs=[]
    dis_name_list = []
    starttimes = []
    endtimes =[]
    for l in xrange(numdis):
        line = inputfile.readline().strip().split(',')
        dloc = (N.float(line[0]),N.float(line[1])) # Tuple containing lat/lon
        dname = line[2] # Disdrometer name
        dfile = line[3] # Disdrometer filename
        dftype = line[4] # Disdrometer file type
        if(dftype[0:2] == 'CU'):
            dis_dir_list.append(dis_dir_CU)
            image_dir_list.append(image_dir_CU)
        else:
            dis_dir_list.append(dis_dir_glen)
            image_dir_list.append(image_dir_glen)
        tprh_file = line[5] # Thermo data file
        wind_file = line[6] # Wind data file
        st = map(int,line[7:13])
        et = map(int,line[13:])
        starttime = date2num(datetime(st[0],st[1],st[2],st[3],st[4],st[5]))
        endtime = date2num(datetime(et[0],et[1],et[2],et[3],et[4],et[5]))
        # Append to appropriate lists
        dlocs.append(dloc)
        dis_list.append(dfile)
        dis_name_list.append(dname)
        dis_ftype_list.append(dftype)
        tprh_filenames.append(tprh_file)
        wind_filenames.append(wind_file)
        starttimes.append(starttime)
        endtimes.append(endtime)
    
    inputfile.readline()
    nexrad_dir = inputfile.readline().strip()
print len(dis_dir_list)
print len(image_dir_list)
print len(dlocs)
print len(dis_name_list)
print len(dis_ftype_list)
print len(tprh_filenames)
print len(wind_filenames)
print len(starttimes)
print len(endtimes)

# ------
# Grab radar data for comparison if desired
if(comp_radar):
    if(not loadradopt):
        nexrad_filelist = glob.glob(nexrad_dir+'cfrad*.nc')
        
        raddates = []
        dBZ_D_tlist = []
        
        for path in nexrad_filelist:
        
            print "Opening file: ",path
            sweepfile_netcdf = Nio.open_file(path)
    
            # Grab the time information from the file
            
            sweeptime = sweepfile_netcdf.variables['time_coverage_start'].get_value().tostring()
            #''.join(sweeptime)
            print sweeptime
            
            radyear = int(sweeptime[:4])
            radmonth = int(sweeptime[5:7])
            radday = int(sweeptime[8:10])
            radhour = int(sweeptime[11:13])
            radmin = int(sweeptime[14:16])
            radsec = int(sweeptime[17:19])
            
            radtime = datetime(radyear,radmonth,radday,radhour,radmin,radsec)-timedelta(hours=6) # Get it to CST time
            raddates.append(radtime)
            
            sweepfile_netcdf.close()
            
            print "Time of sweep = ",radtime.strftime(fmt)
            
            # Now compute reflectivity at disdrometer location(s)
            
            dBZ_D_arr = dis.rad2DD(path,dlocs)
            
            dBZ_D_tlist.append(dBZ_D_arr)
        
        dBZ_D_tarr = N.array(dBZ_D_tlist)
        
        if(saveradopt):
            raddate_file=open('raddates','w')
            #dBZ_D_tlist_file=open('dBZ_D_tlist','w')
            pickle.dump(raddates,raddate_file)
            #pickle.dump(dBZ_D_tlist,dBZ_D_tlist_file)  # This apparently doesn't work... bug in numpy?
            savevars={}
            savevars['dBZ_D_tarr'] = dBZ_D_tarr
            N.savez('dBZ_D_tarr',**savevars)
        
    if(loadradopt):
        raddate_file=open('raddates','r')
        #dBZ_D_tlist_file=open('dBZ_D_tlist','r')
        raddates = pickle.load(raddate_file)
        #dBZ_D_tlist = pickle.load(dBZ_D_tlist_file)
        dBZ_D_tarr_file = N.load('dBZ_D_tarr.npz')
        dBZ_D_tarr = dBZ_D_tarr_file['dBZ_D_tarr']

    # Quick plot to test reflectivity extraction
    
    #fig1 = plt.figure()
    #ax1 = fig1.add_subplot(111)
    #ax1.plot_date(date2num(raddates),dBZ_D_plt,ls='-')
    #plt.title('Observed reflectivity at disdrometer')

# Outer disdrometer (and deployment) loop

for index,dis_dir,image_dir,dis_filename,dis_ftype,tprh_filename,wind_filename,dis_name,starttime,endtime in \
    zip(xrange(0,len(dis_list)),dis_dir_list,image_dir_list,dis_list,dis_ftype_list,tprh_filenames,wind_filenames,dis_name_list,starttimes,endtimes):


    print "Here?"

    if(timetospace):
        centertime = centertimes[index]

    if (not os.path.exists(image_dir)):
        os.makedirs(image_dir)

    # First, read in the thermodynamic data file
    if(plot_thermo):
        if(dis_ftype == 'Romine_netcdf'):
            tprh_file_name = dis_dir+'/'+tprh_filename
            thermodates,temps,rhs,pressures = dis.readthermoGlen(tprh_file_name)
        
            # Also read in wind data, which is in a separate file
            if(wind_filename != 'None'):
                wind_file_name = dis_dir+'/'+wind_filename
                winddates,winds,winddirs,windgusts = dis.readwindGlen(wind_file_name)
        
        elif(dis_ftype == 'CU_txt'):
        
            tprh_file_name = dis_dir+'/'+tprh_filename
            thermodates,temps,rhs,pressures,winds,winddirs,sws,swd = dis.readthermoCU(tprh_file_name)
        
            winddates=thermodates
        
            # Convert to CST from UTC
            for i,date in enumerate(thermodates):
                date=date-timedelta(hours=6)
                thermodates[i] = date
                
    
    # Find index in thermo data file corresponding to chosen start and end times
    
    
    if(plot_thermo):
        thermo_startindex = 0
        thermo_endindex = len(thermodates)
        
        for i,date in enumerate(thermodates):
            datenum = date2num(date)
            #print "starttime,datenum = ",starttime,datenum
            if(datenum == starttime):
                thermo_startindex = i
            if(datenum == endtime):
                thermo_endindex = i
        
        #print thermo_startindex,thermo_endindex
     
        # Do the same for the wind data file (for Glen's data only)
        
        if(wind_filename != 'None'):
            wind_startindex = 0
            wind_endindex = len(winddates)
            
            for i,date in enumerate(winddates):
                datenum = date2num(date)
                #print "starttime,datenum = ",starttime,datenum
                if(datenum == starttime):
                    wind_startindex = i
                if(datenum == endtime):
                    wind_endindex = i
        else:
            wind_startindex = thermo_startindex
            wind_endindex = thermo_endindex
        
        #print wind_startindex,wind_endindex
     
        #N.set_printoptions(threshold='nan')
        #print sws[thermo_startindex:thermo_endindex]
 
    # Now read in the disdrometer data file
    
    dis_file_name = dis_dir+'/'+dis_filename
    
    if(dis_ftype == 'Romine_netcdf'):
                
        pcount,counts_bin,Nc_bin,logNc_bin,min_size,avg_size,max_size,bin_width,disdates = \
            dis.readdataGlennetCDF(starttime,endtime,dis_file_name)
        
        pcount2 = pcount
        
        if(windQC):
            windmask = N.ones((N.size(avg_size)))
        
            # Compute average wind speed for each 1-min DSD interval
            winds_1minavg = N.zeros_like(disdates)
            #print winds_1minavg.shape
            #print winds.shape
            for i in xrange(N.size(winds_1minavg)-1):
                #print i*6,i*6+6
                winds_1minavg[i+1] = N.mean(winds[i*6:i*6+6])   
            
            for i in range(N.size(winds_1minavg)):
                if(winds_1minavg[i] > wind_thresh):
                    Nc_bin[:,i] = N.ma.masked_array(Nc_bin[:,i],mask=windmask)
                    logNc_bin[:,i] = N.ma.masked_array(logNc_bin[:,i],mask=windmask)
                
    elif(dis_ftype == 'CU_txt'):
        
        dates,times,intensities,preciptots,weathercodes,reflectivities,visibilities,pcount,pcount2, \
               sensortemps,min_size,max_size,avg_size,Nc_bin = dis.readdataCU(dis_file_name)
        
        year_end   = [int(date[6:10]) for date in dates]
        month_end  = [int(date[3:5]) for date in dates]
        day_end    = [int(date[:2]) for date in dates]
        hour_end   = [int(time[:2]) for time in times]
        min_end    = [int(time[3:5]) for time in times]
        sec_end    = [int(time[6:8]) for time in times]
        
        # Construct the datetime objects
        
        disdate_end = []
        dis_startindex = 0
        dis_endindex = len(year_end)
        j = 0        
        Nc_bin_temp = Nc_bin
        pcount_temp = pcount
        pcount2_temp = pcount2
        for i in range(len(year_end)):
            distimeend = datetime(year_end[i],month_end[i],day_end[i],hour_end[i],min_end[i],sec_end[i])-timedelta(hours=6) # Get to CST time
            # There appear to be duplicate time records in some files, so check for them and throw them out
            # Although some of them seem to contain different data.  I need to check with Katja about this.
            if(i < len(year_end)-1):
                distimenext = datetime(year_end[i+1],month_end[i+1],day_end[i+1],hour_end[i+1],min_end[i+1],sec_end[i+1])-timedelta(hours=6) # Get to CST time
            else:
                distimenext = None
            if(distimeend != distimenext):
                # Find index corresponding to chosen start and endtimes
                if(date2num(distimeend) == starttime):
                    dis_startindex = j
                if(date2num(distimeend) == endtime):
                    dis_endindex = j
                disdate_end.append(distimeend)
                j = j + 1
            else:
                # Remove duplicate record from number concentration and pcount lists
                Nc_bin_temp = N.delete(Nc_bin,[i+1],axis=0)
                pcount_temp = N.delete(pcount,[i+1])
                pcount2_temp = N.delete(pcount2,[i+1])
        
        Nc_bin = Nc_bin_temp
        pcount = pcount_temp
        pcount2 = pcount2_temp
        
        disdates = disdate_end[dis_startindex:dis_endindex+1]
        
        #print "CU disdrometer times: ",disdates
        
        bin_width = max_size-min_size
        
        # Enhance small drops for testing if desired
        if(enhancesmalldrops):
            Nc_bin[:,0:7] = Nc_bin[:,0:7]*enhancefactor
        
        Nc_bin = Nc_bin.swapaxes(0,1)
        logNc_bin = N.ma.log10(Nc_bin)
        
        #print Nc_bin
        #print logNc_bin
        
        #Restrict time dimension to start and end times
        Nc_bin = Nc_bin[:,dis_startindex:dis_endindex+1]
        logNc_bin = logNc_bin[:,dis_startindex:dis_endindex+1]
        pcount = pcount[dis_startindex:dis_endindex+1]
        pcount2 = pcount2[dis_startindex:dis_endindex+1]
        
        #Nc_bin = N.ma.masked_array(Nc_bin,mask=N.where(Nc_bin == -999., True, False))
        #logNc_bin = N.ma.masked_array(logNc_bin,mask=N.where(Nc_bin == -999., True, False))
        logNc_bin = N.ma.masked_where(Nc_bin <= 0.0, logNc_bin)
        
        #print Nc_bin.shape
    
    elif(dis_ftype == 'CU_txt_2010'):
        
        dates,times,intensities,preciptots,weathercodes,reflectivities,visibilities,pcount,pcount2, \
               sensortemps,min_size,max_size,avg_size,Nc_bin = dis.readdataCU_2010(dis_file_name)
        
        year_end   = [int(date[6:10]) for date in dates]
        month_end  = [int(date[3:5]) for date in dates]
        day_end    = [int(date[:2]) for date in dates]
        hour_end   = [int(time[:2]) for time in times]
        min_end    = [int(time[3:5]) for time in times]
        sec_end    = [int(time[6:8]) for time in times]
        
        # Construct the datetime objects
        
        disdate_end = []
        dis_startindex = 0
        dis_endindex = len(year_end)
        j = 0        
        Nc_bin_temp = Nc_bin
        pcount_temp = pcount
        pcount2_temp = pcount2
        for i in range(len(year_end)):
            distimeend = datetime(year_end[i],month_end[i],day_end[i],hour_end[i],min_end[i],sec_end[i])-timedelta(hours=6) # Get to CST time
            # There appear to be duplicate time records in some files, so check for them and throw them out
            # Although some of them seem to contain different data.  I need to check with Katja about this.
            if(i < len(year_end)-1):
                distimenext = datetime(year_end[i+1],month_end[i+1],day_end[i+1],hour_end[i+1],min_end[i+1],sec_end[i+1])-timedelta(hours=6) # Get to CST time
            else:
                distimenext = None
            if(distimeend != distimenext):
                # Find index corresponding to chosen start and endtimes
                if(date2num(distimeend) == starttime):
                    dis_startindex = j
                if(date2num(distimeend) == endtime):
                    dis_endindex = j
                disdate_end.append(distimeend)
                j = j + 1
            else:
                # Remove duplicate record from number concentration and pcount lists
                Nc_bin_temp = N.delete(Nc_bin,[i+1],axis=0)
                pcount_temp = N.delete(pcount,[i+1])
                pcount2_temp = N.delete(pcount2,[i+1])
                
        
        Nc_bin = Nc_bin_temp
        pcount = pcount_temp
        pcount2 = pcount2_temp
        
        disdates = disdate_end[dis_startindex:dis_endindex+1]
        
        #print "CU disdrometer times: ",disdates
        
        bin_width = max_size-min_size
        
        # Enhance small drops for testing if desired
        if(enhancesmalldrops):
            Nc_bin[:,0:7] = Nc_bin[:,0:7]*enhancefactor
        
        Nc_bin = Nc_bin.swapaxes(0,1)
        logNc_bin = N.ma.log10(Nc_bin)
        
        #print Nc_bin
        #print logNc_bin
        
        #Restrict time dimension to start and end times
        Nc_bin = Nc_bin[:,dis_startindex:dis_endindex+1]
        logNc_bin = logNc_bin[:,dis_startindex:dis_endindex+1]
        pcount = pcount[dis_startindex:dis_endindex+1]
        pcount2 = pcount2[dis_startindex:dis_endindex+1]
        
        #Nc_bin = N.ma.masked_array(Nc_bin,mask=N.where(Nc_bin == -999., True, False))
        #logNc_bin = N.ma.masked_array(logNc_bin,mask=N.where(Nc_bin == -999., True, False))
        logNc_bin = N.ma.masked_where(Nc_bin <= 0.0, logNc_bin)
    
    if(dis_ftype == 'CU_txt_2010'):
        disdates_start = [x-timedelta10s for x in disdates]
        disdates_start.append(disdates_start[-1]+timedelta10s)    # Add an extra 10 sec for the last time bin boundary
        disdates_avg = [x-timedelta5s for x in disdates]
    else:
        disdates_start = [x-timedelta1min for x in disdates]
        disdates_start.append(disdates_start[-1]+timedelta1min)    # Add an extra 60 sec for the last time bin boundary
        disdates_avg = [x-timedelta30s for x in disdates]
 
    # Thin the thermodynamic data (which is output every 10 or 1 s) to every 1 min
    
    if(dis_ftype == 'Romine_netcdf'):
        stride = 6
    elif(dis_ftype == 'CU_txt'):
        stride = 60
    if(plot_thermo):
        thermodates_1min = thermodates[thermo_startindex:thermo_endindex+1:stride]
        winddates_1min = winddates[wind_startindex:wind_startindex+1:stride]
        pressures_1min = pressures[thermo_startindex:thermo_endindex+1:stride]
        temps_1min = temps[thermo_startindex:thermo_endindex+1:stride]
        rhs_1min = rhs[thermo_startindex:thermo_endindex+1:stride]
    
        # Sometimes there are missing records in the disdrometer data file (only found one affected file so far)
        # Deal with this by removing the thermodynamic records that correspond with the missing disdrometer records
        
        for i,thermodate in enumerate(thermodates_1min[:]):      # Need a copy here because thermodates_1min may be modified
            count = disdates.count(thermodate)
            if(count == 0):
                thermodates_1min.remove(thermodate)
                pressures_1min = N.delete(pressures_1min,[i])
                temps_1min = N.delete(temps_1min,[i])
                rhs_1min = N.delete(rhs_1min,[i])
    
    # For time-to-space conversion, compute the positions corresponding to the times above
    
    if(timetospace and comp_radar):
        delta_t_rad = [(x-centertime).total_seconds() for x in raddates]
        delta_t_dis = [(x-centertime).total_seconds() for x in disdates]
        delta_t_dis_start = [(x-centertime).total_seconds() for x in disdates_start]
        delta_t_dis_avg = [(x-centertime).total_seconds() for x in disdates_avg]
        if(plot_thermo):
            delta_t_thermo = [(x-centertime).total_seconds() for x in thermodates]
            delta_t_thermo_1min = [(x-centertime).total_seconds() for x in thermodates_1min]
            delta_t_wind = [(x-centertime).total_seconds() for x in winddates]
            delta_t_wind_1min = [(x-centertime).total_seconds() for x in winddates_1min]
        
        xloc_rad = [(x*stormmotion)/1000.0 for x in delta_t_rad]
        xloc_dis = [(x*stormmotion)/1000.0 for x in delta_t_dis]
        xloc_dis_start = [(x*stormmotion)/1000.0 for x in delta_t_dis_start]
        xloc_dis_avg = [(x*stormmotion)/1000.0 for x in delta_t_dis_avg]
        if(plot_thermo):
            xloc_thermo = [(x*stormmotion)/1000.0 for x in delta_t_thermo]
            xloc_thermo_1min = [(x*stormmotion)/1000.0 for x in delta_t_thermo_1min]
            xloc_wind = [(x*stormmotion)/1000.0 for x in delta_t_wind]
            xloc_wind_1min = [(x*stormmotion)/1000.0 for x in delta_t_wind_1min]
        
        plotx_rad = N.array(xloc_rad)
        plotx_dis = N.array(xloc_dis)
        plotx_dis_start = N.array(xloc_dis_start)
        plotx_dis_avg = N.array(xloc_dis_avg)
        if(plot_thermo):
            plotx_thermo = N.array(xloc_thermo)
            plotx_thermo_1min = N.array(xloc_thermo_1min)
            plotx_wind = N.array(xloc_wind)
            plotx_wind_1min = N.array(xloc_wind_1min)
    else:
        if(comp_radar):
            plotx_rad = date2num(raddates)
        plotx_dis = date2num(disdates)
        plotx_dis_start = date2num(disdates_start)
        plotx_dis_avg = date2num(disdates_avg)
        if(plot_thermo):
            plotx_thermo = date2num(thermodates)
            plotx_thermo_1min = date2num(thermodates_1min)
            plotx_wind = date2num(winddates)
            plotx_wind_1min = date2num(winddates_1min)
    
    # Next, calculate air density from temperature and pressure (read in above)
    
    if(plot_thermo):
        P_Pa = pressures_1min*100.0
        T_K = temps_1min+273.16
        rho = P_Pa/(Rd*T_K)
    else:
        # For now just assume a density of 1.204 kg/m^3 if we don't have the info
        rho = N.empty_like(disdates,dtype=N.float)
        rho[:] = 1.204
    
    # Now for the fun part.  Calculate exponential and gamma distribution size distribution parameters
    # using the method of moments, after Zhang et al. 2008 and Tokay and Short 1996
    
    synthbins,exp_DSD,gam_DSD,dis_DSD = dis.calc_DSD(min_size,avg_size,max_size,bin_width,Nc_bin,logNc_bin,rho,qrQC,qr_thresh)
    
    # Unpack needed values from returned tuples
    
    N_expDSD,N0_exp,lamda_exp,mu_exp,qr_exp,Ntr_exp,refl_DSD_exp,D_med_exp,D_m_exp = exp_DSD
    N_gamDSD,N0_gam,lamda_gam,mu_gam,qr_gam,Ntr_gam,refl_DSD_gam,D_med_gam,D_m_gam = gam_DSD
    Nc_bin,logNc_bin,D_med_disd,D_m_disd,D_mv_disd,D_ref_disd,QR_disd,refl_disd = dis_DSD
        
    if(avg_data):
        # Also calculate quantities for the average DSD
        pcount_sum = pcount.sum()
        pcount2_sum = pcount2.sum()
        Nc_bin_avg = N.ma.sum(Nc_bin,axis=1)/N.ma.count(Nc_bin,axis=1)
        
        #print N.ma.sum(Nc_bin,axis=1)
        print N.ma.count(Nc_bin,axis=1)
        logNc_bin_avg = N.ma.log10(Nc_bin_avg)
        logNc_bin_avg = N.ma.masked_where(Nc_bin_avg <= 0.0, logNc_bin_avg)
        #print Nc_bin_avg
        #print logNc_bin_avg
    
        rho_avg = N.mean(rho)
        #print "rho_avg = ",rho_avg
        dummy,exp_DSD_avg,gam_DSD_avg,dis_DSD_avg = dis.calc_DSD(min_size,avg_size,max_size,bin_width,Nc_bin_avg,logNc_bin_avg,rho_avg,False,qr_thresh)
        
        N_expDSD_avg,N0_exp_avg,lamda_exp_avg,mu_exp_avg,qr_exp_avg,Ntr_exp_avg,refl_DSD_exp_avg,D_med_exp_avg,D_m_exp_avg = exp_DSD_avg
        N_gamDSD_avg,N0_gam_avg,lamda_gam_avg,mu_gam_avg,qr_gam_avg,Ntr_gam_avg,refl_DSD_gam_avg,D_med_gam_avg,D_m_gam_avg = gam_DSD_avg
        Nc_bin_avg,logNc_bin_avg,D_med_disd_avg,D_m_disd_avg,D_mv_disd_avg,D_ref_disd_avg,QR_disd_avg,refl_disd_avg = dis_DSD_avg
    
    if(calc_evap):
    
        # Ok, now that we have the exponential and gamma distribution fits, we can do some other cool stuff
        # like derive bulk evaporative cooling rates, and compare with the observed cooling rate.
        # Let's do that now, shall we?

        # Evaporation rate based on exponential and gamma distribution (may also want to try to compute
        # evaporation directly from bins, which will require a separate calculation).
        
        QREVP_exp,COOL_RATE_exp = dsd.calc_evap(rho,T_K,P_Pa,rhs_1min,N0_exp,lamda_exp,mu_exp)
        QREVP_gam,COOL_RATE_gam = dsd.calc_evap(rho,T_K,P_Pa,rhs_1min,N0_gam,lamda_gam,mu_gam)
            
        # Now we can compare this with the observed cooling rate
    
        # Compute 1-min averages of original temperature data
        
        if(dis_ftype == 'Romine_netcdf'):
            filter_length = 7
            timestep = 10.0
        else:
            filter_length = 61
            timestep = 1.0
    
        temps_1minavg = ndimage.uniform_filter1d(temps,filter_length)
        temps_1minavg = temps_1minavg[thermo_startindex:thermo_endindex+1:stride]
    
        #OBS_COOL_RATE = (temps_1minavg[1:]-temps_1minavg[:-1])/timestep       # Centered time difference on original temp data
        #OBS_COOL_RATE = N.insert(OBS_COOL_RATE,[0],OBS_COOL_RATE[0])    # Just fill in the value at the beginning and end
        #OBS_COOL_RATE = N.append(OBS_COOL_RATE,OBS_COOL_RATE[-1])
        
        #OBS_COOL_RATE = ndimage.uniform_filter1d(OBS_COOL_RATE,filter_length)
        
        OBS_COOL_RATE = N.gradient(temps_1minavg,1.0)
        
        #OBS_COOL_RATE_1min = OBS_COOL_RATE[thermo_startindex:thermo_endindex+1:stride]
        if(plot_opt):
            figevap = plt.figure(figsize=(8,6))
            axtemp = figevap.add_subplot(311)
            axrh = axtemp.twinx()
            #axtemp.plot(date2num(thermodates),temps,ls='-',c='r',marker='None',label='Temp')
            axtemp.plot(plotx_thermo,temps,ls='-',c='r',marker='None',label='Temp')
            axrh.plot(plotx_thermo,rhs,ls='-',c='b',marker='None',label='RH (%)')
            if(not timetospace):
                axtemp.xaxis.set_major_locator(MinuteLocator(interval=5))
                axtemp.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                axtemp.set_xlim(min(starttimes),max(endtimes))
                if(reverse_times):
                    axtemp.invert_xaxis()
            else:
                axtemp.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
                axtemp.set_xlim(xstart,xstop)
                axtemp.invert_xaxis()
            axtemp.set_ylim(10.0,30.0)
            axrh.set_ylim(0.0,100.0)
            axtemp.set_ylabel(r'Temperature ($^\circ$C)')
            #plt.legend(loc='best')
            
            axevap = figevap.add_subplot(312)
            axevap.plot(plotx_thermo[thermo_startindex:thermo_endindex+1:stride],OBS_COOL_RATE,ls='-',c='k',marker='None',label='LCR (obs)')
            axevap.plot(plotx_dis,COOL_RATE_exp*60.0,ls='-',c='b',label='LCR (exp)')
            axevap.plot(plotx_dis,COOL_RATE_gam*60.0,ls='-',c='g',label='LCR (gam)')
            if(not timetospace):
                axevap.xaxis.set_major_locator(MinuteLocator(interval=5))
                axevap.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                axevap.set_xlim(min(starttimes),max(endtimes))
                if(reverse_times):
                    axevap.invert_xaxis()
            else:
                axevap.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
                axevap.set_xlim(xstart,xstop)
                axevap.invert_xaxis()
            #axevap.set_ylim(-0.002,0.001)
            axevap.set_ylabel(r'Cooling rate ($K min^{-1}$)')
            
            axwind = figevap.add_subplot(313)
            axwind.plot(plotx_wind,winds,ls='-',c='k',marker='None',label='wind speed')
            if(dis_ftype == 'Romine_netcdf'):
                axwind.plot(plotx_wind,windgusts,ls='--',c='k',marker='None',label='wind gust')
            if(windQC):
                #print plotx_dis.shape,winds_1minavg.shape
                axwind.plot(plotx_dis,winds_1minavg,ls='--',c='k',marker='None',label='wind speed (1-min avg)')
            if(not timetospace):
                axwind.xaxis.set_major_locator(MinuteLocator(interval=5))
                axwind.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                axwind.set_xlim(min(starttimes),max(endtimes))
                if(reverse_times):
                    axwind.invert_xaxis()
            else:
                axwind.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
                axwind.set_xlim(xstart,xstop)
                axwind.invert_xaxis()
            axwind.set_ylim(0.0,20.0)
            axwind.yaxis.set_major_locator(MultipleLocator(base=5.0))
            axwind.set_ylabel(r'Wind speed ($m s^{-1}$)')
        
            plt.setp(axtemp.get_xticklabels(),visible=False)
            plt.setp(axrh.get_xticklabels(),visible=False)
            plt.setp(axevap.get_xticklabels(),visible=False)
            if(not timetospace):
                figevap.autofmt_xdate()
            plt.savefig(image_dir+dis_name+'_evap.eps',dpi=200)
        
    if (plot_opt):
        # Plot the log of number concentration for each size bin
        # Compute the DSD for a single drop in each bin, for comparison
        counts_1drop = N.ones_like(avg_size)
        vt = dis.assignfallspeed(avg_size)
        if(dis_ftype == 'CU_txt_2010'):
            sampling_period = 10.
        else:
            sampling_period = 60.
        Nc_bin_1drop = counts_1drop/(vt*sampling_period*dis.eff_sensor_area*(max_size-min_size))
        if (plot_DSDs):
            for t in range(len(disdates)):
                #print "Here! disdates(t)",disdates[t],N.ma.sum(Nc_bin[:,t])
                #print Nc_bin[:,t]
                sum = N.ma.sum(Nc_bin[:,t])
                # Only plot if there are actual data
                if(sum != 0.0):
                    fig1=plt.figure(figsize=(8,6))
                    ax1=fig1.add_subplot(111)
                    if(dis_ftype == 'CU_txt_2010'):
                        plt.title('Disdrometer 10-sec DSD Fits for Time '+disdates[t].strftime(fmt2)+' CST')
                    else:
                        plt.title('Disdrometer 1-min DSD Fits for Time '+disdates[t].strftime(fmt2)+' CST')
                    ax1.bar(min_size,Nc_bin[:,t]*1000.0,max_size-min_size,log=True,color='tan')
                    ax1.bar(min_size,Nc_bin_1drop*1000.0,max_size-min_size,log=True,color='None',edgecolor='k')
                    # Plot derived exponential DSD
                    ax1.plot(synthbins,N_expDSD[t],lw=2)
                    #ax1.plot(synthbins,N_gamDSD[t],lw=2)
                    ax1.plot(avg_size,N_gamDSD[t],lw=2) # Use disdrometer bins
                    ax1.set_yscale('log')
                    ax1.set_ylim(10.**2.0,10.**8.5)
                    ax1.set_xlim(0.0,9.0)
                    #ax1.set_xlim(0.0,16.0)
                    ax1.set_xlabel('Drop Diameter (mm)')
                    ax1.set_ylabel(r'N(D) $(m^{-4})$')
                    N0_disp = N0_exp[t]/1.0e6
                    ax1.text(0.50,0.8,r'$N_0$ (exp) = %2.2f'%N0_disp+r' x $10^6 m^{-4}$',transform=ax1.transAxes)
                    ax1.text(0.50,0.75,'Shape parameter (gamma) = %2.2f'%mu_gam[t],transform=ax1.transAxes)
                    ax1.text(0.50,0.7,r'$D_{m43}$ = %2.2f'%D_m_gam[t],transform=ax1.transAxes)
                    ax1.text(0.50,0.65,'Mixing ratio = %2.2f'%QR_disd[t]+' g/kg',transform=ax1.transAxes)
                    ax1.text(0.50,0.60,'Reflectivity =%2.2f'%refl_disd[t]+' dbZ',transform=ax1.transAxes)
                    ax1.text(0.50,0.55,'Particle count = '+str(pcount[t]),transform=ax1.transAxes)
                    ax1.text(0.50,0.50,'Particle count (QC) = '+str(pcount2[t]),transform=ax1.transAxes)
#                     for bin in range(N.size(avg_size)):
#                         if(counts_bin[bin,t] >= 0.0):
#                             ax1.annotate(N.int(counts_bin[bin,t]),(avg_size[bin],10.**8.))
                    #ax1.xaxis.set_major_locator(MinuteLocator(byminute=range(0,60,2)))
                    #ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                    plt.savefig(image_dir+dis_name+'_t'+str(t)+'DSD_plot.eps',dpi=200,bbox_inches='tight')
        
        if (avg_data):
        
            # Compute averaged DSD assuming a single drop in each bin, for comparison
            Nc_bin_1drop_avg = Nc_bin_1drop/N.max(N.ma.count(Nc_bin,axis=1))
        
            fig1=plt.figure(figsize=(8,6))
            ax1=fig1.add_subplot(111)
            plt.title('Disdrometer full-deployment DSD')
            ax1.bar(min_size,Nc_bin_avg*1000.0,max_size-min_size,log=True,color='tan')
            ax1.bar(min_size,Nc_bin_1drop_avg*1000.0,max_size-min_size,log=True,color='None',edgecolor='k')
            # Plot derived exponential DSD
            ax1.plot(synthbins,N_expDSD_avg[0],lw=2)
            #ax1.plot(synthbins,N_gamDSD_avg[0],lw=2)
            ax1.plot(avg_size,N_gamDSD_avg[0],lw=2)
            ax1.set_yscale('log')
            ax1.set_ylim(10.**2.0,10.**8.5)
            ax1.set_xlim(0.0,9.0)
            #ax1.set_xlim(0.0,16.0)
            ax1.set_xlabel('Drop Diameter (mm)')
            ax1.set_ylabel(r'N(D) $(m^{-4})$')
            N0_disp = N0_exp_avg[0]/1.0e6
            ax1.text(0.50,0.8,r'$N_0$ (exp) = %2.2f'%N0_disp+r' x $10^6 m^{-4}$',transform=ax1.transAxes)
            ax1.text(0.50,0.75,'Shape parameter (gamma) = %2.2f'%mu_gam_avg[0],transform=ax1.transAxes)
            ax1.text(0.50,0.7,r'$D_{m43}$ = %2.2f'%D_m_gam_avg[0],transform=ax1.transAxes)
            ax1.text(0.50,0.65,'Mixing ratio = %2.2f'%QR_disd_avg[0]+' g/kg',transform=ax1.transAxes)
            ax1.text(0.50,0.60,'Reflectivity =%2.2f'%refl_disd_avg[0]+' dbZ',transform=ax1.transAxes)
            ax1.text(0.50,0.55,'Particle count (QC) = '+str(pcount2_sum),transform=ax1.transAxes)
            ax1.text(0.50,0.50,'Duration of DSD = '+str(N.max(N.ma.count(Nc_bin,axis=1)))+' min',transform=ax1.transAxes)
            #ax1.xaxis.set_major_locator(MinuteLocator(byminute=range(0,60,2)))
            #ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            plt.savefig(image_dir+dis_name+'_DSD_avg.eps',dpi=200,bbox_inches='tight')
                    
        
        # Plot a timeseries of reflectivity, shape parameter, and intercept parameter (exponential)
        
        if(comp_radar):
            #dBZ_D_plt = [x[index] for x in dBZ_D_tlist]
            dBZ_D_plt = dBZ_D_tarr[:,index]
            #print dBZ_D_plt
        
        if(False):
            fontP = FontProperties()
            fontP.set_size('small')
            
            fig1 = plt.figure(figsize=(8,6))
            ax1 = fig1.add_subplot(411)
            ax1.plot_date(date2num(disdates),refl_disd,ls='-',c='b')
            #ax1.plot_date(date2num(disdates),refl_expDSD,ls='--')
            #ax1.plot_date(date2num(disdates),refl_gamDSD,ls=':')
            if(comp_radar):
                ax1.plot_date(date2num(raddates),dBZ_D_plt,ls='-',c='g')
            ax2 = fig1.add_subplot(412)
            ax2.plot_date(date2num(disdates),N0_exp/1.0e6,ls='-')
            ax3 = fig1.add_subplot(413)
            ax3.plot_date(date2num(disdates),mu_gam,ls='-')
            ax4 = fig1.add_subplot(414)
            ax4.plot_date(date2num(disdates),D_m_disd,ls='-',c='k',label = "D (M-weight mean)")
            ax4.plot_date(date2num(disdates),D_mv_disd,ls='-',c='b',label = "D (mean V)")
            ax4.plot_date(date2num(disdates),D_med_disd,ls='-',c='g',label = "D (med V)")
            ax4.plot_date(date2num(disdates),D_ref_disd,ls='-',c='purple',label = "D (Z-weight mean)")
            ax4.legend(bbox_to_anchor=(0.0,-0.25,1.0,-0.2),ncol=4,mode="expand",borderaxespad=0,prop = fontP)
            ax1.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax2.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax3.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax4.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax3.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax4.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            #plt.setp(ax1.get_xticklabels(),visible=False)
            #plt.setp(ax2.get_xticklabels(),visible=False)
            #plt.setp(ax3.get_xticklabels(),visible=False)
            ax1.set_xlim(starttime,endtime)
            ax1.set_ylim(0.0,80.0)
            ax1.set_ylabel('Reflectivity (dBZ)')
            ax2.set_xlim(starttime,endtime)
            ax2.set_ylim(0.0,15.0)
            ax2.set_ylabel(r'$N_0 x 10^6$ (exp.)')
            ax3.set_xlim(starttime,endtime)
            ax3.set_ylim(0.0,20.0)
            ax3.set_ylabel(r'$\mu$')
            ax4.set_xlim(starttime,endtime)
            ax4.set_ylim(0.0,10.0)
            ax4.set_ylabel('Diameter (mm)')
            
            if(not timetospace):
                fig1.autofmt_xdate()
            
            plt.savefig(image_dir+dis_name+'_dBZ_comp.eps',dpi=200)
        
        # Plot a timeseries pcolor plot of log of number concentration
        
        fig2 = plt.figure(figsize=(8,6))
        #ax2 = fig2.add_subplot(111)
        ax1 = host_subplot(111)
        if(not timetospace):
            fig2.autofmt_xdate()
        ax2 = ax1.twinx()
        divider = make_axes_locatable(ax1)
        
        #print "logNc_bin",logNc_bin
        # Times are valid at end of DSD intervals
        print "plotx_dis_start",plotx_dis_start
        print "plotx_dis_avg",plotx_dis_avg
        C = ax1.pcolor(plotx_dis_start,min_size,logNc_bin,vmin=-1.0,vmax=3.0)
        #Overlay mass-weighted mean diameter
        #ax1.plot(plotx_dis_avg,D_m_disd,ls=':',c='k',lw=3,label=r'$D_{mr43} (mm)$')
        #Overlay median volume diameter
        ax1.plot(plotx_dis_avg,D_med_disd,ls=':',c='k',lw=3,label=r'$D_{0} (mm)$')
        ax1.plot(plotx_dis_avg,D_med_gam,ls='--',c='k',lw=3,label=r'$D_{0,gam} (mm)$')
        if(not timetospace):
            ax1.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax1.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax1.invert_xaxis()
        else:
            ax1.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax1.set_xlim(xstart,xstop)
            ax1.invert_xaxis()
        ax1.yaxis.set_major_locator(MultipleLocator(base=1.0))
        ax1.set_ylim(0.0,9.0)
        ax1.set_ylabel('D (mm)')
        ax1.legend(loc='upper right')
        cax = divider.append_axes("bottom", size="5%", pad=0.35)
        cb = fig2.colorbar(C, cax=cax,orientation='horizontal')
        cb.set_label(r'log[N ($m^{-3} mm^{-1}$)]')
        
        #ax2.set_aspect(5000.0)
        #ax2.set_aspect(0.0001)
            
        #divider.set_aspect(True)
        
        ax2.plot(plotx_dis_avg,refl_disd,ls='-',c='r',marker='o',label=r'$Z_{dis} (dBZ)$')
        if(comp_radar):
            ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',marker='o',label=r'$Z_{rad} (dBZ)$')
        if(not timetospace):
            ax2.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax2.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax2.invert_xaxis()
        elif(comp_radar):
            ax2.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax2.set_xlim(xstart,xstop)
            ax2.invert_xaxis()
        ax2.yaxis.set_major_locator(MultipleLocator(base=5.0))
        ax2.set_ylim(-10.0,80.0)
        ax2.set_ylabel('Reflectivity (dBZ)')
        ax2.legend(loc='upper left')
    
#         pcountax = divider.append_axes("top",size="50%",pad=0.2,sharex=ax1)
#         pcountax.plot(plotx_dis_avg,pcount,ls='-',marker='o',c='k')
#         pcountax.plot(plotx_dis_avg,pcount2,ls='-',marker='o',c='b')
#         if(not timetospace):
#             pcountax.xaxis.set_major_locator(MinuteLocator(interval=5))
#             pcountax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
#             pcountax.set_xlim(starttime,endtime)
#             if(reverse_times):
#                 pcountax.invert_xaxis()
#         else:
#             pcountax.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
#             pcountax.set_xlim(xstart,xstop)
#             pcountax.invert_xaxis()
#         pcountax.yaxis.set_major_locator(MultipleLocator(base=5.0))
#         pcountax.set_yscale('log')
#         pcountax.set_ylim(1.0,6000.0)
#         pcountax.set_ylabel('# of particles')
#         pcountax.set_title('Disdrometer Log(N) and Reflectivity')
#         plt.xticks(visible=False)
        
        
        plt.savefig(image_dir+dis_name+'_logNc.eps',dpi=300)
        
        # Plot a timeseries pcolor plot of log of number concentration as above, but
        # using the gamma distribution fits instead of the raw data
        
        fig2 = plt.figure(figsize=(8,6))
        #ax2 = fig2.add_subplot(111)
        ax1 = host_subplot(111)
        if(not timetospace):
            fig2.autofmt_xdate()
        ax2 = ax1.twinx()
        divider = make_axes_locatable(ax1)
        
        logN_gamDSD = N.log10(N_gamDSD/1000.) # Get to log(m^-3 mm^-1)
       
        logN_gamDSD = logN_gamDSD.swapaxes(0,1)
        
        #print "logNc_bin",logNc_bin
        # Times are valid at end of DSD intervals
        C = ax1.pcolor(plotx_dis_start,min_size,logN_gamDSD,vmin=-1.0,vmax=3.0)
        #Overlay mass-weighted mean diameter
        #ax1.plot(plotx_dis_avg,D_m_gam,ls=':',c='k',lw=2,label=r'$D_{mr43} (mm)$')
        #Overlay median volume diameter
        ax1.plot(plotx_dis_avg,D_med_gam,ls=':',c='k',lw=3,label=r'$D_{0} (mm)$')
        if(not timetospace):
            ax1.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax1.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax1.invert_xaxis()
        else:
            ax1.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax1.set_xlim(xstart,xstop)
            ax1.invert_xaxis()
        ax1.yaxis.set_major_locator(MultipleLocator(base=1.0))
        ax1.set_ylim(0.0,9.0)
        ax1.set_ylabel('D (mm)')
        ax1.legend(loc='upper right')
        cax = divider.append_axes("bottom", size="5%", pad=0.35)
        cb = fig2.colorbar(C, cax=cax,orientation='horizontal')
        cb.set_label(r'log[N ($m^{-4} mm^{-1}$)]')
        
        #ax2.set_aspect(5000.0)
        #ax2.set_aspect(0.0001)
            
        #divider.set_aspect(True)
        
        ax2.plot(plotx_dis_avg,refl_DSD_gam,ls='-',c='r',marker='o',label=r'$Z_{dis} (dBZ)$')
        if(comp_radar):
            ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',marker='o',label=r'$Z_{rad} (dBZ)$')
        if(not timetospace):
            ax2.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax2.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax2.invert_xaxis()
        elif(comp_radar):
            ax2.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax2.set_xlim(xstart,xstop)
            ax2.invert_xaxis()
        ax2.yaxis.set_major_locator(MultipleLocator(base=5.0))
        ax2.set_ylim(-10.0,80.0)
        ax2.set_ylabel('Reflectivity (dBZ)')
        ax2.legend(loc='upper left')
    
#         pcountax = divider.append_axes("top",size="50%",pad=0.2,sharex=ax1)
#         pcountax.plot(plotx_dis_avg,pcount,ls='-',marker='o',c='k')
#         pcountax.plot(plotx_dis_avg,pcount2,ls='-',marker='o',c='b')
#         if(not timetospace):
#             pcountax.xaxis.set_major_locator(MinuteLocator(interval=5))
#             pcountax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
#             pcountax.set_xlim(starttime,endtime)
#             if(reverse_times):
#                 pcountax.invert_xaxis()
#         else:
#             pcountax.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
#             pcountax.set_xlim(xstart,xstop)
#             pcountax.invert_xaxis()
#         pcountax.yaxis.set_major_locator(MultipleLocator(base=5.0))
#         pcountax.set_yscale('log')
#         pcountax.set_ylim(1.0,6000.0)
#         pcountax.set_ylabel('# of particles')
#         pcountax.set_title('Disdrometer Log(N) and Reflectivity')
#         plt.xticks(visible=False)
        
        
        plt.savefig(image_dir+dis_name+'_logNc_gam.eps',dpi=300)
        
        # This clunky code allows for setting the aspect ratio of the plot+colorbar correctly...
        # I really wish there was a better way to do this and that matplotlib colorbars would respect
        # the aspect ratio of the plots automatically...
        #ax2.set_aspect("auto")
        #divider.set_aspect(True)
        #divider.get_horizontal()[0]._aspect=1000.0

        # Plot exponential/gamma distribution parameters
        
        fig1 = plt.figure(figsize=(8,6))
        ax1 = fig1.add_subplot(211)
        ax1.plot(plotx_dis,N0_exp,ls='-',c='k')
        ax1b = ax1.twinx()
        ax1b.plot(plotx_dis,mu_gam,ls='--',c='k')
        ax1.axhline(y=8.0e6,c='b')  # Plot horizontal line corresponding to M-P value
        if(not timetospace):
            ax1.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax1.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax1.invert_xaxis()
        else:
            ax1.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax1.set_xlim(xstart,xstop)
            ax1.invert_xaxis()
        ax1.set_yscale('log')
        ax1.set_ylim(1.0e3,1.0e8)
        ax1.set_ylabel(r'$N_0$ (exp.)')
        ax1b.set_ylim(-5.0,20.0)
        ax1b.set_ylabel(r'$\mu$ (gamma)')
#       ax2 = fig1.add_subplot(312)
#       ax2.plot(plotx_dis,QR_disd,ls='-',c='k')
#       #ax2.plot_date(date2num(disdates),qr_exp*1000.0,ls='--',c='k')
#       if(not timetospace):
#           ax2.xaxis.set_major_locator(MinuteLocator(interval=5))
#           ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
#           ax2.set_xlim(min(starttimes),max(endtimes))
#           if(reverse_times):
#               ax2.invert_xaxis()
#       else:
#           ax2.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
#           ax2.set_xlim(xstart,xstop)
#           ax2.invert_xaxis()
#       ax2.set_ylim(0.0,2.0)
#       ax2.yaxis.set_major_locator(MultipleLocator(base=0.2))
#       ax2.set_ylabel(r'Mixing ratio ($g kg^{-1}$)')
        ax3 = fig1.add_subplot(212)
        ax3.plot(plotx_dis,D_med_disd,ls='-',c='k')
        ax3.plot(plotx_dis,D_med_exp,ls='--',c='k')
        ax3.plot(plotx_dis,D_med_gam,ls=':',c='k')
        if(not timetospace):
            ax3.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax3.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax3.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax3.invert_xaxis()
        else:
            ax3.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax3.set_xlim(xstart,xstop)
            ax3.invert_xaxis()
        ax3.set_ylim(0.0,6.0)
        ax3.yaxis.set_major_locator(MultipleLocator(base=0.5))
        ax3.set_ylabel(r'$D_{0}$ ($mm$)')
        plt.setp(ax1.get_xticklabels(),visible=False)
        plt.setp(ax2.get_xticklabels(),visible=False)
        
        if(not timetospace):
            fig1.autofmt_xdate()
        
        plt.savefig(image_dir+dis_name+'_exp_dis.eps',dpi=200)
        
        # Plot gamma distribution parameters
        
        fig1 = plt.figure(figsize=(8,6))
        ax1 = fig1.add_subplot(311)
        ax1.plot(plotx_dis,mu_gam,ls='-',c='k')
        if(not timetospace):
            ax1.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax1.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax1.invert_xaxis()
        else:
            ax1.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax1.set_xlim(xstart,xstop)
            ax1.invert_xaxis()
        ax1.set_ylim(-5.0,20.0)
        ax1.set_ylabel(r'$\mu$ (gamma)')
        ax2 = fig1.add_subplot(312)
        ax2.plot(plotx_dis,QR_disd,ls='-',c='k')
        #ax2.plot_date(date2num(disdates),qr_gam*1000.0,ls='--',c='k')
        if(not timetospace):
            ax2.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax2.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax2.invert_xaxis()
        else:
            ax2.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax2.set_xlim(xstart,xstop)
            ax2.invert_xaxis()
        ax2.set_ylim(0.0,10.0)
        ax2.yaxis.set_major_locator(MultipleLocator(base=1.0))
        ax2.set_ylabel(r'Mixing ratio ($g kg^{-1}$)')
        ax3 = fig1.add_subplot(313)
        ax3.plot(plotx_dis,D_med_disd,ls='-',c='k')
        ax3.plot(plotx_dis,D_med_gam,ls='--',c='k')
        if(not timetospace):
            ax3.xaxis.set_major_locator(MinuteLocator(interval=5))
            ax3.xaxis.set_major_formatter(DateFormatter('%H:%M'))
            ax3.set_xlim(min(starttimes),max(endtimes))
            if(reverse_times):
                ax3.invert_xaxis()
        else:
            ax3.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
            ax3.set_xlim(xstart,xstop)
            ax3.invert_xaxis()
        ax3.set_ylim(0.0,15.0)
        ax3.yaxis.set_major_locator(MultipleLocator(base=2.0))
        ax3.set_ylabel(r'$D_{0}$ ($mm$)')
        plt.setp(ax1.get_xticklabels(),visible=False)
        plt.setp(ax2.get_xticklabels(),visible=False)
        
        if(not timetospace):
            fig1.autofmt_xdate()
        
        plt.savefig(image_dir+dis_name+'_gam_dis.eps',dpi=200)

# if(plot_opt):
#     plt.show()


    
