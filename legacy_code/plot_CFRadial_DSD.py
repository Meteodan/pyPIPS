# Plot_CFRadial.py
#
# Plots data from a CFRadial netCDF NEXRAD radar file
# Eventually I will merge part of this script with Plot_Disdrometer


#import netCDF4 as netcdf
#from pylab import *
import Nio
import numpy as N
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import *
from matplotlib.ticker import *
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
import pytz as pytz
from scipy import special
import thermolib as thermo
import obanmodule as oban
import math as M
from scipy.interpolate import interp2d
from scipy.interpolate import Rbf
import glob
import os

# This is a nifty little function picked up from:
# http://stackoverflow.com/questions/5328128/scipy-interpolation-of-large-matrix
# It helps avoid the MemoryError I keep getting using the entire grid in interp2d by only picking
# points on the grid close to the requested points

# def my_interp(X, Y, Z, x, y, spn=33):
#     xs,ys = map(N.array,(x,y))
#     z = N.zeros(xs.shape)
#     for i,(x,y) in enumerate(zip(xs,ys)):
#         # get the indices of the nearest x,y
#         xi = N.argmin(N.abs(X[0,:]-x))
#         yi = N.argmin(N.abs(Y[:,0]-y))
#         xlo = max(xi-spn, 0)
#         ylo = max(yi-spn, 0)
#         xhi = min(xi+spn, X[0,:].size)
#         yhi = min(yi+spn, Y[:,0].size)
#         print "xlo,ylo,xhi,yhi",xlo,ylo,xhi,yhi
#         # make slices of X,Y,Z that are only a few items wide
#         nX = X[xlo:xhi, ylo:yhi]
#         nY = Y[xlo:xhi, ylo:yhi]
#         nZ = Z[xlo:xhi, ylo:yhi]
#         intp = interp2d(nX, nY, nZ,kind='linear')
#         z[i] = intp(x,y)[0]
#     return nX,nY,nZ,z

# def my_interp(X, Y, Z, x, y, spn=21):
#     xs,ys = map(N.array,(x,y))
#     z = N.zeros(xs.shape)
#     for i,(x,y) in enumerate(zip(xs,ys)):
#         # get the indices of the nearest x,y
#         xi = N.argmin(N.abs(X[:,0]-x))
#         yi = N.argmin(N.abs(Y[0,:]-y))
#         xlo = max(xi-spn, 0)
#         ylo = max(yi-spn, 0)
#         xhi = min(xi+spn, X[:,0].size)
#         yhi = min(yi+spn, Y[0,:].size)
#         print "xlo,ylo,xhi,yhi",xlo,ylo,xhi,yhi
#         # make slices of X,Y,Z that are only a few items wide
#         nX = X[ylo:yhi, xlo:xhi]
#         nY = Y[ylo:yhi, xlo:xhi]
#         nZ = Z[ylo:yhi, xlo:xhi]
#         #intp = interp2d(nX, nY, nZ,kind='linear')
#         #z[i] = intp(x,y)[0]
#     return nX,nY,nZ,z

# The following code is taken from the following URL: 
# http://stackoverflow.com/questions/2417794/how-to-make-the-angles-in-a-matplotlib-polar-plot-go-clockwise-with-0-at-the-to
# It defines a Polar projection with 0 degrees at the top and angles increasing clockwise

from matplotlib.projections import PolarAxes, register_projection
from matplotlib.transforms import Affine2D, Bbox, IdentityTransform

class NorthPolarAxes(PolarAxes):
    '''
    A variant of PolarAxes where theta starts pointing north and goes
    clockwise.
    '''
    name = 'northpolar'

    class NorthPolarTransform(PolarAxes.PolarTransform):
        def transform(self, tr):
            xy   = N.zeros(tr.shape, N.float_)
            t    = tr[:, 0:1]
            r    = tr[:, 1:2]
            x    = xy[:, 0:1]
            y    = xy[:, 1:2]
            x[:] = r * N.sin(t)
            y[:] = r * N.cos(t)
            return xy

        transform_non_affine = transform

        def inverted(self):
            return NorthPolarAxes.InvertedNorthPolarTransform()

    class InvertedNorthPolarTransform(PolarAxes.InvertedPolarTransform):
        def transform(self, xy):
            x = xy[:, 0:1]
            y = xy[:, 1:]
            r = N.sqrt(x*x + y*y)
            theta = N.arctan2(y, x)
            return N.concatenate((theta, r), 1)

        def inverted(self):
            return NorthPolarAxes.NorthPolarTransform()

    def _set_lim_and_transforms(self):
        PolarAxes._set_lim_and_transforms(self)
        self.transProjection = self.NorthPolarTransform()
        self.transData = (
            self.transScale + 
            self.transProjection + 
            (self.transProjectionAffine + self.transAxes))
        self._xaxis_transform = (
            self.transProjection +
            self.PolarAffine(IdentityTransform(), Bbox.unit()) +
            self.transAxes)
        self._xaxis_text1_transform = (
            self._theta_label1_position +
            self._xaxis_transform)
        self._yaxis_transform = (
            Affine2D().scale(N.pi * 2.0, 1.0) +
            self.transData)
        self._yaxis_text1_transform = (
            self._r_label1_position +
            Affine2D().scale(1.0 / 360.0, 1.0) +
            self._yaxis_transform)

register_projection(NorthPolarAxes)

deg2rad = N.pi/180.

# Time zone stuff

fmt = '%Y-%m-%d %H:%M:%S %Z%z'
fmt2 = '%Y-%m-%d %H:%M:%S'

# Disdrometer lat/lons

# ------ 2009/05/15 Oklahoma squall line ------
# 
# # Katja Deployment 1,2 and Glen deployment 1
# dlocs = [(36.55714,-98.42660),(36.66674,-98.35547),(36.50682,-98.24802),(36.47796,-98.24765)]
# dis_name_list = ['CU_D1_P2','CU_D2_P1','D1_P1','D1_P2']
# 
# starttime = datetime(2009,5,15,17,20,0)
# endtime = datetime(2009,5,15,19,25,0)
# 
# plotxmin =  -60000.0
# plotxmax =   20000.0
# plotymin =  -60000.0
# plotymax =   20000.0
# 
# # Radar name
# 
# radar_name = 'KVNX'
# 
# # Directory containing NEXRAD CFRadial files
# nexrad_dir = '/Users/ddawson/NSF_V2_study/Disdrometer_study/NEXRAD_CFradial/KVNX/20090515/'
# full_nexrad_filelist = glob.glob(nexrad_dir+'cfrad*.nc')
# 
# # Min range to plot (m)
# minrange = 0.0
# # Max range to plot (m)
# maxrange = 200000.0
# # Min azimuth to plot (deg)
# minazim = 0.
# # Max azimuth to plot (deg)
# maxazim = 360.

# ------------

# ------ 2009/06/05 LaGrange, WY supercell ------
# 
# # Glen Deployment 1
# dlocs = [(41.65121,-104.363),(41.64994,-104.3924)]
# dis_name_list = ['D1_P1','D1_P2'] # i.e. deployment #, probe #
# 
# starttime = datetime(2009,6,5,15,56,00)
# endtime = datetime(2009,6,5,16,21,00)
# 
# plotxmin = 20000.0
# plotxmax = 50000.0
# plotymin = 45000.0
# plotymax = 65000.0
# 
# # Glen Deployment 2
# dlocs = [(41.5054,-103.6582),(41.5215,-103.6596)]
# dis_name_list = ['D2_P1','D2_P2']
# 
# starttime = datetime(2009,6,5,17,51,00)
# endtime = datetime(2009,6,5,18,15,00)
# 
# plotxmin = 70000.0
# plotxmax = 110000.0
# plotymin = 30000.0
# plotymax = 70000.0
# 
# # Katja Deployment 1
# dlocs = [(41.64579,-104.24197)]
# dis_name_list = ['CU_D1_P1']
# 
# starttime = datetime(2009,6,5,15,43,0)
# endtime = datetime(2009,6,5,16,48,0)
# 
# plotxmin = 20000.0
# plotxmax = 60000.0
# plotymin = 40000.0
# plotymax = 80000.0
# 
# # Radar name
# 
# radar_name = 'KCYS'
# 
# # Directory containing NEXRAD CFRadial files
# nexrad_dir = '/Users/ddawson/NSF_V2_study/Disdrometer_study/NEXRAD_CFradial/KCYS/20090605/'
# full_nexrad_filelist = glob.glob(nexrad_dir+'cfrad*.nc')
# 
# # Min range to plot (m)
# minrange = 0.0
# # Max range to plot (m)
# maxrange = 200000.0
# # Min azimuth to plot (deg)
# minazim = 0.
# # Max azimuth to plot (deg)
# maxazim = 360.

# ------------

# ------ 2009/06/07 Oregon, MO supercell ------
# 
# # Glen Deployment 1
# dlocs = [(39.986568,-95.204717),(40.012383,-95.219788)]
# dis_name_list = ['D1_P1','D1_P2'] # i.e. deployment #, probe #
# 
# starttime = datetime(2009,6,7,17,29,00)
# endtime = datetime(2009,6,7,18,04,00)
# 
# plotxmin = -100000.0
# plotxmax = -50000.0
# plotymin = 110000.0
# plotymax = 150000.0

# Glen Deployment 2
dlocs = [(39.8786967,-94.3589417),(39.895135,-94.3574367)]
dis_name_list = ['D2_P1','D2_P2']

starttime = datetime(2009,6,7,19,10,00)
endtime = datetime(2009,6,7,19,53,00)

plotxmin = -30000.0
plotxmax = 10000.0
plotymin = 100000.0
plotymax = 140000.0

# # Katja Deployment 1
# dlocs = [(40.18403,-95.32748)]
# dis_name_list = ['CU_D1_P2']
# 
# starttime = datetime(2009,6,7,17,17,00)
# endtime =   datetime(2009,6,7,17,29,00)
# 
# plotxmin = -100000.0
# plotxmax = -50000.0
# plotymin = 110000.0
# plotymax = 160000.0

# Katja Deployment 2
# dlocs = [(40.09222,-94.87020)]
# dis_name_list = ['CU_D2_P2']
# 
# starttime = datetime(2009,6,7,18,9,00)
# endtime =   datetime(2009,6,7,18,20,00)
# 
# plotxmin = -70000.0
# plotxmax = -20000.0
# plotymin = 110000.0
# plotymax = 160000.0

# Katja Deployment 3
# dlocs = [(40.04708,-94.35077)]
# dis_name_list = ['CU_D2_P2']
# 
# starttime = datetime(2009,6,7,18,59,00)
# endtime =   datetime(2009,6,7,19,04,00)
# 
# plotxmin = -70000.0
# plotxmax = -20000.0
# plotymin = 110000.0
# plotymax = 160000.0

# Katja Deployment 4
# dlocs = [(39.97626,-94.15823)]
# dis_name_list = ['CU_D2_P2']
# 
# starttime = datetime(2009,6,7,19,41,00)
# endtime =   datetime(2009,6,7,20,01,00)
# 
# plotxmin = -10000.0
# plotxmax = 40000.0
# plotymin = 100000.0
# plotymax = 150000.0

# Radar name

radar_name = 'KEAX'

# Directory containing NEXRAD CFRadial files
nexrad_dir = '/Users/ddawson/NSF_V2_study/Disdrometer_study/NEXRAD_CFradial/KEAX/20090607/'
full_nexrad_filelist = glob.glob(nexrad_dir+'cfrad*.nc')

# Min range to plot (m)
minrange = 0.0
# Max range to plot (m)
maxrange = 200000.0
# Min azimuth to plot (deg)
minazim = 0.
# Max azimuth to plot (deg)
maxazim = 360.

# ------------

# ------ 2009/06/09 Ford-Greensburg, KS supercell ------

# dlocs = [(37.66112,-99.64175),(37.66858,-99.63229),(37.66395,-99.63215),(37.67288,-99.63229)]
# dis_name_list = ['D1_P1','D1_P2','CU_D1_P1','CU_D1_P2'] # i.e. deployment #, probe #
# 
# starttime = datetime(2009,6,9,17,23,00)
# endtime = datetime(2009,6,9,17,47,00)
# 
# plotxmin = 20000.0
# plotxmax = 50000.0
# plotymin = -20000.0
# plotymax = 10000.0
# 
# # Radar name
# 
# radar_name = 'KDDC'
# 
# # Directory containing NEXRAD CFRadial files
# nexrad_dir = '/Users/ddawson/NSF_V2_study/Disdrometer_study/NEXRAD_CFradial/KDDC/20090609/'
# full_nexrad_filelist = glob.glob(nexrad_dir+'cfrad*.nc')
# 
# # Min range to plot (m)
# minrange = 0.0
# # Max range to plot (m)
# maxrange = 100000.0
# # Min azimuth to plot (deg)
# minazim = 0.
# # Max azimuth to plot (deg)
# maxazim = 360.

# ------------

# Plot approximate locations of storm at certain time intervals based on assume storm motion?
plot_interval = False
storm_u = 12.55
storm_v = 0.0
time_interval = 30.0
num_intervals = 2   #Number of intervals (plus and minus) to plot around each center time and associated scan

raddates = []
nexrad_filelist = []
# Only choose volumes within desired timerange
for path in full_nexrad_filelist:
    
    dummy,filename = os.path.split(path)
    
    print filename
    
    radyear = int(filename[6:10])
    radmonth = int(filename[10:12])
    radday = int(filename[12:14])
    radhour = int(filename[15:17])
    radmin = int(filename[17:19])
    radsec = int(filename[19:21])
    
    sweeptime = datetime(radyear,radmonth,radday,radhour,radmin,radsec)-timedelta(hours=6) # Get it to CST time
    
    if(sweeptime >= starttime and sweeptime <= endtime):
        raddates.append(sweeptime)
        nexrad_filelist.append(path)
    
dBZ_D_tlist = []

average_gates = False        # True if reflectivity should be averaged in the closest 9 gates
                            # False if just picking the value from the closest gate
                                
Cressman = True             # Perform a Cressman analysis on nearby radar gates to determine
                            # reflectivity value

roi = 750.                   # Radius of influence of Cressman analysis in m

image_dir = nexrad_dir+'/images/'
if(not os.path.exists(image_dir)):
    os.makedirs(image_dir)

for path,sweeptime in zip(nexrad_filelist,raddates):

    print "Opening file: ",path
    sweepfile_netcdf = Nio.open_file(path)

    # Grab the time information from the file
     
    print "Time of sweep = ",sweeptime.strftime(fmt)
        
    # Grab some needed variables from the netCDF file
    
    numgates = sweepfile_netcdf.dimensions['range']     # Number of gates
    numtimes = sweepfile_netcdf.dimensions['time']     # Number of times
    numsweeps = sweepfile_netcdf.dimensions['sweep']     # Number of sweeps
    sweep_start_ray_index = sweepfile_netcdf.variables['sweep_start_ray_index'][:]
    sweep_end_ray_index = sweepfile_netcdf.variables['sweep_end_ray_index'][:]
    ray_n_gates = sweepfile_netcdf.variables['ray_n_gates'][:]
            
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
    
    beam_width = sweepfile_netcdf.variables['radar_beam_width_h'].get_value()
    #print "Radar beam width (degrees): "+str(beam_width)
    
    azimuth_rad = sweepfile_netcdf.variables['azimuth'][:]  # Azimuth of center of each gate
    azimuth_rad = azimuth_rad[:sweep_end_ray_index[0]] # Just grab first sweep for now
    
    #print "number of azimuths in sweep ",N.size(azimuth_rad)
    
    # Roll the azimuth dimension around so that it starts at 0 and ends at 360
    
    shift = N.where(azimuth_rad < azimuth_rad[0])[0][0]
    
    #print "shift is ",shift
    azimuth_rad = N.roll(azimuth_rad,shift=-shift)
    beamwidth = azimuth_rad[1:]-azimuth_rad[0:-1]
        
    # For plotting we need the azimuth at the borders of each ray (start of each ray plus one more on the end)
    azimuth_start_rad = N.zeros((N.size(azimuth_rad)+1))
    #azimuth_start_rad = N.zeros_like(azimuth_rad)
    # Find azimuth of "start" of each gate (in azimuth) -- approximate
    azimuth_start_rad[1:-1] = azimuth_rad[1:]-0.5*beamwidth[:]
    azimuth_start_rad[0] = azimuth_rad[0]-0.5*beamwidth[0]
    azimuth_start_rad[-1] = azimuth_rad[-1]+0.5*beamwidth[-1]
    
    #print azimuth_start_rad
    #print azimuth_rad
    
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
    numpointsinsweep = ray_n_gates[0]*sweep_end_ray_index[0]
    
    #print "numpointsinsweep",numpointsinsweep
    #print N.size(dBZ[:numpointsinsweep])
    dBZ = dBZ[:numpointsinsweep]
    dBZ = dBZ.reshape((sweep_end_ray_index[0],-1))
    #print "dBZ.shape = ",dBZ.shape
    
    #print "scale_factor,add_offset",tempdBZ.scale_factor[0],tempdBZ.add_offset[0]
    
    #Unpack values
    dBZ = dBZ*tempdBZ.scale_factor[0]+tempdBZ.add_offset[0] 
    
    # Adjust azimuth axis
    
    dBZ = N.roll(dBZ,shift=-shift,axis=0)
        
    # Quick plot for testing
    
    #fig1 = plt.figure()
    #ax1 = fig1.add_subplot(111,projection='northpolar')
    #ax1 = fig1.add_subplot(111,polar=True)
    #ax1 = fig1.add_subplot(111)
    
    #print azimuth_rad.shape
    #print range.shape
    #print dBZ.shape
    
    minrangeindex = N.where(range_start > minrange)[0][0]
    maxrangeindex = N.where(range_start > maxrange)[0][0]
    #print minrangeindex,maxrangeindex
    
    minazimindex = N.where(azimuth_start_rad/deg2rad >= min(azimuth_start_rad[0]/deg2rad,minazim))[0][0]
    print "azimuth_start_rad[minazimindex] = ",azimuth_start_rad[minazimindex]/deg2rad
    try:
        maxazimindex = N.where(azimuth_start_rad/deg2rad >= maxazim)[0][0]
    except:
        maxazimindex = N.size(azimuth_start_rad)-1
    
    print "azimuth_start_rad[maxazimindex] = ",azimuth_start_rad[maxazimindex]
    
    #print minazimindex,maxazimindex
    
    # if(maxazim < minazim):  # Sector wraps around end of array, roll it around till it lines up properly
    #     azimuth_rad_hi = N.roll(azimuth_rad_hi,shift=maxazimindex)
    #     Zh_hi = N.roll(Zh_hi,shift=maxazimindex,axis=0)
    #     minazimindex = minazimindex+maxazimindex
    #     maxazimindex = N.size(azimuth_rad_hi)
    
    #print minazimindex,maxazimindex
    
    theta,rad = N.meshgrid(azimuth_start_rad[minazimindex:maxazimindex+1],range_start[minrangeindex:maxrangeindex])
    theta_c,rad_c = N.meshgrid(azimuth_rad[minazimindex:maxazimindex+1],range[minrangeindex:maxrangeindex])
    dBZplt = dBZ.swapaxes(0,1)[minrangeindex:maxrangeindex,minazimindex:maxazimindex+1]
    
    dBZplt = N.ma.masked_invalid(dBZplt)
    
    #print theta.shape
    #print rad.shape
    
    #ax1.pcolormesh(theta,rad,dBZplt)
    
    # Plot on rectangular mesh using x,y points derived from r,theta
    # The points are valid on the edges of the gates, while the reflectivity values
    # are valid at the centers.  This is so that pcolor will plot the gates in the correct
    # locations (hopefully).
    
    xplt,yplt = oban.xyloc(rad,theta,el,rlat,rlon,ralt,1)
    xplt_c,yplt_c = oban.xyloc(rad_c,theta_c,el,rlat,rlon,ralt,1)
    
    #print "xplt,yplt,dBZplt shape",xplt.shape,yplt.shape,dBZplt.shape
    
    # Decide whether to perform a shift in time and space 
    subtimes = []
    if(plot_interval):
        for dt in N.arange(-num_intervals*time_interval,num_intervals*time_interval+time_interval,time_interval):
            subtime = sweeptime+timedelta(seconds=dt)
            subtimes.append(subtime)
    else:
        subtimes.append(sweeptime)
    
    for dt,time in zip(N.arange(-num_intervals*time_interval,num_intervals*time_interval+time_interval,time_interval),subtimes):   
        if(plot_interval):
            xplttmp = xplt+dt*storm_u
            xplt_ctmp = xplt_c+dt*storm_u
        else:
            xplttmp = xplt
            xplt_ctmp = xplt_c
        
        yplttmp = yplt
        yplt_ctmp = yplt_c
    
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        plot1 = ax2.pcolormesh(xplttmp,yplttmp,dBZplt)
        #plot2 = ax2.contour(xplt_c,yplt_c,dBZplt,levels=[40.0],c='k')
        
        for dname,dloc in zip(dis_name_list,dlocs):
            Dx,Dy = oban.ll_to_xy(dloc[0]*deg2rad,dloc[1]*deg2rad,rlat,rlon,1)
        
            #print "Dx,Dy",Dx,Dy
            
            # Find closest gate to disdrometer location
            # First, find index along range; first compute slant range
            
    #         srange_dis = oban.computeslantrange([Dx],[Dy],el)[0][0]
    #         rad_diff = srange_dis-rad_c[:,0]    # Can use 0 index here; range should be independent of azimuth (mostly)
    #         srange_index = N.argmin(N.abs(rad_diff))
    #         
    #         print "range of disdrometer "+dname+": "+str(srange_dis)
    #         print srange_index
    #         
    #         # Next find the closest azimuth
    #         
    #         theta_dis = -N.arctan2(Dy,Dx)+N.pi/2.       # I *think* this is correct
    #         
    #         print "azimuth of disdrometer "+dname+": "+str(theta_dis/deg2rad)
    #         #Find closest index of azimuth
    #         
    #         theta_diff = theta_dis-theta_c[srange_index,:]
    #         theta_index = N.argmin(N.abs(theta_diff))
    #         print theta_index
            
            # Try another way to get the indices (the above seems to have problems)
            # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why the 
            # above doesn't work in all situations, but it sometimes picks a gate adjacent to the one
            # we want...
            
            # First, compute the euclidian distance of the x,y location of the disdrometer to each of the
            # radar gate centers
            distance = N.sqrt(((Dx-xplt_ctmp)**2.+(Dy-yplt_ctmp)**2.))
            
            # Now, find the index of the closest radar gate
            srange_index,theta_index = N.unravel_index(distance.argmin(),distance.shape)
            print "srange_index,theta_index",srange_index,theta_index
            
            # Finally, grab reflectivity at closest gate to disdrometer
            # Average the reflectivity in the closest gate and surrounding 8 gates if desired
            if(Cressman):
                #print xplt_c.shape,yplt_c.shape
                dBZ_D = oban.Cresmn(xplt_c,yplt_c,dBZplt,Dx,Dy,roi)
                #print dBZ_D,dBZ_D.shape
            elif(not average_gates):
                dBZ_D = dBZplt[srange_index,theta_index]
            else:
                dBZ_D = (1./9.)*(dBZplt[srange_index-1,theta_index-1]+dBZplt[srange_index,theta_index-1]+
                              dBZplt[srange_index-1,theta_index]+dBZplt[srange_index,theta_index]+
                              dBZplt[srange_index+1,theta_index]+dBZplt[srange_index,theta_index+1]+
                              dBZplt[srange_index+1,theta_index+1]+dBZplt[srange_index+1,theta_index-1]+
                              dBZplt[srange_index-1,theta_index+1])
            
            print "Reflectivity at disdrometer "+dname+": "+"%.2f"%dBZ_D
            
            ax2.plot(Dx,Dy,'y*',ms=5)
            
            #ax2.plot(xplt[srange_index,theta_index],yplt[srange_index,theta_index],'rx',ms=10)
            #ax2.plot(xplt_c[srange_index,theta_index],yplt_c[srange_index,theta_index],'bx',ms=10)
            
            #ax1.plot(theta[srange_index,theta_index],rad[srange_index,theta_index],'rx',ms=10)
            
            plt.annotate(dname+" %.2f"%dBZ_D,(Dx,Dy))
        
        # Determine plot bounds
        if(plotxmin == -1):
            plotxmin = N.min(xplt)
        if(plotxmax == -1):
            plotxmax = N.max(xplt)
        if(plotymin == -1):
            plotymin = N.min(yplt)
        if(plotymax == -1):
            plotymax = N.max(yplt)
        
        ax2.set_xlim(plotxmin,plotxmax)
        ax2.set_ylim(plotymin,plotymax)
        
        plt.colorbar(plot1)
        plt.title('dBZ at el = %.2f'%(el/deg2rad)+'and time '+time.strftime(fmt))
        ax2.set_aspect('equal')
    
        plt.savefig(image_dir+'dbZ_disd_t'+time.strftime(fmt).strip()+'.png',dpi=200)
    
plt.show()

# # Get rid of negative infinity values (set them to 0)
# 
# Zh = N.where(N.isneginf(Zh),0.0,Zh)
# Zh = N.where(N.isposinf(Zh),0.0,Zh)
# Zh = N.where(N.isnan(Zh),0.0,Zh)
# 
# # Quick plot for testing
# 
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111,polar=True)
# 
# ax1.pcolor(azimuth_rad,gates,Zh.swapaxes(0,1))
# 
# plt.show()



