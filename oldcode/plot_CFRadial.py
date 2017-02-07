# Plot_CFRadial.py
#
# Plots data from a CFRadial netCDF radar file
# Eventually I will merge part of this script with Plot_Disdrometer


#import netCDF4 as netcdf
#from pylab import *
import Nio
import numpy as N
import matplotlib
matplotlib.use('TkAgg')
#matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import *
from matplotlib.ticker import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import ImageGrid
import ctablesfrompyesviewer as ctables
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
import sys
import radarmodule as radar

fmt = '%Y-%m-%d %H:%M:%S %Z%z'
fmt2 = '%Y-%m-%d %H:%M:%S'
fmt3 = '%Y%m%d%H%M%S'

def mtokm(val,pos):
    """Convert m to km for formatting axes tick labels"""
    val=val/1000.0
    return '%i' % val

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

# Plot approximate locations of storm at certain time intervals based on assume storm motion?
plot_interval = False
storm_u = 12.55
storm_v = 0.0
time_interval = 30.0
num_intervals = 2   #Number of intervals (plus and minus) to plot around each center time and associated scan

# Parse command line argument and read in disdrometer/radar information from text file
# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
# version first, and if not, reading the textfile and pickling the read-in variables for future read-ins.  Think about this.

with open(sys.argv[1],'r') as inputfile:
    # Read and discard header
    inputfile.readline()
    inputfile.readline()
    inputfile.readline()
    numdis = N.int(inputfile.readline().strip().split()[0]) # Read in number of disdrometers
    
    # Read in disdrometer locations and names
    dlocs=[]
    dis_name_list = []
    for l in xrange(numdis):
        line = inputfile.readline().strip().split(',')
        dloc = (N.float(line[0]),N.float(line[1])) # Tuple containing lat/lon
        dname = line[2] # Disdrometer name
        # Append disdrometer locations and names to lists
        dlocs.append(dloc)
        dis_name_list.append(dname)
    
    # Read in start and end times
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_int = map(int,line)
    starttime = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    line = inputfile.readline().strip().split(',')
    line_int = map(int,line)
    endtime = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
    
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
    
    print line[4]
    print N.float(line[4])
    try:
        el = N.float(line[4])
        print "elevation angle",el
    except:
        el = 0.5    # Default to 0.5 degrees
    
    # Read in directory containing radar files
    inputfile.readline()
    radar_dir = inputfile.readline().strip()
    
    # Read in min range,max range, min azimuth, and max azimuth for radar plotting (may deprecate this)
    inputfile.readline()
    line = inputfile.readline().strip().split(',')
    line_float = map(float,line)
    minrange = line_float[0]
    maxrange = line_float[1]
    minazim = line_float[2]
    maxazim = line_float[3]

full_radar_filelist = glob.glob(radar_dir+'cfrad*.nc')

radlims = [minrange,maxrange,minazim,maxazim]
plotlims = [plotxmin,plotxmax,plotymin,plotymax]

raddates = []
radar_filelist = []
# Only choose volumes within desired timerange
for path in full_radar_filelist:
    
    dummy,filename = os.path.split(path)
    
    print filename
    
    radyear = int(filename[6:10])
    radmonth = int(filename[10:12])
    radday = int(filename[12:14])
    radhour = int(filename[15:17])
    radmin = int(filename[17:19])
    radsec = int(filename[19:21])
    
    #el = float(filename[40:44])
    #print el
    
    sweeptime = datetime(radyear,radmonth,radday,radhour,radmin,radsec)
    #sweeptime = datetime(radyear,radmonth,radday,radhour,radmin,radsec)-timedelta(hours=6) # Get it to CST time
    
    if(sweeptime >= starttime and sweeptime <= endtime): # and el == 0.5):
        raddates.append(sweeptime)
        radar_filelist.append(path)
    
dBZ_D_tlist = []

ovrdis = False               # Overlay disdrometer locations on radar plot?

average_gates = False        # True if reflectivity should be averaged in the closest 9 gates
                            # False if just picking the value from the closest gate
                                
Cressman = True             # Perform a Cressman analysis on nearby radar gates to determine
                            # reflectivity value

roi = 750.                   # Radius of influence of Cressman analysis in m

clevels_ref = N.arange(5.0,85.0,5.0)          # Contour levels for reflectivity (dBZ)
clevels_zdr = N.arange(0.0,6.25,0.25)         # Contour levels for Zdr (dB)
clevels_vr  = N.arange(-40.0,41.0,1.0)        # Contour levels for Vr (m/s)
cmapdBZ = ctables.__getattribute__('REF_default')
cmapzdr = cm.Reds
cmapvr = cm.RdBu_r

fieldnames = ['dBZ','Zdr','Vr']

image_dir = radar_dir+'/images/'
if(not os.path.exists(image_dir)):
    os.makedirs(image_dir)

for path,sweeptime in zip(radar_filelist,raddates):

    fieldlist,range_start,range,azimuth_start_rad,azimuth_rad,rlat,rlon,ralt,el = \
        radar.readCFRadial(True,el,rlat,rlon,ralt,path,sweeptime,fieldnames)
    
    print "rlat,rlon,ralt,el",rlat,rlon,ralt,el
    
    # Prepare masks for fields by reflectivity < some threshold
    
    for field,fieldname in zip(fieldlist,fieldnames):   # Probably should do this using a dictionary
        if(fieldname == 'dBZ'):
            mask = N.where(field > 5.0,False,True)
            
    masklist = [mask,mask,mask]
    
    # Compute x,y locations of disdrometer (assumes lambert conformal like plotsweep for now)
    dxy_list = []
    if(ovrdis):
        dxy_list = []
        for dloc,dname in zip(dlocs,dis_name_list):
            Dx,Dy = oban.ll_to_xy(dloc[0]*deg2rad,dloc[1]*deg2rad,rlat*deg2rad,rlon*deg2rad,1)
            dxy_list.append((Dx,Dy))
            h,r = oban.computeheightrangesingle(Dx,Dy,el*deg2rad)
            print "Disdrometer name,x,y,radar elevation angle,slant range, beam height:"
            print dname,Dx,Dy,el,r,h
            if(dloc == dlocs[0] and plotxmin == -1):
                plotlims = [Dx-20000.,Dx+20000.,Dy-20000.,Dy+20000.]
            
    figlist,gridlist = radar.plotsweep(radlims,plotlims,fieldnames,fieldlist,masklist,
                    range_start,range,azimuth_start_rad,azimuth_rad,rlat,rlon,ralt,el,False,
                    ovrdis,dis_name_list,dxy_list)
    
    # Save figures
    for fieldname,fig in zip(fieldnames,figlist):
        plt.figure(fig.number)   
        plt.savefig(image_dir+fieldname+sweeptime.strftime(fmt3).strip()+'el'+str(el)+'.png',dpi=200,bbox_inches='tight')
     
#plt.show()
