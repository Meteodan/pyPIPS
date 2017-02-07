# obanmodule.py
# A set of helpful functions to facilitate reading output from OPAWS and comparing with COMMAS
# output

import netCDF4 as netcdf
import numpy as N
from netcdftime import utime
import sys
import os
from . import markup
from datetime import datetime,timedelta
from scipy.interpolate import interp1d
from scipy import ndimage


re = 6367.0             # Radius of Earth in km
rem = re*1000.0
ere = (4./3.)*re        # Effective radius of earth (for typical atmospheric refraction of radar beam)
erem = ere*1000.0
deg2rad = N.pi/180.

def vintrp(model_z,swp_height,var_col):
    """Given a model column, the heights of each point in the column, and the height of the sweep
       surface valid for that column, vertically interpolate to the sweep height"""
    
    print "model_z.shape",model_z.shape
    
    intrpz = interp1d(model_z,var_col)
    swpvar = interpz(swp_height)
    return swpvar

def intrp2swp(swp_height,model_z,var):
    """Given a 2D array of heights of the radar sweep surface at each x,y grid location, and a 3D
       model field, interpolate (in the vertical) from the 3D model field to the sweep surface"""
    swpvar = N.zeros_like(var[...,0])
    
    #vec_vintrp = N.vectorize(vintrp)
    
    #print "model_z.shape",model_z.shape
    
    #swpvar[:] = vec_vintrp(model_z[:],swp_height[:,:],var[:,:,:])
    
    # Find the indices of the model level above and below the sweep surface for each column
    
    for index,height in N.ndenumerate(swp_height):
        k = N.searchsorted(model_z,height)
        h = (height-model_z[k-1])/(model_z[k]-model_z[k-1])
        swpvar[index] = (1.0-h)*var[index[0],index[1],k-1]+h*var[index[0],index[1],k]
    
    return swpvar

def findswpindices(swp_height,model_z):
    """Given a 2D array of heights of the radar sweep surface at each x,y grid location,
       and an array of model heights, find the index of the model heights below
       the sweep surface for each x,y location."""
    
    z_index = N.zeros_like(swp_height,dtype='int32')
    
    # Find the vertical index of the model level below the sweep surface for each column
    
    for index,height in N.ndenumerate(swp_height):
        if(model_z.ndim == 1):
            z_column = model_z
        else:
            z_column = model_z[index]
        z_index[index] = max(N.searchsorted(z_column,height)-1,0)
    
    return z_index
    
def intrp2swp2(swp_height,model_z,var):
    """Given a 2D array of heights of the radar sweep surface at each x,y grid location, and a 3D
       model field, interpolate (in the vertical) from the 3D model field to the sweep surface"""
    swpvar = N.zeros_like(var[...,0])
    
    for index,height in N.ndenumerate(swp_height):
        if(model_z.ndim == 1):
            z_column = model_z
        else:
            z_column = model_z[index]
        #intrpz = interp1d(z_column,var[index[0],index[1],:],bounds_error=False)
        #swpvar[index] = intrpz(swp_height[index])
        # EDIT 05/22/2013: Use numpy's interp routine instead
        swpvar[index] = N.interp(swp_height[index],z_column,var[index[0],index[1],:])
    
    return swpvar
    
def intrp2swp3(swp_height,model_z,var):
    """Given a 2D array of heights of the radar sweep surface at each x,y grid location, and a 3D
       model field, interpolate (in the vertical) from the 3D model field to the sweep surface.
       This version handles multiple tilts at once."""
       
    swp_var = N.zeros_like(swp_height)
        
    for index,dummy in N.ndenumerate(swp_height[...,0]):
        intrpz = interp1d(model_z,var[index[0],index[1],:],bounds_error=False)
        swp_var[index[0],index[1],:] = intrpz(swp_height[index[0],index[1],:])
    
    return swp_var
    
def computeheightrangesingle(x,y,el):
    """Given x,y coordinates relative to radar location, and the radar elevation, computes the
       height and slant path range of the radar beam at a single x,y location.  Modified from OPAWS compute_height function"""
    r = N.sqrt(x*x+y*y)/N.cos(el)     #Slant range
    h = N.sqrt(r*r+erem*erem+2.0*r*erem*N.sin(el))-erem
    return h,r

def computeheight(x,y,el):
    """Given x,y coordinates relative to radar location, and the radar elevation, computes the
       height of the radar beam at each x,y point.  Modified from OPAWS compute_height function"""
    X,Y = N.meshgrid(x,y)
    r = N.sqrt(X*X+Y*Y)/N.cos(el)     #Slant range
    h = N.sqrt(r*r+erem*erem+2.0*r*erem*N.sin(el))-erem
    return h.swapaxes(0,1)

def computeslantrange(x,y,el):
    """Given x,y coordinates relative to radar location, and the radar elevation angle, compute 
       the slant range to the radar gate."""
    X,Y = N.meshgrid(x,y)
    r = N.sqrt(X*X+Y*Y)/N.cos(el)
    return r

def computeheight2(x,y,elev_list):
    """Given x,y coordinates relative to radar location, and the radar elevation, computes the
       height of the radar beam at each x,y point.  Modified from OPAWS compute_height function.
       This version computes for multiple elevation angles at once."""
    
    X,Y = N.meshgrid(x,y)
    
    swp_height = N.zeros((N.size(X,0),N.size(X,1),len(elev_list)))
    
    for k,el in enumerate(elev_list):
        r = N.sqrt(X*X+Y*Y)/N.cos(el)   # Slant range
        swp_height[...,k] = N.sqrt(r*r+erem*erem+2.0*r*erem*N.sin(el))-erem
    
    return swp_height

def llobs(r,az,el,rlat,rlon,ralt):
    """Given slant-path range, azimuth, and elevation angle of a radar observation, computes the 
       lat/lon of the observation.  Modified from David Dowell's xyzloc subroutine in OPAWS"""
    
    # Note r is assumed to be in meters
    
    # First compute height of radar observation
    h = N.sqrt(r*r+erem*erem + 2.0*r*erem*N.sin(el))-erem
    h = h + ralt
    s = erem*N.arcsin((r*N.cos(el))/(erem+h))
    d = s/rem
    
    # Now compute the observation lat/lon
    lat = N.arcsin(N.sin(rlat)*N.cos(d) + N.cos(rlat)*N.sin(d)*N.cos(az))
    dlon = N.arctan2(N.sin(az)*N.sin(d)*N.cos(rlat), N.cos(d)-N.sin(rlat)*N.sin(lat))
    lon = N.mod(rlon+dlon+N.pi, 2.0*N.pi) - N.pi
    
    return lat,lon
    
def ll_to_xy(lat,lon,rlat,rlon,map_proj=0):
    """Given a lat/lon pair, compute the x,y location relative to the radar lat/lon using the given
       map projection.  Based on David Dowell's ll_to_xy function in OPAWS as modified by Mike 
       Coniglio."""
    
    tlat1 = 30.0    # True latitude 1 for LC
    tlat2 = 60.0    # True latitude 2 for LC
    clon = rlon     # Just choose radar longitude for center longitude
    
    if(map_proj == 0):                              # Flat Earth approximation
        x = rem*N.cos(0.5*(rlat+lat))*(lon-rlon)
        y = rem*(lat-rlat)
    elif(map_proj == 1):                            # Lambert conformal
        tlat   = deg2rad*tlat1
        blat   = deg2rad*tlat2
        cenlon = clon
        b      = (N.log(N.sin(tlat))-N.log(N.sin(blat)))/(N.log(N.tan(tlat/2.))-N.log(N.tan(blat/2.)))
        c1     = rem*N.sin(tlat)/(b*N.tan(tlat/2.)**b)
        # first find the x and y of the first point
        xlat   = (N.pi/2.)-rlat
        xlon   = rlon
        d      = c1*(N.sin(xlat)/(1.+N.cos(xlat)))**b
        x1     = d*N.sin(b*(xlon-cenlon))
        y1     = d*N.cos(b*(xlon-cenlon))
        # now find the x and y of the second point
        xlat   = (N.pi/2.)-lat
        xlon   = lon
        d      = c1*(N.sin(xlat)/(1.+N.cos(xlat)))**b
        x2     = d*N.sin(b*(xlon-cenlon))
        y2     = d*N.cos(b*(xlon-cenlon))
        # subtract to find x, y distances (km)
        x      = x2-x1
        y      = y1-y2  # makes positive values point to north instead of south
    return x,y


def xyloc(r,az,el,rlat,rlon,ralt,map_proj=0):
    """Given slant-path range, azimuth, and elevation angle of a radar observation, computes the x,y 
       location relative to the radar location of the observation.  Modified from David Dowell's
       xyzloc subroutine in OPAWS"""
       
    lat,lon = llobs(r,az,el,rlat,rlon,ralt)
       
    x,y = ll_to_xy(lat,lon,rlat,rlon,map_proj)
       
    return x,y


def oban2D2model(sweepdir,obantemplatefile,timecorrection,stime,target_elev,timetol):
    """Create an oban on a single sweep on the model grid closest to the
       requested time."""
    

def readobangrid(obanfilename):
    """Read grid information from oban file"""
    
    root = netcdf.Dataset(obanfilename, "r")
    
    ref_lat = root.variables['grid_latitude'][:]
    ref_lon = root.variables['grid_longitude'][:]
    
    root.close()
    
    return ref_lat,ref_lon

def readsweep(obanfilename,stime,target_elev,timetol):
    """Read a single sweep from the 2D oban module corresponding to the requested elevation and 
       requested time.  Returns dBZ, Zdr, Kdp, and Vr if they exist, otherwise just return False booleans"""

    # Open oban file
    
    root  = netcdf.Dataset(obanfilename, "r")

    dBZ_found = False
    Zdr_found = False
    Kdp_found = False
    rhv_found = False
    Vr_found = False
    
    x     = root.variables['x'][:]
    y     = root.variables['y'][:]

    # Attempt to read dBZ
    try:
        dBZ = root.variables['DZ']
        dBZ_found = True
    except: 
        print "\n Cannot find variable DZ"
    
    if(not dBZ_found):  
        try:
            dBZ = root.variables['DBZ']
            dBZ_found = True
        except:
            print "\n Cannot find variable DBZ"
    
    # Attempt to read Zdr
    try:
        Zdr = root.variables['DR']
        Zdr_found = True
    except:
        print "\n Cannot find variable DR"
    
    # Attempt to read Kdp
    try:
        Kdp = root.variables['KD']
        Kdp_found = True
    except:
        print "\n Cannot find variable KD"
    
    # Attempt to read rho_hv
    try:
        rhv = root.variables['RH']
        rhv_found = True
    except:
        print "\n Cannot find variable RH"
    
    # Attempt to read Vr

    try:
        Vr = root.variables['VEL']
        Vr_found = True
    except:
            print "\n Cannot find variable VEL"

    if(not Vr_found):
        try:
            Vr = root.variables['VR']
            Vr_found = True
        except:
            print "\n Cannot find variable VR"
    
    if(not Vr_found):    
        try:
            Vr = root.variables['VT']
            Vr_found = True
        except:
            print "\n Cannot find variable VT"
    
    if(not Vr_found):
        try:
            Vr = root.variables['VE']
            Vr_found = True
        except:
            print "\n Cannot find variable VE"
    
    if(not Vr_found):
        try:
            Vr = root.variables['VU']
            Vr_found = True
        except:
            print "\n Cannot find variable VU"

    height_array   = root.variables['HEIGHT']
    elev_array     = root.variables['EL']
    azim_array     = root.variables['AZ']
    time_array     = root.variables['TIME']
    date  = netcdf.chartostring(root.variables['start_date'][:])
    #print "start_date = "+str(date)
    time  = netcdf.chartostring(root.variables['start_time'][:])
    #print "start_time = "+str(time
    date2 = str.replace(date.tostring(),"/","-")

    date_string = "seconds since "+ date2[6:10]+"-"+date2[0:5] + " " + time.tostring()
    print
    print "COARDS string from file:  ",date_string

    sec_utime = utime(date_string)

    model_time = datetime(int(stime[0:4]),int(stime[4:6]),int(stime[6:8]),int(stime[8:10]),int(stime[10:12]),int(stime[12:14]))

    model_time_sec_from_oban_time = sec_utime.date2num(model_time)

    print "Difference in seconds of requested model time and oban start time = ",model_time_sec_from_oban_time

    if(dBZ_found):
        var_shape = dBZ.shape
    elif(Zdr_found):
        dBZ = False
        var_shape = Zdr.shape
    elif(Kdp_found):
        Zdr = False
        var_shape = Kdp.shape
    elif(rhv_found):
        Kdp = False
        var_shape = rhv.shape
    elif(Vr_found):
        rhv = False
        var_shape = Vr.shape
    else:
        print "\nNo radar data found in file! Exiting!"
        return x,y,False,dBZ_found,Zdr_found,Kdp_found,rhv_found,Vr_found
 
    # Now, search through the file and find the sweeps closest to the target sweep for the given variable
    # Store the data and times of each sweep in lists, which will be used later to attempt to match times from the
    # given model netCDF file as well as the target time

    times = []
    files = []

    model_time_diff = 100000 # Set initial model time difference to some large number

    foundtime = False

    for p in [var_shape[1]-1]:
    #for p in range(0,1):
        print "Getting data for pass "+str(p)
        #for k in N.arange(0, var_shape[2], options.skip):
        for k in N.arange(0, var_shape[2]):
            elev_mask = N.ma.masked_array(elev_array[0,k,:,:], elev_array[0,k,:,:] <= -32767.)
            mean_elev=elev_mask.mean()
            time_mask = N.ma.masked_array(time_array[0,k,:,:],   time_array[0,k,:,:]   <= -32767.)
            mean_time=time_mask.mean()
            if(abs(mean_elev-target_elev) <= 0.1):
                # Calculate difference in current sweep time from requested model time
                model_time_diff_old = model_time_diff
                model_time_diff_min = min(model_time_diff,model_time_diff_old)
                model_time_diff = abs(sec_utime.date2num(model_time)-mean_time)
                # Look for a sweep time that is within the time tolerance of the requested model time
                # and corresponds to the target elevation
                if(model_time_diff <= timetol and model_time_diff < model_time_diff_min):
                    foundtime = True
                    #print "Elevation = ",mean_elev
                    #print "Sweep Time = ",sec_utime.num2date(mean_time)
                    #height_swp = N.ma.masked_array(height_array[0,k,...],height_array[0,k,...] <= -32767.)
                    height_swp = height_array[0,k,...]
                    if(dBZ_found):
                        dBZ_swp = N.ma.masked_array(dBZ[0,p,k,...],dBZ[0,p,k,...] <= -32767.)
                    if(Zdr_found):
                        Zdr_swp = N.ma.masked_array(Zdr[0,p,k,...],Zdr[0,p,k,...] <= -32767.)
                    if(Kdp_found):
                        Kdp_swp = N.ma.masked_array(Kdp[0,p,k,...],Kdp[0,p,k,...] <= -32767.)
                    if(rhv_found):
                        rhv_swp = N.ma.masked_array(rhv[0,p,k,...],rhv[0,p,k,...] <= -32767.)
                    if(Vr_found):
                        Vr_swp = N.ma.masked_array(Vr[0,p,k,...],Vr[0,p,k,...] <= -32767.)
                        
                    print "Slice: %3d  Sweep elevation:  %4.2f  Sweep time %s " % (k,mean_elev, sec_utime.num2date(mean_time))
                    print "Requested model time %s " % sec_utime.num2date(model_time_sec_from_oban_time)
                    print "Time difference in seconds is %s " % model_time_diff
                    sweeptime = sec_utime.num2date(mean_time)
                #print "sec_utime, mean_time",sec_utime, mean_time

                #if mean_time != False and N.size(N.where(obanvar2.mask != False)) > 0:
                #files.append(plotoban(dd, vv, x, y, el.mean(), az, sec_utime.num2date(ti.mean()), pass_no = p+1, directory = options.dir))
    
    root.close()
    
    if(not foundtime):
        print "No matching times found in file! Exiting!"
        return False,x,y,False,False,False,False
    else:
        return sweeptime,x,y,height_swp,dBZ_swp,Zdr_swp,Kdp_swp,rhv_swp,Vr_swp
    
def Cresmn(x, y, obs, x0, y0, roi, missing=0.0):
    """ Returns a data value for the point
      Arguments: x/y/obs:  1D arrays of location
                 x0, y0:   point to analyze to
                 roi:      radius of influence
                 missing:  value to assign if no data, default = 0.0
    """
    # Create distance array
    
    dis = N.sqrt( (x-x0)**2 + (y-y0)**2 )
    #print dis.shape
    
    # Cut the size o n by choosing a smart threshold (2 km)
    
    indices = N.where(dis <= roi)
    #print indices
    
    # Check to see if there are any data ponts that are within search radius
    
    size = N.size(indices)
    
    if size != 0:
        R2    = roi**2.0
        w_sum = 0.0
        top   = 0.0
        #     print size, x, y, obs, x0, y0
    
    
        # go thru values w/in radius: calculate weight, mult by value and sum also sum weights
#        for n, value in N.ndenumerate(dis[indices]):
        for i,j in zip(indices[0],indices[1]):
            n = (i,j)
            value = dis[n]
            if(N.isfinite(obs[n])):
                #print value
                #print n
                rk2 = value**2.0
                wk = (R2-rk2) / (R2+rk2)
                top = top + wk*obs[n]
                w_sum = w_sum + wk
                #print n, value, R2, w_sum, top
        
        if w_sum < 0.01:
            return missing
        else:
            return top/w_sum
    
    else:
        return missing