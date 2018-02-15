#A module containing some functions to read ARPS history files and
#a few other auxilliary functions
#This is an updated module that uses the PyNio libraries to open HDF4 files instead of the older
#pyhdf

import numpy as N

try:
    import Nio
except ImportError:
    print "PyNIO not available, manipulating ARPS data will fail ..."

from mpl_toolkits.basemap import Basemap
import thermolib as thermo
import dualpara as dualpol
import DSDlib as dsd
import sys
import weave
from weave import converters
from scipy.special import gammaln
from datahandler import DataHandler
from matplotlib.colors import LinearSegmentedColormap
import datetime

from utils import log, warning, fatal

# New blue-to-red colormap kindly provided by Robin T.
cdict1 = {'red':  ((0.0, 0.0, 0.0),
                   (0.3, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.7, 0.9, 0.9),
                   (1.0, 0.4, 0.4)),

         'green': ((0.0, 0.0, 0.0),
                   (0.3, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.7, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.4, 0.4),
                   (0.3, 0.9, 0.9),
                   (0.5, 1.0, 1.0),
                   (0.7, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

class ARPSDataHandler(DataHandler):
    def __init__(self, base_dir, times, microphys=None):
        self._base_dir = base_dir
        self._time_list = times
        self._file = None
        self._file_name = ""
        self._grdbas_name = ""
        self._cur_time = -1
        self._microphys = microphys

        super(ARPSDataHandler, self).__init__('ARPS')
        return

    def setRun(self, file_base, member, time=None):
        if time is None:
            time = self._time_list[0]

        file_name = self.fileNameBuilder(self._base_dir + file_base, time, member)
        grdbas_name = self.fileNameBuilder(self._base_dir + file_base, None, member)
        if self._file_name != file_name:
            if self._file:
                self._file.close()
            if(time == -1): # Try to open grdbas file instead of time-dependent file
                self._file = Nio.open_file(grdbas_name,format='hdf') # Need to add support for netCDF later.
            else:
                self._file = Nio.open_file(file_name,format='hdf') # Need to add support for netCDF later.
            self._file_name = file_name
        return

    def setTime(self, time):
        self._cur_time = time
        try:
            self._timeindex = N.where(self._filetimelist==time)[0][0]
        except IndexError:
            fatal("Time '%d' not found!" % time) # TODO: Since ARPS stores times in separate files, need
                                                 # to do different error checking here.  Probably
                                                 # by globbing of directory and making sure file
                                                 # corresponding to chosen time exists...
        return self._timeobjlist[self._timeindex].strftime("%d %b %Y %H%M UTC")

    def loadGrid(self):
        log("Reading grid and map information from ARPS HDF file!")
        file_obj = None
        # Read grid variables
        nx,ny,nz,dx,dy,xe,ye,ze,xc,yc,zc,zeagl,zcagl = self._readarpsgrid()
        map_proj,ctrlat,ctrlon,trulat1,trulat2,trulon = self._readarpsmap()
        nxm = nx-1
        nym = ny-1
        nxe = nx
        nye = ny
        zc = zc.T
        ze = ze.T
        zcagl = zcagl.T
        zeagl = zeagl.T
        zc1d = N.mean(N.mean(zcagl,axis=1),axis=1) # Quick and dirty approximate heights AGL
                                                    # by taking mean of each column in domain.
        ze1d = N.mean(N.mean(zeagl,axis=1),axis=1)

        xcor, ycor = N.meshgrid(xe, ye)
        xe, _ = N.meshgrid(xe, yc)
        _, ye = N.meshgrid(xc, ye)
        xc, yc = N.meshgrid(xc, yc)

        # Create a basemap projection using the ARPS coordinates and map projection
        # Set up the map projection based on the info in the ARPS data file

        if(map_proj == 2): # Only support lambert conformal for now
            width = nx*dx
            height = ny*dy

            bgmap = Basemap(projection='lcc', width=width, height=height, lat_1=trulat1, lat_2=trulat2, \
                    lat_0=ctrlat, lon_0=ctrlon, resolution='h', area_thresh=10.,suppress_ticks=False)
        else:
            bgmap = None


        # Now compute x/y coordinates relative to radar location if desired (TODO).

#         if pc.plot_slice != 4:
#             xc_rel = xc
#             yc_rel = yc
#             xe_rel = xe
#             ye_rel = ye
#         else:
#             # Figure out x/y location of radar on ARPS grid
#             x_rad = bgmap(pc.rlon,pc.rlat)[0]
#             y_rad = bgmap(pc.rlon,pc.rlat)[1]
#
#             xc_rel = xc-x_rad
#             yc_rel = yc-y_rad
#             xe_rel = xe-x_rad
#             ye_rel = ye-y_rad
#
#         # A little trick to duplicate the xs and ys over the vertical coordinate without actually using any more memory.
#         xc_rel = N.lib.stride_tricks.as_strided(xc_rel, strides=((0, ) + xc_rel.strides), shape=((nz-1, ) + xc_rel.shape))
#         yc_rel = N.lib.stride_tricks.as_strided(yc_rel, strides=((0, ) + yc_rel.strides), shape=((nz-1, ) + yc_rel.shape))
#         xe_rel = N.lib.stride_tricks.as_strided(xe_rel, strides=((0, ) + xe_rel.strides), shape=((nz-1, ) + xe_rel.shape))
#         ye_rel = N.lib.stride_tricks.as_strided(ye_rel, strides=((0, ) + ye_rel.strides), shape=((nz-1, ) + ye_rel.shape))

        xc = N.lib.stride_tricks.as_strided(xc, strides=((0, ) + xc.strides), shape=((nz-1, ) + xc.shape))
        yc = N.lib.stride_tricks.as_strided(yc, strides=((0, ) + yc.strides), shape=((nz-1, ) + yc.shape))
        xe = N.lib.stride_tricks.as_strided(xe, strides=((0, ) + xe.strides), shape=((nz-1, ) + xe.shape))
        ye = N.lib.stride_tricks.as_strided(ye, strides=((0, ) + ye.strides), shape=((nz-1, ) + ye.shape))
        xcor = N.lib.stride_tricks.as_strided(xcor, strides=((0, ) + xcor.strides), shape=((nz-1, ) + xcor.shape))
        ycor = N.lib.stride_tricks.as_strided(ycor, strides=((0, ) + ycor.strides), shape=((nz-1, ) + ycor.shape))

        return xc, yc, zc, zc1d, xe, ye, ze, ze1d, bgmap

    def loadTimes(self):
        # First time setup only
        filetimes = N.array(self._time_list)    # for ARPS, each time is in a separate file, so we just use
                             # the requested time list and assume that these files
                             # are all present.  There is thus no timeindex.

        # Read in initial time from file
        day,year,month,hour,minute,second = self._readinitime()

        datetimestart = datetime.datetime(year,month,day,hour,minute,second)

        # Construct datetime objects
        filetimelist = []   # Will contain list of times in seconds since model start time in file.
        timeobjlist = []    # Will contain list of corresponding datetime objects
        timeindices = []

        for i,filetime in enumerate(filetimes):
            time_dt = datetime.timedelta(seconds=filetimes[i])
            timeobj = datetimestart+time_dt
            if i == 0:
                timeobj0 = timeobj
                inittime = timeobj
            filetimelist.append(time_dt.seconds)
            timeobjlist.append(timeobj)
            timeindices.append(i)

        self._filetimelist = N.array(filetimelist)
        self._timeobjlist = timeobjlist
        self._inittime = inittime
        return

    def loadMicrophysics(self):
        # Right now, only read MY or ZVD scheme.  Support for other schemes to be
        # added later
        if self._microphys[0:2] == 'WR':
            log("Reading microphysics information for "+self._microphys+" scheme.")
            micro, consts = self._readarpsmicroWR()
            micro['qs'] = N.zeros_like(micro['qr'])
            micro['qg']= N.zeros_like(micro['qr'])
            micro['qh'] = N.zeros_like(micro['qr'])
            micro['ntr'] = N.zeros_like(micro['qr'])
            micro['nts'] = N.zeros_like(micro['qr'])
            micro['ntg'] = N.zeros_like(micro['qr'])
            micro['nth'] = N.zeros_like(micro['qr'])
            micro['zr'] = N.zeros_like(micro['qr'])
            micro['zs'] = N.zeros_like(micro['qr'])
            micro['zg'] = N.zeros_like(micro['qr'])
            micro['zh'] = N.zeros_like(micro['qr'])
            micro['alphar'] = N.zeros_like(micro['qr'])
            micro['alphas'] = N.zeros_like(micro['qr'])
            micro['alphag'] = N.zeros_like(micro['qr'])
            micro['alphah'] = N.zeros_like(micro['qr'])
            micro['rhos'] = consts['rhoscst']
            micro['rhog'] = consts['rhogcst']
            micro['rhoh'] = consts['rhohcst']
            consts['MPflg'] = 5
            consts['microphys'] = 'WR'
        elif self._microphys[0:3] == 'LFO':
            log("Reading microphysics information for "+self._microphys+" scheme.")
            micro, consts = self._readarpsmicroLFO()
            micro['qg']= N.zeros_like(micro['qr'])
            micro['ntg'] = N.zeros_like(micro['qr'])
            micro['zg'] = N.zeros_like(micro['qr'])
            micro['rhos'] = consts['rhoscst']
            micro['rhog'] = consts['rhogcst']
            micro['rhoh'] = consts['rhohcst']
            consts['MPflg'] = 5
            consts['microphys'] = 'LFO'
        elif self._microphys[0:2] == 'MY':
            log("Reading microphysics information for "+self._microphys+" scheme.")
            micro, consts = self._readarpsmicroMY()
            consts['MPflg'] = 0
            consts['microphys'] = 'MY'
        elif self._microphys[0:3] == 'ZVD':
            micro, consts = self._readarpsmicroZVD()
            micro['rhos'] = consts['rhoscst']
            consts['imurain'] = 1
            consts['MPflg'] = 2
            consts['microphys'] = 'ZVD'
        else:
            fatal("This option is not supported for now! Exiting!")

        if 'MPflg' not in consts:
            consts['MPflg'] = 0

        if 'imurain' not in consts:
            consts['imurain'] = 1

        # Also initialize water fraction on snow, graupel, and hail
        micro['qsw'] = N.asfortranarray(N.zeros_like(micro['qs']))
        micro['qgw'] = N.asfortranarray(N.zeros_like(micro['qg']))
        micro['qhw'] = N.asfortranarray(N.zeros_like(micro['qh']))

        return micro, consts

    def loadMeanWind(self):
        self._abstract('loadMeanWind')
        return

    def loadHorizWind(self, average=True):
        u = self._hdfreadvar3d('u',stag='u')
        v = self._hdfreadvar3d('v',stag='v')

        #Average u,v to the scalar points

        if average:
            uc = 0.5*(u[:-1,:,:]+u[1:,:,:])
            vc = 0.5*(v[:,:-1,:]+v[:,1:,:])
        else:
            uc = u
            vc = v

        return uc,vc
        #return N.transpose(uc,(1,0,2)),N.transpose(vc,(1,0,2))

    def loadVertWind(self, average=True):
        w = self._hdfreadvar3d('w',stag='w')
        if average:
            wc = 0.5*(w[:,:,:-1]+w[:,:,1:])
        else:
            wc = w
        return wc.T

    def loadModelReflectivity(self):
        return

    def getTimeObjs(self):
        return self._inittime, self._timeobjlist

    def fileNameBuilder(self, file_base, time, member):
        if(time is None):
            file_base = file_base+'.hdfgrdbas'
        else:
            file_base = file_base+'.hdf'+'%06d' % time
        return file_base

    def _readinitime(self):
        """This function reads the initial model time from an ARPS HDF file"""
        day = self._file.day[0]
        year = self._file.year[0]
        month = self._file.month[0]
        hour = self._file.hour[0]
        minute = self._file.minute[0]
        second = self._file.second[0]

        return day,year,month,hour,minute,second

    def _readarpsgrid(self):
        """This function reads grid information from an ARPS HDF file"""

        nx = self._file.nx[0]
        ny = self._file.ny[0]
        nz = self._file.nz[0]

        xe=self._file.variables['x']
        xe=xe[:]
        dx=xe[1]-xe[0]

        xc=N.empty_like(xe)
        xc[:-1] = 0.5*(xe[:-1]+xe[1:])
        xc=xc[:-1]
        ye=self._file.variables['y']
        ye=ye[:]
        dy=ye[1]-ye[0]

        yc=N.empty_like(ye)
        yc[:-1] = 0.5*(ye[:-1]+ye[1:])
        yc=yc[:-1]

        ze=self._hdfreadvar3d('zp',stag='w')

        zc=N.empty_like(ze)
        zc[:,:,:-1] = 0.5*(ze[:,:,:-1]+ze[:,:,1:])
        zcagl=N.empty_like(zc)
        zeagl=N.empty_like(ze)
        for k in range(nz):
            zcagl[...,k] = zc[...,k]-ze[:,:,1]
            zeagl[...,k] = ze[...,k]-ze[:,:,1]

        zc=zc[:,:,:-1]
        zcagl=zcagl[:,:,:-1]

        return nx,ny,nz,dx,dy,xe,ye,ze,xc,yc,zc,zeagl,zcagl

    def _readarpsmap(self):
        """This function reads the arps map projection information from the given file"""
        map_proj=self._file.mapproj[0]
        ctrlat=self._file.ctrlat[0]
        ctrlon=self._file.ctrlon[0]
        trulat1=self._file.trulat1[0]
        trulat2=self._file.trulat2[0]
        trulon=self._file.trulon[0]
        return map_proj,ctrlat,ctrlon,trulat1,trulat2,trulon

    def _readarpsmicroparam(self):
        """This function reads the arps microphsyics parameter information"""

        n0rain=self._file.n0rain
        n0snow=self._file.n0snow
        n0hail=self._file.n0hail
        rhosnow=self._file.rhosnow
        rhohail=self._file.rhohail

        try:
            ntcloud=self._file.ntcloud
            n0grpl=self._file.n0grpl
            rhoice=self._file.rhoice
            rhogrpl=self._file.rhogrpl
            alpharain=self._file.alpharain
            alphaice=self._file.alphaice
            alphasnow=self._file.alphasnow
            alphagrpl=self._file.alphagrpl
            alphahail=self._file.alphahail

        # For files that were converted from WRF. This is a temporary fix. (Y. Jung 3/23/2017)
            if rhoice == 0.0:
               rhoice = 500.
            if rhosnow == 0.0:
               rhosnow = 100.
            if rhogrpl == 0.0:
               rhogrpl = 400.
            if rhohail == 0.0:
               rhohail = 913.
        except:
            ntcloud=-1.0
            n0grpl=n0hail
            rhoice=rhohail
            rhogrpl=rhohail
            alpharain=0.0
            alphaice=0.0
            alphasnow=0.0
            alphagrpl=0.0
            alphahail=0.0

        return n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,ntcloud,alpharain,alphaice,alphasnow,alphagrpl,alphahail

    def readarpsmicroWR(self):
        """Reads in several microphysical variables for the ARPS Kessler warm-rain scheme and derives others from an ARPS
           HDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        verbose = True
        fortran = True

        micro, consts = {}, {}

        # Calculate moist air density

        pt = self._hdfreadvar3d('pt',stag='s')
        qv = self._hdfreadvar3d('qv',stag='s')
        p = self._hdfreadvar3d('p',stag='s')
        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)
        # Read in microphysics parameter information from file

        n0rain,n0snow,n0grpl,n0hail,rhoicst,rhoscst,rhogcst, \
        rhohcst,ntcloud,alpharcst,alphaicst,alphascst,alphagcst,alphahcst = self._readarpsmicroparam()
        if(verbose):
            print 'n0rain = ',n0rain
            print 'n0snow = ',n0snow
            print 'n0grpl = ',n0grpl
            print 'n0hail = ',n0hail
            print 'rhoice = ',rhoicst
            print 'rhosnow = ',rhoscst
            print 'rhogrpl = ',rhogcst
            print 'rhohail = ',rhohcst
            print 'ntcloud = ',ntcloud
            print 'alpharain = ',alpharcst
            print 'alphaice = ',alphaicst
            print 'alphasnow = ',alphascst
            print 'alphagrpl = ',alphagcst
            print 'alphahail = ',alphahcst

        consts['n0rain'] = n0rain
        consts['n0snow'] = n0snow
        consts['n0grpl'] = n0grpl
        consts['n0hail'] = n0hail
        consts['rhoicst'] = rhoicst
        consts['rhoscst'] = rhoscst
        consts['rhogcst'] = rhogcst
        consts['rhohcst'] = rhohcst
        consts['alpharcst'] = alpharcst
        consts['alphascst'] = alphascst
        consts['alphagcst'] = alphagcst
        consts['alphahcst'] = alphahcst

        micro['rhor'] = 1000.*N.ones_like(pt)
        micro['rhos'] = rhoscst*N.ones_like(pt)
        micro['rhog'] = rhogcst*N.ones_like(pt)
        micro['rhoh'] = rhohcst*N.ones_like(pt)

        # Read in mixing ratios

        micro['qr'] = self._hdfreadvar3d('qr',stag='s')
        consts['graupel_on'] = 0
        consts['hail_on'] = 0

        for var, data in micro.iteritems():
            if type(data) == N.ndarray:
                micro[var] = N.asfortranarray(data)

        return micro,consts

    def _readarpsmicroLFO(self):
        """Reads in several microphysical variables for the ARPS LFO scheme and derives others from an ARPS
           HDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        verbose = True
        fortran = True

        micro, consts = {}, {}

        # Calculate moist air density

        pt = self._hdfreadvar3d('pt',stag='s')
        qv = self._hdfreadvar3d('qv',stag='s')
        p = self._hdfreadvar3d('p',stag='s')
        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)
        # Read in microphysics parameter information from file

        n0rain,n0snow,n0grpl,n0hail,rhoicst,rhoscst,rhogcst, \
        rhohcst,ntcloud,alpharcst,alphaicst,alphascst,alphagcst,alphahcst = self._readarpsmicroparam()
        if(verbose):
            print 'n0rain = ',n0rain
            print 'n0snow = ',n0snow
            print 'n0grpl = ',n0grpl
            print 'n0hail = ',n0hail
            print 'rhoice = ',rhoicst
            print 'rhosnow = ',rhoscst
            print 'rhogrpl = ',rhogcst
            print 'rhohail = ',rhohcst
            print 'ntcloud = ',ntcloud
            print 'alpharain = ',alpharcst
            print 'alphaice = ',alphaicst
            print 'alphasnow = ',alphascst
            print 'alphagrpl = ',alphagcst
            print 'alphahail = ',alphahcst

        consts['n0rain'] = n0rain
        consts['n0snow'] = n0snow
        consts['n0grpl'] = n0grpl
        consts['n0hail'] = n0hail
        consts['rhoicst'] = rhoicst
        consts['rhoscst'] = rhoscst
        consts['rhogcst'] = rhogcst
        consts['rhohcst'] = rhohcst
        consts['alpharcst'] = alpharcst
        consts['alphascst'] = alphascst
        consts['alphagcst'] = alphagcst
        consts['alphahcst'] = alphahcst

        micro['rhor'] = 1000.*N.ones_like(pt)
        micro['rhos'] = rhoscst*N.ones_like(pt)
        micro['rhog'] = rhogcst*N.ones_like(pt)
        micro['rhoh'] = rhohcst*N.ones_like(pt)

        # Read in mixing ratios

        micro['qr'] = self._hdfreadvar3d('qr',stag='s')
        micro['qs'] = self._hdfreadvar3d('qs',stag='s')
        micro['qh'] = self._hdfreadvar3d('qh',stag='s')
        consts['graupel_on'] = 0
        consts['hail_on'] = 1

        micro['ntr'] = N.zeros_like(pt)
        micro['nts'] = N.zeros_like(pt)
        micro['nth'] = N.zeros_like(pt)

        micro['zr'] = N.zeros_like(pt)
        micro['zs'] = N.zeros_like(pt)
        micro['zh'] = N.zeros_like(pt)

        micro['alphar'] = alpharcst*N.ones_like(pt)
        micro['alphas'] = alphascst*N.ones_like(pt)
        micro['alphah'] = alphahcst*N.ones_like(pt)

        for var, data in micro.iteritems():
            if type(data) == N.ndarray:
                micro[var] = N.asfortranarray(data)

        return micro,consts

    def _readarpsmicroMY(self):
        """Reads in several microphysical variables for the MY scheme and derives others from an ARPS
           HDF file.  Also checks if the sixth moments are present and computes the shape parameter as appropriate.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""

        verbose = True
        fortran = True

        micro, consts = {}, {}

        # Constants in diagnostic alpha relations (for mphyopt = 10)

        c1r = 19.0
        c2r = 0.6
        c3r = 1.8
        c4r = 17.0
        c1i = 12.0
        c2i = 0.7
        c3i = 1.7
        c4i = 11.0
        c1s = 4.5
        c2s = 0.5
        c3s = 5.0
        c4s = 5.5
        c1g = 5.5
        c2g = 0.7
        c3g = 4.5
        c4g = 8.5
        c1h = 3.7
        c2h = 0.3
        c3h = 9.0
        c4h = 6.5
        c5h = 1.0
        c6h = 6.5

        # Calculate moist air density

        pt = self._hdfreadvar3d('pt',stag='s')
        qv = self._hdfreadvar3d('qv',stag='s')
        p = self._hdfreadvar3d('p',stag='s')
        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)
        # Read in microphysics parameter information from file

        n0rain,n0snow,n0grpl,n0hail,rhoicst,rhoscst,rhogcst, \
        rhohcst,ntcloud,alpharcst,alphaicst,alphascst,alphagcst,alphahcst = self._readarpsmicroparam()

        if(verbose):
            print 'n0rain = ',n0rain
            print 'n0snow = ',n0snow
            print 'n0grpl = ',n0grpl
            print 'n0hail = ',n0hail
            print 'rhoice = ',rhoicst
            print 'rhosnow = ',rhoscst
            print 'rhogrpl = ',rhogcst
            print 'rhohail = ',rhohcst
            print 'ntcloud = ',ntcloud
            print 'alpharain = ',alpharcst
            print 'alphaice = ',alphaicst
            print 'alphasnow = ',alphascst
            print 'alphagrpl = ',alphagcst
            print 'alphahail = ',alphahcst

        consts['n0rain'] = n0rain
        consts['n0snow'] = n0snow
        consts['n0grpl'] = n0grpl
        consts['n0hail'] = n0hail
        consts['rhoicst'] = rhoicst
        consts['rhoscst'] = rhoscst
        consts['rhogcst'] = rhogcst
        consts['rhohcst'] = rhohcst
        consts['alpharcst'] = alpharcst
        consts['alphascst'] = alphascst
        consts['alphagcst'] = alphagcst
        consts['alphahcst'] = alphahcst

        micro['rhor'] = 1000.*N.ones_like(pt)
        micro['rhos'] = rhoscst*N.ones_like(pt)
        micro['rhog'] = rhogcst*N.ones_like(pt)
        micro['rhoh'] = rhohcst*N.ones_like(pt)

        # Read in mixing ratios

        micro['qr'] = self._hdfreadvar3d('qr',stag='s')
        micro['qs'] = self._hdfreadvar3d('qs',stag='s')
        micro['qg'] = self._hdfreadvar3d('qg',stag='s')
        micro['qh'] = self._hdfreadvar3d('qh',stag='s')
        consts['graupel_on'] = 1  # TODO: This should be present in the file now, so should update this to read from file...
        consts['hail_on'] = 1 # TODO: This should be present in the file now, so should update this to read from file...

        # Try to read in number concentration variables

        try:
            micro['ntr'] = self._hdfreadvar3d('nr',stag='s')
            micro['nts'] = self._hdfreadvar3d('ns',stag='s')
            micro['ntg'] = self._hdfreadvar3d('ng',stag='s')
            micro['nth'] = self._hdfreadvar3d('nh',stag='s')
        except:
            micro['ntr'] = N.zeros_like(pt)
            micro['nts'] = N.zeros_like(pt)
            micro['ntg'] = N.zeros_like(pt)
            micro['nth'] = N.zeros_like(pt)

        # Now try to read 6th moment (reflectivity) variables
        # and compute the corresponding shape parameters
        # If they don't exist, use the constant alpha values
        # or diagnose them (if microphys == 'MY2DA')

        try:
            micro['zr'] = self._hdfreadvar3d('zr',stag='s')
            micro['zs'] = self._hdfreadvar3d('zs',stag='s')
            micro['zg'] = self._hdfreadvar3d('zg',stag='s')
            micro['zh'] = self._hdfreadvar3d('zh',stag='s')
            mu = 1./3.
            micro['alphar'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qr'],micro['ntr'],micro['zr'],micro['rhor'])
            micro['alphas'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qs'],micro['nts'],micro['zs'],micro['rhos'])
            micro['alphag'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qg'],micro['ntg'],micro['zg'],micro['rhog'])
            micro['alphah'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qh'],micro['nth'],micro['zh'],micro['rhoh'])
            print "Found Z arrays! Computed shape parameters!"
        except:
            if(self._microphys == 'MY2DA'):
                print "Did not find Z arrays! Computing shape parameters for MY2DA."
                dmr = (rhoa*micro['qr']/((N.pi/6.)*micro['rhor']*micro['nr']))**(1./3.)
                micro['alphar'] = c1r*N.tanh(c2r*(dmr-c3r))+c4r
                dmi = (rhoa*micro['qi']/((N.pi/6.)*micro['rhoi']*micro['ni']))**(1./3.)
                micro['alphai'] = c1i*N.tanh(c2i*(dmr-c3i))+c4i
                dms = (rhoa*micro['qs']/((N.pi/6.)*micro['rhos']*micro['ns']))**(1./3.)
                micro['alphas'] = c1s*N.tanh(c2s*(dms-c3s))+c4s
                dmg = (rhoa*micro['qg']/((N.pi/6.)*micro['rhog']*micro['ng']))**(1./3.)
                micro['alphag'] = c1g*N.tanh(c2g*(dmg-c3g))+c4g
                dmh = (rhoa*micro['qh']/((N.pi/6.)*micro['rhoh']*micro['nh']))**(1./3.)
                micro['alphah'] = N.where(dmh < 8e-3, c1h*N.tanh(c2h*(dmh-c3h))+c4h, c5h*dmh-c6h)
            else:
                print "Did not find Z arrays! Assuming 2-moment with constant shape parameter."
            micro['zr'] = N.zeros_like(pt)
            micro['zs'] = N.zeros_like(pt)
            micro['zg'] = N.zeros_like(pt)
            micro['zh'] = N.zeros_like(pt)
            micro['alphar'] = alpharcst*N.ones_like(pt)
            micro['alphas'] = alphascst*N.ones_like(pt)
            micro['alphag'] = alphagcst*N.ones_like(pt)
            micro['alphah'] = alphahcst*N.ones_like(pt)

        for var, data in micro.iteritems():
            if type(data) == N.ndarray:
                micro[var] = N.asfortranarray(data)

        return micro, consts

    def _readarpsmicroZVD(self):
        """Reads in several microphysical variables for the ZVD scheme and derives others from an ARPS hdf file
           Determines if graupel and hail are present and returns the appropriate flags.  Also checks if sixth moments are present
           and computes the shape parameter as appropriate.Note, for now qi and qc are ignored
           since they aren't used in the polarimetric emulator, but obviously need to add these
           for use by other applications."""

        verbose = True
        fortran = True

        micro, consts = {}, {}

        # Calculate moist air density

        pt = self._hdfreadvar3d('pt',stag='s')
        qv = self._hdfreadvar3d('qv',stag='s')
        p = self._hdfreadvar3d('p',stag='s')
        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)
        # Read in microphysics parameter information from file

        n0rain,n0snow,n0grpl,n0hail,rhoicst,rhoscst,rhogcst, \
        rhohcst,ntcloud,alpharcst,alphaicst,alphascst,alphagcst,alphahcst = self._readarpsmicroparam()

        if(verbose):
            print 'n0rain = ',n0rain
            print 'n0snow = ',n0snow
            print 'n0grpl = ',n0grpl
            print 'n0hail = ',n0hail
            print 'rhoice = ',rhoicst
            print 'rhosnow = ',rhoscst
            print 'rhogrpl = ',rhogcst
            print 'rhohail = ',rhohcst
            print 'ntcloud = ',ntcloud
            print 'alpharain = ',alpharcst
            print 'alphaice = ',alphaicst
            print 'alphasnow = ',alphascst
            print 'alphagrpl = ',alphagcst
            print 'alphahail = ',alphahcst

        consts['n0rain'] = n0rain
        consts['n0snow'] = n0snow
        consts['n0grpl'] = n0grpl
        consts['n0hail'] = n0hail
        consts['rhoicst'] = rhoicst
        consts['rhoscst'] = rhoscst
        consts['rhogcst'] = rhogcst
        consts['rhohcst'] = rhohcst
        consts['alpharcst'] = alpharcst
        consts['alphascst'] = alphascst
        consts['alphagcst'] = alphagcst
        consts['alphahcst'] = alphahcst

        micro['rhor'] = 1000.*N.ones_like(pt)
        micro['rhos'] = rhoscst*N.ones_like(pt)

        micro['qr'] = self._hdfreadvar3d('qr',stag='s')
        micro['qs'] = self._hdfreadvar3d('qs',stag='s')
        micro['qg'] = self._hdfreadvar3d('qg',stag='s')
        micro['qh'] = self._hdfreadvar3d('qh',stag='s')
        consts['graupel_on'] = 1  # TODO: This should be present in the file now, so should update this to read from file...
        consts['hail_on'] = 1 # TODO: This should be present in the file now, so should update this to read from file...

        # Try to read in graupel and hail volumes and compute bulk densities

        try:
            micro['vg'] = self._hdfreadvar3d('vg',stag='s')
            micro['rhog'] = micro['qg']*rhoa/micro['vg']
        except:
            micro['rhog'] = rhogcst*N.ones_like(pt)

        try:
            micro['vh'] = self._hdfreadvar3d('vh',stag='s')
            micro['rhoh'] = micro['qh']*rhoa/micro['vh']
        except:
            micro['rhoh'] = rhohcst*N.ones_like(pt)

        # Try to read in number concentration variables

        try:
            micro['ntr'] = self._hdfreadvar3d('nr',stag='s')
            micro['nts'] = self._hdfreadvar3d('ns',stag='s')
            micro['ntg'] = self._hdfreadvar3d('ng',stag='s')
            micro['nth'] = self._hdfreadvar3d('nh',stag='s')
        except:
            micro['ntr'] = N.zeros_like(pt)
            micro['nts'] = N.zeros_like(pt)
            micro['ntg'] = N.zeros_like(pt)
            micro['nth'] = N.zeros_like(pt)

        # Now try to read 6th moment (reflectivity) variables
        # and compute the corresponding shape parameters
        # If they don't exist, use the constant alpha values
        # or diagnose them (if microphys == 'MY2DA')

        try:
            micro['zr'] = self._hdfreadvar3d('zr',stag='s')
            micro['zg'] = self._hdfreadvar3d('zg',stag='s')
            micro['zh'] = self._hdfreadvar3d('zh',stag='s')
            mu = 1./3.
            micro['alphar'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qr'],micro['ntr'],micro['zr'],micro['rhor'])
            micro['alphag'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qg'],micro['ntg'],micro['zg'],micro['rhog'])
            micro['alphah'] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['qh'],micro['nth'],micro['zh'],micro['rhoh'])
            print "Found Z arrays! Computed shape parameters!"
        except:
            print "Did not find Z arrays! Assuming 2-moment with constant shape parameter."
            micro['zr'] = N.zeros_like(pt)
            micro['zg'] = N.zeros_like(pt)
            micro['zh'] = N.zeros_like(pt)
            micro['alphar'] = alpharcst*N.ones_like(pt)
            micro['alphag'] = alphagcst*N.ones_like(pt)
            micro['alphah'] = alphahcst*N.ones_like(pt)

        # Snow is never 3-moment in ZVD
        micro['zs'] = N.zeros_like(qs)
        micro['alphas'] = alphascst*N.ones_like(pt)

        for var, data in micro.iteritems():
            if type(data) == N.ndarray:
                micro[var] = N.asfortranarray(data)

        return micro, consts

    def _hdfreadvar3d(self,varname,stag=None):
        """Read and return a 3d variable from an ARPS HDF file"""
        varname = varname.rstrip('_')
        vards=self._file.variables[varname]
        #Check to see if the variable is packed to 16 bit integers
        try:
            packed=vards.packed16
        except:
            packed=0
        #extract the variable to a numpy array and convert to float
        var=vards[:].astype('float32')
        #Put in ijk indexing order
        var=var.swapaxes(0,1)
        var=var.swapaxes(1,2)
        var=var.swapaxes(0,1)
        #Unpack the data from 16-bit integers if needed
        if (packed == 1):
            #Get the max and min values for each level
            hmax=N.array(vards.max[:])
            hmin=N.array(vards.min[:])
            #Perform the unpacking
            for k in xrange(0,len(var[0,0,:])):
                scalef = (hmax[k]-hmin[k])/65534.0
                var[:,:,k] = scalef*(var[:,:,k] + 32767)+hmin[k]
        if(stag=='s'):
            var=var[:-1,:-1,:-1]
        elif(stag == 'u'):
            var=var[:,:-1,:-1]
        elif(stag == 'v'):
            var=var[:-1,:,:-1]
        elif(stag == 'w'):
            var=var[:-1,:-1,:]
        return var
