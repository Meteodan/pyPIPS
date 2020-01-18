# commasmodule.py
# A collection of functions related to reading and manipulating commas netCDF files

import netCDF4 as netcdf
import numpy as N
from mpl_toolkits.basemap import Basemap
from . import thermolib as thermo
from . import DSDlib as dsd
from scipy.special import gammaln
from .datahandler import DataHandler
import datetime

from .utils import log, fatal

# Gravitational acceleration
g = 9.80665   # m s^-2


class WRFDataHandler(DataHandler):
    def __init__(self, base_dir, times):
        self._base_dir = base_dir
        self._time_list = times
        self._file = None
        self._file_name = ""
        self._cur_time = -1

        super(WRFDataHandler, self).__init__('WRF')
        return

    def setRun(self, file_base, member, time=None):
        if time is None:
            time = self._time_list[0]

        file_name = self.fileNameBuilder(self._base_dir + file_base, time, member)
        if self._file_name != file_name:
            if self._file:
                self._file.close()

            self._file = netcdf.Dataset(file_name)
            self._file_name = file_name
        return

    def setTime(self, time):
        self._cur_time = time
        try:
            self._timeindex = N.where(self._filetimelist == time)[0][0]
        except IndexError:
            fatal("Time '%d' not found in WRF input file." % time)
        return self._timeobjlist[self._timeindex].strftime("%d %b %Y %H%M UTC")

    def loadGrid(self):
        log("Reading grid and map information from WRF netCDF file!")
        # Read grid variables
        xlon, xlon_u, xlon_v, xlat, xlat_u, xlat_v, zc, ze = self._readwrfgrid()
        map_proj, ctrlat, ctrlon, trulat1, trulat2, trulon = self._readwrfmap()

        zc1d = N.mean(N.mean(zc, axis=1), axis=1)  # Quick and dirty approximate heights AGL
        # by taking mean of each column in domain.
        ze1d = N.mean(N.mean(ze, axis=1), axis=1)

        nxm = N.size(xlon, axis=1)
        nym = N.size(xlat, axis=0)
        nxe = N.size(xlon_u, axis=1)
        nye = N.size(xlat_v, axis=0)
        nzm = N.size(zc, axis=0)
        # nze = N.size(ze, axis=0)

        dx = self._file.DX
        dy = self._file.DY

        if(map_proj == 0):  # No map projection (probably idealized run)
            ctrlat = 36.0   # Are these even needed?
            ctrlon = -100.0
            trulat1 = 30.0
            trulat2 = 60.0

            # Create idealized grids
            xc, yc = N.meshgrid(N.arange(nxm) * dx, N.arange(nym) * dy)
            # Note, grid edges are offset by dx/2.
            xe, dummy = N.meshgrid(N.arange(nxe) * dx - dx / 2., N.arange(nym) * dy)
            dummy, ye = N.meshgrid(N.arange(nxm) * dx, N.arange(nye) * dy - dy / 2.)
            # Also compute the corner points for pcolormesh
            xcor, ycor = N.meshgrid(N.arange(nxe) * dx - dx / 2., N.arange(nye) * dy - dy / 2.)

            bgmap = None

        else:   # Assume lambert conformal for now
            # Create a basemap projection using the WRF coordinates and map projection
            # Set up the map projection based on the info in the WRF netCDF file

            width = nxe * dx
            height = nye * dy

            bgmap = Basemap(projection='lcc', width=width, height=height, lat_1=trulat1,
                            lat_2=trulat2, lat_0=ctrlat, lon_0=ctrlon, resolution='h',
                            area_thresh=10., suppress_ticks=False)

            # Now find the x/y coordinates of the grid points.
            xc, yc = bgmap(xlon, xlat)
            # xe, ye = bgmap(xlon_u,xlat_u) # I don't think this is correct...
            # Need to come back to this (DTD: 06/10/2014): pcolormesh expects grid corners
            xe, yxe = bgmap(xlon_u, xlat_u)
            xye, ye = bgmap(xlon_v, xlat_v)
            # Now compute x/y coordinates relative to radar location if desired.

#           # A little trick to duplicate the xs and ys over the vertical coordinate without
#           # actually using any more memory.
#           xc_rel = N.lib.stride_tricks.as_strided(xc_rel, strides=((0, ) + xc_rel.strides),
#                                                   shape=((nzm, ) + xc_rel.shape))
#           yc_rel = N.lib.stride_tricks.as_strided(yc_rel, strides=((0, ) + yc_rel.strides),
#                                                   shape=((nzm, ) + yc_rel.shape))
#           xe_rel = N.lib.stride_tricks.as_strided(xe_rel, strides=((0, ) + xe_rel.strides),
#                                                   shape=((nzm, ) + xe_rel.shape))
#           ye_rel = N.lib.stride_tricks.as_strided(ye_rel, strides=((0, ) + ye_rel.strides),
#                                                   shape=((nzm, ) + ye_rel.shape))

        xc = N.lib.stride_tricks.as_strided(xc, strides=(
            (0, ) + xc.strides), shape=((nzm, ) + xc.shape))
        yc = N.lib.stride_tricks.as_strided(yc, strides=(
            (0, ) + yc.strides), shape=((nzm, ) + yc.shape))
        xe = N.lib.stride_tricks.as_strided(xe, strides=(
            (0, ) + xe.strides), shape=((nzm, ) + xe.shape))
        ye = N.lib.stride_tricks.as_strided(ye, strides=(
            (0, ) + ye.strides), shape=((nzm, ) + ye.shape))

        return xc, yc, zc, zc1d, xe, ye, ze, ze1d, bgmap

    def loadTimes(self):
        # First time setup only
        filetimes = self._file.variables['Times'][:]
        # Convert times to seconds since start time in file
        # First construct datetime objects
        filetimelist = []   # Will contain list of times in seconds since model start time in file.
        timeobjlist = []    # Will contain list of corresponding datetime objects
        for i, filetime in enumerate(filetimes):
            timeobj = datetime.datetime.strptime("".join(filetime), "%Y-%m-%d_%H:%M:%S")

            if i == 0:
                timeobj0 = timeobj
            # timedelta object representing difference in time from start time
            time_dt = timeobj - timeobj0
            filetimelist.append(time_dt.seconds)
            timeobjlist.append(timeobj)

        self._filetimelist = N.array(filetimelist)
        self._timeobjlist = timeobjlist
        return

    def loadMicrophysics(self):
        # Support for other schemes is easily added later.
        if(self._file.MP_PHYSICS == 16):
            log("Reading microphysics information for WDM6 in WRF")
            micro, consts = self._readwrfmicroWDM6()
            consts['microphys'] = 'WDM6'
            # WDM6 doesn't have a separate hail category, but set these arrays to zero
            # since the code below requires they be present anyway.
            micro['qh'] = N.asfortranarray(N.zeros_like(micro['qg']))
            micro['nth'] = N.asfortranarray(N.zeros_like(micro['ntg']))
            micro['alphah'] = micro['alphag']
            micro['rhoh'] = micro['rhog']
            consts['rhohcst'] = consts['rhogcst']
            consts['n0hail'] = consts['n0grpl']
            consts['MPflg'] = 5
        elif(self._file.MP_PHYSICS == 9):
            log("Reading microphysics information for MY in WRF")
            micro, consts = self._readwrfmicroMY()
            consts['microphys'] = 'MY'
        elif(self._file.MP_PHYSICS == 18):
            log("Reading microphysics information for NSSL in WRF")
            micro, consts = self._readwrfmicroNSSL()
            consts['microphys'] = 'NSSL'
            wrf_ver_str = self._file.TITLE.split()
            # wrf_ver = wrf_ver_str[3][1:]
            if wrf_ver_str[3][1:] > '3.5.1':
                consts['MPflg'] = 2
                micro['alphar'] = N.zeros_like(micro['qr'])
            else:
                consts['MPflg'] = 1
                micro['alphar'] = -0.4 * N.ones_like(micro['qr'])
            print('MPflg is ', consts['MPflg'])
        elif(self._file.MP_PHYSICS == 10):
            log("Reading microphysics information for Morrison in WRF")
            micro, consts = self._readwrfmicroMorrison()
            consts['microphys'] = 'Morrison'
            micro['qh'] = N.asfortranarray(N.zeros_like(micro['qg']))
            micro['nth'] = N.asfortranarray(N.zeros_like(micro['ntg']))
            consts['n0hail'] = consts['n0grpl']
            micro['alphah'] = micro['alphag']
            micro['rhoh'] = micro['rhog']
        elif(self._file.MP_PHYSICS == 8):
            log("Reading microphysics information for Thompson in WRF")
            micro, consts = self._readwrfmicroThomp()
            consts['microphys'] = 'Thompson'
            micro['qh'] = N.asfortranarray(N.zeros_like(micro['qg']))
            micro['nth'] = N.asfortranarray(N.zeros_like(micro['ntg']))
            consts['n0hail'] = consts['n0grpl']
            micro['alphah'] = micro['alphag']
            micro['rhoh'] = micro['rhog']
            consts['MPflg'] = 4
        elif(self._file.MP_PHYSICS == 6):
            log("Reading microphysics information for WSM6 in WRF")
            micro, consts = self._readwrfmicroWSM6()
            consts['microphys'] = 'WSM6'
            micro['qh'] = N.asfortranarray(N.zeros_like(micro['qg']))
            micro['nth'] = N.asfortranarray(N.zeros_like(micro['ntg']))
            consts['n0hail'] = consts['n0grpl']
            micro['alphah'] = micro['alphag']
            micro['rhoh'] = micro['rhog']
            consts['MPflg'] = 5
        elif(self._file.MP_PHYSICS == 2):
            log("Reading microphysics information for LIN in WRF")
            micro, consts = self._readwrfmicroLIN()
            consts['microphys'] = 'LIN'
            micro['qg'] = N.asfortranarray(N.zeros_like(micro['qh']))
            micro['ntg'] = N.asfortranarray(N.zeros_like(micro['nth']))
            consts['n0grpl'] = consts['n0hail']
            micro['alphag'] = micro['alphah']
            micro['rhog'] = micro['rhoh']
            consts['MPflg'] = 5
        else:
            # Check microphysics flag in WRF netcdf file.  If it isn't one of those above, exit.
            fatal("This option is not supported for now! Exiting!")

        if 'MPflg' not in consts:
            consts['MPflg'] = 0

        # Also initialize water fraction on snow, graupel, and hail
        micro['qsw'] = N.asfortranarray(N.zeros_like(micro['qs']))
        micro['qgw'] = N.asfortranarray(N.zeros_like(micro['qg']))
        micro['qhw'] = N.asfortranarray(N.zeros_like(micro['qh']))

        return micro, consts

    def loadMeanWind(self):
        self._abstract('loadMeanWind')
        return

    def loadHorizWind(self, average=True):
        # Compute mean storm-relative wind in 0-3km and 0-6km layer
        u = self._file.variables['U'][self._timeindex, :]
        v = self._file.variables['V'][self._timeindex, :]

        # Average u,v to the scalar points

        if average:
            us = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
            vs = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
        else:
            us = u
            vs = v

        return us.T, vs.T

    def loadVertWind(self, average=True):
        w = self._file.variables['W'][self._timeindex, :]
        if average:
            ws = 0.5 * (w[:-1, :, :] + w[1:, :, :])
        else:
            ws = w
        return ws.T

    def loadModelReflectivity(self):
        return

    def getTimeObjs(self):
        return self._timeobjlist

    def fileNameBuilder(self, file_base, time, member):
        return file_base

    def _readwrfgrid(self):
        """This function reads grid information from a WRF netCDF file"""

        xlon = self._file.variables['XLONG'][0, ...]
        xlon_u = self._file.variables['XLONG_U'][0, ...]
        xlon_v = self._file.variables['XLONG_V'][0, ...]
        xlat = self._file.variables['XLAT'][0, ...]
        xlat_u = self._file.variables['XLAT_U'][0, ...]
        xlat_v = self._file.variables['XLAT_V'][0, ...]
        phb = self._file.variables['PHB'][0, ...]
        hgt = self._file.variables['HGT'][0, ...]
        zc = 0.5 * (phb[:-1, :, :] + phb[1:, :, :]) / g - hgt[:, :]
        ze = phb / g - hgt

        return xlon, xlon_u, xlon_v, xlat, xlat_u, xlat_v, zc, ze

    def _readwrfmap(self):
        """This function reads map projection information from a WRF netCDF file"""
        map_proj = self._file.MAP_PROJ
        ctrlat = self._file.CEN_LAT
        ctrlon = self._file.CEN_LON
        trulat1 = self._file.TRUELAT1
        trulat2 = self._file.TRUELAT2
        trulon = self._file.STAND_LON
        return map_proj, ctrlat, ctrlon, trulat1, trulat2, trulon

    def _readwrfmicroLIN(self):
        """This function reads Lin microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # TODO: find out if we need this
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        # scale_rho = True

        micro, consts = {}, {}
        # Many parameters are hardcoded in the LIN WRF code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        # Constant intercept parameter for rain -- This isn't actually used
        consts['n0rain'] = 8.0e6
        # in the WDM6 scheme since rain is double-moment, but is set here
        # for safety/completeness
        # Constant intercept parameter for graupel -- used below to compute Ntg
        consts['n0hail'] = 4.0e4
        # Constant intercept parameter for graupel -- used below to compute Ntg
        consts['n0snow'] = 3.0e6

        consts['rhorcst'] = 1000.    # Bulk raindrop density (kt/m^3)
        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhohcst'] = 913.     # Bulk hail density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        micro['rhod'] = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        consts['graupel_on'] = 0  # No separate graupel category in LIN

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qh'] = self._file.variables['QHAIL'][self._timeindex, :, :, :]
            consts['hail_on'] = 1
        except BaseException:
            micro['qh'] = N.zeros_like(micro['qr'])
            consts['hail_on'] = 0

        micro['rhor'] = consts['rhorcst'] * N.ones_like(micro['qr'])
        micro['rhos'] = consts['rhoscst'] * N.ones_like(micro['qs'])
        micro['rhoh'] = consts['rhohcst'] * N.ones_like(micro['qh'])

        micro['alphar'] = N.ones_like(micro['qr'])
        micro['alphas'] = N.zeros_like(micro['qs'])
        micro['alphah'] = N.zeros_like(micro['qh'])

        # Compute Nts
        cs = (N.pi / 6.) * micro['rhos']
        micro['nts'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qs'],
            consts['n0snow'],
            cs,
            micro['alphas'])

        # Compute Nth
        cg = (N.pi / 6.) * micro['rhoh']
        micro['nth'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qh'],
            consts['n0hail'],
            cg,
            micro['alphah'])

        # Compute Ntr
        cr = (N.pi / 6.) * micro['rhor']
        micro['ntr'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qr'],
            consts['n0rain'],
            cr,
            micro['alphar'])

        for var, data in micro.items():
            if isinstance(data, N.ndarray):
                micro[var] = N.asfortranarray(data.T)

        return micro, consts

    def _readwrfmicroWSM6(self):
        """This function reads WSM6 microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # # Set to true if the number concentration and reflectivity moments need
        # # to be scaled by air density upon read-in
        # scale_rho = True

        micro, consts = {}, {}
        # Many parameters are hardcoded in the WSM6 WRF code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        # Constant intercept parameter for rain -- This isn't actually used
        consts['n0rain'] = 8.0e6
        # in the WSM6 scheme since rain is double-moment, but is set here
        # for safety/completeness
        # Constant intercept parameter for graupel -- used below to compute Ntg
        consts['n0grpl'] = 4.0e6

        consts['rhorcst'] = 1000.    # Bulk raindrop density (kt/m^3)
        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 500.     # Bulk graupel density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        micro['rhod'] = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        consts['hail_on'] = 0  # No separate hail category in WSM6

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0

        micro['rhor'] = consts['rhorcst'] * N.ones_like(micro['qr'])
        micro['rhos'] = consts['rhoscst'] * N.ones_like(micro['qs'])
        micro['rhog'] = consts['rhogcst'] * N.ones_like(micro['qg'])

        micro['alphar'] = N.ones_like(micro['qr'])
        micro['alphas'] = N.zeros_like(micro['qs'])
        micro['alphag'] = N.zeros_like(micro['qg'])

        # Snow in WSM6 uses a temperature-dependent intercept parameter
        consts['n0snow'] = 2.e6 * N.exp(0.12 * (273.15 - micro['tair']))
        # Compute Nts
        cs = (N.pi / 6.) * micro['rhos']
        micro['nts'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qs'],
            consts['n0snow'],
            cs,
            micro['alphas'])

        # Compute Ntg
        cg = (N.pi / 6.) * micro['rhog']
        micro['ntg'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qg'],
            consts['n0grpl'],
            cg,
            micro['alphag'])

        # Compute Ntr
        cr = (N.pi / 6.) * micro['rhor']
        micro['ntr'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qr'],
            consts['n0rain'],
            cr,
            micro['alphar'])

        for var, data in micro.items():
            if isinstance(data, N.ndarray):
                micro[var] = N.asfortranarray(data.T)

        return micro, consts

    def _readwrfmicroWDM6(self):
        """This function reads WDM6 microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        scale_rho = True

        micro, consts = {}, {}
        # Many parameters are hardcoded in the WDM6 WRF code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        # Constant intercept parameter for rain -- This isn't actually used
        consts['n0rain'] = 8.0e6
        # in the WDM6 scheme since rain is double-moment, but is set here
        # for safety/completeness
        # Constant intercept parameter for graupel -- used below to compute Ntg
        consts['n0grpl'] = 4.0e6

        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 500.     # Bulk graupel density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        rhod = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        consts['hail_on'] = 0  # No separate hail category in WDM6

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0

        # Read in rain number concentration field
        try:
            micro['ntr'] = self._file.variables['QNRAIN'][self._timeindex, :, :, :]
            if(scale_rho):
                micro['ntr'] = micro['ntr'] * rhod  # Convert from #/kg to #/m^3
        except BaseException:
            micro['ntr'] = N.zeros_like(micro['qr'])

        micro['rhos'] = consts['rhoscst'] * N.ones_like(micro['qs'])
        micro['rhog'] = consts['rhogcst'] * N.ones_like(micro['qg'])

        micro['alphar'] = N.ones_like(micro['qr'])
        micro['alphas'] = N.zeros_like(micro['qs'])
        micro['alphag'] = N.zeros_like(micro['qg'])

        # Snow in WDM6 uses a temperature-dependent intercept parameter
        consts['n0snow'] = 2.e6 * N.exp(0.12 * (273.15 - micro['tair']))
        # Compute Nts
        cs = (N.pi / 6.) * micro['rhos']
        micro['nts'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qs'],
            consts['n0snow'],
            cs,
            micro['alphas'])

        # Compute Ntg
        cg = (N.pi / 6.) * micro['rhog']
        micro['ntg'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qg'],
            consts['n0grpl'],
            cg,
            micro['alphag'])

        # swap the axes of the arrays
#       rhoa = reshape_array(rhoa,fortran)
#       qr = reshape_array(qr,fortran)
#       qs = reshape_array(qs,fortran)
#       qg = reshape_array(qg,fortran)
#       Ntr = reshape_array(Ntr,fortran)
#       Nts = reshape_array(Nts,fortran)
#       Ntg = reshape_array(Ntg,fortran)
#       alphar = reshape_array(alphar,fortran)
#       alphas = reshape_array(alphas,fortran)
#       alphag = reshape_array(alphag,fortran)
#       rhos = reshape_array(rhos,fortran)
#       rhog = reshape_array(rhog,fortran)
#       tair = reshape_array(tair,fortran)

        for var, data in micro.items():
            if isinstance(data, N.ndarray):
                micro[var] = N.asfortranarray(data.T)

#       return (graupel_on,hail_on,rhoa,qr,qs,qg,Ntr,Nts,Ntg,alphar,alphas,alphag,n0rain,n0snow,
#               n0grpl,rhos,rhog,tair)
        return micro, consts

    def _readwrfmicroMY(self):
        """This function reads MY2 microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        scale_rho = True
        # Many parameters are hardcoded in the WRF microphysics code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.
        micro, consts = {}, {}

        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 400.     # Bulk graupel density (kt/m^3)
        consts['rhohcst'] = 900.     # Bulk hail density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        rhod = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0
        try:
            micro['qh'] = self._file.variables['QHAIL'][self._timeindex, :, :, :]
            consts['hail_on'] = 1
        except BaseException:
            micro['qh'] = N.zeros_like(micro['qr'])
            consts['hail_on'] = 0

        for spec, ncname in [('r', 'RAIN'), ('s', 'SNOW'), ('g', 'GRAUPEL'), ('h', 'HAIL')]:
            # Read in number concentration fields
            try:
                micro['nt' + spec] = self._file.variables['QN' + ncname][self._timeindex, :, :, :]
                if(scale_rho):
                    micro['nt' + spec] *= rhod  # Convert from #/kg to #/m^3
            except BaseException:
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])

            if spec != 'r':
                micro['rho' + spec] = consts["rho%scst" % spec] * N.ones_like(micro['q' + spec])

            micro['alpha' + spec] = N.zeros_like(micro['q' + spec])

        # swap the axes of the arrays
        for var, dat in micro.items():
            micro[var] = N.asfortranarray(dat.T)

        # YJ: Snow in MY uses a temperature-dependent intercept parameter
        # I am using a constant here for simplicity as it is not being used.
        consts['n0rain'] = 1.e7
        consts['n0snow'] = 8.e7
        consts['n0grpl'] = 4.e6
        consts['n0hail'] = 1.e5

        return micro, consts

    def _readwrfmicroMorrison(self):
        """This function reads Morrison microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        scale_rho = True

        # Many parameters are hardcoded in the WRF microphysics code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        micro, consts = {}, {}

        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 400.     # Bulk graupel density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        rhod = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        consts['hail_on'] = 0  # No separate hail category in Morrison

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0

        for spec, ncname in [('r', 'RAIN'), ('s', 'SNOW'), ('g', 'GRAUPEL')]:
            # Read in number concentration fields
            try:
                micro['nt' + spec] = self._file.variables['QN' + ncname][self._timeindex, :, :, :]
                if(scale_rho):
                    micro['nt' + spec] *= rhod  # Convert from #/kg to #/m^3
            except BaseException:
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])

            if spec != 'r':
                micro['rho' + spec] = consts["rho%scst" % spec] * N.ones_like(micro['q' + spec])

            micro['alpha' + spec] = N.zeros_like(micro['q' + spec])

        # swap the axes of the arrays
        for var, dat in micro.items():
            micro[var] = N.asfortranarray(dat.T)

        # YJ: Intercept parameters are never defined in Morrison.
        #     These valuse are defined for computational reason but never be used.
        consts['n0rain'] = 1.e7
        consts['n0snow'] = 8.e7
        consts['n0grpl'] = 4.e6

        return micro, consts

    def _readwrfmicroThomp(self):
        """This function reads Thompson microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        scale_rho = True

        # Many parameters are hardcoded in the WRF code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        micro, consts = {}, {}

        # Constant intercept parameter for rain -- This isn't actually used
        consts['n0rain'] = 8.0e6
        # in the Thompson scheme since rain is double-moment, but is set here
        # for safety/completeness

        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 500.     # Bulk graupel density (kg/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        rhod = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        consts['hail_on'] = 0  # No separate hail category in WDM6

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0

        # Note that rhos is not constant in Thompson scheme but
        micro['rhos'] = consts['rhoscst'] * N.ones_like(micro['qs'])
        # varies as a function of diameter.  This is handled automatically
        # within the emulator code via the assumed mass/diameter relation.
        micro['rhog'] = consts['rhogcst'] * N.ones_like(micro['qg'])

        # N.ones_like(qr), looks like Thompson has 0 for alphar by default
        micro['alphar'] = N.zeros_like(micro['qr'])
        micro['alphas'] = N.zeros_like(micro['qs'])
        micro['alphag'] = N.zeros_like(micro['qg'])

        # Read in rain number concentration field
        try:
            micro['ntr'] = self._file.variables['QNRAIN'][self._timeindex, :, :, :]
            if(scale_rho):
                micro['ntr'] *= rhod  # Convert from #/kg to #/m^3
        except BaseException:
            micro['ntr'] = N.zeros_like(micro['qr'])

        cr = (N.pi / 6.) * 1000.

        # Snow and graupel intercept parameters in Thompson are not constants.
        # Borrowed the code of Bryan Putnam
#       alphas2 = 0.
#       thom_lam0 = 20.78
#       thom_lam1 = 3.29
#       thom_k0 = 490.6
#       thom_k1 = 17.46

        cs = (N.pi / 6.) * micro['rhos']
        cg = (N.pi / 6.) * micro['rhog']
#       moma = dsd.power_mom(2,cs,tair,qs)
#       momb = dsd.power_mom(3,cs,tair,qs)

#       N0s = ((dble(moma)**4)/(dble(momb)**3))*dble(thom_k0)
#       N0s2 = ((dble(moma)**4)/(dble(momb)**3))*dble(thom_k1)*      \
#                      ((dble(moma)/dble(momb))**dble(alphas2))
#
        # Borrow n0snow for WDM6 for now!!!!
        # DTD Note 09/19/2014: The correct distribution is computed within the emulator code,
        # but we should also compute the correct number concentration here as well for completeness
        # at some point.
        consts['n0snow'] = 2.e6 * N.exp(0.12 * (273.15 - micro['tair']))

        # Compute Nts
        micro['nts'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qs'],
            consts['n0snow'],
            cs,
            micro['alphas'])

        # Compute Ntg
        # n0grpl = N.maximum(1.e4, N.minimum(200./qg, 5.e6)) # Old formula
        # The following calculation of N0_g is transcribed from the Thompson scheme code in WRF
        gonv_min = 1.e4
        gonv_max = 3.e6
        N0_min = N.ones_like(micro['qr']) * gonv_max
        qrmask = N.where(micro['qr'] > 1.e-12, True, False)
        D0r = 50.e-6
        mvd_r = dsd.cal_D0(micro['rhoa'], micro['qr'], micro['ntr'], cr, micro['alphar'])
        mvd_r = N.where(mvd_r > 2.5e-3, 2.5e-3, mvd_r)
        mvd_r = N.where(mvd_r < 0.75 * D0r, 0.75 * D0r, mvd_r)
        xslw1 = N.where(
            (micro['tair'] < 270.65) & qrmask & (
                mvd_r > 100.e-6),
            4.01 + N.log10(mvd_r),
            0.01)
        rg = micro['qg'] * micro['rhoa']
        rg = N.where(rg > 5.e-5, rg, 5.e-5)
        ygra1 = 4.31 + N.log10(rg)
        zans1 = 3.1 + (100. / (300. * xslw1 * ygra1 /
                               (10. / xslw1 + 1. + 0.25 * ygra1) + 30. + 10. * ygra1))
        N0_exp = 10.**(zans1)
        N0_exp = N.where(N0_exp < gonv_max, N0_exp, gonv_max)
        N0_exp = N.where(N0_exp > gonv_min, N0_exp, gonv_min)
        N0_min = N.where(N0_exp < N0_min, N0_exp, N0_min)
        am_g = cg
        bm_g = 3.
        oge1 = 1. / (bm_g + 1.)
        cgg_1 = N.exp(gammaln(bm_g + 1.))
        cgg_2 = N.exp(gammaln(micro['alphag'] + 1.))
        cgg_3 = N.exp(gammaln(bm_g + micro['alphag'] + 1.))
        ogg1 = 1. / cgg_1
        ogg2 = 1. / cgg_2
        obmg = 1. / bm_g
        cge_2 = micro['alphag'] + 1.
        lam_exp = (N0_exp * am_g * cgg_1 / rg)**oge1
        lamg = lam_exp * (cgg_3 * ogg2 * ogg1)**obmg
        # ilamg = 1. / lamg
        consts['n0grpl'] = N0_exp / (cgg_2 * lam_exp) * lamg**cge_2
        micro['ntg'] = dsd.cal_Nt(
            micro['rhoa'],
            micro['qg'],
            consts['n0grpl'],
            cg,
            micro['alphag'])

        # swap the axes of the arrays
        for var, dat in micro.items():
            micro[var] = N.asfortranarray(dat.T)

        return micro, consts

    def _readwrfmicroNSSL(self):
        """This function reads MY2 microphysics scheme information from a WRF netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""
        # fortran = True
        # Set to true if the number concentration and reflectivity moments need
        # to be scaled by air density upon read-in
        scale_rho = True

        # Many parameters are hardcoded in the WRF microphysics code and not saved to the history
        # file.  Set them here.  If the user changes any of these parameters in the code,
        # they should also be updated here otherwise bad things will happen.

        micro, consts = {}, {}

        consts['rhoscst'] = 100.     # Bulk snow density (kg/m^3)
        consts['rhogcst'] = 400.     # Bulk graupel density (kt/m^3)
        consts['rhohcst'] = 900.     # Bulk hail density (kt/m^3)

        # Calculate air density

        pt = self._file.variables['T'][self._timeindex, :, :, :] + 300.  # Weird, I know...
        qv = self._file.variables['QVAPOR'][self._timeindex, :, :, :]
        pprt = self._file.variables['P'][self._timeindex, :, :, :]
        pbar = self._file.variables['PB'][self._timeindex, :, :, :]
        p = pbar + pprt

        micro['rhoa'] = thermo.calrho(p, pt, qv)   # Moist air density
        rhod = thermo.calrhod(p, pt)     # Dry air density
        micro['tair'] = thermo.calT(p, pt)        # Air temperature (K)

        # Read in mixing ratios
        micro['qr'] = self._file.variables['QRAIN'][self._timeindex, :, :, :]
        try:
            micro['qs'] = self._file.variables['QSNOW'][self._timeindex, :, :, :]
        except BaseException:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QGRAUP'][self._timeindex, :, :, :]
            consts['graupel_on'] = 1
        except BaseException:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0
        try:
            micro['qvolg'] = self._file.variables['QVGRAUPEL'][self._timeindex, :, :, :]
            micro['rhog'] = micro['qg'] / micro['qvolg']
        except BaseException:
            micro['rhog'] = consts['rhogcst'] * N.ones_like(micro['qg'])
        try:
            micro['qh'] = self._file.variables['QHAIL'][self._timeindex, :, :, :]
            consts['hail_on'] = 1
        except BaseException:
            micro['qh'] = N.zeros_like(micro['qr'])
            consts['hail_on'] = 0
        try:
            micro['qvolh'] = self._file.variables['QVHAIL'][self._timeindex, :, :, :]
            micro['rhoh'] = micro['qh'] / micro['qvolh']
        except BaseException:
            micro['rhoh'] = consts['rhohcst'] * N.ones_like(micro['qh'])

        for spec, ncname in [('r', 'RAIN'), ('s', 'SNOW'), ('g', 'GRAUPEL'), ('h', 'HAIL')]:
            # Read in number concentration fields
            try:
                micro['nt' + spec] = self._file.variables['QN' + ncname][self._timeindex, :, :, :]
                if(scale_rho):
                    micro['nt' + spec] *= rhod  # Convert from #/kg to #/m^3
            except BaseException:
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])

            if spec == 's':
                micro['rho' + spec] = consts["rho%scst" % spec] * N.ones_like(micro['q' + spec])

        micro['alphas'] = -0.4 * N.ones_like(micro['qs'])
        micro['alphag'] = N.zeros_like(micro['qg'])
        micro['alphah'] = N.ones_like(micro['qh'])

        # swap the axes of the arrays
        for var, dat in micro.items():
            micro[var] = N.asfortranarray(dat.T)

        # YJ: Snow in MY uses a temperature-dependent intercept parameter
        # I am using a constant here for simplicity as it is not being used.
        consts['n0rain'] = 1.e7
        consts['n0snow'] = 8.e7
        consts['n0grpl'] = 4.e6
        consts['n0hail'] = 1.e5

        micro['alphah'] = N.ones_like(micro['qr'])

        return micro, consts


def reshape_array(a, fortran=False):
    """Reshapes a COMMAS array into i,j,k (x,y,z) ordering.
       Also allows for converting to fortran memory order (default no)."""
    a = a.swapaxes(0, 1)
    a = a.swapaxes(1, 2)
    a = a.swapaxes(0, 1)
    if fortran:
        a = N.asfortranarray(a)
    return a


def calp_exner(timeindex, exner):
    """Calculate pressure from exner function"""
    Rd = 287.0
    cp = 1005.6
    p = (exner**(cp / Rd)) * 100000.0
    # p = reshape_array(p)

    return p
