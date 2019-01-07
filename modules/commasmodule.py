# commasmodule.py
# A collection of functions related to reading and manipulating commas netCDF files

import netCDF4 as netcdf
import numpy as N
import thermolib as thermo
import dualpara as dualpol
import sys

from utils import log, warning, fatal
from datahandler import DataHandler

import datetime

class COMMASDataHandler(DataHandler):
    def __init__(self, base_dir, times, multitime=True):
        self._base_dir = base_dir
        self._time_list = times
        self._file = None
        self._file_name = ""
        self._cur_time = -1
        self._multitime = multitime
        self.vardict = None

        super(COMMASDataHandler, self).__init__('COMMAS')
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
        temptimeindex = N.where(self._filetimes == time)[0][0]
        timestring = self._datetimelist[temptimeindex].strftime("%d %B %Y %H%M UTC")
        if(self._multitime):
            self._timeindex = temptimeindex
        else:
            self._timeindex = 0
        return timestring

    def loadGrid(self):
        # Read grid variables
        xc, xe, yc, ye, zc1d, ze1d = self._readcommasgrid()

        nxm = xc.shape[0]
        nxe = xe.shape[0]
        nym = yc.shape[0]
        nye = ye.shape[0]
        nzm = zc1d.shape[0]
        nze = ze1d.shape[0]

        # Recast everything as 3D arrays
        xc = N.lib.stride_tricks.as_strided(xc, strides=((0, 0) + xc.strides), shape=((nzm, nym) + xc.shape))
        xe = N.lib.stride_tricks.as_strided(xe, strides=((0, 0) + xe.strides), shape=((nzm, nym) + xe.shape))

        yc = N.lib.stride_tricks.as_strided(yc, strides=((0, ) + yc.strides + (0, )), shape=((nzm, ) + yc.shape + (nxm, )))
        ye = N.lib.stride_tricks.as_strided(ye, strides=((0, ) + ye.strides + (0, )), shape=((nzm, ) + ye.shape + (nxm, )))

        zc = N.lib.stride_tricks.as_strided(zc1d, strides=(zc1d.strides + (0, 0)), shape=(zc1d.shape + (nym, nxm)))
        ze = N.lib.stride_tricks.as_strided(ze1d, strides=(ze1d.strides + (0, 0)), shape=(ze1d.shape + (nym, nxm)))

        return xc, yc, zc, zc1d, xe, ye, ze, ze1d, None


    def loadVars(self, varnamelist, boxindices=None, squeeze=True):
        vardict = {}
        for varname in varnamelist:
            try:
                if boxindices is not None:
                    var = self._file.variables[varname][self._timeindex,
                                                        boxindices[4]:boxindices[5],
                                                        boxindices[2]:boxindices[3],
                                                        boxindices[0]:boxindices[1]]
                    if squeeze:
                        var = var.squeeze()
                else:
                    var = self._file.variables[varname][self._timeindex,...]
                vardict[varname] = var
            except:
                warning('Variable '+varname+' not found in file!')
        return vardict


    def loadTimes(self):
        if(self._multitime): # All times in one file
            self._filetimes = self._file.variables['TIME'][:]
        else:
            self._filetimes = N.array(self._time_list)

        # Read in time string for start of times in file
        timestringstart = self._file.COARDS
        datetimestart = datetime.datetime.strptime(timestringstart, "seconds since %Y-%m-%d %H:%M:%S")
        # Now construct datetime objects for each time in desired timelist
        time_dt = datetime.timedelta(seconds=int(self._filetimes[0]))
        self._datetimelist = [ datetimestart+time_dt ]

        for i in xrange(1, len(self._filetimes)):
            time_dt = datetime.timedelta(seconds=int(self._filetimes[i] - self._filetimes[i-1]))
            self._datetimelist.append(self._datetimelist[i-1]+time_dt)
        return


    def loadMicrophysics(self):
        microphys = self._file.MICROPHYS

        if microphys == "LFO":
            log("Reading microphysics information for "+microphys+" scheme.")
            micro, consts = self._readcommasmicroLFO()

            micro['qh'] = N.zeros_like(micro['qr'])
            micro['rhoh'] = N.zeros_like(micro['rhor'])
            consts['n0hail'] = 0.0

            for spec in ['r', 's', 'g', 'h']:

                # shape parameters are all assumed 0 for LFO (exponential)
                micro['alpha' + spec] = N.zeros_like(micro['q' + spec])

                # Set to zero unused microphysics variables
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])

                if spec != 'r':
                    consts['rho%scst' % spec] = N.zeros_like(micro['rho' + spec])
                    micro['q%sw' % spec] = N.zeros_like(micro['q' + spec])

            consts['ZVDflg'] = 0
            consts['imurain'] = -1
            consts['MPflg'] = 0

        elif microphys[:2] == "MY":
            log("Reading microphysics information for "+microphys+" scheme.")
            micro, consts = self._readcommasmicroMY()

            for spec in ['s', 'g', 'h']:
                micro['q%sw' % spec] = N.zeros_like(micro['q' + spec])
                consts['rho%scst' % spec] = micro['rho' + spec]

            consts['imurain'] = -1
            consts['ZVDflg'] = 0
            consts['MPflg'] = 0

        elif microphys[:3] == "ZVD" or microphys[:3] == "WAR" or microphys[:3] == "ZIE":
            log("Reading microphysics information for "+microphys+" scheme.")
            micro, consts = self._readcommasmicroZVD()

            if(consts['imurain'] == 3):
                consts['MPflg'] = 1
            else:
                consts['MPflg'] = 2 # New ZVD option for gamma-diameter rain

            if(microphys == "ZVDM" or microphys == "ZVDHM" or microphys == "ZIEGM" or microphys == "ZIEGHM"):
                consts['MFflg'] = 2
            if(microphys == "ZVDM" or microphys == "ZVD" or microphys == "ZIEG"): # No hail in this run, set graupel to off
                warning("No hail in this run, setting hail flag to off.")
                consts['hail_on'] = 0

        elif microphys[:3] == "TAK":
            log("Reading microphysics information for "+microphys+" scheme.")
            graupel_on,hail_on,rhoa,qr,qs,qg,qh,nr_bin,ns_bin,ng_bin,nh_bin,  \
                rhorcst,rhoscst,rhogcst,rhohcst = self._readcommasmicroTAK()
            imurain = -1
            ZVDflg = 0

            # Compute total number concentration variables from bins
            Ntr = N.zeros_like(qr)
            Nts = N.zeros_like(qs)
            Ntg = N.zeros_like(qg)
            Nth = N.zeros_like(qh)

            # Rain total number concentration
            for b in xrange(23,34):
                Ntr = Ntr + nr_bin[b,...]*tak.rbinwidth[b]

            # Snow
            for b in xrange(21):
                for k in range(5):
                    Nts = Nts + ns_bin[k,b,...]*tak.sbinwidth[k,b]

            # Graupel and Hail
            for b in xrange(45):
                Ntg = Ntg + ng_bin[b,...]*tak.gbinwidth[b]
                Nth = Nth + nh_bin[b,...]*tak.hbinwidth[b]

            alphar = N.zeros_like(qr)
            alphas = N.zeros_like(qr)
            alphag = N.zeros_like(qr)
            alphah = N.zeros_like(qr)
            qsw = N.zeros_like(qr)
            qgw = N.zeros_like(qr)
            qhw = N.zeros_like(qr)
            rhos = rhoscst
            rhog = rhogcst
            rhoh = rhohcst

            n0rain = 0.0
            n0snow = 0.0
            n0grpl = 0.0
            n0hail = 0.0
        else:
            fatal("Microphysics scheme '%s' is not supported for COMMAS!" % microphys)

        consts['microphys'] = microphys

        return micro, consts

    def loadVertWind(self):
        w = self._file.variables['W'][self._timeindex,:,:,:]
        wsovr = 0.5*(w[:-1,:,:]+w[1:,:,:])
        return wsovr.T

    def loadHorizWind(self):
        u = self._file.variables['U'][self._timeindex,:] - self._file.variables['UGRID'][0]
        v = self._file.variables['V'][self._timeindex,:] - self._file.variables['VGRID'][0]

        us = 0.5*(u[:,:,:-1]+u[:,:,1:])
        vs = 0.5*(v[:,:-1,:]+v[:,1:,:])

        return us.T, vs.T

    def loadModelReflectivity(self):
        Zmod = self._file.variables['DBZ'][self._timeindex,:,:,:]
        Zmod = reshape_array(Zmod,False)
        return Zmod

    def fileNameBuilder(self, file_base, time, member):
        if self._multitime:
            filename = "%s.%03d.nc" % (file_base, member)
        else:
            # Not sure why COMMAS removes the member number when the times are in separate files?
            filename = "%s.%06d.nc" % (file_base, time)
        return filename

    def getTimeObjs(self):
        return self._datetimelist

    def getGridPosition(self):
        xg_pos = self._file.variables['XG_POS'][self._timeindex]
        yg_pos = self._file.variables['YG_POS'][self._timeindex]

        return xg_pos, yg_pos

    def _readcommasgrid(self):
        """This function reads grid information from a COMMAS netCDF file"""

        xc = self._file.variables['XC'][...]
        yc = self._file.variables['YC'][...]
        xe = self._file.variables['XE'][...]
        ye = self._file.variables['YE'][...]
        zc = self._file.variables['ZC'][...]
        ze = self._file.variables['ZE'][...]

        return xc,xe,yc,ye,zc,ze

    def _readcommasmicroMY(self):
        """Reads in several microphysical variables for the MY scheme and derives others from a COMMAS
           netCDF file.  Determines if graupel and/or hail are present and returns appropriate flags.
           Also checks if the sixth moments are present and computes the shape parameter as appropriate.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""

        fortran = True

        scale_rho = False       # Set to true if the number concentration and reflectivity moments need
                                # to be scaled by air density upon read-in

        # Calculate moist air density

        pt = self._file.variables['TH'][self._timeindex,:,:,:]
        qv = self._file.variables['QV'][self._timeindex,:,:,:]
        zc = self._file.variables['ZC'][...]

        # Retrieve exner function

        exner = N.zeros_like(pt)

        for k in range(0,zc.size):
            exner[k] = self._file.variables['PI'][self._timeindex,k,:,:]+self._file.variables['PIINIT'][self._timeindex, k]

        p = calp_exner(self._timeindex,exner)

        micro, consts = {}, {}

        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density
        # Compute air temperature
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)

        # Read in mixing ratios

        micro['qr'] = self._file.variables['QR'][self._timeindex,:,:,:]
        micro['qs'] = self._file.variables['QS'][self._timeindex,:,:,:]
        try:
            micro['qg'] = self._file.variables['QH'][self._timeindex,:,:,:]
            consts['graupel_on'] = 1
        except:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0
        try:
            micro['qh'] = self._file.variables['QHL'][self._timeindex,:,:,:]
            consts['hail_on'] = 1
        except:
            micro['qh'] = N.zeros_like(micro['qr'])
            consts['hail_on'] = 0

        for spec, name, ncname in [('r', 'rain', 'RW'), ('s', 'snow', 'SW'), ('g', 'grpl', 'HW'), ('h', 'hail', 'HL')]:
            # Try to read in number concentration variables
            try:
                micro['nt' + spec] = self._file.variables['C' + ncname][self._timeindex,:,:,:]
            except:
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])
            consts['n0' + name] = self._file.variables['CNO%s_MY' % spec.upper()][:]

            # Read in particle densities
            micro['rho' + spec] = self._file.variables['RHO_Q%s_MY' % spec.upper()][:] * N.ones_like(micro['q' + spec])

            # Now try to read 6th moment (reflectivity) variables
            # and compute the corresponding shape parameters
            # If they don't exist, use the constant alpha values
            try:
                micro['z' + spec] = self._file.variables['Z' + ncname][self._timeindex,:,:,:]
                if(scale_rho):
                    micro['z' + spec] *= rhod

                mu = 1./3.
                micro['alpha' + spec] = dualpol.solve_alpha_iter(micro['rhoa'], mu, micro['q' + spec], micro['nt' + spec], micro['z' + spec], micro['rho' + spec])
                print "Found Z%s array! Computed shape parameter for %s." % (spec.upper(), name)
            except:
                print "Did not find Z%s array! Assuming 2-moment with constant shape parameter." % spec.upper()
                zr = N.zeros_like(micro['q' + spec])
                try:
                    micro['alpha' + spec] = self._file.variables['ALPHA%s_MY' % spec.upper()][:] * N.ones_like(micro['q' + spec])
                except:
                    micro['alpha' + spec] = N.zeros_like(micro['q' + spec])

            # Scale N and Z by dry air density if needed
            if(scale_rho):
                micro['nt' + spec] *= rhod

        # swap the axes of the arrays
        for var, dat in micro.iteritems():
            micro[var] = N.asfortranarray(dat.T)

        return micro, consts


    def _readcommasmicroLFO(self):
        """Reads in several microphysical variables for the LFO scheme and derives others from a COMMAS
           netCDF file.
           Note, for now qi and qc are ignored since they aren't used in the polarimetric
           emulator, but obviously need to add these for use by other applications."""

        fortran = True

        # Calculate moist air density

        pt = self._file.variables['TH'][self._timeindex,:,:,:]
        qv = self._file.variables['QV'][self._timeindex,:,:,:]
        zc = self._file.variables['ZC'][...]

        # Retrieve exner function

        exner = N.zeros_like(pt)

        micro, consts = {}, {}

        # Compute air temperature
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)

        for k in range(0,zc.size):
            exner[k] = self._file.variables['PI'][self._timeindex,k,:,:]+self._file.variables['PIINIT'][self._timeindex, k]

        p = calp_exner(self._timeindex,exner)

        micro['rhoa'] = thermo.calrho(p,pt,qv)   # Moist air density
        rhod = thermo.calrhod(p,pt)     # Dry air density

        consts['graupel_on'] = 1
        consts['hail_on'] = 0 # No separate hail category in LFO

        for spec, name, ncname in [ ('r', 'rain', 'R'), ('s', 'snow', 'S'), ('g', 'grpl', 'H') ]:
            # Read in mixing ratios
            micro['q' + spec] = self._file.variables['Q' + ncname][self._timeindex,:,:,:]

            # read in intercept parameters
            consts['n0' + name] = self._file.variables['CNO' + ncname][:]

            # Read in particle densities
            micro['rho' + spec] = self._file.variables['RHO_Q' + ncname][:]*N.ones_like(micro['q' + spec])

        # swap the axes of the arrays
        for var, dat in micro.iteritems():
            micro[var] = N.asfortranarray(dat.T)

        return micro, consts

    def _readcommasmicroZVD(self):
        """Reads in several microphysical variables for the ZVD scheme and derives others from a COMMAS netCDF file
           Determines if graupel and hail are present and returns the appropriate flags.  Also
           checks whether liquid water fraction, variable density, and sixth moments are present
           and computes the shape parameter as appropriate.Note, for now qi and qc are ignored
           since they aren't used in the polarimetric emulator, but obviously need to add these
           for use by other applications."""

        fortran = True
        # Calculate moist air density, needed for computation of graupel and hail density below

        pt = self._file.variables['TH'][self._timeindex,:,:,:]
        qv = self._file.variables['QV'][self._timeindex,:,:,:]
        zc = self._file.variables['ZC'][...]

        # Retrieve exner function

        exner = N.zeros_like(pt)

        for k in range(0,zc.size):
            exner[k] = self._file.variables['PI'][self._timeindex,k,:,:].astype(N.float64)+self._file.variables['PIINIT'][self._timeindex, k]

        p = calp_exner(self._timeindex,exner)

        micro, consts = {}, {}

        micro['rhoa'] = thermo.calrho(p,pt,qv)

        # Compute air temperature
        micro['tair'] = thermo.calT(p,pt)        # Air temperature (K)

        # Read in mixing ratios

        micro['qr'] = self._file.variables['QR'][self._timeindex,:,:,:]
        try:
            micro['qs'] = self._file.variables['QS'][self._timeindex,:,:,:]
        except:
            micro['qs'] = N.zeros_like(micro['qr'])
        try:
            micro['qg'] = self._file.variables['QH'][self._timeindex,:,:,:]
            consts['graupel_on'] = 1
        except:
            micro['qg'] = N.zeros_like(micro['qr'])
            consts['graupel_on'] = 0
        try:
            micro['qh'] = self._file.variables['QHL'][self._timeindex,:,:,:]
            consts['hail_on'] = 1
        except:
            micro['qh'] = N.zeros_like(micro['qr'])
            consts['hail_on'] = 0

        try:
            volg = self._file.variables['VHW'][self._timeindex,:,:,:]
            micro['rhog'] = micro['qg'] * micro['rhoa'] / volg
        except:
            micro['rhog'] = self._file.variables['RHO_QH'][:]*N.ones_like(micro['qg'])
        try:
            volh = self._file.variables['VHL'][self._timeindex,:,:,:]
            micro['rhoh'] = micro['qh'] * micro['rhoa'] / volh
        except:
            micro['rhoh'] = self._file.variables['RHO_QHL'][:] * N.ones_like(micro['qh'])

        #With new gamma-diameter option, we need to find out which one we are using
        try:
            consts['imurain'] = self._file.variables['IMURAIN'].getValue()
        except:
            consts['imurain'] = 3 # Default is gamma-volume

        for spec, name, ncname in [('r', 'rain', 'RW'), ('s', 'snow', 'SW'), ('g', 'grpl', 'HW'), ('h', 'hail', 'HL') ]:
            # Try to read in number concentration variables
            try:
                micro['nt' + spec] = self._file.variables['C' + ncname][self._timeindex,:,:,:]
            except:
                micro['nt' + spec] = N.zeros_like(micro['q' + spec])

            try:
                consts['n0' + name] = self._file.variables['CNO' + ncname.strip('W')][:].astype(N.float64)
            except:
                consts['n0' + name] = 0.0

            # Read in particle densities
            try:
                micro["rho%scst" % spec] = self._file.variables['RHO_Q' + spec.upper()][:] * N.ones_like(micro['q' + spec])
            except:
                micro["rho%scst" % spec] = N.zeros_like(micro['q' + spec])

            if spec in ['r', 's']:
                micro['rho' + spec] = micro["rho%scst" % spec]


            if spec != 'r':
                # Try to read melted water fraction variables
                # If they don't exist, just set them to zero
                try:
                    micro["q%sw" % spec] = self._file.variables["Q%sW" % ncname.strip('W')][self._timeindex,:,:,:]
                except:
                    micro["q%sw" % spec] = N.zeros_like(micro['q' + spec])

            # Now try to read 6th moment (reflectivity) variables
            # and compute the corresponding shape parameters
            # If they don't exist, use the constant alpha values

            try:
                micro['z' + spec] = self._file.variables['Z' + ncname][self._timeindex,:,:,:]

                if spec == 'r' and consts['imurain'] == 3:
            	    print "Rain is gamma-volume."
                    mu = 1.0
                elif spec == 's':
                    mu = 1.0
                else:
                    if spec == 'r': print "Rain is gamma-diameter."
                    mu = 1./3.

                micro['alpha' + spec] = dualpol.solve_alpha_iter(micro['rhoa'],mu,micro['q' + spec], micro['nt' + spec], micro['z' + spec], micro['rho' + spec])
                print "Found Z%s array! Computed shape parameter for %s." % (spec.upper(), name)
            except:
                print "Did not find Z%s array! Assuming 2-moment with constant shape parameter." % spec.upper()
                micro['z' + spec] = N.zeros_like(micro['q' + spec])
                try:
                    micro['alpha' + spec] = self._file.variables['ALPHA' + ncname.strip('W')][:]*N.ones_like(qr)
                except:
                    if spec in ['r', 's']:
                        micro['alpha' + spec] = -0.4*N.ones_like(micro['q' + spec])
                    else:
                        micro['alpha' + spec] = N.zeros_like(micro['q' + spec])


        for var, dat in micro.iteritems():
            micro[var] = N.asfortranarray(dat.T)

        return micro, consts

    def _readcommasmicroTAK(self):
        """Reads in several microphysical variables for the TAK scheme and derives others from a COMMAS netCDF file"""

        # Set up arrays for bin radius and volume to be used later for TAK scheme
        rr_bin = 2.0e-6*N.exp(N.arange(0.,34.,1.)/4.329)
        rbinwidth = 2.*rr_bin/4.329
        rs_bin = N.zeros((5,21))
        for i in xrange(5):
            rs_bin[i,:] = 2.0e-4*N.exp(N.arange(0.,21.,1.)/2.885)
        #print rs_bin
        hs_bin = 1.0e-4*N.exp(N.tile(N.arange(0.,5.,1.),(21,1)).transpose()/(4./N.log(2.*rs_bin/2.0e-3)))    # Snow/Ice thickness
        #print hs_bin
        sbinwidth = 2.*rs_bin/2.885
        rg_bin = 2.0e-6*N.exp(N.arange(0.,45.,1.)/4.329)
        rh_bin = rg_bin
        gbinwidth = 2.*rg_bin/4.329
        hbinwidth = 2.*rh_bin/4.329
        vr_bin = (4.0*N.pi/3.0)*rr_bin**3.
        ms_bin = (N.pi*(2.0*rs_bin)**2./4.0)*((hs_bin-1.0e-5)+1.0e-5*100.)  # Mass of ice crystal
        vg_bin = (4.0*N.pi/3.0)*rg_bin**3.
        vh_bin = vg_bin


        fortran = True

        if(fortran):
        	order='F'
        else:
        	order='C'

        nxe = self._file.NXEND
        nye = self._file.NYEND
        nze = self._file.NZEND

        # Calculate moist air density, needed for computation of graupel and hail density below

        pt = self._file.variables['TH'][self._timeindex,:,:,:]
        qv = self._file.variables['QV'][self._timeindex,:,:,:]
        zc = self._file.variables['ZC'][...]

        # Retrieve exner function

        exner = N.zeros_like(pt)

        for k in range(0,zc.size):
            exner[k] = self._file.variables['PI'][self._timeindex,k,:,:]+self._file.variables['PIINIT'][self._timeindex, k]

        p = calp_exner(self._timeindex,exner)

        rhoa = thermo.calrho(p,pt,qv)
        # Compute air temperature
        tair = thermo.calT(p,pt)        # Air temperature (K)

        # Read in concentration bins for cloud/rain,ice/snow,graupel, and hail
        # They are stored in separate 4D variables for each bin

        # Cloud/rain

        nr_bin = N.zeros((34,nxe-1,nye-1,nze-1),order=order)

        for i in xrange(34):
        	if(i == 0):
        		varname = 'CRW'
        	else:
        		varname = 'CRW%02d' % (i+1)
        	nrtmp = self._file.variables[varname][self._timeindex,:,:,:]*1.e6/rbinwidth[i] # Convert from cm^-3 to m^-3/m
        	nrtmp = reshape_array(nrtmp,fortran)
        	nr_bin[i,:] = nrtmp


        # Ice/snow

        ns_bin = N.zeros((5,21,nxe-1,nye-1,nze-1),order=order)

        for k in xrange(5):
        	for i in xrange(21):
        		if(k == 0 and i == 0):
        			varname = 'CCI'
        		else:
        			varname = 'CCI'+str(k+1)+'_%02d' % (i+1)
        		nstmp = self._file.variables[varname][self._timeindex,:,:,:]*1.e6/sbinwidth[k,i]
        		nstmp = reshape_array(nstmp,fortran)
        		ns_bin[k,i,:] = nstmp

        # Graupel

        ng_bin = N.zeros((45,nxe-1,nye-1,nze-1),order=order)

        for i in xrange(45):
        	if(i == 0):
        		varname = 'CHW'
        	else:
        		varname = 'CHW%02d' % (i+1)
        	ngtmp = self._file.variables[varname][self._timeindex,:,:,:]*1.e6/gbinwidth[i] # Convert from cm^-3 to m^-3
        	ngtmp = reshape_array(ngtmp,fortran)
        	ng_bin[i,:] = ngtmp

        # Hail

        nh_bin = N.zeros_like(ng_bin)

        for i in xrange(45):
        	if(i == 0):
        		varname = 'CHL'
        	else:
        		varname = 'CHL%02d' % (i+1)
        	nhtmp = self._file.variables[varname][self._timeindex,:,:,:]*1.e6/hbinwidth[i] # Convert from cm^-3 to m^-3
        	nhtmp = reshape_array(nhtmp,fortran)
        	nh_bin[i,:] = nhtmp

        # Read in mixing ratios (summed from Takahashi bins in COMMAS)

        qr = self._file.variables['QR'][self._timeindex,:,:,:]
        try:
            qs = self._file.variables['QS'][self._timeindex,:,:,:]
        except:
            qs = N.zeros_like(qr)
        try:
            qg = self._file.variables['QH'][self._timeindex,:,:,:]
            graupel_on = 1
        except:
            qg = N.zeros_like(qr)
            graupel_on = 0
        try:
            qh = self._file.variables['QHL'][self._timeindex,:,:,:]
            hail_on = 1
        except:
            qh = N.zeros_like(qr)
            hail_on = 0

        # Read rain and snow densities

        rhorcst = self._file.variables['RHO_QR'][:]*N.ones_like(qr)
        try:
            rhosread = self._file.variables['RHO_QS'][:]
        except:
            rhosread = N.zeros_like(rhorcst)
        rhoscst = rhosread*N.ones_like(qs)

        # Read in (constant) densities for graupel and hail

        try:
            rhogcst = self._file.variables['RHO_QH'][:]*N.ones_like(qg)
        except:
            rhogcst = N.zeros_like(rhorcst)
        try:
            rhohcst = self._file.variables['RHO_QHL'][:]*N.ones_like(qh)
        except:
            rhohcst = N.zeros_like(rhohcst)

     	# TODO below: compute 6th moments from Takahashi bins and return them?

        # swap the axes of the arrays
        rhoa = reshape_array(rhoa,fortran)
        qr = reshape_array(qr,fortran)
        qs = reshape_array(qs,fortran)
        qg = reshape_array(qg,fortran)
        qh = reshape_array(qh,fortran)
        rhorcst = reshape_array(rhorcst,fortran)
        rhoscst = reshape_array(rhoscst,fortran)
        rhogcst = reshape_array(rhogcst,fortran)
        rhohcst = reshape_array(rhohcst,fortran)
        tair = reshape_array(tair,fortran)

        # Hack to hard code constant rhog and rhoh for now, since they are not saved in the netCDF
        # file properly for the TAK scheme

        rhogcst[:] = 300.
        rhohcst[:] = 900.

        # Now return all the variables

        return graupel_on,hail_on,rhoa,qr,qs,qg,qh,nr_bin,ns_bin,ng_bin,nh_bin,  \
        rhorcst,rhoscst,rhogcst,rhohcst,tair

def commascalradTmat(microphysics,MFflg,timeindex,filename,wavelen,dirscatt,fortran=False):
    """Calculates dual-pol radar variables for COMMAS model output.
       Uses Youngsun Jung's T-matrix scattering-based code."""

    # Determine whether to read ZVD or MY info (LFO support to be added later)

    if (microphysics[:2] == "MY"):
        MPflg = 0
        MFflg = 0
        print "Reading microphysics information for MY scheme"
        graupel_on,hail_on,rhoa,qr,qs,qg,qh,Ntr,Nts,Ntg,Nth,zr,zs,zg,zh,n0rain,n0snow,n0grpl,n0hail,rhor,rhos,rhog,rhoh = \
        readcommasmicroMY(timeindex,filename,fortran)
    elif (microphysics[:3] == "ZVD"):
        MPflg = 1
        print "Reading microphysics information for ZVD scheme"
        graupel_on,hail_on,rhoa,qr,qs,qg,qh,Ntr,Nts,Ntg,Nth,zr,zs,zg,zh,n0rain,n0snow,n0grpl,n0hail,rhosnow,  \
        rhogrpl,rhohail,qsw,qgw,qhw,rhor,rhos,rhog,rhoh = readcommasmicroZVD(timeindex,filename,fortran)
    else:
        print "Microphysics not supported yet!"
        sys.exit()

    dualpol.dualpara.setgrplhl(graupel_on,hail_on)

    print "Initializing DSD parameters"
    #dualpol.dualpara.init_dsd()
    dualpol.dualpara.model_dsd(n0rain,n0snow,n0grpl,n0hail,rhosnow,rhogrpl,rhohail)

    print "Calling dualpol subroutine"
    logz,sumzh,sumzv,logzdr,sumzhv,kdp,ahh,avv = dualpol.refl_rsa_array(MPflg,MFflg,dirscatt,wavelen,rhoa,qr,qs,qg,qh,
                                                                    Ntr,Nts,Ntg,Nth,alphar,alphas,alphag,alphah,
                                                                    qsw,qgw,qhw,rhog,rhoh)

def dxy_2_ll(x, y, lat1, lon1):
    """dxy_2_ll returns the approximate lat/lon between an x,y coordinate and
       a reference lat/lon point.  Assumes a flat earth approximation (which is
       sufficient for radar data) and should be good out to distances of ~200 km.

       INPUTS:  x,y in meters, lat1, lon1 in radians.  Returns degrees

       if x > 0, lon > lon1

       if y > 0, lat > lat1

       OUTPUTS:  lat, lon in radians
    """
    rearth = 1000.0 * 6367.0

    if lon1 < 0.0:
        lon1p = lon1+2.0*N.pi
    else:
        lon1p = lon1

    lat = lat1 + y / rearth
    lon = lon1 + x / ( rearth * N.cos(0.5*(lat1+lat)) )

    lat = N.rad2deg(lat)
    lon = N.rad2deg(lon)
    lon = N.where( lon < 180., lon, 180 - lon)

    return lat, lon

# Comment 01/18/2012: Need to revisit this fortran memory order issue. It may not be doing what
# I think it is doing.
def reshape_array(a,fortran=False):
    """Reshapes a 3D COMMAS array into i,j,k (x,y,z) ordering.
       Also allows for converting to fortran memory order (default no)."""
    a = a.swapaxes(0,1)
    a = a.swapaxes(1,2)
    a = a.swapaxes(0,1)
    if fortran:
        a = N.asfortranarray(a)
    return a

def calp_exner(timeindex,exner):
    """Calculate pressure from exner function"""
    Rd = 287.0
    cp = 1005.6
    p = (exner**(cp/Rd))*100000.0
    #p = reshape_array(p)

    return p
def my_colormaps(name="jet"):
    """ A few colormaps stolen from one of Lou's old scripts.
    They apparently don't work anymore"""
    if name == "dbz_15_75":
        return[(0.000,0.000,0.804), \
            (0.482,0.988,0.000), \
            (0.000,0.933,0.000), \
            (0.000,0.545,0.000), \
            (1.000,1.000,0.000), \
            (0.726,0.525,0.043), \
            (1.000,0.549,0.000), \
            (1.000,0.000,0.000), \
            (0.804,0.000,0.000), \
            (0.549,0.000,0.000), \
            (0.933,0.070,0.537), \
            (0.604,0.196,0.078) ]
    if name == "dbz_0_75":
        return [ (0.000,1.000,1.000), \
             (0.118,0.566,1.000), \
             (0.000,0.000,0.804), \
             (0.482,0.988,0.000), \
             (0.000,0.933,0.000), \
             (0.000,0.545,0.000), \
             (1.000,1.000,0.000), \
             (0.726,0.525,0.043), \
             (1.000,0.549,0.000), \
             (1.000,0.000,0.000), \
             (0.804,0.000,0.000), \
             (0.549,0.000,0.000), \
             (0.933,0.070,0.537), \
             (0.604,0.196,0.078) ]

    if name == "jet":
        return "jet"

