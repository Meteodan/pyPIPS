# radarmodule.py: A collection of functions to read and plot radar data

import Nio
import numpy as np
import matplotlib
# matplotlib.use('TkAgg')
# matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from matplotlib.projections import PolarAxes, register_projection
from matplotlib.transforms import Affine2D, Bbox, IdentityTransform
from mpl_toolkits.axes_grid1 import ImageGrid
# import ctablesfrompyesviewer as ctables
from metpy.plots import ctables
from datetime import datetime, timedelta
from . import obanmodule as oban
import glob
import os
import pyart
from . import utils
from . import timemodule as tm
from . import disdrometer_module as dis
import pickle

clevels_ref = np.arange(0.0, 85.0, 5.0)          # Contour levels for reflectivity (dBZ)
clevels_zdr = np.arange(0.0, 6.25, 0.25)         # Contour levels for Zdr (dB)
clevels_kdp = np.arange(0.0, 12.5, 0.5)            # Contour levels for Kdp (deg km^-1)
clevels_rhv = np.arange(0.9, 1.01, 0.01)           # Contour levels for rhv
clevels_vr = np.arange(-40.0, 41.0, 1.0)         # Contour levels for Vr (m/s)
clevels_kdp = np.arange(-1.0, 7.0, 0.5)		  # Contour levels for Kdp (deg/km)
clevels_rhv = np.arange(0.0, 1.0, 0.1)			  # Contour levels for Rhv
clevels_rain = np.arange(0.1, 200.0)        # Contour levels for rainrate
clevels_d0 = np.arange(0.0, 4.0, 0.1)             # Contour levels for median drop size
clevels_lwc = np.arange(0.01, 10.0)        # Contour levels for liquid water content
clevels_nt = np.arange(-0.5, 4.5, 0.5)         # Contour levels for log total number concentration
clevels_sigm = np.arange(0.0, 2.0, 0.05)        # Contour levels for sigma
clevels_mu = np.arange(-2.0, 20.0, 1.0)        # Contour levels for shape parameter
clevels_lam = np.arange(0.0, 20.0, 1.0)        # Contour levels for slope parameter
# cmapdBZ = ctables.__getattribute__('REF_default')
normdBZ, cmapdBZ = ctables.registry.get_with_steps('NWSReflectivity', 5., 5.)
cmapzdr = cm.Reds
cmapvr = cm.RdBu_r
cmapkdp = cm.Set1
cmaprhv = cm.Spectral
cmapretrievals = cm.viridis


def mtokm(val, pos):
    """Convert m to km for formatting axes tick labels"""
    val = val / 1000.0
    return '%i' % val

# The following code is taken from the following URL:
# http://stackoverflow.com/questions/2417794/
# how-to-make-the-angles-in-a-matplotlib-polar-plot-go-clockwise-with-0-at-the-to
# It defines a Polar projection with 0 degrees at the top and angles increasing clockwise
class NorthPolarAxes(PolarAxes):
    '''
    A variant of PolarAxes where theta starts pointing north and goes
    clockwise.
    '''
    name = 'northpolar'

    class NorthPolarTransform(PolarAxes.PolarTransform):
        def transform(self, tr):
            xy = np.zeros(tr.shape, np.float_)
            t = tr[:, 0:1]
            r = tr[:, 1:2]
            x = xy[:, 0:1]
            y = xy[:, 1:2]
            x[:] = r * np.sin(t)
            y[:] = r * np.cos(t)
            return xy

        transform_non_affine = transform

        def inverted(self):
            return NorthPolarAxes.InvertedNorthPolarTransform()

    class InvertedNorthPolarTransform(PolarAxes.InvertedPolarTransform):
        def transform(self, xy):
            x = xy[:, 0:1]
            y = xy[:, 1:]
            r = np.sqrt(x * x + y * y)
            theta = np.arctan2(y, x)
            return np.concatenate((theta, r), 1)

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
            Affine2D().scale(np.pi * 2.0, 1.0) +
            self.transData)
        self._yaxis_text1_transform = (
            self._r_label1_position +
            Affine2D().scale(1.0 / 360.0, 1.0) +
            self._yaxis_transform)


register_projection(NorthPolarAxes)

deg2rad = np.pi / 180.


def _getsweeptime(path, CFRadial=True):
    """Attempts to read sweep time from the radar file or construct it from the file name."""
    if(CFRadial):
        # print "Opening file: ",path
        sweepfile_netcdf = Nio.open_file(path)
        sweeptime = sweepfile_netcdf.variables['time_coverage_start'].get_value().tostring()
        radyear = int(sweeptime[:4])
        radmonth = int(sweeptime[5:7])
        radday = int(sweeptime[8:10])
        radhour = int(sweeptime[11:13])
        radmin = int(sweeptime[14:16])
        radsec = int(sweeptime[17:19])
        sweepfile_netcdf.close()
    else:
        filename = os.path.split(path)[1]   # Get time from file name for now
        sweeptime = filename[1:-4]
        radyear = int(sweeptime[:4])
        radmonth = int(sweeptime[4:6])
        radday = int(sweeptime[6:8])
        radhour = int(sweeptime[8:10])
        radmin = int(sweeptime[10:12])
        radsec = int(sweeptime[12:14])

    # -timedelta(hours=5) # Get it to CST time
    return datetime(radyear, radmonth, radday, radhour, radmin, radsec)


def getradarfilelist(radar_dir, radar_save_dir, starttime=None, stoptime=None, platform='NEXRAD',
                     radar_name=None, el_req=0.5):
    """
    Given a directory of radar data in netCDF format, return a list of the files containing
    successive sweeps at the desired elevation angle.  Currently only CFRadial files for NEXRAD
    and (proprietary) netCDF files for UMASS X-Pol are supported.
    """

    if(platform not in ['NEXRAD', 'UMXP', 'SMARTR']):
        print("Sorry, platform: ", platform, " not supported.  Please try again!")
        return

    if(platform == 'NEXRAD' or platform == 'SMARTR'):
        CFRadial = True
    else:
        CFRadial = False

    if(platform == 'UMXP' or platform == 'SMARTR'):
        radar_name = platform

    radar_filelist = glob.glob(radar_dir + '*.nc')

    # Filter only those files containing sweeps closest to the desired elevation angle
    # Only need to do this for UMXP or SMARTR, which has separate files for each elevation
    if(platform == 'UMXP' or platform == 'SMARTR'):
        # First check to see if text file containing file list of requested
        # elevation already exists (saves loads of time)
        radar_filelist_filename = (radar_name + '_' + starttime.strftime(tm.timefmt3).strip() +
                                   '_' + stoptime.strftime(tm.timefmt3).strip() + '_el' +
                                   str(el_req) + '_filelist.txt')
        radar_filelist_path = os.path.join(radar_save_dir, radar_filelist_filename)
        if(os.path.exists(radar_filelist_path)):
            print(len(radar_filelist))
            with open(radar_filelist_path) as f:
                radar_filelist = [line.rstrip('\n') for line in f]

        else:
            radar_filelist = getncfilelist(platform, radar_filelist, el_req, tolerance=0.5)
            # Save the file list with requested elevation to text file for later retrieval
            radar_filelist_fileout = open(radar_filelist_path, 'w')
            for file in radar_filelist:
                print(file, file=radar_filelist_fileout)

    radtimes = []

    for path in radar_filelist:
        # Grab the time information from the file
        radtime = _getsweeptime(path, CFRadial)
        radtimes.append(radtime)

        # print "Time of sweep = ",radtime.strftime(fmt)

    # Sort the list of radar files using the sweep times
    radar_filelist_sorted = [path for (radtime, path) in sorted(zip(radtimes, radar_filelist),
                                                                key=lambda pair: pair[0])]
    radtimes.sort()

    if(not starttime):
        starttime = radtimes[0]
    if(not stoptime):
        stoptime = radtimes[-1]

    # Find times within desired range
    try:
        radstartindex = next(i for i, radtime in enumerate(radtimes) if radtime >= starttime)
    except Exception:
        radstartindex = 0
    try:
        radstopindex = next(i for i, radtime in enumerate(radtimes) if radtime >= stoptime)
    except Exception:
        radstopindex = len(radtimes)

    radar_filelist = radar_filelist_sorted[radstartindex:radstopindex + 1]
    radtimes = radtimes[radstartindex:radstopindex + 1]

    return radar_filelist, radtimes


def readCFRadial_pyART(el, filename, sweeptime, fieldnames, radlat=None, radlon=None,
                       radalt=None, compute_kdp=True):
    """Reads radar data from a CFRadial netCDF file.  Attempts to extract fields given by
       the input list "fieldnames".  For nexrad files, which contain an entire volume
       (apparently), only the desired elevation angle will be returned.  The fields are
       returned in numpy arrays dimensioned by (range,azimuth), along with other details
       about the radar. This version utilizes the pyART libraries."""

    fieldlist = []

    print("Opening file: ", filename)
    radarobj = pyart.io.read_cfradial(filename)

    elevs = radarobj.elevation['data']

    # Try to match up requested elevation angle with one in the file
    # Pick closest elevation angle to that requested

    elev_list = []
    for elev in radarobj.iter_elevation():
        elev_list.append(elev[0])  # Just grab the elevation of the first ray in each sweep

    # print "elev_list",elev_list

    elev_list = np.array(elev_list)

    print("Requested elevation angle ", el)

    # Returns index of sweep with closest elevation angle to that requested
    sweepindex = (np.abs(elev_list - el)).argmin()
    radarsweep = radarobj.extract_sweeps([sweepindex])

    el = elevs[sweepindex]
    # print elevs[swp_start_index:swp_end_index]
    print("Actual elevation angle at start of sweep: ", el)
    # print "swp_start_index,swp_end_index",swp_start_index,swp_end_index

    # Grab the time information from the file
    if(not sweeptime):
        sweeptime = _getsweeptime(filename, True)

    print("Time of sweep = ", sweeptime.strftime(tm.timefmt))

    # Assign some needed variables from the PyART radar object

    numgates = radarsweep.ngates              # Number of gates
    # sweep_start_ray_index and sweep_end_ray_index contain the time index of the start and
    # end of each sweep (i.e. at one elevation angle) of radar data

#     try:
#         ray_n_gates = sweepfile_netcdf.variables['ray_n_gates'][:]
#         twoDarray = False
#         # ray_start_index contains the data indices of the start of each ray of radar data
#         # needed below to reshape the 1D data arrays into 2D arrays by azimuth and range
#         ray_start_index = sweepfile_netcdf.variables['ray_start_index'][:]
#     except:
#         print "No ray_n_gates in file, assuming 2D arrays (azimuth,range)"
#         twoDarray = True

    if(radlat is None):
        rlat = radarsweep.latitude['data']
        rlon = radarsweep.longitude['data']
        ralt = radarsweep.altitude['data']

        if(np.size(rlat) > 1):
            rlat = rlat[0]
        if(np.size(rlon) > 1):
            rlon = rlon[0]
        if(np.size(ralt) > 1):
            ralt = ralt[0]
    else:
        rlat = radlat
        rlon = radlon
        ralt = radalt

    print("Number of gates: ", numgates)
    print("Radar lat,lon,alt", rlat, rlon, ralt)

    range = radarsweep.range['data']
    gatewidth = range[1] - range[0]  # gate spacing

    print("Gatewidth ", gatewidth)

    # Get the Azimuth info

    beam_width = radarsweep.instrument_parameters['radar_beam_width_h']['data'][0]
    print("Radar beam width (degrees): " + str(beam_width))

    num_azim = radarsweep.nrays

    print("Number of azimuths in sweep ", num_azim)
    # print azimuth_rad

    # Roll the azimuth dimension around so that it starts at 0 and ends at 360

#     try:
#         shift = np.where(azimuth_rad < azimuth_rad[0])[0][0]
#     except:
#         shift = 0   # Not sure if this is correct
#
#     #print "shift = ",shift
#
#     azimuth_rad = np.roll(azimuth_rad,shift=-shift)

#     beamwidth = azimuth_rad[1:]-azimuth_rad[0:-1]
#     # For plotting we need the azimuth at the borders of each ray (start of each ray plus one
#     # more on the end)
#     azimuth_start = np.zeros((np.size(azimuth)+1))
#
#     # Find azimuth of "start" of each gate (in azimuth) -- approximate
#     azimuth_start[1:-1] = azimuth[1:]-0.5*beamwidth[:]
#     azimuth_start_rad[0] = azimuth_rad[0]-0.5*beamwidth[0]
#     azimuth_start_rad[-1] = azimuth_rad[-1]+0.5*beamwidth[-1]
#
#     azimuth_rad = azimuth_rad*deg2rad
#     azimuth_start_rad = azimuth_start_rad*deg2rad

    # Read the desired fields from the file, if they exist
    # Need to find better solution than nested try/excepts

    # print fieldnames

    # outfieldnames = []
    outfieldnames = {}

    for fieldname in fieldnames:
        f = None
        fieldtype = 'unknown'
        # Reflectivity
        if fieldname in ['dBZ', 'DBZ', 'Z']:
            fieldtype = 'reflectivity'
            f = next((f for f in list(radarsweep.fields.items())
                      if f[0] in ['REF', 'DZ']), None)
        if fieldname in ['Zdr', 'ZDR']:
            fieldtype = 'differential reflectivity'
            f = next((f for f in list(radarsweep.fields.items())
                      if f[0] in ['ZDR', 'DB_ZDR']), None)
        if fieldname in ['Kdp', 'KDP']:
            fieldtype = 'specific differential phase'
            f = next((f for f in list(radarsweep.fields.items())
                      if f[0] in ['KD']), None)
            # In this case, if Kdp isn't found, we may want to additionally check if PhiDP is in
            # the file, and if it is, compute Kdp using the PyART built-in linear programming method
            if(not f and compute_kdp):
                tempfield = next(
                    (field_data['data'] for field_key,
                     field_data in list(radarsweep.fields.items()) if field_key in ['PHI']),
                    np.empty(
                        (0)))
                if(tempfield.size):
                    # First check that normalized coherent power is in the file.  If not, set to
                    # all ones (Scott Collis, personal communication 2016). It seems the level-II
                    # data doesn't have this field.
                    ncp_exists = any(
                        x in list(radarsweep.fields.keys()) for x in [
                            'normalized_coherent_power', 'SQI', 'SQI2', 'NCP', 'NCP_F'])
                    if(not ncp_exists):
                        ncp_field = pyart.config.get_field_name('normalized_coherent_power')
                        ncp = pyart.config.get_metadata(ncp_field)
                        ncp['data'] = np.ones_like(tempfield)
                        radarsweep.add_field('normalized_coherent_power', ncp)
                    print("Computing specific differential phase from differential phase.")
                    refl_field_name = next((field_key for field_key in radarsweep.fields
                                            if field_key in ['REF', 'DZ']), None)
                    rhv_field_name = next((field_key for field_key in radarsweep.fields
                                           if field_key in ['RHO']), None)
                    if(refl_field_name and rhv_field_name):
                        phidp, kdp = pyart.correct.phase_proc_lp(
                            radarsweep, 0.0, refl_field=refl_field_name, rhv_field=rhv_field_name,
                            phidp_field='PHI', debug=True)
                        f = (fieldname, kdp)
                        radarsweep.add_field(fieldname, kdp)
                    else:
                        print("Sorry, couldn't compute KDP!")
                        filefieldname = ''
                        filefield = np.empty((0))
        if fieldname in ['rhv', 'RHV']:
            fieldtype = 'cross correlation coefficient'
            f = next((f for f in list(radarsweep.fields.items())
                      if f[0] in ['RHO']), None)
        if fieldname in ['vr', 'VR', 'Vr']:
            fieldtype = 'radial velocity'
            f = next((f for f in list(radarsweep.fields.items())
                      if f[0] in ['VEL', 'VR']), None)
        if(not f):
            fieldname = ''
            filefield = np.empty((0))
        else:
            filefieldname = f[0]
            filefield = f[1]['data']
        if(not filefield.size and fieldtype != 'unknown'):
            print("Cannot find " + fieldtype + " field in file, setting to to empty")
        elif(fieldtype == 'unknown'):
            print("Unknown or unimplemented field: " + fieldname)
        if(filefield.size):
            outfieldnames[fieldname] = filefieldname
            fieldlist.append(filefield)

    # Return stuff

    # return outfieldnames,fieldlist,range,azimuth,rlat_rad,rlon_rad,ralt,el_rad
    return outfieldnames, radarsweep


def readCFRadial(nexrad, el, radlat, radlon, radalt, file, sweeptime, fieldnames):
    """Reads radar data from a CFRadial netCDF file.  Attempts to extract fields given by
       the input list "fieldnames".  For nexrad files, which contain an entire volume
       (apparently), only the desired elevation angle will be returned.  The fields are
       returned in numpy arrays dimensioned by (range,azimuth), along with other details
       about the radar."""

    fieldlist = []

    print("Opening file: ", file)
    sweepfile_netcdf = Nio.open_file(file)

    # Grab the time information from the file

    print("Time of sweep = ", sweeptime.strftime(tm.timefmt))

    # Grab some needed variables from the netCDF file

    numgates = sweepfile_netcdf.dimensions['range']     # Number of gates
    numtimes = sweepfile_netcdf.dimensions['time']     # Number of times
    # numsweeps = sweepfile_netcdf.dimensions['sweep']     # Number of sweeps
    # sweep_start_ray_index and sweep_end_ray_index contain the time index of the start and
    # end of each sweep (i.e. at one elevation angle) of radar data
    sweep_start_ray_index = sweepfile_netcdf.variables['sweep_start_ray_index'][:]
    sweep_end_ray_index = sweepfile_netcdf.variables['sweep_end_ray_index'][:]

    try:
        # ray_n_gates = sweepfile_netcdf.variables['ray_n_gates'][:]
        twoDarray = False
        # ray_start_index contains the data indices of the start of each ray of radar data
        # needed below to reshape the 1D data arrays into 2D arrays by azimuth and range
        ray_start_index = sweepfile_netcdf.variables['ray_start_index'][:]
    except Exception:
        print("No ray_n_gates in file, assuming 2D arrays (azimuth,range)")
        twoDarray = True

    if(radlat is None):
        rlat = sweepfile_netcdf.variables['latitude'].get_value()
        rlon = sweepfile_netcdf.variables['longitude'].get_value()
        ralt = sweepfile_netcdf.variables['altitude'].get_value()

        if(np.size(rlat) > 1):
            rlat = rlat[0]
        if(np.size(rlon) > 1):
            rlon = rlon[0]
        if(np.size(ralt) > 1):
            ralt = ralt[0]
    else:
        rlat = radlat
        rlon = radlon
        ralt = radalt

    print("Number of gates: ", numgates)
    print("Radar lat,lon,alt", rlat, rlon, ralt)

    rlat_rad = rlat * deg2rad
    rlon_rad = rlon * deg2rad

    elevs = sweepfile_netcdf.variables['elevation'][:]

    # Try to match up requested elevation angle with one in the file
    # Pick closest elevation angle to that requested

    elev_list = []
    for swp_index in sweep_start_ray_index:
        elev_list.append(elevs[swp_index])

    # print "elev_list",elev_list

    elev_list = np.array(elev_list)

    print("Requested elevation angle ", el)

    index = (np.abs(elev_list - el)).argmin()    # Returns closest elevation angle to that requested
    # print "index = ",index
    swp_start_index = sweep_start_ray_index[index]
    swp_end_index = sweep_end_ray_index[index] + 1
    # DTD, some sweeps seem to have problems where there aren't enough times in the file so limit
    # the swp_end_index here to the number of times available.
    swp_end_index = min(swp_end_index, numtimes - 1)
    el = elevs[swp_start_index]
    # print elevs[swp_start_index:swp_end_index]
    print("Actual elevation angle at start of sweep: ", el)
    # print "swp_start_index,swp_end_index",swp_start_index,swp_end_index
    el_rad = el * deg2rad

    temp = sweepfile_netcdf.variables['range']    # range to center of each gate
    gatewidth = temp.meters_between_gates[0]  # gate spacing for each gate
    range = temp[:]
    range_start = temp[:] - gatewidth / 2.  # We also want range to start of each gate

    print("Gatewidth ", gatewidth)

    # Get the Azimuth info

    beam_width = sweepfile_netcdf.variables['radar_beam_width_h'].get_value()
    print("Radar beam width (degrees): " + str(beam_width))

    azimuth_rad = sweepfile_netcdf.variables['azimuth'][:]  # Azimuth of center of each gate
    azimuth_rad = azimuth_rad[swp_start_index:swp_end_index]
    num_azim = np.size(azimuth_rad)

    print("Number of azimuths in sweep ", num_azim)
    # print azimuth_rad

    # Roll the azimuth dimension around so that it starts at 0 and ends at 360

    try:
        shift = np.where(azimuth_rad < azimuth_rad[0])[0][0]
    except Exception:
        shift = 0   # Not sure if this is correct

    # print "shift = ",shift

    azimuth_rad = np.roll(azimuth_rad, shift=-shift)
    beamwidth = azimuth_rad[1:] - azimuth_rad[0:-1]

    # For plotting we need the azimuth at the borders of each ray (start of
    # each ray plus one more on the end)
    azimuth_start_rad = np.zeros((np.size(azimuth_rad) + 1))
    # Find azimuth of "start" of each gate (in azimuth) -- approximate
    azimuth_start_rad[1:-1] = azimuth_rad[1:] - 0.5 * beamwidth[:]
    azimuth_start_rad[0] = azimuth_rad[0] - 0.5 * beamwidth[0]
    azimuth_start_rad[-1] = azimuth_rad[-1] + 0.5 * beamwidth[-1]

    azimuth_rad = azimuth_rad * deg2rad
    azimuth_start_rad = azimuth_start_rad * deg2rad

    # Read the desired fields from the file, if they exist
    # Need to find better solution than nested try/excepts

    # print fieldnames

    outfieldnames = []
    for field in fieldnames:
        tempfield = None
        if(field == 'dBZ' or field == 'Z'):
            try:
                tempfield = sweepfile_netcdf.variables['REF']
            except Exception:
                try:
                    tempfield = sweepfile_netcdf.variables['DZ']
                except Exception:
                    print("Cannot find reflectivity field in file, setting to to empty")
                    tempfield = np.empty((0))
        if(field == 'Zdr' or field == 'ZDR'):
            try:
                tempfield = sweepfile_netcdf.variables['ZDR']
            except Exception:
                try:
                    tempfield = sweepfile_netcdf.variables['DB_ZDR']
                except Exception:
                    print("Cannot find differential reflectivity field in file, ",
                          "setting to to empty")
                    tempfield = np.empty((0))
        if(field == 'Kdp' or field == 'KDP'):
            try:
                tempfield = sweepfile_netcdf.variables['KD']
            except Exception:
                print("Cannot find Kdp field in file, setting to to empty")
                tempfield = np.empty((0))
        if(field == 'rhv' or field == 'RHV'):
            try:
                tempfield = sweepfile_netcdf.variables['RHO']
            except Exception:
                print("Cannot find rHV field in file, setting to to empty")
                tempfield = np.empty((0))
        if(field == 'Vr' or field == 'VR'):
            try:
                tempfield = sweepfile_netcdf.variables['VEL']
            except Exception:
                try:
                    tempfield = sweepfile_netcdf.variables['VR']
                except Exception:
                    print("Cannot find radial velocity field in file, setting to to empty")
                    tempfield = np.empty((0))

        if(tempfield):
            outfieldnames.append(field)

        # Grab requested sweep out of dBZ array
        # Need to check if the array is dimensioned by total number of gates or just
        # number of gates for each sweep

        if(tempfield):
            if(not twoDarray):  # Note, this needs to be tested!
                # numpointsinsweep = ray_n_gates[swp_start_index]*num_azim
                # print numpointsinsweep
                field_start_index = ray_start_index[swp_start_index]
                field_end_index = ray_start_index[swp_end_index]
                # field = tempfield[swp_start_index:numpointsinsweep+swp_start_index]
                field = tempfield[field_start_index:field_end_index]
                # print field.shape
                field = field.reshape((num_azim, -1))
                # print field.shape
            else:
                field = tempfield[swp_start_index:swp_end_index, :]

            # Unpack values
            # print tempfield.scale_factor[0]
            # print tempfield.add_offset[0]
            field = field * tempfield.scale_factor[0] + tempfield.add_offset[0]

            # Adjust azimuth axis
            field = np.roll(field, shift=-shift, axis=0)
            # Add field to list for return
            fieldlist.append(field)

    # Return stuff

    return (outfieldnames, fieldlist, range_start, range, azimuth_start_rad, azimuth_rad,
            rlat_rad, rlon_rad, ralt, el_rad)


def getncfilelist(platform, filelist, el_req, tolerance=0.5):
    """Given a list of radar netcdf files that contain one elevation angle apiece,
       find those that contain the elevation angle closest
       to that requested and return a list of them"""

    if(platform == 'UMXP'):
        elvarname = 'Elevation'
    else:
        elvarname = 'elevation'

    elevations = []

    for index, file in enumerate(filelist):
        sweepfile_netcdf = Nio.open_file(file)
        elevations.append(sweepfile_netcdf.variables[elvarname][0])

    elevations = np.array(elevations)
    print(elevations)
    diffs = np.abs(elevations - el_req)
    indices = np.where(diffs < tolerance)[0].tolist()
    newfilelist = [filelist[i] for i in indices]

    return newfilelist


def readUMXPnc(file, sweeptime, fieldnames, heading=None, correct_attenuation=True):
    """Reads a single sweep from a UMASS X-pol netCDF file"""

    debug = False

    fieldlist = []

    print("Opening file: ", file)
    sweepfile_netcdf = Nio.open_file(file)

    # print "Time of sweep = ",sweeptime.strftime(tm.fmt)

    # Grab some needed variables from the netCDF file

    numgates = sweepfile_netcdf.dimensions['Gate']
    if not heading:
        heading = sweepfile_netcdf.Heading  # -90.0
    print("Heading: ", heading)
    heading = heading * np.pi / 180.
    el = sweepfile_netcdf.variables['Elevation'][0]
    print("Elevation angle: ", el)
    el_rad = el * deg2rad

    rlat = sweepfile_netcdf.Latitude[0]
    rlon = sweepfile_netcdf.Longitude[0]
    ralt = sweepfile_netcdf.Altitude[0]

    print("Number of gates: ", numgates)
    print("Radar lat,lon,alt", rlat, rlon, ralt)

    rlat_rad = rlat * deg2rad
    rlon_rad = rlon * deg2rad

    gatewidth = sweepfile_netcdf.variables['GateWidth'][:]

    # Gatewidth appears to be the same everywhere, so only need first element
    gatewidth = gatewidth[0]
    print("Gatewidth ", gatewidth)

    # Set up the range 1-D array

    range_start = np.linspace(0.0, (numgates - 1) * gatewidth, num=numgates)
    range = range_start + gatewidth / 2.

    # Get the Azimuth info

    beam_width = sweepfile_netcdf.AntennaBeamwidth[0]
    print("Radar beam width (degrees): " + str(beam_width))

    azimuth = sweepfile_netcdf.variables['Azimuth'][:]

    # Some of the UMXP files have extra sweeps in them (or have bad azimuths as the azimuth wraps
    # around 360 deg).  This is a stopgap measure to remove the extra sweeps while I get Will's
    # code to work correctly to split the netCDF files that have extra sweeps up, or Will uploads
    # the reprocessed data, whichever comes first.
    # UPDATE 08/18/2017: Something isn't working correctly with this code at
    # the moment. Disabled for now.

#     diffazim = azimuth[1:]-azimuth[:-1]
#     stopindices = np.where(diffazim < -90.)[0]
#     if(stopindices.size):
#         if(stopindices[0] == 0):
#             stopindex = stopindices[1]
#         else:
#             stopindex = stopindices[0]
#     else:
#         stopindex = azimuth.size-1
#
#     print "stopindex = ",stopindex
#     azimuth = azimuth[:stopindex+1]
    # Add heading of truck to azimuth
    azimuth_rad = azimuth * np.pi / 180
    num_azim = np.size(azimuth_rad)
    print("Number of azimuths in sweep ", num_azim)

    # Add heading of truck to azimuth
    azimuth_rad = azimuth_rad + heading
    # Actually difference between successive azimuths: may be different from half-power beam width
    beamwidth = azimuth_rad[1:] - azimuth_rad[0:-1]
    azimuth_rad = azimuth_rad % (2. * np.pi)  # wraps angles greater than 2.*pi back around 0

    # For plotting we need the azimuth at the borders of each ray (start of
    # each ray plus one more on the end)
    azimuth_start_rad = np.zeros((np.size(azimuth_rad) + 1))
    # Find azimuth of "start" of each gate (in azimuth) -- approximate
    azimuth_start_rad[1:-1] = azimuth_rad[1:] - 0.5 * beamwidth[:]
    azimuth_start_rad[0] = azimuth_rad[0] - 0.5 * beamwidth[0]
    azimuth_start_rad[-1] = azimuth_rad[-1] + 0.5 * beamwidth[-1]
    # wraps angles greater than 2.*pi back around 0
    azimuth_start_rad = azimuth_start_rad % (2. * np.pi)

    # Read the desired fields from the file, if they exist
    # print fieldnames

    outfieldnames = []
    for fieldname in fieldnames:
        tempfield = None
        if(fieldname == 'dBZ' or fieldname == 'Z'):
            try:
                tempfield = sweepfile_netcdf.variables['Zh'][:]
                if(correct_attenuation):
                    try:
                        tempfield = tempfield + sweepfile_netcdf.variables['Ah'][:]
                    except Exception:
                        print("No attenuation correction information in the file!")
            except Exception:
                print("Cannot find reflectivity field in file, setting to to empty")
                tempfield = np.empty((0))
        if(fieldname == 'Zdr' or fieldname == 'ZDR'):
            try:
                tempfield = sweepfile_netcdf.variables['ZDR'][:]
                if(correct_attenuation):
                    try:
                        tempfield = tempfield + sweepfile_netcdf.variables['Ad'][:]
                    except Exception:
                        print("No attenuation correction information in the file!")
            except Exception:
                print("Cannot find differential reflectivity field in file, setting to to empty")
                tempfield = np.empty((0))
        if(fieldname == 'Kdp' or fieldname == 'KDP'):
            try:
                tempfield = sweepfile_netcdf.variables['KD'][:]
            except Exception:
                print("Cannot find Kdp field in file, setting to to empty")
                tempfield = np.empty((0))
        if(fieldname == 'rhv' or fieldname == 'RHV'):
            try:
                tempfield = sweepfile_netcdf.variables['RhoHV'][:]
            except Exception:
                print("Cannot find rHV field in file, setting to to empty")
                tempfield = np.empty((0))
        if(fieldname == 'Vr' or fieldname == 'VR'):
            try:
                tempfield = sweepfile_netcdf.variables['VE'][:]
            except Exception:
                print("Cannot find radial velocity field in file, setting to to empty")
                tempfield = np.empty((0))

        if(tempfield.size):
            outfieldnames.append(fieldname)

        # Grab requested sweep out of dBZ array
        # Need to check if the array is dimensioned by total number of gates or just number of
        # gates for each sweep

        if(tempfield.size):
            field = tempfield
            # Set missing values to np.nan
            field = np.where(field == -99, np.nan, field)
            # Add field to list for return
            fieldlist.append(field)

        if(debug):
            # Quick plot for testing

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, polar=True)
            ax1.set_theta_zero_location('N')
            ax1.set_theta_direction(-1)

            ax1.pcolor(azimuth_start_rad, range_start, field.T)

    # Get rid of negative infinity values (set them to 0)
#
#     Zh = np.where(np.isneginf(Zh),0.0,Zh)
#     Zh = np.where(np.isposinf(Zh),0.0,Zh)
#     Zh = np.where(np.isnan(Zh),0.0,Zh)
#
#     Zdr = np.where(np.isneginf(Zdr),0.0,Zdr)
#     Zdr = np.where(np.isposinf(Zdr),0.0,Zdr)
#     Zdr = np.where(np.isnan(Zdr),0.0,Zdr)
#

    # Return stuff

    return (outfieldnames, fieldlist, range_start, range, azimuth_start_rad, azimuth_rad, rlat_rad,
            rlon_rad, ralt, el_rad)


def sweep2xy(azimuth_start_rad, azimuth_rad, range_start, range, el_rad, rlat_rad, rlon_rad, ralt,
             map_proj=1):
    """Brings a radar sweep from r,theta to x,y coordinates with the desired map projection."""
    theta, rad = np.meshgrid(azimuth_start_rad, range_start)
    theta_c, rad_c = np.meshgrid(azimuth_rad, range)
    xrad, yrad = oban.xyloc(rad, theta, el_rad, rlat_rad, rlon_rad, ralt, map_proj)
    xrad_c, yrad_c = oban.xyloc(rad_c, theta_c, el_rad, rlat_rad, rlon_rad, ralt, map_proj)

    return xrad, yrad, xrad_c, yrad_c


def plotsweep(radlims, plotlims, fieldnames, fieldlist, masklist, range_start, range,
              azimuth_start_rad, azimuth_rad, rlat_rad, rlon_rad, ralt, el_rad, ovrmap, ovrdis,
              dis_name_list, dxy_list, fields_D_list):
    """Plots an individual radar sweep, with an option to overlay a map.  This assumes the
       sweep has been read in and adjusted using the readCFRadial function"""

    # Plot approximate locations of storm at certain time intervals based on
    # assume storm motion? (no effect for now)
    # plot_interval = False
    # storm_u = 12.55
    # storm_v = 0.0
    # time_interval = 30.0
    # Number of intervals (plus and minus) to plot around each center time and associated scan
    # num_intervals = 2

    # Low reflectivity threshold
    # dBZthresh = 5.0

    # Unpack range and azimuth limits for plot
    minrange = radlims[0]
    maxrange = radlims[1]
    minazim = radlims[2]
    maxazim = radlims[3]

    # Unpack plotting limits
    plotxmin = plotlims[0]
    plotxmax = plotlims[1]
    plotymin = plotlims[2]
    plotymax = plotlims[3]

    # Find the indices corresponding to the plot limits
    minrangeindex = np.where(range_start > minrange)[0][0]
    try:
        maxrangeindex = np.where(range_start > maxrange)[0][0]
    except Exception:
        maxrangeindex = np.size(range_start)
    # print "minrangeindex,maxrangeindex",minrangeindex,maxrangeindex

    minazimindex = np.where(
        azimuth_start_rad /
        deg2rad >= min(
            azimuth_start_rad[0] /
            deg2rad,
            minazim))[0][0]
    # print "minazimindex,azimuth_start_rad[minazimindex] =
    # ",minazimindex,azimuth_start_rad[minazimindex]/deg2rad
    try:
        maxazimindex = np.where(azimuth_start_rad / deg2rad >= maxazim)[0][0]
    except Exception:
        maxazimindex = np.size(azimuth_start_rad) - 1

    # print "maxazimindex,azimuth_start_rad[maxazimindex] =
    # ",maxazimindex,azimuth_start_rad[maxazimindex]/deg2rad

    # print azimuth_start_rad/deg2rad
    # Set up the coordinate arrays
    theta, rad = np.meshgrid(
        azimuth_start_rad[minazimindex:maxazimindex + 1], range_start[minrangeindex:maxrangeindex])
    theta_c, rad_c = np.meshgrid(
        azimuth_rad[minazimindex:maxazimindex + 1], range[minrangeindex:maxrangeindex])

    # Plot on rectangular mesh using x,y points derived from r,theta
    # The points are valid on the edges of the gates, while the reflectivity values
    # are valid at the centers.  This is so that pcolor will plot the gates in the correct
    # locations (hopefully).

    xplt, yplt = oban.xyloc(rad, theta, el_rad, rlat_rad, rlon_rad, ralt, 1)
    xplt_c, yplt_c = oban.xyloc(rad_c, theta_c, el_rad, rlat_rad, rlon_rad, ralt, 1)

    figlist = []
    gridlist = []

    # Loop through the provided fields and prepare them for plotting
    for fieldname, field, fields_D, mask in zip(fieldnames, fieldlist, fields_D_list, masklist):
        if(field is None):
            continue

        # mask array by provided mask, if desired (usually some reflectivity threshold)
        if(mask is not None):
            field = np.ma.masked_array(field, mask=mask)

        fieldplt = field.swapaxes(0, 1)[minrangeindex:maxrangeindex, minazimindex:maxazimindex + 1]
        fieldplt = np.ma.masked_invalid(fieldplt)

        # Decide whether to perform a shift in time and space (disabled for now)
#         subtimes = []
#         if(plot_interval):
#             for dt in np.arange(-num_intervals*time_interval,
#                                 num_intervals*time_interval+time_interval,time_interval):
#                 subtime = sweeptime+timedelta(seconds=dt)
#                 subtimes.append(subtime)
#         else:
#             subtimes.append(sweeptime)
#
#         for dt,time in zip(np.arange(-num_intervals*time_interval,
#                                      num_intervals*time_interval+time_interval,time_interval),
#                                      subtimes):
#             if(plot_interval):
#                 xplttmp = xplt+dt*storm_u
#                 xplt_ctmp = xplt_c+dt*storm_u
#             else:
#                 xplttmp = xplt
#                 xplt_ctmp = xplt_c

        xplttmp = xplt
        yplttmp = yplt

        if(fieldname == 'dBZ'):
            clevels = clevels_ref
            cmap = cmapdBZ
            clvls = 5.0
            disfmtstr = "{:3.1f} dBZ"
        elif(fieldname == 'ZDR'):
            clevels = clevels_zdr
            cmap = cmapzdr
            clvls = 0.5
            disfmtstr = "{:3.1f} dB"
        elif(fieldname == 'KDP'):
            clevels = clevels_kdp
            cmap = cmapkdp
            clvls = 0.5
            disfmtstr = "{:3.1f} deg/km"
        elif(fieldname == 'RHV'):
            clevels = clevels_rhv
            cmap = cmaprhv
            clvls = 0.1
            disfmtstr = "{:3.1f}"
        elif(fieldname == 'Vr'):
            clevels = clevels_vr
            cmap = cmapvr
            clvls = 5.0
            disfmtstr = "{:3.1f} m/s"
        elif(fieldname == 'Rain'):
            clevels = clevels_rain
            cmap = cmapretrievals
            clvls = [0.1, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0]
            disfmtstr = "{:3.1f} mm/hr"
        elif(fieldname == 'D0'):
            clevels = clevels_d0
            cmap = cmapretrievals
            clvls = 0.5
            disfmtstr = "{:3.1f} mm"
        elif(fieldname == 'W'):
            clevels = clevels_lwc
            cmap = cmapretrievals
            clvls = [0.01, 0.1, 0.3, 0.5, 1.0, 1.5, 2., 3., 5., 7.5, 10.]
            disfmtstr = "{:3.1f} g m^-3"
        elif(fieldname == 'Nt'):
            clevels = clevels_nt
            cmap = cmapretrievals
            clvls = 0.5
            disfmtstr = "{:3.1f} # m^3"
        elif(fieldname == 'sigm'):
            clevels = clevels_sigm
            cmap = cmapretrievals
            clvls = 0.2
            disfmtstr = "{:3.1f} mm"
        elif(fieldname == 'mu'):
            clevels = clevels_mu
            cmap = cmapretrievals
            clvls = 2.0
            disfmtstr = "{:3.1f}"
        elif(fieldname == 'lam'):
            clevels = clevels_lam
            cmap = cmapretrievals
            clvls = 2.0
            disfmtstr = "{:3.1f}"

        fig = plt.figure()
        grid = ImageGrid(
            fig,
            111,
            nrows_ncols=(
                1,
                1),
            axes_pad=0.0,
            cbar_mode="single",
            cbar_location="right",
            aspect=True,
            label_mode="1")
        ax = grid[0]
        if(fieldname == 'Rain'):
            norm = matplotlib.colors.LogNorm(vmin=0.1, vmax=200.0)
        elif(fieldname == 'W'):
            norm = matplotlib.colors.LogNorm(vmin=0.01, vmax=10.)
        elif(fieldname == 'dBZ'):
            norm = normdBZ
        else:
            norm = matplotlib.colors.BoundaryNorm(clevels, cmap.N)
        plot1 = ax.pcolormesh(xplttmp,
                              yplttmp,
                              fieldplt,
                              vmin=clevels[0],
                              vmax=clevels[-1],
                              cmap=cmap,
                              norm=norm,
                              edgecolors='None',
                              antialiased=False)
        # plot2 = ax2.contour(xplt_c,yplt_c,dBZplt,levels=[40.0],c='k')

        # DTD: originally did this outside the function, but have to do it here
        # because of some weird problem
        if(ovrdis):
            for dname, dloc, field_D in zip(dis_name_list, dxy_list, fields_D):
                Dx = dloc[0]
                Dy = dloc[1]
                ax.plot(Dx, Dy, 'k*', ms=8)
                ax.annotate(dname, (Dx + 550., Dy - 100.), clip_on=True)
                if(field_D is not None):
                    ax.annotate(disfmtstr.format(field_D), (Dx + 550., Dy - 2100.), clip_on=True)

        # Determine plot bounds
        if(plotxmin == -1):
            plotxmin = np.min(xplt)
        if(plotxmax == -1):
            plotxmax = np.max(xplt)
        if(plotymin == -1):
            plotymin = np.min(yplt)
        if(plotymax == -1):
            plotymax = np.max(yplt)

        # FIXME: something is screwed up with ImageGrid colorbar ticks, but it looks like the
        # correct behavior can be recovered by reverting to the matplotlib version. (Already
        # corrected in plotsweep_pyART)
        if(fieldname == 'Rain'):
            clvllocator = ticker.LogLocator(base=10.0, subs=(0.1, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0,
                                                             30.0, 50.0, 75.0, 100.0, 150.0, 200.0))
            grid.cbar_axes[0].colorbar(plot1)
            grid.cbar_axes[0].toggle_label(True)
            grid.cbar_axes[0].yaxis.set_ticks(
                [0.1, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0])
            grid.cbar_axes[0].yaxis.set_ticklabels(['0.1', '1.0', '3.0', '5.0', '10.0', '15.0',
                                                    '20.0', '30.0', '50.0', '75.0', '100.0',
                                                    '150.0', '200.0'])
        elif(fieldname == 'W'):
            clvllocator = ticker.LogLocator(base=10.0, subs=(0.01, 0.1, 0.3, 0.5, 1.0, 1.5, 2., 3.,
                                                             5., 7.5, 10.))
            grid.cbar_axes[0].colorbar(plot1)
            grid.cbar_axes[0].toggle_label(True)
            grid.cbar_axes[0].yaxis.set_ticks(
                [0.01, 0.1, 0.3, 0.5, 1.0, 1.5, 2., 3., 5., 7.5, 10.])
            grid.cbar_axes[0].yaxis.set_ticklabels(
                ['0.01', '0.1', '0.3', '0.5', '1.0', '1.5', '2.', '3.', '5.', '7.5', '10.'])
        else:
            clvllocator = ticker.MultipleLocator(base=clvls)
            grid.cbar_axes[0].colorbar(plot1)
            grid.cbar_axes[0].toggle_label(True)
            grid.cbar_axes[0].yaxis.set_major_locator(clvllocator)
        # plt.title('dBZ at el = %.2f'%(el/deg2rad)+'and time '+time.strftime(fmt))
        ax.set_xlim(plotxmin, plotxmax)
        ax.set_ylim(plotymin, plotymax)
        formatter = ticker.FuncFormatter(mtokm)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=10000.0))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=10000.0))
        ax.set_aspect('equal')

        figlist.append(fig)
        gridlist.append(grid)

    return figlist, gridlist


def plotsweep_pyART(radlims, plotlims, fieldnames, radarsweep, ovrmap, ovrdis, dis_name_list,
                    dxy_list, fields_D_list, xoffset=0.0, yoffset=0.0, alpha=0.5):
    """
    Plots an individual radar sweep, with an option to overlay a map.  This version uses pyART
    routines and assumes the radar sweep has been read in using readCFRadial_pyART.
    """

    # STOPPED HERE!

    # Plot approximate locations of storm at certain time intervals based on
    # assume storm motion? (no effect for now)
    # plot_interval = False
    # storm_u = 12.55
    # storm_v = 0.0
    # time_interval = 30.0
    # Number of intervals (plus and minus) to plot around each center time and associated scan
    # num_intervals = 2

    # Low reflectivity threshold
    # dBZthresh = 5.0

    # Unpack range and azimuth limits for plot
    minrange = radlims[0]
    maxrange = radlims[1]
    minazim = radlims[2]
    maxazim = radlims[3]

    # Unpack plotting limits
    plotxmin = plotlims[0]
    plotxmax = plotlims[1]
    plotymin = plotlims[2]
    plotymax = plotlims[3]

    # Find the indices corresponding to the plot limits
    range = radarsweep.range['data']
    minrangeindex = np.where(range > minrange)[0][0]
    try:
        maxrangeindex = np.where(range > maxrange)[0][0]
    except Exception:
        maxrangeindex = np.size(range)
    # print "minrangeindex,maxrangeindex",minrangeindex,maxrangeindex

    azimuth = radarsweep.azimuth['data']
    minazimindex = np.where(azimuth >= min(azimuth[0], minazim))[0][0]
    # print "minazimindex,azimuth_start_rad[minazimindex] =
    # ",minazimindex,azimuth_start_rad[minazimindex]/deg2rad
    try:
        maxazimindex = np.where(azimuth >= maxazim)[0][0]
    except Exception:
        maxazimindex = np.size(azimuth) - 1

    # print "maxazimindex,azimuth_start_rad[maxazimindex] =
    # ",maxazimindex,azimuth_start_rad[maxazimindex]/deg2rad

    # print azimuth_start_rad/deg2rad
    # Set up the coordinate arrays
    # theta,rad = np.meshgrid(azimuth_start_rad[minazimindex:maxazimindex+1],
    #                         range_start[minrangeindex:maxrangeindex])
    # theta_c,rad_c = np.meshgrid(azimuth_rad[minazimindex:maxazimindex+1],
    #                             range[minrangeindex:maxrangeindex])

    # Plot on rectangular mesh using x,y points derived from r,theta
    # The points are valid on the edges of the gates, while the reflectivity values
    # are valid at the centers.  This is so that pcolor will plot the gates in the correct
    # locations (hopefully).

    xplt, yplt, zplt = radarsweep.get_gate_x_y_z(0, edges=True)
    xplt = xplt[minazimindex:maxazimindex + 1, minrangeindex:maxrangeindex]
    yplt = yplt[minazimindex:maxazimindex + 1, minrangeindex:maxrangeindex]
    xplt = xplt + xoffset
    yplt = yplt + yoffset

    figlist = []
    gridlist = []

    # Loop through the provided fields and prepare them for plotting

    for fieldname, filefieldname in list(fieldnames.items()):
        #
        #         # mask array by provided mask, if desired (usually some reflectivity threshold)
        #         if(mask != None):
        #             field = np.ma.masked_array(field,mask=mask)

        field = radarsweep.fields[filefieldname]['data']
        fieldplt = field[minazimindex:maxazimindex + 1, minrangeindex:maxrangeindex]
        fieldplt = np.ma.masked_invalid(fieldplt)

        # Decide whether to perform a shift in time and space (disabled for now)
#         subtimes = []
#         if(plot_interval):
#             for dt in np.arange(-num_intervals*time_interval,
#                                 num_intervals*time_interval+time_interval,time_interval):
#                 subtime = sweeptime+timedelta(seconds=dt)
#                 subtimes.append(subtime)
#         else:
#             subtimes.append(sweeptime)
#
#         for dt,time in zip(np.arange(-num_intervals*time_interval,
#                                      num_intervals*time_interval+time_interval,
#                                      time_interval),subtimes):
#             if(plot_interval):
#                 xplttmp = xplt+dt*storm_u
#                 xplt_ctmp = xplt_c+dt*storm_u
#             else:
#                 xplttmp = xplt
#                 xplt_ctmp = xplt_c

        norm = None
        if fieldname in ['dBZ', 'DBZ', 'Z']:
            clevels = clevels_ref
            cmap = cmapdBZ
            norm = normdBZ
            clvls = 5.0
            # disfmtstr = "{:3.1f} dBZ"
        if fieldname in ['Zdr', 'ZDR']:
            clevels = clevels_zdr
            cmap = cmapzdr
            clvls = 0.5
            # disfmtstr = "{:3.1f} dB"
        if fieldname in ['Kdp', 'KDP']:
            clevels = clevels_kdp
            cmap = cmapkdp
            clvls = 0.5
            # disfmtstr = "{:3.1f} deg/km"
        if fieldname in ['rhv', 'RHV']:
            clevels = clevels_rhv
            cmap = cmaprhv
            clvls = 0.1
            # disfmtstr = "{:3.1f}"
        if fieldname in ['vr', 'VR', 'Vr', 'VEL']:
            clevels = clevels_vr
            cmap = cmapvr
            clvls = 5.0
            # disfmtstr = "{:3.1f} m/s"

        fig = plt.figure()
        grid = ImageGrid(
            fig,
            111,
            nrows_ncols=(
                1,
                1),
            axes_pad=0.0,
            cbar_mode="single",
            cbar_location="right",
            aspect=True,
            label_mode="1")
        ax = grid[0]
        if norm is None:
            norm = matplotlib.colors.BoundaryNorm(clevels, cmap.N)
        plot1 = ax.pcolormesh(xplt, yplt, fieldplt, vmin=clevels[0], vmax=clevels[-1], cmap=cmap,
                              norm=norm, edgecolors='None', antialiased=False, alpha=alpha)

        # DTD: originally did this outside the function, but have to do it here
        # because of some weird problem
        # FIXME
        # if(ovrdis):
        #     fields_D = fields_D_dict[fieldname]
        #     for dname, dloc, field_D in zip(dis_name_list, dxy_list, fields_D):
        #         Dx = dloc[0] + xoffset
        #         Dy = dloc[1] + yoffset
        #         ax.plot(Dx, Dy, 'k*', ms=8)
        #         ax.annotate(dname, (Dx + 250., Dy + 250.), clip_on=True)
        #         if(field_D is not None):
        #             ax.annotate(disfmtstr.format(field_D), (Dx + 250., Dy - 250.), clip_on=True)

        # Determine plot bounds
        if(plotxmin == -1):
            plotxmin = np.min(xplt)
        if(plotxmax == -1):
            plotxmax = np.max(xplt)
        if(plotymin == -1):
            plotymin = np.min(yplt)
        if(plotymax == -1):
            plotymax = np.max(yplt)

        clvllocator = ticker.MultipleLocator(base=clvls)
        # grid.cbar_axes[0].colorbar(plot1)
        plt.colorbar(plot1, orientation='vertical', ticks=clvllocator, cax=grid.cbar_axes[0])
        # FIXME: something is screwed up with ImageGrid colorbar ticks, but it looks like the
        # correct behavior can be recovered by reverting to the matplotlib version above
#         grid.cbar_axes[0].toggle_label(True)
#         grid.cbar_axes[0].yaxis.set_major_locator(clvllocator)
        # plt.title('dBZ at el = %.2f'%(el/deg2rad)+'and time '+time.strftime(fmt))
        ax.set_xlim(plotxmin, plotxmax)
        ax.set_ylim(plotymin, plotymax)
        formatter = ticker.FuncFormatter(mtokm)
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=10000.0))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(base=10000.0))
        ax.set_aspect('equal')

        figlist.append(fig)
        gridlist.append(grid)

    return figlist, gridlist
#        plt.savefig(image_dir+fieldname+time.strftime(fmt3).strip()+'.png',dpi=200)


def readsweeps2PIPS(fieldnames, pc, ib):
    """Reads in a list of radar sweeps and interpolates fields to disdrometer locations.
       Uses the plotcontrol and text input dictionary bunches. Returns a dictionary containing
       arrays of sweeps and values of the radar fields at the disdrometer locations. Optionally
       also saves or loads the dictionary to or from an npz file."""

    sweepdict = {}

    raddate_path = os.path.join(pc.radar_save_dir, ib.radar_name + '_'
                                + ib.starttimerad.strftime(tm.timefmt3).strip()
                                + '_' + ib.stoptimerad.strftime(tm.timefmt3).strip()
                                + '.txt')
    radnpz_filename = (ib.radar_name + '_' + ib.starttimerad.strftime(tm.timefmt3).strip()
                       + '_' + ib.stoptimerad.strftime(tm.timefmt3).strip()
                       + str(ib.el_req) + '.npz')
    radnpz_path = os.path.join(pc.radar_save_dir, radnpz_filename)

    if(not pc.loadradopt):

        radar_filelist, radtimes = getradarfilelist(ib.radar_dir, pc.radar_save_dir,
                                                    starttime=ib.starttimerad,
                                                    stoptime=ib.stoptimerad, platform=ib.platform,
                                                    el_req=ib.el_req)

        print("radar_filelist = ", radar_filelist)

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

        newradtimes = []
        newradarfilelist = []
        outfieldnames = []

        for index, path, sweeptime in zip(
                range(len(radar_filelist)), radar_filelist, radtimes):
            print("Processing file: " + path)

            if(ib.platform == 'NEXRAD'):
                (outfieldnames, fields, range_start, radar_range, azimuth_start_rad,
                    azimuth_rad, rlat_rad, rlon_rad, ralt, el_rad) = \
                    readCFRadial(True, ib.el_req, ib.rlat,
                                 ib.rlon, ib.ralt, path, sweeptime, fieldnames)
            elif(ib.platform == 'SMARTR'):
                (outfieldnames, fields, range_start, radar_range, azimuth_start_rad,
                    azimuth_rad, rlat_rad, rlon_rad, ralt, el_rad) = \
                    readCFRadial(False, ib.el_req, ib.rlat,
                                 ib.rlon, ib.ralt, path, sweeptime, fieldnames)
            elif(ib.platform == 'UMXP'):
                (outfieldnames, fields, range_start, radar_range, azimuth_start_rad,
                    azimuth_rad, rlat_rad, rlon_rad, ralt, el_rad) = \
                    readUMXPnc(path, sweeptime, fieldnames, heading=ib.heading)

            # Some NEXRAD files appear to only have a couple fields in them. This creates problems
            # with casting to arrays below. For now, we'll just
            # throw out any sweeps that don't have at least 4 ('dBZ', 'ZDR', 'RHV', 'VR') but
            # will need to come up with a better solution in the future.

            if(len(outfieldnames) >= len(fieldnames)):
                fields_tlist.append(fields)
                range_start_tlist.append(range_start)
                range_tlist.append(radar_range)
                azimuth_start_tlist.append(azimuth_start_rad)
                azimuth_tlist.append(azimuth_rad)
                rlat_tlist.append(rlat_rad)
                rlon_tlist.append(rlon_rad)
                ralt_tlist.append(ralt)
                el_tlist.append(el_rad)

                dxy_list, fields_D = dis.rad2DD2(fields, range_start, radar_range,
                                                 azimuth_start_rad, azimuth_rad, rlat_rad,
                                                 rlon_rad, ralt, el_rad, ib.dlocs,
                                                 average_gates=False, Cressman=True,
                                                 roi=1000.)
                print("fields_D.shape", fields_D, fields_D.shape)

                fields_D_tlist.append(fields_D)
                dxy_tlist.append(dxy_list)
                newradtimes.append(sweeptime)
                newradarfilelist.append(path)
            else:
                print("This sweep only has ", len(outfieldnames), " fields. Throwing out!")

        sweepdict['radtimes'] = newradtimes
        sweepdict['fields_tarr'] = np.array(fields_tlist)
        sweepdict['range_start_tarr'] = np.array(range_start_tlist)
        sweepdict['range_tarr'] = np.array(range_tlist)
        sweepdict['azimuth_start_tarr'] = np.array(azimuth_start_tlist)
        sweepdict['azimuth_tarr'] = np.array(azimuth_tlist)
        sweepdict['rlat_tarr'] = np.array(rlat_tlist)
        sweepdict['rlon_tarr'] = np.array(rlon_tlist)
        sweepdict['ralt_tarr'] = np.array(ralt_tlist)
        sweepdict['el_tarr'] = np.array(el_tlist)
        sweepdict['fields_D_tarr'] = np.array(fields_D_tlist)
        sweepdict['dxy_tarr'] = np.array(dxy_tlist)
        sweepdict['outfieldnames'] = outfieldnames
        sweepdict['radar_filelist'] = newradarfilelist

        if(pc.saveradopt):
            if (not os.path.exists(pc.radar_save_dir)):
                os.makedirs(pc.radar_save_dir)

            raddate_file = open(raddate_path, 'wb')
            pickle.dump(radtimes, raddate_file)
            pickle.dump(radar_filelist, raddate_file)
            pickle.dump(outfieldnames, raddate_file)
            raddate_file.close()

            np.savez_compressed(radnpz_path, **sweepdict)

    if(pc.loadradopt):
        raddate_file = open(raddate_path, 'rb')
        radtimes = pickle.load(raddate_file)
        radar_filelist = pickle.load(raddate_file)
        outfieldnames = pickle.load(raddate_file)
        raddate_file.close()

        radnpz_file = np.load(radnpz_path)

        sweepdict = dict(radnpz_file)

        sweepdict['radtimes'] = radtimes
        sweepdict['radar_filelist'] = radar_filelist
        sweepdict['outfieldnames'] = outfieldnames

#         fields_tarr = radnpz_file['fields_tarr']
#         range_start_tarr = radnpz_file['range_start_tarr']
#         range_tarr = radnpz_file['range_tarr']
#         azimuth_start_tarr = radnpz_file['azimuth_start_tarr']
#         azimuth_tarr = radnpz_file['azimuth_tarr']
#         rlat_tarr = radnpz_file['rlat_tarr']
#         rlon_tarr = radnpz_file['rlon_tarr']
#         ralt_tarr = radnpz_file['ralt_tarr']
#         el_tarr = radnpz_file['el_tarr']
#         fields_D_tarr = radnpz_file['fields_D_tarr']
#         dxy_tarr = radnpz_file['dxy_tarr']

    return utils.Bunch(sweepdict)


def plotsweeps(pc, ib, sb, sweepstart=-1, sweepstop=-1):
    """Plots radar sweeps with overlaid disdrometer locations"""

    radar_image_dir = os.path.join(ib.image_dir, 'radar_PPI/')
    if (not os.path.exists(radar_image_dir)):
        os.makedirs(radar_image_dir)

    if(sweepstart == -1):
        sweepstart = sb.radtimes[0]
    if(sweepstop == -1):
        sweepstop = sb.radtimes[-1]

    print("Plotting radar sweeps with overlaid disdrometer locations and data.")
    for index, path, sweeptime in zip(range(len(sb.radar_filelist)), sb.radar_filelist,
                                      sb.radtimes):
        if(sweepstart - timedelta(hours=1) <= sweeptime <= sweepstop + timedelta(hours=1)):

            fields_arr = sb.fields_tarr[index]
            fields_D_arr = sb.fields_D_tarr[index]
            rlat_rad = sb.rlat_tarr[index]
            rlon_rad = sb.rlon_tarr[index]
            ralt = sb.ralt_tarr[index]
            el_rad = sb.el_tarr[index]
            range_start = sb.range_start_tarr[index]
            radar_range = sb.range_tarr[index]
            azimuth_start_rad = sb.azimuth_start_tarr[index]
            azimuth_rad = sb.azimuth_tarr[index]
            dxy_arr = sb.dxy_tarr[index]

            # list(array) where array is 2D or higher yields a list of arrays!
            fields_list = list(fields_arr)
            dxy_list = dxy_arr.tolist()
            # array.tolist() yields a nested list for a 2D+ array!
            # Note, need to transpose here because plotsweep expects the first axis to be the
            # field (i.e. dBZ, ZDR, etc.) and the second axis to be the disdrometer.  What a
            # tangled web I weave!
            fields_D_list = fields_D_arr.T.tolist()

            # In most of our cases the radar location isn't going to change with time,
            # but in the more general case, this may not be true (i.e. if we are dealing with
            # a mobile radar).
            print("Radar sweep time: ", sweeptime.strftime(tm.timefmt))
            # Prepare masks for fields by reflectivity < some threshold

            # Probably should do this using a dictionary
            for field, fieldname in zip(fields_list, sb.outfieldnames):
                if(fieldname == 'dBZ'):
                    mask = np.where(field > 5.0, False, True)

            masklist = [mask, mask, mask]

            # Compute height and range of radar beam above disdrometer (assumes lambert
            # conformal like plotsweep for now)
            for dloc, dname, dxy in zip(ib.dlocs, ib.dis_name_list, dxy_list):
                Dx, Dy = dxy
                print(Dx, Dy)
                # FIXME: below calculation of beam height doesn't take into account disdrometer
                # elevation.
                h, r = oban.computeheightrangesingle(Dx, Dy, np.deg2rad(ib.el_req))
                # In most of our cases the radar location isn't going to change with time,
                # but in the more general case, this may not be true (i.e. if we are dealing
                # with a mobile radar).
                print("Disdrometer name,x,y,radar elevation angle,slant range, "
                      "approximate beam height:")
                print(dname, Dx, Dy, el_rad / deg2rad, r, h)
                if(dloc == ib.dlocs[0] and ib.plotlims[0] == -1):
                    ib.plotlims = [Dx - 25000., Dx + 25000., Dy - 30000., Dy + 20000.]

            figlist, gridlist = plotsweep(ib.radlims, ib.plotlims, sb.outfieldnames,
                                          fields_list, masklist, range_start,
                                          radar_range, azimuth_start_rad,
                                          azimuth_rad, rlat_rad, rlon_rad, ralt, el_rad,
                                          False, pc.plot_radar, ib.dis_name_list,
                                          dxy_list, fields_D_list)

            # Save figures
            for fieldname, fig in zip(sb.outfieldnames, figlist):
                plt.figure(fig.number)
                plt.title(sweeptime.strftime(tm.timefmt2).strip())
                plt.savefig(radar_image_dir + fieldname + sweeptime.strftime(tm.timefmt3).strip()
                            + 'el' + str(ib.el_req) + '.png', dpi=200, bbox_inches='tight')
                plt.close(fig)