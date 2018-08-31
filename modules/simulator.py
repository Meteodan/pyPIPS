# simulator.py: a collection of functions related to the Parsivel simulator

import numpy as np
from scipy.stats import gamma, uniform
from . import disdrometer_module as dis
from . import plotmodule as pm
from shapely.geometry import MultiLineString, LineString
from matplotlib.collections import LineCollection
from datetime import datetime, timedelta
import matplotlib.dates as dates
import pandas as pd
import glob
import os
from . import radarmodule as radar
from . import thermolib as thermo
from .datahandler import getDataHandler
import pyart as pyart
import matplotlib.pyplot as plt
from shapely.geometry import MultiLineString, LineString
from matplotlib.collections import LineCollection
from metpy.plots import ctables

def samplegammaDSD(Nt, lamda, alpha, bins=None):
    """Randomly samples a gamma DSD given Nt, lamda, and alpha."""
    scale = 1. / lamda
    shape = alpha + 1.
    s = gamma.rvs(shape, scale=scale, size=int(Nt))
    # If desired, bin up the resulting samples into the given diameter ranges
    if(bins is None):
        return s
    else:
        ND_sample, _ = np.histogram(s, bins)
        return ND_sample / (bins[1:] - bins[:-1])


def create_random_gamma_DSD(Nt, lamda, alpha, Vt, sampling_length, sampling_width, Dl, Dmid, Dr,
                            Dmin=0, Dmax=20.0, sampling_interval=10., remove_margins=False,
                            verbose=False, rhocorrect=False, rho=None, mask_lowest=True):
    """Given Nt, lamda, alpha, create a spatial distribution in a volume"""
    # First, determine the sampling volume. Use the sampling area A multiplied by the
    # depth that the fasted falling particle would cover in the time given by sampling_interval
    Dmax_index = np.searchsorted(Dr, Dmax / 1000.)
    if verbose:
        print "Dmax_index = ", Dmax_index
    Vt = dis.assignfallspeed(Dmid*1000., rhocorrect=rhocorrect, rho=rho)
    Vtmax = Vt[Dmax_index]
    sampling_height = Vtmax * sampling_interval
    sampling_area = sampling_length * sampling_width
    sampling_volume = sampling_area * sampling_height  # Maximum sampling volume
    # Sampling volumes as a function of diameter D
    sampling_volumes_D = calc_sampling_volumes_D(Vt, Dr, Dmax, sampling_interval, sampling_area)

    if verbose:
        print "sampling height = ", sampling_height
        print "sampling volume = ", sampling_volume

    # Next, create a uniform distribution of n=Nt*sampling_volume drops within
    # the volume given by sampling_volume
    n = int(Nt * sampling_volume)
    if verbose:
        print "number concentration = ", Nt
        print "number of particles in sampling volume = ", n
    xpos = uniform.rvs(0., sampling_length, ((n, 1)))
    ypos = uniform.rvs(0., sampling_width, ((n, 1)))
    zpos = uniform.rvs(0., sampling_height, ((n, 1)))

    # Next, determine the sizes of the drops by drawing the diameters randomly from a gamma
    # distribution given by n, lamda, and alpha

    diameters = samplegammaDSD(n, lamda, alpha)
    if verbose:
        print "minimum, maximum diameter in sample = ", diameters.min(), diameters.max()
        print "maximumm allowed diameter = ", Dmax / 1000.
    # Restrict diameters to be less than Dmax
    diameter_mask = diameters <= Dmax / 1000.
    if verbose:
        print "number of particles less than Dmax = ", diameter_mask.sum()
    # Mask the lowest two diameter bins by default (the real Parsivel does this owing to low SNR)
    if mask_lowest:
        low_mask = diameters > dis.max_diameter[1] / 1000.
        if verbose:
            print "number of particles above the lowest two bins = ", low_mask.sum()
        diameter_mask = diameter_mask & low_mask
    diameters = diameters[diameter_mask]
    xpos = xpos[diameter_mask]
    ypos = ypos[diameter_mask]
    zpos = zpos[diameter_mask]

    if verbose:
        print "number of particles within allowable diameter range = ", diameter_mask.sum()
        print "minimum, maximum particle diameter in truncated sample = ", diameters.min(), diameters.max()

    # Now, figure out which drops in the volume won't fall through the sensor
    # area in the given time, and remove them

    velocities = dis.assignfallspeed(diameters * 1000., rhocorrect=rhocorrect, rho=rho)
    depths = velocities * sampling_interval
    keepers = np.where(zpos.squeeze() - depths <= 0.)

    xpos = xpos[keepers]
    ypos = ypos[keepers]
    zpos = zpos[keepers]

    diameters = diameters[keepers]
    velocities = velocities[keepers]

    # Now, figure out which drops are margin fallers and flag them if desired
    # Below also masks out the margins on the ends of the beam, which may not be correct. One can
    # envision a drop falling and partially striking the top of the Parsivel, shearing it apart with
    # only part of the drop falling through the beam. OTOH, this should show up as a small drop
    # falling too fast, just like an "edge" margin faller, though...

    margins_xl = xpos.squeeze() - diameters / 2. < 0.
    margins_xr = xpos.squeeze() + diameters / 2. > sampling_length
    margins_yl = ypos.squeeze() - diameters / 2. < 0.
    margins_yr = ypos.squeeze() + diameters / 2. > sampling_width
    margin_mask = margins_xl | margins_xr | margins_yl | margins_yr

    if verbose:
        print "number of particles that fall through sampling area = ", xpos.size
        print "number of these that are margin fallers = ", margin_mask.sum()

    if remove_margins:
        if verbose:
            print "Removing margin fallers!"
        xpos = xpos[~margin_mask]
        ypos = ypos[~margin_mask]
        zpos = zpos[~margin_mask]
        diameters = diameters[~margin_mask]
        velocities = velocities[~margin_mask]

    # Bin up particles into the Parsivel diameter bins (later will add velocity too; right now all
    # particles are assumed to fall strictly along theoretical/empirical fall speed curve)
    Dedges = np.append(Dl, Dr[-1])
    Dedges = Dedges[:Dmax_index + 2]
    pcount_binned, _ = np.histogram(diameters, Dedges)
    # Compute ND of sample in original Parsivel diameter and velocity bins
    # Gets particle count/unit volume/diameter interval
    ND = calc_ND(pcount_binned, sampling_volumes_D, Dr, Dl, Dmax)

    positions = np.hstack((xpos, ypos, zpos))
    sample_dict = {'positions': positions, 'diameters': diameters, 'velocities': velocities,
                   'ND': ND, 'margin_mask': margin_mask, 'pcount_binned': pcount_binned,
                   'sampling_volumes_D': sampling_volumes_D}
    return sample_dict


def calc_ND(pcount_binned, sampling_volumes_D, Dr, Dl, Dmax):
    """Calculate the number density ND for diameter bins bounded by Dr, Dl
       given the particle counts in each bin (pcount_binned)"""
    Dmax_index = np.searchsorted(Dr, Dmax / 1000.)
    return pcount_binned / (sampling_volumes_D * (Dr[:Dmax_index + 1] - Dl[:Dmax_index + 1]))


def calc_sampling_volumes_D(Vt, Dr, Dmax, sampling_interval, sampling_area):
    """Calculate the sampling volumes as a function of terminal velocity."""
    Dmax_index = np.searchsorted(Dr, Dmax / 1000.)
    return Vt[:Dmax_index + 1] * sampling_interval * sampling_area

def uniquetol(a, tol=1.e-3):
    """Returns an array with duplicates removed within a given tolerance.
       This one assumes the array is already sorted in either
       increasing or decreasing order."""
    # from https://stackoverflow.com/questions/5426908/
    # find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
    b = np.array(a)
    d = np.append(True, np.abs(np.diff(b)))
    return b[d > tol]

def find_transect_grid_intersections(casedate, grid_dict, dis_dict, model_dict, radar_dict,
                                     vardict, plot_locations=True,  debug=False):
    """Find intersections of probe transects with the model grid"""


    dx = grid_dict['dx']
    dy = grid_dict['dy']
    xe1d = grid_dict['xe1d']
    ye1d = grid_dict['ye1d']
    xcplot = grid_dict['xcplot']
    ycplot = grid_dict['ycplot']
    xeplot = grid_dict['xeplot']
    yeplot = grid_dict['yeplot']
    xcorplot = grid_dict['xcorplot']
    ycorplot = grid_dict['ycorplot']

    lines=[]
    for x in xe1d:
        lines.append(((x, ye1d[0]), (x, ye1d[-1])))

    for y in ye1d:
        lines.append(((xe1d[0], y), (xe1d[-1], y)))

    grid = MultiLineString(lines)

    dh = model_dict[casedate]['DataHandler']
    model_times = model_dict[casedate]['model_times']
    reference_time_sec = model_dict[casedate]['modeltimesec_ref']
    reference_time = model_dict[casedate]['modeltime_ref']
    umove, vmove = radar_dict[casedate]['feature_motion']
    # Grab locations of disdrometers relative to radars and use the
    # last disdrometer in the list as the "reference"
    dxlist = [i[0] for i in dis_dict[casedate]['dradloclist']]
    dylist = [i[1] for i in dis_dict[casedate]['dradloclist']]
    dis_x_rad = dxlist[-1]
    dis_y_rad = dylist[-1]

    if model_dict[casedate]['composite']:

        dis_x, dis_y = model_dict[casedate]['ref_coords_comp']
        model_times_rel = [0.]
        fixed_time = True
    else:
        dis_x, dis_y = model_dict[casedate]['ref_coords']
        model_times_rel = np.array(model_times - reference_time_sec)
        if model_times_rel.size == 1:
            fixed_time = True
        else:
            fixed_time = False

    # Remove element from single-element list
    if fixed_time:
        try:
            vardict = vardict[0]
        except:
            pass

    # Below, we figure out where the disdrometers are relative to the model grid
    xshift = dis_x - dis_x_rad
    yshift = dis_y - dis_y_rad
    dxmodlist = [d_x + xshift for d_x in dxlist]
    dymodlist = [d_y + yshift for d_y in dylist]

    PSDtimes = dis_dict[casedate]['timeseries']['PSDtimes']
#     print PSDtimes
#     print reference_time
    # Loop through each disdrometer

    dis_ts_xylocs = []
    dis_ts_xyslocs = []
    dis_ts_ijlocs = []
    dis_ts_tlocs = []
    dis_ts_vars = []
    dis_ts_vars_points = []
    dis_ts_times = []
    dis_ts_stimes = []

    for dis in xrange(len(dxlist)):
        # Calculate the range of times and the x and y locations at the end of each sampling window
        # of the virtual disdrometer
        # sampling_times contains the times relative to the reference time using the actual
        # disdrometer sampling times

        sampling_times = np.array([(PSDtime - reference_time).total_seconds() for PSDtime in PSDtimes[dis]])
        if fixed_time:
            combined_times = sampling_times
            sample_xlocs = dxmodlist[dis] - umove*sampling_times
            sample_ylocs = dymodlist[dis] - vmove*sampling_times
            dis_xlocs = sample_xlocs
            dis_ylocs = sample_ylocs
        else:
            # Combine the sampling times and model times into a single array
            combined_times = np.concatenate((sampling_times, model_times_rel))
            combined_times = np.unique(combined_times)
            sample_xlocs = dxmodlist[dis] - umove*sampling_times
            sample_ylocs = dymodlist[dis] - vmove*sampling_times
            dis_xlocs = dxmodlist[dis] - umove*combined_times
            dis_ylocs = dymodlist[dis] - vmove*combined_times
#
#     # Sampling interval of the virtual disdrometer
#     sampling_interval = 60.
#     # Seconds before reference time to start the sampling
#     sampling_start = -600.
#     # Seconds after the reference time to stop the sampling
#     sampling_stop = 1500.

        if debug:
            print "sampling_times = ", sampling_times
            print "combined_times = ", combined_times
            print "dis_xlocs = ", dis_xlocs
            print "dis_ylocs = ", dis_ylocs

        line = LineString(np.c_[dis_xlocs, dis_ylocs])

        # fig = plt.figure(figsize=(8,2))
        # ax = fig.add_subplot(111)

        all_times = []
        xlocs = []
        ylocs = []
        ilocs = []
        jlocs = []

        for i, segment in enumerate(line.difference(grid)):
            x, y = segment.xy
        #     print "i = ", i
        #     print "x = ", np.array(x)
        #     print "y = ", np.array(y)
        #     plt.plot(x, y)
        #     plt.text(np.mean(x), np.mean(y), str(i))

            # compute times of each crossing (relative to reference time)
            t = (dxmodlist[dis] - np.array(x))/umove
        #     print "t = ", t
            xlocs = xlocs + list(x)
            ylocs = ylocs + list(y)
            all_times = all_times + list(t)

        # lc = LineCollection(lines, color="gray", lw=1, alpha=0.5)
        # ax.add_collection(lc)
        # ax.set_aspect('equal')

        if debug:
            print "xlocs = ", xlocs
            print "ylocs = ", ylocs
            print "all_times = ", all_times

        # In some cases, just using pandas unique can lead to problems when one list has two values that are
        # very close but not quite the same, and the other has precisely the same. I encountered
        # this problem with P2 for the June 5th case. Need to figure out how to handle it. Does
        # pd.unique have a tolerance parameter?
        # Answer, use this approach from https://stackoverflow.com/questions/5426908/
        # find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
#         xlocs = pd.unique(xlocs)
#         ylocs = pd.unique(ylocs)
#         all_times = pd.unique(all_times)
        xlocs = uniquetol(xlocs)
        ylocs = uniquetol(ylocs)
        all_times = uniquetol(all_times)
        # This creates another problem. Sometimes the above operation removes the wrong "duplicates".
        # We need it to remove those values that are not the ones that coincide with the ones in
        # sampling_times. So we'll go through the all_times array and set any times that are very
        # close to times in sampling_times equal to those corresponding values in sampling_times
        # FIXME: maybe try to improve efficiency of algorithm? (probably not worth it)
        for sidx, stime in np.ndenumerate(sampling_times):
            for aidx, atime in np.ndenumerate(all_times):
                if np.abs(atime-stime) <= 1.e-3:
                    all_times[aidx] = stime

        if debug:
            print "xlocs = ", xlocs
            print "ylocs = ", ylocs
            print "all_times = ", all_times

        # Calculate the fractional i and j indices relative to grid edges corresponding to the x, y locations
        ilocs = (xlocs-xeplot[0,0])/dx
        jlocs = (ylocs-yeplot[0,0])/dy
        if debug:
            print "ilocs = ", ilocs
            print "jlocs = ", jlocs
        # Calculate the fractional time indices corresponding to each point
        if not fixed_time:
            tlocs = (all_times-model_times_rel[0])/model_dt
            tlocs = np.where(tlocs < 0.0, 0.0, tlocs)
        else:
            tlocs = np.array(len(all_times) * [0.0])
        if debug:
            print "tlocs = ", tlocs
        # Find the indices of the edges of the grid boxes that the disdrometer traverses during the sampling time
        igllocs = ilocs.astype(int)
        jgllocs = jlocs.astype(int)
        if debug:
            print "igllocs = ", igllocs
            print "jgllocs = ", jgllocs
        # The -1 below is because we want to identify the right edges with the next lowest grid centers
        igrlocs = np.ceil(ilocs).astype(int)-1
        jgrlocs = np.ceil(jlocs).astype(int)-1
        if debug:
            print "igrlocs = ", igrlocs
            print "jgrlocs = ", jgrlocs
    #     if not fixed_time:
        tgblocs = tlocs.astype(int)
        if debug:
            print "tgblocs = ", tgblocs

        # The grid edge indices we want to use for calculating the sampling depend on the grid motion
        if umove > 0 and vmove > 0:
            locs = [jgrlocs, igrlocs]
        elif umove < 0 and vmove > 0:
            locs = [jgrlocs, igllocs]
        elif umove < 0 and vmove < 0:
            locs = [jgllocs, igllocs]
        else:
            locs = [jgllocs, igrlocs]

        # For plotting we want the indices of the west and south edges. EDIT: Or do we? I don't
        # think so...This way shows the actual grid boxes that are used in the sampling.
        plocs = locs # [jgllocs, igllocs]
        if debug:
            print zip(plocs[0], plocs[1])
            print zip(xcplot[plocs], ycplot[plocs])

        vars_at_ts = []
        time = None
        # TODO: Fix convoluted logic below. Just doing it this way right now for backwards-compatibility
        # with original notebook cell
        for t, tloc in enumerate(tgblocs):
            if model_dict[casedate]['composite']:
                vars_at_ts.append(vardict)
            elif fixed_time:
                if t == 0:
                    time = modeltimeref_sec
                    vars_at_ts.append(vardict)
                else:
                    vars_at_ts.append(vars_at_ts[-1])
            else:
                if(model_times[tloc] != time):
                    time = model_times[tloc]
                    vars_at_ts.append(vardict[tloc])
                else:
                    vars_at_ts.append(vars_at_ts[-1])

            # qrplot = dp_data_2D['qr'].T*1000.

            if plot_locations and dis == 0:
                Zmodplot = vars_at_ts[tloc]['DBZ']
                fig = None
                ax = None
                ptype = 2
                xlim = [xlocs.min() - 5000., xlocs.max() + 10000.]
                ylim = [ylocs.min() - 10000., ylocs.max() + 10000.]
                clevels = np.arange(0., 85., 5.)
                norm, cmap = ctables.registry.get_with_steps('NWSReflectivity', 0., 5.)
                clabel = 'Z (dBZ)'
                cformat = None
                ovrmap = False
                gis_info = None
                numovr = 0
                axesticks = [2000., 2000.]

                fig, ax = pm.plotsingle(fig, ax, ptype, xcplot, ycplot, xcorplot, ycorplot, xlim, ylim, Zmodplot, clevels, cmap, norm,
                                    clevels, clabel, cformat, ovrmap, gis_info, numovr, None, None, None, None, None,
                                    axesticks)

                # Plot the locations of the disdrometer at each time in the sampling period
                #ax.plot(dis_xlocs, dis_ylocs, 'ko', ms=2)
                print "x, y = ", xlocs[t], ylocs[t]
                ax.plot(xlocs[t], ylocs[t], 'ko', ms=4)
                ax.plot(dxmodlist[dis], dymodlist[dis], 'bo', ms=4)
                ax.plot(xcplot[plocs], ycplot[plocs], 'kx', ms=5, alpha=0.5)

                if(all_times[t] == model_times_rel[tloc]):
                    titlestring = 'Time = %06d' % all_times[t] + ' (model time)'
                else:
                    titlestring = 'Time = %06d' % all_times[t] + ' (intermediate)'

                ax.set_title(titlestring)

        # For each variable in the dictionary for each transect time, find the value of that
        # variable at the nearest model grid point to the transect point
        vars_at_ts_points = []
        for t, vard in enumerate(vars_at_ts):
            vard_points = {}
            loc = (locs[0][t], locs[1][t])
            for key, var in vard.iteritems():
                vard_points[key] = var[loc]
            vars_at_ts_points.append(vard_points)

        dis_ts_xylocs.append(zip(xlocs, ylocs))
        dis_ts_xyslocs.append(zip(sample_xlocs, sample_ylocs))
        dis_ts_ijlocs.append(locs)
        dis_ts_times.append(all_times)
        dis_ts_stimes.append(sampling_times)
        dis_ts_tlocs.append(tgblocs)
        dis_ts_vars.append(vars_at_ts)
        dis_ts_vars_points.append(vars_at_ts_points)

    dis_ts_model_dict = {'dis_ts_xylocs': dis_ts_xylocs, 'dis_ts_ijlocs': dis_ts_ijlocs,
                         'dis_ts_times': dis_ts_times, 'dis_ts_stimes': dis_ts_stimes,
                         'dis_ts_tlocs': dis_ts_tlocs, 'dis_ts_vars': dis_ts_vars,
                         'dis_ts_vars_points': dis_ts_vars_points,
                         'dis_ts_xyslocs': dis_ts_xyslocs}
    return dis_ts_model_dict

def init_composite(compositedict, grid_dict):
    # Extract stuff from dictionaries
    xc1d = grid_dict['xc1d']
    yc1d = grid_dict['yc1d']
    xe1d = grid_dict['xe1d']
    ye1d = grid_dict['ye1d']
    zc1d = grid_dict['zc1d']
    ze1d = grid_dict['ze1d']
    dx = grid_dict['dx']
    dy = grid_dict['dy']
    compositewidthx, compositewidthy = compositedict['compositewidth']
    searchboxwidthx, searchboxwidthy = compositedict['searchboxwidth']
    tracking_level = compositedict['tracking_level']
    tracking_varname = compositedict['tracking_varname']
    tracking_extremum = compositedict['tracking_extremum']
    gridlims = compositedict['gridlims']

    igbgn = np.searchsorted(xc1d, gridlims[0], side='left')
    igend = np.searchsorted(xc1d, gridlims[1], side='right')
    jgbgn = np.searchsorted(yc1d, gridlims[2], side='left')
    jgend = np.searchsorted(yc1d, gridlims[3], side='right')

    gridlimindices = [igbgn, igend+1, jgbgn, jgend+1]

    xe1dg = xe1d[igbgn:igend]
    ye1dg = ye1d[jgbgn:jgend]
    xc1dg = xc1d[igbgn:igend]
    yc1dg = yc1d[jgbgn:jgend]

    xckm = xc1dg/1000.
    yckm = yc1dg/1000.

    nxg = igend - igbgn + 1
    nyg = jgend - jgbgn + 1

    print """Starting and ending grid coordinates are:
             igbgn = {:d}, igend = {:d}, jgbgn = {:d}, jgend = {:d}""".format(igbgn, igend, jgbgn, jgend)

    zeagl1d = ze1d[:] - ze1d[0]
    zcagl1d = zc1d[:] - zc1d[0]

    print zeagl1d, zcagl1d

    # Halfwidths of composite box in grid points
    ichw = int(compositewidthx/(2.*dx))
    jchw = int(compositewidthy/(2.*dy))

    xckm_comp = np.arange(-ichw, ichw+1) * dx / 1000.
    yckm_comp = np.arange(-jchw, jchw+1) * dy / 1000.

    # Halfwidths of search box in grid points
    ishw = int(searchboxwidthx/(2.*dx))
    jshw = int(searchboxwidthy/(2.*dy))

    if(tracking_varname == 'w'):
        tracking_klvl = np.where(zeagl1d >= tracking_level)[0][0]
    else:
        tracking_klvl = np.where(zcagl1d >= tracking_level)[0][0]

    print "Tracking "+tracking_extremum+" "+ \
          tracking_varname+" at height {:.1f} m".format(tracking_level)+" (k = {:d})".format(tracking_klvl)

    # Pack new stuff into compositedict
    compositedict['gridlimindices'] = gridlimindices
    compositedict['gcoords'] = (xe1dg, ye1dg, xc1dg, yc1dg, xckm, yckm, zeagl1d, zcagl1d)
    compositedict['ccoords'] = (xckm_comp, yckm_comp)
    compositedict['compositehw'] = (ichw, jchw)
    compositedict['searchhw'] = (ishw, jshw)
    compositedict['tracking_klvl'] = tracking_klvl

    return compositedict


def read_probe_time_series(casedate, dis_dict, radar_dict):
    """Reads in needed timeseries data for the probes for a given case date. May replace
       read_convdata_at_sweeptimes since it is more general."""

    convtimeslist = []
    PSDtimeslist = []
    windspdlist = []
    windspdavgveclist = []
    winddirabslist = []
    winddiravgveclist = []
    templist = []
    dewpointlist = []
    pressurelist = []
    NDlist = []
    intensitylist = []
    pcountlist = []

    # Extract stuff from dictionary
    dis_types = dis_dict[casedate]['dis_types']
    disfilenames = dis_dict[casedate]['disfilenames']
    convfilenames = dis_dict[casedate]['convfilenames']
    starttimes = dis_dict[casedate]['starttimes']
    stoptimes = dis_dict[casedate]['stoptimes']
    dis_dir = dis_dict[casedate]['dis_dir']

    for dis_type, disfilename, convfilename, starttime, stoptime in zip(dis_types, disfilenames, convfilenames,
                                                                        starttimes, stoptimes):

        dis_filepath = os.path.join(dis_dir, disfilename)
        conv_filepath = os.path.join(dis_dir, convfilename)

        if dis_type in 'CU':

            PSD_dict = dis.readCU(conv_filepath, dis_filepath, requested_interval=60.0,
                                   starttime=starttime, stoptime=stoptime)
            windspdstr = 'bwindspd'
            winddirstr = 'bwinddirabs'
        elif dis_type in 'NV2':
            PSD_dict = dis.readNV2netCDF(conv_filepath, dis_filepath, requested_interval=60.0,
                                          starttime=starttime, stoptime=stoptime)
            windspdstr = 'swindspd'
            winddirstr = 'swinddirabs'

        # Extract conventional data timeseries
        conv_df = PSD_dict['conv_df']
        windspds = conv_df[windspdstr].values
        winddirabss = conv_df[winddirstr].values
        temps = conv_df['slowtemp'].values
        dewpoints = conv_df['dewpoint'].values
        pressures = conv_df['pressure'].values

        # Extract PSD data timeseries
        ND = PSD_dict['ND']
        PSD_df = PSD_dict['PSD_df']
        intensity = PSD_df['intensity'].values
        pcount = PSD_df['pcount'].values

        # Plot wind meteogram
        windavgintv = 60
        windgustintv = 3

        # Compute wind speed and direction, and wind gusts
        if dis_type in 'NV2':
            windspdsavg = windspds
            windspdsavgvec = windspds
            windgustsavg = conv_df['swindgust'].values
            winddirsavgvec = winddirabss
        elif dis_type in 'CU':
            windspdsavg,windspdsavgvec,winddirsavgvec,windgusts,windgustsavg = dis.avgwind(winddirabss,
                windspds,windavgintv,gusts=True,gustintv=windgustintv,center=False)
    #     offset = 0
    #     windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg,usavg,vsavg,unit_usavg,unit_vsavg = \
    #     dis.resamplewind(datetimesUTC,offset,winddirabss,windspds,'60S',gusts=True,gustintvstr='3S',center=False)

        convtimes = PSD_dict['convtimestamps']
        PSDtimes = PSD_dict['PSDtimestamps']

        convtimeslist.append(convtimes)
        PSDtimeslist.append(PSDtimes)

        windspdlist.append(windspds)
        windspdavgveclist.append(windspdsavgvec)
        winddirabslist.append(winddirabss)
        winddiravgveclist.append(winddirsavgvec)
        templist.append(temps)
        dewpointlist.append(dewpoints)
        pressurelist.append(pressures)
        NDlist.append(ND)
        intensitylist.append(intensity)
        pcountlist.append(pcount)

    # Stuff all the data into the dis_dict dictionary
    dis_dict[casedate]['timeseries'] = {'windspd': windspdlist, 'windspdavgvec': windspdavgveclist,
                                                    'winddirabs': winddirabslist, 'winddiravgvec': winddiravgveclist,
                                                    'temp': templist, 'dewpoint': dewpointlist,
                                                    'pressure': pressurelist, 'ND': NDlist,
                                                    'intensity': intensitylist, 'pcount': pcountlist,
                                                    'convtimes': convtimeslist, 'PSDtimes': PSDtimeslist}

    return dis_dict


def read_convdata_at_sweeptimes(casedate, dis_dict, radar_dict):
    """Reads in the conventional data valid at the radar times from the probes for the given case
       and stuffs it into the dis_dict dictionary"""
    deployedlist = []
    windspdlist = []
    windspdavgveclist = []
    winddirabslist = []
    winddiravgveclist = []
    templist = []
    dewpointlist = []
    pressurelist = []
    NDlist = []

    # Extract stuff from dictionary
    dis_types = dis_dict[casedate]['dis_types']
    disfilenames = dis_dict[casedate]['disfilenames']
    convfilenames = dis_dict[casedate]['convfilenames']
    starttimes = dis_dict[casedate]['starttimes']
    stoptimes = dis_dict[casedate]['stoptimes']
    dis_dir = dis_dict[casedate]['dis_dir']

    for dis_type, disfilename, convfilename, starttime, stoptime in zip(dis_types, disfilenames, convfilenames,
                                                                        starttimes, stoptimes):

        dis_filepath = os.path.join(dis_dir, disfilename)
        conv_filepath = os.path.join(dis_dir, convfilename)

        if dis_type in 'CU':

            DSD_dict = dis.readCU(conv_filepath, dis_filepath, requested_interval=60.0,
                                   starttime=starttime, stoptime=stoptime)
            windspdstr = 'bwindspd'
            winddirstr = 'bwinddirabs'
        elif dis_type in 'NV2':
            DSD_dict = dis.readNV2netCDF(conv_filepath, dis_filepath, requested_interval=60.0,
                                          starttime=starttime, stoptime=stoptime)
            windspdstr = 'swindspd'
            winddirstr = 'swinddirabs'

        # Extract conventional data timeseries
        conv_df = DSD_dict['conv_df']
        windspds = conv_df[windspdstr].values
        winddirabss = conv_df[winddirstr].values
        temps = conv_df['slowtemp'].values
        dewpoints = conv_df['dewpoint'].values
        pressures = conv_df['pressure'].values

        # Plot wind meteogram
        windavgintv = 60
        windgustintv = 3

        # Compute wind speed and direction, and wind gusts
        if dis_type in 'NV2':
            windspdsavg = windspds
            windspdsavgvec = windspds
            windgustsavg = conv_df['swindgust'].values
            winddirsavgvec = winddirabss
        elif dis_type in 'CU':
            windspdsavg,windspdsavgvec,winddirsavgvec,windgusts,windgustsavg = dis.avgwind(winddirabss,
                windspds,windavgintv,gusts=True,gustintv=windgustintv,center=False)
    #     offset = 0
    #     windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg,usavg,vsavg,unit_usavg,unit_vsavg = \
    #     dis.resamplewind(datetimesUTC,offset,winddirabss,windspds,'60S',gusts=True,gustintvstr='3S',center=False)

        datetimesUTC = DSD_dict['convtimestamps']

        deployedtlist = []
        windspdtlist = []
        windspdavgvectlist = []
        winddirabstlist = []
        winddiravgvectlist = []
        temptlist = []
        dewpointtlist = []
        pressuretlist = []

        # Find the times closest to the sweeptimes
        for sweeptime in radar_dict[casedate]['sweeptimelist']:
            try:
                index = next(i for i, t in enumerate(datetimesUTC) if np.abs((t-sweeptime).total_seconds()) <= 10.)
                deployedtlist.append(True)
            except:
                index = None
                deployedtlist.append(False)

            if index is not None:
                windspdtlist.append(windspds[index])
                winddirabstlist.append(winddirabss[index])
                temptlist.append(temps[index])
                dewpointtlist.append(dewpoints[index])
                pressuretlist.append(pressures[index])
                windspdavgvectlist.append(windspdsavgvec[index])
                winddiravgvectlist.append(winddirsavgvec[index])
            else:
                windspdtlist.append(np.nan)
                winddirabstlist.append(np.nan)
                temptlist.append(np.nan)
                dewpointtlist.append(np.nan)
                pressuretlist.append(np.nan)
                windspdavgvectlist.append(np.nan)
                winddiravgvectlist.append(np.nan)

        deployedarr = np.array(deployedtlist)
        windspdarr = np.array(windspdtlist)
        winddirabsarr = np.array(winddirabstlist)
        temparr = np.array(temptlist)
        dewpointarr = np.array(dewpointtlist)
        pressurearr = np.array(pressuretlist)
        windspdavgvecarr = np.array(windspdavgvectlist)
        winddiravgvecarr = np.array(winddiravgvectlist)

        deployedlist.append(deployedarr)
        windspdlist.append(windspdarr)
        windspdavgveclist.append(windspdavgvecarr)
        winddirabslist.append(winddirabsarr)
        winddiravgveclist.append(winddiravgvecarr)
        templist.append(temparr)
        dewpointlist.append(dewpointarr)
        pressurelist.append(pressurearr)

    # Stuff all the data into the dis_dict dictionary
    dis_dict[casedate]['convdata_at_sweeptimes'] = {'windspd': windspdlist, 'windspdavgvec': windspdavgveclist,
                                                    'winddirabs': winddirabslist, 'winddiravgvec': winddiravgveclist,
                                                    'temp': templist, 'dewpoint': dewpointlist,
                                                    'pressure': pressurelist, 'deployed': deployedlist}

    return dis_dict


def read_sweeps(casedate, radar_dict):
    """Reads sweeps from CFRadial files for a case as defined in the radar_dict dictionary.
       Stuffs these sweeps into the dictionary"""
    # Initial latitude, longitude, altitude set to None (will be read from files)
    radlat = None
    radlon = None
    radalt = None
    radardir = radar_dict[casedate]['radardir']
    radpathlist = glob.glob(radardir+'/*nc')
    radstarttime = radar_dict[casedate]['radstarttime']
    radstoptime = radar_dict[casedate]['radstoptime']
    fieldnames = radar_dict[casedate]['fieldnames']
    el_req = radar_dict[casedate]['el_req']

    # Now read in all the sweeps between radstarttime and radstoptime closest to the requested
    # elevation angle
    radstarttimedt = datetime(np.int(radstarttime[:4]), np.int(radstarttime[4:6]), np.int(
            radstarttime[6:8]), np.int(radstarttime[8:10]), np.int(radstarttime[10:12]))
    radstoptimedt = datetime(np.int(radstoptime[:4]), np.int(radstoptime[4:6]), np.int(
            radstoptime[6:8]), np.int(radstoptime[8:10]), np.int(radstoptime[10:12]))

    outfieldnameslist = []
    radarsweeplist = []
    sweeptimelist = []

    for radpath in radpathlist:
        sweeptime = radar._getsweeptime(radpath)

        if(radstarttimedt <= sweeptime and sweeptime <= radstoptimedt):
            outfieldnames, radarsweep = radar.readCFRadial_pyART(False, el_req, radlat, radlon, radalt, radpath,
                                                                 None, fieldnames, compute_kdp=False)
            outfieldnameslist.append(outfieldnames)
            radarsweeplist.append(radarsweep)
            sweeptimelist.append(sweeptime)

    # Stuff the lists into the dictionary
    radar_dict[casedate]['outfieldnameslist'] = outfieldnameslist
    radar_dict[casedate]['radarsweeplist'] = radarsweeplist
    radar_dict[casedate]['sweeptimelist'] = sweeptimelist

    return radar_dict


def compute_storm_motion(radar_dict):
    """Computes storm motion from start and end points of a given feature.
       Adds a new key to the given dictionary with the storm motion vector
       in a tuple."""
    for casedate, case in radar_dict.iteritems():
        deltat = (case['feature_end_time'] - case['feature_start_time']).total_seconds()
        ustorm = (case['feature_end_loc'][0] - case['feature_start_loc'][0]) * 1000. / deltat
        vstorm = (case['feature_end_loc'][1] - case['feature_start_loc'][1]) * 1000. / deltat
        case['feature_motion'] = (ustorm, vstorm)

    return radar_dict


def get_dis_locs_relative_to_radar(casedate, dis_dict, radar_dict):
    """Gets the disdrometer locations relative to the radar."""
    dgeoloclist = []
    dradloclist = []
    # Extract parameters from dictionary
    dis_names = dis_dict[casedate]['dis_names']
    dis_types = dis_dict[casedate]['dis_types']
    disfilenames = dis_dict[casedate]['disfilenames']
    convfilenames = dis_dict[casedate]['convfilenames']
    starttimes = dis_dict[casedate]['starttimes']
    stoptimes = dis_dict[casedate]['stoptimes']
    dis_dir = dis_dict[casedate]['dis_dir']

    # Grab radar location information from the first sweep in radar_dict
    radarsweep = radar_dict[casedate]['radarsweeplist'][0]
    rlat = radarsweep.latitude['data'][0]
    rlon = radarsweep.longitude['data'][0]

    for dis_name, dis_filename, conv_filename, starttime, stoptime, dis_type in \
            zip(dis_names, disfilenames, convfilenames, starttimes, stoptimes, dis_types):
        filepath = os.path.join(dis_dir, dis_filename)
        if(dis_type == 'PIPS'):
            GPS_lats, GPS_lons, GPS_stats, GPS_alts, dgeoloc = dis.readPIPSloc(filepath)
        elif(dis_type == 'CU'):
            convfilepath = os.path.join(dis_dir, conv_filename)
            GPS_lats, GPS_lons, GPS_stats, GPS_alts, dgeoloc = dis.readCUloc(convfilepath,
                                                                          starttime=starttime,
                                                                          stoptime=stoptime)
        elif(dis_type == 'NV2'):
            dgeoloc = dis.readNV2loc(filepath)

        dradx, drady = pyart.core.geographic_to_cartesian_aeqd(dgeoloc[1], dgeoloc[0], rlon, rlat)
        dx = dradx[0]
        dy = drady[0]
        print dx, dy

        dgeoloclist.append(dgeoloc)
        dradloclist.append((dx, dy))

    # Stuff the locations into the dis_dict dictionary
    dis_dict[casedate]['dgeoloclist'] = dgeoloclist
    dis_dict[casedate]['dradloclist'] = dradloclist

    return dis_dict



def set_dh(casedate, model_dict, radar_dict, fixed_time=True, multitime=True):
    """Reads a DataHandler instance for COMMAS simulation output given information in model_dict, using
       casedate as the key. Puts the DataHandler instance into model_dict. Also associates
       the reference time of the model simulation with that of the reference sweep in radar_dict"""
    modelname = 'COMMAS'
    runname = model_dict[casedate]['runname']
    dirname = model_dict[casedate]['dirname']
    model_dt = model_dict[casedate]['model_dt']
    modeltimesec_ref = model_dict[casedate]['modeltimesec_ref']
    if model_dict[casedate]['composite']:
        model_times = model_dict[casedate]['model_times']
    else:
        model_times = [modeltimesec_ref]
    microphys = model_dict[casedate]['microphys']
    dh = getDataHandler(modelname, dirname, model_times, microphys, multitime=multitime)
    dh.setRun(runname, 0)
    dh.loadTimes()
    dh.setTime(modeltimesec_ref)
    modeltime_ref = radar_dict[casedate]['sweeptime_ref']
    model_dict[casedate]['DataHandler'] = dh
    model_dict[casedate]['modeltime_ref'] = modeltime_ref

    return model_dict


def read_model_grid(dh):
    grid_dict = {}
    gridvarnames = ['xc', 'yc', 'zc', 'zc1d', 'xe', 'ye', 'ze', 'ze1d', 'bgmap']
    gridvars = dh.loadGrid()
    for gridvarname, gridvar in zip(gridvarnames, gridvars):
        grid_dict[gridvarname] = gridvar

    grid_dict['dx'] = grid_dict['xc'][0, 0, 1] - grid_dict['xc'][0, 0, 0]
    grid_dict['dy'] = grid_dict['yc'][0, 1, 0] - grid_dict['yc'][0, 0, 0]
    grid_dict['xe1d'] = grid_dict['xe'][0, 0, :]
    grid_dict['ye1d'] = grid_dict['ye'][0, :, 0]
    grid_dict['xc1d'] = grid_dict['xc'][0, 0, :]
    grid_dict['yc1d'] = grid_dict['yc'][0, :, 0]
    grid_dict['xcplot'] = grid_dict['xc'][0, :, :]
    grid_dict['ycplot'] = grid_dict['yc'][0, :, :]
    grid_dict['xeplot'] = grid_dict['xe'][0, :, :]
    grid_dict['yeplot'] = grid_dict['ye'][0, :, :]
    xcorplot, ycorplot = pm.computecorners(grid_dict['xeplot'], grid_dict['yeplot'])
    grid_dict['xcorplot'] = xcorplot
    grid_dict['ycorplot'] = ycorplot
    return grid_dict


def trackfeature(var, searchboxlims=None, guesscoords=None, extremum='max', debug=False):
    """Tracks a feature based on the extremum (max or min) of a given 2D field. Returns the coordinates
       of the extremum relative to the initial shape, as well as the value of the extremum"""
    if(searchboxlims == None):
        ibgn = 0
        iend = var.shape[1]-1
        jbgn = 0
        jend = var.shape[0]-1
    else:
        ibgn = searchboxlims[0]
        iend = searchboxlims[1]
        jbgn = searchboxlims[2]
        jend = searchboxlims[3]

    if debug:
        print ibgn, iend, jbgn, jend, var.shape

        fig = plt.figure()
        ax = fig.add_subplot(111)

    varsearch = var[jbgn:jend+1, ibgn:iend+1]
    print varsearch.shape

    if debug:
        ax.imshow(varsearch, origin='lower')

    if(extremum == 'max'):
        var_extremum = varsearch.max()
        flatindex = np.argmax(varsearch)
    else:
        var_extremum = varsearch.min()
        flatindex = np.argmin(varsearch)

    #print varsearch.shape
    jrel, irel = np.unravel_index(flatindex, varsearch.shape)

    if debug:
        print "irel, jrel = ", irel, jrel
        ax.plot([irel], [jrel], marker='o', color='red')

    if(guesscoords == None):
        iref = irel
        jref = jrel
    else:
        iref = ibgn+irel
        jref = jbgn+jrel

    if debug:
        print "iref, jref = ", iref, jref

    return var_extremum, iref, jref


def get_composite_grid(grid_dict, compositedict):
    composite_grid_dict = {}
    ichw, jchw = compositedict['compositehw']
    dx = grid_dict['dx']
    dy = grid_dict['dy']
    xc_comp = np.arange(-ichw, ichw+1) * dx
    yc_comp = np.arange(-jchw, jchw+1) * dy
    xe_comp = np.arange(-ichw-0.5, ichw+1.5) * dx
    ye_comp = np.arange(-jchw-0.5, jchw+1.5) * dy

    nxm = xc_comp.shape[0]
    nxe = xe_comp.shape[0]
    nym = yc_comp.shape[0]
    nye = ye_comp.shape[0]
    nzm = grid_dict['zc1d'].shape[0]
    nze = grid_dict['ze1d'].shape[0]

    # Recast xc,xe,yc,ye as 3D arrays
    composite_grid_dict['xc'] = np.lib.stride_tricks.as_strided(xc_comp, strides=((0, 0) + xc_comp.strides), shape=((nzm, nym) + xc_comp.shape))
    composite_grid_dict['xe'] = np.lib.stride_tricks.as_strided(xe_comp, strides=((0, 0) + xe_comp.strides), shape=((nzm, nym) + xe_comp.shape))
    composite_grid_dict['yc'] = np.lib.stride_tricks.as_strided(yc_comp, strides=((0, ) + yc_comp.strides + (0, )), shape=((nzm, ) + yc_comp.shape + (nxm, )))
    composite_grid_dict['ye'] = np.lib.stride_tricks.as_strided(ye_comp, strides=((0, ) + ye_comp.strides + (0, )), shape=((nzm, ) + ye_comp.shape + (nxm, )))
    composite_grid_dict['dx'] = dx
    composite_grid_dict['dy'] = dy
    composite_grid_dict['xe1d'] = composite_grid_dict['xe'][0, 0, :]
    composite_grid_dict['ye1d'] = composite_grid_dict['ye'][0, :, 0]
    composite_grid_dict['xc1d'] = composite_grid_dict['xc'][0, 0, :]
    composite_grid_dict['yc1d'] = composite_grid_dict['yc'][0, :, 0]
    composite_grid_dict['xcplot'] = composite_grid_dict['xc'][0, :, :]
    composite_grid_dict['ycplot'] = composite_grid_dict['yc'][0, :, :]
    composite_grid_dict['xeplot'] = composite_grid_dict['xe'][0, :, :]
    composite_grid_dict['yeplot'] = composite_grid_dict['ye'][0, :, :]
    xcorplot, ycorplot = pm.computecorners(composite_grid_dict['xeplot'], composite_grid_dict['yeplot'])
    composite_grid_dict['xcorplot'] = xcorplot
    composite_grid_dict['ycorplot'] = ycorplot

    return composite_grid_dict


def read_vardict(casedate, model_dict, varlists, varlistv, varlist_derived):
    """Reads in variables from the model files listed in model_dict and returns a list of
       dictionaries containing the variables for each time."""
    model_times = model_dict[casedate]['model_times']
    vardictlist = []
    for t, time in enumerate(model_times):
        vardict = dh.loadVars(varlists)
        temp = dh.loadVars(varlistv)
        vardict.update(temp)
        if 'U' in varlistv:
            vardict['UC'] = 0.5*(vardict['U'][:-1, :-1] + vardict['U'][:-1, 1:])
        if 'V' in varlistv:
            vardict['VC'] = 0.5*(vardict['V'][:-1, :-1] + vardict['V'][1:, :-1])
        if 'PTE' in varlist_derived:
            try:
                vardict['PTE'] = thermo.calpte(vardict['P'], vardict['TH'], vardict['QV'])
            except:
                print "Cannot calculate PTE!"
        vardictlist.append(vardict)
    return vardictlist


def build_composite(casedate, model_dict, compositedict, dh, plotlocs=True):
    # Extract stuff from dictionaries
    ichw, jchw = compositedict['compositehw']
    ishw, jshw = compositedict['searchhw']
    tracking_varname = compositedict['tracking_varname']
    tracking_thresh = compositedict['tracking_thresh']
    tracking_extremum = compositedict['tracking_extremum']
    tracking_klvl = compositedict['tracking_klvl']
    tracking_level = compositedict['tracking_level']
    xe1dg, ye1dg, xc1dg, yc1dg, xckm, yckm, zeagl1d, zcagl1d = compositedict['gcoords']
    gridlimindices = compositedict['gridlimindices']
    igbgn = gridlimindices[0]
    igend = gridlimindices[1]-1
    jgbgn = gridlimindices[2]
    jgend = gridlimindices[3]-1

    model_times = model_dict[casedate]['model_times']
    ntimes = model_times.size

    ireflist = []
    jreflist = []
    varlists = ['DBZ', 'TH', 'QV', 'P']
    varlists_derived = ['PTE', 'UC', 'VC']
    varlistv = ['U', 'V']
    varcompdict = {}
    for varname in varlists:
        varcompdict[varname] = np.zeros((ichw*2+1, jchw*2+1), dtype=np.float)
    for varname in varlists_derived:
        varcompdict[varname] = np.zeros((ichw*2+1, jchw*2+1), dtype=np.float)
    for varname in varlistv:
        varcompdict[varname] = np.zeros((ichw*2+2, jchw*2+2), dtype=np.float)

    if plotlocs:
        figcl = plt.figure()
        axcl = figcl.add_subplot(111)
        axclcolorcycle = [plt.cm.plasma(i) for i in np.linspace(0, 1, ntimes)]

    nt = 0
    tflag = True
    for t, time in enumerate(model_times):
        if t == 0:
            searchboxlims = None
            guesscoords = None

        #print "guesscoords = ",guesscoords
        #print "searchboxlims = ",searchboxlims

        print "The model time is {:d} s".format(int(time))
        dh.setTime(time)
        # Right now, only vortz or w for tracking variable
        if tracking_varname == 'w':
            print "Reading variable " + tracking_varname
            trackvardict = dh.loadVars([var])
            trackvar = trackvardict[tracking_varname][tracking_klvl, ...]
            trackvar = trackvar.squeeze()
        elif tracking_varname == 'vortz':
            print "Deriving variable " + tracking_varname
            # TODO: Make this less brittle (maybe create relational dictionaries for variable names
            # for each model [ARPS, CM1, COMMAS], so that we can use the same variable names for each
            # and it automatically gets translated within the DataHandler code)
            trackvardict = dh.loadVars(['WZ'])
            trackvar = trackvardict['WZ'][tracking_klvl, ...]
            trackvar = trackvar.squeeze()

        var_extremum, iref, jref = trackfeature(trackvar, searchboxlims=searchboxlims,
                                              guesscoords=guesscoords, extremum=tracking_extremum)
        guesscoords = [iref, jref]
        searchboxlims = [iref-ishw, iref+ishw+1, jref-jshw, jref+jshw+1]
        print tracking_extremum + ' ' + tracking_varname + ' at height z = {:.1f} m'.format(tracking_level) + '\n' + \
              'is {:.2f}'.format(var_extremum)
        print "The location is x,y = {:.2f},{:.2f} (i,j = {:d},{:d})".format(xckm[iref], yckm[jref], iref, jref)

        ireflist.append(iref)
        jreflist.append(jref)

        xref = xckm[iref]
        yref = yckm[jref]

        if np.abs(var_extremum) >= tracking_thresh: # Restrict what goes into composite above a certain threshold
            if plotlocs:
                axcl.plot([xref] , [yref], marker='o', color=axclcolorcycle[t])
            nt = nt + 1
            # Read in variables to composite (separate grid dimensions for scalars vs. vector wind components)
            compositeboxindices = [iref-ichw+igbgn, iref+ichw+1+igbgn, jref-jchw+jgbgn, jref+jchw+1+jgbgn, 0, 1]
            compositeboxindices_uv = [iref-ichw+igbgn, iref+ichw+1+igbgn+1, jref-jchw+jgbgn, jref+jchw+1+jgbgn+1, 0, 1]

            vardict = dh.loadVars(varlists, compositeboxindices)
            temp = dh.loadVars(varlistv, compositeboxindices_uv)
            vardict.update(temp)

            # Compute velocity components at scalar points
            vardict['UC'] = 0.5*(vardict['U'][:-1, :-1] + vardict['U'][:-1, 1:])
            vardict['VC'] = 0.5*(vardict['V'][:-1, :-1] + vardict['V'][1:, :-1])

            # Compute equivalent potential temperature
            vardict['PTE'] = thermo.calpte(vardict['P'], vardict['TH'], vardict['QV'])

            # Now read in the microphysics scalars
            mp_data, consts = dh.loadMicrophysics()
            # Extract the lowest model level and store in mp_data_2D
            mp_data_2D = {}
            for key, dat in mp_data.iteritems():
                mp_data_2D[key] = dat.T[0, compositeboxindices[2]:compositeboxindices[3],
                                           compositeboxindices[0]:compositeboxindices[1]]
                if tflag:
                    varcompdict[key] = np.zeros((ichw*2+1, jchw*2+1), dtype=np.float)

            print mp_data_2D.keys()
            print varcompdict.keys()

            # Cumulative sum for each variable for each time
            for varname, var in vardict.iteritems():
                varcompdict[varname] = varcompdict[varname] + var
            for varname, var in mp_data_2D.iteritems():
                varcompdict[varname] = varcompdict[varname] + var
            tflag = False
    # Compute composite averages for each variable (nt is the total number of times)
    for varname, var in varcompdict.iteritems():
        varcompdict[varname] = varcompdict[varname]/nt

    if plotlocs:
        axcl.set_aspect('equal')
        figcl.savefig(model_dict[casedate]['runname'] + '_sfcvortloc_{:06d}_{:06d}.png'.format(int(model_times[0]), int(model_times[-1])), dpi=200)

    return varcompdict