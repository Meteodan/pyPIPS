# simulator.py: a collection of functions related to the Parsivel simulator

import inspect
import numpy as np
from scipy.stats import gamma, uniform
# from . import disdrometer_module as dis
from . import PIPS as pips
from . import parsivel_params as pp
from . import plotmodule as pm
from . import DSDlib as dsd
from .legacy import disdrometer_module as dis
from shapely.geometry import MultiLineString, LineString
from datetime import datetime
import glob
import os
from . import radarmodule as radar
from . import thermolib as thermo
from .legacy.datahandler import getDataHandler
from . import dualpara as dualpol
import pyart as pyart
import matplotlib.pyplot as plt
from metpy.plots import ctables
import matplotlib.ticker as ticker
import xarray as xr
from pyCRMtools.modules import utils as CRMutils
from pyCRMtools.pycaps import arps_read
from pyCRMtools.pycaps import pycaps_fields
from joblib import Parallel, delayed

# Some global parameters to make my life easier
rhoacst = 1.0  # kg m^-3
rhorcst = 1000.  # kg m^-3
cr = rhorcst * np.pi / 6.
mur = 1. / 3.  # FIXME: Assume rain is gamma-diameter for now!

sampling_area = pp.parsivel_parameters['sensor_area_mm2']
sampling_width = pp.parsivel_parameters['sensor_width_mm']
sampling_length = pp.parsivel_parameters['sensor_length_mm']

D = pips.parsivel_parameters['avg_diameter_bins_mm']
Dl = pips.parsivel_parameters['min_diameter_bins_mm']
Dr = pips.parsivel_parameters['max_diameter_bins_mm']
Dedges = np.append(Dl, Dr[-1])
bin_width = Dr - Dl


def get_Dmax_index(Dr, Dmax):
    if Dmax is not None:
        Dmax_index = np.searchsorted(Dr, Dmax / 1000., side='right')
    else:
        Dmax_index = np.size(Dr) - 1
        Dmax = Dr[Dmax_index] * 1000.
    return Dmax, Dmax_index


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
                            Dmin=0, Dmax=None, sampling_interval=10., remove_margins=False,
                            verbose=False, rhocorrect=False, rho=None, mask_lowest=True,
                            perturb_vel=True):
    """Given Nt, lamda, alpha, create a spatial distribution in a volume"""
    # First, determine the sampling volume. Use the sampling area A multiplied by the
    # depth that the fasted falling particle would cover in the time given by sampling_interval
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    if verbose:
        print("Dmax_index = ", Dmax_index)
    # Vt = dis.assignfallspeed(Dmid * 1000., rhocorrect=rhocorrect, rho=rho)
    Vtmax = Vt[Dmax_index]
    sampling_height = Vtmax * sampling_interval
    sampling_area = sampling_length * sampling_width
    sampling_volume = sampling_area * sampling_height  # Maximum sampling volume
    # Sampling volumes as a function of diameter D
    sampling_volumes_D = calc_sampling_volumes_D(Vt, Dr, Dmax, sampling_interval, sampling_area)

    if verbose:
        print("sampling height = ", sampling_height)
        print("sampling volume = ", sampling_volume)

    # Next, create a uniform distribution of n=Nt*sampling_volume drops within
    # the volume given by sampling_volume
    n = int(Nt * sampling_volume)
    if verbose:
        print("number concentration = ", Nt)
        print("number of particles in sampling volume = ", n)
    xpos = uniform.rvs(0., sampling_length, ((n, 1)))
    ypos = uniform.rvs(0., sampling_width, ((n, 1)))
    zpos = uniform.rvs(0., sampling_height, ((n, 1)))

    # Next, determine the sizes of the drops by drawing the diameters randomly from a gamma
    # distribution given by n, lamda, and alpha

    diameters = samplegammaDSD(n, lamda, alpha)
    if verbose and diameters.size:
        print("minimum, maximum diameter in sample = ", diameters.min(), diameters.max())
        print("maximum allowed diameter = ", Dmax / 1000.)
    # Restrict diameters to be less than Dmax
    diameter_mask = diameters <= Dmax / 1000.
    if verbose:
        print("number of particles less than Dmax = ", diameter_mask.sum())
    # Mask the lowest two diameter bins by default (the real Parsivel does this owing to low SNR)
    if mask_lowest:
        low_mask = diameters > Dr[1] / 1000.
        if verbose:
            print("number of particles above the lowest two bins = ", low_mask.sum())
        diameter_mask = diameter_mask & low_mask
    diameters = diameters[diameter_mask]
    xpos = xpos[diameter_mask]
    ypos = ypos[diameter_mask]
    zpos = zpos[diameter_mask]

    if verbose:
        print("number of particles within allowable diameter range = ", diameter_mask.sum())
        if diameter_mask.sum() > 0.:
            print(
                "minimum, maximum particle diameter in truncated sample = ",
                diameters.min(),
                diameters.max())

    # Now, figure out which drops in the volume won't fall through the sensor
    # area in the given time, and remove them
    velocities = pips.calc_empirical_fallspeed(diameters * 1000., correct_rho=rhocorrect, rho=rho)
    # TODO: add gaussian perturbations to velocities as an option
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
        print("number of particles that fall through sampling area = ", xpos.size)
        print("number of these that are margin fallers = ", margin_mask.sum())

    if remove_margins:
        if verbose:
            print("Removing margin fallers!")
        xpos = xpos[~margin_mask]
        ypos = ypos[~margin_mask]
        zpos = zpos[~margin_mask]
        diameters = diameters[~margin_mask]
        velocities = velocities[~margin_mask]

    # Bin up particles into the Parsivel diameter bins (later will add velocity too; right now all
    # particles are assumed to fall strictly along theoretical/empirical fall speed curve)
    Dedges = np.append(Dl, Dr[-1])
    print(Dedges.shape)
    # Dedges = Dedges[:Dmax_index + 1]
    pcount_binned, _ = np.histogram(diameters, Dedges)
    # print('diameters shape, dedges shape, pcount_binned shape', diameters.shape, Dedges.shape,
    # pcount_binned.shape)
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
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    return pcount_binned / (sampling_volumes_D * (Dr[:Dmax_index+1] - Dl[:Dmax_index+1]))


def calc_sampling_volumes_D(Vt, Dr, Dmax, sampling_interval, sampling_area):
    """Calculate the sampling volumes as a function of terminal velocity."""
    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)
    return Vt[:Dmax_index+1] * sampling_interval * sampling_area


def uniquetol(a, tol=1.e-3):
    """Returns an array with duplicates removed within a given tolerance.
       This one assumes the array is already sorted in either
       increasing or decreasing order."""
    # from https://stackoverflow.com/questions/5426908/
    # find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
    b = np.array(a)
    d = np.append(True, np.abs(np.diff(b)))
    return b[d > tol]


def combine_sample_and_model_times(model_datetimes, probe_datetimes):
    """Given a range of model datetimes and a range of probe datetimes, combine them
       into a single range"""


def find_transect_grid_intersections(grid_dict, dis_dict, model_dict, radar_dict,
                                     vardict, plot_locations=True, debug=False):
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

    lines = []
    for x in xe1d:
        lines.append(((x, ye1d[0]), (x, ye1d[-1])))

    for y in ye1d:
        lines.append(((xe1d[0], y), (xe1d[-1], y)))

    grid = MultiLineString(lines)

    # dh = model_dict['DataHandler']
    model_times = model_dict['model_times']
    reference_time_sec = model_dict['modeltimesec_ref']
    reference_time = model_dict['modeltime_ref']
    # TODO: Make sure that there is a "model_dt" key in model_dict
    model_dt = model_dict['model_dt']
    umove, vmove = radar_dict['feature_motion']
    # Grab locations of disdrometers relative to radars and use the
    # last disdrometer in the list as the "reference"
    dxlist = [i[0] for i in dis_dict['dradloclist']]
    dylist = [i[1] for i in dis_dict['dradloclist']]
    dis_x_rad = dxlist[-1]
    dis_y_rad = dylist[-1]

    if model_dict['composite']:

        dis_x, dis_y = model_dict['ref_coords_comp']
        model_times_rel = [0.]
        fixed_time = True
    else:
        dis_x, dis_y = model_dict['ref_coords']
        model_times_rel = np.array(model_times) - reference_time_sec
        if model_times_rel.size == 1:
            fixed_time = True
        else:
            fixed_time = False

    # Remove element from single-element list
    if fixed_time:
        try:
            vardict = vardict[0]
        except BaseException:
            pass

    # Below, we figure out where the disdrometers are relative to the model grid
    xshift = dis_x - dis_x_rad
    yshift = dis_y - dis_y_rad
    dxmodlist = [d_x + xshift for d_x in dxlist]
    dymodlist = [d_y + yshift for d_y in dylist]

    PSDtimes = dis_dict['timeseries']['times']
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

    for disnum in range(len(dxlist)):
        # Calculate the range of times and the x and y locations at the end of each sampling window
        # of the virtual disdrometer
        # sampling_times contains the times relative to the reference time using the actual
        # disdrometer sampling times

        sampling_times = np.array([(PSDtime - reference_time).total_seconds()
                                   for PSDtime in PSDtimes[dis]])
        if fixed_time:
            combined_times = sampling_times
            sample_xlocs = dxmodlist[dis] - umove * sampling_times
            sample_ylocs = dymodlist[dis] - vmove * sampling_times
            dis_xlocs = sample_xlocs
            dis_ylocs = sample_ylocs
        else:
            # Combine the sampling times and model times into a single array
            combined_times = np.concatenate((sampling_times, model_times_rel))
            combined_times = np.unique(combined_times)
            sample_xlocs = dxmodlist[dis] - umove * sampling_times
            sample_ylocs = dymodlist[dis] - vmove * sampling_times
            dis_xlocs = dxmodlist[dis] - umove * combined_times
            dis_ylocs = dymodlist[dis] - vmove * combined_times
#
#     # Sampling interval of the virtual disdrometer
#     sampling_interval = 60.
#     # Seconds before reference time to start the sampling
#     sampling_start = -600.
#     # Seconds after the reference time to stop the sampling
#     sampling_stop = 1500.

        if debug:
            print("sampling_times = ", sampling_times)
            print("combined_times = ", combined_times)
            print("dis_xlocs = ", dis_xlocs)
            print("dis_ylocs = ", dis_ylocs)

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
            t = (dxmodlist[dis] - np.array(x)) / umove
        #     print "t = ", t
            xlocs = xlocs + list(x)
            ylocs = ylocs + list(y)
            all_times = all_times + list(t)

        # lc = LineCollection(lines, color="gray", lw=1, alpha=0.5)
        # ax.add_collection(lc)
        # ax.set_aspect('equal')

        if debug:
            print("xlocs = ", xlocs)
            print("ylocs = ", ylocs)
            print("all_times = ", all_times)

        # In some cases, just using pandas unique can lead to problems when one list has two values
        # that are very close but not quite the same, and the other has precisely the same.
        # I encountered this problem with P2 for the June 5th case. Need to figure out how to
        # handle it. Does pd.unique have a tolerance parameter?
        # Answer, use this approach from https://stackoverflow.com/questions/5426908/
        # find-unique-elements-of-floating-point-array-in-numpy-with-comparison-using-a-d
#         xlocs = pd.unique(xlocs)
#         ylocs = pd.unique(ylocs)
#         all_times = pd.unique(all_times)
        xlocs = uniquetol(xlocs)
        ylocs = uniquetol(ylocs)
        all_times = uniquetol(all_times)
        # This creates another problem. Sometimes the above operation removes the wrong
        # "duplicates". We need it to remove those values that are not the ones that coincide with
        # the ones in sampling_times. So we'll go through the all_times array and set any times
        # that are very close to times in sampling_times equal to those corresponding values in
        # sampling_times
        # FIXME: maybe try to improve efficiency of algorithm? (probably not worth it)
        for sidx, stime in np.ndenumerate(sampling_times):
            for aidx, atime in np.ndenumerate(all_times):
                if np.abs(atime - stime) <= 1.e-3:
                    all_times[aidx] = stime

        if debug:
            print("xlocs = ", xlocs)
            print("ylocs = ", ylocs)
            print("all_times = ", all_times)

        # Calculate the fractional i and j indices relative to grid edges
        # corresponding to the x, y locations
        ilocs = (xlocs - xeplot[0, 0]) / dx
        jlocs = (ylocs - yeplot[0, 0]) / dy
        if debug:
            print("ilocs = ", ilocs)
            print("jlocs = ", jlocs)
        # Calculate the fractional time indices corresponding to each point
        if not fixed_time:
            tlocs = (all_times - model_times_rel[0]) / model_dt
            tlocs = np.where(tlocs < 0.0, 0.0, tlocs)
        else:
            tlocs = np.array(len(all_times) * [0.0])
        if debug:
            print("tlocs = ", tlocs)
        # Find the indices of the edges of the grid boxes that the disdrometer
        # traverses during the sampling time
        igllocs = ilocs.astype(int)
        jgllocs = jlocs.astype(int)
        if debug:
            print("igllocs = ", igllocs)
            print("jgllocs = ", jgllocs)
        # The -1 below is because we want to identify the right edges with the
        # next lowest grid centers
        igrlocs = np.ceil(ilocs).astype(int) - 1
        jgrlocs = np.ceil(jlocs).astype(int) - 1
        if debug:
            print("igrlocs = ", igrlocs)
            print("jgrlocs = ", jgrlocs)
    #     if not fixed_time:
        tgblocs = tlocs.astype(int)
        if debug:
            print("tgblocs = ", tgblocs)

        # The grid edge indices we want to use for calculating the sampling depend
        # on the grid motion
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
        plocs = locs  # [jgllocs, igllocs]
        if debug:
            print(list(zip(plocs[0], plocs[1])))
            print(list(zip(xcplot[plocs], ycplot[plocs])))

        vars_at_ts = []
        time = None
        # TODO: Fix convoluted logic below. Just doing it this way right now for
        # backwards-compatibility with original notebook cell
        for t, tloc in enumerate(tgblocs):
            if model_dict['composite']:
                vars_at_ts.append(vardict)
            elif fixed_time:
                if t == 0:
                    time = reference_time_sec
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

                fig, ax = pm.plotsingle(fig, ax, ptype, xcplot, ycplot, xcorplot, ycorplot, xlim,
                                        ylim, Zmodplot, clevels, cmap, norm, clevels, clabel,
                                        cformat, ovrmap, gis_info, numovr, None, None, None, None,
                                        None, axesticks)

                # Plot the locations of the disdrometer at each time in the sampling period
                # ax.plot(dis_xlocs, dis_ylocs, 'ko', ms=2)
                print("x, y = ", xlocs[t], ylocs[t])
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
            for key, var in vard.items():
                vard_points[key] = var[loc]
            vars_at_ts_points.append(vard_points)

        dis_ts_xylocs.append(list(zip(xlocs, ylocs)))
        dis_ts_xyslocs.append(list(zip(sample_xlocs, sample_ylocs)))
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

    gridlimindices = [igbgn, igend + 1, jgbgn, jgend + 1]

    xe1dg = xe1d[igbgn:igend]
    ye1dg = ye1d[jgbgn:jgend]
    xc1dg = xc1d[igbgn:igend]
    yc1dg = yc1d[jgbgn:jgend]

    xckm = xc1dg / 1000.
    yckm = yc1dg / 1000.

    # nxg = igend - igbgn + 1
    # nyg = jgend - jgbgn + 1

    print("""Starting and ending grid coordinates are:
             igbgn = {:d}, igend = {:d}, jgbgn = {:d}, jgend = {:d}""".format(igbgn, igend, jgbgn,
                                                                              jgend))

    zeagl1d = ze1d[:] - ze1d[0]
    zcagl1d = zc1d[:] - zc1d[0]

    print(zeagl1d, zcagl1d)

    # Halfwidths of composite box in grid points
    ichw = int(compositewidthx / (2. * dx))
    jchw = int(compositewidthy / (2. * dy))

    xckm_comp = np.arange(-ichw, ichw + 1) * dx / 1000.
    yckm_comp = np.arange(-jchw, jchw + 1) * dy / 1000.

    # Halfwidths of search box in grid points
    ishw = int(searchboxwidthx / (2. * dx))
    jshw = int(searchboxwidthy / (2. * dy))

    if(tracking_varname == 'w'):
        tracking_klvl = np.where(zeagl1d >= tracking_level)[0][0]
    else:
        tracking_klvl = np.where(zcagl1d >= tracking_level)[0][0]

    print("Tracking " + tracking_extremum + " " +
          tracking_varname + " at height {:.1f} m".format(tracking_level) +
          " (k = {:d})".format(tracking_klvl))

    # Pack new stuff into compositedict
    compositedict['gridlimindices'] = gridlimindices
    compositedict['gcoords'] = (xe1dg, ye1dg, xc1dg, yc1dg, xckm, yckm, zeagl1d, zcagl1d)
    compositedict['ccoords'] = (xckm_comp, yckm_comp)
    compositedict['compositehw'] = (ichw, jchw)
    compositedict['searchhw'] = (ishw, jshw)
    compositedict['tracking_klvl'] = tracking_klvl

    return compositedict


def read_probe_time_series(dis_dict, radar_dict):
    """Reads in needed timeseries data for the probes for a given case date. May replace
       read_convdata_at_sweeptimes since it is more general."""

    timeslist = []
    # windspdlist = []
    windspdavgveclist = []
    # winddirabslist = []
    winddiravgveclist = []
    templist = []
    dewpointlist = []
    pressurelist = []
    rholist = []
    NDlist = []
    intensitylist = []
    pcountlist = []

    # Extract stuff from dictionary
    dis_types = dis_dict['dis_types']
    disfilenames = dis_dict['disfilenames']
    convfilenames = dis_dict['convfilenames']
    starttimes = dis_dict['starttimes']
    stoptimes = dis_dict['stoptimes']
    dis_dir = dis_dict['dis_dir']

    for dis_type, disfilename, convfilename, starttime, stoptime in zip(dis_types, disfilenames,
                                                                        convfilenames, starttimes,
                                                                        stoptimes):

        dis_filepath = os.path.join(dis_dir, disfilename)
        conv_filepath = os.path.join(dis_dir, convfilename)

        if dis_type in 'CU':

            PSD_dict = dis.readCU(conv_filepath, dis_filepath, requested_interval=60.0,
                                  starttime=starttime, stoptime=stoptime)
            # windspdstr = 'bwindspd'
            # winddirstr = 'bwinddirabs'
        elif dis_type in 'NV2':
            PSD_dict = dis.readNV2netCDF(conv_filepath, dis_filepath, requested_interval=60.0,
                                         starttime=starttime, stoptime=stoptime)
            # windspdstr = 'swindspd'
            # winddirstr = 'swinddirabs'

        # Extract conventional data timeseries and resample to PSD times
        conv_df = PSD_dict['conv_df']
        PSDtimes = PSD_dict['PSDtimestamps']
        sec_offset = PSDtimes[0].second
        DSD_interval = PSD_dict['DSD_interval']
        DSD_index = PSD_dict['DSD_index']
        conv_resampled_df = dis.resampleconv(dis_type, DSD_interval, sec_offset, conv_df)
        conv_resampled_df = conv_resampled_df.loc[conv_resampled_df.index.intersection(DSD_index)]

        windspdsavgvec = conv_resampled_df['windspdavgvec'].values
        winddirsavgvec = conv_resampled_df['winddiravgvec'].values
        temps = conv_resampled_df['slowtemp'].values
        dewpoints = conv_resampled_df['dewpoint'].values
        pressures = conv_resampled_df['pressure'].values
        rhos = conv_resampled_df['rho'].values

        # Extract PSD data timeseries
        ND = PSD_dict['ND']
        PSD_df = PSD_dict['PSD_df']
        intensity = PSD_df['intensity'].values
        pcount = PSD_df['pcount'].values
#
#         # Plot wind meteogram
#         windavgintv = 60
#         windgustintv = 3
#
#         # Compute wind speed and direction, and wind gusts
#         if dis_type in 'NV2':
#             windspdsavg = windspds
#             windspdsavgvec = windspds
#             windgustsavg = conv_resampled_df['swindgust'].values
#             winddirsavgvec = winddirabss
#         elif dis_type in 'CU':
#             windspdsavg,windspdsavgvec,winddirsavgvec,windgusts,windgustsavg = \
#               dis.avgwind(winddirabss,
#                 windspds,windavgintv,gusts=True,gustintv=windgustintv,center=False)
    #     offset = 0
    #     windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg, \
    #       usavg,vsavg,unit_usavg,unit_vsavg = \
    #     dis.resamplewind(datetimesUTC,offset,winddirabss,windspds,'60S',gusts=True,
    #                      gustintvstr='3S',center=False)

        timeslist.append(PSDtimes)

        windspdavgveclist.append(windspdsavgvec)
        winddiravgveclist.append(winddirsavgvec)
        templist.append(temps)
        dewpointlist.append(dewpoints)
        pressurelist.append(pressures)
        rholist.append(rhos)
        NDlist.append(ND)
        intensitylist.append(intensity)
        pcountlist.append(pcount)

    # Stuff all the data into the dis_dict dictionary
    dis_dict['timeseries'] = {'windspdavgvec': windspdavgveclist,
                              'winddiravgvec': winddiravgveclist,
                              'temp': templist, 'dewpoint': dewpointlist,
                              'pressure': pressurelist, 'rho': rholist, 'ND': NDlist,
                              'intensity': intensitylist, 'pcount': pcountlist,
                              'times': timeslist}

    return dis_dict


def read_convdata_at_modeltimes(dis_dict, model_dict):
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
    # NDlist = []

    # Extract stuff from dictionary
    dis_types = dis_dict['dis_types']
    disfilenames = dis_dict['disfilenames']
    try:
        convfilenames = dis_dict['convfilenames']
    except KeyError:
        convfilenames = [None] * len(dis_types)
    starttimes = dis_dict['starttimes']
    stoptimes = dis_dict['stoptimes']
    dis_dir = dis_dict['dis_dir']
    interval = dis_dict['interval']

    for dis_type, disfilename, convfilename, starttime, stoptime in zip(dis_types, disfilenames,
                                                                        convfilenames, starttimes,
                                                                        stoptimes):

        dis_filepath = os.path.join(dis_dir, disfilename)
        if convfilename is not None:
            conv_filepath = os.path.join(dis_dir, convfilename)

        if dis_type in 'CU':

            DSD_dict = dis.readCU(conv_filepath, dis_filepath, requested_interval=interval,
                                  starttime=starttime, stoptime=stoptime)
            windspdstr = 'bwindspd'
            winddirstr = 'bwinddirabs'
            tempstr = 'slowtemp'
        elif dis_type in 'NV2':
            DSD_dict = dis.readNV2netCDF(conv_filepath, dis_filepath, requested_interval=interval,
                                         starttime=starttime, stoptime=stoptime)
            windspdstr = 'swindspd'
            winddirstr = 'swinddirabs'
            tempstr = 'slowtemp'
        elif dis_type in 'PIPS':
            DSD_dict = dis.readPIPS(dis_filepath, requested_interval=interval,
                                    starttime=starttime, stoptime=stoptime)
            windspdstr = 'windspd'
            winddirstr = 'winddirabs'
            tempstr = 'fasttemp'

        # Extract conventional data timeseries
        conv_df = DSD_dict['conv_df']
        windspds = conv_df[windspdstr].values
        winddirabss = conv_df[winddirstr].values
        temps = conv_df[tempstr].values
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
            windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg = \
                dis.avgwind(winddirabss, windspds, windavgintv, gusts=True, gustintv=windgustintv,
                            center=False)
        elif dis_type in 'PIPS':
            windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg = dis.avgwind(
                winddirabss, windspds, windavgintv, gusts=True, gustintv=windgustintv, center=False)
    #     offset = 0
    #     windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg, \
    #       usavg,vsavg,unit_usavg,unit_vsavg = \
    #     dis.resamplewind(datetimesUTC,offset,winddirabss,windspds,'60S',gusts=True,
    #                      gustintvstr='3S',center=False)

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
        for modeltime in model_dict['datetime_range']:
            try:
                index = next(
                    i for i, t in enumerate(datetimesUTC) if np.abs(
                        (t - modeltime).total_seconds()) <= 10.)
                deployedtlist.append(True)
            except Exception:
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
    dis_dict['convdata_at_modeltimes'] = {
        'windspd': windspdlist, 'windspdavgvec': windspdavgveclist, 'winddirabs': winddirabslist,
        'winddiravgvec': winddiravgveclist, 'temp': templist, 'dewpoint': dewpointlist,
        'pressure': pressurelist, 'deployed': deployedlist
    }

    return dis_dict


def read_convdata_at_sweeptimes(dis_dict, radar_dict, resample_interval=60.):
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
    # NDlist = []

    # Extract stuff from dictionary
    dis_types = dis_dict['dis_types']
    disfilenames = dis_dict['disfilenames']
    try:
        convfilenames = dis_dict['convfilenames']
    except KeyError:
        convfilenames = [None] * len(dis_types)
    starttimes = dis_dict['starttimes']
    stoptimes = dis_dict['stoptimes']
    dis_dir = dis_dict['dis_dir']
    interval = dis_dict['interval']

    for dis_type, disfilename, convfilename, starttime, stoptime in zip(dis_types, disfilenames,
                                                                        convfilenames, starttimes,
                                                                        stoptimes):

        dis_filepath = os.path.join(dis_dir, disfilename)
        if convfilename is not None:
            conv_filepath = os.path.join(dis_dir, convfilename)

        if dis_type in 'CU':

            DSD_dict = dis.readCU(conv_filepath, dis_filepath, requested_interval=interval,
                                  starttime=starttime, stoptime=stoptime)
            windspdstr = 'bwindspd'
            winddirstr = 'bwinddirabs'
            tempstr = 'slowtemp'
        elif dis_type in 'NV2':
            DSD_dict = dis.readNV2netCDF(conv_filepath, dis_filepath, requested_interval=interval,
                                         starttime=starttime, stoptime=stoptime)
            windspdstr = 'swindspd'
            winddirstr = 'swinddirabs'
            tempstr = 'slowtemp'
        elif dis_type in 'PIPS':
            DSD_dict = dis.readPIPS(dis_filepath, requested_interval=interval,
                                    starttime=starttime, stoptime=stoptime)
            windspdstr = 'windspd'
            winddirstr = 'winddirabs'
            tempstr = 'fasttemp'

        # Extract conventional data timeseries
        conv_df = DSD_dict['conv_df']
        windspds = conv_df[windspdstr].values
        winddirabss = conv_df[winddirstr].values
        temps = conv_df[tempstr].values
        dewpoints = conv_df['dewpoint'].values
        pressures = conv_df['pressure'].values

        sec_offset = conv_df.index.to_pydatetime()[0].second
        conv_df = pips.resample_conv(dis_type, resample_interval, sec_offset, conv_df, gusts=True,
                                     gustintvstr='3S', center=False)

        windspdsavgvec = conv_df['windspdavgvec'].values
        winddirsavgvec = conv_df['winddiravgvec'].values

        # Plot wind meteogram
        # windavgintv = 60
        # windgustintv = 3

    #     # Compute wind speed and direction, and wind gusts
    #     if dis_type in 'NV2':
    #         windspdsavg = windspds
    #         windspdsavgvec = windspds
    #         windgustsavg = conv_df['swindgust'].values
    #         winddirsavgvec = winddirabss
    #     elif dis_type in 'CU':
    #         windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg = \
    #             dis.avgwind(winddirabss, windspds, windavgintv, gusts=True, gustintv=windgustintv,
    #                         center=False)
    #     elif dis_type in 'PIPS':
    #         windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg = dis.avgwind(
    #             winddirabss, windspds, windavgintv, gusts=True, gustintv=windgustintv,
    # center=False)
    # #     offset = 0
    # #     windspdsavg,windspdsavgvec,winddirsavgvec,winddirsunitavgvec,windgusts,windgustsavg, \
    # #       usavg,vsavg,unit_usavg,unit_vsavg = \
    # #     dis.resamplewind(datetimesUTC,offset,winddirabss,windspds,'60S',gusts=True,
    # #                      gustintvstr='3S',center=False)

        datetimesUTC = conv_df.index.to_pydatetime()  # DSD_dict['convtimestamps']

        deployedtlist = []
        windspdtlist = []
        windspdavgvectlist = []
        winddirabstlist = []
        winddiravgvectlist = []
        temptlist = []
        dewpointtlist = []
        pressuretlist = []

        # Find the times closest to the sweeptimes
        # TODO: just resmpale to the sweeptimes for crying out loud.
        for sweeptime in radar_dict['sweeptimelist']:
            try:
                index = next(
                    i for i, t in enumerate(datetimesUTC) if np.abs(
                        (t - sweeptime).total_seconds()) <= resample_interval)
                deployedtlist.append(True)
            except Exception:
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
    dis_dict['convdata_at_sweeptimes'] = {
        'windspd': windspdlist, 'windspdavgvec': windspdavgveclist, 'winddirabs': winddirabslist,
        'winddiravgvec': winddiravgveclist, 'temp': templist, 'dewpoint': dewpointlist,
        'pressure': pressurelist, 'deployed': deployedlist
    }

    return dis_dict


def read_sweeps(radar_dict):
    """Reads sweeps from CFRadial files for a case as defined in the radar_dict dictionary.
       Stuffs these sweeps into the dictionary"""
    # Initial latitude, longitude, altitude set to None (will be read from files)
    radardir = radar_dict['radardir']
    radpathlist = glob.glob(radardir + '/*{}*nc'.format(radar_dict['radname']))
    radstarttime = radar_dict['radstarttimestamp']
    radstoptime = radar_dict['radstoptimestamp']
    fieldnames = radar_dict['fieldnames']
    el_req = radar_dict['el_req']

    # Now read in all the sweeps between radstarttime and radstoptime closest to the requested
    # elevation angle
    radstarttimedt = datetime.strptime(radstarttime, '%Y%m%d%H%M%S')
    radstoptimedt = datetime.strptime(radstoptime, '%Y%m%d%H%M%S')

    # outfieldnameslist = []
    radarsweeplist = []
    sweeptimelist = []

    for radpath in radpathlist:
        sweeptime = radar._getsweeptime(radpath)

        if radstarttimedt <= sweeptime and sweeptime <= radstoptimedt:
            radarsweep = radar.readCFRadial_pyART(el_req, radpath, sweeptime,
                                                  fieldnames, compute_kdp=False)
            radarsweeplist.append(radarsweep)
            sweeptimelist.append(sweeptime)

    # Sort the lists by increasing time since glob doesn't sort in any particular order
    sorted_sweeptimelist = sorted(sweeptimelist)
    sorted_radarsweeplist = [x for _, x in sorted(zip(sweeptimelist, radarsweeplist),
                                                  key=lambda pair: pair[0])]

    # Stuff the lists into the dictionary
    radar_dict['radarsweeplist'] = sorted_radarsweeplist
    radar_dict['sweeptimelist'] = sorted_sweeptimelist

    return radar_dict


def compute_storm_motion(radar_dict):
    """Computes storm motion from start and end points of a given feature.
       Adds a new key to the given dictionary with the storm motion vector
       in a tuple."""

    deltat = (radar_dict['feature_end_time'] - radar_dict['feature_start_time']).total_seconds()
    ustorm = (radar_dict['feature_end_loc'][0] - radar_dict['feature_start_loc'][0]) * 1000. / deltat
    vstorm = (radar_dict['feature_end_loc'][1] - radar_dict['feature_start_loc'][1]) * 1000. / deltat
    radar_dict['feature_motion'] = (ustorm, vstorm)

    return radar_dict


def get_dis_locs_relative_to_radar(dis_dict, radar_dict):
    """Gets the disdrometer locations relative to the radar."""
    dradloclist = []
    try:
        dgeoloclist = dis_dict['dgeoloclist']
    except KeyError:
        dis_dict = get_dis_geoloc(dis_dict)
        dgeoloclist = dis_dict['dgeoloclist']

    # Grab radar location information from the first sweep in radar_dict
    radarsweep = radar_dict['radarsweeplist'][0]
    rlat = radarsweep.latitude['data'][0]
    rlon = radarsweep.longitude['data'][0]

    for dgeoloc in dgeoloclist:
        dradx, drady = pyart.core.geographic_to_cartesian_aeqd(dgeoloc[1], dgeoloc[0], rlon, rlat)
        dx = dradx[0]
        dy = drady[0]
        print(dx, dy)
        dradloclist.append((dx, dy))

    # Stuff the locations into the dis_dict dictionary
    dis_dict['dradloclist'] = dradloclist

    return dis_dict


def get_dis_geoloc(dis_dict):
    """Get geographic coordinates of probes"""
    dgeoloclist = []
    # Extract parameters from dictionary
    dis_names = dis_dict['dis_names']
    dis_types = dis_dict['dis_types']
    disfilenames = dis_dict['disfilenames']
    try:
        convfilenames = dis_dict['convfilenames']
    except KeyError:
        convfilenames = [None] * len(dis_names)
    starttimes = dis_dict['starttimes']
    stoptimes = dis_dict['stoptimes']
    dis_dir = dis_dict['dis_dir']

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
        dgeoloclist.append(dgeoloc)

    dis_dict['dgeoloclist'] = dgeoloclist
    return dis_dict


def get_dis_locs_arps_real_grid(grid_dict, PIPS_geo_locs):
    """Returns the locations of the disdrometers within the ARPS grid. Only works for real cases
       with a map projection."""
    # Get model coordinates
    xc = grid_dict['xs']
    yc = grid_dict['ys']
    # Get basemap instance
    bgmap = grid_dict['bgmap']

    # Find coordinates of PIPS stations in the model
    dx = grid_dict['dx']
    dy = grid_dict['dy']
    modloc_list = []
    coord_list = []
    for i, PIPS_loc in enumerate(PIPS_geo_locs):
        xloc, yloc = bgmap(PIPS_loc[1], PIPS_loc[0])
        modloc_list.append(np.array([xloc, yloc]))
        iloc = (xloc - xc[1]) / dx
        jloc = (yloc - yc[1]) / dy
        coord_list.append(np.array([iloc, jloc]))

    return modloc_list, coord_list


def set_dh(model_dict, radar_dict, fixed_time=True):
    """Reads a DataHandler instance for COMMAS simulation output given information in model_dict.
       Puts the DataHandler instance into model_dict. Also associates
       the reference time of the model simulation with that of the reference sweep in radar_dict"""
    modelname = 'COMMAS'
    runname = model_dict['runname']
    dirname = model_dict['dirname']
    # model_dt = model_dict['model_dt']
    modeltimesec_ref = model_dict['modeltimesec_ref']
    if model_dict['composite']:
        model_times = model_dict['model_times']
    else:
        model_times = [modeltimesec_ref]
    multitime = model_dict.get('multitime', True)
    microphys = model_dict['microphys']
    dh = getDataHandler(modelname, dirname, model_times, microphys, multitime=multitime)
    dh.setRun(runname, 0)
    dh.loadTimes()
    dh.setTime(modeltimesec_ref)
    modeltime_ref = radar_dict['sweeptime_ref']
    model_dict['DataHandler'] = dh
    model_dict['modeltime_ref'] = modeltime_ref

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
    if(searchboxlims is None):
        ibgn = 0
        iend = var.shape[1] - 1
        jbgn = 0
        jend = var.shape[0] - 1
    else:
        ibgn = searchboxlims[0]
        iend = searchboxlims[1]
        jbgn = searchboxlims[2]
        jend = searchboxlims[3]

    if debug:
        print(ibgn, iend, jbgn, jend, var.shape)

        fig = plt.figure()
        ax = fig.add_subplot(111)

    varsearch = var[jbgn:jend + 1, ibgn:iend + 1]
    print(varsearch.shape)

    if debug:
        ax.imshow(varsearch, origin='lower')

    if(extremum == 'max'):
        var_extremum = varsearch.max()
        flatindex = np.argmax(varsearch)
    else:
        var_extremum = varsearch.min()
        flatindex = np.argmin(varsearch)

    # print varsearch.shape
    jrel, irel = np.unravel_index(flatindex, varsearch.shape)

    if debug:
        print("irel, jrel = ", irel, jrel)
        ax.plot([irel], [jrel], marker='o', color='red')

    if(guesscoords is None):
        iref = irel
        jref = jrel
    else:
        iref = ibgn + irel
        jref = jbgn + jrel

    if debug:
        print("iref, jref = ", iref, jref)

    return var_extremum, iref, jref


def get_composite_grid(grid_dict, compositedict):
    composite_grid_dict = {}
    ichw, jchw = compositedict['compositehw']
    dx = grid_dict['dx']
    dy = grid_dict['dy']
    xc_comp = np.arange(-ichw, ichw + 1) * dx
    yc_comp = np.arange(-jchw, jchw + 1) * dy
    xe_comp = np.arange(-ichw - 0.5, ichw + 1.5) * dx
    ye_comp = np.arange(-jchw - 0.5, jchw + 1.5) * dy

    nxm = xc_comp.shape[0]
    # nxe = xe_comp.shape[0]
    nym = yc_comp.shape[0]
    # nye = ye_comp.shape[0]
    nzm = grid_dict['zc1d'].shape[0]
    # nze = grid_dict['ze1d'].shape[0]

    # Recast xc,xe,yc,ye as 3D arrays
    composite_grid_dict['xc'] = np.lib.stride_tricks.as_strided(
        xc_comp, strides=((0, 0) + xc_comp.strides), shape=((nzm, nym) + xc_comp.shape))
    composite_grid_dict['xe'] = np.lib.stride_tricks.as_strided(
        xe_comp, strides=((0, 0) + xe_comp.strides), shape=((nzm, nym) + xe_comp.shape))
    composite_grid_dict['yc'] = np.lib.stride_tricks.as_strided(yc_comp, strides=(
        (0, ) + yc_comp.strides + (0, )), shape=((nzm, ) + yc_comp.shape + (nxm, )))
    composite_grid_dict['ye'] = np.lib.stride_tricks.as_strided(ye_comp, strides=(
        (0, ) + ye_comp.strides + (0, )), shape=((nzm, ) + ye_comp.shape + (nxm, )))
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
    xcorplot, ycorplot = pm.computecorners(
        composite_grid_dict['xeplot'], composite_grid_dict['yeplot'])
    composite_grid_dict['xcorplot'] = xcorplot
    composite_grid_dict['ycorplot'] = ycorplot

    return composite_grid_dict


def read_vardict(model_dict, varlists, varlistv, varlist_derived):
    """Reads in variables from the model files listed in model_dict and returns a list of
       dictionaries containing the variables for each time."""
    model_times = model_dict['model_times']
    dh = model_dict['DataHandler']
    vardictlist = []
    for t, time in enumerate(model_times):
        vardict = dh.loadVars(varlists)
        temp = dh.loadVars(varlistv)
        vardict.update(temp)
        # Extract only lowest model level here
        for key, var in vardict.items():
            vardict[key] = var[0, ...]
        # Now read in the microphysics scalars
        mp_data, consts = dh.loadMicrophysics()
        # Extract the lowest model level and store in mp_data_2D
        # mp_data_2D = {}
        for key, dat in mp_data.items():
            vardict[key] = dat.T[0, ...]
        if 'U' in varlistv:
            vardict['UC'] = 0.5 * (vardict['U'][:-1, :-1] + vardict['U'][:-1, 1:])
        if 'V' in varlistv:
            vardict['VC'] = 0.5 * (vardict['V'][:-1, :-1] + vardict['V'][1:, :-1])
        if 'PTE' in varlist_derived:
            try:
                vardict['PTE'] = thermo.calpte(vardict['P'], vardict['TH'], vardict['QV'])
            except BaseException:
                print("Cannot calculate PTE!")
        vardictlist.append(vardict)
    return vardictlist


def build_composite(model_dict, compositedict, dh, plotlocs=True):
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
    igend = gridlimindices[1] - 1
    jgbgn = gridlimindices[2]
    jgend = gridlimindices[3] - 1

    model_times = model_dict['model_times']
    multitime = model_dict.get('multitime', True)
    ntimes = model_times.size

    ireflist = []
    jreflist = []
    varlists = ['DBZ', 'TH', 'QV', 'P']
    varlists_derived = ['PTE', 'UC', 'VC']
    varlistv = ['U', 'V']
    varcompdict = {}
    for varname in varlists:
        varcompdict[varname] = np.zeros((ichw * 2 + 1, jchw * 2 + 1), dtype=np.float)
    for varname in varlists_derived:
        varcompdict[varname] = np.zeros((ichw * 2 + 1, jchw * 2 + 1), dtype=np.float)
    for varname in varlistv:
        varcompdict[varname] = np.zeros((ichw * 2 + 2, jchw * 2 + 2), dtype=np.float)

    if plotlocs:
        figcl = plt.figure()
        axcl = figcl.add_subplot(111)
        axclcolorcycle = [plt.cm.plasma(i) for i in np.linspace(0, 1, ntimes)]

    nt = 0
    tflag = True
    for t, time in enumerate(model_times):
        if t == 0:
            # Choose center of domain to start
            iref = int((igend - igbgn) / 2.)
            jref = int((jgend - jgbgn) / 2.)
            searchboxlims = [
                iref - ishw * 2,
                iref + ishw * 2 + 1,
                jref - jshw * 2,
                jref + jshw * 2 + 1]
            guesscoords = [iref, jref]

        # print "guesscoords = ",guesscoords
        # print "searchboxlims = ",searchboxlims

        print("The model time is {:d} s".format(int(time)))
        dh.setTime(time)
        if not multitime:
            dh.setRun(model_dict['runname'], 0, time=time)
        # Right now, only vortz or w for tracking variable
        if tracking_varname == 'w':
            print("Reading variable " + tracking_varname)
            trackvardict = dh.loadVars([tracking_varname])
            trackvar = trackvardict[tracking_varname][tracking_klvl, ...]
            trackvar = trackvar.squeeze()
        elif tracking_varname == 'vortz':
            print("Deriving variable " + tracking_varname)
            # TODO: Make this less brittle (maybe create relational dictionaries for variable names
            # for each model [ARPS, CM1, COMMAS], so that we can use the same variable names for
            # each and it automatically gets translated within the DataHandler code)
            trackvardict = dh.loadVars([tracking_varname])
            trackvar = trackvardict[tracking_varname][tracking_klvl, ...]
            trackvar = trackvar.squeeze()

        var_extremum, iref, jref = trackfeature(trackvar, searchboxlims=searchboxlims,
                                                guesscoords=guesscoords, extremum=tracking_extremum)
        guesscoords = [iref, jref]
        searchboxlims = [iref - ishw, iref + ishw + 1, jref - jshw, jref + jshw + 1]
        print(tracking_extremum + ' ' + tracking_varname +
              ' at height z = {:.1f} m'.format(tracking_level) + '\n' +
              'is {:.2f}'.format(var_extremum))
        print(
            "The location is x,y = {:.2f},{:.2f} (i,j = {:d},{:d})".format(
                xckm[iref],
                yckm[jref],
                iref,
                jref))

        ireflist.append(iref)
        jreflist.append(jref)

        xref = xckm[iref]
        yref = yckm[jref]
        # Restrict what goes into composite above a certain threshold
        if np.abs(var_extremum) >= tracking_thresh:
            if plotlocs:
                axcl.plot([xref], [yref], marker='o', color=axclcolorcycle[t])
            nt = nt + 1
            # Read in variables to composite (separate grid dimensions for scalars vs.
            # vector wind components)
            compositeboxindices = [iref - ichw + igbgn, iref + ichw + 1 +
                                   igbgn, jref - jchw + jgbgn, jref + jchw + 1 + jgbgn, 0, 1]
            compositeboxindices_uv = [
                iref - ichw + igbgn,
                iref + ichw + 1 + igbgn + 1,
                jref - jchw + jgbgn,
                jref + jchw + 1 + jgbgn + 1,
                0,
                1]
            print(compositeboxindices)
            print(compositeboxindices_uv)

            vardict = dh.loadVars(varlists, compositeboxindices)
            temp = dh.loadVars(varlistv, compositeboxindices_uv)
            vardict.update(temp)

            # Compute velocity components at scalar points
            vardict['UC'] = 0.5 * (vardict['U'][:-1, :-1] + vardict['U'][:-1, 1:])
            vardict['VC'] = 0.5 * (vardict['V'][:-1, :-1] + vardict['V'][1:, :-1])

            # Compute equivalent potential temperature
            vardict['PTE'] = thermo.calpte(vardict['P'], vardict['TH'], vardict['QV'])

            # Now read in the microphysics scalars
            mp_data, consts = dh.loadMicrophysics()
            # Extract the lowest model level and store in mp_data_2D
            mp_data_2D = {}
            for key, dat in mp_data.items():
                mp_data_2D[key] = dat.T[0, compositeboxindices[2]:compositeboxindices[3],
                                        compositeboxindices[0]:compositeboxindices[1]]
                if tflag:
                    varcompdict[key] = np.zeros((ichw * 2 + 1, jchw * 2 + 1), dtype=np.float)

            print(list(mp_data_2D.keys()))
            print(list(varcompdict.keys()))

            # Cumulative sum for each variable for each time
            for varname, var in vardict.items():
                varcompdict[varname] = varcompdict[varname] + var
            for varname, var in mp_data_2D.items():
                varcompdict[varname] = varcompdict[varname] + var
            tflag = False
    # Compute composite averages for each variable (nt is the total number of times)
    for varname, var in varcompdict.items():
        if model_times.size > 1:
            varcompdict[varname] = varcompdict[varname] / nt

    if plotlocs:
        axcl.set_aspect('equal')
        figcl.savefig(model_dict['runname'] + '_sfcvortloc_{:06d}_{:06d}.png'.format(
            int(model_times[0]), int(model_times[-1])), dpi=200)

    return varcompdict


def calc_obs_transect(dis_dict, dis_ts_model_dict, Dr, Dmax=None, calc_fits=True,
                      plot_transects=False):
    """Gathers and computes DSD quantities based on disdrometer observations along transects.
       Optionally computes exponential and gamma fits based on MofM and MofTM and plots them."""

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    D0r_obs = []
    ND_obs = []
    if calc_fits:
        #         D0r_obs_exp = []
        D0r_obs_gam = []
        D0r_obs_tmf = []
        ND_obs_gam = []
        ND_obs_tmf = []

    for d, dis_name in enumerate(dis_dict['dis_names']):
        sample_xlocs = np.array([xylocs[0] for xylocs in dis_ts_model_dict['dis_ts_xyslocs'][d]])
        ND = dis_dict['timeseries']['ND'][d]
        if calc_fits:
            rho = dis_dict['timeseries']['rho'][d]
            synthbins, exp_DSD, gam_DSD, tmf_DSD, dis_DSD = dis.calc_DSD(ND.T, rho)
            (ND_expDSD, N0_exp, lamda_exp, mu_exp, qr_exp, Ntr_exp, refl_DSD_exp, D_med_exp,
             D_m_exp) = exp_DSD
            (ND_gamDSD, N0_gam, lamda_gam, mu_gam, qr_gam, Ntr_gam, refl_DSD_gam, D_med_gam,
             D_m_gam, LWC_gam, rainrate_gam) = gam_DSD
            (ND_tmfDSD, N0_tmf, lamda_tmf, mu_tmf, qr_tmf, Ntr_tmf, refl_DSD_tmf, LWC_tmf,
             rainrate_tmf) = tmf_DSD
#             ND_expDSD = ND_expDSD.T
#             logND_expDSD = np.ma.log10(ND_expDSD / 1000.)  # Get to log(m^-3 mm^-1)
#             logND_expDSD = np.ma.masked_where(logND_expDSD <= -1.0, logND_expDSD)
            ND_gamDSD = ND_gamDSD[:, :Dmax_index]
            logND_gamDSD = np.ma.log10(ND_gamDSD / 1000.)  # Get to log(m^-3 mm^-1)
            logND_gamDSD = np.ma.masked_where(logND_gamDSD <= -1.0, logND_gamDSD)
            ND_tmfDSD = ND_tmfDSD[:, :Dmax_index]
            logND_tmfDSD = np.ma.log10(ND_tmfDSD / 1000.)  # Get to log(m^-3 mm^-1)
            logND_tmfDSD = np.ma.masked_where(logND_tmfDSD <= -1.0, logND_tmfDSD)

#             D0r_exp = np.zeros((np.size(np.array(sample_xlocs))))
            D0r_gam = np.zeros((np.size(np.array(sample_xlocs))))
            D0r_tmf = np.zeros_like(D0r_gam)
            for t in range(sample_xlocs.size):
                # D0r_exp[t] = dis.calc_D0_bin(D, Dl, Dr, ND_expDSD[t, :], bin_width)
                D0r_gam[t] = dis.calc_D0_bin(D[:Dmax_index], Dl[:Dmax_index], Dr[:Dmax_index],
                                             ND_gamDSD[t, :], bin_width[:Dmax_index])
                D0r_tmf[t] = dis.calc_D0_bin(D[:Dmax_index], Dl[:Dmax_index], Dr[:Dmax_index],
                                             ND_tmfDSD[t, :], bin_width[:Dmax_index])
#             D0r_obs_exp.append(D0r_exp)
            D0r_obs_gam.append(D0r_gam)
            D0r_obs_tmf.append(D0r_tmf)

        # FIXME: Have to do this afterwards because the call to calc_DSD uses the full diameter
        # range internally.
        ND = ND[:, :Dmax_index]
        logND = np.log10(ND)
        logND = np.ma.masked_where(logND <= -1.0, logND)
        # TODO: can simplify this by improving vectorization of calc_D0_bin
        D0r = np.zeros((np.size(np.array(sample_xlocs))))
        # Calculate median drop diameter
        for t in range(sample_xlocs.size):
            D0r[t] = dis.calc_D0_bin(D[:Dmax_index], Dl[:Dmax_index],
                                     Dr[:Dmax_index], ND[t, :], bin_width[:Dmax_index])
        D0r_obs.append(D0r)

        if plot_transects:
            # Set up the figure for whether or not we are calculating fits
            if calc_fits:
                plotcbar = False
                fig = plt.figure(figsize=(8, 9))
                ax = fig.add_subplot(311)
            else:
                plotcbar = True
                fig = plt.figure(figsize=(8, 3))
                ax = fig.add_subplot(111)

            # Plot the raw obs
            plotparamdicts = [{
                'type': 'pcolor',
                'vlimits': (-1.0, 3.0),
                'clabel': r'log[N ($m^{-3} mm^{-1}$)]',
                'plotcbar': plotcbar
            }]
            xvals = [sample_xlocs / 1000.]
        #     print "xvals = ", xvals
            yvals = [Dl[:Dmax_index] * 1000.]
            zvals = [logND.T]
            ax = pm.plotmeteogram(ax, xvals, zvals, plotparamdicts, yvals=yvals)
            axparamdicts = [{'majorxlocator': ticker.MultipleLocator(base=1.0),
                             'majorylocator': ticker.MultipleLocator(base=2.0),
                             'axeslimits': [None, (0.0, Dmax)],
                             'axeslabels': ['x position', 'D (mm)'],
                             'axesautofmt': False}]
            # FIXME
            # if not calc_fits:
            #     axlist = pm.set_meteogram_axes([ax], axparamdicts)
            ax.plot(xvals[0], D0r * 1000., c='k', ls='-', lw=1)

            # Plot additional panels for the fitted distributions (just gamma and
            # gamma with TMM for now)
            if calc_fits:
                ax2 = fig.add_subplot(312)
                plotparamdicts = [{
                    'type': 'pcolor',
                    'vlimits': (-1.0, 3.0),
                    'clabel': r'log[N ($m^{-3} mm^{-1}$)]',
                    'plotcbar': False
                }]
                xvals = [sample_xlocs / 1000.]
            #     print "xvals = ", xvals
                yvals = [Dl[:Dmax_index] * 1000.]
                zvals = [logND_gamDSD.T]
                ax2 = pm.plotmeteogram(ax2, xvals, zvals, plotparamdicts, yvals=yvals)
                axparamdicts.append({'majorxlocator': ticker.MultipleLocator(base=1.0),
                                     'majorylocator': ticker.MultipleLocator(base=2.0),
                                     'axeslimits': [None, (0.0, Dmax)],
                                     'axeslabels': ['x position', 'D (mm)'],
                                     'axesautofmt': False})
                ax2.plot(xvals[0], D0r_gam * 1000., c='k', ls='-', lw=1)

                ax3 = fig.add_subplot(313)
                plotparamdicts = [{
                    'type': 'pcolor',
                    'vlimits': (-1.0, 3.0),
                    'clabel': r'log[N ($m^{-3} mm^{-1}$)]',
                    'plotcbar': True
                }]
                xvals = [sample_xlocs / 1000.]
            #     print "xvals = ", xvals
                yvals = [Dl[:Dmax_index] * 1000.]
                zvals = [logND_tmfDSD.T]
                ax3 = pm.plotmeteogram(ax3, xvals, zvals, plotparamdicts, yvals=yvals)
                axparamdicts.append({'majorxlocator': ticker.MultipleLocator(base=1.0),
                                     'majorylocator': ticker.MultipleLocator(base=2.0),
                                     'axeslimits': [None, (0.0, Dmax)],
                                     'axeslabels': ['x position', 'D (mm)'],
                                     'axesautofmt': False})
                # axlist = pm.set_meteogram_axes([ax, ax2, ax3], axparamdicts)
                ax3.plot(xvals[0], D0r_tmf * 1000., c='k', ls='-', lw=1)

        ND_obs.append(ND)
        if calc_fits:
            ND_obs_gam.append(ND_gamDSD)
            ND_obs_tmf.append(ND_tmfDSD)

    # Compile all return items into a convenient dictionary
    transect_DSD_obs_dict = {'ND': ND_obs, 'D0r_obs': D0r_obs}
    if calc_fits:
        transect_DSD_obs_dict.update({'ND_gam': ND_obs_gam, 'D0r_gam': D0r_obs_gam,
                                      'ND_tmf': ND_obs_tmf, 'D0r_tmf': D0r_obs_tmf})

    return transect_DSD_obs_dict


def interp_model_to_transect(dis_dict, model_dict, dis_ts_model_dict, Dr,
                             sampling_interval=60., add_hail=False, use_bins_for_interp=False,
                             use_Parsivel_simulator=False, Dmax=None, plot_transects=False):
    """Interpolates model variables, including bulk DSDs, to a simulated disdrometer transect.
       Also bins the resulting DSDs using the Parsivel diameter bins. Options:
       add_hail: if True, add hail and graupel to the rain distribution
       use_bins_for_resample: if True, discretize the model gamma DSDs to the Parsivel bins
       first, and then interpolate (at sampling_interval) along the transect. If False, directly
       interpolate the moments (using weighted averages), and then discretize.
       use_Parsivel_simulator: if True, sample the model DSDs using the Parsivel simulator. If
       False, keep the raw DSDs."""

    Dmax, Dmax_index = get_Dmax_index(Dr, Dmax)

    # TODO: allow for adding hail but still
    if not use_bins_for_interp and add_hail:
        print("Directly using moments for resampling (instead of discretizing to bins first) \n"
              "is currently not working with the add hail option, sorry!"
              "Setting add_hail to False!")
        add_hail = False

    # FIXME: should probably rearrange things so that the timeseries for each variable is contained
    # in the dictionary of variables rather than the other way around
    D0r_mod = []
    ND_list = []
    if use_Parsivel_simulator:
        D0r_mod_ps = []
        ND_ps_list = []

    for d, dis_name in enumerate(dis_dict['dis_names']):
        sampling_times = dis_ts_model_dict['dis_ts_stimes'][d]
        all_times = dis_ts_model_dict['dis_ts_times'][d]
        dt = all_times[1:] - all_times[:-1]
        sample_xlocs = np.array([xylocs[0] for xylocs in dis_ts_model_dict['dis_ts_xyslocs'][d]])
        # sample_ylocs = np.array([xylocs[1] for xylocs in dis_ts_model_dict['dis_ts_xyslocs'][d]])
        vars_at_ts_points = dis_ts_model_dict['dis_ts_vars_points'][d]
        ntimes = len(vars_at_ts_points)

        rhoa = np.array([vars_at_ts_points[t]['rhoa'] for t in range(ntimes)])
        qr = np.array([vars_at_ts_points[t]['qr'] for t in range(ntimes)])
        ntr = np.array([vars_at_ts_points[t]['ntr'] for t in range(ntimes)])
        zr = np.array([vars_at_ts_points[t]['zr'] for t in range(ntimes)])
        # TODO: maybe change arguments of calls to cal_N0, cal_lambda, etc. to optionally directly
        # read in Z moments and calculate alpha within
        if use_bins_for_interp:
            alphar = np.array([vars_at_ts_points[t]['alphar'] for t in range(ntimes)])
            N0r, _ = dsd.cal_N0(rhoa, qr, ntr, cr, alphar)
            lamdar = dsd.cal_lamda(rhoa, qr, ntr, cr, alphar)
        # dBZ = np.array([vars_at_ts_points[t]['DBZ'] for t in range(ntimes)])

        if add_hail:
            qh = np.array([vars_at_ts_points[t]['qh'] for t in range(ntimes)])
            nth = np.array([vars_at_ts_points[t]['nth'] for t in range(ntimes)])
            # zh = np.array([vars_at_ts_points[t]['zh'] for t in range(ntimes)])
            alphah = np.array([vars_at_ts_points[t]['alphah'] for t in range(ntimes)])
            rhoh = np.array([vars_at_ts_points[t]['rhoh'] for t in range(ntimes)])
            ch = rhoh * np.pi / 6.

            qg = np.array([vars_at_ts_points[t]['qg'] for t in range(ntimes)])
            ntg = np.array([vars_at_ts_points[t]['ntg'] for t in range(ntimes)])
            # zg = np.array([vars_at_ts_points[t]['zg'] for t in range(ntimes)])
            alphag = np.array([vars_at_ts_points[t]['alphag'] for t in range(ntimes)])
            rhog = np.array([vars_at_ts_points[t]['rhog'] for t in range(ntimes)])
            cg = rhog * np.pi / 6.

            N0h, _ = dsd.cal_N0(rhoa, qh, nth, ch, alphah)
            lamdah = dsd.cal_lamda(rhoa, qh, nth, ch, alphah)
            N0g, _ = dsd.cal_N0(rhoa, qg, ntg, cg, alphag)
            lamdag = dsd.cal_lamda(rhoa, qg, ntg, cg, alphag)

        # If sampling the model DSD with the Parsivel simulator, we need the array of assumed
        # fall speeds vs. diameter, with the appropriate correction for air density
        if use_Parsivel_simulator:
            Vtr = []
            for t in range(ntimes):
                Vtr.append(dis.assignfallspeed(dis.avg_diameter, rhocorrect=True, rho=rhoa[t]))
            Vtr = np.array(Vtr)
            Vtr = Vtr[:, :Dmax_index]

        Nc_bin_tmp = np.empty((np.size(N0r), np.size(D[:Dmax_index])))
        Nc_bin = np.zeros((np.size(np.array(sampling_times)), np.size(D[:Dmax_index])))

        if use_bins_for_interp:
            for index, _ in np.ndenumerate(N0r):
                Nc_bin_tmp[index, :] = 1.e-3 * N0r[index] * \
                    (D[:Dmax_index])**alphar[index] * np.exp(-lamdar[index] * (D[:Dmax_index]))
                if add_hail:
                    Nc_bin_tmp[index, :] = Nc_bin_tmp[index, :] + 1.e-3 * N0h[index] * \
                        (D[:Dmax_index])**alphah[index] * np.exp(-lamdah[index] * (D[:Dmax_index]))
                    Nc_bin_tmp[index, :] = Nc_bin_tmp[index, :] + 1.e-3 * N0g[index] * \
                        (D[:Dmax_index])**alphag[index] * np.exp(-lamdag[index] * (D[:Dmax_index]))
        else:
            qr_stimes = np.zeros((np.size(np.array(sampling_times))))
            ntr_stimes = np.zeros_like(qr_stimes)
            zr_stimes = np.zeros_like(qr_stimes)
            N0r_stimes = np.zeros_like(qr_stimes)
            lamdar_stimes = np.zeros_like(qr_stimes)
            alphar_stimes = np.zeros_like(qr_stimes)
            rhoa_stimes = np.zeros_like(qr_stimes)

            # Just first sampling time
            Nc_bin_tmp[0, :] = 1.e-3 * N0r[0] * \
                (D[:Dmax_index])**alphar[0] * np.exp(-lamdar[0] * (D[:Dmax_index]))
            if add_hail:  # Not working right now for this option
                Nc_bin_tmp[0, :] = Nc_bin_tmp[0, :] + 1.e-3 * N0h[0] * \
                    (D[:Dmax_index])**alphah[0] * np.exp(-lamdah[0] * (D[:Dmax_index]))
                Nc_bin_tmp[0, :] = Nc_bin_tmp[0, :] + 1.e-3 * N0g[0] * \
                    (D[:Dmax_index])**alphag[0] * np.exp(-lamdag[0] * (D[:Dmax_index]))

        Nc_bin_tmp = np.ma.masked_invalid(Nc_bin_tmp)

        # Now loop through and "resample" the number concentrations as weighted averages
        # for the case that the disdrometer crosses a grid edge or the model DSD
        # changes during the sampling interval.
        sample_indices = np.searchsorted(all_times, sampling_times, side='left')
        if use_Parsivel_simulator:
            Nc_bin_tmp_ps = np.empty((np.size(lamdar), np.size(D[:Dmax_index])))
            Nc_bin_ps = np.zeros((np.size(np.array(sampling_times)), np.size(D[:Dmax_index])))
            # Special treatment for first sampling time. Just assume DSD valid at that time was
            # constant for the previous sampling interval
            sample_dict = create_random_gamma_DSD(ntr[0], lamdar[0],
                                                  alphar[0], Vtr[0], sampling_length,
                                                  sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                  sampling_interval=sampling_interval,
                                                  remove_margins=True, rhocorrect=True, rho=rhoa[0])
            ND_sample = sample_dict['ND']
            pcount_binned_sample = sample_dict['pcount_binned']
            if add_hail:
                sample_dict_h = create_random_gamma_DSD(nth[0], lamdah[0],
                                                        alphah[0], Vtr[0], sampling_length,
                                                        sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                        sampling_interval=sampling_interval,
                                                        remove_margins=True, rhocorrect=True,
                                                        rho=rhoa[0])

                sample_dict_g = create_random_gamma_DSD(ntg[0], lamdag[0],
                                                        alphag[0], Vtr[0], sampling_length,
                                                        sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                        sampling_interval=sampling_interval,
                                                        remove_margins=True, rhocorrect=True,
                                                        rho=rhoa[0])

                ND_sample = ND_sample + sample_dict_h['ND'] + sample_dict_g['ND']
                pcount_binned_sample = pcount_binned_sample + \
                    sample_dict_h['pcount_binned'] + sample_dict_g['pcount_binned']

            Nc_bin_tmp_ps[0, :] = 1.e-3 * ND_sample
            Nc_bin_ps[0, :] = Nc_bin_tmp_ps[0, :]

            pcount_binned_samples = []
            for index, _ in np.ndenumerate(lamdar[:-1]):
                sample_dict = create_random_gamma_DSD(ntr[index], lamdar[index],
                                                      alphar[index], Vtr[index], sampling_length,
                                                      sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                      sampling_interval=dt[index],
                                                      remove_margins=True, rhocorrect=True,
                                                      rho=rhoa[index])
                ND_sample = sample_dict['ND']
                pcount_binned_sample = sample_dict['pcount_binned']

                if add_hail:
                    sample_dict_h = create_random_gamma_DSD(nth[index], lamdah[index],
                                                            alphah[index], Vtr[index],
                                                            sampling_length,
                                                            sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                            sampling_interval=dt[index],
                                                            remove_margins=True, rhocorrect=True,
                                                            rho=rhoa[index])

                    sample_dict_g = create_random_gamma_DSD(ntg[index], lamdag[index],
                                                            alphag[index], Vtr[index],
                                                            sampling_length,
                                                            sampling_width, Dl, D, Dr, Dmax=Dmax,
                                                            sampling_interval=sampling_interval,
                                                            remove_margins=True, rhocorrect=True,
                                                            rho=rhoa[index])
                    ND_sample = ND_sample + sample_dict_h['ND'] + sample_dict_g['ND']
                    pcount_binned_sample = pcount_binned_sample + \
                        sample_dict_h['pcount_binned'] + sample_dict_g['pcount_binned']

                pcount_binned_samples.append(sample_dict['pcount_binned'])
                Nc_bin_tmp_ps[index, :] = 1.e-3 * ND_sample

            pcount_binned_samples = np.array(pcount_binned_samples)
            Nc_bin_tmp_ps = np.ma.masked_invalid(Nc_bin_tmp_ps)
        else:
            # Special treatment for first sampling time. Just assume DSD valid at that time was
            # constant for the previous sampling interval
            Nc_bin[0, :] = Nc_bin_tmp[sample_indices[0], :]

        for s, sample_index in enumerate(sample_indices[:-1]):
            sample_index_end = sample_indices[s + 1]
            current_sample_indices = slice(sample_index, sample_index_end, None)

            if use_bins_for_interp:
                Nc_bin[s + 1, :] = np.sum(Nc_bin_tmp[current_sample_indices, :] *
                                          dt[current_sample_indices, None], axis=0) / \
                                          sampling_interval
            else:
                # Weighted average of moments
                qr_stimes[s + 1] = np.sum(qr[current_sample_indices] *
                                          dt[current_sample_indices]) / sampling_interval
                ntr_stimes[s + 1] = np.sum(ntr[current_sample_indices] *
                                           dt[current_sample_indices]) / sampling_interval
                zr_stimes[s + 1] = np.sum(zr[current_sample_indices] *
                                          dt[current_sample_indices]) / sampling_interval
                rhoa_stimes[s + 1] = np.sum(rhoa[current_sample_indices] *
                                            dt[current_sample_indices]) / sampling_interval
                # Re-compute N0r, lamdar, and alphar for weighted average DSD
                alphar_stimes[s + 1] = dualpol.solve_alpha_iter(rhoa_stimes[s + 1], mur,
                                                                qr_stimes[s + 1], ntr_stimes[s + 1],
                                                                zr_stimes[s + 1], rhorcst)
                N0r_stimes[s + 1], _ = dsd.cal_N0(rhoa_stimes[s + 1], qr_stimes[s + 1],
                                                  ntr_stimes[s + 1], cr, alphar_stimes[s + 1])
                lamdar_stimes[s + 1] = dsd.cal_lamda(rhoa_stimes[s + 1], qr_stimes[s + 1],
                                                     ntr_stimes[s + 1], alphar_stimes[s + 1])
                # Now compute discretized distribution
                Nc_bin[s + 1, :] = 1.e-3 * N0r_stimes[s + 1] * \
                    (D[:Dmax_index])**alphar_stimes[s + 1] * \
                    np.exp(-lamdar_stimes[s + 1] * (D[:Dmax_index]))
            if use_Parsivel_simulator:
                # Take weighted average of Vtr over each sample interval for now and use that to
                # calculate the sampling volumes
                Vtr_mean = np.sum(Vtr[current_sample_indices, :] *
                                  dt[current_sample_indices, None], axis=0) / sampling_interval
                sampling_volumes_D = calc_sampling_volumes_D(Vtr_mean, Dr, Dmax, sampling_interval,
                                                             sampling_area)

                pcount_binned = np.sum(pcount_binned_samples[current_sample_indices], axis=0)
                Nc_bin_ps[s + 1, :] = 1.e-3 * \
                    calc_ND(pcount_binned, sampling_volumes_D, Dr, Dl, Dmax)

        logNc_bin = np.log10(Nc_bin)
        logNc_bin = np.ma.masked_where(logNc_bin <= -1.0, logNc_bin)

        if use_Parsivel_simulator:
            Nc_bin_ps = np.ma.masked_invalid(Nc_bin_ps)
            logNc_bin_ps = np.log10(Nc_bin_ps)
            logNc_bin_ps = np.ma.masked_where(logNc_bin_ps <= -1.0, logNc_bin_ps)

        # Calculate median volume diameter
        D0r = np.zeros((np.size(np.array(sampling_times))))
        if use_Parsivel_simulator:
            D0r_ps = np.zeros_like(D0r)
        for t in range(sampling_times.size):
            D0r[t] = dis.calc_D0_bin(D[:Dmax_index], Dl[:Dmax_index], Dr[:Dmax_index], Nc_bin[t, :],
                                     bin_width[:Dmax_index])
            if use_Parsivel_simulator:
                D0r_ps[t] = dis.calc_D0_bin(D[:Dmax_index], Dl[:Dmax_index], Dr[:Dmax_index],
                                            Nc_bin_ps[t, :], bin_width[:Dmax_index])

        D0r_mod.append(D0r)
        if use_Parsivel_simulator:
            D0r_mod_ps.append(D0r_ps)

        if plot_transects:
            # Set up the figure for whether or not we are using the Parsivel simulator
            if use_Parsivel_simulator:
                plotcbar = False
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(211)
            else:
                plotcbar = True
                fig = plt.figure(figsize=(8, 3))
                ax = fig.add_subplot(111)

            plotparamdicts = [{
                'type': 'pcolor',
                'vlimits': (-1.0, 3.0),
                'clabel': r'log[N ($m^{-3} mm^{-1}$)]',
                'plotcbar': plotcbar
            }]
            xvals = [sample_xlocs / 1000.]
        #     print "xvals = ", xvals
            yvals = [Dl[:Dmax_index] * 1000.]
            zvals = [logNc_bin.T]
            ax = pm.plotmeteogram(ax, xvals, zvals, plotparamdicts, yvals=yvals)
            axparamdicts = [{'majorxlocator': ticker.MultipleLocator(base=1.0),
                             'majorylocator': ticker.MultipleLocator(base=2.0),
                             'axeslimits': [None, (0.0, Dmax)],
                             'axeslabels': ['x position', 'D (mm)'],
                             'axesautofmt': False}]
            # FIXME
            # if not use_Parsivel_simulator:
            #     axlist = pm.set_meteogram_axes([ax], axparamdicts)
            ax.plot(xvals[0], D0r * 1000., c='k', ls='-', lw=1)

            # Plot an additional panel for the sampled DSD from the Parsivel simulator
            if use_Parsivel_simulator:
                ax2 = fig.add_subplot(212)
                plotparamdicts = [{
                    'type': 'pcolor',
                    'vlimits': (-1.0, 3.0),
                    'clabel': r'log[N ($m^{-3} mm^{-1}$)]',
                    'plotcbar': True
                }]
                xvals = [sample_xlocs / 1000.]
            #     print "xvals = ", xvals
                yvals = [Dl[:Dmax_index] * 1000.]
                zvals = [logNc_bin_ps.T]
                ax2 = pm.plotmeteogram(ax2, xvals, zvals, plotparamdicts, yvals=yvals)
                axparamdicts.append({'majorxlocator': ticker.MultipleLocator(base=1.0),
                                     'majorylocator': ticker.MultipleLocator(base=2.0),
                                     'axeslimits': [None, (0.0, Dmax)],
                                     'axeslabels': ['x position', 'D (mm)'],
                                     'axesautofmt': False})
                # axlist = pm.set_meteogram_axes([ax, ax2], axparamdicts)
                ax2.plot(xvals[0], D0r_ps * 1000., c='k', ls='-', lw=1)
        #     for ax in axlist:
        #         ax.invert_xaxis()
            # plt.savefig('raw_modelDSD.png',dpi=300)

        ND_list.append(Nc_bin)
        if use_Parsivel_simulator:
            ND_ps_list.append(Nc_bin_ps)

    # Compile all return items into a convenient dictionary
    transect_DSD_dict = {'ND': ND_list, 'D0r': D0r_mod}
    if use_Parsivel_simulator:
        transect_DSD_dict.update({'ND_ps': ND_ps_list, 'D0r_ps': D0r_mod_ps})

    return transect_DSD_dict


def get_ARPS_member_dir_and_prefix(member, cycle):
    """
    Gets the proper form for the subdirectory and file prefix name given a member number
    and cycle type (either 'posterior' or 'prior'). member number 0 is interpreted as the mean.
    """
    if member == 0:
        if cycle in 'posterior':
            member_dir = 'ENamean'
            member_prefix = 'enmean'
        elif cycle in 'prior':
            member_dir = 'ENfmean'
            member_prefix = 'efmean'
    else:
        if cycle in 'posterior':
            member_dir = 'EN{:03d}'.format(int(member))
            member_prefix = 'ena{:03d}'.format(int(member))
        elif cycle in 'prior':
            member_dir = 'ENF{:03d}'.format(int(member))
            member_prefix = 'enf{:03d}'.format(int(member))

    return member_dir, member_prefix


def read_ARPS_ensemble(f, member_list, f_args=None, f_kwargs=None, iterate_over='member_list',
                       process_parallel=True, n_jobs=5, verbose=0):
    """
    Reads multiple runs either in serial or parallel (using joblib). TODO. This is duplicating some
    functionality written by T. Supinie for pycaps (util.run_concurrently). Maybe just use that
    instead?

    Parameters
    ----------
    f : callable
        A function used to read in data from a particular run. It must have a positional argument
        that will be iterated over, the name of which is given by "iterate_over"
    member_list : list
        List of member numbers to read
    args : tuple
        List of positional arguments to pass to f
    iterate_over : str
        The name of the positional argument in args to iterate over
    process_parallel : bool
        Whether to process the runs in parallel or not (using joblib)
    n_jobs : int
        Number of jobs to run in parallel
    verbose : int
        verbosity level passed on to joblib.Parallel
    kwargs : dict
        List of keyword arguments to pass to f
    Returns
    -------
    list
        A list of updated run dictionaries

    """
    # The following avoids this gotcha:
    # https://pythonconquerstheuniverse.wordpress.com/2012/02/15/mutable-default-arguments/
    if f_args is None:
        f_args = ()
    if f_kwargs is None:
        f_kwargs = {}
    # Get a list of the function argument names, and find the one that we want to iterate over
    f_argnames = inspect.getargspec(f).args
    try:
        arg_to_iterate = f_argnames.index(iterate_over)
    except ValueError:
        print("{} not found in the function argument list! Stopping!".format(iterate_over))
        return

    # Define a wrapper function to replace the value of the iterated argument in f with a new one
    # Not sure if this is the most elegant way to do this, but it works and lets us use the
    # joblib.Parallel function with list comprehension for the "val" argument
    def g(f, val):
        new_args = tuple(val if i == arg_to_iterate else f_args[i] for i in range(len(f_args)))
        return f(*new_args, **f_kwargs)

    if process_parallel:
        runs = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(g)(f, val) for val in member_list)
    else:
        runs = []
        for val in member_list:
            print("Reading files for member {:d}".format(val))
            runs.append(g(f, val))
    return runs


def read_ARPS_member_data(basedir, expname, member, cycle, fileformat, time_range, varnames,
                          filetype='history', ibgn=None, iend=None, jbgn=None, jend=None, klvls=[1],
                          nproc_x=1, nproc_y=1, dump_nc_output=True, ncdir=None, fileprefix='',
                          grid_dict=None, datetime_range=None, x=None, y=None, mid_diameters=None):
    """[summary]

    Parameters
    ----------
    basedir : [type]
        [description]
    expname : [type]
        [description]
    member : [type]
        [description]
    cycle : [type]
        [description]
    fileformat : [type]
        [description]
    time_range : [type]
        [description]
    varnames : [type]
        [description]
    filetype : str, optional
        [description], by default 'history'
    ibgn : [type], optional
        [description], by default None
    iend : [type], optional
        [description], by default None
    jbgn : [type], optional
        [description], by default None
    jend : [type], optional
        [description], by default None
    klvls : list, optional
        [description], by default [1]
    nproc_x : int, optional
        [description], by default 1
    nproc_y : int, optional
        [description], by default 1
    dump_nc_output : bool, optional
        [description], by default True
    ncdir : [type], optional
        [description], by default None
    datetime_range : [type], optional
        [description], by default None
    x : [type], optional
        [description], by default None
    y : [type], optional
        [description], by default None
    mid_diameters : [type], optional
        [description], by default None
    """
    print("Loading member #{:d}".format(member))
    # Get member run prefix
    member_dir, member_prefix = get_ARPS_member_dir_and_prefix(member, cycle)
    member_absdir = os.path.join(basedir, expname, member_dir)
    vardict_list = []
    # all_files_exist_list = []
    for time in time_range:
        print("Loading time ", time)
        # Get member file path
        filepath = arps_read.get_file_path(member_absdir, member_prefix, fileformat, time=time,
                                           filetype='history')
#         all_files_exist = True
#         for ip in range(1, nproc_x + 1):
#             for jp in range(1, nproc_y + 1):
#                 filenamepatch = filepath + "_%03d" % ip + "%03d" % jp
#                 if not os.path.exists(filenamepatch):
#                     print("File {} does not exist!".format(filenamepatch))
#                     all_files_exist = False
#         all_files_exist_list.append(all_files_exist)
        # Get variables from file
        print(filepath)
        vardict = arps_read.read_hdfvars(filepath, varnames, ibgn=ibgn, jbgn=jbgn, iend=iend,
                                         jend=jend, klvls=klvls, nproc_x=nproc_x, nproc_y=nproc_y)
        vardict_list.append(vardict)

    if dump_nc_output:
        output_prefix = fileprefix + member_prefix
        dump_ARPS_xyslice_nc(ncdir=ncdir, vardict_list=vardict_list, member_prefix=output_prefix,
                             time=datetime_range, x=x, y=y, mid_diameters=mid_diameters, ibgn=ibgn,
                             iend=iend, jbgn=jbgn, jend=jend, grid_dict=grid_dict)

    return vardict_list  # all_files_exist_list


def dump_ARPS_xyslice_nc(ncdir=None, vardict_list=None, member_prefix=None, time=None, x=None,
                         y=None, mid_diameters=None, ibgn=None, iend=None, jbgn=None, jend=None,
                         grid_dict=None):
    """[summary]

    Parameters
    ----------
    ncdir : [type], optional
        [description], by default None
    vardict_list : [type], optional
        [description], by default None
    member_prefix : [type], optional
        [description], by default None
    time : [type], optional
        [description], by default None
    x : [type], optional
        [description], by default None
    y : [type], optional
        [description], by default None
    mid_diameters : [type], optional
        [description], by default None
    ibgn : [type], optional
        [description], by default None
    iend : [type], optional
        [description], by default None
    jbgn : [type], optional
        [description], by default None
    jend : [type], optional
        [description], by default None
    """

    xc = grid_dict['xc']
    yc = grid_dict['yc']
    xe = grid_dict['xe']
    ye = grid_dict['ye']
    coord_dict = {
        'time': time,
        'yc': ('yc', yc),
        'xc': ('xc', xc),
        'ye': ('ye', ye),
        'xe': ('xe', xe)
        }
    # First, create a dict of lists out of the above list of dicts
    vardict_combined = CRMutils.make_dict_of_lists(vardict_list)
    # Set things up for creating the xr Dataset
    for varname, var in vardict_combined.items():
        var_arr = np.array(var).T.squeeze()
        var_arr = np.rollaxis(var_arr, 2, 0)
        # Trim variables down to just the patch we want to work with
        if varname == 'u':
            var_arr_patch = var_arr[:, jbgn:jend+1, ibgn:iend+2]
            vardict_combined[varname] = (['time', 'yc', 'xe'], var_arr_patch)
        elif varname == 'v':
            var_arr_patch = var_arr[:, jbgn:jend+2, ibgn:iend+1]
            vardict_combined[varname] = (['time', 'ye', 'xc'], var_arr_patch)
        else:
            var_arr_patch = var_arr[:, jbgn:jend+1, ibgn:iend+1]
            vardict_combined[varname] = (['time', 'yc', 'xc'], var_arr_patch)

    # Create an xarray Dataset out of the variable dictionary
    var_ds = xr.Dataset(vardict_combined, coords=coord_dict)

    # Compute raw model DSD paramters and add them to the model Dataset
    rhor = 1000.
    cr = np.pi / 6. * rhor
    var_ds['rho'] = thermo.calrho(var_ds['p'], var_ds['pt'], var_ds['qv'])

    # Shape parameter
    var_ds['alphar'] = dsd.solve_alpha(var_ds['rho'], cr, var_ds['qr'], var_ds['nr'], var_ds['zr'])
    # var_ds['alphar'] = var_ds['alphar'].interpolate_na()
    # Intercept parameter
    var_ds['N0r'] = dsd.calc_N0_gamma(var_ds['rho'], var_ds['qr'], var_ds['nr'], cr, var_ds['alphar'])
    # Slope parameter
    var_ds['lamdar'] = dsd.calc_lamda_gamma(var_ds['rho'], var_ds['qr'], var_ds['nr'], cr, var_ds['alphar'])

    # Try computing ND and logND here

    # Broadcast DSD parameter DataArrays to get everyone on the same dimensional page
    mid_diameters, N0r_model_da, lamdar_model_da, alphar_model_da = \
        xr.broadcast(mid_diameters, var_ds['N0r'], var_ds['lamdar'], var_ds['alphar'])

    # Transpose these DataArrays to get time as the first dimension
    mid_diameters = mid_diameters.transpose('time', 'diameter_bin', 'yc', 'xc')
    N0r_model_da = N0r_model_da.transpose('time', 'diameter_bin', 'yc', 'xc')
    lamdar_model_da = lamdar_model_da.transpose('time', 'diameter_bin', 'yc', 'xc')
    alphar_model_da = alphar_model_da.transpose('time', 'diameter_bin', 'yc', 'xc')

    ND_model = dsd.calc_binned_DSD_from_params(N0r_model_da, lamdar_model_da, alphar_model_da,
                                               mid_diameters)
    ND_model = ND_model.fillna(0.0)
    logND_model = np.log10(ND_model)
    logND_model = logND_model.where(logND_model > -np.inf)
    #logND_model = logND_model.fillna(0.0)
    #logND_model = logND_model.where(logND_model > -1.0)

    var_ds['ND'] = ND_model
    var_ds['logND'] = logND_model

    var_ds.attrs['nx_full'] = grid_dict['nx_full']
    var_ds.attrs['ny_full'] = grid_dict['ny_full']
    var_ds.attrs['dx'] = grid_dict['dx']
    var_ds.attrs['dy'] = grid_dict['dy']
    var_ds.attrs['ctrlat'] = grid_dict['ctrlat']
    var_ds.attrs['ctrlon'] = grid_dict['ctrlon']
    var_ds.attrs['trulat1'] = grid_dict['trulat1']
    var_ds.attrs['trulat2'] = grid_dict['trulat2']
    var_ds.attrs['trulon'] = grid_dict['trulon']

    # print(var_ds)
    # Save Dataset to nc file
    filename = "{}_fields.nc".format(member_prefix)
    filepath = os.path.join(ncdir, filename)
    var_ds.to_netcdf(filepath)
