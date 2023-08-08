# plotmodule.py: A module containing some functions related to plotting model output
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable, host_subplot
import matplotlib.ticker as ticker
from matplotlib.collections import Collection, LineCollection
from matplotlib.artist import allow_rasterization
from matplotlib.font_manager import FontProperties
# import ctablesfrompyesviewer as ctables
from metpy.plots import ctables
from . import timemodule as tm
from . import PIPS as pips
from itertools import cycle
import pandas as pd
import xarray.plot as xrplot

# Set global font size for axes and colorbar labels, etc.

font = {'size': 9}
matplotlib.rc('font', **font)

fontP = FontProperties()
fontP.set_size('x-small')  # 'small'

# Contour levels for reflectivity (dBZ)
clevels_ref = np.arange(0.0, 85.0, 5.0)
clevels_zdr = np.arange(0.0, 6.25, 0.25)         # Contour levels for Zdr (dB)
clevels_vr = np.arange(-40.0, 41.0, 1.0)        # Contour levels for Vr (m/s)
# cmapdBZ = ctables.__getattribute__('REF_default')
normdBZ, cmapdBZ = ctables.registry.get_with_steps('NWSReflectivity', 5., 5.)
cmapzdr = cm.Reds
cmapvr = cm.RdBu_r

RdBu1_colors = ("#84001D", "#890025", "#8D002D", "#920034", "#97103C", "#9C2143",
                "#A12E4B", "#A63A53", "#AA455B", "#AF5063", "#B45A6C", "#B96575",
                "#BE707E", "#C37B88", "#C88792", "#CD949D", "#D2A1A9", "#D7B0B6",
                "#DDC1C5", "#E3D6D8", "#DAD8E2", "#CAC5DA", "#BDB6D4", "#B2A9CF",
                "#A89DCA", "#9F91C6", "#9687C1", "#8E7DBE", "#8773BA", "#7F69B7",
                "#7960B4", "#7257B1", "#6D4DAF", "#6743AE", "#6239AE", "#5E2CAE",
                "#5A1AB1", "#5800B7", "#5800C2", "#5C00D9")

# Define meteogram plot parameters below. These are aggregated and used in the meteogram
# plotting functions below. In the future, it may be a good idea to make a "Meteogram"
# class to encapsulate this stuff.

# Windspeed, gusts, and direction

windspeed_params = {'type': 'fill_between', 'linestyle': '-',
                    'color': 'b', 'alpha': 0.5, 'plotmin': 0}

windgust_params = {'type': 'line', 'linestyle': 'none',
                   'marker': 'o', 'color': 'r', 'ms': 2}

winddir_params = {'type': 'line', 'linestyle': 'none',
                  'marker': 'o', 'color': 'k', 'ms': 2}

compass_dir_params = {'type': 'line', 'linestyle': 'none',
                      'marker': 'o', 'color': 'k', 'ms': 2}

winddiag_params = {'type': 'vertical line',
                   'linestyle': '-', 'color': 'g', 'linewidth': 0.5}

# Temperature, dewpoint, and relative humidity

temp_params = {'type': 'fill_between', 'linestyle': '-',
               'color': 'r', 'alpha': 0.5, 'plotmin': -10}

dewpoint_params = {'type': 'fill_between', 'linestyle': '-',
                   'color': 'g', 'alpha': 0.5, 'plotmin': -10}

RH_params = {'type': 'fill_between', 'linestyle': '-',
             'color': 'g', 'alpha': 0.5, 'plotmin': 0}

# Battery voltage

battery_params = {'type': 'fill_between', 'linestyle': '-',
                  'color': 'brown', 'alpha': 1.0, 'plotmin': 8}

# GPS speed

GPS_speed_params = {'type': 'line', 'linestyle': '-', 'color': 'b'}

# Pressure

pressure_params = {'type': 'fill_between', 'linestyle': '-',
                   'color': 'purple', 'alpha': 0.5, 'plotmin': 0}

# Rain rate
rainrate_params = {'type': 'fill_between', 'linestyle': '-',
                   'color': 'g', 'alpha': 1.0, 'plotmin': 0}

# Reflectivity
reflectivity_params = {'type': 'fill_between', 'linestyle': '-',
                       'color': 'b', 'alpha': 1.0, 'plotmin': 0}

# Particle counts
pcount_params = {'type': 'fill_between', 'linestyle': '-',
                 'color': 'r', 'alpha': 0.5, 'plotmin': 0, 'label': r'Particle count (internal)'}

pcount2_params = {'type': 'fill_between', 'linestyle': '-',
                  'color': 'b', 'alpha': 0.5, 'plotmin': 0, 'label': r'Particle count (calculated)'}

# Signal amplitude
amplitude_params = {'type': 'fill_between', 'linestyle': '-',
                    'color': 'b', 'alpha': 1.0, 'plotmin': 0}


# The following Class and method are from
# http://stackoverflow.com/questions/12583970/matplotlib-contour-plots-as-postscript
# and allow rasterization of a contour or filled contour plot

class ListCollection(Collection):
    def __init__(self, collections, **kwargs):
        Collection.__init__(self, **kwargs)
        self.set_collections(collections)

    def set_collections(self, collections):
        self._collections = collections

    def get_collections(self):
        return self._collections

    @allow_rasterization
    def draw(self, renderer):
        for _c in self._collections:
            _c.draw(renderer)


def insert_rasterized_contour_plot(c):
    collections = c.collections
    for _c in collections:
        _c.remove()
    cc = ListCollection(collections, rasterized=True)
    ax = plt.gca()
    ax.add_artist(cc)
    return cc


def mtokm(val, pos):
    """Convert m to km for formatting axes tick labels"""
    val = val / 1000.0
    return '%i' % val


def mtokmr1(val, pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val = val / 1000.0
    # return '{:2.{prec}f}'.format(val,prec=prec)
    return '%2.1f' % val


def mtokmr2(val, pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val = val / 1000.0
    # return '{:2.{prec}f}'.format(val,prec=prec)
    return '%2.2f' % val


'''
Defines a function colorline that draws a (multi-)colored 2D line with coordinates x and y.
The color is taken from optional data in z, and creates a LineCollection.

z can be:
- empty, in which case a default coloring will be used based on the position along the input arrays
- a single number, for a uniform color [this can also be accomplished with the usual plt.plot]
- an array of the length of at least the same length as x, to color according to this data
- an array of a smaller length, in which case the colors are repeated along the curve

The function colorline returns the LineCollection created, which can be modified afterwards.

See also: plt.streamplot
'''

# Data manipulation:


def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    return segments


# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3,
              alpha=1.0, ax=None):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    if ax is None:
        ax = plt.gca()
    ax.add_collection(lc)

    return lc


def clear_frame(ax=None):
    # Taken from a post by Tony S Yu
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plotsingle(fig, axes, ptype, xs, ys, x, y, xlim, ylim, field, clevels, cmap, fieldnorm,
               cbarlevels, clabel, cformat, ovrmap, gis_info, numovr, ovrx, ovry, ovrfields,
               ovrfieldlvls, ovrfieldcolors, axesticks, rasterized=True):
    """Plot a single-paneled figure from model output (possibly staggered grid)"""
    if fig is None:
        fig = plt.figure()
    if axes is None:
        axes = fig.add_subplot(111)
    if rasterized:
        axes.set_rasterization_zorder(2)
    norm = fieldnorm
    if ptype == 1:  # Contour plot
        plot = axes.contourf(xs, ys, field, levels=clevels, cmap=cmap, norm=norm, zorder=1)
    elif ptype == 2:  # pcolor plot
        # Mask values below lower bounds
        field = np.ma.masked_where(field < clevels[0], field)
        # norm = matplotlib.colors.BoundaryNorm(clevels,cmap.N)
        # This is a hack to force masked areas to be white.  By default masked areas are
        # transparent, but the Agg backends apparently have problems when converting pcolor images
        # to eps and changes the transparent color to black.
        cmap.set_bad('white', alpha=None)
        plot = axes.pcolormesh(x, y, field, vmin=clevels[0], vmax=clevels[-1], cmap=cmap, norm=norm,
                               edgecolors='None', antialiased=False, rasterized=rasterized)
    # Find a nice number of ticks for the colorbar

    if cbarlevels is None:
        cintv = clevels[1] - clevels[0]
        cintvs = np.arange(clevels[0], clevels[-1], cintv)
        while True:
            if cintvs.size > 20:
                cintv = (cintvs[1] - cintvs[0]) * 2.
                cintvs = np.arange(cintvs[0], cintvs[-1], cintv)
            else:
                break
        cbarlevels = ticker.MultipleLocator(base=cintv)

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(plot, orientation='vertical', ticks=cbarlevels, cax=cax, format=cformat)
    if clabel is not None:
        cax.set_ylabel(clabel)
    if numovr > 0:
        for i in range(numovr):
            #             if(i != 1): # Hack!
            axes.contour(ovrx[i], ovry[i], ovrfields[i], levels=ovrfieldlvls[i],
                         colors=ovrfieldcolors[i], lw=2)
    if gis_info is not None:
        axes.plot([gis_info[1]], [gis_info[2]], 'ko')
    if xlim is None:
        xlim = [x.min(), x.max()]
    if ylim is None:
        ylim = [y.min(), y.max()]
    print(xlim, ylim)
    axes.set_xlim(xlim[0], xlim[1])
    axes.set_ylim(ylim[0], ylim[1])
    if axesticks[0] < 1000.:
        if axesticks[0] < 500.:
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.xaxis.set_major_formatter(formatter)
    if axesticks[1] < 1000.:
        if axesticks[0] < 500.:
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.yaxis.set_major_formatter(formatter)
    print(axesticks[0], axesticks[1])
    axes.xaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[0]))
    axes.yaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[1]))
    axes.set_xlabel('km')
    axes.set_ylabel('km')
    if axesticks[1] == axesticks[0]:
        axes.set_aspect('equal')
    else:
        axes.set_aspect('auto')

#   if(ovrmap): # Overlay map
#       readshapefile(track_shapefile_location,'track',drawbounds=True,linewidth=0.5,color='black'
#                     ax=axes)
# readshapefile(county_shapefile_location,'counties',drawbounds=True,
# linewidth=0.5, color='gray',ax=axes)  #Draws US county boundaries.

    return fig, axes


def plotsingle2(fig, axes, ptype, xs, ys, x, y, xlim, ylim, field, clevels, cmap, fieldnorm,
                cbarlevels, clabel, cformat, ovrmap, gis_info, numovr, ovrx, ovry, ovrfields,
                ovrfieldlvls, ovrfieldcolors, axesticks, rasterized=True):
    """Plot a single-paneled figure from model output (possibly staggered grid).
       This version is set up to use ImageGrid with a single panel and is provided
       only to attempt to match the dimensions of the resulting figure with that
       produced by plotsweep_PyArt in radarmodule.py. Clearly need to refactor some
       of this code..."""
    if fig is None:
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
        axes = grid[0]
    if rasterized:
        axes.set_rasterization_zorder(2)
    norm = fieldnorm
    if ptype == 1:  # Contour plot
        plot = axes.contourf(xs, ys, field, levels=clevels, cmap=cmap, norm=norm, zorder=1)
    elif ptype == 2:  # pcolor plot
        # Mask values below lower bounds
        field = np.ma.masked_where(field < clevels[0], field)
        # norm = matplotlib.colors.BoundaryNorm(clevels,cmap.N)
        # This is a hack to force masked areas to be white.  By default masked areas are
        # transparent, but the Agg backends apparently have problems when converting pcolor images
        # to eps and changes the transparent color to black.
        cmap.set_bad('white', alpha=None)
        plot = axes.pcolormesh(x, y, field, vmin=clevels[0], vmax=clevels[-1], cmap=cmap, norm=norm,
                               edgecolors='None', antialiased=False, rasterized=rasterized)
    # Find a nice number of ticks for the colorbar

    if cbarlevels is None:
        cintv = clevels[1] - clevels[0]
        cintvs = np.arange(clevels[0], clevels[-1], cintv)
        while True:
            if cintvs.size > 20:
                cintv = (cintvs[1] - cintvs[0]) * 2.
                cintvs = np.arange(cintvs[0], cintvs[-1], cintv)
            else:
                break
        cbarlevels = ticker.MultipleLocator(base=cintv)

    grid.cbar_axes[0].colorbar(plot)
    grid.cbar_axes[0].toggle_label(True)
    grid.cbar_axes[0].yaxis.set_major_locator(cbarlevels)
    if clabel is not None:
        grid.cbar_axes[0].set_ylabel(clabel)
    if numovr > 0:
        for i in range(numovr):
            #             if(i != 1): # Hack!
            axes.contour(ovrx[i], ovry[i], ovrfields[i], levels=ovrfieldlvls[i],
                         colors=ovrfieldcolors[i], lw=2)
    if gis_info is not None:
        axes.plot([gis_info[1]], [gis_info[2]], 'ko')
    if xlim is None:
        xlim = [x.min(), x.max()]
    if ylim is None:
        ylim = [y.min(), y.max()]
    print(xlim, ylim)
    axes.set_xlim(xlim[0], xlim[1])
    axes.set_ylim(ylim[0], ylim[1])
    if axesticks[0] < 1000.:
        if axesticks[0] < 500.:
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.xaxis.set_major_formatter(formatter)
    if axesticks[1] < 1000.:
        if axesticks[0] < 500.:
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.yaxis.set_major_formatter(formatter)
    print(axesticks[0], axesticks[1])
    axes.xaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[0]))
    axes.yaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[1]))
    axes.set_xlabel('km')
    axes.set_ylabel('km')
    if axesticks[1] == axesticks[0]:
        axes.set_aspect('equal')
    else:
        axes.set_aspect('auto')

#     if(ovrmap): # Overlay map
#         readshapefile(track_shapefile_location,'track',drawbounds=True,linewidth=0.5,color='black',ax=axes)
# readshapefile(county_shapefile_location,'counties',drawbounds=True,
# linewidth=0.5, color='gray',ax=axes)  #Draws US county boundaries.

    return fig, axes, grid


def plot_wind_meteogram(plottimes, conv_plot_ds, global_plot_config_dict, avgwind=True,
                        windavgintv=60, windgustintv=3, xlimits=None, ptype='PIPS'):
    """[summary]

    Parameters
    ----------
    plottimes : [type]
        [description]
    conv_plot_ds : [type]
        [description]
    global_plot_config_dict : [type]
        [description]
    windavgintv : int, optional
        [description], by default 60
    windgustintv : int, optional
        [description], by default 3
    xlimits : [type], optional
        [description], by default None
    ptype : str, optional
        [description], by default 'PIPS'

    Returns
    -------
    [type]
        [description]
    """
    # Plot wind meteogram
    if ptype == 'PIPS':
        winddirstr = 'winddirabs'
        windspdstr = 'windspd'
        windguststr = 'windgust'
    elif ptype == 'CU':
        winddirstr = 'bwinddirabs'
        windspdstr = 'bwindspd'
    elif ptype == 'NV2':
        winddirstr = 'swinddirabs'
        windspdstr = 'swindspd'

    # Note that we are extracting the underlying numpy arrays here. In the future will try to
    # make everything work with xarray dataarrays/datasets
    winddirabs = conv_plot_ds[winddirstr].values
    windspd = conv_plot_ds[windspdstr].values

    # Handle bad wind values
    maskbadwind = global_plot_config_dict.get('maskbadwind', False)
    plot_diagnostics = global_plot_config_dict['plot_diagnostics']

    if (maskbadwind or plot_diagnostics) and ptype == 'PIPS':
        winddiag = conv_plot_ds['winddiag']
        # Extract indices for "bad" wind data
        winddiag_index = np.where(np.any([winddiag > 0, np.isnan(winddiag)], axis=0))[0]
        if maskbadwind:
            windspd = np.where(winddiag > 0, np.nan, windspd)
            winddirabs = np.where(winddiag > 0, np.nan, winddirabs)

    # Compute average wind speed and direction, and wind gusts
    # For the NV2 probes, the data are already at 60-s intervals and gust information is
    # provided, so does not need to be computed. In the future, need to make a more general
    # interface to allow for averaging of NV2 data at longer intervals.
    if ptype == 'NV2':
        windspdavg = windspd
        windspdavgvec = windspd
        windgustavg = conv_plot_ds['swindgust'].values
        # TODO: Check that wind directions from NV2 probes are vector averages
        winddiravgvec = winddirabs
    elif avgwind:
        windspdavg, windspdavgvec, winddiravgvec, windgust, windgustavg = \
            pips.avgwind(winddirabs, windspd, windavgintv, gusts=True, gustintv=windgustintv,
                         center=False)
    else:   # No additional averaging/fields already present in file
        windspdavg = windspd
        try:
            windgustavg = conv_plot_ds[windguststr]
        except KeyError:
            windgustavg = None
        winddiravgvec = conv_plot_ds[winddirstr]

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()

    if avgwind:
        titlestring = ("Wind direction, wind speed ({:d}-s mean), "
                       "and gust ({:d}-s max of {:d}-s mean)").format(windavgintv, windgustintv,
                                                                      windavgintv)
    else:
        titlestring = 'Wind direction and wind speed'
    plt.title(titlestring)

    fields = [windspdavg]
    fieldparamdicts = [windspeed_params]

    if windgustavg is not None:
        fields.append(windgustavg)
        fieldparamdicts.append(windgust_params)

    # Indicate bad wind times with vertical green lines if desired
    if plot_diagnostics and ptype == 'PIPS':
        # These are the times with bad wind data
        winddiag_plot = plottimes[winddiag_index]
        fields.append(winddiag_plot)
        fieldparamdicts.append(winddiag_params)

    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    fields = [winddiravgvec]
    fieldparamdicts = [winddir_params]
    ax2 = plotmeteogram(ax2, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': global_plot_config_dict['majorxlocator'],
                    'majorxformatter': global_plot_config_dict['majorxformatter'],
                    'minorxlocator': global_plot_config_dict['minorxlocator'],
                    'axeslimits': [xlimits, global_plot_config_dict['ws_range']],
                    'axeslabels': [global_plot_config_dict['xlabel'],
                                   r'wind speed (m s$^{-1}$)']}
    axparamdict2 = {'majorylocator': ticker.MultipleLocator(45.),
                    'axeslimits': [None, [0.0, 360.0]],
                    'axeslabels': [None, r'Wind direction ($^{\circ}$)']}
    axparamdicts = [axparamdict1, axparamdict2]
    ax1, ax2 = set_meteogram_axes([ax1, ax2], axparamdicts)

    return fig, ax1, ax2


def plot_compass_dir_meteogram(plottimes, conv_plot_ds, global_plot_config_dict, xlimits=None,
                               ptype='PIPS'):
    """_summary_

    Parameters
    ----------
    plottimes : _type_
        _description_
    conv_plot_ds : _type_
        _description_
    global_plot_config_dict : _type_
        _description_
    xlimits : _type_, optional
        _description_, by default None
    ptype : str, optional
        _description_, by default 'PIPS'
    """

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [conv_plot_ds['compass_dir'].values]
    fieldparamdicts = [compass_dir_params]

    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': global_plot_config_dict['majorxlocator'],
                    'majorxformatter': global_plot_config_dict['majorxformatter'],
                    'minorxlocator': global_plot_config_dict['minorxlocator'],
                    'majorylocator': ticker.MultipleLocator(45.),
                    'axeslimits': [xlimits, [0.0, 360.0]],
                    'axeslabels': [global_plot_config_dict['xlabel'],
                                   r'Compass heading ($^{\circ}$)']}

    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    return fig, ax1


def plot_temperature_dewpoint_meteogram(plottimes, conv_plot_ds, global_plot_config_dict,
                                        xlimits=None, ptype='PIPS'):
    # Plot temperature and dewpoint
    # tavgintv = 10  # Currently not used

    if ptype == 'PIPS':
        tempstr = 'fasttemp'
    elif ptype == 'CU':
        tempstr = 'slowtemp'
    elif ptype == 'NV2':
        tempstr = 'slowtemp'

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [conv_plot_ds[tempstr].values, conv_plot_ds['dewpoint'].values]
    temp_params['plotmin'] = global_plot_config_dict['T_Td_range'][0]
    dewpoint_params['plotmin'] = global_plot_config_dict['T_Td_range'][0]
    fieldparamdicts = [temp_params, dewpoint_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    ax1.axhline(0.0, ls=':', color='k')

    axparamdict1 = {
        'majorxlocator': global_plot_config_dict['majorxlocator'],
        'majorxformatter': global_plot_config_dict['majorxformatter'],
        'minorxlocator': global_plot_config_dict['minorxlocator'],
        'axeslimits': [xlimits, global_plot_config_dict['T_Td_range']],
        'axeslabels': [global_plot_config_dict['xlabel'], r'Temperature ($^{\circ}$C)']
    }
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    return fig, ax1


def plot_RH_meteogram(plottimes, conv_plot_ds, global_plot_config_dict, xlimits=None,
                      ptype='PIPS'):
    # Plot relative humidity
    # avgintv = 10  # Currently not used

    if ptype == 'PIPS':
        RHstr = 'RH_derived'
    elif ptype == 'CU':
        RHstr = 'RH'
    elif ptype == 'NV2':
        RHstr = 'RH'

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [conv_plot_ds[RHstr].values]
    fieldparamdicts = [RH_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {
        'majorxlocator': global_plot_config_dict['majorxlocator'],
        'majorxformatter': global_plot_config_dict['majorxformatter'],
        'minorxlocator': global_plot_config_dict['minorxlocator'],
        'axeslimits': [xlimits, [0., 100.]],
        'axeslabels': [global_plot_config_dict['xlabel'], 'Relative Humidity (%)']
    }

    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    return fig, ax1


def plot_pressure_meteogram(plottimes, conv_plot_ds, global_plot_config_dict,
                            xlimits=None, ptype='PIPS'):

    pmin = np.nanmin(conv_plot_ds['pressure'].values)
    pmax = np.nanmax(conv_plot_ds['pressure'].values)

    # pmean = conv_plot_ds['pressure'].values.mean()
    # avgintv = 1  # Currently not used

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [conv_plot_ds['pressure'].values]
    fieldparamdicts = [pressure_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {
        'majorxlocator': global_plot_config_dict['majorxlocator'],
        'majorxformatter': global_plot_config_dict['majorxformatter'],
        'minorxlocator': global_plot_config_dict['minorxlocator'],
        'axeslimits': [xlimits, [pmin - 2.5, pmax + 2.5]],
        'axeslabels': [global_plot_config_dict['xlabel'], 'Station pressure (hPa)']
    }
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    return fig, ax1


def plot_voltage_meteogram(plottimes, conv_plot_ds, global_plot_config_dict,
                           xlimits=None, ptype='PIPS'):

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [conv_plot_ds['voltage'].values]
    fieldparamdicts = [battery_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {
        'majorxlocator': global_plot_config_dict['majorxlocator'],
        'majorxformatter': global_plot_config_dict['majorxformatter'],
        'minorxlocator': global_plot_config_dict['minorxlocator'],
        'axeslimits': [xlimits, [8.0, 16.0]],
        'axeslabels': [global_plot_config_dict['xlabel'], r'Battery Voltage (V)']
    }
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    return fig, ax1


def plot_GPS_speed_meteogram(plottimes, conv_plot_ds, global_plot_config_dict,
                             xlimits=None, ptype='PIPS'):
    try:
        # GPS-derived speed
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(111)

        fields = [conv_plot_ds['GPS_speed'].values]
        np.set_printoptions(threshold=np.inf)
        fieldparamdicts = [GPS_speed_params]
        ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

        axparamdict1 = {
            'majorxlocator': global_plot_config_dict['majorxlocator'],
            'majorxformatter': global_plot_config_dict['majorxformatter'],
            'minorxlocator': global_plot_config_dict['minorxlocator'],
            'axeslimits': [xlimits, [0.0, 20.0]],
            'axeslabels': [global_plot_config_dict['xlabel'], r'GPS speed (m s$^{-1}$)']}
        axparamdicts = [axparamdict1]
        ax1, = set_meteogram_axes([ax1], axparamdicts)

        return fig, ax1
    except KeyError:
        print("No GPS Speed information in file!")
        return


def plotDSDderivedmeteograms(PIPS_index, pc, ib, **PSDderiveddict):
    """Plots meteograms of the various derived DSD quantities from the PIPS"""
    PSDmidtimes = PSDderiveddict.get('PSDmidtimes')
    PSD_plot_df = PSDderiveddict.get('PSD_plot_df')

    xaxislimits = [PSDmidtimes[0], PSDmidtimes[-1]]
    dis_name = ib.dis_name_list[PIPS_index]

    # Rain rates (intensities)
    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [PSD_plot_df['intensity'].values]
    fieldparamdicts = [rainrate_params]
    ax1 = plotmeteogram(ax1, [PSDmidtimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [0.0, 150.0]],
                    'axeslabels': [pc.timelabel, r'Rain rate (mm hr$^{-1}$)']}
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_rainrate.png', dpi=300)
    plt.close(fig)

    # Reflectivity
    if ib.type[PIPS_index] != 'NV2':
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(111)

        fields = [PSD_plot_df['reflectivity'].values]
        fieldparamdicts = [reflectivity_params]
        ax1 = plotmeteogram(ax1, [PSDmidtimes], fields, fieldparamdicts)

        axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                        'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [0.0, 80.0]],
                        'axeslabels': [pc.timelabel, r'Radar reflectivity (dBZ)']}
        axparamdicts = [axparamdict1]
        ax1, = set_meteogram_axes([ax1], axparamdicts)

        plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_reflectivity.png', dpi=300)
        plt.close(fig)

    # Particle counts
    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(111)

    fields = [PSD_plot_df['pcount'].values]
    fieldparamdicts = [pcount_params]
    if ib.type[PIPS_index] != 'NV2':
        fields.append(PSD_plot_df['pcount2'].values)
        fieldparamdicts.append(pcount2_params)

    ax1 = plotmeteogram(ax1, [PSDmidtimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [0.0, 1000.0]],
                    'axeslabels': [pc.timelabel, r'Particle count']}
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)
    ax1.legend(bbox_to_anchor=(1., 1.), loc='upper right', ncol=1, fancybox=True, shadow=False,
               prop=fontP)
    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_pcounts.png', dpi=300)
    plt.close(fig)

    # Signal amplitude
    if ib.type[PIPS_index] == 'PIPS':
        fig = plt.figure(figsize=(8, 3))
        ax1 = fig.add_subplot(111)

        fields = [PSD_plot_df['amplitude'].values]
        fieldparamdicts = [amplitude_params]
        ax1 = plotmeteogram(ax1, [PSDmidtimes], fields, fieldparamdicts)

        axparamdict1 = {
            'majorxlocator': pc.locator,
            'majorxformatter': pc.formatter,
            'minorxlocator': pc.minorlocator,
            'axeslimits': [xaxislimits, [0.0, 30000.0]],
            'axeslabels': [pc.timelabel, r'Signal amplitude']
        }
        axparamdicts = [axparamdict1]
        ax1, = set_meteogram_axes([ax1], axparamdicts)

        plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_amplitude.png', dpi=300)
        plt.close(fig)


def plotDSDmeteograms(dis_name, image_dir, axparams, disvars, radvars=None, close_fig=True):
    """Plots one or more meteograms of disdrometer number concentrations vs. diameter bins,
       along with one or more derived variables and optionally radar variables for comparison.
       One meteogram is plotted per dualpol variable (i.e. Z,ZDR,KDP,RHV)"""
    diameter_bin_edges = disvars.get('diameter_bin_edges', np.empty((0)))
    PSDstarttimes = disvars.get('PSDstarttimes', np.empty((0)))
    PSDmidtimes = disvars.get('PSDmidtimes', np.empty((0)))
    logND = disvars.get('logND', np.empty((0)))
    if not logND.size or not PSDstarttimes.size or not PSDmidtimes.size:
        print("No DSD info to plot! Quitting!")
        return
    # D_0_dis = disvars.get('D_0', np.empty((0)))
    D_m_dis = disvars.get('D_m', np.empty((0)))
    dBZ_ray_dis = disvars.get('dBZ_ray', np.empty((0)))
    flaggedtimes = disvars.get('flaggedtimes', np.empty((0)))
    hailflag = disvars.get('hailflag', np.empty((0)))
    if radvars is not None:
        # D_0_rad = radvars.get('D_0_rad', np.empty((0)))
        D_m_keys = [k for k, v in radvars.items() if 'D_m' in k]
        radmidtimes = radvars.get('radmidtimes', np.empty((0)))
    else:
        # D_0_rad = np.empty((0))
        D_m_keys = []
        radmidtimes = np.empty((0))

    # Try to find the desired dualpol variables for plotting in the provided dictionary

    dualpol_dis_varnames = []
    dualpol_dis_vars = []
    for key, value in disvars.items():
        if key in ['REF', 'ZDR', 'KDP', 'RHO']:
            dualpol_dis_varnames.append(key)
            dualpol_dis_vars.append(value)

    # If no dualpol variables are in the dictionary, we'll still go ahead and plot the
    # meteogram, with just the number concentrations and possibly median volume diameters, etc.
    if not dualpol_dis_varnames:
        dualpol_dis_varnames.append(None)
        dualpol_dis_vars.append(np.empty((0)))

    # Start the plotting loop
    for dualpol_dis_varname, dualpol_dis_var in zip(dualpol_dis_varnames, dualpol_dis_vars):
        # See if the variable is also provided in the radvars dictionary
        if radvars is not None:
            dualpol_rad_var = radvars.get(dualpol_dis_varname, np.empty((0)))
        else:
            dualpol_rad_var = None

        # Create the figure
        fig = plt.figure(figsize=(8, 3))
        ax1 = host_subplot(111)
        fig.autofmt_xdate()
        if dualpol_dis_var.size:
            ax2 = ax1.twinx()
        # divider = make_axes_locatable(ax1)

        xvals = [PSDstarttimes]
        plotvars = [logND]
        plotparamdict1 = {'type': 'pcolor', 'vlimits': [-1.0, 3.0],
                          'clabel': r'log[N ($m^{-3} mm^{-1}$)]'}
        plotparamdicts = [plotparamdict1]

        # # Median volume diameter
        # if(D_0_dis.size):
        #     xvals.append(PSDmidtimes)
        #     plotvars.append(D_0_dis)
        #     plotparamdict2 = {'type': 'line', 'linestyle': ':', 'color': 'r', 'linewidth': 1.0,
        #                       'label': r'$D_{0,dis} (mm)$'}
        #     plotparamdicts.append(plotparamdict2)

        # # Retrieved median volume diameter
        # if D_0_rad.size:
        #     xvals.append(PSDmidtimes)
        #     plotvars.append(D_0_rad)
        #     plotparamdict2 = {'type': 'line', 'linestyle': ':', 'color': 'purple',
        #                       'linewidth': 1.0, 'label': r'$D_{0,rad} (mm)$'}
        #     plotparamdicts.append(plotparamdict2)

        # Mass-weighted mean diameter
        if D_m_dis.size:
            xvals.append(PSDmidtimes)
            plotvars.append(D_m_dis)
            plotparamdict2 = {'type': 'line', 'linestyle': ':', 'color': 'r', 'linewidth': 1.0,
                              'label': r'$D_{m,dis} (mm)$'}
            plotparamdicts.append(plotparamdict2)

        # Retrieved Mass-weighted mean diameter

        if len(D_m_keys) > 0:
            colors = ['purple', 'navy', 'navy', 'green']
            linestyles = [':', '--', '-.', '-']
            for col, ls, D_m_rad_key in zip(colors, linestyles, D_m_keys):
                plot_key = D_m_rad_key.replace('D_m_rad_', '')
                xvals.append(PSDmidtimes)
                plotvars.append(radvars[D_m_rad_key])

                plotparamdict2 = {'type': 'line', 'linestyle': ls, 'color': col, 'linewidth': 1.0,
                                  'label': r'$D_{{m}}$ {} (mm)'.format(plot_key)}
                plotparamdicts.append(plotparamdict2)

        # Vertical lines for flagged times (such as from wind contamination).
        if flaggedtimes.size:
            xvals.append(PSDmidtimes)
            plotvars.append(flaggedtimes)
            plotparamdict = {'type': 'vertical line', 'linestyle': '-', 'color': 'r',
                             'linewidth': 0.5}
            plotparamdicts.append(plotparamdict)

        # Mark times with hail detected with a vertical purple line
        if hailflag.size:
            xvals.append(PSDmidtimes)
            plotvars.append(hailflag)
            plotparamdict = {'type': 'vertical line', 'linestyle': '-', 'color': 'purple',
                             'linewidth': 0.5}
            plotparamdicts.append(plotparamdict)

        ax1 = plotmeteogram(ax1, xvals, plotvars, plotparamdicts,
                            yvals=[diameter_bin_edges] * len(plotvars))

        axparamdict1 = axparams
        axes = [ax1]
        axparamdicts = [axparamdict1]
        ax1.legend(bbox_to_anchor=(0., 1.), loc='upper left',
                   ncol=1, fancybox=True, shadow=False, prop=fontP)

        # Plot the disdrometer-derived dualpol variable
        if dualpol_dis_varname is not None:
            xvals = [PSDmidtimes]
            plotvars = [dualpol_dis_var]
            if dualpol_dis_varname == 'REF':
                dualpol_dis_varlabel = r'$Z_{dis} (dBZ)$'
                axis_label = 'Reflectivity (dBZ)'
                axis_limits = [0.0, 80.0]
                axis_intv = 10.0
                # Also plot Rayleigh version
                if dBZ_ray_dis.size:
                    xvals.append(PSDmidtimes)
                    plotvars.append(dBZ_ray_dis)
                    dBZ_ray_dis_varlabel = r'$Z_{dis,ray} (dBZ)$'
            elif dualpol_dis_varname == 'ZDR':
                dualpol_dis_varlabel = r'$Z_{DR,dis} (dB)$'
                axis_label = 'Differential Reflectivity (dBZ)'
                axis_limits = [0.0, 6.0]
                axis_intv = 0.5
            elif dualpol_dis_varname == 'KDP':
                dualpol_dis_varlabel = r'$K_{DP,dis} (deg km^{-1})$'
                axis_label = r'Specific Differential Phase (deg km$^{-1}$)'
                axis_limits = [0.0, 12.0]
                axis_intv = 1.0
            elif dualpol_dis_varname == 'RHO':
                dualpol_dis_varlabel = r'$\sigma_{HV,dis}$'
                axis_label = r'Cross-correlation Coefficient'
                axis_limits = [0.0, 1.05]
                axis_intv = 0.1

            plotparamdict1 = {'type': 'line', 'linestyle': '-', 'color': 'r', 'linewidth': 1.5,
                              'label': dualpol_dis_varlabel}
            plotparamdicts = [plotparamdict1]
            if dualpol_dis_varname == 'REF' and dBZ_ray_dis.size:
                plotparamdict = {'type': 'line', 'linestyle': '--', 'color': 'r', 'linewidth': 1.5,
                                 'label': dBZ_ray_dis_varlabel}
                plotparamdicts.append(plotparamdict)

        # Now add the radar dualpol variables if desired
        if dualpol_rad_var is not None:

            if dualpol_dis_varname == 'REF':
                dualpol_rad_varlabel = r'$Z_{rad} (dBZ)$'
            elif dualpol_dis_varname == 'ZDR':
                dualpol_rad_varlabel = r'$Z_{DR,rad} (dB)$'
            elif dualpol_dis_varname == 'KDP':
                dualpol_rad_varlabel = r'$K_{DP,rad} (deg km^{-1})$'
            elif dualpol_dis_varname == 'RHO':
                dualpol_rad_varlabel = r'$\sigma_{HV,rad}$'

            if dualpol_rad_var.size:
                xvals.append(radmidtimes)
                plotvars.append(dualpol_rad_var)
                plotparamdict2 = {'type': 'line', 'linestyle': '-', 'color': 'purple',
                                  'linewidth': 1.5, 'label': dualpol_rad_varlabel}
                plotparamdicts.append(plotparamdict2)

        # print xvals,plotvars,plotparamdicts
        if dualpol_dis_varname is not None:
            ax2 = plotmeteogram(ax2, xvals, plotvars, plotparamdicts)
            axparamdict2 = {'majorylocator': ticker.MultipleLocator(base=axis_intv),
                            'axeslimits': [axparams['axeslimits'][0], axis_limits],
                            'axeslabels': [None, axis_label]}
            axes.append(ax2)
            axparamdicts.append(axparamdict2)
            ax2.legend(bbox_to_anchor=(1., 1.), loc='upper right',
                       ncol=1, fancybox=True, shadow=False, prop=fontP)

        axes = set_meteogram_axes(axes, axparamdicts)
        if dualpol_dis_varname is not None:
            figpath = os.path.join(image_dir, dis_name + '_' + dualpol_dis_varname + '_logNc.png')
            plt.savefig(figpath, dpi=300)
        else:
            figpath = os.path.join(image_dir, dis_name + '_logNc.png')
            plt.savefig(figpath, dpi=300)
        if close_fig:
            plt.close(fig)


def plotmeteogram(ax, xvals, zvals, plotparamdicts, yvals=None, plot_data_bounds=True):
    """Plots a meteogram (time series) of one or more meteorological variables"""
    ax = ax or plt.figure().add_subplot(111)
    for idx, xval, zval, plotparamdict in zip(range(len(zvals)), cycle(xvals), zvals,
                                              plotparamdicts):
        mtype = plotparamdict.get('type', 'fill_between')
        linestyle = plotparamdict.get('linestyle', '-')
        linewidth = plotparamdict.get('linewidth', 1.0)
        marker = plotparamdict.get('marker', None)
        color = plotparamdict.get('color', None)
        alpha = plotparamdict.get('alpha', 0.5)
        markeredgecolor = plotparamdict.get('markeredgecolor', 'none')
        ms = plotparamdict.get('ms', 0)
        plotmin = plotparamdict.get('plotmin', 0)
        plotlabel = plotparamdict.get('label', None)
        vlimits = plotparamdict.get('vlimits', None)
        clabel = plotparamdict.get('clabel', None)
        plotcbar = plotparamdict.get('plotcbar', True)

        if mtype == 'fill_between':
            ax.plot_date(xval, zval, ls=linestyle, lw=linewidth, marker=marker, color=color,
                         markeredgecolor=markeredgecolor, ms=ms, label=plotlabel, fmt="")
            ax.fill_between(xval, zval, plotmin, facecolor=color, alpha=alpha)
        elif mtype == 'pcolor':
            divider = make_axes_locatable(ax)
            C = ax.pcolormesh(xval, yvals[idx], zval, vmin=vlimits[0], vmax=vlimits[1])
            cax = divider.append_axes("bottom", size="5%", pad=0.4)
            # Go ahead and create the appended axes, but simply shut it off if we don't
            # want it. This is so that if we are plotting subplots and only want one colorbar
            # between them, they all end up the same size.
            if not plotcbar:
                cax.axis('off')
            else:
                cb = ax.get_figure().colorbar(C, cax=cax, orientation='horizontal')
                if clabel:
                    cb.set_label(clabel)
        # For flagging times with bad data, etc.
        # zval is interpreted as a list of x-indices
        elif mtype == 'vertical line':
            for x in zval:
                ax.axvline(x=x, ls=linestyle, lw=linewidth, color=color)
        else:
            ax.plot_date(xval, zval, ls=linestyle, lw=linewidth, marker=marker, color=color,
                         markeredgecolor=markeredgecolor, ms=ms, label=plotlabel, fmt="")

    # If desired, plot the start and end times of data collection as vertical black dashed lines
    if plot_data_bounds:
        ax.axvline(x=xvals[0][0], ls='--', lw=0.5, color='k', alpha=0.75)
        ax.axvline(x=xvals[0][-1], ls='--', lw=0.5, color='k', alpha=0.75)

    return ax


def set_meteogram_axes(axes, axparamdicts):
    """Sets up and formats the axes for a meteogram plot"""
    for ax, axparamdict in zip(axes, axparamdicts):
        ax = ax or plt.figure().add_subplot(111)
        majorxlocator = axparamdict.get('majorxlocator', None)
        majorxformatter = axparamdict.get('majorxformatter', None)
        majorylocator = axparamdict.get('majorylocator', None)
        # majoryformatter = axparamdict.get('majoryformatter', None)
        minorxlocator = axparamdict.get('minorxlocator', None)
        axeslimits = axparamdict.get('axeslimits', [None, None])
        axeslabels = axparamdict.get('axeslabels', [None, None])
        axesautofmt = axparamdict.get('axesautofmt', True)

        if majorxlocator:
            ax.xaxis.set_major_locator(majorxlocator)
        if majorxformatter:
            ax.xaxis.set_major_formatter(majorxformatter)
        if minorxlocator:
            ax.xaxis.set_minor_locator(minorxlocator)
        if axeslimits[0]:
            ax.set_xlim(axeslimits[0][0], axeslimits[0][1])
        if axeslabels[0]:
            ax.set_xlabel(axeslabels[0])
        if axeslabels[1]:
            ax.set_ylabel(axeslabels[1])
        if axeslimits[1]:
            ax.set_ylim(axeslimits[1][0], axeslimits[1][1])
        if majorylocator:
            ax.yaxis.set_major_locator(majorylocator)

    if axes[0] and axesautofmt:
        axes[0].get_figure().autofmt_xdate()

    return axes


def plot_DSD(axdict, PSDdict, PSDfitdict, PSDparamdict):
    """Plots an individual measured PSD on a semilog plot, along with optional exponential and gamma
       fits and associated parameters."""

    time_to_plot = axdict.get('time', None)
    if time_to_plot is not None:
        time_to_plot_datetime = pd.to_datetime(time_to_plot).to_pydatetime()
    xbin_left = axdict.get('xbin_left', np.empty((0)))
    xbin_right = axdict.get('xbin_right', np.empty((0)))
    xbin_mid = axdict.get('xbin_mid', np.empty((0)))
    ND = PSDdict.get('ND', np.empty((0)))
    # FIXME
    ND_onedrop = PSDdict.get('ND_onedrop', np.empty((0)))
    interval = axdict.get('interval', 10)

    fig1 = plt.figure(figsize=(8, 6))
    ax1 = fig1.add_subplot(111)
    # TODO:
    if time_to_plot is not None:
        titlestring = '{0:d}-s DSD fits for time {1} EST'
        plt.title(titlestring.format(interval, time_to_plot_datetime.strftime(tm.timefmt2)))
    else:
        titlestring = 'Full deployment DSD ({:d} s)'
        plt.title(titlestring.format(interval))

    # print('ND', ND)
    ax1.bar(xbin_left, ND * 1000.0, xbin_right - xbin_left, 10.**2., align='edge', log=True,
            color='tan', edgecolor='k')
    # print('ND_onedrop', ND_onedrop)
    ax1.bar(xbin_left, ND_onedrop * 1000.0, xbin_right - xbin_left, 10.**2., align='edge',
            log=True, fill=False, edgecolor='k')

    # ax1.plot(xbin_mid, ND * 1000.0, lw=2, label='obs')

    for fitname, ND_tuple in PSDfitdict.items():
        if (fitname == 'Dis Retr') | (fitname == 'Rad Retr'):
            ND_fit = ND_tuple[0] * 1000.
        else:
            ND_fit = ND_tuple[0]
#         print ND_fit
        label = ND_tuple[1]
        if ND_fit.size:
            ax1.plot(xbin_mid, ND_fit, lw=2, label=label)

    ax1.set_yscale('log')
    ax1.set_ylim(10.**2.0, 10.**8.5)
    ax1.set_xlim(0.0, 9.0)
    ax1.set_xlabel('D (mm)')
    ax1.set_ylabel(r'N(D) $(m^{-4})$')

    # FIXME
    ypos = 0.95
    for paramname, paramtuple in PSDparamdict.items():
        ax1.text(0.50, ypos, paramtuple[1] + ' = {:2.2f}'.format(float(paramtuple[0])),
                 transform=ax1.transAxes)
        ypos = ypos - 0.05

    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=8)

    return fig1, ax1

    # plt.savefig(ib.image_dir + 'DSDs/' + dis_name + '/' + dis_name + '_DSD_t{:04d}.png'.format(t),
    #             dpi=200, bbox_inches='tight')
    # plt.close(fig1)


def plot_vel_D(axdict, PSDdict, rho):
    """Plots the terminal velocity vs. diameter matrix for a given DSD"""

    time = axdict.get('time', None)
    xlim = axdict.get('xlim', (0.0, 9.0))
    ylim = axdict.get('ylim', (0.0, 15.0))
    diameter_bin_edges = axdict.get('diameter_bin_edges', None)
    fallspeed_bin_edges = axdict.get('fallspeed_bin_edges', None)
    avg_diameter = axdict.get('avg_diameter', None)
    cblim = axdict.get('cblim', (1, 50))

    vd_matrix_da = PSDdict.get('vd_matrix_da', None)
    flaggedtime = PSDdict.get('flaggedtime', 0)
    DSD_interval = PSDdict.get('DSD_interval', 10.)

    fig1 = plt.figure(figsize=(8, 6))
    ax1 = fig1.add_subplot(111)
    if time:
        titlestring = 'Fall speed vs. diameter for time {} and interval {:d} s'
        plt.title(titlestring.format(time.strftime(tm.timefmt2), int(DSD_interval)))
    else:
        titlestring = 'Fall speed vs. diameter for full deployment ({:d} s)'
        plt.title(titlestring.format(int(DSD_interval)))

    countsplot = vd_matrix_da
    # countsplot = np.ma.masked_where(vd_matrix_da <= 0, vd_matrix_da)
    C = ax1.pcolormesh(diameter_bin_edges, fallspeed_bin_edges, countsplot, vmin=cblim[0],
                       vmax=cblim[1], edgecolors='w', linewidths=0.2, cmap=cm.plasma)
    rainvd = pips.calc_empirical_fallspeed(avg_diameter, correct_rho=True, rho=rho)

    ax1.plot(avg_diameter, rainvd, c='r')
    # ax1.scatter(X[0:10,20:31],Y[0:10,20:31],c='r',marker='x')
    fig1.colorbar(C)

    # FIXME
    if flaggedtime > 1:
        ax1.text(0.5, 0.5, 'Flagged for strong wind contamination!',
                 horizontalalignment='center',
                 verticalalignment='center', color='y',
                 transform=ax1.transAxes)
    # FIXME
    # if(dis.plot_strongwindQC):
    #     ax1.scatter(X[dis.strongwindmask], Y[dis.strongwindmask], c='r', marker='x', alpha=1.0)
    # if(dis.plot_splashingQC):
    #     ax1.scatter(X[dis.splashmask], Y[dis.splashmask], c='w', marker='o', alpha=0.75)
    #     # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(splashmask,mask=-splashmask),
    #                  cmap=cm.Reds,alpha=0.1)
    # if(dis.plot_marginQC):
    #     ax1.scatter(X[dis.marginmask], Y[dis.marginmask], c='g', marker='x', alpha=0.1)
    #     # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(marginmask,mask=-marginmask),
    #                  cmap=cm.Reds,alpha=0.1)
    # if(dis.plot_rainfallspeedQC):
    #     ax1.scatter(X[dis.fallspeedmask], Y[dis.fallspeedmask], c='k', marker='x', alpha=0.5)
    #     # ax1.pcolor(min_diameter,min_fall_bins,
    #                  ma.masked_array(fallspeedmask,mask=-fallspeedmask),cmap=cm.gray,alpha=0.1)
    # if(dis.plot_rainonlyQC):
    #     ax1.scatter(X[dis.rainonlymask], Y[dis.rainonlymask], c='g', marker='x', alpha=0.5)

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax1.set_xlabel('diameter (mm)')
    ax1.set_ylim(ylim[0], ylim[1])
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax1.set_ylabel('fall speed (m/s)')

    return fig1, ax1


def computecorners(xe, ye, UM=False):
    """
    Given 2D meshes of the coordinates of the grid edges (i.e. xe(Nys,Nxe), ye(Nye,Nxs)
    where Nxs,Nxe,Nys,Nye are the number of grid points in each direction for the scalar (s)
    and staggered (e) grid), compute the grid corners as simple horizontal averages of adjacent
    edge points in each direction.  For the edges of the domain, an optional simple linear
    extrapolation is performed. Returns xcor and ycor of dimensions (Nye,Nxe).
    """

    if not UM:
        cor_shp = xe.shape[:-2] + (ye.shape[-2], xe.shape[-1])
        xcor = np.zeros(cor_shp)
        ycor = np.zeros(cor_shp)

#       print xcor.shape, ycor.shape, xe.shape, ye.shape

        xcor[..., 1:-1, :] = 0.5 * (xe[..., 1:, :] + xe[..., :-1, :]
                                    )  # All but south and north edges
        ycor[..., :, 1:-1] = 0.5 * (ye[..., :, 1:] + ye[..., :, :-1]
                                    )  # All but west and east edges
        xcor[..., 0, :] = xe[..., 2, :] + 1.5 * \
            (xe[..., 1, :] - xe[..., 2, :])  # Extrapolate for south edge
        xcor[..., -1, :] = xe[..., -2, :] + 1.5 * \
            (xe[..., -1, :] - xe[..., -2, :])  # Extrapolate for north edge
        ycor[..., :, 0] = ye[..., :, 2] + 1.5 * \
            (ye[..., :, 1] - ye[..., :, 2])  # Extrapolate for west edge
        ycor[..., :, -1] = ye[..., :, -2] + 1.5 * \
            (ye[..., :, -1] - ye[..., :, -2])  # Extrapolate for east edge
    else:   # Case for UM grid
        xcor = 0.5 * (xe[..., 1:, :] + xe[..., :-1, :])
        ycor = 0.5 * (ye[..., :, 1:] + ye[..., :, :-1])

    return xcor, ycor


# From https://gist.github.com/syrte/592a062c562cd2a98a83
def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, )
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)
        `c` can be a 2-D array in which the rows are RGB or RGBA, however.
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls),
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection
# Below is some extra code for DSD plotting with no home right now

    #
#         #print "logND",logND
#         # Times are valid at end of DSD intervals
#         C = ax1.pcolor(plotx_dis_start[pstartindex:pstopindex+1],min_diameter,
#                        logND[:,pstartindex:pstopindex+1],vmin=-1.0,vmax=3.0)
#
#         if(calc_DSD):
#             #Overlay mass-weighted mean diameter
#             #ax1.plot(plotx_dis_avg,D_m_disd,ls=':',c='k',lw=2,label=r'$D_{mr43} (mm)$')
#             #Overlay median volume diameter (centered running mean)
#             D_med_disd_avg = pd.rolling_mean(D_med_disd,12,center=True,win_type='triang',
#                                              min_periods=1)
#             ax1.plot(plotx_dis_avg,D_med_disd_avg,ls=':',c='b',lw=0.5,label=r'$D_{0} (mm)$')
#             #ax1.plot(plotx_dis_avg,D_med_gam,ls='--',c='k',lw=0.5,label=r'$D_{0,gam} (mm)$')
#
#         ax1.xaxis.set_major_locator(locator)
#         ax1.xaxis.set_major_formatter(DateFormatter(dateformat))
#         ax1.xaxis.set_minor_locator(minorlocator)
#         ax1.set_xlim(starttime,stoptime)
#
#         ax1.yaxis.set_major_locator(MultipleLocator(base=1.0))
#         ax1.set_ylim(0.0,9.0)
#         ax1.set_ylabel('D (mm)')
#         #ax1.legend(loc='upper right')
#         cax = divider.append_axes("bottom", size="5%", pad=0.35)
#         cb = fig2.colorbar(C, cax=cax,orientation='horizontal')
#         cb.set_label(r'log[N ($m^{-3} mm^{-1}$)]')
#
#         #ax2.set_aspect(5000.0)
#         #ax2.set_aspect(0.0001)
#
#         #divider.set_aspect(True)
#
#         if(calc_DSD):
#             # Average reflectivity from disdrometer using a centered running mean
#             refl_disd_avg = pd.rolling_mean(refl_disd,12,center=True,win_type='triang',
#                                             min_periods=1)
#             #ax2.plot(plotx_dis_avg,refl_disd,ls='-',c='r',marker='o',label=r'$Z_{dis} (dBZ)$')
#             ax2.plot(plotx_dis_avg,refl_disd_avg,ls='-',c='green',lw=0.5,marker=None,
#                      label=r'$Z_{dis} (dBZ)$')
#             if(comp_radar):
#                 # Assumes reflectivity is the first.  Need to come back and work on this.
#                 dBZ_D_plt = fields_D_tarr[:,index,0]
#                 #ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',marker='o',label=r'$Z_{rad} (dBZ)$')
#                 ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',lw=0.5,marker=None,
#                          label=r'$Z_{rad} (dBZ)$')
# #             if(not timetospace):
# #                 ax2.xaxis.set_major_locator(locator)
# #                 ax2.xaxis.set_minor_locator(minorlocator)
# #                 ax2.xaxis.set_major_formatter(DateFormatter(dateformat))
# #                 ax2.set_xlim(starttime,stoptime)
# #                 if(reverse_times):
# #                     ax2.invert_xaxis()
# #             elif(comp_radar):
# #                 ax2.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
# #                 ax2.set_xlim(xstart,xstop)
# #                 ax2.invert_xaxis()
#             ax2.yaxis.set_major_locator(MultipleLocator(base=5.0))
#             ax2.set_ylim(-10.0,80.0)
#             ax2.set_ylabel('Reflectivity (dBZ)')
#             ax2.legend(loc='best')
# #             pcountax = divider.append_axes("top",size="50%",pad=0.2,sharex=ax1)
# #             pcountax.plot(plotx_dis_avg,pcount,ls='-',marker='o',c='k')
# #             pcountax.plot(plotx_dis_avg,pcount2,ls='-',marker='o',c='b')
# #             if(not timetospace):
# #                 pcountax.xaxis.set_major_locator(MinuteLocator(interval=5))
# #                 pcountax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
# #                 pcountax.set_xlim(starttime,endtime)
# #                 if(reverse_times):
# #                     pcountax.invert_xaxis()
# #             else:
# #                 pcountax.xaxis.set_major_locator(MultipleLocator(base=xtickintv))
# #                 pcountax.set_xlim(xstart,xstop)
# #                 pcountax.invert_xaxis()
# #             pcountax.yaxis.set_major_locator(MultipleLocator(base=5.0))
# #             pcountax.set_yscale('log')
# #             pcountax.set_ylim(1.0,6000.0)
# #             pcountax.set_ylabel('# of particles')
# #             pcountax.set_title('Disdrometer Log(N) and Reflectivity')
# #             plt.xticks(visible=False)


def plot_mu_lamda(poly_coeff=None, poly=None, lamda=None, mu=None, ax=None, title=None,
                  plot_Z01_C08=True, plot_legend=True, plot_poly_coeff=True):

    xx = np.linspace(0.0, 30.0)
    if plot_Z01_C08:
        y_Cao = -0.0201 * xx**2. + 0.902 * xx - 1.718
        y_Zhang = -0.016 * xx**2. + 1.213 * xx - 1.957
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = plt.gcf()
    if title:
        plt.title(title)
    if lamda is not None:
        ax.scatter(lamda, mu, color='k', marker='.')
        if poly is not None:
            yy = poly(xx)
            ax.plot(xx, yy, label='Present study')
    if plot_Z01_C08:
        ax.plot(xx, y_Cao, label='C08')
        ax.plot(xx, y_Zhang, label='Z01')
    ax.set_xlim(0.0, 20.0)
    ax.set_ylim(-5.0, 32.0)
    ax.set_xlabel(r'$\lambda$ (mm$^{-1}$)')
    ax.set_ylabel(r'$\mu$')
    if lamda is not None:
        ax.text(0.05, 0.85, '# of Points: {:d}'.format(len(lamda)),
                transform=ax.transAxes, fontsize=12.)
    if plot_poly_coeff and poly_coeff is not None:
        op1 = '+' if np.sign(poly_coeff[1]) == 1 else '-'
        op2 = '+' if np.sign(poly_coeff[0]) == 1 else '-'
        polytext = r'$\mu = {0:2.4f}\lambda^{{2}} {1} {2:2.4f}\lambda {3} {4:2.4f}$'

        polytext = polytext.format(poly_coeff[2], op1, np.abs(poly_coeff[1]), op2,
                                   np.abs(poly_coeff[0]))
        ax.text(0.05, 0.80, polytext, transform=ax.transAxes, fontsize=12.)
    if plot_legend:
        plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=12.)

    return fig, ax


def plot_scatter(ds, var_x, var_y, axparams, fig=None, ax=None, add_colorbar=True):
    label_x = axparams.get('label_x', var_x)
    label_y = axparams.get('label_y', var_y)
    var_lims = axparams.get('var_lims', [[0., 1.], [0., 1.]])
    plot_log = axparams.get('plot_log', [False, False])
    col_field = axparams.get('col_field', None)
    label_cb = axparams.get('label_cb', col_field)
    norm = axparams.get('norm', None)
    col_field_lims = axparams.get('col_field_lims', [0., 1.])
    color = axparams.get('color', 'k')
    alpha = axparams.get('alpha', 1.)
    markersize = axparams.get('markersize', 10)
    markerstyle = axparams.get('markerstyle', 'o')

    if plot_log[0]:
        x_lims = [10.**v for v in var_lims[0]]
    else:
        x_lims = var_lims[0]
    if plot_log[1]:
        y_lims = [10.**v for v in var_lims[1]]
    else:
        y_lims = var_lims[1]

    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif not fig:
        fig = plt.gcf()

    if col_field:
        xrplot.scatter(ds, var_x, var_y, ax=ax, hue=col_field, hue_style='continuous', norm=norm,
                       vmin=col_field_lims[0], vmax=col_field_lims[1], alpha=alpha, s=markersize,
                       marker=markerstyle, cbar_kwargs={'label': label_cb},
                       add_guide=add_colorbar)
    else:
        xrplot.scatter(ds, var_x, var_y, ax=ax, colors=color)

    ax.set_xlabel(label_x)
    ax.set_ylabel(label_y)
    if plot_log[0]:
        ax.set_xscale('log')
    if plot_log[1]:
        ax.set_yscale('log')
    ax.set_xlim(x_lims)
    ax.set_ylim(y_lims)
    ax.set_aspect('equal')

    return fig, ax


def plot_one2one(ds, var_x, var_y, axparams, fig=None, ax=None, compute_stats=True,
                 add_colorbar=True):

    var_lims = axparams.get('var_lims', [[0., 1.], [0., 1.]])
    plot_log = axparams.get('plot_log', [False, False])
    if compute_stats:
        stat_text_loc = axparams.get('stat_text_loc', [(0.1, 0.9), (0.1, 0.85)])
        stat_labels = axparams.get('stat_labels', [r'$\rho_{{rd}}$: {:.3f}',
                                                   r'Bias$_{{rd}}$: {:.3f}'])

    fig, ax = plot_scatter(ds, var_x, var_y, axparams, fig=fig, ax=ax, add_colorbar=add_colorbar)
    if plot_log[0]:
        x_lims = [10.**v for v in var_lims[0]]
    else:
        x_lims = var_lims[0]
    if plot_log[1]:
        y_lims = [10.**v for v in var_lims[1]]
    else:
        y_lims = var_lims[1]
    ax.plot(x_lims, y_lims, color='k')

    if compute_stats:
        cc, bias = pips.calc_stats(ds, var_x, var_y)
        ax.text(*stat_text_loc[0], stat_labels[0].format(cc.iloc[0, 1]), transform=ax.transAxes)
        ax.text(*stat_text_loc[1], stat_labels[1].format(bias), transform=ax.transAxes)

    return fig, ax


def plot_retr_timeseries(obs_dict, retr_dis_dict, retr_rad_dict, DSDmidtimes, axparams, fig=None,
                         ax=None, compute_stats=True, name=None):

    obs = obs_dict.get('field', None)
    retr_dis = retr_dis_dict.get('field', None)
    retr_rad = retr_rad_dict.get('field', None)

    obs_param_dict = obs_dict.get('plotparams', {
        'linestyle': '-',
        'color': 'k',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'observed'
    })

    retr_dis_param_dict = retr_dis_dict.get('plotparams', {
        'linestyle': '-',
        'color': 'c',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'dis retrieved'
    })

    retr_rad_param_dict = retr_rad_dict.get('plotparams', {
        'linestyle': '-',
        'color': 'g',
        'alpha': 0.25,
        'plotmin': 0,
        'label': r'rad retrieved'
    })

    if compute_stats:
        bias_dis = (100. * (retr_dis - obs).mean() / obs.mean()).values
        bias_rad = (100. * (retr_rad - obs).mean() / obs.mean()).values
        cc_dis = pd.DataFrame({'x': obs, 'y': retr_dis}).corr()
        cc_rad = pd.DataFrame({'x': obs, 'y': retr_rad}).corr()

    if not ax:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif not fig:
        fig = plt.gcf()

    fields = [obs, retr_dis, retr_rad]
    fieldparamdicts = [obs_param_dict, retr_dis_param_dict, retr_rad_param_dict]
    xvals = [DSDmidtimes] * len(fields)
    ax = plotmeteogram(ax, xvals, fields, fieldparamdicts)
    axparamdicts = [axparams]
    ax, = set_meteogram_axes([ax], axparamdicts)
    if (name == 'Nt' or name == 'R'):
        ax.set_yscale('log')
    ax.text(0.05, 0.93, 'Dis Retr. Bias =%2.2f' % bias_dis + '%', transform=ax.transAxes)
    ax.text(0.05, 0.86, 'Rad Retr. Bias =%2.2f' % bias_rad + '%', transform=ax.transAxes)
    ax.text(0.05, 0.79, 'Dis Retr. Corr Coeff =%2.3f' %
            cc_dis.iloc[0, 1], transform=ax.transAxes)
    ax.text(0.05, 0.72, 'Rad Retr. Corr Coeff =%2.3f' %
            cc_rad.iloc[0, 1], transform=ax.transAxes)
    ax.legend(
        bbox_to_anchor=(
            1.,
            1.),
        loc='upper right',
        ncol=1,
        fancybox=True,
        shadow=False,
        prop=fontP)
    # plt.savefig(image_dir + 'meteograms/' + dis_name + '_' + name + '.png', dpi=300)
    # plt.close(fig)
    return fig, ax


def plot_animation(xplt, yplt, field_da, clevels, cbarlabel=None, cbarintv=None,
                   cmap='pyart_HomeyerRainbow', norm=None, PIPS_list=None, PIPS_xy_list=None,
                   ax=None, ptype='pcolor', axestickintv=10000., axeslimits=None):

    if norm is None:
        norm = cm.colors.Normalize(vmin=clevels[0], vmax=clevels[-1])
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    else:
        fig = ax.get_figure()
    ims = []
    for i, var in enumerate(field_da):
        plotdata = []
        time = np.datetime_as_string(var.coords['time'].values, unit='m')  # Ugly, but whatever

        title = ax.text(0.5, 1.05, f"Time: {time}",
                        size=plt.rcParams["axes.titlesize"],
                        ha="center", transform=ax.transAxes)
        plotdata.append(title)

        if ptype == 'pcolor':
            ci = ax.pcolormesh(xplt, yplt, var.squeeze(), vmin=clevels[0], vmax=clevels[-1],
                               cmap=cmap, norm=norm)
            plotdata.append(ci)
        else:
            ci = ax.contourf(xplt, yplt, var.squeeze(), levels=clevels,
                             cmap=cmap, norm=norm)
            plotdata.extend(ci.collections)

        if PIPS_list is not None and PIPS_xy_list is not None:
            # Plot PIPS locations
            for PIPS, PIPS_xy in zip(PIPS_list, PIPS_xy_list):
                PIPS_x = PIPS_xy[0]
                PIPS_y = PIPS_xy[1]
                ax.plot([PIPS_x], [PIPS_y], 'k*')
        if i == 0.:
            if cbarintv is None:
                cbarintv = clevels[1] - clevels[0]
            cbarlevels = ticker.MultipleLocator(base=cbarintv)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(ci, orientation='vertical', ticks=cbarlevels, cax=cax)
            if cbarlabel is not None:
                cax.set_ylabel(cbarlabel)
            formatter = ticker.FuncFormatter(mtokm)
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(base=axestickintv))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(base=axestickintv))
            ax.set_xlabel('km')
            ax.set_ylabel('km')
            if axeslimits is None:
                xmin = xplt[0]
                xmax = xplt[-1]
                ymin = yplt[0]
                ymax = yplt[-1]
            else:
                xmin, xmax, ymin, ymax = axeslimits
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect('equal')

        ims.append(plotdata)

    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    plt.close()
    return ani