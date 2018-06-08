# plotmodule.py: A module containing some functions related to plotting model output
import numpy as N
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable, host_subplot
import matplotlib.ticker as ticker
from matplotlib.collections import Collection, LineCollection
from matplotlib.artist import allow_rasterization
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.transforms as mtransforms
from matplotlib.font_manager import FontProperties
import ctablesfrompyesviewer as ctables
import disdrometer_module as dis
import timemodule as tm
from itertools import cycle

# Set global font size for axes and colorbar labels, etc.

font = {'size': 10}
matplotlib.rc('font', **font)

fontP = FontProperties()
fontP.set_size('small')

# Contour levels for reflectivity (dBZ)
clevels_ref = N.arange(5.0, 85.0, 5.0)
clevels_zdr = N.arange(0.0, 6.25, 0.25)         # Contour levels for Zdr (dB)
clevels_vr = N.arange(-40.0, 41.0, 1.0)        # Contour levels for Vr (m/s)
cmapdBZ = ctables.__getattribute__('REF_default')
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

    points = N.array([x, y]).T.reshape(-1, 1, 2)
    segments = N.concatenate([points[:-1], points[1:]], axis=1)

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
        z = N.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = N.array([z])

    z = N.asarray(z)

    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    if(ax is None):
        ax = plt.gca()
    ax.add_collection(lc)

    return lc


def clear_frame(ax=None):
    # Taken from a post by Tony S Yu
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.itervalues():
        spine.set_visible(False)


def plotsingle(fig, axes, ptype, xs, ys, x, y, xlim, ylim, field, clevels, cmap, fieldnorm,
               cbarlevels, clabel, cformat, ovrmap, gis_info, numovr, ovrx, ovry, ovrfields,
               ovrfieldlvls, ovrfieldcolors, axesticks, rasterized=True):
    """Plot a single-paneled figure from model output (possibly staggered grid)"""
    if(fig is None):
        fig = plt.figure()
    if(axes is None):
        axes = fig.add_subplot(111)
    if(rasterized):
        axes.set_rasterization_zorder(2)
    norm = fieldnorm
    if(ptype == 1):  # Contour plot
        plot = axes.contourf(xs, ys, field, levels=clevels, cmap=cmap, norm=norm, zorder=1)
    elif(ptype == 2):  # pcolor plot
        # Mask values below lower bounds
        field = N.ma.masked_where(field < clevels[0], field)
        # norm = matplotlib.colors.BoundaryNorm(clevels,cmap.N)
        # This is a hack to force masked areas to be white.  By default masked areas are
        # transparent, but the Agg backends apparently have problems when converting pcolor images
        # to eps and changes the transparent color to black.
        cmap.set_bad('white', alpha=None)
        plot = axes.pcolormesh(x, y, field, vmin=clevels[0], vmax=clevels[-1], cmap=cmap, norm=norm,
                               edgecolors='None', antialiased=False, rasterized=rasterized)
    # Find a nice number of ticks for the colorbar

    if(cbarlevels is None):
        cintv = clevels[1] - clevels[0]
        cintvs = N.arange(clevels[0], clevels[-1], cintv)
        while True:
            if(cintvs.size > 20):
                cintv = (cintvs[1] - cintvs[0]) * 2.
                cintvs = N.arange(cintvs[0], cintvs[-1], cintv)
            else:
                break
        cbarlevels = ticker.MultipleLocator(base=cintv)

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(plot, orientation='vertical', ticks=cbarlevels, cax=cax, format=cformat)
    if(clabel is not None):
        cax.set_ylabel(clabel)
    if(numovr > 0):
        for i in xrange(numovr):
            #             if(i != 1): # Hack!
            plotovr = axes.contour(ovrx[i], ovry[i], ovrfields[i], levels=ovrfieldlvls[i],
                                   colors=ovrfieldcolors[i], lw=2)
    if(gis_info is not None):
        axes.plot([gis_info[1]], [gis_info[2]], 'ko')
    axes.set_xlim(xlim[0], xlim[1])
    axes.set_ylim(ylim[0], ylim[1])
    if(axesticks[0] < 1000.):
        if(axesticks[0] < 500.):
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.xaxis.set_major_formatter(formatter)
    if(axesticks[1] < 1000.):
        if(axesticks[0] < 500.):
            formatter = ticker.FuncFormatter(mtokmr2)
        else:
            formatter = ticker.FuncFormatter(mtokmr1)
    else:
        formatter = ticker.FuncFormatter(mtokm)
    axes.yaxis.set_major_formatter(formatter)
    print axesticks[0], axesticks[1]
    axes.xaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[0]))
    axes.yaxis.set_major_locator(ticker.MultipleLocator(base=axesticks[1]))
    axes.set_xlabel('km')
    axes.set_ylabel('km')
    if(axesticks[1] == axesticks[0]):
        axes.set_aspect('equal')
    else:
        axes.set_aspect('auto')

#     if(ovrmap): # Overlay map
#         readshapefile(track_shapefile_location,'track',drawbounds=True,linewidth=0.5,color='black',ax=axes)
# readshapefile(county_shapefile_location,'counties',drawbounds=True,
# linewidth=0.5, color='gray',ax=axes)  #Draws US county boundaries.

    return fig, axes


def plotonesecmeteograms(dis_index, pc, ib, convmeteodict):
    """Plots meteograms of the one-second PIPS data"""
    plottimes = convmeteodict.get('plottimes')
    onesec_plot_df = convmeteodict.get('onesec_plot_df')
    windavgintv = convmeteodict.get('windavgintv')
    windgustintv = convmeteodict.get('windgustintv')

    xaxislimits = [plottimes[0], plottimes[-1]]
    dis_name = ib.dis_name_list[dis_index]

    # Plot wind meteogram

    winddirabs = onesec_plot_df['winddirabs'].values
    windspd = onesec_plot_df['windspd'].values

    # Compute wind speed and direction, and wind gusts
    windspdavg, windspdavgvec, winddiravgvec, windgust, windgustavg = dis.avgwind(
        winddirabs, windspd, windavgintv, gusts=True, gustintv=windgustintv, center=False)

    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    # plt.title('Wind speed (5-min mean) and gust (5-min max of 3-s mean)')

    fields = [windspdavg, windgustavg]
    fieldparamdicts = [windspeed_params, windgust_params]

    # Add vertical lines to indicate bad wind values if desired
    if(pc.plot_diagnostics):
        winddiag = onesec_plot_df['winddiag']
        # Extract indices for "bad" wind data
        winddiag_index = N.where(N.any([winddiag > 0, N.isnan(winddiag)], axis=0))[0]
        # These are the times with bad wind data
        winddiag_plot = plottimes[winddiag_index]
        fields.append(winddiag_plot)
        fieldparamdicts.append(winddiag_params)

    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    fields = [winddiravgvec]
    fieldparamdicts = [winddir_params]
    ax2 = plotmeteogram(ax2, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, pc.meteo_ws_range],
                    'axeslabels': [pc.timelabel, r'wind speed (m s$^{-1}$)']}
    axparamdict2 = {'majorylocator': ticker.MultipleLocator(45.),
                    'axeslimits': [None, [0.0, 360.0]],
                    'axeslabels': [None, r'Wind direction ($^{\circ}$C)']}
    axparamdicts = [axparamdict1, axparamdict2]
    ax1, ax2 = set_meteogram_axes([ax1, ax2], axparamdicts)

    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_wind.png', dpi=300)
    plt.close(fig)

    # Plot temperature and dewpoint
    tavgintv = 10  # Currently not used

    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)

    fields = [onesec_plot_df['fasttemp'].values, onesec_plot_df['dewpoint'].values]
    fieldparamdicts = [temp_params, dewpoint_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    ax1.axhline(0.0, ls=':', color='k')

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, pc.meteo_T_Td_range],
                    'axeslabels': [pc.timelabel, r'Temperature ($^{\circ}$C)']}
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_T_Td.png', dpi=300)
    plt.close(fig)

    # Plot relative humidity
    avgintv = 10  # Currently not used

    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)

    fields = [onesec_plot_df['RH_derived'].values]
    fieldparamdicts = [RH_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [0., 100.]],
                    'axeslabels': [pc.timelabel, 'Relative Humidity (%)']}

    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_RH.png', dpi=300)
    plt.close(fig)

    # Plot station pressure

    pmin = onesec_plot_df['pressure'].values.min()
    pmax = onesec_plot_df['pressure'].values.max()
    pmean = onesec_plot_df['pressure'].values.mean()
    avgintv = 1  # Currently not used

    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)

    fields = [onesec_plot_df['pressure'].values]
    fieldparamdicts = [pressure_params]
    ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator,
                    'axeslimits': [xaxislimits, [pmin - 2.5, pmax + 2.5]],
                    'axeslabels': [pc.timelabel, 'Station pressure (hPa)']}
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)
    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_pres.png', dpi=300)
    plt.close(fig)

    # Plot some additional diagnostic time series
    if(pc.plot_diagnostics):
        # Battery voltages
        fig = plt.figure(figsize=(5, 3))
        ax1 = fig.add_subplot(111)

        fields = [onesec_plot_df['voltage'].values]
        fieldparamdicts = [battery_params]
        ax1 = plotmeteogram(ax1, [plottimes], fields, fieldparamdicts)

        axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                        'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [8.0, 16.0]],
                        'axeslabels': [pc.timelabel, r'Battery Voltage (V)']}
        axparamdicts = [axparamdict1]
        ax1, = set_meteogram_axes([ax1], axparamdicts)

        plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_voltage.png', dpi=300)


def plotDSDderivedmeteograms(dis_index, pc, ib, **PSDderiveddict):
    """Plots meteograms of the various derived DSD quantities from the PIPS"""
    PSDmidtimes = PSDderiveddict.get('PSDmidtimes')
    PSD_plot_df = PSDderiveddict.get('PSD_plot_df')

    xaxislimits = [PSDmidtimes[0], PSDmidtimes[-1]]
    dis_name = ib.dis_name_list[dis_index]

    # Rain rates (intensities)
    fig = plt.figure(figsize=(5, 3))
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
    fig = plt.figure(figsize=(5, 3))
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
    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)

    fields = [PSD_plot_df['pcount'].values, PSD_plot_df['pcount2'].values]
    fieldparamdicts = [pcount_params, pcount2_params]
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
    fig = plt.figure(figsize=(5, 3))
    ax1 = fig.add_subplot(111)

    fields = [PSD_plot_df['amplitude'].values]
    fieldparamdicts = [amplitude_params]
    ax1 = plotmeteogram(ax1, [PSDmidtimes], fields, fieldparamdicts)

    axparamdict1 = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                    'minorxlocator': pc.minorlocator, 'axeslimits': [xaxislimits, [0.0,  30000.0]],
                    'axeslabels': [pc.timelabel, r'Signal amplitude']}
    axparamdicts = [axparamdict1]
    ax1, = set_meteogram_axes([ax1], axparamdicts)

    plt.savefig(ib.image_dir + 'meteograms/' + dis_name + '_amplitude.png', dpi=300)
    plt.close(fig)


def plotDSDmeteograms(dis_name, image_dir, axparams, disvars, radvars):
    """Plots one or more meteograms of disdrometer number concentrations vs. diameter bins,
       along with one or more derived variables and optionally radar variables for comparison.
       One meteogram is plotted per dualpol variable (i.e. Z,ZDR,KDP,RHV)"""
    min_diameter = disvars.get('min_diameter', N.empty((0)))
    PSDstarttimes = disvars.get('PSDstarttimes', N.empty((0)))
    PSDmidtimes = disvars.get('PSDmidtimes', N.empty((0)))
    logND = disvars.get('logND', N.empty((0)))
    if(not logND.size or not PSDstarttimes.size or not PSDmidtimes.size):
        print "No DSD info to plot! Quitting!"
        return
    D_0_dis = disvars.get('D_0', N.empty((0)))
    D_0_rad = radvars.get('D_0_rad', N.empty((0)))
    radmidtimes = radvars.get('radmidtimes', N.empty((0)))
    dBZ_ray_dis = disvars.get('dBZ_ray', N.empty((0)))
    flaggedtimes = disvars.get('flaggedtimes', N.empty((0)))
    hailflag = disvars.get('hailflag', N.empty((0)))

    # Try to find the desired dualpol variables for plotting in the provided dictionary

    dualpol_dis_varnames = []
    dualpol_dis_vars = []
    for key, value in disvars.iteritems():
        if key in ['dBZ', 'ZDR', 'KDP', 'RHV']:
            dualpol_dis_varnames.append(key)
            dualpol_dis_vars.append(value)

    # If no dualpol variables are in the dictionary, we'll still go ahead and plot the
    # meteogram, with just the number concentrations and possibly median volume diameters, etc.
    if(not dualpol_dis_varnames):
        dualpol_dis_varnames.append(None)
        dualpol_dis_vars.append(N.empty((0)))

    # Start the plotting loop
    for dualpol_dis_varname, dualpol_dis_var in zip(dualpol_dis_varnames, dualpol_dis_vars):
        # See if the variable is also provided in the radvars dictionary
        dualpol_rad_var = radvars.get(dualpol_dis_varname, N.empty((0)))

        # Create the figure
        fig = plt.figure(figsize=(8, 3))
        ax1 = host_subplot(111)
        fig.autofmt_xdate()
        if(dualpol_dis_var.size):
            ax2 = ax1.twinx()
        divider = make_axes_locatable(ax1)

        xvals = [PSDstarttimes]
        plotvars = [logND]
        plotparamdict1 = {'type': 'pcolor', 'vlimits': [-1.0, 3.0],
                          'clabel': r'log[N ($m^{-3} mm^{-1}$)]'}
        plotparamdicts = [plotparamdict1]

        # Median volume diameter
        if(D_0_dis.size):
            xvals.append(PSDmidtimes)
            plotvars.append(D_0_dis)
            plotparamdict2 = {'type': 'line', 'linestyle': ':', 'color': 'b', 'linewidth': 0.5,
                              'label': r'$D_{0,dis} (mm)$'}
            plotparamdicts.append(plotparamdict2)

        # Vertical lines for flagged times (such as from wind contamination).
        if(flaggedtimes.size):
            xvals.append(PSDmidtimes)
            plotvars.append(flaggedtimes)
            print "flagged times", flaggedtimes
            plotparamdict = {'type': 'vertical line', 'linestyle': '-', 'color': 'r',
                             'linewidth': 0.5}
            plotparamdicts.append(plotparamdict)

        # Mark times with hail detected with a vertical purple line
        if(hailflag.size):
            xvals.append(PSDmidtimes)
            plotvars.append(hailflag)
            plotparamdict = {'type': 'vertical line', 'linestyle': '-', 'color': 'purple',
                             'linewidth': 0.5}
            plotparamdicts.append(plotparamdict)

        ax1 = plotmeteogram(ax1, xvals, plotvars, plotparamdicts,
                            yvals=[min_diameter] * len(plotvars))

        axparamdict1 = axparams
        axes = [ax1]
        axparamdicts = [axparamdict1]

        # Now plot the dualpol variable
        if(dualpol_dis_var.size):
            xvals = [PSDmidtimes]
            plotvars = [dualpol_dis_var]
            if(dualpol_dis_varname == 'dBZ'):
                dualpol_dis_varlabel = r'$Z_{dis} (dBZ)$'
                dualpol_rad_varlabel = r'$Z_{rad} (dBZ)$'
                axis_label = 'Reflectivity (dBZ)'
                axis_limits = [0.0, 80.0]
                axis_intv = 5.0

                # Also plot Rayleigh version
                if(dBZ_ray_dis.size):
                    xvals.append(PSDmidtimes)
                    plotvars.append(dBZ_ray_dis)
                    dBZ_ray_dis_varlabel = r'$Z_{dis,ray} (dBZ)$'
            elif(dualpol_dis_varname == 'ZDR'):
                dualpol_dis_varlabel = r'$Z_{DR,dis} (dB)$'
                dualpol_rad_varlabel = r'$Z_{DR,rad} (dB)$'
                axis_label = 'Differential Reflectivity (dBZ)'
                axis_limits = [0.0, 6.0]
                axis_intv = 0.5
            elif(dualpol_dis_varname == 'KDP'):
                dualpol_dis_varlabel = r'$K_{DP,dis} (deg km^{-1})$'
                dualpol_rad_varlabel = r'$K_{DP,rad} (deg km^{-1})$'
                axis_label = r'Specific Differential Phase (deg km$^{-1}$)'
                axis_limits = [0.0, 12.0]
                axis_intv = 1.0
            elif(dualpol_dis_varname == 'RHV'):
                dualpol_dis_varlabel = r'$\sigma_{HV,dis}$'
                dualpol_rad_varlabel = r'$\sigma_{HV,rad}$'
                axis_label = r'Cross-correlation Coefficient'
                axis_limits = [0.0, 1.05]
                axis_intv = 0.1

            plotparamdict1 = {'type': 'line', 'linestyle': '-', 'color': 'b', 'linewidth': 1.5,
                              'label': dualpol_dis_varlabel}
            plotparamdicts = [plotparamdict1]
            if(dualpol_dis_varname == 'dBZ' and dBZ_ray_dis.size):
                plotparamdict = {'type': 'line', 'linestyle': '--', 'color': 'b', 'linewidth': 1.5,
                                 'label': dBZ_ray_dis_varlabel}
                plotparamdicts.append(plotparamdict)

            if(dualpol_rad_var.size):
                xvals.append(radmidtimes)
                plotvars.append(dualpol_rad_var)
                plotparamdict2 = {'type': 'line', 'linestyle': '-', 'color': 'purple',
                                  'linewidth': 1.5, 'label': dualpol_rad_varlabel}
                plotparamdicts.append(plotparamdict2)

            # print xvals,plotvars,plotparamdicts
            ax2 = plotmeteogram(ax2, xvals, plotvars, plotparamdicts)
            axparamdict2 = {'majorylocator': ticker.MultipleLocator(base=axis_intv),
                            'axeslimits': [axparams['axeslimits'][0], axis_limits], 'axeslabels': [None, axis_label]}
            axes.append(ax2)
            axparamdicts.append(axparamdict2)
            ax2.legend(bbox_to_anchor=(1., 1.), loc='upper right',
                       ncol=1, fancybox=True, shadow=False, prop=fontP)

        axes = set_meteogram_axes(axes, axparamdicts)
        if(dualpol_dis_varname):
            plt.savefig(image_dir + dis_name + '_' + dualpol_dis_varname + '_logNc.png', dpi=300)
        else:
            plt.savefig(image_dir + dis_name + '_logNc.png', dpi=300)
        plt.close(fig)


def plotmeteogram(ax, xvals, zvals, plotparamdicts, yvals=None):
    """Plots a meteogram (time series) of one or more meteorological variables"""
    ax = ax or plt.figure().add_subplot(111)
    for idx, xval, zval, plotparamdict in zip(xrange(len(zvals)), cycle(xvals), zvals,
                                              plotparamdicts):
        type = plotparamdict.get('type', 'fill_between')
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

        if(type == 'fill_between'):
            ax.plot_date(xval, zval, ls=linestyle, lw=linewidth, marker=marker, color=color,
                         markeredgecolor=markeredgecolor, ms=ms, label=plotlabel)
            ax.fill_between(xval, zval, plotmin, facecolor=color, alpha=alpha)
        elif(type == 'pcolor'):
            divider = make_axes_locatable(ax)
            C = ax.pcolor(xval, yvals[idx], zval, vmin=vlimits[0], vmax=vlimits[1])
            cax = divider.append_axes("bottom", size="5%", pad=0.40)
            cb = ax.get_figure().colorbar(C, cax=cax, orientation='horizontal')
            if(clabel):
                cb.set_label(clabel)
        elif(type == 'vertical line'):  # For flagging times with bad data, etc.
                                        # zval is interpreted as a list of x-indices
            for x in zval:
                ax.axvline(x=x, ls=linestyle, lw=linewidth, color=color)
        else:
            ax.plot_date(xval, zval, ls=linestyle, lw=linewidth, marker=marker, color=color,
                         markeredgecolor=markeredgecolor, ms=ms, label=plotlabel)
    return ax


def set_meteogram_axes(axes, axparamdicts):
    """Sets up and formats the axes for a meteogram plot"""
    for ax, axparamdict in zip(axes, axparamdicts):
        ax = ax or plt.figure().add_subplot(111)
        majorxlocator = axparamdict.get('majorxlocator', None)
        majorxformatter = axparamdict.get('majorxformatter', None)
        majorylocator = axparamdict.get('majorylocator', None)
        majoryformatter = axparamdict.get('majoryformatter', None)
        minorxlocator = axparamdict.get('minorxlocator', None)
        axeslimits = axparamdict.get('axeslimits', [None, None])
        axeslabels = axparamdict.get('axeslabels', [None, None])

        if(majorxlocator):
            ax.xaxis.set_major_locator(majorxlocator)
        if(majorxformatter):
            ax.xaxis.set_major_formatter(majorxformatter)
        if(minorxlocator):
            ax.xaxis.set_minor_locator(minorxlocator)
        if(axeslimits[0]):
            ax.set_xlim(axeslimits[0][0], axeslimits[0][1])
        if(axeslabels[0]):
            ax.set_xlabel(axeslabels[0])
        if(axeslabels[1]):
            ax.set_ylabel(axeslabels[1])
        if(axeslimits[1]):
            ax.set_ylim(axeslimits[1][0], axeslimits[1][1])
        if(majorylocator):
            ax.yaxis.set_major_locator(majorylocator)

    if(axes[0]):
        axes[0].get_figure().autofmt_xdate()

    return axes


def plot_DSD(ib, axdict, PSDdict, PSDfitdict, PSDparamdict):
    """Plots an individual measured PSD on a semilog plot, along with optional exponential and gamma
       fits and associated parameters."""

    times = axdict.get('times', N.empty((0)))
    xbin_left = axdict.get('xbin_left', N.empty((0)))
    xbin_right = axdict.get('xbin_right', N.empty((0)))
    xbin_mid = axdict.get('xbin_mid', N.empty((0)))
    ND = PSDdict.get('ND', N.empty((0)))
    ND_onedrop = PSDdict.get('ND_onedrop', N.empty((0)))
    interval = axdict.get('interval', 10)
    dis_name = axdict.get('dis_name', None)
    t = axdict.get('time', None)

#     print ND*1000.
#     print xbin_left
#     print xbin_right

    fig1 = plt.figure(figsize=(8, 6))
    ax1 = fig1.add_subplot(111)
    plt.title('{0:d}-s DSD fits for time {1} EST'.format(interval, times[t].strftime(tm.timefmt2)))

    ax1.bar(xbin_left, ND * 1000.0,
            xbin_right - xbin_left, 10.**2., align='edge', log=True, color='tan', edgecolor='k')

    # ax1.plot(xbin_mid, ND * 1000.0, lw=2, label='obs')

    for fitname, ND_tuple in PSDfitdict.iteritems():
        if ((fitname == 'Dis Retr') | (fitname == 'Rad Retr')):
            ND_fit = ND_tuple[0]*1000.
        else:
            ND_fit = ND_tuple[0]
#         print ND_fit
        label = ND_tuple[1]
        if(ND_fit.size):
            ax1.plot(xbin_mid, ND_fit, lw=2, label=label)

    ax1.set_yscale('log')
    ax1.set_ylim(10.**2.0, 10.**8.5)
    ax1.set_xlim(0.0, 9.0)
    ax1.set_xlabel('D (mm)')
    ax1.set_ylabel(r'N(D) $(m^{-4})$')
    ypos = 0.95
    for paramname, paramtuple in PSDparamdict.iteritems():
        ax1.text(0.50, ypos, paramtuple[1] + ' = %2.2f'%paramtuple[0],
                 transform=ax1.transAxes)
        ypos = ypos - 0.05

    plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=8)
    plt.savefig(ib.image_dir + 'DSDs/' + dis_name + '/' + dis_name + '_DSD_t{:04d}.png'.format(t),
                dpi=200, bbox_inches='tight')
    plt.close(fig1)

def plot_vel_D(ib, axdict, PSDdict, rho):
    """Plots the terminal velocity vs. diameter matrix for a given DSD"""

    times = axdict.get('times', N.empty((0)))
    dis_name = axdict.get('dis_name', None)
    t = axdict.get('time', None)
    xlim = axdict.get('xlim', (0.0, 9.0))
    ylim = axdict.get('ylim', (0.0, 15.0))
    min_diameter = axdict.get('min_diameter', None)
    min_fall_bins = axdict.get('min_fall_bins', None)
    avg_diameter = axdict.get('avg_diameter', None)

    countsMatrix = PSDdict.get('countsMatrix', None)
    flaggedtime = PSDdict.get('flaggedtime', 0)


    fig1 = plt.figure(figsize=(8, 6))
    ax1 = fig1.add_subplot(111)
    plt.title('Fall speed vs. diameter for time {0}'.format(times[t].strftime(tm.timefmt2)))

    countsplot = N.ma.masked_where(countsMatrix[:] <= 0, countsMatrix[:])
    C = ax1.pcolor(min_diameter, min_fall_bins, countsplot, vmin=1, vmax=50, edgecolors='w')
    rainvd = dis.assignfallspeed(avg_diameter, rhocorrect=True, rho=rho)

    ax1.plot(avg_diameter, rainvd, c='r')
    # ax1.scatter(X[0:10,20:31],Y[0:10,20:31],c='r',marker='x')
    fig1.colorbar(C)

    if(flaggedtime > 1):
        ax1.text(0.5, 0.5, 'Flagged for strong wind contamination!',
                 horizontalalignment='center',
                 verticalalignment='center', color='y',
                 transform=ax1.transAxes)
    if(dis.plot_strongwindQC):
        ax1.scatter(X[dis.strongwindmask], Y[dis.strongwindmask], c='r', marker='x', alpha=1.0)
    if(dis.plot_splashingQC):
        ax1.scatter(X[dis.splashmask], Y[dis.splashmask], c='w', marker='o', alpha=0.75)
        # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(splashmask,mask=-splashmask),cmap=cm.Reds,alpha=0.1)
    if(dis.plot_marginQC):
        ax1.scatter(X[dis.marginmask], Y[dis.marginmask], c='g', marker='x', alpha=0.1)
        # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(marginmask,mask=-marginmask),cmap=cm.Reds,alpha=0.1)
    if(dis.plot_rainfallspeedQC):
        ax1.scatter(X[dis.fallspeedmask], Y[dis.fallspeedmask], c='k', marker='x', alpha=0.5)
        # ax1.pcolor(min_diameter,min_fall_bins,ma.masked_array(fallspeedmask,mask=-fallspeedmask),cmap=cm.gray,alpha=0.1)
    if(dis.plot_rainonlyQC):
        ax1.scatter(X[dis.rainonlymask], Y[dis.rainonlymask], c='g', marker='x', alpha=0.5)

    ax1.set_xlim(xlim[0], xlim[1])
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax1.set_xlabel('diameter (mm)')
    ax1.set_ylim(ylim[0], ylim[1])
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax1.set_ylabel('fall speed (m/s)')

    plt.savefig(ib.image_dir + 'vel_D/' + dis_name + '/' + dis_name + '_vel_D_t{:04d}.png'.format(t),
                dpi=200, bbox_inches='tight')
    plt.close(fig1)

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
#             D_med_disd_avg = pd.rolling_mean(D_med_disd,12,center=True,win_type='triang',min_periods=1)
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
#             refl_disd_avg = pd.rolling_mean(refl_disd,12,center=True,win_type='triang',min_periods=1)
#             #ax2.plot(plotx_dis_avg,refl_disd,ls='-',c='r',marker='o',label=r'$Z_{dis} (dBZ)$')
#             ax2.plot(plotx_dis_avg,refl_disd_avg,ls='-',c='green',lw=0.5,marker=None,label=r'$Z_{dis} (dBZ)$')
#             if(comp_radar):
#                 dBZ_D_plt = fields_D_tarr[:,index,0] # Assumes reflectivity is the first.  Need to come back and work on this.
#                 #ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',marker='o',label=r'$Z_{rad} (dBZ)$')
#                 ax2.plot(plotx_rad,dBZ_D_plt,ls='-',c='k',lw=0.5,marker=None,label=r'$Z_{rad} (dBZ)$')
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
