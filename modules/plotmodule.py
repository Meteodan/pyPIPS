#plotmodule.py: A module containing some functions related to plotting model output
import numpy as N
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import ImageGrid,make_axes_locatable,host_subplot
from matplotlib.ticker import *
from matplotlib.collections import Collection,LineCollection
from matplotlib.artist import allow_rasterization
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.transforms as mtransforms
from matplotlib.font_manager import FontProperties

fontP = FontProperties()
fontP.set_size('small')


RdBu1_colors = ("#84001D","#890025","#8D002D","#920034","#97103C","#9C2143",
          "#A12E4B","#A63A53","#AA455B","#AF5063","#B45A6C","#B96575",
          "#BE707E","#C37B88","#C88792","#CD949D","#D2A1A9","#D7B0B6",
          "#DDC1C5","#E3D6D8","#DAD8E2","#CAC5DA","#BDB6D4","#B2A9CF",
          "#A89DCA","#9F91C6","#9687C1","#8E7DBE","#8773BA","#7F69B7",
          "#7960B4","#7257B1","#6D4DAF","#6743AE","#6239AE","#5E2CAE",
          "#5A1AB1","#5800B7","#5800C2","#5C00D9")

# The following Class and method are from http://stackoverflow.com/questions/12583970/matplotlib-contour-plots-as-postscript
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

def mtokm(val,pos):
    """Convert m to km for formatting axes tick labels"""
    val=val/1000.0
    return '%i' % val

def mtokmr1(val,pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val=val/1000.0
    #return '{:2.{prec}f}'.format(val,prec=prec)
    return '%2.1f' % val

def mtokmr2(val,pos):
    """Convert m to km for formatting axes tick labels.
       Floating point version."""
    val=val/1000.0
    #return '{:2.{prec}f}'.format(val,prec=prec)
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

def colorline(x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0,ax=None):
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
    
    if(ax == None):
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
        
def plotsingle(fig,axes,ptype,xs,ys,x,y,xlim,ylim,field,clevels,cmap,fieldnorm,cbarlevels,clabel,cformat,ovrmap,gis_info,numovr,ovrx,ovry,ovrfields,ovrfieldlvls,ovrfieldcolors,axesticks,rasterized=True):
    """Plot a single-paneled figure from model output (possibly staggered grid)"""
    if(fig==None):
        fig = plt.figure()
    if(axes==None):
        axes = fig.add_subplot(111)
    if(rasterized):
        axes.set_rasterization_zorder(2)
    norm = fieldnorm
    if(ptype == 1): # Contour plot
#         if(True): # Hack!
#             plot = axes.contourf(ovrx[1],ovry[1],ovrfields[1],levels=ovrfieldlvls[1],cmap=cm.Greys_r,zorder=0,extend='both')
        plot = axes.contourf(xs,ys,field,levels=clevels,cmap=cmap,norm=norm,zorder=1)
    elif(ptype == 2): # pcolor plot
        # Mask values below lower bounds
        field = N.ma.masked_where(field < clevels[0],field)
        #norm = matplotlib.colors.BoundaryNorm(clevels,cmap.N)
        cmap.set_bad('white',alpha=None) # This is a hack to force masked areas to be white.  By default masked areas are
                                         # transparent, but the Agg backends apparently have problems when converting pcolor images
                                         # to eps and changes the transparent color to black.
        plot = axes.pcolormesh(x,y,field,vmin=clevels[0],vmax=clevels[-1],cmap=cmap,norm=norm,edgecolors='None',antialiased=False,rasterized=rasterized)
    # Find a nice number of ticks for the colorbar
    
    if(cbarlevels == None):
        cintv = clevels[1]-clevels[0]
        cintvs = N.arange(clevels[0],clevels[-1],cintv)
        while True:
            if(cintvs.size > 20):
                cintv = (cintvs[1]-cintvs[0])*2.
                cintvs = N.arange(cintvs[0],cintvs[-1],cintv)
            else:
                break
        cbarlevels = MultipleLocator(base=cintv)

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right",size="5%",pad=0.05)
    plt.colorbar(plot,orientation='vertical',ticks=cbarlevels,cax=cax,format=cformat)
    if(clabel != None):
        cax.set_ylabel(clabel)
    if(numovr > 0):
        for i in xrange(numovr):
#             if(i != 1): # Hack!
            plotovr = axes.contour(ovrx[i],ovry[i],ovrfields[i],levels=ovrfieldlvls[i],colors=ovrfieldcolors[i],lw=2)
    if(gis_info != None):
        axes.plot([gis_info[1]],[gis_info[2]],'ko')
    axes.set_xlim(xlim[0],xlim[1])
    axes.set_ylim(ylim[0],ylim[1])
    if(axesticks[0] < 1000.):
        if(axesticks[0] < 500.):
            formatter = FuncFormatter(mtokmr2)
        else:
            formatter = FuncFormatter(mtokmr1)
    else:
        formatter = FuncFormatter(mtokm)
    axes.xaxis.set_major_formatter(formatter)
    if(axesticks[1] < 1000.):
        if(axesticks[0] < 500.):
            formatter = FuncFormatter(mtokmr2)
        else:
            formatter = FuncFormatter(mtokmr1)
    else:
        formatter = FuncFormatter(mtokm)
    axes.yaxis.set_major_formatter(formatter)
    print axesticks[0],axesticks[1]
    axes.xaxis.set_major_locator(MultipleLocator(base=axesticks[0]))
    axes.yaxis.set_major_locator(MultipleLocator(base=axesticks[1]))
    axes.set_xlabel('km')
    axes.set_ylabel('km')
    if(axesticks[1] == axesticks[0]):
        axes.set_aspect('equal')
    else:
        axes.set_aspect('auto')
    
#     if(ovrmap): # Overlay map
#         readshapefile(track_shapefile_location,'track',drawbounds=True,linewidth=0.5,color='black',ax=axes)
#         readshapefile(county_shapefile_location,'counties',drawbounds=True, linewidth=0.5, color='gray',ax=axes)  #Draws US county boundaries.
    
    return fig,axes

def plotDSDmeteograms(dis_name,image_dir,axparams,disvars,radvars):
    """Plots one or more meteograms of disdrometer number concentrations vs. diameter bins,
       along with one or more derived variables and optionally radar variables for comparison.
       One meteogram is plotted per dualpol variable (i.e. Z,ZDR,KDP,RHV)"""
    min_diameter = disvars.pop('min_diameter',N.empty((0)))
    DSDstarttimes = disvars.pop('DSDstarttimes',N.empty((0)))
    DSDmidtimes = disvars.pop('DSDmidtimes',N.empty((0)))
    logNc_bin = disvars.pop('logNc_bin',N.empty((0)))
    if(not logNc_bin.size or not DSDstarttimes.size or not DSDmidtimes.size):
        print "No DSD info to plot! Quitting!"
        return
    D_0_dis = disvars.pop('D_0_dis',N.empty((0)))
    D_0_rad = radvars.pop('D_0_rad',N.empty((0)))
    radmidtimes = radvars.pop('radmidtimes',N.empty((0)))
    dBZ_ray_dis = disvars.pop('dBZ_ray',N.empty((0)))
    flaggedtimes = disvars.pop('flaggedtimes',N.empty((0)))
    hailflag = disvars.pop('hailflag',N.empty((0)))
    
    # Try to find the desired dualpol variables for plotting in the provided dictionary
    
    dualpol_dis_varnames = []
    dualpol_dis_vars = []
    for key,value in disvars.iteritems():
        if key in ['dBZ','ZDR','KDP','RHV']:
            dualpol_dis_varnames.append(key)
            dualpol_dis_vars.append(value)
    
    # If no dualpol variables are in the dictionary, we'll still go ahead and plot the 
    # meteogram, with just the number concentrations and possibly median volume diameters, etc.
    if(not dualpol_dis_varnames):
        dualpol_dis_varnames.append(None)
        dualpol_dis_vars.append(N.empty((0)))
        
    # Start the plotting loop
    for dualpol_dis_varname,dualpol_dis_var in zip(dualpol_dis_varnames,dualpol_dis_vars):
        # See if the variable is also provided in the radvars dictionary
        dualpol_rad_var = radvars.pop(dualpol_dis_varname,N.empty((0)))
    
        # Create the figure
        fig = plt.figure(figsize=(8,3))
        ax1 = host_subplot(111)
        fig.autofmt_xdate()
        if(dualpol_dis_var.size):
            ax2 = ax1.twinx()
        divider = make_axes_locatable(ax1)

        xvals = [DSDstarttimes]
        plotvars = [logNc_bin]
        plotparamdict1 = {'type':'pcolor','vlimits':[-1.0,3.0],'clabel':r'log[N ($m^{-3} mm^{-1}$)]'}
        plotparamdicts = [plotparamdict1]
        
        # Median volume diameter
        if(D_0_dis.size):
            xvals.append(DSDmidtimes)
            plotvars.append(D_0_dis)
            plotparamdict2 = {'type':'line','linestyle':':','color':'b','linewidth':0.5,'label':r'$D_{0,dis} (mm)$'}
            plotparamdicts.append(plotparamdict2)
        
        # Vertical lines for flagged times (such as from wind contamination).
        if(flaggedtimes.size):
            xvals.append(DSDmidtimes)
            plotvars.append(flaggedtimes)
            print "flagged times",flaggedtimes
            plotparamdict = {'type':'vertical line','linestyle':'-','color':'r','linewidth':0.5}
            plotparamdicts.append(plotparamdict)
        
        # Mark times with hail detected with a vertical purple line
        if(hailflag.size):
            xvals.append(DSDmidtimes)
            plotvars.append(hailflag)
            plotparamdict = {'type':'vertical line','linestyle':'-','color':'purple','linewidth':0.5}
            plotparamdicts.append(plotparamdict)

        ax1 = plotmeteogram(ax1,xvals,plotvars,plotparamdicts,yvals=[min_diameter]*len(plotvars))

        axparamdict1 = axparams
        axes = [ax1]
        axparamdicts = [axparamdict1]

        # Now plot the dualpol variable
        if(dualpol_dis_var.size):
            xvals = [DSDmidtimes]
            plotvars = [dualpol_dis_var]
            if(dualpol_dis_varname == 'dBZ'):
                dualpol_dis_varlabel = r'$Z_{dis} (dBZ)$'
                dualpol_rad_varlabel = r'$Z_{rad} (dBZ)$'
                axis_label = 'Reflectivity (dBZ)'
                axis_limits = [0.0,80.0]
                axis_intv = 5.0
                
                # Also plot Rayleigh version
                if(dBZ_ray_dis.size):
                    xvals.append(DSDmidtimes)
                    plotvars.append(dBZ_ray_dis)
                    dBZ_ray_dis_varlabel = r'$Z_{dis,ray} (dBZ)$'
            elif(dualpol_dis_varname == 'ZDR'):
                dualpol_dis_varlabel = r'$Z_{DR,dis} (dB)$'
                dualpol_rad_varlabel = r'$Z_{DR,rad} (dB)$'
                axis_label = 'Differential Reflectivity (dBZ)'
                axis_limits = [0.0,6.0]
                axis_intv = 0.5
            elif(dualpol_dis_varname == 'KDP'):
                dualpol_dis_varlabel = r'$K_{DP,dis} (deg km^{-1})$'
                dualpol_rad_varlabel = r'$K_{DP,rad} (deg km^{-1})$'
                axis_label = r'Specific Differential Phase (deg km$^{-1}$)'
                axis_limits = [0.0,12.0]
                axis_intv = 1.0
            elif(dualpol_dis_varname == 'RHV'):
                dualpol_dis_varlabel = r'$\sigma_{HV,dis}$'
                dualpol_rad_varlabel = r'$\sigma_{HV,rad}$'
                axis_label = r'Cross-correlation Coefficient'
                axis_limits = [0.0,1.05]
                axis_intv = 0.1
            
            plotparamdict1 = {'type':'line','linestyle':'-','color':'g','linewidth':0.5,'label':dualpol_dis_varlabel}
            plotparamdicts = [plotparamdict1]
            if(dualpol_dis_varname == 'dBZ' and dBZ_ray_dis.size):
                plotparamdict = {'type':'line','linestyle':'--','color':'g','linewidth':0.5,'label':dBZ_ray_dis_varlabel}
                plotparamdicts.append(plotparamdict)
            
            if(dualpol_rad_var.size):
                xvals.append(radmidtimes)
                plotvars.append(dualpol_rad_var)
                plotparamdict2 = {'type':'line','linestyle':'-','color':'k','linewidth':0.5,'label':dualpol_rad_varlabel}
                plotparamdicts.append(plotparamdict2)
            
            #print xvals,plotvars,plotparamdicts
            ax2 = plotmeteogram(ax2,xvals,plotvars,plotparamdicts)
            axparamdict2 = {'majorylocator':MultipleLocator(base=axis_intv),'axeslimits':[None,axis_limits],
                            'axeslabels':[None,axis_label]}
            axes.append(ax2)
            axparamdicts.append(axparamdict2)
            ax2.legend(bbox_to_anchor=(1.,1.), loc='upper right',
                                    ncol=1, fancybox=True, shadow=False, prop = fontP)
                
        axes = set_meteogram_axes(axes,axparamdicts)
        if(dualpol_dis_varname):
            plt.savefig(image_dir+dis_name+'_'+dualpol_dis_varname+'_logNc.png',dpi=300)
        else:
            plt.savefig(image_dir+dis_name+'_logNc.png',dpi=300)
        plt.close(fig)
            
def plotmeteogram(ax,xvals,zvals,plotparamdicts,yvals=None):
    """Plots a meteogram (time series) of one or more meteorological variables"""
    ax = ax or plt.figure().add_subplot(111)
    for idx,xval,zval,plotparamdict in zip(xrange(len(xvals)),xvals,zvals,plotparamdicts):
        type = plotparamdict.pop('type','fill_between')
        linestyle = plotparamdict.pop('linestyle','-')
        linewidth = plotparamdict.pop('linewidth',1.0)
        marker = plotparamdict.pop('marker',None)
        color = plotparamdict.pop('color',None)
        alpha = plotparamdict.pop('alpha',0.5)
        markeredgecolor = plotparamdict.pop('markeredgecolor','none')
        ms = plotparamdict.pop('ms',0)
        plotmin = plotparamdict.pop('plotmin',0)
        plotlabel = plotparamdict.pop('label',None)
        vlimits = plotparamdict.pop('vlimits',None)
        clabel = plotparamdict.pop('clabel',None)
        
        if(type == 'fill_between'):
            ax.plot_date(xval,zval,ls=linestyle,lw=linewidth,marker=marker,color=color,
                         markeredgecolor=markeredgecolor,ms=ms,label=plotlabel)
            ax.fill_between(xval,zval,plotmin,facecolor=color,alpha=alpha)
        elif(type == 'pcolor'):
            divider = make_axes_locatable(ax)
            C = ax.pcolor(xval,yvals[idx],zval,vmin=vlimits[0],vmax=vlimits[1])
            cax = divider.append_axes("bottom", size="5%", pad=0.40)
            cb = ax.get_figure().colorbar(C, cax=cax,orientation='horizontal')
            if(clabel):
                cb.set_label(clabel)
        elif(type == 'vertical line'):  # For flagging times with bad data, etc. zval is interpreted as a list of x-indices
            for x in zval:
                ax.axvline(x=x,ls=linestyle,lw=linewidth,color=color)
        else:
            ax.plot_date(xval,zval,ls=linestyle,lw=linewidth,marker=marker,color=color,
                         markeredgecolor=markeredgecolor,ms=ms,label=plotlabel)
    return ax
    
def set_meteogram_axes(axes,axparamdicts):
    """Sets up and formats the axes for a meteogram plot"""
    for ax,axparamdict in zip(axes,axparamdicts):
        ax = ax or plt.figure().add_subplot(111)
        majorxlocator = axparamdict.get('majorxlocator',None)
        majorxformatter = axparamdict.get('majorxformatter',None)
        majorylocator = axparamdict.get('majorylocator',None)
        majoryformatter = axparamdict.get('majoryformatter',None)
        minorxlocator = axparamdict.get('minorxlocator',None)
        axeslimits = axparamdict.get('axeslimits',[None,None])
        axeslabels = axparamdict.get('axeslabels',[None,None])
        
        if(majorxlocator):
            ax.xaxis.set_major_locator(majorxlocator)
        if(majorxformatter):
            ax.xaxis.set_major_formatter(majorxformatter)
        if(minorxlocator):
            ax.xaxis.set_minor_locator(minorxlocator)
        if(axeslimits[0]):
            ax.set_xlim(axeslimits[0][0],axeslimits[0][1])
        if(axeslabels[0]):
            ax.set_xlabel(axeslabels[0])
        if(axeslabels[1]):
            ax.set_ylabel(axeslabels[1])
        if(axeslimits[1]):
            ax.set_ylim(axeslimits[1][0],axeslimits[1][1])
        if(majorylocator):
            ax.yaxis.set_major_locator(majorylocator)
    
    if(axes[0]):
        axes[0].get_figure().autofmt_xdate()
    
    return axes

# Below is some extra code for DSD plotting with no home right now
    
    # 
#         #print "logNc_bin",logNc_bin
#         # Times are valid at end of DSD intervals
#         C = ax1.pcolor(plotx_dis_start[pstartindex:pstopindex+1],min_diameter,
#                        logNc_bin[:,pstartindex:pstopindex+1],vmin=-1.0,vmax=3.0)
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
      
