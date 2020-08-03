# pyXTRRA.py
# A module containing some basic utilities for reading and plotting XTRRA data.
# Credit to S. Harrell, D. Dawson, and M. Baldwin for some of the base code.
# RLT 20181011


import datetime
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.io.shapereader import Reader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pyart
from pyart.config import get_metadata
from pyart.core import Radar
import glob
import os
from matplotlib.colors import LinearSegmentedColormap
import pytz

pyart.config.load_config('../configs/xtrra_config.py')

figure_directory = "../figures"
#shapefile_directory = '/home/rtanama/public/shapefiles'
#countyshapefile = os.path.join(shapefile_directory, 'countyl010g.shp')
#counties = cfeat.ShapelyFeature(Reader(countyshapefile).geometries(), ccrs.PlateCarree(),
#                                facecolor='none')

def CreateRadarFromRXM25netCDF(netcdf_file, instrument_name, cfradial_outfile=None, heading=None):
    data  = netCDF4.Dataset(netcdf_file, 'r')
    #fileTime = datetime.datetime.strptime(netcdf_file[-22:-7],'%Y%m%d-%H%M%S')


    ngates = data.dimensions['Gate'].size
    rays_per_sweep = data.dimensions['Radial'].size
    radar = pyart.testing.make_empty_ppi_radar(ngates, rays_per_sweep, 1)

    # Time needs to be converted from nss1970 to nss1980 and added to radar object
    nineteen89 = datetime.datetime(1989, 1, 1, 0, 0, 1, tzinfo = pytz.utc)
    baseTime = np.array([datetime.datetime.fromtimestamp(t,tz=pytz.UTC) for t in data.variables['Time'][:]])
    radar.time['data'] = np.array([t.total_seconds() for t in baseTime - nineteen89])

    if heading is not None:
        radar.heading = heading
        radar.azimuth['data'] = np.mod(data['Azimuth'][:] - radar.heading, 360.)
    else:
        radar.azimuth['data'] = data['Azimuth'][:]

    radar.longitude['data'] = np.array([data.Longitude], dtype='float64')
    radar.latitude['data'] = np.array([data.Latitude], dtype='float64')
    radar.elevation['data'] = data['Elevation'][:]
    radar.altitude['data'] = np.array([data.Height], dtype='float64')

    fixed_agl_data = np.empty((1, ), dtype='float32')
    fixed_agl_data[:] = np.mean(
        radar.elevation['data'][:rays_per_sweep])

    radar.fixed_angle['data'] = fixed_agl_data

    radar.range['data'] = np.linspace(
        data['StartRange'][0]/1000,
       (ngates - 1)*data['GateWidth'][0]/1000 + data['StartRange'][0]/1000,
       ngates)

    ref = data['Reflectivity'][:]
    norm_pow = data['NormalizedCoherentPower'][:]
    spec_w = data['SpectralWidth'][:]
    vel = data['Velocity'][:]
    corr_ref = data['CorrectedReflectivity'][:]
    diff_ref = data['DifferentialReflectivity'][:]
    diff_phase = data['DifferentialPhase'][:]
    spec_phase = data['SpecificPhase'][:]
    corr_diff_ref = data['CorrectedDifferentialReflectivity'][:]
    sig_noise = data['SignalToNoiseRatio'][:]
    rain_rate = data['RainfallRate'][:]
    cross_ra = data['CrossPolCorrelation'][:]

    fields = {
        'reflectivity': get_metadata('reflectivity'),
        'normalized_coherent_power': get_metadata('normalized_coherent_power'),
        'spectral_width': get_metadata('spectral_width'),
        'velocity': get_metadata('velocity'),
        'corrected_reflectivity': get_metadata('correct_reflectivity'),
        'differential_reflectivity': get_metadata('differential_reflectivity'),
        'differential_phase': get_metadata('differential_phase'),
        'specific_differential_phase': get_metadata('specific_differential_phase'),
        'corrected_differential_reflectivity': get_metadata('corrected_differential_reflectivity'),
        'signal_to_noise_ratio': get_metadata('signal_to_noise_ratio'),
        'rain_rate': get_metadata('rain_rate'),
        'cross_correlation_ratio': get_metadata('cross_correlation_ratio')}

    radar.fields = fields
#    radar.fields['reflectivity']['data'] = ref
#    radar.fields['normalized_coherent_power']['data'] = norm_pow
#    radar.fields['spectral_width']['data'] = spec_w
#    radar.fields['velocity']['data'] = vel
#    radar.fields['corrected_reflectivity']['data'] = corr_ref
#    radar.fields['differential_reflectivity']['data'] = diff_ref
#    radar.fields['differential_phase']['data'] = diff_phase
#    radar.fields['specific_differential_phase']['data'] = spec_phase
#    radar.fields['corrected_differential_reflectivity']['data'] = corr_diff_ref
    radar.fields['reflectivity']['data'] = np.ma.masked_where(sig_noise < 3.,ref)
    radar.fields['normalized_coherent_power']['data'] = norm_pow
    radar.fields['spectral_width']['data'] = np.ma.masked_where(sig_noise < 3.,spec_w)
    radar.fields['velocity']['data'] = np.ma.masked_where(sig_noise < 3.,vel)
    radar.fields['corrected_reflectivity']['data'] = np.ma.masked_where(sig_noise < 3.,corr_ref)
    radar.fields['differential_reflectivity']['data'] = np.ma.masked_where(sig_noise < 3.,diff_ref)
    radar.fields['differential_phase']['data'] = np.ma.masked_where(sig_noise < 3.,diff_phase)
    radar.fields['specific_differential_phase']['data'] = np.ma.masked_where(sig_noise < 3.,spec_phase)
    radar.fields['corrected_differential_reflectivity']['data'] = np.ma.masked_where(sig_noise < 3.,corr_diff_ref)
    radar.fields['signal_to_noise_ratio']['data'] = sig_noise
    radar.fields['rain_rate']['data'] = rain_rate
    radar.fields['cross_correlation_ratio']['data'] = cross_ra

    radar.metadata['instrument_name'] = instrument_name
    if cfradial_outfile is not None:
        pyart.io.write_cfradial(cfradial_outfile, radar, arm_time_variables=True)

    return radar

def CreatePlot(radar, variable, fname, elevation, vmin, vmax, timestamp=None, title=False,
               ovrmap=False, verbose=False, separate_colorbar=True):
    ###### VERY IMPORTANT - ELEVATION IS NOT USED IN ANY WAY YET ####
    if ovrmap:
        display = pyart.graph.RadarMapDisplayCartopy(radar)
    else:
        display = pyart.graph.RadarDisplay(radar)
    #plt.axis('off')
    f = plt.figure(figsize = [8, 8])
#     projection = ccrs.AzimuthalEquidistant(central_longitude=radar.longitude['data'][0],
#                                            central_latitude=radar.latitude['data'][0])
#     ax = plt.axes(projection=projection)
#     # Gotcha. Matplotlib transparency doesn't work with Cartopy. Need to use the following
#     # from https://stackoverflow.com/questions/30076560/make-a-transparent-plot-with-cartopy
#     ax.outline_patch.set_visible(False)
#     ax.background_patch.set_visible(False)
#     ax.axis('off')
    if ovrmap:
        display.plot_ppi_map(variable['radar_name'], sweep=0, vmin=vmin, vmax=vmax,
                             title_flag=False, colorbar_flag=False, mask_outside=True, alpha=0.7)
        display.ax.axis('off')
        projection = display.grid_projection
        display.ax.add_feature(counties, edgecolor='brown', alpha=0.7, lw=0.75)
        display.ax.add_feature(cfeat.NaturalEarthFeature('cultural', 'roads', '10m'),
                               facecolor='none', edgecolor='b')
    else:
        # Have to do this here because we still need map projection information in this case
        # even though we are using the regular plot_ppi that doesn't care about it. The reason
        # we want to use the regular plot_ppi in the first place is that plot_ppi_map does stuff
        # like plot lat/lon lines when we don't want them, etc. This way we recover the original
        # behavior of a transparent background that we can overlay on the Google Maps but still
        # have the projection information we can use to get the correct horizontal extents for the
        # overlay.
        projection = ccrs.AzimuthalEquidistant(central_longitude=radar.longitude['data'][0],
                                               central_latitude=radar.latitude['data'][0])
        display.ax = plt.axes(projection=projection)
        # Gotcha. Matplotlib transparency doesn't work with Cartopy. Need to use the following
        # from https://stackoverflow.com/questions/30076560/make-a-transparent-plot-with-cartopy
        display.ax.outline_patch.set_visible(False)
        display.ax.background_patch.set_visible(False)
        display.ax.axis('off')
        display.plot_ppi(variable['radar_name'], sweep=0, vmin=vmin, vmax=vmax, title_flag=False,
                         axislabels_flag=False, colorbar_flag=False, mask_outside=True, alpha=0.7)

    display.plot_range_rings(np.arange(10., 60., 10.), ls = '--', lw=0.5) # RLT: Stopped rings at 50 km (where data end)
    # Something isn't working properly with the axes transforms below. I get a bizarre error:
    # TypeError: descriptor '_as_mpl_transform' requires a 'cartopy._crs.CRS' object but received a 'GeoAxesSubplot'
    # Possibly a bug? EDIT 08/24/18: I think I know the problem. I need to put () after PlateCarree.
    # Will revisit this
    # So for now I'm just going to use what I know is the correct transformed Azimuthal Equidistant projection
#     display.ax.plot([-86.912196], [40.430134], color='k', marker='*', transform=ccrs.PlateCarree)
#     display.ax.text(-86.912196+0.02, 40.430134, 'XTRRA', fontsize=7, transform=ccrs.PlateCarree)
    # Using 0., 0. as the position since the projection is centered on the radar location by default
    display.ax.plot([0.], [0.], color='k', marker='*')
    display.ax.text(1000., 0., 'XTRRA', fontsize=7)
    if title:
        timestamptitle = '{:d}-{:d}-{:d} {}:{}:{}'.format(int(timestamp[4:6]), int(timestamp[6:8]),
                                                        int(timestamp[:4]), timestamp[9:11],
                                                        timestamp[11:13], timestamp[13:15])
        title_1 = str('XTRRA')+' : Lat '+str(radar.latitude['data'][0])+' deg, Lon '+str(radar.longitude['data'][0])+' deg'
        title_2 = variable['short_name'] +' ('+variable['units']+')'
        plot_title = title_1+'\n'+title_2+'\n'+timestamptitle+' UTC \n'+'Elevation angle: {:.1f} (deg)'.format(elevation)
        plt.title(plot_title, loc='left')
    # Can set the axes extents to whatever we want here. Here we use the native Azimuthal Equidistant
    # projection of the geoaxes associated with the pyart display (i.e. display.ax)). The projection
    # instance that PyART is using is stored in display.grid_projection. So, let's set it to go out to
    # 70 km in each direction from the radar.
    if ovrmap:
        display.ax.set_extent((-70000., 70000., -70000., 70000.), crs=projection)
    else:
        # There's a bug (feature?) in Cartopy that appears to have units of km being used for
        # the regular plot_ppi/RadarDisplay vs. units of m being used for the plot_ppi_map/RadarDisplayCartopy!
        # Probably should open up an issue on pyArt github. Better yet, see if there's a way we
        # can get better consistency overall between the various plotting methods so that we don't
        # have to finagle around with projections like above.
        display.ax.set_extent((-70., 70., -70., 70.), crs=projection)
    # Now, if we call get_extent() and feed it the PlateCarree transform, we should get the lat/lon
    # points of the corners back. This information can then (presumably) be used to locate the
    # domain in another map if desired.

    # latlonbounds = display.ax.get_extent(crs=ccrs.PlateCarree())
    if verbose:
        print("Extent of domain in lat-lon coordinates: ", latlonbounds)
    # Try the reverse transformation just to make sure we get the expected result back
    # Answer: yes!
    # testbounds = display.ax.get_extent(crs=projection)
    # print testbounds
    extent = display.ax.get_window_extent().transformed(f.dpi_scale_trans.inverted())
    if timestamp is not None:
        fname = os.path.splitext(fname)[0]+'_'+timestamp

#     print display.ax.collections
    if separate_colorbar:
        cbar_name = 'scale-%s.png' % variable['website_name']
        colorbar_label = '%s' % variable['short_name']
        fcb = plt.figure(figsize = [4,1])
        axcb = fcb.add_subplot(111)
        display.plot_colorbar(mappable=display.ax.collections[0], label=colorbar_label,
                      orient='horizontal', ax=axcb, fig=fcb)
        fcb.patch.set_facecolor('white')
        fcb.patch.set_alpha(0.7)
        axcb.set_visible(False)
        fcb.savefig('%s/%s' % (figure_directory, cbar_name), bbox_inches='tight', dpi=100)
        plt.close(fcb)

    f.savefig('%s/%s' % (figure_directory, fname), transparent=True, dpi=600, bbox_inches=extent,
              pad_inches=0)
    plt.close(f)

# DTD: Added "units" key to variable dictionaries below
# RLT: adjusted some variable max/min values to align with typical observations
variables = [
        {'short_name': 'Z',
         'radar_name': 'reflectivity',
         'website_name': 'z',
         'max': 60,
         'min': 0,
         'units': 'dBZ' },
        {'short_name': 'V',
         'radar_name': 'velocity',
         'website_name': 'v',
         'max': 35,
         'min': -35,
         'units': r'm s$^{-1}$' },
        {'short_name': 'ZDR',
         'radar_name': 'differential_reflectivity',
         'website_name': 'zdr',
         'max': 8,
         'min': -2,
         'units': 'dB' },
        {'short_name': 'PDP',
         'radar_name': 'differential_phase',
         'website_name': 'pdp',
         'max': 360,
         'min': 0,
         'units': r'$^{\circ}$' },
        {'short_name': 'KDP',
         'radar_name': 'specific_differential_phase',
         'website_name': 'kdp',
         'max': 10,
         'min': -2,
         'units': r'$^{\circ}$ km$^{-1}$' },
        {'short_name': 'RHV',
         'radar_name': 'cross_correlation_ratio',
         'website_name': 'rhv',
         'max': 1,
         'min': 0,
         'units': '' },
        {'short_name': 'W',
         'radar_name': 'spectral_width',
         'website_name': 'w',
         'max': 25,
         'min': 0,
         'units': '' },
        ]