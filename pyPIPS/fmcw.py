"""Contains functions for reading and analyzing data from the UMass FMCW radar."""
import os
from glob import glob
import itertools
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
from skimage.restoration import unwrap_phase
from scipy import ndimage
from scipy.constants import c
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib import cm
from .utils import mtokm, interp_along_1D


# Some small functions needed throughout
# Source:
# https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
def nearest(items, pivot):
    """ Finds closest object to 'pivot' in list 'items' """
    return min(items, key=lambda x: abs(x - pivot))


# This function formatter will replace integers with target names
formatter = plt.FuncFormatter(lambda val, loc: categories[val+1])

# My favorite blue-to-red colormap
cdict1 = {
    'red':   ((0.0, 0.0, 0.0),
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

# Custom black-blue-red colormap for hydrometeor classification
hydroclass3 = LinearSegmentedColormap.from_list('hydroclass3',
                                                [(1, 0, 0), (0, 0, 1), (1, 1, 1)], N=3)

# Colormap for reflectivity that goes to blue instead of black
plasma_white = LinearSegmentedColormap.from_list('plasma_white',
                                                 ['white', 'purple', 'red', 'yellow'],
                                                 N=256)


def read_fmcw(filename):
    """
    Reads basic (moment) variables from a UMass FMCW NetCDF file.

    Input:
        fileName : string
           Name of the UMass FMCW data NetCDF file
    Output:
        A dictionary with the following arrays keyed by the array name:

        height : array
            Height (m AGL)
        vbins : array
            Spectral Velocity Bins (m/s)
        UTCtime : array
            Time in UTC (24-hr)
        Ze : array
            Logarithmic reflectivity factor (dBZ)
        vel : array
            Mean Radial Velocity (m/s)
        wid : array
            Spectrum Width (m/s)
        snr : array
            Signal-to-noise ratio (dB)
        PRF : float
            Pulse repetition frequency (Hz)
        Frequency : float
            Radar frequency (Hz)
        Wavelength : float
            Radar wavelength (m)
        Vmax : float
            Maximum unambiguous (Nyquist) velocity (m/s)
    """
    data = netcdf.Dataset(filename, "r")
    height = data.variables['height'][:]
    vbins = data.variables['vbins'][:]
    secs = data.variables['secs'][:]
    Ze = data.variables['Ze'][:]/10.
    ZeName = data.variables['Ze'].Name
    vel = data.variables['vel'][:]
    velName = data.variables['vel'].Name
    wid = data.variables['wid'][:]
    widName = data.variables['wid'].Name
    snr = data.variables['snr'][:]/10.
    snrName = data.variables['snr'].Name
    PRF = data.PRF
    Frequency = data.Frequency

    data.close()
    # Convert frequency to wavelength
    Wavelength = c/Frequency
    # Nyquist velocity
    Vmax = PRF * Wavelength / 4.
    # Convert secs (number of seconds since epoch) into UTC time
    UTCtime = []
    for s in secs:
        UTCtime.append(datetime(1970, 1, 1, 0, 0, 0) + timedelta(0, int(s)))
    UTCtime = np.array(UTCtime)
    return {'height': height, 'vbins': vbins, 'UTCtime': UTCtime, 'Ze': (ZeName, Ze),
            'vel': (velName, vel), 'wid': (widName, wid), 'snr': (snrName, snr),
            'PRF': PRF, 'Frequency': Frequency, 'Wavelength': Wavelength, 'Vmax': Vmax}


def read_fmcw_xarray(filename):
    """Reads an FMCW netCDF file into an xarray DataSet"""

    fmcw_dataset = xr.open_dataset(filename)
    # Rename a couple of the variables to the same name as the dimensions
    # and assign them as coordinates
    fmcw_dataset = fmcw_dataset.rename({'secs': 'time', 'gate': 'height'})
    fmcw_dataset = fmcw_dataset.set_coords(['time', 'height'])

    # Fix the 'units' attribute of the new time coordinate so that it can be
    # decoded into datetime objects
    fmcw_dataset['time'].attrs['units'] = fmcw_dataset['time'].attrs.pop('Units')
    fmcw_dataset['time'].attrs['units'] = 'seconds since 1-1-1970'
    fmcw_dataset = xr.decode_cf(fmcw_dataset)

    fmcw_dataset['Ze'] = fmcw_dataset['Ze']/10.
    fmcw_dataset['Zef'] = fmcw_dataset['Zef']/10.
    fmcw_dataset['snr'] = fmcw_dataset['snr']/10.

    # Add new Wavelength and Nyquist velocity attributes
    fmcw_dataset.attrs['Wavelength'] = c/fmcw_dataset.Frequency
    fmcw_dataset.attrs['Vmax'] = fmcw_dataset.PRF * fmcw_dataset.Wavelength / 4.

    return fmcw_dataset


def preprocess_fmcw_xarray(ds):
    ds = ds.rename({'secs': 'time', 'gate': 'height', 'vbins': 'vels'})
    ds = ds.set_coords(['time', 'height', 'vels'])
    # ds = ds.set_index(time='time', height='height', vels='vels')
    # print(ds)
    # ds = ds.set_index({'time': 'time', 'height': 'height', 'vels': 'vels'})
    # for dim in ds.dims:
    #     print(dim)
    #     indexes = ds.indexes.get(dim)
    #     print(indexes)
    # print(ds)
    return ds


def read_fmcw_multifile_xarray(filelist):
    """Reads a series of FMCW netCDF file into an xarray DataSet"""

    # For some reason, no matter what I do, combine='by_coords' refuses to work here. So we need
    # to sort the dataset by time after we open it.
    fmcw_dataset = xr.open_mfdataset(filelist, preprocess=preprocess_fmcw_xarray)
    fmcw_dataset = fmcw_dataset.sortby('time')
    # Rename a couple of the variables to the same name as the dimensions
    # and assign them as coordinates
    # fmcw_dataset = fmcw_dataset.rename({'secs': 'time', 'gate': 'height'})
    # fmcw_dataset = fmcw_dataset.set_coords(['time', 'height'])

    # Fix the 'units' attribute of the new time coordinate so that it can be
    # decoded into datetime objects
    fmcw_dataset['time'].attrs['units'] = fmcw_dataset['time'].attrs.pop('Units')
    fmcw_dataset['time'].attrs['units'] = 'seconds since 1-1-1970'
    fmcw_dataset = xr.decode_cf(fmcw_dataset)

    fmcw_dataset['Ze'] = fmcw_dataset['Ze']/10.
    fmcw_dataset['Zef'] = fmcw_dataset['Zef']/10.
    fmcw_dataset['snr'] = fmcw_dataset['snr']/10.

    # Add new Wavelength and Nyquist velocity attributes
    fmcw_dataset.attrs['Wavelength'] = c/fmcw_dataset.Frequency
    fmcw_dataset.attrs['Vmax'] = fmcw_dataset.PRF * fmcw_dataset.Wavelength / 4.

    return fmcw_dataset


def get_fmcw_filepaths(basedir, starttimestring, endtimestring, subdirs=False):
    """Gets a list of FMCW file paths given a base directory and a starttime and endtime"""
    timeformat = '%Y%m%d%H%M'

    starttime = datetime.strptime(starttimestring, timeformat)
    endtime = datetime.strptime(endtimestring, timeformat)

    delta = endtime - starttime
    totalhours = int(delta.days * 24 + delta.seconds/3600)

    delta_day_bounds = (endtime.date() - starttime.date()).days

    datetimerange = [starttime + timedelta(hours=x) for x in range(0, totalhours + 1)]
    daterange = [starttime.date() + timedelta(days=x) for x in range(0, delta_day_bounds + 1)]

    if subdirs:
        subdir_template = "{:02d}{:02d}"
        subdirs = [subdir_template.format(x.month, x.day) for x in daterange]

    filename_format = 'S{:04d}{:02d}{:02d}T{:02d}'
    filename_templates = [filename_format.format(x.year, x.month, x.day, x.hour)
                          for x in datetimerange]

    ncfiles2d = []
    if subdirs:
        for subdir in subdirs:
            absdir = os.path.join(basedir, subdir)
            ncfiles2d.append(glob(absdir + '/*nc'))
        ncfiles = list(itertools.chain.from_iterable(ncfiles2d))
    else:
        ncfiles = glob(basedir + '/*nc')

    pathlist = [f for f in ncfiles if any(m in f for m in filename_templates)]

    return pathlist


def plot_fmcw_4panel_xarray(dataset, figname='FMCW', dBZ_var='Ze'):
    """Basic four-panel plot of file contents."""

#     height, vbins, UTCtime, Ze, vel, wid, snr, PRF, Frequency, Wavelength, Vmax, \
#         ZeName, velName, widName, snrName = readfmcw_more(fileName)
    if dataset['time'].dt.year[0].item() == 2017: Vmax = 7.3  # 2017 data were collected with a
    # solid state amplifier that had a larger Nyquist interval than the default
    # +/-4.9 m/s.

    # 4-panel plot of output
    fig, ax = plt.subplots(2, 2, sharex=True, sharey=True)
    fig.set_size_inches(8.5,4.75)
    fig.subplots_adjust(bottom=0.15) # Keeps xlabel from getting chopped off
    fig.subplots_adjust(wspace=0.05) # Expand subplots slightly in the horizontal

    # Reflectivity (Ze)
#     m1 = ax[0, 0].pcolormesh(UTCtime, height/1000., Ze.T,cmap=plasma_white, \
#                              vmin=-30., vmax=30.)

    m1 = dataset[dBZ_var].plot(x='time', y='height', ax=ax[0, 0], cmap=plasma_white, vmin=-30., vmax=30.)

    ax[0, 0].set_title("(a) " + dataset[dBZ_var].name + " (dBZ)")
    ax[0, 0].set_ylim(dataset['height'].min(), dataset['height'].max())
    # ax[0, 0].xaxis.set_major_locator(md.MinuteLocator(byminute=range(0, 60, 10)))
    # ax[0, 0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
    ax[0, 0].yaxis.set_major_formatter(ticker.FuncFormatter(mtokm))
    #plt.xticks(rotation=30)
    #ax[0, 0].set_xlabel("Time (UTC)")
    ax[0, 0].set_ylabel("Height (km AGL)")
    # plt.colorbar(m1,ax=ax[0, 0])

    # Dealiased Doppler (vertical) velocity (vel_da)

    m2 = dataset['vel_da'].plot(x='time', y='height', ax=ax[0, 1], cmap=cm.RdBu_r,
                                vmin=-2.*np.ceil(dataset.Vmax), vmax=2.*np.ceil(dataset.Vmax))

#     m2 = ax[0, 1].pcolormesh(UTCtime,height/1000., vel.T, cmap=cm.RdBu_r, \
#                              vmin=-1.*np.ceil(Vmax), vmax=np.ceil(Vmax))

    ax[0, 1].set_title("(b) " + dataset['vel_da'].name + r" (m s$^{-1}$)")
    ax[0, 1].set_ylim(dataset['height'].min(), dataset['height'].max())
    # ax[0, 1].xaxis.set_major_locator(md.MinuteLocator(byminute=range(0, 60, 10)))
    # ax[0, 1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
    ax[0, 1].yaxis.set_major_formatter(ticker.FuncFormatter(mtokm))
    #plt.xticks(rotation=30)
    #ax[0, 1].set_xlabel("Time (UTC)")
    #ax[0, 1].set_ylabel("Height (km AGL)")
    # plt.colorbar(m2,ax=ax[0, 1])

    # Spectrum width (wid)

    m3 = dataset['wid'].plot(x='time', y='height', ax=ax[1, 0], cmap=cm.BuPu,
                             vmin=0., vmax=5.)

#     m3 = ax[1, 0].pcolormesh(UTCtime, height/1000., wid.T, cmap=cm.BuPu, \
#                              vmin=0., vmax=5.)
    ax[1, 0].set_title("(c) " + dataset['wid'].name + r" (m s$^{-1}$)")
    ax[1, 0].set_ylim(dataset['height'].min(), dataset['height'].max())
    # ax[1, 0].xaxis.set_major_locator(md.MinuteLocator(byminute=range(0, 60, 10)))
    # ax[1, 0].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
    ax[1, 0].yaxis.set_major_formatter(ticker.FuncFormatter(mtokm))
    #plt.xticks(rotation=30)
    ax[1, 0].set_xlabel("Time (UTC)")
    ax[1, 0].set_ylabel("Height (km AGL)")
    # plt.colorbar(m3,ax=ax[1, 0])

    # SNR (db)

    m4 = dataset['snr'].plot(x='time', y='height', ax=ax[1, 1], cmap=cm.viridis,
                             vmin=-20., vmax=30.)

#     m4 = ax[1, 1].pcolormesh(UTCtime, height/1000., snr.T, cmap=cm.viridis, \
#                              vmin=-20., vmax=30.)
    ax[1, 1].set_title("(d) " + dataset['snr'].name + " (dB)")
    ax[1, 1].set_ylim(dataset['height'].min(), dataset['height'].max())
    # ax[1, 1].xaxis.set_major_locator(md.MinuteLocator(byminute=range(0, 60, 10)))
    # ax[1, 1].xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
    ax[1, 1].yaxis.set_major_formatter(ticker.FuncFormatter(mtokm))
    #plt.xticks(rotation=30)
    ax[1, 1].set_xlabel("Time (UTC)")
    #ax[1, 1].set_ylabel("Height (km AGL)")
    # plt.colorbar(m4,ax=ax[1, 1])

    figname = figname + "_4panel.png"
    plt.savefig(figname, dpi=300)

    # Need to have ImageMagick installed on your system for the following command to
    # work. It deletes white space around the edges of the plot (thus saving disk space).
    try:
        os.system("convert -trim " + figname + " " + figname)
    except OSError:
        print("Can't find ImageMagick! No conversion will be done.")


def linear_texture_1d(image, n, axis=0, boundary_condition='reflect'):
    """
    Compute the 1D linear texture of an image given an n x ! window.
    This will only do the texture along one axis (useful for vertically
    pointing radar data)

    Parameters
    ----------
    image: 2D array of floats
        The image to calculate the texture of
    n: int
        The size of the window to calculate the standard deviation
        over.
    axis:
        Axis to calculate standard deviation over.
    boundary_condition:
        This determines how the edges are handled when calculating texture.
        'reflect' =

    Returns
    -------
    std_dev: 2D array of floats (n x n)
        The linear texture of the image

    Written by Bobby Jackson, Argonne National Laboratory, Feb 2018
    Tweaked by Robin Tanamachi, Purdue University, Feb 2018
    """

    mean_kernel = np.ones(n)/float(n)
    sum_kernel = np.ones(n)
    image = image.astype(float)
    image_squared = image**2
    sum_squares = ndimage.convolve1d(
        image_squared, sum_kernel, axis=axis, mode='reflect')
    sum_array = ndimage.convolve1d(
        image, sum_kernel, axis=axis, mode='reflect')
    # print(sum_array)

    variance = 1.0/float(n+1)*(sum_squares.astype(float) - sum_array.astype(float)**2/float(n))
    # print(variance)
    # print(np.max(variance),np.min(variance))

    return np.sqrt(variance)


def unwrap_fmcw_vel(vel, Vmax, wid_texture, snr, snrThreshold=-50.):
    """
    Dealiases Doppler velocities using the unwrap_phase method of scikit-image. The
    fields snr and wid_texture are combined to create a mask for the raw Doppler
    velocities; masked elements are not dealiased.

    Input:
        vel : 2D array of floats
            Doppler velocity array (vel) from an FMCW NetCDF file. (m/s)
        Vmax : float
            Maximum unambiguous velocity of the radar system. (m/s)
        wid_texture : 2D array of floats
            Output from linear_texture_1d applied to the wid field, used to construct
            a mask for the dealiasing. Masked elements won't be dealiased. (m/s)
        snr : 2D array of floats
            Signal-to-noise ratio array (snr) from an FMCW NetCDF file. (dB)
        snrThreshold : float
            Mask vel on snr below this value (in dB). Masked elements won't be
            dealiased.
    Output:
        vel_masked_unwrap : 2D array of floats (masked)
            Dealiased Doppler velocity array of the same size and shape as vel. (m/s)
    """
    wid_textureThreshold = 0.2  # Determined through experimentation. Anything above this
                                # value is likely to be nonmeteorological or attenuated.
    # Scale velocities to phase in range (-pi,pi]. Required by unwrap_phase.
    psi = vel / Vmax * np.pi
    # Combine SNR threshold and spectrum width texture threshold into a mask.
    # Note that the spectrum width texture mask is dilated by 1 pixel to mitigate edge
    # effects.
    psi_masked = np.ma.masked_array(data=psi, mask=np.any(np.dstack([snr <= snrThreshold,
                                    ndimage.morphology.binary_dilation(wid_texture > 0.2)]),
                                    axis=2))
    psi_masked_unwrap = unwrap_phase(psi_masked)
    # Convert phase back to velocity
    vel_masked_unwrap = psi_masked_unwrap * Vmax / np.pi

    return vel_masked_unwrap


def dealias(dataset, texture_kernel=3):
    """
    Given an FMCW xarray DataSet, dealias the velocity data and create a new variable to store it
    """
    # texture_kernel must be an odd, positive integer
    # Calculate spectrum width texture
    wid_texture = linear_texture_1d(dataset['wid'], texture_kernel, axis=1)
    # Mask the velocity on spectrum width and SNR, then dealias
    vel_masked_unwrap = unwrap_fmcw_vel(dataset['vel'], dataset.Vmax, wid_texture, dataset['snr'])
    dataset['vel_da'] = xr.DataArray(data=vel_masked_unwrap, dims=('time', 'height'),
                                     attrs={'name': 'Mean Radial Velocity (dealiased)',
                                            'units': 'm/s'})
    return dataset


def correct_fmcw_with_nexrad(fmcw_ds, PIPS_ds, radar_name='KHTX', plot=True):
    beam_height = PIPS_ds['{}_beam_height_at_PIPS'.format(radar_name)]
    beam_height = beam_height.interpolate_na(dim='time')
    beam_height = beam_height.ffill(dim='time')
    beam_height = beam_height.interp_like(fmcw_ds)
    beam_height = beam_height.ffill(dim='time')
    fmcw_ds_interp = interp_along_1D(fmcw_ds, beam_height, 'height', 'time')
    fmcw_dBZ_at_nexrad_beam = fmcw_ds_interp['Ze']
    nexrad_dBZ_at_fmcw = PIPS_ds['{}_at_PIPS'.format(radar_name)].sel(
        {'fields_{}'.format(radar_name): 'REF_filtered'})
    nexrad_dBZ_at_fmcw = nexrad_dBZ_at_fmcw.interp_like(fmcw_dBZ_at_nexrad_beam)

    if plot:
        fig1, ax1 = plt.subplots()
        fmcw_dBZ_at_nexrad_beam.plot(ax=ax1)
        nexrad_dBZ_at_fmcw.plot(ax=ax1)

    # Compute difference between FMCW and NEXRAD
    diff_dBZ = nexrad_dBZ_at_fmcw - fmcw_dBZ_at_nexrad_beam
    # Compute mean difference across time
    diff_dBZ_mean = diff_dBZ.mean(dim='time')

    if plot:
        fig2, ax2 = plt.subplots()
        diff_dBZ.plot(ax=ax2)
        ax2.axhline(y=diff_dBZ_mean)

    # Finally, correct the FMCW reflectivity
    fmcw_ds['Ze_corr'] = fmcw_ds['Ze'] + diff_dBZ_mean.values
    # set -99. values to missing in new Zef_corr variable
    fmcw_ds['Zef_corr'] = fmcw_ds['Zef'].copy()
    fmcw_ds['Zef_corr'] = fmcw_ds['Zef_corr'].where(fmcw_ds['Zef_corr'] > -99.)
    fmcw_ds['Zef_corr'] = fmcw_ds['Zef_corr'] + diff_dBZ_mean.values
    fmcw_ds['Ze_corr'].attrs['bias'] = diff_dBZ_mean.values
    fmcw_ds['Zef_corr'].attrs['bias'] = diff_dBZ_mean.values

    return fmcw_ds
