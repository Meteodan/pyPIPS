"""PIPS.py: a collection of functions to work with data from the PIPS and predecessors."""
import numpy as np
import pandas as pd
import xarray as xr
from . import thermolib as thermo
from . import utils
import pyPIPS.parsivel_params as pp
from .parsivel_params import parsivel_parameters
from .pips_io import combine_parsivel_data
from numba import jit

deg2rad = np.pi / 180.


min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
diameter_edges = np.append(min_diameter, max_diameter[-1])
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']
fall_bins_edges = np.append(min_fall_bins, max_fall_bins[-1])


def calc_thermo(conv_df):
    """[summary]

    Parameters
    ----------
    conv_df : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    p_Pa = conv_df['pressure'] * 100.
    T_K = conv_df['fasttemp'] + 273.15
    RH = conv_df['RH_derived'] / 100.

    pt = thermo.caltheta(p_Pa, T_K)
    qv = thermo.calqv(RH, p_Pa, T_K)
    rho = thermo.calrho(p_Pa, pt, qv)

    conv_df['pt'] = pt
    conv_df['qv'] = qv
    conv_df['rho'] = rho

    return conv_df


def get_offset_seconds(datetime_range):
    """[summary]

    Parameters
    ----------
    datetime_range : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    return datetime_range[0].second


# TODO: Move this somewhere else?
def wind_dir_and_speed_from_u_and_v(u, v):
    """[summary]

    Arguments:
        u {[type]} -- [description]
        v {[type]} -- [description]
    """
    windspdvec = np.sqrt(u**2. + v**2.)
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    winddirvec = (270.0 - (180. / np.pi) * np.arctan2(v, u)) % 360.

    return windspdvec, winddirvec


# TODO: not used right now, but will eventually refactor everything to use xarray DataArrays
# instead of pandas DataFrames
def resample_wind_da(wind_dir, wind_spd, intervalstr, offset, gusts=True, gustintervalstr='3S',
                     center=False):
    wind_spd_avg = wind_spd.resample(time=intervalstr, label='right', closed='right',
                                     base=offset).mean()
    if gusts:
        wind_gust = wind_spd.resample(time=gustintervalstr, label='right', closed='right',
                                      base=offset).mean()
        wind_gust_avg = wind_gust.resample(time=intervalstr, label='right', closed='right',
                                           base=offset).max()
    else:
        wind_gust = None
        wind_gust_avg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    u = wind_spd * np.cos(np.deg2rad(-wind_dir + 270.))
    v = wind_spd * np.sin(np.deg2rad(-wind_dir + 270.))

    # Linearly interpolate for bad values of us,vs
    # u = interpnan1D(u)
    # v = interpnan1D(v)
    # Compute averages of wind components
    u_avg = u.resample(time=intervalstr, label='right', closed='right', base=offset).mean()
    v_avg = v.resample(time=intervalstr, label='right', closed='right', base=offset).mean()

    wind_spd_vec_avg = np.sqrt(u_avg**2. + v_avg**2.)
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    wind_dir_vec_avg = (270.0 - (180. / np.pi) * np.arctan2(v_avg, u_avg)) % 360.

    # unit average wind direction
    unit_u = u / wind_spd
    unit_v = v / wind_spd
    unit_u_avg = unit_u.resample(time=intervalstr, label='right', closed='right',
                                 base=offset).mean()
    unit_v_avg = unit_v.resample(time=intervalstr, label='right', closed='right',
                                 base=offset).mean()

    # Need to use %360 to keep wind dir between 0 and 360 degrees
    wind_dir_unit_vec_avg = (270.0 - (180. / np.pi) * np.arctan2(unit_u_avg, unit_v_avg)) % 360.

    # Pack everything into a dictionary and return
    wind_dict = {'windspd': wind_spd_avg, 'windspdavgvec': wind_spd_vec_avg,
                 'winddirabs': wind_dir_vec_avg, 'winddirunitavgvec': wind_dir_unit_vec_avg,
                 'windgust': wind_gust_avg, 'uavg': u_avg, 'vavg': v_avg,
                 'unit_uavg': unit_u_avg, 'unit_vavg': unit_v_avg}

    return wind_dict


def resample_wind(datetimes, offset, winddirs, windspds, intervalstr, gusts=True, gustintvstr='3S',
                  center=False):
    """Given a timeseries of wind directions and speeds, and an interval for resampling,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""

    windspdsavg = pd.Series(data=windspds, index=datetimes).resample(intervalstr, label='right',
                                                                     closed='right',
                                                                     base=offset).mean()
    if gusts:
        windgusts = pd.Series(data=windspds, index=datetimes).resample(gustintvstr, label='right',
                                                                       closed='right',
                                                                       base=offset).mean()
        windgustsavg = windgusts.resample(intervalstr, label='right', closed='right',
                                          base=offset).max()
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds * np.cos(np.deg2rad(-winddirs + 270.))
    vs = windspds * np.sin(np.deg2rad(-winddirs + 270.))

    # Linearly interpolate for bad values of us,vs
    # us = interpnan1D(us)
    # vs = interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.Series(
        data=us,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    vsavg = pd.Series(
        data=vs,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()

    windspdsavgvec = np.sqrt(usavg**2. + vsavg**2.)
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    winddirsavgvec = (270.0 - (180. / np.pi) * np.arctan2(vsavg, usavg)) % 360.

    # unit average wind direction
    unit_us = us / windspds
    unit_vs = vs / windspds
    unit_usavg = pd.Series(
        data=unit_us,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    unit_vsavg = pd.Series(
        data=unit_vs,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    # Need to use %360 to keep wind dir between 0 and 360 degrees
    winddirsunitavgvec = (270.0 - (180. / np.pi) * np.arctan2(unit_vsavg, unit_usavg)) % 360.

    # Pack everything into a dictionary and then into a pd.DataFrame
    wind_dict = {'windspd': windspdsavg, 'windspdavgvec': windspdsavgvec,
                 'winddirabs': winddirsavgvec, 'winddirunitavgvec': winddirsunitavgvec,
                 'windgust': windgustsavg, 'uavg': usavg, 'vavg': vsavg,
                 'unit_uavg': unit_usavg, 'unit_vavg': unit_vsavg}

    wind_df = pd.DataFrame(wind_dict)

    return wind_df


def resample_compass(datetimes, compass_dir, offset, intervalstr):

    # Compute x- and y-components of compass direction
    x = np.cos(np.deg2rad(-compass_dir + 270.))
    y = np.sin(np.deg2rad(-compass_dir + 270.))

    # Compute averages of components
    x_avg = pd.Series(
        data=x,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()
    y_avg = pd.Series(
        data=y,
        index=datetimes).resample(
        intervalstr,
        label='right',
        closed='right',
        base=offset).mean()

    # Need to use %360 to keep direction between 0 and 360 degrees
    compass_dir_avg = (270.0 - (180. / np.pi) * np.arctan2(y_avg, x_avg)) % 360.
    return compass_dir_avg


def resample_conv(probe_type, resample_interval, sec_offset, conv_df, gusts=False, gustintvstr='3S',
                  center=False):
    """Resamples the conventional data to a longer interval"""

    if probe_type == 'PIPS':
        winddirkey = 'winddirabs'
        windspdkey = 'windspd'
        other_data = ['fasttemp', 'slowtemp', 'RH', 'RH_derived', 'pressure',
                      'GPS_lat', 'GPS_lon', 'GPS_alt', 'voltage', 'dewpoint', 'rho']
    elif probe_type == 'TriPIPS':
        winddirkey = 'winddirabs'
        windspdkey = 'windspd'
        other_data = ['slowtemp', 'RH', 'RH_derived', 'pressure',
                      'GPS_lat', 'GPS_lon', 'GPS_alt', 'voltage', 'dewpoint', 'rho']
    elif probe_type == 'CU':
        winddirkey = 'bwinddirabs'
        windspdkey = 'bwindspd'
        other_data = ['slowtemp', 'RH', 'pressure',
                      'GPS_lat', 'GPS_lon', 'GPS_alt', 'voltage', 'dewpoint', 'rho']
    elif probe_type == 'NV2':
        winddirkey = 'swinddirabs'
        windspdkey = 'swindspd'
        other_data = ['slowtemp', 'RH', 'pressure', 'dewpoint', 'rho']

    intervalstr = '{:d}S'.format(int(resample_interval))

    # First, resample the winds

    conv_resampled_df = resample_wind(conv_df.index, sec_offset, conv_df[winddirkey],
                                      conv_df[windspdkey], intervalstr, gusts=gusts,
                                      gustintvstr=gustintvstr, center=center)

    # Resample compass direction
    if probe_type == 'PIPS':
        compass_dir_avg = resample_compass(conv_df.index, conv_df['compass_dir'], sec_offset,
                                           intervalstr)
        conv_resampled_df['compass_dir'] = compass_dir_avg

    # Special treatment for wind diagnostic flags
    # Note, have to use numpy max() as a lambda function because the
    # pandas resample.max() does not propagate NaN's!
    # TODO: Find out if this is still true since pandas has been updated several times since then
    if probe_type == 'PIPS':
        winddiags_rs = pd.Series(data=conv_df['winddiag'], index=conv_df.index).resample(
            intervalstr, label='right', closed='right',
            base=sec_offset).apply(lambda x: utils.trymax(x.values))

        conv_resampled_df['winddiag'] = winddiags_rs

    # Now, resample the other one-sec data
    other_resampled_df = \
        conv_df[other_data].resample(intervalstr, label='right', closed='right',
                                     base=sec_offset).mean()

    conv_resampled_df = conv_resampled_df.join(other_resampled_df)

    return conv_resampled_df


def check_requested_resampling_interval(requested_interval, sampling_interval):
    """[summary]

    Parameters
    ----------
    requested_interval : [type]
        [description]
    sampling_interval : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    resample_index_interval = int(requested_interval / sampling_interval)
    print("Requested DSD interval: {:.1f}. Actual DSD interval: {:.1f}".format(
        requested_interval, resample_index_interval * sampling_interval))
    resample_interval = resample_index_interval * sampling_interval

    return resample_interval


def resample_vd_matrix(resample_interval, vd_matrix):
    """[summary]

    Parameters
    ----------
    resample_interval : [type]
        [description]
    vd_matrix : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    intervalstr = '{:d}S'.format(int(resample_interval))

    # We need to find the offset corresponding to the starting second and then
    # generate the frequency string. Seems like there should be an easier way...
    sec_offset = pd.to_datetime(vd_matrix['time'].values)[0].second
    # Resample the vd_matrix in time, filling missing values with zero
    vd_matrix = vd_matrix.resample(time=intervalstr, label='right', closed='right',
                                   base=sec_offset).sum(dim='time').fillna(0)

    return vd_matrix


def resample_parsivel(resample_interval, parsivel_df):
    """[summary]

    Parameters
    ----------
    resample_interval : [type]
        [description]
    parsivel_df : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    intervalstr = '{:d}S'.format(int(resample_interval))
    # We need to find the offset corresponding to the starting second and then
    # generate the frequency string. Seems like there should be an easier way...
    sec_offset = parsivel_df.index.to_pydatetime()[0].second
    # Resample parsivel_df in time, filling missing values with zero: TEMPORARY FIX. later will deal
    # with missing values differently
    parsivel_rs = parsivel_df.resample(intervalstr, label='right', closed='right',
                                       base=sec_offset)

    # Each column of the dataframe needs to downsample differently. For example, the
    # precipitation totals need to be summed, while the reflectivity should be averaged,
    # etc. We can do this with the Pandas "agg" function on the resampler. Neat, huh?
    # parsivel_df = parsivel_rs.agg({'intensity': np.mean, 'preciptot': np.sum,
    #                                'reflectivity': np.mean, 'pcount': np.sum, 'pcount2': np.sum,
    #                                'amplitude': np.mean, 'flaggedtimes': np.any,
    #                                'hailflag': np.any}).fillna(0)

    parsivel_df = parsivel_rs.agg({'precipintensity': np.mean, 'precipaccum': np.sum,
                                   'parsivel_dBZ': np.mean, 'pcount': np.sum,
                                   'signal_amplitude': np.mean, 'pvoltage': np.mean,
                                   'sensor_temp': np.mean,
                                   'sample_interval': np.mean}).fillna(0)

    return parsivel_df


def resample_ND(resample_interval, ND_df):
    """[summary]

    Parameters
    ----------
    resample_interval : [type]
        [description]
    ND_df : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    intervalstr = '{:d}S'.format(int(resample_interval))
    # We need to find the offset corresponding to the starting second and then
    # generate the frequency string. Seems like there should be an easier way...
    sec_offset = ND_df.index.to_pydatetime()[0].second
    # Resample parsivel_df in time, filling missing values with zero: TEMPORARY FIX. later will deal
    # with missing values differently
    ND_df = ND_df.resample(intervalstr, label='right', closed='right',
                           base=sec_offset).mean().fillna(0)

    return ND_df


def calc_ND(vd_matrix, fallspeed_spectrum, sample_interval, use_measured_fallspeed=True):
    """Computes the binned number density from the disdrometer (m^-3 mm^-1)

    Parameters
    ----------
    vd_matrix : array_like
        Velocity-diameter count matrix from PIPS (32 x 32)
    fallspeed_spectrum : array_like
        Velocity spectrum
    sample_interval : array_like
        Time interval of DSD integration

    Returns
    -------
    array_like
        Binned number density (m^-3 mm^-1)
    """
    eff_sensor_area = xr.DataArray(pp.parsivel_parameters['eff_sensor_area_mm2'],
                                   dims=['diameter_bin']) * 1.e-6  # Get to m2
    bin_width = xr.DataArray((pp.parsivel_parameters['max_diameter_bins_mm'] -
                              pp.parsivel_parameters['min_diameter_bins_mm']), dims=['diameter_bin'])

    if not use_measured_fallspeed:
        vd_matrix = vd_matrix.sum(dim='fallspeed_bin')
    ND = vd_matrix / (fallspeed_spectrum * sample_interval * eff_sensor_area * bin_width)
    if use_measured_fallspeed:
        ND = ND.sum(dim='fallspeed_bin')

    return ND


def calc_ND_onedrop(sample_interval, correct_rho=False, rho=None):

    diameter_bins = pp.parsivel_parameters['avg_diameter_bins_mm']
    fallspeed_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

    fallspeed_spectrum = calc_fallspeed_spectrum(diameter_bins, fallspeed_bins,
                                                 correct_rho=correct_rho, rho=rho,
                                                 use_measured_fallspeed=False)
    if 'time' in rho.dims:
        vd_matrix_onedrop = np.ones((rho.size, 1, diameter_bins.size))
        coords = {
            'time': rho['time'],
            'diameter': ('diameter_bin', diameter_bins),
        }
        dims = ['time', 'fallspeed_bin', 'diameter_bin']
    else:
        vd_matrix_onedrop = np.ones((1, diameter_bins.size))
        coords = {
            'diameter': ('diameter_bin', diameter_bins)
        }
        dims = ['fallspeed_bin', 'diameter_bin']
    vd_matrix_onedrop_da = \
        xr.DataArray(vd_matrix_onedrop,
                     name='velocity_diameter_onedrop_per_bin',
                     coords=coords,
                     dims=dims)

    ND_onedrop = calc_ND(vd_matrix_onedrop_da, fallspeed_spectrum, sample_interval)
    return ND_onedrop


def calc_fallspeed_spectrum(diameter_bins, fallspeed_bins,
                            correct_rho=False, rho=None, use_measured_fallspeed=True):
    """[summary]

    Parameters
    ----------
    diameter_bins : [type]
        [description]
    fallspeed_bins : [type]
        [description]
    correct_rho : bool, optional
        [description], by default False
    rho : [type], optional
        [description], by default None
    use_measured_fallspeed : bool, optional
        [description], by default True

    Returns
    -------
    [type]
        [description]
    """
    if not use_measured_fallspeed:
        fallspeed_spectrum = calc_empirical_fallspeed(diameter_bins, correct_rho=correct_rho,
                                                      rho=rho)
        # print('rho', rho)
        # exit
        if 'time' in rho.dims:
            coords = {
                'time': rho['time'],
                'diameter': ('diameter_bin', diameter_bins),
            }
            dims = ['time', 'diameter_bin']
        else:
            coords = {
                'diameter': ('diameter_bin', diameter_bins),
            }
            dims = ['diameter_bin']
        fallspeed_da = xr.DataArray(fallspeed_spectrum, coords=coords, dims=dims)
    else:
        _, fallspeed_spectrum = np.meshgrid(diameter_bins, fallspeed_bins)
        # TODO: STOPPED HERE! 09/03/19. Need to fix the dimension names here to match
        # up with those of the spectrum dataarray, so that broadcasting works
        # NOTE: (01/28/2020) I think I must have fixed this?
        fallspeed_da = xr.DataArray(fallspeed_spectrum,
                                    coords={
                                        'fallspeed': ('fallspeed_bin', fallspeed_bins),
                                        'diameter': ('diameter_bin', diameter_bins)
                                    },
                                    dims=['fallspeed_bin', 'diameter_bin'])

    return fallspeed_da


@jit
def calc_empirical_fallspeed(d, correct_rho=False, rho=None):
    """Assigns a fall speed for a range of diameters based on code
       from David Dowell (originally from Terry Schuur).  It appears that
       the formulas originate from Atlas et al. (1973), but this took a bit of sleuthing!"""

    v = np.where(d < 3.0, 3.78 * d**0.67, 9.65 - 10.3 * np.exp(-0.6 * d))
    v = np.where(d < 0.0, 0.0, v)

    # Correct fall speed based on air density
    # Based on Foote and duToit (1969): v = v0*(rho0/rho)^(0.4)
    # where rho0 = 1.204 kg/m^3 -- that corresponding to a T of 20 C and pressure of 1013 mb.

    if correct_rho and rho is not None:
        v = v[:, None] * (1.204 / rho.values)**(0.4)
        v = v.squeeze()
        v = np.atleast_1d(v)
        v = v.T
    return v


def avgwind(winddirs, windspds, avgintv, gusts=True, gustintv=3, center=True):
    """Given a timeseries of wind directions and speeds, and an interval for averaging,
       compute the vector and scalar average wind speed, and vector average wind direction.
       Optionally also compute gusts."""

    windspdsavg = pd.Series(windspds).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    if gusts:
        windgusts = pd.Series(windspds).rolling(
            window=gustintv,
            center=center,
            min_periods=1).mean().values
        windgustsavg = pd.Series(windgusts).rolling(
            window=avgintv,
            center=center,
            min_periods=1).max().values
    else:
        windgusts = None
        windgustsavg = None

    # Compute vector average wind speed and direction
    # First compute the u and v wind components
    us = windspds * np.cos(np.deg2rad(-winddirs + 270.))
    vs = windspds * np.sin(np.deg2rad(-winddirs + 270.))

    # Linearly interpolate for bad values of us,vs
    us = utils.interpnan1D(us)
    vs = utils.interpnan1D(vs)
    # Compute averages of wind components
    usavg = pd.Series(us).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    vsavg = pd.Series(vs).rolling(
        window=avgintv,
        center=center,
        min_periods=1).mean().values
    windspdsavgvec = np.sqrt(usavg**2. + vsavg**2.)
    winddirsavgvec = (270.0 - (180. / np.pi) * np.arctan2(vsavg, usavg)
                      ) % 360.  # Need to use %360 to keep wind dir between 0 and 360 degrees

    # TODO: get all this in a pandas DataFrame once all calling routines have been checked.

    return windspdsavg, windspdsavgvec, winddirsavgvec, windgusts, windgustsavg


def get_PSD_datetimes(vd_matrix, dim_name='time'):
    return pd.to_datetime(vd_matrix[dim_name].values).to_pydatetime()


def get_datetimes(PIPS_ds, dim_name='time'):
    return pd.to_datetime(PIPS_ds[dim_name].values).to_pydatetime()


def get_conv_datetimes(conv_df, dim_name='time'):
    return conv_df.index.to_pydatetime()


def get_PSD_time_bins(PSD_datetimes):
    PSD_interval_td = PSD_datetimes[1] - PSD_datetimes[0]
    PSD_half_interval_td = PSD_interval_td / 2.

    PSD_datetimes_edges = PSD_datetimes - PSD_interval_td
    last_edge = np.array(PSD_datetimes_edges[-1] + PSD_interval_td)
    PSD_datetimes_edges = np.append(PSD_datetimes_edges, last_edge)
    PSD_datetimes_centers = PSD_datetimes - PSD_half_interval_td

    return {'PSD_datetimes': PSD_datetimes, 'PSD_datetimes_edges': PSD_datetimes_edges,
            'PSD_datetimes_centers': PSD_datetimes_centers}


def get_interval_str(interval):
    return '{:d}S'.format(int(interval))


def calc_parsivel_wind_angle(wind_dir, compass_dir, parsivel_angle):
    parsivel_true = np.mod(compass_dir + parsivel_angle, 360.)
    parsivel_wind_diff = wind_dir - parsivel_true
    parsivel_wind_diff = np.abs(np.rad2deg(np.arcsin(np.sin(np.deg2rad(parsivel_wind_diff)))))
    return parsivel_wind_diff


def reindex_velocity_bins(VD_matrix, interval):
    # Set up new regularly spaced bin edges
    min_fallspeed = VD_matrix.coords['min_fallspeeds'][0]
    max_fallspeed = VD_matrix.coords['max_fallspeeds'][-1]
    new_vt_bin_edges = np.arange(min_fallspeed, max_fallspeed+interval, interval)
    new_min_vt_bins = np.delete(new_vt_bin_edges, -1)
    new_max_vt_bins = np.delete(new_vt_bin_edges, 0)
    new_center_vt_bins = 0.5 * (new_min_vt_bins + new_max_vt_bins)
    # Scale velocity counts by velocity bin width
    weights = VD_matrix['max_fallspeeds'] - VD_matrix['min_fallspeeds']
    VD_matrix_scaled = VD_matrix / weights
    # Set the minimum velocity for each bin as the new index. Needed so that the reindexing will
    # properly center the scaled counts in each sub-bin.
    VD_matrix_scaled = VD_matrix_scaled.set_index(fallspeed_bin='min_fallspeeds')
    # Now re-bin the scaled velocity counts into the new regularly spaced bins
    VD_matrix_scaled_rebinned = VD_matrix_scaled.reindex({'fallspeed_bin': new_min_vt_bins},
                                                         method='pad')
    # Recover original counts. NOTE: in general there may now be fractional counts in the new bins,
    # but the total number of drops is preserved by this procedure.
    VD_matrix_rebinned = VD_matrix_scaled_rebinned * interval
    # Set index back to fallspeed bin centers and add min_fallspeed coordinate back
    VD_matrix_rebinned.coords['fallspeed'] = ('fallspeed_bin', new_center_vt_bins)
    VD_matrix_rebinned = VD_matrix_rebinned.set_index(fallspeed_bin='fallspeed')
    VD_matrix_rebinned.coords['fallspeed'] = ('fallspeed_bin', new_center_vt_bins)
    VD_matrix_rebinned.coords['min_fallspeeds'] = ('fallspeed_bin', new_min_vt_bins)
    VD_matrix_rebinned.coords['max_fallspeeds'] = ('fallspeed_bin', new_max_vt_bins)
    # VD_matrix_rebinned.coords['diameter_bin'] = ('diameter_bin', VD_matrix_rebinned['diameter'])
    return VD_matrix_rebinned


def calc_mean_velocity(vd_matrix_rebinned):
    # Need to iterate across diameter dimension, which is annoying.
    mean_vels = []
    for vel_da in vd_matrix_rebinned.transpose():
        weights = vel_da.fillna(0)
        velocities = vel_da.coords['fallspeed']
        vel_weighted = velocities.weighted(weights)
        mean_vel_d = vel_weighted.mean('fallspeed_bin')
        mean_vels.append(mean_vel_d)
    mean_vel = xr.concat(mean_vels, dim='diameter_bin')
    return mean_vel


def shift_mean_velocity(vd_matrix_rebinned, vt_rain):
    # We have to iterate over both time and diameter dimensions, because the shift needed is a
    # function of both
    # TODO: This is painfully slow and memory-intensive for large timeseries of DSDs.
    # Need to optimize!
    vel_interval = (vd_matrix_rebinned['fallspeed_bin'][1] -
                    vd_matrix_rebinned['fallspeed_bin'][0].values)
    # print(vel_interval)
    mean_vel = calc_mean_velocity(vd_matrix_rebinned)
    time_flag = 'time' in vd_matrix_rebinned.dims
    if time_flag:
        vt_rain = vt_rain.stack(stack_dim=['time', 'diameter_bin'])
        mean_vel = mean_vel.stack(stack_dim=['time', 'diameter_bin'])
        vd_matrix_rebinned = vd_matrix_rebinned.stack(stack_dim=['time', 'diameter_bin'])
        groupby_dims = 'stack_dim'
    else:
        groupby_dims = 'diameter_bin'
    vt_d_groups = vt_rain.groupby(groupby_dims, squeeze=False)
    mean_vel_d_groups = mean_vel.groupby(groupby_dims, squeeze=False)
    vel_da_groups = vd_matrix_rebinned.groupby(groupby_dims, squeeze=False)
    if time_flag:
        vd_matrix_rebinned = vd_matrix_rebinned.unstack('stack_dim')
    for vt_l, mean_vel_l, vel_da_l in zip(list(vt_d_groups), list(mean_vel_d_groups),
                                          list(vel_da_groups)):
        vt_d = vt_l[1]
        mean_vel_d = mean_vel_l[1]
        vel_da = vel_da_l[1]
        vt_diff = vt_d - mean_vel_d
        if np.isfinite(vt_diff.values):
            vt_shift = int(vt_diff / vel_interval)
        else:
            vt_shift = 0
        vel_da = vel_da.shift(fallspeed_bin=vt_shift)
        if time_flag:
            vel_da = vel_da.unstack('stack_dim')
            vd_matrix_rebinned.loc[dict(time=vel_da['time'].values,
                                        diameter_bin=vel_da['diameter_bin'].values)] = vel_da
        else:
            vd_matrix_rebinned.loc[dict(diameter_bin=vel_da['diameter_bin'])] = vel_da
    return vd_matrix_rebinned


def rebin_to_parsivel(vd_matrix_rebinned):
    # Regroup into original fallspeed bins and get the dimension and coordinate names back in order.
    # This is rather clunky, but it works
    vd_matrix_groups = vd_matrix_rebinned.groupby_bins('fallspeed_bin', fall_bins_edges)
    vd_matrix = vd_matrix_groups.sum()
    vd_matrix = vd_matrix.swap_dims({'fallspeed_bin_bins': 'fallspeed_bin'})
    vd_matrix = vd_matrix.rename({'fallspeed_bin_bins': 'fallspeed_bin'})
    vd_matrix = vd_matrix.reindex({'fallspeed_bin': avg_fall_bins})
    vd_matrix.coords['min_fallspeeds'] = ('fallspeed_bin', min_fall_bins)
    vd_matrix.coords['max_fallspeeds'] = ('fallspeed_bin', max_fall_bins)
    vd_matrix.coords['fallspeed'] = ('fallspeed_bin', avg_fall_bins)
    # Get the dimensions back into their original order
    vd_matrix = vd_matrix.transpose('time', 'fallspeed_bin', 'diameter_bin')
    return vd_matrix


def correct_ND_RB15(parsivel_ds, ND_name='ND_RB15_vshift_qc'):
    ND = parsivel_ds[ND_name]
    ND_RB15 = ND.copy()
    RR_edges = np.append(pp.RB15_RR_min, pp.RB15_RR_max[-1])
    # According to RB15, should use internal parsivel rainrate ('precipintensity') to categorize
    parsivel_ds_RR_groups = parsivel_ds.groupby_bins('precipintensity', RR_edges, right=False,
                                                     labels=pp.RB15_RR_min)
    # For some reason the group labels aren't sorted in increasing order...

    # Get indices in original array of each member of each group
    group_indices = parsivel_ds_RR_groups.groups
    # Perform correction to ND using the appropriate correction factors for each rainrate category
    for group in parsivel_ds_RR_groups:
        group_idx = group_indices[group[0]]
        ND_temp = ND[group_idx]
        RB15_correction_factor = pp.RB15_correction_factors.sel(rainrate=group[0])
        ND_temp2 = RB15_correction_factor * ND_temp
        ND_temp2 = ND_temp2.transpose()
        ND_RB15[group_idx] = ND_temp2

    # Add the corrected ND to the parsivel_ds
    return combine_parsivel_data(parsivel_ds, ND_RB15, name='ND_RB15_qc')


def calc_stats(ds, var_x, var_y):
    bias = (100. * (ds[var_y] - ds[var_x]).mean() / ds[var_x].mean()).values
    cc = pd.DataFrame({'x': ds[var_x], 'y': ds[var_y]}).corr()
    return cc, bias
