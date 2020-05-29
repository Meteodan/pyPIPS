"""PIPS.py: a collection of functions to work with data from the PIPS and predecessors."""
import numpy as np
import pandas as pd
import xarray as xr
from . import obanmodule as oban
from . import radarmodule as radar
from . import thermolib as thermo
from . import utils
from .parsivel_params import parsivel_parameters
from numba import jit

deg2rad = np.pi / 180.


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
        windgusts = windgusts
        windgustsavg = windgustsavg
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
                                   'signal_amplitude': np.mean, 'pvoltage': np.mean, 'sensor_temp': np.mean,
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
    eff_sensor_area = xr.DataArray(parsivel_parameters['eff_sensor_area_mm2'],
                                   dims=['diameter_bin']) * 1.e-6  # Get to m2
    bin_width = xr.DataArray((parsivel_parameters['max_diameter_bins_mm'] -
                              parsivel_parameters['min_diameter_bins_mm']), dims=['diameter_bin'])

    if not use_measured_fallspeed:
        vd_matrix = vd_matrix.sum(dim='fallspeed_bin')
    ND = vd_matrix / (fallspeed_spectrum * sample_interval * eff_sensor_area * bin_width)
    if use_measured_fallspeed:
        ND = ND.sum(dim='fallspeed_bin')

    return ND


def calc_ND_onedrop(sample_interval, correct_rho=False, rho=None):

    diameter_bins = parsivel_parameters['avg_diameter_bins_mm']
    fallspeed_bins = parsivel_parameters['avg_fallspeed_bins_mps']

    fallspeed_spectrum = calc_fallspeed_spectrum(diameter_bins, fallspeed_bins,
                                                 correct_rho=correct_rho, rho=rho,
                                                 use_measured_fallspeed=False)
    vd_matrix_onedrop = np.ones((rho.size, 1, diameter_bins.size))
    vd_matrix_onedrop_da = \
        xr.DataArray(vd_matrix_onedrop,
                     name='velocity_diameter_onedrop_per_bin',
                     coords={
                         'time': rho['time'],
                         'diameter': ('diameter_bin', diameter_bins)
                     },
                     dims=['time', 'fallspeed_bin', 'diameter_bin'])

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
        fallspeed_da = xr.DataArray(fallspeed_spectrum,
                                    coords={
                                        'time': rho['time'],
                                        'diameter': ('diameter_bin', diameter_bins),
                                    },
                                    dims=['time', 'diameter_bin'])
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

@jit(parallel=True)
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


# TODO: refactor this function to use xarray
def rad2DD2(fieldlist, range_start, rrange, azimuth_start_rad, azimuth_rad, rlat, rlon, ralt, el,
            dlocs, average_gates=True, Cressman=False, roi=750., map_proj=1):
    """
    Another version of rad2DD: assumes radar sweep has been read in and computes values of fields in
    fieldlist at the disdrometer location(s).  Eventually will allow for optional
    advection/sedimentation correction of radar reflectivity to account for height above surface of
    radar scan.  Returns list of lists of field values at each disdrometer location
    """

    # First find x,y locations of radar gates

    xrad, yrad, xrad_c, yrad_c = radar.sweep2xy(
        azimuth_start_rad, azimuth_rad, range_start, rrange, el, rlat, rlon, ralt, map_proj)

    field_D_list = []
    dxy_list = []

    for dloc in dlocs:
        Dx, Dy = oban.ll_to_xy(dloc[0] * deg2rad, dloc[1] * deg2rad, rlat, rlon, map_proj)
        dxy_list.append((Dx, Dy))
        # print "Dx,Dy",Dx,Dy

        # Find closest gate to disdrometer location
        # First find the closest azimuth

#         theta_dis = -np.arctan2(Dy,Dx)+np.pi/2.       # I *think* this is correct
#
#         #print "azimuth of disdrometer: ",theta_dis/deg2rad
#         #Find closest index of azimuth
#
#         theta_diff = theta_dis-theta_c[0,:]
#         theta_index = np.argmin(np.abs(theta_diff))
#
#         # Now find index along range; first compute slant range
#
#         srange_dis = oban.computeslantrange([Dx],[Dy],el)
#         rad_diff = srange_dis-rad_c[:,0]
#         srange_index = np.argmin(np.abs(rad_diff))
#         #print "range of disdrometer: ",srange_dis

        # Try another way to get the indices (the above seems to have problems).
        # EDIT: 03/18/2012 -- the following way is more reliable but is also slower.  Not sure why
        # the above doesn't work in all situations, but it sometimes picks a gate adjacent to the
        # onewe want...

        field_list = []

        for field in fieldlist:
            field = field.swapaxes(0, 1)
            field = np.ma.masked_invalid(field)

            if Cressman:
                # print xplt_c.shape,yplt_c.shape
                field_D = oban.Cresmn(xrad_c, yrad_c, field, Dx, Dy, roi)
                # print dBZ_D,dBZ_D.shape
                # First, compute the euclidian distance of the x,y location of the disdrometer to
                # each of the radar gate centers
            else:
                distance = np.sqrt(((Dx - xrad_c)**2. + (Dy - yrad_c)**2.))

                # Now, find the index of the closest radar gate
                srange_index, theta_index = np.unravel_index(distance.argmin(), distance.shape)
                # print "srange_index,theta_index",srange_index,theta_index

                print("Distance to closest gate: ", distance[srange_index, theta_index])

                # Finally, grab field at closest gate to disdrometer
                # Average the field in the closest gate and surrounding 8 gates if desired
                # Disdrometer is probably not covered by sweep, set value to np.nan
                if(distance[srange_index, theta_index] > 3000.):
                    field_D = np.nan
                else:
                    if not average_gates:
                        field_D = field[srange_index, theta_index]
                    else:
                        # field_D = (1./9.)*(field[srange_index-1,theta_index-1]+
                        #                    field[srange_index,theta_index-1]+
                        #                    field[srange_index-1,theta_index]+
                        #                    field[srange_index,theta_index]+
                        #                    field[srange_index+1,theta_index]+
                        #                    field[srange_index,theta_index+1]+
                        #                    field[srange_index+1,theta_index+1]+
                        #                    field[srange_index+1,theta_index-1]+
                        #                    field[srange_index-1,theta_index+1])
                        field_D = np.nanmean(field[srange_index - 1:srange_index + 2,
                                                   theta_index - 1:theta_index + 2])

            # print "Value of field at disdrometer = ","%.2f"%field_D

            field_list.append(field_D)

        field_D_list.append(field_list)

    field_D_arr = np.array(field_D_list)

    return dxy_list, field_D_arr


def get_PSD_datetimes(vd_matrix, dim_name='time'):
    return pd.to_datetime(vd_matrix[dim_name].values).to_pydatetime()


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
