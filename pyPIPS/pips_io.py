"""pips_io.py: a collection of functions to read and write disdrometer and other surface probe
data.
"""
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
import numpy as np
from . import thermolib as thermo
from . import utils
from . import parsivel_params


PIPS_timestamp_format1 = '%Y-%m-%d %H:%M:%S'
PIPS_timestamp_format2 = '%Y-%m-%d %H:%M:%S.000'


def get_field_indices(header, tripips=False):

    header_fields = header.strip().split(',')
    if 'TIMESTAMP' not in header_fields:
        if not tripips:
            header_fields = parsivel_params.parsivel_parameters['PIPS_file_field_names']
        else:
            header_fields = parsivel_params.parsivel_parameters['TriPIPS_file_field_names']
    return {field: header_fields.index(field) for field in header_fields}


def parse_PIPS_record(record, field_indices, tripips=False):
    """[summary]

    Parameters
    ----------
    record : [type]
        [description]
    tripips : bool, optional
        [description], by default False

    Returns
    -------
    [type]
        [description]
    """
    # TODO: refactor this function to do more efficient parsing?
    token_dict = {}

    tokens = record.strip().split(',')
    timestamp = tokens[field_indices['TIMESTAMP']]
    if timestamp[-4:] == '.000':
        token_dict['logger_datetime'] = datetime.strptime(timestamp, PIPS_timestamp_format2)
    else:
        token_dict['logger_datetime'] = datetime.strptime(timestamp, PIPS_timestamp_format1)
    token_dict['voltage'] = np.float(tokens[field_indices['BattV']])
    token_dict['winddirrel'] = np.float(tokens[field_indices['WindDir']])
    token_dict['windspd'] = np.float(tokens[field_indices['WS_ms']])
    token_dict['winddiag'] = np.float(tokens[field_indices['WSDiag']])
    if not tripips:
        token_dict['fasttemp'] = np.float(tokens[field_indices['FastTemp']])
    token_dict['slowtemp'] = np.float(tokens[field_indices['SlowTemp']])
    if tripips:
        token_dict['fasttemp'] = token_dict['slowtemp']
    token_dict['RH'] = np.float(tokens[field_indices['RH']])
    token_dict['pressure'] = np.float(tokens[field_indices['Pressure']])
    token_dict['compass_dir'] = np.float(tokens[field_indices['FluxDirection']])
    token_dict['GPS_time'] = tokens[field_indices['GPSTime']]
    token_dict['GPS_status'] = tokens[field_indices['GPSStatus']]
    GPS_lat = np.float(tokens[field_indices['GPSLat']])
    GPS_lat_hem = tokens[field_indices['GPSLatHem']]
    token_dict['GPS_lat'] = utils.DDMtoDD(GPS_lat, GPS_lat_hem)
    GPS_lon = np.float(tokens[field_indices['GPSLon']])
    GPS_lon_hem = tokens[field_indices['GPSLonHem']]
    token_dict['GPS_lon'] = utils.DDMtoDD(GPS_lon, GPS_lon_hem)
    token_dict['GPS_spd'] = np.float(tokens[field_indices['GPSSpd']])
    token_dict['GPS_dir'] = np.float(tokens[field_indices['GPSDir']])
    token_dict['GPS_date'] = tokens[field_indices['GPSDate']]
    try:
        token_dict['GPS_magvar'] = np.float(tokens[field_indices['GPSMagVar']])
    except ValueError:
        token_dict['GPS_magvar'] = np.nan
    try:
        token_dict['GPS_alt'] = np.float(tokens[field_indices['GPSAlt']])
    except ValueError:
        token_dict['GPS_alt'] = np.nan
    try:
        winddirabs = np.float(tokens[field_indices['WindDirAbs']])
        if np.isnan(winddirabs):
            winddirabs = token_dict['winddirrel']
    except ValueError:
        winddirabs = np.nan
    token_dict['winddirabs'] = winddirabs
    try:
        dewpoint = np.float(tokens[field_indices['Dewpoint']])
        if np.isnan(dewpoint):
            dewpoint = (thermo.calTdfromRH(token_dict['pressure'] * 100.,
                                           token_dict['fasttemp'] + 273.15,
                                           token_dict['RH'] / 100.) - 273.15)
    except ValueError:
        dewpoint = thermo.calTdfromRH(token_dict['pressure'] * 100.,
                                      token_dict['fasttemp'] + 273.15,
                                      token_dict['RH'] / 100.) - 273.15
    token_dict['dewpoint'] = dewpoint
    try:
        RH_derived = np.float(tokens[field_indices['RHDer']])
        if np.isnan(RH_derived):
            RH_derived = token_dict['RH']
    except ValueError:
        RH_derived = token_dict['RH']
    token_dict['RH_derived'] = RH_derived
    token_dict['parsivel_telegram'] = tokens[field_indices['ParsivelStr']]

    return token_dict


def get_PIPS_GPS_offset(filename, field_indices, tripips=False):
    """[summary]

    Parameters
    ----------
    filename : [type]
        [description]
    tripips : bool, optional
        [description], by default False

    Returns
    -------
    [type]
        [description]
    """

    firstgoodGPS = False
    with open(filename, 'r') as disfile:
        for line in disfile:
            if firstgoodGPS:
                break
            tokens = line.strip().split(',')
            timestamp = tokens[field_indices['TIMESTAMP']]
            try:
                if '.000' in timestamp:
                    logger_datetime = datetime.strptime(timestamp, PIPS_timestamp_format2)
                else:
                    logger_datetime = datetime.strptime(timestamp, PIPS_timestamp_format1)
            except:
                continue
            GPS_time = tokens[field_indices['GPSTime']]
            GPS_status = tokens[field_indices['GPSStatus']]
            GPS_date = tokens[field_indices['GPSDate']]
            try:
                GPS_alt = np.float(tokens[field_indices['GPSAlt']])
            except ValueError:
                GPS_alt = np.nan

            # print("GPS_status: {}, GPS_alt: {:g}".format(GPS_status, GPS_alt))

            if not np.isnan(GPS_alt) and not firstgoodGPS and GPS_status == 'A':
                firstgoodGPS = True

                # Construct datetime object
                gyear = np.int('20' + GPS_date[4:])
                gmonth = np.int(GPS_date[2:4])
                gday = np.int(GPS_date[:2])
                ghour = np.int(GPS_time[:2])
                gmin = np.int(GPS_time[2:4])
                gsec = np.int(GPS_time[4:])

                GPS_datetime = datetime(gyear, gmonth, gday, ghour, gmin, gsec)
                GPS_offset = GPS_datetime - logger_datetime
                print('GPS time: {}, Logger time: {}'.format(GPS_datetime.ctime(),
                                                             logger_datetime.ctime()))
                print('GPS Offset: {}'.format(str(GPS_offset)))

    if not firstgoodGPS:
        print('No GPS time information in the file! Setting offset to zero!')
        GPS_offset = timedelta(seconds=0)

    return GPS_offset


def parse_parsivel_telegram(parsivel_telegram, logger_datetime):
    """[summary]

    Parameters
    ----------
    parsivel_telegram : [type]
        [description]
    logger_datetime : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """
    parsivel_tokens = parsivel_telegram.strip().split(';')
    serialnum = parsivel_tokens[0]
    valid_serial_nums = [v['serialnum'] for k, v in parsivel_params.probe_info.items()]
    if serialnum in valid_serial_nums and len(parsivel_tokens) >= 11:
        precipintensity = np.float(parsivel_tokens[1])
        precipaccum = np.float(parsivel_tokens[2])
        parsivel_dBZ = np.float(parsivel_tokens[3])
        sample_interval = np.float(parsivel_tokens[4])
        signal_amplitude = np.float(parsivel_tokens[5])
        pcount = np.int(parsivel_tokens[6])
        sensor_temp = np.float(parsivel_tokens[7])
        pvoltage = np.float(parsivel_tokens[8])
        try:
            vd_matrix = [float(x) if x != '' else 0 for x in parsivel_tokens[11:]]
        except Exception:  # TODO: find out what exception to catch here. ValueError?
            vd_matrix = np.array([np.nan for i in range(1025)])
        if len(vd_matrix) < 1024:
            print("Problem with Parsivel vd_matrix.  Flagging as bad!")
            print("Time: {}".format(logger_datetime.ctime()))
            vd_matrix = np.array([np.nan for i in range(1025)])
        else:
            if len(vd_matrix) == 1025:
                vd_matrix = vd_matrix[:-1]  # Strip off bogus last value
        # Now create an array out of the vd_matrix and reshape it to 32x32
        vd_matrix = np.array(vd_matrix, dtype='int')
        if vd_matrix.size == 1024:
            vd_matrix = vd_matrix.reshape((32, 32))
        else:
            vd_matrix = np.empty((32, 32), dtype='int')
            vd_matrix[:] = np.nan
        parsivel_dict = {
            'parsivel_datetime': logger_datetime, 'precipintensity': precipintensity,
            'precipaccum': precipaccum, 'parsivel_dBZ': parsivel_dBZ,
            'sample_interval': sample_interval, 'signal_amplitude': signal_amplitude,
            'pcount': pcount, 'sensor_temp': sensor_temp, 'pvoltage': pvoltage
        }
    else:
        parsivel_dict = None
        vd_matrix = None

    return parsivel_dict, vd_matrix


def read_PIPS(filename, start_timestamp=None, end_timestamp=None, tripips=False,
              correct_logger_time=True, sort=False, check_order=False):
    """[summary]

    Parameters
    ----------
    filename : [type]
        [description]
    start_timestamp : [type], optional
        [description], by default None
    end_timestamp : [type], optional
        [description], by default None
    tripips : bool, optional
        [description], by default False
    correct_logger_time : bool, optional
        [description], by default True

    Returns
    -------
    [type]
        [description]
    """
    if start_timestamp is not None:
        starttime = datetime.strptime(start_timestamp, '%Y%m%d%H%M%S')
    else:
        starttime = None
    if end_timestamp is not None:
        stoptime = datetime.strptime(end_timestamp, '%Y%m%d%H%M%S')
    else:
        stoptime = None

    conv_dict_list = []
    parsivel_dict_list = []
    vd_matrix_list = []
    # Open the file and attempt to read the header
    with open(filename, 'r') as disfile:
        header = disfile.readline()
        field_indices = get_field_indices(header, tripips=tripips)

    # Open the file again and start parsing the records
    with open(filename, 'r') as disfile:
        for record in disfile:
            try:
                record_dict = parse_PIPS_record(record, field_indices, tripips=tripips)
            except:
                continue
            # Skip this record if it lies before or after the desired period
            if starttime is not None and record_dict['logger_datetime'] < starttime:
                continue
            if stoptime is not None and record_dict['logger_datetime'] > stoptime:
                continue

            parsivel_dict, vd_matrix = parse_parsivel_telegram(record_dict['parsivel_telegram'],
                                                               record_dict['logger_datetime'])
            conv_dict = record_dict
            conv_dict.pop('parsivel_telegram')
            conv_dict_list.append(conv_dict)
            if parsivel_dict is not None:
                parsivel_dict_list.append(parsivel_dict)
                vd_matrix_list.append(vd_matrix)

    # Create pandas dataframes out of these lists
    conv_df = pd.DataFrame(conv_dict_list)
    parsivel_df = pd.DataFrame(parsivel_dict_list)

    if correct_logger_time:
        GPS_offset = get_PIPS_GPS_offset(filename, field_indices, tripips=tripips)
        conv_df['logger_datetime'] += GPS_offset
        parsivel_df['parsivel_datetime'] += GPS_offset

    # Have to do this because of some weird issue where the DataArray constructor below errors out
    # if you try to give the 'time' coordinate the actual parsivel_df.index as its values.
    # Instead, we save the pd.Series of datetimes ina  separate variable here *before* setting it
    # as the index to parsivel_df, and then use it as the values for the DataArray
    # time coordinates. Bizarre.
    parsivel_datetime = parsivel_df['parsivel_datetime']

    # Set the pd.Series of datetimes as the new index for both the conv_df and parsivel_df
    # DataFrames
    conv_df = conv_df.rename(columns={'logger_datetime': 'time'})
    conv_df = conv_df.set_index('time')
    parsivel_df = parsivel_df.rename(columns={'parsivel_datetime': 'time'})
    parsivel_df = parsivel_df.set_index('time')


    # Create an xarray DataArray for the vd_matrix
    vd_matrix_arr = np.dstack(vd_matrix_list)
    vd_matrix_arr = np.rollaxis(vd_matrix_arr, 2, 0)

    diameters = parsivel_params.parsivel_parameters['avg_diameter_bins_mm']
    min_diameters = parsivel_params.parsivel_parameters['min_diameter_bins_mm']
    max_diameters = parsivel_params.parsivel_parameters['max_diameter_bins_mm']
    fallspeeds = parsivel_params.parsivel_parameters['avg_fallspeed_bins_mps']
    min_fallspeeds = parsivel_params.parsivel_parameters['min_fallspeed_bins_mps']
    max_fallspeeds = parsivel_params.parsivel_parameters['max_fallspeed_bins_mps']
    # Note, on below, I get an error if I try to use the parsivel_df.index as the value of the
    # coordinate 'time'.
    vd_matrix_da = \
        xr.DataArray(vd_matrix_arr,
                     name='VD_matrix',
                     coords={'time': parsivel_datetime,
                             'fallspeed': ('fallspeed_bin', fallspeeds),
                             'diameter': ('diameter_bin', diameters),
                             'min_diameter': ('diameter_bin', min_diameters),
                             'max_diameter': ('diameter_bin', max_diameters),
                             'min_fallspeeds': ('fallspeed_bin', min_fallspeeds),
                             'max_fallspeeds': ('fallspeed_bin', max_fallspeeds)
                             },
                     dims=['time', 'fallspeed_bin', 'diameter_bin'])

    if check_order:
        time_diff = conv_df.to_xarray()['time'].diff('time').astype(np.float)*1.e-9
        print(time_diff)
        out_of_order_times = time_diff.where(time_diff < 0, drop=True)['time']
        print(out_of_order_times)
        if out_of_order_times.sizes['time'] > 0:
            print("Data contains times that are out of order! You might want to sort them...")
            print(out_of_order_times)

    if sort:
        conv_df = conv_df.sort_index()
        parsivel_df = parsivel_df.sort_index()
        vd_matrix_da = vd_matrix_da.sortby('time')

    return conv_df, parsivel_df, vd_matrix_da


def get_PIPS_loc(GPS_stats, GPS_lats, GPS_lons, GPS_alts):
    """[summary]

    Parameters
    ----------
    GPS_stats : [type]
        [description]
    GPS_lats : [type]
        [description]
    GPS_lons : [type]
        [description]
    GPS_alts : [type]
        [description]
    """

    # Find the disdrometer location by averaging the valid GPS lats,lons, and alts
    GPSnonvalid = (np.array(GPS_stats) != 'A')
    GPS_lats_masked = np.ma.masked_where(GPSnonvalid, np.array(GPS_lats))
    GPS_lons_masked = np.ma.masked_where(GPSnonvalid, np.array(GPS_lons))
    GPS_alts_masked = np.ma.masked_where(GPSnonvalid, np.array(GPS_alts))

    GPS_lats_masked = np.ma.masked_invalid(GPS_lats_masked)
    GPS_lons_masked = np.ma.masked_invalid(GPS_lons_masked)
    GPS_alts_masked = np.ma.masked_invalid(GPS_alts_masked)

    lat = GPS_lats_masked.mean()
    lon = GPS_lons_masked.mean()
    alt = GPS_alts_masked.mean()

    # There's a bug in numpy.ma that sometimes causes operations such as mean() to return
    # a 0d array instead of a scalar.  To remedy this, explicitly cast them as scalars
    # here.

    lat = lat.item()
    lon = lon.item()
    alt = alt.item()

    return lat, lon, alt


def correct_PIPS(serialnum, infile, outfile):
    """Corrects offset Parsivel strings in a PIPS data file"""

    with open(infile, 'r') as disfile:

        lines_out = []
        parsivel_string_out_list = []
        parsivel_string_out = ''
        truncated = False
        first = True
        # Get the header line
        header_line = next(disfile)
        for line in disfile:
            line = line.rstrip('\r\n')
            tokens = line.strip().split(',')
            header_info = ",".join(tokens[:-1])
            parsivel_string = tokens[-1]
            parsivel_tokens = parsivel_string.strip().split(';')
            if len(parsivel_tokens) > 1:
                if truncated:  # Previous record was truncated
                    # Look for the parsivel serial number somewhere in the next record
                    # and find the index if it exists
                    try:
                        line_out = header_info
                        sindex = [i for i, x in enumerate(parsivel_tokens) if serialnum in x][0]
                        # Concatenate portion of string before serial number onto end of previous
                        # record
                        parsivel_string_out = (parsivel_string_out +
                                               ";".join(parsivel_tokens[:sindex]))
                        lines_out.append(line_out + ',' + parsivel_string_out)
                        parsivel_string_out_list.append(parsivel_string_out)
                        parsivel_string_out = ";".join(parsivel_tokens[sindex:]).lstrip()
                        truncated = False
                    except Exception:
                        print("Something went wrong!")
                elif not first:
                    # Should be able to just rip through the rest of the records and piece
                    # them together
                    sindex = [i for i, x in enumerate(parsivel_tokens) if serialnum in x][0]
                    parsivel_string_out = parsivel_string_out + ";".join(parsivel_tokens[:sindex])
                    parsivel_string_out_list.append(parsivel_string_out)
                    lines_out.append(line_out + ',' + parsivel_string_out)
                    line_out = header_info
                    parsivel_string_out = ";".join(parsivel_tokens[sindex:]).lstrip()
                elif first:
                    if len(parsivel_tokens) < 1036:    # Likely a truncated record
                        truncated = True
                        first = False
                        parsivel_string_out = parsivel_string
                        line_out = header_info
                    else:
                        lines_out.append(line)
            else:
                lines_out.append(line)

        # # Sort the output by time again
        # numrecords = len(lines_out)
        # datetime_onesec_list = []
        # for i in range(numrecords):  # Loop through records
        #     # Parse the timestamp for the records and create a datetime object out of it
        #     timestring = dict_onesec['TIMESTAMP'][i].strip().split()
        #     # print('i', 'timestring_onesec', i, timestring_onesec)
        #     datetime_onesec = parseTimeStamp(timestring_onesec)
        #     datetime_onesec_list.append(datetime_onesec)

        # Sort the records in increasing order of time. (Some files are corrupted and have some
        # of their records out of order)
        indices = list(range(numrecords))
        indices_sorted = sortby(indices, datetime_onesec_list)
        if not indices == indices_sorted:
            print("Times out of order! Need to sort!")
            datetime_onesec_list.sort()
            for field, records in dict_onesec.copy().items():
                dict_onesec[field] = [records[i] for i in indices_sorted]

        # Sort the output lines by record #
        # sorted_lines_out = sorted(lines_out, key=lambda record: int(record.strip().split(',')[1]))
        sorted_lines_out = lines_out

    with open(outfile, 'w') as outdisfile:
        outdisfile.write(header_line)
        for line in sorted_lines_out:
            outdisfile.write(line + '\n')


def conv_df_to_ds(conv_df):
    """Converts the conventional data from the PIPS from a pandas DataFrame to an xarray Dataset.
    Adds units as metadata

    Parameters
    ----------
    conv_df : pandas.DataFrame
        pandas DataFrame containing the 1-s "conventional" observations from the PIPS

    Returns
    -------
    xarray.Dataset
        Dataset version with additional metadata
    """
    conv_ds = conv_df.to_xarray()

    # These are not present when conv_df is resampled to longer intervals
    # TODO: update pips.resample_conv to properly handle these variables
    # NOTE: compass_dir should be working now
    try:
        conv_ds['compass_dir'].attrs['units'] = 'degrees'
    except KeyError:
        pass
    try:
        conv_ds['GPS_date'].attrs['units'] = 'DDMMYY'
        conv_ds['GPS_dir'].attrs['units'] = 'degrees'
        conv_ds['GPS_spd'].attrs['units'] = 'meters per second'
        conv_ds['GPS_time'].attrs['units'] = 'HHMMSS'
        conv_ds['winddirrel'].attrs['units'] = 'degrees'
    except KeyError:
        pass

    conv_ds['GPS_alt'].attrs['units'] = 'meters'
    conv_ds['GPS_lat'].attrs['units'] = 'degrees N'
    conv_ds['GPS_lon'].attrs['units'] = 'degrees E'
    conv_ds['GPS_lon'].attrs['description'] = 'west is negative'
    conv_ds['RH'].attrs['units'] = 'percent'
    conv_ds['RH_derived'].attrs['units'] = 'percent'
    conv_ds['dewpoint'].attrs['units'] = 'degrees Celsius'
    try:
        conv_ds['fasttemp'].attrs['units'] = 'degrees Celsius'
    except KeyError:  # TriPIPS doesn't yet have fasttemp
        pass
    conv_ds['pressure'].attrs['units'] = 'hectoPascals'
    conv_ds['slowtemp'].attrs['units'] = 'degrees Celsius'
    conv_ds['voltage'].attrs['units'] = 'volts'
    # TODO: add attributes to signify if the winds below are raw or resampled (i.e. averaged
    # over a period)
    conv_ds['winddirabs'].attrs['units'] = 'degrees'
    conv_ds['windspd'].attrs['units'] = 'meters per second'
    try:
        conv_ds['windgust'].attrs['units'] = 'meters per second'
        conv_ds['windgust'].attrs['description'] = \
            'Max 3-s wind over interval (given by DSD_interval)'
    except KeyError:
        pass

    return conv_ds


def parsivel_df_to_ds(parsivel_df):
    """Converts the parsivel derived data from the PIPS from a pandas DataFrame to an xarray
    Dataset. Adds units as metadata

    Parameters
    ----------
    parsivel_df : pandas.DataFrame
        pandas DataFrame containing the 10-s "derived" observations from the PIPS (i.e. everything
        in the telegram *but* the velocity-diameter matrix).

    Returns
    -------
    xarray.Dataset
        Dataset version with additional metadata
    """
    parsivel_ds = parsivel_df.to_xarray()
    parsivel_ds['parsivel_dBZ'].attrs['units'] = 'dBZ'
    parsivel_ds['pcount'].attrs['units'] = 'count'
    parsivel_ds['precipaccum'].attrs['units'] = 'millimeters'
    parsivel_ds['precipintensity'].attrs['units'] = 'millimeters per hour'
    parsivel_ds['pvoltage'].attrs['units'] = 'volts'
    parsivel_ds['sample_interval'].attrs['units'] = 'seconds'
    parsivel_ds['sensor_temp'].attrs['units'] = 'degrees Celsius'
    parsivel_ds['signal_amplitude'].attrs['units'] = 'unknown'
    parsivel_ds.attrs['nominal sample interval'] = '10 seconds'

    return parsivel_ds


def combine_parsivel_data(parsivel_ds, data_array, name=None, coord='time'):
    """Adds a new DataArray to an existing Parsivel Dataset.

    Parameters
    ----------
    parsivel_ds : xarray.Dataset
        xarray Dataset containing data from the Parsivel disdrometer
    data_array : xarray.DataArray
        xarray DataArray containing data to add to the Dataset. Must have the same times.
    name : str, optional
        name of the new DataArray, by default None. If None, will attempt to extract the name
        from the DataArray.name attribute.
    coord : str, optional
        coordinate name to check for consistency, by default 'time'

    Returns
    -------
    xarray.Dataset or None
        Dataset with the new DataArray added, or None if the operation fails.
    """

    try:
        if not np.array_equal(parsivel_ds[coord].values, data_array[coord].values):
            message = "Coordinate values {} for parsivel dataset and new data array do not match!"
            print(message.format(coord))
            return
        else:
            if name is None:
                if data_array.name is not None:
                    name = data_array.name
                else:
                    print("""Please provide a name for the new array with the "name" keyword
                          argument""")
                    return
            parsivel_ds[name] = data_array
            return parsivel_ds
    except KeyError:
        print("Problem with the time dimension! Aborting!")
        return


def reconstruct_MultiIndex(da, index_level_names, MultiIndex_name):
    """Reconstructs a MultiIndex from the given list of indices in a DataArray. Needed because
    serializing a MultiIndex to a netCDF file doesn't work, so if we need the MultiIndex (i.e. such
    as for performing grouping and averaging on it), after reading the DataArray back in from
    the netCDF file, we need to reconstruct it. This function does that.

    Parameters
    ----------
    da : xr.DataArray
        The DataArray
    index_level_names : list of str
        List of names of the index coordinates we want to recreate the MultiIndex from
    MultiIndex_name : str
        The name of the new MultiIndex

    Returns
    -------
    xr.DataArray
        The DataArray with the MultiIndex reconstructed
    """
    dim_name = da[index_level_names[0]].dims[0]  # NOTE: assumes only one dimension!
    index_values = []
    for index_level_name in index_level_names:
        index_values.append(da[index_level_name].values)
    da = da.reset_coords(names=index_level_names, drop=True)
    da.coords[MultiIndex_name] = (dim_name,
                                  pd.MultiIndex.from_arrays(index_values, names=index_level_names))
    return da


def remove_unneeded(ds):
    dim_names = [dim for dim in ds.dims if 'fields' in dim]
    ds = ds.drop_dims(dim_names, errors='ignore')
    # TODO: fix this!
    ds = ds.drop_vars(['flagged times'], errors='ignore')
    return ds


def reset_time_index(ds):
    return ds.reset_index('time')


def preprocess_dataset(ds):
    ds = remove_unneeded(ds)
    ds = reset_time_index(ds)
    return ds