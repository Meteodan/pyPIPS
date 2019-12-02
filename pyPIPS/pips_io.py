"""pyPIPS.io: a collection of functions to read and write disdrometer and other surface probe
data.
"""
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
import numpy as np
from . import thermolib as thermo
from . import utils
from . import parsivel_params


PIPS_timestamp_format = '%Y-%m-%d %H:%M:%S'


def parse_PIPS_record(record, tripips=False):
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

    # Figure out which version we are reading in
    if not tripips:
        curfieldnames = parsivel_params.parsivel_parameters['PIPS_file_field_names']
    else:
        curfieldnames = parsivel_params.parsivel_parameters['TriPIPS_file_field_names']

    tokens = record.strip().split(',')
    timestamp = tokens[curfieldnames.index('TIMESTAMP')]
    token_dict['logger_datetime'] = datetime.strptime(timestamp, PIPS_timestamp_format)
    token_dict['voltage'] = np.float(tokens[curfieldnames.index('BattV')])
    token_dict['winddirrel'] = np.float(tokens[curfieldnames.index('WindDir')])
    token_dict['windspd'] = np.float(tokens[curfieldnames.index('WS_ms')])
    token_dict['winddiag'] = np.float(tokens[curfieldnames.index('WSDiag')])
    if not tripips:
        token_dict['fasttemp'] = np.float(tokens[curfieldnames.index('FastTemp')])
    token_dict['slowtemp'] = np.float(tokens[curfieldnames.index('SlowTemp')])
    if tripips:
        token_dict['fasttemp'] = token_dict['slowtemp']
    token_dict['RH'] = np.float(tokens[curfieldnames.index('RH')])
    token_dict['pressure'] = np.float(tokens[curfieldnames.index('Pressure')])
    token_dict['compass_dir'] = np.float(tokens[curfieldnames.index('FluxDirection')])
    token_dict['GPS_time'] = tokens[curfieldnames.index('GPSTime')]
    token_dict['GPS_status'] = tokens[curfieldnames.index('GPSStatus')]
    GPS_lat = np.float(tokens[curfieldnames.index('GPSLat')])
    GPS_lat_hem = tokens[curfieldnames.index('GPSLatHem')]
    token_dict['GPS_lat'] = utils.DDMtoDD(GPS_lat, GPS_lat_hem)
    GPS_lon = np.float(tokens[curfieldnames.index('GPSLon')])
    GPS_lon_hem = tokens[curfieldnames.index('GPSLonHem')]
    token_dict['GPS_lon'] = utils.DDMtoDD(GPS_lon, GPS_lon_hem)
    token_dict['GPS_spd'] = np.float(tokens[curfieldnames.index('GPSSpd')])
    token_dict['GPS_dir'] = np.float(tokens[curfieldnames.index('GPSDir')])
    token_dict['GPS_date'] = tokens[curfieldnames.index('GPSDate')]
    try:
        token_dict['GPS_magvar'] = np.float(tokens[curfieldnames.index('GPSMagVar')])
    except ValueError:
        token_dict['GPS_magvar'] = np.nan
    try:
        token_dict['GPS_alt'] = np.float(tokens[curfieldnames.index('GPSAlt')])
    except ValueError:
        token_dict['GPS_alt'] = np.nan
    try:
        winddirabs = np.float(tokens[curfieldnames.index('WindDirAbs')])
        if np.isnan(winddirabs):
            winddirabs = token_dict['winddirrel']
    except ValueError:
        winddirabs = np.nan
    token_dict['winddirabs'] = winddirabs
    try:
        dewpoint = np.float(tokens[curfieldnames.index('Dewpoint')])
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
        RH_derived = np.float(tokens[curfieldnames.index('RHDer')])
        if np.isnan(RH_derived):
            RH_derived = token_dict['RH']
    except ValueError:
        RH_derived = token_dict['RH']
    token_dict['RH_derived'] = RH_derived
    token_dict['parsivel_telegram'] = tokens[curfieldnames.index('ParsivelStr')]

    return token_dict


def get_PIPS_GPS_offset(filename, tripips=False):
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
    # Figure out which version we are reading in
    if not tripips:
        curfieldnames = parsivel_params.parsivel_parameters['PIPS_file_field_names']
    else:
        curfieldnames = parsivel_params.parsivel_parameters['TriPIPS_file_field_names']

    firstgoodGPS = False
    with open(filename, 'r') as disfile:
        for line in disfile:
            if firstgoodGPS:
                break
            tokens = line.strip().split(',')
            timestamp = tokens[curfieldnames.index('TIMESTAMP')]
            logger_datetime = datetime.strptime(timestamp, PIPS_timestamp_format)
            GPS_time = tokens[curfieldnames.index('GPSTime')]
            GPS_status = tokens[curfieldnames.index('GPSStatus')]
            GPS_date = tokens[curfieldnames.index('GPSDate')]
            try:
                GPS_alt = np.float(tokens[curfieldnames.index('GPSAlt')])
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


def read_PIPS(filename, starttimestamp=None, stoptimestamp=None, tripips=False,
              correct_logger_time=True):
    """[summary]

    Parameters
    ----------
    filename : [type]
        [description]
    starttimestamp : [type], optional
        [description], by default None
    stoptimestamp : [type], optional
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
    if starttimestamp is not None:
        starttime = datetime.strptime(starttimestamp, '%Y%m%d%H%M%S')
    else:
        starttime = None
    if stoptimestamp is not None:
        stoptime = datetime.strptime(stoptimestamp, '%Y%m%d%H%M%S')
    else:
        stoptime = None

    conv_dict_list = []
    parsivel_dict_list = []
    vd_matrix_list = []
    # Open the file and start parsing the records
    with open(filename, 'r') as disfile:
        for record in disfile:
            record_dict = parse_PIPS_record(record, tripips=tripips)
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
        GPS_offset = get_PIPS_GPS_offset(filename, tripips=tripips)
        conv_df['logger_datetime'] += GPS_offset
        parsivel_df['parsivel_datetime'] += GPS_offset

    # Have to do this because of some weird issue where the DataArray constructor below errors out
    # if you try to give the 'time_10s' coordinate the actual parsivel_df.index as its values.
    # Instead, we save the pd.Series of datetimes ina  separate variable here *before* setting it
    # as the index to parsivel_df, and then use it as the values for the DataArray
    # time_10s coordinates. Bizarre.
    parsivel_datetime = parsivel_df['parsivel_datetime']

    # Set the pd.Series of datetimes as the new index for both the conv_df and parsivel_df
    # DataFrames
    conv_df = conv_df.set_index('logger_datetime')
    parsivel_df = parsivel_df.set_index('parsivel_datetime')

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
    # coordinate 'time_10s'.
    vd_matrix_da = \
        xr.DataArray(vd_matrix_arr,
                     coords={'time_10s': parsivel_datetime,
                             'fallspeed': ('fallspeed_bin', fallspeeds),
                             'diameter': ('diameter_bin', diameters),
                             'min_diameter': ('diameter_bin', min_diameters),
                             'max_diameter': ('diameter_bin', max_diameters),
                             'min_fallspeeds': ('fallspeed_bin', min_fallspeeds),
                             'max_fallspeeds': ('fallspeed_bin', max_fallspeeds)
                             },
                     dims=['time_10s', 'fallspeed_bin', 'diameter_bin'])

    return conv_df, parsivel_df, vd_matrix_da


def correct_PIPS(serialnum, infile, outfile):
    """Corrects offset Parsivel strings in a PIPS data file"""

    with open(infile, 'r') as disfile:

        lines_out = []
        parsivel_string_out_list = []
        parsivel_string_out = ''
        truncated = False
        first = True

        for line in disfile:
            line = line.rstrip('\r\n')
            tokens = line.strip().split(',')
            header_info = ",".join(tokens[:26])
            parsivel_string = tokens[26]
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

        # Sort the output lines by record #
        sorted_lines_out = sorted(lines_out, key=lambda record: int(record.strip().split(',')[1]))
        # sorted_lines_out = lines_out

    with open(outfile, 'w') as outdisfile:
        for line in sorted_lines_out:
            outdisfile.write(line + '\n')