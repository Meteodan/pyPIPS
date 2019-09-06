"""pyPIPS.parsivel_params: contains several fixed parameters related to the Parsivel disdrometers
"""
import numpy as np
# Dictionaries containing some metadata for the various probes

# length units: mm
# velocity units: m/s
parsivel_parameters = {
    'PIPS_file_field_names': ['TIMESTAMP', 'RECORD', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms',
                              'WSDiag', 'FastTemp', 'SlowTemp', 'RH', 'Pressure', 'FluxDirection',
                              'GPSTime', 'GPSStatus', 'GPSLat', 'GPSLatHem', 'GPSLon', 'GPSLonHem',
                              'GPSSpd', 'GPSDir', 'GPSDate', 'GPSMagVar', 'GPSAlt', 'WindDirAbs',
                              'Dewpoint', 'RHDer', 'ParsivelStr'],
    'TriPIPS_file_field_names': ['TIMESTAMP', 'BattV', 'PTemp_C', 'WindDir', 'WS_ms', 'WSDiag',
                                 'SlowTemp', 'RH', 'Pressure', 'FluxDirection', 'GPSTime',
                                 'GPSStatus', 'GPSLat', 'GPSLatHem', 'GPSLon', 'GPSLonHem',
                                 'GPSSpd', 'GPSDir', 'GPSDate', 'GPSMagVar', 'GPSAlt', 'WindDirAbs',
                                 'Dewpoint', 'RHDer', 'ParsivelStr'],
    'min_diameter_bins_mm': np.array([0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000,
                                      1.125, 1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 3.000, 3.500,
                                      4.000, 4.500, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000,
                                      12.000, 14.000, 16.000, 18.000, 20.000, 23.000]),
    'max_diameter_bins_mm': np.array([0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000, 1.125,
                                      1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 3.000, 3.500, 4.000,
                                      4.500, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000, 12.000,
                                      14.000, 16.000, 18.000, 20.000, 23.000, 26.000]),
    'min_fallspeed_bins_mps': np.array([0.000, 0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700,
                                        0.800, 0.900, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000,
                                        2.400, 2.800, 3.200, 3.600, 4.000, 4.800, 5.600, 6.400,
                                        7.200, 8.000, 9.600, 11.200, 12.800, 14.400, 16.000,
                                        19.200]),
    'max_fallspeed_bins_mps': np.array([0.100, 0.200, 0.300, 0.400, 0.500, 0.600, 0.700, 0.800,
                                        0.900, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.400,
                                        2.800, 3.200, 3.600, 4.000, 4.800, 5.600, 6.400, 7.200,
                                        8.000, 9.600, 11.200, 12.800, 14.400, 16.000, 19.200,
                                        22.400]),
    'sensor_length_mm': 180.,
    'sensor_width_mm': 30.,
    'sensor_area_mm2': 5400.
}

parsivel_parameters['avg_diameter_bins_mm'] = (parsivel_parameters['min_diameter_bins_mm'] +
                                               parsivel_parameters['max_diameter_bins_mm']) / 2.

parsivel_parameters['avg_fallspeed_bins_mps'] = (parsivel_parameters['min_fallspeed_bins_mps'] +
                                                 parsivel_parameters['max_fallspeed_bins_mps']) / 2.

# Jaffrain and Berne (2011), Tokay et al. (2011). Units are mm^2
parsivel_parameters['eff_sensor_area_mm2'] = \
    (180. * (30. - parsivel_parameters['avg_diameter_bins_mm'] / 2.))

probe_info = {
    'PIPS1A': {
        'serialnum': '304545'
    },
    'PIPS1B': {
        'serialnum': '295153'
    },
    'PIPS2A': {
        'serialnum': '295166'
    },
    'PIPS2B': {
        'serialnum': '304543'
    },
    'TriPIPS': {
        'serialnum': '390654'
    }
}
