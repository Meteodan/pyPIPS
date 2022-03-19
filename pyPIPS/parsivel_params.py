"""pyPIPS.parsivel_params: contains several fixed parameters related to the Parsivel disdrometers
"""
import numpy as np
import xarray as xr
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

parsivel_parameters['diameter_bin_edges_mm'] = np.append(
    parsivel_parameters['min_diameter_bins_mm'],
    parsivel_parameters['max_diameter_bins_mm'][-1])
parsivel_parameters['fallspeed_bin_edges_mps'] = np.append(
    parsivel_parameters['min_fallspeed_bins_mps'],
    parsivel_parameters['max_fallspeed_bins_mps'][-1])

parsivel_parameters['avg_diameter_bins_mm'] = (parsivel_parameters['min_diameter_bins_mm'] +
                                               parsivel_parameters['max_diameter_bins_mm']) / 2.

parsivel_parameters['avg_fallspeed_bins_mps'] = (parsivel_parameters['min_fallspeed_bins_mps'] +
                                                 parsivel_parameters['max_fallspeed_bins_mps']) / 2.

# Jaffrain and Berne (2011), Tokay et al. (2011). Units are mm^2
parsivel_parameters['eff_sensor_area_mm2'] = \
    (180. * (30. - parsivel_parameters['avg_diameter_bins_mm'] / 2.))

probe_info = {
    'PIPS1A': {
        'serialnum': '304545',
        'parsivel_angle': -45.,
    },
    'PIPS1B': {
        'serialnum': '295153',
        'parsivel_angle': 45.,
    },
    'PIPS2A': {
        'serialnum': '295166',
        'parsivel_angle': -45.,
    },
    'PIPS2B': {
        'serialnum': '304543',
        'parsivel_angle': 45.,
    },
    'TriPIPS': {
        'serialnum': '390654',
        'parsivel_angle': 90.,
    }
}

RB15_RR_min = np.array([0., 0.1, 0.25, 0.5, 1., 2., 200.])
RB15_RR_max = np.array([0.1, 0.25, 0.5, 1., 2., 200., np.inf])

RB15_correction_array = np.array([[1.00, 1.00, 0.02, 0.03, 0.11, 0.20, 0.36, 0.55, 0.86, 0.74, 1.04,
                                   1.10, 1.14, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00, 1.00, 0.04, 0.05, 0.16, 0.26, 0.47, 0.55, 0.85, 0.84, 1.13,
                                   1.20, 0.97, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00, 1.00, 0.04, 0.05, 0.19, 0.29, 0.52, 0.67, 0.94, 1.08, 1.22,
                                   1.35, 1.34, 1.25, 1.29, 1.43, 0.51, 1.00, 1.00, 1.00, 1.00, 1.00,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00, 1.00, 0.05, 0.07, 0.22, 0.36, 0.53, 0.67, 0.89, 0.90, 1.12,
                                   1.19, 1.17, 1.17, 1.17, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00, 1.00, 0.06, 0.11, 0.30, 0.45, 0.71, 0.80, 1.01, 1.17, 1.36,
                                   1.37, 1.41, 1.22, 1.43, 1.37, 1.31, 1.00, 1.00, 1.00, 1.00, 1.00,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00, 1.00, 0.07, 0.16, 0.36, 0.54, 0.78, 0.86, 1.03, 1.03, 1.12,
                                   1.10, 1.04, 0.97, 1.06, 1.07, 1.02, 0.97, 0.73, 0.58, 0.45, 0.32,
                                   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00],
                                  [1.00] * 32])
avg_diameter = parsivel_parameters['avg_diameter_bins_mm']
RB15_correction_factors = xr.DataArray(RB15_correction_array, name='RB15_correction_factors',
                                       coords={'rainrate': ('rainrate', RB15_RR_min),
                                               'RR_min': ('rainrate', RB15_RR_min),
                                               'RR_max': ('rainrate', RB15_RR_max),
                                               'diameter_bin': ('diameter_bin', avg_diameter),
                                               },
                                       dims=['rainrate', 'diameter_bin'])
