# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Calculates various DSD fits and parameters from PIPS data"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')

args = parser.parse_args()

# Dynamically import the case configuration file
utils.log("Case config file is {}".format(args.case_config_path))
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Dynamically import the plotting configuration file
utils.log("Plotting configuration file is {}".format(args.plot_config_path))
try:
    pc = utils.import_all_from(args.plot_config_path)
    utils.log("Successfully imported pyPIPS control parameters!")
except Exception:
    utils.warning(
        "Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
    import configs.plot_config_default as pc

# Extract needed lists and variables from PIPS_IO_dict configuration dictionary
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Create the directory for the meteogram plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Read in the PIPS data for the deployment
conv_df_list = []
conv_resampled_df_list = []
parsivel_df_list = []
vd_matrix_da_list = []

for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype in zip(range(
        0, len(PIPS_filenames)), PIPS_filenames, PIPS_names, start_times, end_times, geo_locs,
        PIPS_types):

    tripips = (ptype == 'TriPIPS')
    PIPS_data_file_path = os.path.join(PIPS_dir, PIPS_filename)
    conv_df, parsivel_df, vd_matrix_da = pipsio.read_PIPS(PIPS_data_file_path,
                                                          start_timestamp=start_time,
                                                          end_timestamp=end_time, tripips=tripips)

    # We need the disdrometer locations. If they aren't supplied in the input control file, find
    # them from the GPS data
    if not geo_loc:
        geo_locs[index] = pipsio.get_PIPS_loc(conv_df['GPS_status'], conv_df['GPS_lat'],
                                              conv_df['GPS_lon'], conv_df['GPS_alt'])
    print("Lat/Lon/alt of {}: {}".format(PIPS_name, str(geo_loc)))

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Do some QC on the V-D matrix
    strongwindQC = pc.PIPS_qc_dict['strongwindQC']
    splashingQC = pc.PIPS_qc_dict['splashingQC']
    marginQC = pc.PIPS_qc_dict['marginQC']
    rainfallQC = pc.PIPS_qc_dict['rainfallQC']
    rainonlyQC = pc.PIPS_qc_dict['rainonlyQC']
    hailonlyQC = pc.PIPS_qc_dict['hailonlyQC']
    graupelonlyQC = pc.PIPS_qc_dict['graupelonlyQC']
    basicQC = pc.PIPS_qc_dict['basicQC']

    if basicQC:
        strongwindQC = True
        splashingQC = True
        marginQC = True

    if strongwindQC:
        vd_matrix_da = pqc.strongwindQC(vd_matrix_da)
    if splashingQC:
        vd_matrix_da = pqc.splashingQC(vd_matrix_da)
    if marginQC:
        vd_matrix_da = pqc.marginQC(vd_matrix_da)
    if rainfallQC:
        fallspeedmask = pqc.get_fallspeed_mask(avg_diameter, avg_fall_bins)
        vd_matrix_da = pqc.rainfallspeedQC(vd_matrix_da, fallspeedmask)
    if rainonlyQC:
        vd_matrix_da = pqc.rainonlyQC(vd_matrix_da)
    if hailonlyQC:
        vd_matrix_da = pqc.hailonlyQC(vd_matrix_da)
    # if graupelonlyQC:
    #     vd_matrix_da = pqc.graupelonlyQC(vd_matrix_da)

    # Compute the number density ND from the V-D matrix
    # Use measured fallspeed by default

    # Resample the parsivel data to a longer interval if desired
    if requested_interval > 10.:
        DSD_interval = pips.check_requested_resampling_interval(requested_interval, 10.)
        vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
        parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
    else:
        DSD_interval = 10.

    PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
    # Resample conventional data to the parsivel times
    sec_offset = PSD_datetimes[0].second
    conv_resampled_df = pips.resample_conv(ptype, DSD_interval, sec_offset, conv_df)
    conv_resampled_df_index = conv_resampled_df.index.intersection(parsivel_df.index)
    conv_resampled_df = conv_resampled_df.loc[conv_resampled_df_index]

    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_resampled_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    logND = np.log10(ND)

    # Compute various fits using the MM and TMM

    M2, _ = dsd.calc_moment_bin(ND, moment=2)
    M3, _ = dsd.calc_moment_bin(ND, moment=3)
    M4, _ = dsd.calc_moment_bin(ND, moment=4)
    M6, _ = dsd.calc_moment_bin(ND, moment=6)

    DSD_MM24 = dsd.fit_DSD_MM24(M2, M4)
    DSD_MM36 = dsd.fit_DSD_MM36(M3, M6)
    DSD_MM346 = dsd.fit_DSD_MM346(M3, M4, M6)
    DSD_MM246 = dsd.fit_DSD_MM246(M2, M4, M6)
    DSD_MM234 = dsd.fit_DSD_MM234(M2, M3, M4)

    D_min, D_max = dsd.get_max_min_diameters(ND)
    DSD_TMM246 = dsd.fit_DSD_TMM_xr(M2, M4, M6, D_min, D_max)

    # Wrap fits into a DataSet and dump to netCDF file
    data_arrays = []
    for da_tuple in [DSD_MM24, DSD_MM36, DSD_MM346, DSD_MM246, DSD_MM234, DSD_TMM246]:
        da_concat = xr.concat(da_tuple, pd.Index(['N0', 'lamda', 'alpha'], name='parameter'))
        data_arrays.append(da_concat)
    names = ['DSD_MM24', 'DSD_MM36', 'DSD_MM346', 'DSD_MM246', 'DSD_MM234', 'DSD_TMM246']
    fits_ds = xr.Dataset({name: da for name, da in zip(names, data_arrays)})

    # Add ND to the Dataset
    fits_ds = pipsio.combine_parsivel_data(fits_ds, ND, name='ND')

    # Compute sigma and Dm and add to Dataset
    D = ND['diameter']
    dD = ND['max_diameter'] - ND['min_diameter']

    Dm = dsd.calc_Dmpq_binned(4, 3, ND)
    sigma = dsd.calc_sigma(D, dD, ND)

    fits_ds = pipsio.combine_parsivel_data(fits_ds, Dm, name='Dm43')
    fits_ds = pipsio.combine_parsivel_data(fits_ds, sigma, name='sigma')

    # Add some metadata
    fits_ds.attrs['DSD_interval'] = DSD_interval
    fits_ds.attrs['strongwindQC'] = int(strongwindQC)
    fits_ds.attrs['splashingQC'] = int(splashingQC)
    fits_ds.attrs['marginQC'] = int(marginQC)
    fits_ds.attrs['rainfallQC'] = int(rainfallQC)
    fits_ds.attrs['rainonlyQC'] = int(rainonlyQC)
    fits_ds.attrs['hailonlyQC'] = int(hailonlyQC)
    fits_ds.attrs['graupelonlyQC'] = int(graupelonlyQC)
    fits_ds.attrs['basicQC'] = int(basicQC)

    ncfile_name = 'DSD_fits_params_{}_{:d}s_{}_{}.nc'.format(PIPS_name, int(DSD_interval), start_time,
                                                             end_time)
    ncfile_path = os.path.join(PIPS_dir, ncfile_name)
    print("Dumping {}".format(ncfile_path))
    fits_ds.to_netcdf(ncfile_path)
