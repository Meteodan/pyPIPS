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
description = "Performs SATP on DSDs from PIPS data"
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
dataset_name = config.PIPS_IO_dict.get('dataset_name', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Set up D0 and RR bins for the SATP procedure
D0_bins = np.arange(0.05, 6.05, 0.05)
# Following is from
# https://stackoverflow.com/questions/45234987/numpy-range-created-using-percentage-increment
# TODO: make a function out of this and put in pyPIPS.utils
RR_start = 0.1
RR_stop = 250.
RR_pct_incr = 10.
RR_incr = (100. + RR_pct_incr) / 100.
RR_bins = RR_start * np.full(int(np.log(RR_stop / RR_start) / np.log(RR_incr)),
                                 RR_incr).cumprod()

# Read in the PIPS data for the deployment
conv_df_list = []
conv_resampled_df_list = []
parsivel_df_list = []
vd_matrix_da_list = []
ND_list = []
D0_list = []
RR_list = []

for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype in zip(range(
        0, len(PIPS_filenames)), PIPS_filenames, PIPS_names, start_times, end_times, geo_locs,
        PIPS_types):

    tripips = (ptype == 'TriPIPS')
    PIPS_data_file_path = os.path.join(PIPS_dir, PIPS_filename)
    conv_df, parsivel_df, vd_matrix_da = pipsio.read_PIPS(PIPS_data_file_path,
                                                          start_timestamp=start_time,
                                                          end_timestamp=end_time, tripips=tripips)
    print("Read {}".format(PIPS_filename))
    # We need the disdrometer locations. If they aren't supplied in the input control file, find
    # them from the GPS data
    if not geo_loc:
        geo_locs[index] = pipsio.get_PIPS_loc(conv_df['GPS_status'], conv_df['GPS_lat'],
                                              conv_df['GPS_lon'], conv_df['GPS_alt'])
    print("Lat/Lon/alt of {}: {}".format(PIPS_name, str(geo_locs[index])))

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

    # Compute binned number concentration
    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_resampled_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    ND = ND.where(ND > 0.)
    # Drop all entries without DSDs
    ND = ND.dropna(dim='time', how='all')
    conv_resampled_df_index = conv_resampled_df.index.intersection(ND.indexes['time'])
    conv_resampled_df = conv_resampled_df.loc[conv_resampled_df_index]
    # Compute rainrate using empirical fallspeed curve
    # TODO: allow for the use of the measured fallspeeds in the rainrate calculation.
    fallspeeds_emp = pips.calc_empirical_fallspeed(avg_diameter, correct_rho=True,
                                                   rho=conv_resampled_df['rho'])
    rainrate_bin = (6. * 10.**-4.) * np.pi * fallspeeds_emp * avg_diameter**3. * ND * bin_width
    rainrate = rainrate_bin.sum(dim='diameter_bin')
    rainrate = rainrate.where(rainrate <= RR_bins[-1])
    rainrate = rainrate.where(rainrate > 0.)

    # Compute D0 (in mm)
    D0 = dsd.calc_D0_bin(ND) * 1000.
    D0 = D0.where(D0 <= D0_bins[-1])

    # Add D0 and RR coordinates to the ND DataArray
    ND.coords['D0'] = ('time', D0)
    ND.coords['RR'] = ('time', rainrate)
    # Digitize the D0 and RR using the bins computed earlier to get the indices
    # of the bins for each D0/RR pair for each DSD and make a new MultiIndex out of it
    D0_indices = np.digitize(D0, D0_bins)
    RR_indices = np.digitize(rainrate, RR_bins)
    ND.coords['D0_RR'] = ('time', pd.MultiIndex.from_arrays([D0_indices, RR_indices],
                                                            names=['D0_idx', 'RR_idx']))
    # Change the name of dimension 'time' to 'D0_RR' since we don't care about the timestamps
    # here. Also, this allows us to concatenate each
    # deployment's data into a single DataArray for later grouping and averaging by RR-D0 bin
    ND = ND.swap_dims({'time': 'D0_RR'})
    ND_list.append(ND)

# Ok, now combine the list of DSD DataArrays into a single DataArray. This may take a while...
print("Combining ND data")
ND_combined = xr.concat(ND_list, dim='D0_RR')
ND_combined.name = 'ND_combined_{}'.format(dataset_name)
# Add some metadata
ND_combined.attrs['DSD_interval'] = DSD_interval
ND_combined.attrs['strongwindQC'] = int(strongwindQC)
ND_combined.attrs['splashingQC'] = int(splashingQC)
ND_combined.attrs['marginQC'] = int(marginQC)
ND_combined.attrs['rainfallQC'] = int(rainfallQC)
ND_combined.attrs['rainonlyQC'] = int(rainonlyQC)
ND_combined.attrs['hailonlyQC'] = int(hailonlyQC)
ND_combined.attrs['graupelonlyQC'] = int(graupelonlyQC)
ND_combined.attrs['basicQC'] = int(basicQC)

print("Grouping by D0-RR and averaging")
# Group by D0-RR pairs:
ND_groups = ND_combined.groupby('D0_RR')
# Now average in each RR-D0 bin
ND_avg = ND_groups.mean(dim='D0_RR')
# Unstack the array so that D0 and RR are recovered as individual dimensions
ND_avg = ND_avg.unstack(dim='D0_RR')
# Rename the new dimensions back to what they were before because for some reason doing averages on a groupby object
# resets the names of the multiindex levels. Also do some reindexing and naming of coords.
# TODO change coordinates to midpoints of bins instead of edges
ND_avg = ND_avg.rename({'D0_RR_level_0': 'D0_idx', 'D0_RR_level_1': 'RR_idx'})
ND_avg = ND_avg.reindex({'D0_idx': range(D0_bins.size), 'RR_idx': range(RR_bins.size)})
ND_avg.coords['D0'] = ('D0_idx', D0_bins)
ND_avg.coords['RR'] = ('RR_idx', RR_bins)
ND_avg.name = 'SATP_ND_{}'.format(dataset_name)
# Add some metadata
ND_avg.attrs['DSD_interval'] = DSD_interval
ND_avg.attrs['strongwindQC'] = int(strongwindQC)
ND_avg.attrs['splashingQC'] = int(splashingQC)
ND_avg.attrs['marginQC'] = int(marginQC)
ND_avg.attrs['rainfallQC'] = int(rainfallQC)
ND_avg.attrs['rainonlyQC'] = int(rainonlyQC)
ND_avg.attrs['hailonlyQC'] = int(hailonlyQC)
ND_avg.attrs['graupelonlyQC'] = int(graupelonlyQC)
ND_avg.attrs['basicQC'] = int(basicQC)

# Dump SATP dataset to netCDF files with appropriate metadata
# Before dumping have to reset the MultiIndex because xarray doesn't yet allow for
# serialization of MultiIndex to netCDF files. This means we have to reconstruct the
# MultiIndex upon reading it back from the file if we want to use it to group and average
# as above. There's a function "reconstruct_MultiIndex" in pips_io.py for this purpose
ND_combined = ND_combined.reset_index('D0_RR')
ND_combined_ncfile_name = 'ND_combined_{}_{:d}s.nc'.format(dataset_name, int(DSD_interval))
ND_combined_ncfile_path = os.path.join(PIPS_dir, ND_combined_ncfile_name)
print("Dumping {}".format(ND_combined_ncfile_path))
ND_combined.to_netcdf(ND_combined_ncfile_path)

ND_avg_ncfile_name = 'ND_avg_{}_{:d}s.nc'.format(dataset_name, int(DSD_interval))
ND_avg_ncfile_path = os.path.join(PIPS_dir, ND_avg_ncfile_name)
print("Dumping {}".format(ND_avg_ncfile_path))
ND_avg.to_netcdf(ND_avg_ncfile_path)

# TODO: Add an "SATP" dictionary to the config file to control things like the file name, etc.


