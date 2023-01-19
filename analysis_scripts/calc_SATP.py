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
parser.add_argument('--compute-rr', action='store_true', dest='compute_rr',
                    help='Compute rain-rate from ND bins instead of using internal parsivel values')
parser.add_argument('--ND-tag', dest='ND_tag', default=None,
                    help='Tag for ND variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')

args = parser.parse_args()
if not args.ND_tag:
    ND_tag = ''
else:
    ND_tag = '_{}'.format(args.ND_tag)

# Dynamically import the case configuration file
utils.log("Case config file is {}".format(args.case_config_path))
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Extract needed lists and variables from PIPS_IO_dict configuration dictionary
dataset_name = config.PIPS_IO_dict.get('dataset_name', None)
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
parsivel_combined_filenames = config.PIPS_IO_dict['PIPS_filenames_nc']
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

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

ND_list = []

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    ND = parsivel_combined_ds['ND{}'.format(ND_tag)]
    ND = ND.where(ND > 0.)
    # Drop all entries without DSDs
    ND = ND.dropna(dim='time', how='all')
    # # Compute binned number concentration
    # fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
    #                                                   rho=parsivel_combined_ds['rho'])
    # vd_matrix_da = parsivel_combined_ds['VD_matrix_qc']
    # vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    # ND2 = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    # ND2 = ND2.where(ND2 > 0.)
    # # Drop all entries without DSDs
    # ND2 = ND2.dropna(dim='time', how='all')
    # print("ND2 = ", ND2)

    if args.compute_rr:
        # Compute rainrate using empirical fallspeed curve
        # TODO: allow for the use of the measured fallspeeds in the rainrate calculation.
        # First, see if this has already been computed
        try:
            rainrate = parsivel_combined_ds['rainrate_derived{}'.format(ND_tag)]
        except KeyError:
            fallspeeds_emp = pips.calc_empirical_fallspeed(avg_diameter, correct_rho=True,
                                                           rho=parsivel_combined_ds['rho'])
            rainrate_bin = ((6. * 10.**-4.) * np.pi * fallspeeds_emp * avg_diameter**3. * ND *
                            bin_width)
            rainrate = rainrate_bin.sum(dim='diameter_bin')
    else:
        rainrate = parsivel_combined_ds['precipintensity']

    rainrate = rainrate.loc[ND.indexes['time']]
    RR_ind = (rainrate <= RR_bins[-1]) & (rainrate >= RR_bins[0])
    # rainrate = rainrate.where(RR_ind)

    # Compute D0 (in mm)
    # TODO: change this to only calculate if it is not already present in the file
    D0 = dsd.calc_D0_bin(ND) * 1000.
    D0_ind = (D0 <= D0_bins[-1]) & (D0 >= D0_bins[0])
    # D0 = D0.where(D0_ind)

    # Also mask out ND for D0 and RR outside of range
    D0 = D0.where((D0_ind & RR_ind), drop=True)
    rainrate = rainrate.where((D0_ind & RR_ind), drop=True)
    ND = ND.where((D0_ind & RR_ind), drop=True)

    if ND.sizes['time'] > 0:
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
# print(ND_combined)
ND_combined.name = 'ND_combined_{}{}'.format(dataset_name, ND_tag)

print("Grouping by D0-RR and averaging")
# Group by D0-RR pairs:
ND_groups = ND_combined.groupby('D0_RR')
# Now average in each RR-D0 bin
ND_avg = ND_groups.mean(dim='D0_RR')
# print(ND_avg)
# TODO: Modify the following to keep D0_RR index but to also add the D0_idx and RR_idx dimensions
# Along with coordinates. Can't add new dimensions to a DataArray so need to make this a Dataset
# For now, uncomment this out so that we are just keeping the combined D0_RR index.

# print(ND_avg)
# Unstack the array so that D0 and RR are recovered as individual dimensions
# ND_avg = ND_avg.unstack(dim='D0_RR')
# print(ND_avg)
# Rename the new dimensions back to what they were before because for some reason doing averages on
# a groupby object
# resets the names of the multiindex levels. Also do some reindexing and naming of coords.
# TODO change coordinates to midpoints of bins instead of edges
# ND_avg = ND_avg.rename({'D0_RR_level_0': 'D0_idx', 'D0_RR_level_1': 'RR_idx'})
# Never mind, don't reindex, keep the combined D0_RR index
# ND_avg = ND_avg.reindex({'D0_idx': range(D0_bins.size), 'RR_idx': range(RR_bins.size)})
# ND_avg.coords['D0_idx'] = ('D0_idx', range(D0_bins.size))
# ND_avg.coords['RR_idx'] = ('RR_idx', range(RR_bins.size))

ND_avg.name = 'SATP_ND_{}{}'.format(dataset_name, ND_tag)
# Dump SATP dataset to netCDF files with appropriate metadata
# Before dumping have to reset the MultiIndex because xarray doesn't yet allow for
# serialization of MultiIndex to netCDF files. This means we have to reconstruct the
# MultiIndex upon reading it back from the file if we want to use it to group and average
# as above. There's a function "reconstruct_MultiIndex" in pips_io.py for this purpose
ND_combined = ND_combined.reset_index('D0_RR')
ND_combined = ND_combined.to_dataset()
ND_combined.attrs['DSD_interval'] = DSD_interval
for attr_key, attr_val in parsivel_combined_ds.attrs.items():
    if 'QC' in attr_key:
        ND_combined.attrs[attr_key] = attr_val
ND_combined.coords['D0_bins'] = D0_bins
ND_combined.coords['RR_bins'] = RR_bins
ND_combined_ncfile_name = 'ND_combined_{}{}_{:d}s.nc'.format(dataset_name, ND_tag,
                                                             int(DSD_interval))
ND_combined_ncfile_path = os.path.join(PIPS_dir, ND_combined_ncfile_name)
print("Dumping {}".format(ND_combined_ncfile_path))
ND_combined.to_netcdf(ND_combined_ncfile_path)

ND_avg = ND_avg.reset_index('D0_RR')
ND_avg = ND_avg.to_dataset()
ND_avg.attrs['DSD_interval'] = DSD_interval
for attr_key, attr_val in parsivel_combined_ds.attrs.items():
    if 'QC' in attr_key:
        ND_avg.attrs[attr_key] = attr_val
ND_avg.coords['D0_bins'] = D0_bins
ND_avg.coords['RR_bins'] = RR_bins
ND_avg_ncfile_name = 'ND_avg_{}{}_{:d}s.nc'.format(dataset_name, ND_tag, int(DSD_interval))
ND_avg_ncfile_path = os.path.join(PIPS_dir, ND_avg_ncfile_name)
print("Dumping {}".format(ND_avg_ncfile_path))
ND_avg.to_netcdf(ND_avg_ncfile_path)
