# pyPIPS_meteograms.py
#
# This script plots a Dm-sigma scatterplot from PIPS data
# STOPPED HERE!
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
description = "Plots Dm-sigma scatter and CG prediction from PIPS data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--coefficients', nargs=3, metavar=('c1', 'c2', 'c3'), type=float,
                    dest='coefficients',
                    help='coefficients of mu-lambda polynomial (in decreasing order of exponent')
parser.add_argument('--plot-SATP', action='store_true', dest='plot_SATP',
                    help='plot for the SATP-filtered dataset')
parser.add_argument('--filter-RR', dest='filter_RR', type=float, default=None,
                    help='filter rainrate < # (mm)')
parser.add_argument('--filter-counts', dest='filter_counts', type=int, default=None,
                    help='filter particle counts < #')

args = parser.parse_args()
if args.plot_SATP:
    filter_RR = None
    filter_counts = None
else:
    filter_RR = args.filter_RR
    filter_counts = args.filter_counts

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
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)


if not args.plot_SATP:
    # Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
    parsivel_combined_filenames = [
        'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name,
                                                  int(requested_interval))
        for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
    parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf)
                                  for pcf in parsivel_combined_filenames]

    # print(DSD_fits_filenames)

    # Open entire multi-file dataset
    print("Opening dataset {}".format(dataset_name))
    parsivel_combined_ds = xr.open_mfdataset(parsivel_combined_filelist)
    plot_tag = 'no_SATP'
    dim = 'time'
else:
    parsivel_combined_filename = 'ND_avg_{}_{:d}s.nc'.format(dataset_name, int(requested_interval))
    parsivel_combined_filepath = os.path.join(PIPS_dir, parsivel_combined_filename)
    parsivel_combined_ds = xr.open_dataset(parsivel_combined_filepath)
    plot_tag = 'SATP'
    dim = 'D0_RR'

# Filter by RR and pcount if desired
filter_RR_tag = 'no_RR_filter'
filter_counts_tag = 'no_counts_filter'
if filter_RR:
    parsivel_combined_ds = parsivel_combined_ds.where(
        parsivel_combined_ds['precipintensity'] >= filter_RR, drop=True)
    filter_RR_tag = 'RR_filter_{:d}mm'.format(int(filter_RR))
if filter_counts:
    parsivel_combined_ds = parsivel_combined_ds.where(
        parsivel_combined_ds['pcount'] >= filter_counts, drop=True)
    filter_counts_tag = 'counts_filter_{:d}'.format(filter_counts)

DSD_interval = parsivel_combined_ds.DSD_interval

# Drop points where mu > 30 or lambda > 20000
# parsivel_combined_ds = parsivel_combined_ds.where(
#     parsivel_combined_ds['DSD_TMM246'].sel(parameter='lamda') < 20000., drop=True)
# parsivel_combined_ds = parsivel_combined_ds.where(
#     parsivel_combined_ds['DSD_TMM246'].sel(parameter='alpha') < 30., drop=True)

lamda_CG = np.linspace(0., 20., num=1000.)
mu_CG = dsd.calc_mu_lamda(lamda_CG, args.coefficients)

Dm43_obs = parsivel_combined_ds['Dm43'] * 1000.  # Get to mm
sigma_obs = parsivel_combined_ds['sigma'] * 1000.  # Get to mm

Dm43_CG = (4. + mu_CG) / lamda_CG
sigma_CG = np.sqrt((4. + mu_CG) / lamda_CG**2.)

# Plot the sigma-Dm scatterplot and C-G prediction
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(Dm43_obs, sigma_obs, color='k')
ax.plot(Dm43_CG, sigma_CG, color='b')
ax.set_xlim(0.0, 5.0)
ax.set_ylim(0.0, 2.5)
ax.set_xlabel(r'$D_m$ (mm)')
ax.set_ylabel(r'$\sigma_m$ (mm)')

plt.savefig(plot_dir +
            '/sigma_Dm_{}_{}_{}_{}.png'.format(dataset_name, plot_tag, filter_RR_tag,
                                               filter_counts_tag), dpi=200, bbox_inches='tight')
