# pyPIPS_meteograms.py
#
# This script fits a CG mu-lambda relation to the PIPS data and plots it along with the mu-lambda
# pairs
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
description = "Computes and plots CG relations from PIPS data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--ND-tag', dest='ND_tag', default='qc',
                    help='tag for ND variable in file (either qc or RB15)')

args = parser.parse_args()
ND_tag = args.ND_tag

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

# Get a list of the DSD fit netCDF data files that are present in the PIPS directory

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory

parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# print(DSD_fits_filenames)

# Open entire multi-file dataset
print("Opening dataswet {}".format(dataset_name))
parsivel_combined_ds = xr.open_mfdataset(parsivel_combined_filelist)
DSD_interval = parsivel_combined_ds.DSD_interval
PIPS_name = parsivel_combined_ds.probe_name
deployment_name = parsivel_combined_ds.deployment_name


# Get the untruncated and truncated moment fits MM246.
# TODO: update this for greater flexibility later
DSD_MM246 = parsivel_combined_ds['DSD_MM246']
DSD_TMM246 = parsivel_combined_ds['DSD_TMM246']

# Drop points where mu > 30 or lambda > 20000
DSD_MM246 = DSD_MM246.where(DSD_MM246.sel(parameter='lamda') < 20000.)
DSD_MM246 = DSD_MM246.where(DSD_MM246.sel(parameter='alpha') < 30.)
DSD_MM246 = DSD_MM246.dropna(dim='time', how='any')

DSD_TMM246 = DSD_TMM246.where(DSD_TMM246.sel(parameter='lamda') < 20000.)
DSD_TMM246 = DSD_TMM246.where(DSD_TMM246.sel(parameter='alpha') < 30.)
DSD_TMM246 = DSD_TMM246.dropna(dim='time', how='any')

# Fit polynomials
lamda_MM246 = DSD_MM246.sel(parameter='lamda') / 1000. # Get to mm^-1
mu_MM246 = DSD_MM246.sel(parameter='alpha')

MM246_poly_coeff, MM246_poly = dsd.calc_CG_polynomial(lamda_MM246, mu_MM246)
print("The MM246 polynomial coefficients are: ", MM246_poly_coeff)

# Plot the relationship on a scatterplot along with the Cao and Zhang relations
fig, ax = pm.plot_mu_lamda(lamda_MM246, mu_MM246, MM246_poly_coeff, MM246_poly, title='Untruncated')
plt.savefig(plot_dir + '/Untruncated_MM246_mu_lamda.png', dpi=200, bbox_inches='tight')

# Fit polynomials
lamda_TMM246 = DSD_TMM246.sel(parameter='lamda') / 1000. # Get to mm^-1
mu_TMM246 = DSD_TMM246.sel(parameter='alpha')

TMM246_poly_coeff, TMM246_poly = dsd.calc_CG_polynomial(lamda_TMM246, mu_TMM246)
print("The TMM246 polynomial coefficients are: ", TMM246_poly_coeff)

# Plot the relationship on a scatterplot along with the Cao and Zhang relations
pm.plot_mu_lamda(lamda_TMM246, mu_TMM246, TMM246_poly_coeff, TMM246_poly, title='Truncated')
plt.savefig(plot_dir + '/Truncated_MM246_mu_lamda.png', dpi=200, bbox_inches='tight')
