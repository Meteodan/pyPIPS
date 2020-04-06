# pyPIPS_meteograms.py
#
# This script plots a Dm-sigma scatterplot from PIPS data
# STOPPED HERE!
import os
import argparse
from datetime import datetime, timedelta
import matplotlib.colors as colors
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
description = "Plots one-to-one scatter plots from PIPS/radar data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--filter-RR', dest='filter_RR', type=float, default=None,
                    help='filter rainrate < # (mm)')
parser.add_argument('--filter-counts', dest='filter_counts', type=int, default=None,
                    help='filter particle counts < #')

args = parser.parse_args()
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

# Extract needed lists and variables from the radar_dict configuration dictionary
load_radar_at_PIPS = config.radar_config_dict.get('load_radar_at_PIPS', False)
save_radar_at_PIPS = config.radar_config_dict.get('save_radar_at_PIPS', False)
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
plot_retrieval = config.radar_config_dict.get('plot_retrieval', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_dir = config.radar_config_dict.get('radar_dir', None)
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Create the directory for the one2one plots if it doesn't exist
one2one_image_dir = os.path.join(plot_dir, 'one2one')
if not os.path.exists(one2one_image_dir):
    os.makedirs(one2one_image_dir)

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
radar_fields_at_PIPS_da = parsivel_combined_ds['{}_at_PIPS'.format(radar_name)]
dim_name = 'fields_{}'.format(radar_name)

# Drop points where mu > 30 or lambda > 20000
# parsivel_combined_ds = parsivel_combined_ds.where(
#     parsivel_combined_ds['DSD_TMM246'].sel(parameter='lamda') < 20000., drop=True)
# parsivel_combined_ds = parsivel_combined_ds.where(
#     parsivel_combined_ds['DSD_TMM246'].sel(parameter='alpha') < 30., drop=True)

# Plot one-to-one for various variables

# D0
ND = parsivel_combined_ds['ND_qc'].load()
ND_retr = parsivel_combined_ds['ND_retr'].load()
D0 = dsd.calc_D0_bin(ND) * 1000.
mu_retr = parsivel_combined_ds['mu_retr'].load()
lamda_retr = parsivel_combined_ds['lamda_retr'].load()
D0_retr = (3.67 + mu_retr) / lamda_retr
# D0_retr = parsivel_combined_ds['D0_retr'].load() # dsd.calc_D0_bin(ND_retr) * 1000.
D0_rad = radar_fields_at_PIPS_da.loc[{dim_name: 'D0'}].load()
RR_obs = parsivel_combined_ds['precipintensity'].load()
D0_ds = xr.Dataset({'D0_obs': D0, 'D0_retr': D0_retr, 'RR_obs': RR_obs})
D0_ds_rad = xr.Dataset({'D0_obs': D0, 'D0_rad': D0_rad, 'RR_obs': RR_obs})

axparams1 = {
    'var_lims': [0.0, 4.5],
    'col_field': 'RR_obs',
    'col_field_lims': [0.1, 200.],
    'norm': colors.LogNorm(vmin=0.1, vmax=200.),
    'alpha': 0.75,
    'markersize': 10,
    'markerstyle': 'o',
    'label_x': r'$D_0$ (obs; mm)',
    'label_y': r'$D_0$ (retrieved, disdrometer; mm)',
    'label_cb': r'$RR$ (mm h$^{-1}$)'
}

axparams2 = {
    'var_lims': [0.0, 4.5],
    'col_field': 'RR_obs',
    'col_field_lims': [0.1, 200.],
    'norm': colors.LogNorm(vmin=0.1, vmax=200.),
    'alpha': 0.75,
    'markersize': 10,
    'markerstyle': '+',
    'label_x': r'$D_0$ (obs; mm)',
    'label_y': r'$D_0$ (retrieved, radar; mm)',
    'label_cb': r'$RR$ (mm h$^{-1}$)',
    'stat_text_loc': [(0.1, 0.8), (0.1, 0.75)],
    'stat_labels': [r'$\rho_{{rr}}$: {:.2f}', r'Bias$_{{rr}}$: {:.2f}']
}

fig, ax = pm.plot_one2one(D0_ds, 'D0_obs', 'D0_retr', axparams1, add_colorbar=False)
fig, ax = pm.plot_one2one(D0_ds_rad, 'D0_obs', 'D0_rad', axparams2, fig=fig, ax=ax)
# plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=12.)

plt.savefig(one2one_image_dir +
            '/D0_one2one_{}_{}_{}_{}.png'.format(dataset_name, plot_tag, filter_RR_tag,
                                                 filter_counts_tag), dpi=200, bbox_inches='tight')


