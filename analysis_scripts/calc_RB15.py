# calc_RB15.py
#
# This script applies the correction procedure of Raupach and Berne (2015) to the PIPS parsivel
# observations.
#
# Reference:
# Raupach, T. H., & Berne, A. (2015). Correction of raindrop size distributions measured by
# Parsivel disdrometers, using a two-dimensional video disdrometer as a reference.
# Atmospheric Measurement Techniques, 8(1), 343–365. https://doi.org/10.5194/amt-8-343-2015
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
fall_bins_edges = np.append(min_fall_bins, max_fall_bins[-1])

# Parse the command line options
description = "Corrects observed DSDs using the procedure of Raupach and Berne (2015)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')

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

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in
                              parsivel_combined_filenames]

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_ds = xr.load_dataset(parsivel_combined_file)

    DSD_interval = parsivel_ds.DSD_interval
    PIPS_name = parsivel_ds.probe_name
    deployment_name = parsivel_ds.deployment_name

    # Make some index coordinates out of existing dimensions. This is a limitation of xarray right
    # now that you can't index using coordinates that aren't index coordinates.
    parsivel_ds.coords['fallspeed_bin'] = ('fallspeed_bin', parsivel_ds.coords['fallspeed'].data)
    parsivel_ds.coords['diameter_bin'] = ('diameter_bin', parsivel_ds.coords['diameter'].data)

    # We want the raw VD matrix here.
    vd_matrix = parsivel_ds['VD_matrix']
    # Rebin the velocities into evenly spaced bins of 0.1 m/s
    interval = 0.1
    print("Rebinning velocities.")
    vd_matrix_rebinned = pips.reindex_velocity_bins(vd_matrix, interval)
    # Calculate rain terminal velocity as a function of time and diameter
    # Uses the Atlas (1973) formula with the density correction of Foote and Dutoit
    vt_rain = pips.calc_empirical_fallspeed(vd_matrix_rebinned.coords['diameter'], correct_rho=True,
                                            rho=parsivel_ds['rho'])
    vt_rain_da = xr.DataArray(vt_rain,
                              coords={
                                  'time': ('time', vd_matrix_rebinned.coords['time'].data),
                                  'diameter_bin': ('diameter_bin',
                                                   vd_matrix_rebinned.coords['diameter'].data)
                              },
                              dims=['time', 'diameter_bin'])
    # Shift the velocities in each diameter bin of each DSD such that the mean matches
    # the expected terminal velocity
    print("Shifting velocities.")
    vd_matrix_rebinned_shifted = pips.shift_mean_velocity(vd_matrix_rebinned, vt_rain_da)
    # Then collect the counts into the original parsivel velocity bins
    vd_matrix_shifted = pips.rebin_to_parsivel(vd_matrix_rebinned_shifted)

    print("Computing ND for shifted VD matrix.")
    # Compute ND from the shifted VD matrix
    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=parsivel_ds['rho'])

    vd_matrix_shifted = vd_matrix_shifted.where(vd_matrix_shifted > 0.)
    ND_RB15_vshift = pips.calc_ND(vd_matrix_shifted, fallspeed_spectrum, DSD_interval)

    # Add the new variables to the parsivel_ds
    parsivel_ds = pipsio.combine_parsivel_data(parsivel_ds, vd_matrix_shifted,
                                               name='VD_matrix_RB15_vshift')
    parsivel_ds = pipsio.combine_parsivel_data(parsivel_ds, ND_RB15_vshift,
                                               name='ND_RB15_vshift')

    print("Computing corrected ND.")
    # Compute the ND correction from RB15 and add the new corrected ND array
    parsivel_ds = pips.correct_ND_RB15(parsivel_ds)

    # Save dataset back to file
    print("Dumping {}".format(parsivel_combined_file))
    parsivel_ds.to_netcdf(parsivel_combined_file)
