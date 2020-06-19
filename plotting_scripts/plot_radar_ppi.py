# pyPIPS_meteograms.py
#
# This script plots radar PPI plots overlaid with the locations of the PIPS
import os
import argparse
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp
import pyPIPS.timemodule as tm

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
parser = argparse.ArgumentParser(description="Plots radar PPIs")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--plot-filtered-fields', dest='plot_filtered', default=False,
                    action='store_true',
                    help='Whether to also plot previously filtered dBZ and ZDR fields for the retrieval')
parser.add_argument('--input-tag', dest='input_tag', default=None,
                    help='Input nametag to determine which radar files to read in')

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
comp_radar = config.radar_config_dict.get('comp_radar', False)
clean_radar = config.radar_config_dict.get('comp_radar', False)
calc_dualpol = config.radar_config_dict.get('calc_dualpol', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', None)
radar_dir = config.radar_config_dict.get('radar_dir', None)
field_names = config.radar_config_dict.get('field_names', ['REF'])
if not calc_dualpol:
    field_names = ['REF']
el_req = config.radar_config_dict.get('el_req', 0.5)
radar_start_timestamp = config.radar_config_dict.get('radar_start_timestamp', None)
radar_end_timestamp = config.radar_config_dict.get('radar_end_timestamp', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

# Create the directory for the PPI plots if it doesn't exist
radar_ppi_image_dir = os.path.join(plot_dir, 'radar_ppi')
if not os.path.exists(radar_ppi_image_dir):
    os.makedirs(radar_ppi_image_dir)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# Read radar sweeps
if args.input_tag is None:
    radar_paths = glob(radar_dir + '/*{}*SUR.nc'.format(radar_name))
else:
    radar_paths = glob(radar_dir + '/*{}*SUR_{}.nc'.format(radar_name, args.input_tag))
radar_dict = radar.read_sweeps(radar_paths, radar_start_timestamp,
                               radar_end_timestamp, field_names=field_names, el_req=el_req,
                               radar_type=radar_type)

geo_locs = []
rad_locs = []
PIPS_names = []
for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    PIPS_name = parsivel_combined_ds.probe_name
    PIPS_names.append(PIPS_name)
    deployment_name = parsivel_combined_ds.deployment_name
    image_dir = os.path.join(radar_ppi_image_dir, deployment_name)
    if not os.path.exists(image_dir):
        os.makedirs(image_dir)
    geo_loc_str = parsivel_combined_ds.location
    geo_loc = list(map(np.float, geo_loc_str.strip('()').split(',')))
    geo_locs.append(geo_loc)
    rad_loc = radar.get_PIPS_loc_relative_to_radar(geo_loc, radar_dict['radarsweeplist'][0])
    rad_locs.append(rad_loc)

for radar_obj, sweeptime in zip(radar_dict['radarsweeplist'], radar_dict['sweeptimelist']):
    print(radar_obj.info())
    sweeptime_string = sweeptime.strftime(tm.timefmt3)
    figlist, axlist, fields_plotted = radar.plotsweep_pyART(radar_obj, sweeptime, PIPS_names,
                                                            geo_locs, rad_locs, field_names,
                                                            plot_filtered=args.plot_filtered)

    for fig, ax, field_name in zip(figlist, axlist, fields_plotted):
        PIPS_plot_name = '{}_{}_{}_{}_{}deg.png'.format(field_name, deployment_name,
                                                        sweeptime_string, radar_name, str(el_req))
        PIPS_plot_path = os.path.join(image_dir, PIPS_plot_name)
        fig.savefig(PIPS_plot_path, dpi=200, bbox_inches='tight')
        plt.close(fig)

