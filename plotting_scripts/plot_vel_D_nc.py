# pyPIPS_meteograms.py
#
# This script plots velocity-diameter histograms from the PIPS (netCDF version)
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
parser = argparse.ArgumentParser(description="Plots velocity-diameter histograms from PIPS data")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config.py', help='Location of the plot configuration file')
parser.add_argument('--plot-raw', action='store_true', dest='plot_raw',
                    help='plot raw velocity-diameter matrix')
parser.add_argument('--plot-qc', action='store_true', dest='plot_qc',
                    help='plot QC velocity-diameter matrix')
parser.add_argument('--plot-full', action='store_true', dest='plot_full',
                    help='Plot full-deployment v-d matrix')
parser.add_argument('--normalize', action='store_true', dest='normalize',
                    help='Normalize full-deployment v-d matrix by total drop count')
parser.add_argument('--plot-series', action='store_true', dest='plot_series',
                    help='Plot time series of v-d matrix')

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

# Create the directory for the plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'vel_D')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
parsivel_combined_filenames = [
    'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name, int(requested_interval))
    for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]
print(parsivel_combined_filenames)

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    vd_matrix_list = []
    tag_list = []

    if args.plot_raw:
        vd_matrix_raw = parsivel_combined_ds['VD_matrix']
        vd_matrix_list.append(vd_matrix_raw)
        tag_list.append('raw')
    if args.plot_qc:
        vd_matrix_qc = parsivel_combined_ds['VD_matrix_qc']
        vd_matrix_list.append(vd_matrix_qc)
        tag_list.append('qc')

    for tag, vd_matrix_da in zip(tag_list, vd_matrix_list):
        vel_D_image_dir = os.path.join(plot_dir,
                                       'vel_D/{}_{}_{:d}s'.format(deployment_name, PIPS_name,
                                                                  int(requested_interval)))
        if not os.path.exists(vel_D_image_dir):
            os.makedirs(vel_D_image_dir)

        axdict = {
            'min_diameter': min_diameter,
            'avg_diameter': avg_diameter,
            'min_fall_bins': min_fall_bins,
            'xlim': pc.PIPS_plotting_dict['diameter_range'],
            'ylim': pc.PIPS_plotting_dict['velocity_range'],
            'PIPS_name': PIPS_name
        }

        if args.plot_series:
            for t, time in enumerate(vd_matrix_da['time'].to_index()):
                if parsivel_combined_ds['pcount_derived_qc'].loc[time] > 0:
                    print("Plotting for {} and time {}".format(PIPS_name,
                                                               time.strftime(tm.timefmt3)))
                    axdict['time'] = time
                    axdict['cblim'] = (1, 50)
                    PSDdict = {
                        'vd_matrix_da': vd_matrix_da[t, :],
                        'DSD_interval': DSD_interval
                        # FIXME
                        # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
                    }
                    fig, ax = pm.plot_vel_D(axdict, PSDdict, parsivel_combined_ds['rho'].loc[time])
                    image_name = \
                        '{}_{}_vel_D_{}_{:d}_{}_t{:04d}.png'.format(PIPS_name, deployment_name, tag,
                                                                    int(DSD_interval),
                                                                    time.strftime(tm.timefmt3), t)
                    image_path = os.path.join(vel_D_image_dir, image_name)
                    fig.savefig(image_path, dpi=200, bbox_inches='tight')
                    plt.close(fig)

        if args.plot_full:
            vd_matrix_da_full = vd_matrix_da.sum(dim='time')
            vd_matrix_da_full = vd_matrix_da_full.where(vd_matrix_da_full > 0.)
            print("Plotting full-deployment v-d matrix for {} and {}".format(deployment_name,
                                                                             PIPS_name))
            if args.normalize:
                print("Normalizing by total drop count")
                axdict['cblim'] = (0., 1.)
                vd_matrix_da_full = vd_matrix_da_full / vd_matrix_da_full.sum()
                image_tag = 'full_norm'
            else:
                axdict['cblim'] = (1, 1500)
                image_tag = 'full'
            PSDdict = {
                'vd_matrix_da': vd_matrix_da_full,
                'DSD_interval': len(vd_matrix_da['time']) * DSD_interval
                # FIXME
                # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
            }
            fig, ax = pm.plot_vel_D(axdict, PSDdict, parsivel_combined_ds['rho'].mean(dim='time'))
            image_name =  \
                '{}_{}_vel_D_{}_{}.png'.format(PIPS_name, deployment_name, tag, image_tag)
            image_path = os.path.join(vel_D_image_dir, image_name)
            fig.savefig(image_path, dpi=200, bbox_inches='tight')
            plt.close(fig)

