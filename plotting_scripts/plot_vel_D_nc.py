# pyPIPS_meteograms.py
#
# This script plots velocity-diameter histograms from the PIPS (netCDF version)
from __future__ import annotations

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import xarray as xr

import pyPIPS.parsivel_params as pp
import pyPIPS.plotmodule as pm
import pyPIPS.timemodule as tm
from pyPIPS import utils

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

diameter_bin_edges = pp.parsivel_parameters['diameter_bin_edges_mm']
fallspeed_bin_edges = pp.parsivel_parameters['fallspeed_bin_edges_mps']

# Parse the command line options
parser = argparse.ArgumentParser(description="Plots velocity-diameter histograms from PIPS data")
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--plot-config-path', dest='plot_config_path',
                    default='plot_config_default.py',
                    help='Location of the plot configuration file')
parser.add_argument('--plot-raw', action='store_true', dest='plot_raw',
                    help='plot raw velocity-diameter matrix')
parser.add_argument('--plot-qc', action='store_true', dest='plot_qc',
                    help='plot QC velocity-diameter matrix')
parser.add_argument('--QC-tag', dest='QC_tag', default=None,
                    help='Tag for QC variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--plot-full', action='store_true', dest='plot_full',
                    help='Plot full-deployment v-d matrix')
parser.add_argument('--normalize', action='store_true', dest='normalize',
                    help='Normalize full-deployment v-d matrix by total drop count')
parser.add_argument('--plot-series', action='store_true', dest='plot_series',
                    help='Plot time series of v-d matrix')
parser.add_argument('--cb-limits-full', nargs=2, metavar=('lower', 'upper'), type=float,
                    dest='cb_limits_full', default=[1, 2000],
                    help='color bar limits for drop number for full deployment plot')
parser.add_argument('--cb-limits-indv', nargs=2, metavar=('lower', 'upper'), type=float,
                    dest='cb_limits_indv', default=[1, 50],
                    help='color bar limits for drop number for individual DSDs')
parser.add_argument('--time-dim', dest='time_dim', default='time',
                    help='Name of the time dimension in the netCDF file')

args = parser.parse_args()

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
config = utils.import_all_from(args.case_config_path)
try:
    config = utils.import_all_from(args.case_config_path)
    utils.log("Successfully imported case configuration parameters!")
except Exception:
    utils.fatal(
        "Unable to import case configuration parameters! Aborting!")

# Dynamically import the plotting configuration file
utils.log(f"Plotting configuration file is {args.plot_config_path}")
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
PIPS_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)
conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', None)
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Create the directory for the plots if it doesn't exist
meteogram_image_dir = os.path.join(plot_dir, 'vel_D')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
if PIPS_filenames_nc:
    parsivel_combined_filenames = PIPS_filenames_nc
else:
    parsivel_combined_filenames = [
        f'parsivel_combined_{deployment_name}_{PIPS_name}_{int(requested_interval):d}s.nc'
            for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print(f"Reading {parsivel_combined_file}")  # noqa: T201
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)
    # Set up desired start and end times if they are not "None". Otherwise just use all times in
    # each file
    start_time = start_times[index]
    end_time = end_times[index]
    if start_time is not None and end_time is not None and args.time_dim == 'time':
        start_datetime = datetime.strptime(start_time, tm.timefmt3)
        end_datetime = datetime.strptime(end_time, tm.timefmt3)
        start_time = start_datetime.strftime(tm.timefmt2)
        end_time = end_datetime.strftime(tm.timefmt2)
        print(f"Extracting subset of times between {start_time} and {end_time}")  # noqa: T201
        parsivel_combined_ds = parsivel_combined_ds.sel({
            args.time_dim: slice(start_time, end_time)})

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
        QC_tag = args.QC_tag if args.QC_tag else 'qc'
        vd_matrix_qc = parsivel_combined_ds[f'VD_matrix_{QC_tag}']
        vd_matrix_list.append(vd_matrix_qc)
        tag_list.append(QC_tag)

    for tag, vd_matrix_da in zip(tag_list, vd_matrix_list):
        vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.)
        vel_D_image_dir = \
            os.path.join(plot_dir,
                         f'vel_D/{deployment_name}_{PIPS_name}_{int(requested_interval):d}s')
        if not os.path.exists(vel_D_image_dir):
            os.makedirs(vel_D_image_dir)

        axdict = {
            'diameter_bin_edges': diameter_bin_edges,
            'avg_diameter': avg_diameter,
            'fallspeed_bin_edges': fallspeed_bin_edges,
            'xlim': pc.PIPS_plotting_dict['diameter_range'],
            'ylim': pc.PIPS_plotting_dict['velocity_range'],
            'PIPS_name': PIPS_name
        }

        if args.plot_series:
            for t, ltime in enumerate(vd_matrix_da[args.time_dim].to_index()):
                time = int(ltime) if args.time_dim == 'relative_time' else ltime
                if parsivel_combined_ds['pcount_derived_qc'].loc[time] > 0:
                    if args.time_dim == 'time':
                        time_string = time.strftime(tm.timefmt3)
                    else:
                        time_string = f"{time:04d}"
                    print(f"Plotting for {PIPS_name} and time {time_string}")  # noqa: T201
                    axdict['time'] = time
                    axdict['cblim'] = args.cb_limits_indv
                    PSDdict = {
                        'vd_matrix_da': vd_matrix_da[t, :],
                        'DSD_interval': DSD_interval
                        # FIXME
                        # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
                    }
                    fig, ax = pm.plot_vel_D(axdict, PSDdict, parsivel_combined_ds['rho'].loc[time],
                                            time_dim=args.time_dim)
                    if args.time_dim == 'time':
                        time_string = time.strftime(tm.timefmt3)
                    else:
                        time_string = f'{t:04d}'
                    image_name = \
                        (f'{PIPS_name}_{deployment_name}_vel_D_{tag}_{int(DSD_interval):d}_'
                         f'{time_string}_t{t:04d}.png')
                    image_path = os.path.join(vel_D_image_dir, image_name)
                    fig.savefig(image_path, dpi=200, bbox_inches='tight')
                    plt.close(fig)

        if args.plot_full:
            vd_matrix_da_full = vd_matrix_da.sum(dim=args.time_dim)
            vd_matrix_da_full = vd_matrix_da_full.where(vd_matrix_da_full > 0.)
            print('Maximum bin drop count: ', vd_matrix_da_full.max())  # noqa: T201
            print(f"Plotting full-deployment v-d matrix for {deployment_name} and {PIPS_name}")  # noqa: T201
            if args.normalize:
                print("Normalizing by total drop count")  # noqa: T201
                axdict['cblim'] = (0., 0.05)
                vd_matrix_da_full = vd_matrix_da_full / vd_matrix_da_full.sum()  # noqa: PLR6104
                image_tag = 'full_norm'
            else:
                axdict['cblim'] = args.cb_limits_full
                image_tag = 'full'
            PSDdict = {
                'vd_matrix_da': vd_matrix_da_full,
                'DSD_interval': len(vd_matrix_da[args.time_dim]) * DSD_interval
                # FIXME
                # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
            }
            fig, ax = pm.plot_vel_D(axdict, PSDdict,
                                    parsivel_combined_ds['rho'].mean(dim=args.time_dim),
                                    time_dim=args.time_dim)
            image_name = \
                f'{PIPS_name}_{deployment_name}_vel_D_{tag}_{image_tag}.png'
            image_path = os.path.join(vel_D_image_dir, image_name)
            fig.savefig(image_path, dpi=200, bbox_inches='tight')
            plt.close(fig)

