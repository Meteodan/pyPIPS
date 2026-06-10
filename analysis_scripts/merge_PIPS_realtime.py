#!/usr/bin/env python
"""
merge_PIPS_realtime.py

This script merges real-time PIPS netCDF files with card-based netCDF files to fill data gaps.
It reads real-time data files (onesec, ND, spectrum, telegram) and merges them with
corresponding card-derived data to create complete datasets.

Based on the functionality from merge_PIPS_realtime_nc.ipynb
"""
from __future__ import annotations

import argparse
import glob
import os
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as dates

import pyPIPS.PIPS as pips
import pyPIPS.pips_io as pipsio
import pyPIPS.timemodule as tm
import pyPIPS.plotmodule as pm
from pyPIPS import thermolib as thermo
from pyPIPS import utils


def get_files(file_list, starttime, endtime, ftype='onesec'):
    """
    Find files within specified time range.

    Parameters
    ----------
    file_list : list
        List of filenames to search through
    starttime : datetime
        Start time for file selection
    endtime : datetime
        End time for file selection
    ftype : str
        File type prefix (e.g., 'onesec', 'ND', 'spectrum', 'telegram')

    Returns
    -------
    list
        Filtered and sorted list of files within time range
    """
    if not file_list:
        return []

    len_prefix = 8 + len(ftype)

    day_str = [(starttime + timedelta(days=i)).strftime("%Y%m%d")
               for i in range((endtime - starttime).days + 1)]

    # Find all PIPS files with days between starttime and endtime
    file_list = [file_name for file_name in file_list if any(day in file_name for day in day_str)]

    if file_list:
        # Sort files by date, then find nearest indices for all the dates
        sorted_files = sorted(file_list,
                              key=lambda f: datetime.strptime(f[len_prefix:len_prefix + 14],
                                                              '%Y%m%d%H%M%S'))
        starttimes = [datetime.strptime(f[len_prefix:len_prefix + 14], '%Y%m%d%H%M%S')
                      for f in sorted_files]
        endtimes = [datetime.strptime(f[len_prefix + 15:len_prefix + 29], '%Y%m%d%H%M%S')
                    for f in sorted_files]
        _, idx1 = min((abs(val - starttime), idx) for (idx, val) in enumerate(starttimes))
        _, idx2 = min((abs(val - endtime), idx) for (idx, val) in enumerate(endtimes))
        file_list = sorted_files[idx1:idx2 + 1]

    return file_list


def correct_realtime_dewpoint(onesec_rt_ds):
    """
    Correct dewpoint calculation for real-time data.

    Real-time data incorrectly uses fasttemp instead of slowtemp.
    This function recalculates dewpoint and derived RH.

    Parameters
    ----------
    onesec_rt_ds : xr.Dataset
        Real-time one-second dataset

    Returns
    -------
    xr.Dataset
        Corrected dataset with proper dewpoint and RH_derived
    """
    RH = onesec_rt_ds['RH']
    SlowT = onesec_rt_ds['slowtemp']

    # Calculate correct dewpoint using slowtemp
    dewpoint = 243.04 * (np.log(RH / 100.) + ((17.625 * SlowT) / (243.04 + SlowT))) / \
                        (17.625 - np.log(RH / 100.) - ((17.625 * SlowT) / (243.04 + SlowT)))

    # Recalculate RH_derived using fasttemp and corrected dewpoint
    fasttemp = onesec_rt_ds['fasttemp']
    RH_derived = 100. * (np.exp((17.625 * dewpoint) / (243.04 + dewpoint)) /
                         np.exp((17.625 * fasttemp) / (243.04 + fasttemp)))

    onesec_rt_ds = onesec_rt_ds.copy()
    onesec_rt_ds['dewpoint'] = dewpoint
    onesec_rt_ds['RH_derived'] = RH_derived

    return onesec_rt_ds


def merge_datasets(rt_da, card_da, all_times):
    """
    Merge real-time and card datasets to fill gaps.
    Uses notebook approach: real-time data has priority.

    Parameters
    ----------
    rt_da : xr.DataArray
        Real-time data array
    card_da : xr.DataArray
        Card-derived data array
    all_times : pd.DatetimeIndex
        Complete time index for output

    Returns
    -------
    xr.DataArray
        Merged data array
    """
    # Use real-time data first, fill gaps with card data (notebook approach)
    merged_da = rt_da.combine_first(card_da)

    # Reindex to full time range
    merged_da_full = merged_da.reindex({'time': all_times})

    return merged_da_full


def create_diagnostic_plots(all_onesec_times, onesec_rt_ds, onesec_card_ds,
                          all_tensec_times, parsivel_combined_merged_full_ds,
                          PIPS_name, deployment_name, plot_dir, ND_rt_ds=None,
                          parsivel_combined_card_ds=None):
    """
    Create diagnostic plots showing missing data patterns (notebook approach).

    Parameters
    ----------
    all_onesec_times : pd.DatetimeIndex
        Complete one-second time range
    onesec_rt_ds : xr.Dataset or None
        Real-time one-second dataset
    onesec_card_ds : xr.Dataset
        Card one-second dataset
    all_tensec_times : pd.DatetimeIndex
        Complete ten-second time range
    parsivel_combined_merged_full_ds : xr.Dataset
        Merged parsivel dataset
    PIPS_name : str
        PIPS instrument name
    deployment_name : str
        Deployment name
    plot_dir : str
        Directory to save plots
    ND_rt_ds : xr.Dataset, optional
        Real-time ND dataset
    parsivel_combined_card_ds : xr.Dataset, optional
        Original card parsivel dataset (before merging)
    """
    utils.log("Creating diagnostic plots...")

    # Common DSD plotting parameters
    def setup_dsd_plot_params():
        """Setup common parameters for DSD meteogram plotting."""
        diameter_bins = pips.diameter_edges
        log_ND_params = {
            'type': 'pcolor',
            'vlimits': [-1.0, 3.0],
            'clabel': r'log[N (m$^{-3}$ mm$^{-1}$)]',
            'cmap': 'viridis'
        }
        return diameter_bins, log_ND_params

    def format_time_axis(ax):
        """Format time axis consistently across plots."""
        ax.xaxis.set_major_locator(dates.HourLocator(interval=1))
        ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
        ax.xaxis.set_minor_locator(dates.MinuteLocator(interval=15))
        plt.xticks(rotation=45)

    # Create missing times arrays
    if onesec_rt_ds is not None:
        missing_times_rt = np.array([False if time in onesec_rt_ds.indexes['time']
                                    else True for time in all_onesec_times])
    else:
        missing_times_rt = np.ones(len(all_onesec_times), dtype=bool)

    if 'missing_times' in onesec_card_ds:
        missing_times_card = onesec_card_ds['missing_times'].values.astype('bool')
    else:
        missing_times_card = np.array([False if time in onesec_card_ds.indexes['time']
                                      else True for time in all_onesec_times])

    missing_times_both = missing_times_rt & missing_times_card

    # Plot 1: Missing times comparison
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(all_onesec_times, missing_times_rt, label='Real-time data missing', alpha=0.7)
    ax.plot(all_onesec_times, missing_times_card, label='Card data missing', alpha=0.7)
    ax.set_xlabel('Time')
    ax.set_ylabel('Data Missing')
    ax.set_title(f'Missing Data Comparison - {PIPS_name} ({deployment_name})')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    format_time_axis(ax)

    plot_filename = f'missing_data_comparison_{deployment_name}_{PIPS_name}.png'
    plot_path = os.path.join(plot_dir, plot_filename)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    utils.log(f"Saved diagnostic plot: {plot_path}")

    # Plot 2: Missing times for both datasets
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(all_onesec_times, missing_times_both, label='Both datasets missing', color='red')
    ax.set_xlabel('Time')
    ax.set_ylabel('Data Missing from Both')
    ax.set_title(f'Times Missing from Both Datasets - {PIPS_name} ({deployment_name})')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    format_time_axis(ax)

    plot_filename = f'missing_data_both_{deployment_name}_{PIPS_name}.png'
    plot_path = os.path.join(plot_dir, plot_filename)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    utils.log(f"Saved diagnostic plot: {plot_path}")

    # Plot 3: Temperature and dewpoint from merged data
    if 'slowtemp' in parsivel_combined_merged_full_ds.data_vars and 'dewpoint' in parsivel_combined_merged_full_ds.data_vars:
        fig, ax = plt.subplots(figsize=(12, 6))

        plot_times = parsivel_combined_merged_full_ds['time'].to_index().to_pydatetime()
        temp_data = parsivel_combined_merged_full_ds['slowtemp'].values
        dewpoint_data = parsivel_combined_merged_full_ds['dewpoint'].values

        ax.plot(plot_times, temp_data, label='Temperature', color='red')
        ax.plot(plot_times, dewpoint_data, label='Dewpoint', color='blue')
        ax.set_xlabel('Time')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title(f'Merged Temperature/Dewpoint - {PIPS_name} ({deployment_name})')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        format_time_axis(ax)

        plot_filename = f'merged_temperature_{deployment_name}_{PIPS_name}.png'
        plot_path = os.path.join(plot_dir, plot_filename)
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        utils.log(f"Saved diagnostic plot: {plot_path}")

    # Plot 4: DSD meteogram from merged data
    if 'ND' in parsivel_combined_merged_full_ds.data_vars:
        utils.log("Creating DSD meteogram for merged data...")

        fig, ax = plt.subplots(figsize=(12, 8))

        # Get ND data and time coordinates
        ND_data = parsivel_combined_merged_full_ds['ND']
        plot_times_tmp = ND_data.coords['time'].to_index().to_pydatetime()

        # Prepend additional time for pcolormesh edges (notebook approach)
        time_delta = plot_times_tmp[1] - plot_times_tmp[0]
        plot_times = np.insert(plot_times_tmp, 0, plot_times_tmp[0] - time_delta)
        plot_times = [plot_times]

        # Get diameter bin information and setup plot parameters
        diameter_bins, log_ND_params = setup_dsd_plot_params()

        # Prepare data for plotting
        ND_arr = ND_data.values.T
        log_ND_arr = np.log10(np.maximum(ND_arr, 1e-6))  # Avoid log(0)

        # Plot DSD meteogram
        fields_to_plot = [log_ND_arr]
        field_parameters = [log_ND_params]
        yvals = [diameter_bins] * len(fields_to_plot)

        ax = pm.plotmeteogram(ax, plot_times, fields_to_plot, field_parameters,
                             yvals=yvals, plot_data_bounds=False)

        ax.set_xlabel('Time (UTC)')
        ax.set_ylabel('Diameter (mm)')
        ax.set_title(f'Merged DSD Meteogram - {PIPS_name} ({deployment_name})')
        format_time_axis(ax)

        plot_filename = f'merged_DSD_meteogram_{deployment_name}_{PIPS_name}.png'
        plot_path = os.path.join(plot_dir, plot_filename)
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        utils.log(f"Saved diagnostic plot: {plot_path}")

    # Plot 5: DSD meteogram from real-time data (if available)
    if onesec_rt_ds is not None and ND_rt_ds is not None and 'ND' in ND_rt_ds.data_vars:
        utils.log("Creating DSD meteogram for real-time data...")

        try:
            fig, ax = plt.subplots(figsize=(12, 8))

            # Get real-time ND data
            ND_rt_data = ND_rt_ds['ND']
            plot_times_tmp = ND_rt_data.coords['time'].to_index().to_pydatetime()

            if len(plot_times_tmp) > 1:
                time_delta = plot_times_tmp[1] - plot_times_tmp[0]
                plot_times = np.insert(plot_times_tmp, 0, plot_times_tmp[0] - time_delta)
                plot_times = [plot_times]

                # Get diameter bins and setup plot parameters
                diameter_bins, log_ND_params = setup_dsd_plot_params()

                # Prepare data
                ND_arr = ND_rt_data.values.T
                log_ND_arr = np.log10(np.maximum(ND_arr, 1e-6))

                fields_to_plot = [log_ND_arr]
                field_parameters = [log_ND_params]
                yvals = [diameter_bins] * len(fields_to_plot)

                ax = pm.plotmeteogram(ax, plot_times, fields_to_plot, field_parameters,
                                     yvals=yvals, plot_data_bounds=False)

                ax.set_xlabel('Time (UTC)')
                ax.set_ylabel('Diameter (mm)')
                ax.set_title(f'Real-time DSD Meteogram - {PIPS_name} ({deployment_name})')
                format_time_axis(ax)

                plot_filename = f'realtime_DSD_meteogram_{deployment_name}_{PIPS_name}.png'
                plot_path = os.path.join(plot_dir, plot_filename)
                plt.tight_layout()
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.close()
                utils.log(f"Saved diagnostic plot: {plot_path}")
            else:
                plt.close()
                utils.log("Insufficient real-time data for DSD meteogram")
        except Exception as e:
            plt.close()
            utils.log(f"Error creating real-time DSD meteogram: {e}")

    # Plot 6: DSD meteogram from card data
    if parsivel_combined_card_ds is not None and 'ND' in parsivel_combined_card_ds.data_vars:
        utils.log("Creating DSD meteogram for card data...")

        try:
            fig, ax = plt.subplots(figsize=(12, 8))

            # Get original card ND data (before merging)
            ND_card_data = parsivel_combined_card_ds['ND']
            plot_times_tmp = ND_card_data.coords['time'].to_index().to_pydatetime()

            if len(plot_times_tmp) > 1:
                time_delta = plot_times_tmp[1] - plot_times_tmp[0]
                plot_times = np.insert(plot_times_tmp, 0, plot_times_tmp[0] - time_delta)
                plot_times = [plot_times]

                # Get diameter bins and setup plot parameters
                diameter_bins, log_ND_params = setup_dsd_plot_params()

                # Prepare data
                ND_arr = ND_card_data.values.T
                log_ND_arr = np.log10(np.maximum(ND_arr, 1e-6))

                fields_to_plot = [log_ND_arr]
                field_parameters = [log_ND_params]
                yvals = [diameter_bins] * len(fields_to_plot)

                ax = pm.plotmeteogram(ax, plot_times, fields_to_plot, field_parameters,
                                     yvals=yvals, plot_data_bounds=False)

                ax.set_xlabel('Time (UTC)')
                ax.set_ylabel('Diameter (mm)')
                ax.set_title(f'Card DSD Meteogram - {PIPS_name} ({deployment_name})')
                format_time_axis(ax)

                plot_filename = f'card_DSD_meteogram_{deployment_name}_{PIPS_name}.png'
                plot_path = os.path.join(plot_dir, plot_filename)
                plt.tight_layout()
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.close()
                utils.log(f"Saved diagnostic plot: {plot_path}")
            else:
                plt.close()
                utils.log("Insufficient card data for DSD meteogram")
        except Exception as e:
            plt.close()
            utils.log(f"Error creating card DSD meteogram: {e}")
    else:
        utils.log("No card ND data available for DSD meteogram")


def main():
    # Parse command line arguments
    description = "Merge real-time and card-based PIPS netCDF files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                        help='The path to the case configuration file')
    parser.add_argument('--realtime-dir', dest='realtime_dir', required=True,
                        help='Directory containing real-time netCDF files')
    parser.add_argument('--output-dir', dest='output_dir', default=None,
                        help='Output directory for merged files (defaults to PIPS_dir from config)')
    parser.add_argument('--output-tag', dest='output_tag', default='',
                        help='Tag for output files (default: none)')
    parser.add_argument('--diagnostic-plots', dest='diagnostic_plots', action='store_true',
                        help='Create diagnostic plots showing missing data patterns')
    parser.add_argument('--plot-dir', dest='plot_dir', default=None,
                        help='Directory for diagnostic plots (defaults to output_dir/plots)')

    args = parser.parse_args()

    # Load configuration
    utils.log(f"Case config file is {args.case_config_path}")
    try:
        config = utils.import_all_from(args.case_config_path)
        utils.log("Successfully imported case configuration parameters!")
    except Exception:
        utils.fatal("Unable to import case configuration parameters! Aborting!")

    # Extract configuration parameters
    deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
    PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
    PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
    PIPS_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)
    conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', None)
    requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.0)

    if not all([deployment_names, PIPS_dir, PIPS_names]):
        utils.fatal("Required configuration parameters missing from config file!")

    output_dir = args.output_dir if args.output_dir else PIPS_dir

    # Set up plot directory if diagnostic plots are requested
    plot_dir = None
    if args.diagnostic_plots:
        plot_dir = args.plot_dir if args.plot_dir else os.path.join(output_dir, 'plots')
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
            utils.log(f"Created plot directory: {plot_dir}")
        utils.log(f"Plot directory: {plot_dir}")

    utils.log(f"Output directory: {output_dir}")
    utils.log(f"Processing {len(PIPS_names)} PIPS instruments")

    # Create output directory if needed
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Process each PIPS
    for index, (PIPS_name, deployment_name) in enumerate(zip(PIPS_names, deployment_names)):
        utils.log(f"\n{'='*60}")
        utils.log(f"Processing PIPS {index+1}/{len(PIPS_names)}: {PIPS_name}")
        utils.log(f"Deployment: {deployment_name}")

        try:
            process_single_pips(
                PIPS_name=PIPS_name,
                deployment_name=deployment_name,
                PIPS_dir=PIPS_dir,
                realtime_dir=args.realtime_dir,
                output_dir=output_dir,
                PIPS_filenames_nc=PIPS_filenames_nc,
                conv_filenames_nc=conv_filenames_nc,
                index=index,
                output_tag=args.output_tag,
                requested_interval=requested_interval,
                diagnostic_plots=args.diagnostic_plots,
                plot_dir=plot_dir
            )
            utils.log(f"Successfully processed {PIPS_name}")

        except Exception as e:
            utils.log(f"Error processing {PIPS_name}: {str(e)}")
            continue

    utils.log(f"\n{'='*60}")
    utils.log("Merge processing completed!")


def process_single_pips(PIPS_name, deployment_name, PIPS_dir, realtime_dir, output_dir,
                       PIPS_filenames_nc, conv_filenames_nc, index, output_tag, requested_interval,
                       diagnostic_plots=False, plot_dir=None):
    """
    Process a single PIPS instrument for merging real-time and card-based data.
    """

    # Find and load conventional file first to get time range
    utils.log(f"Loading conventional file from config to get time range...")

    # Get conventional filename from config
    if conv_filenames_nc and index < len(conv_filenames_nc):
        onesec_filename = conv_filenames_nc[index]
        onesec_path = os.path.join(PIPS_dir, onesec_filename)
    else:
        # Fallback to pattern matching if config doesn't have conv filenames
        onesec_pattern = os.path.join(PIPS_dir, f'conventional_raw_*{PIPS_name}*.nc')
        onesec_files = glob.glob(onesec_pattern)
        if not onesec_files:
            utils.log(f"No conventional files found for {PIPS_name}, skipping...")
            return
        onesec_path = onesec_files[0]

    # Check if conventional file exists and load it to get time range
    if not os.path.exists(onesec_path):
        utils.log(f"Conventional file not found: {onesec_path}, skipping {PIPS_name}...")
        return

    utils.log(f"Loading card onesec file: {onesec_path}")
    onesec_card_ds = xr.open_dataset(onesec_path, decode_timedelta=False)
    # Load into memory and close file to prevent permission issues when saving
    onesec_card_ds = onesec_card_ds.load()
    onesec_card_ds.close()

    # Read time range from conventional netCDF attributes
    try:
        starting_time_str = onesec_card_ds.attrs['starting_time']
        ending_time_str = onesec_card_ds.attrs['ending_time']

        # Parse time strings to datetime objects
        starttime_dt = datetime.strptime(starting_time_str, '%Y%m%d%H%M%S')
        endtime_dt = datetime.strptime(ending_time_str, '%Y%m%d%H%M%S')

        utils.log(f"Time range from conventional netCDF: {starttime_dt} to {endtime_dt}")

    except (KeyError, ValueError) as e:
        utils.log(f"Error reading time range from conventional netCDF attributes: {e}")
        utils.log(f"Skipping {PIPS_name}...")
        return

    # Now load parsivel file
    if PIPS_filenames_nc and index < len(PIPS_filenames_nc):
        parsivel_combined_filename = PIPS_filenames_nc[index]
        parsivel_combined_path = os.path.join(PIPS_dir, parsivel_combined_filename)
    else:
        # Fallback to pattern matching if config doesn't have filenames
        parsivel_combined_pattern = os.path.join(PIPS_dir, f'parsivel_combined_*{PIPS_name}*.nc')
        parsivel_combined_files = glob.glob(parsivel_combined_pattern)
        if not parsivel_combined_files:
            utils.log(f"No parsivel combined files found for {PIPS_name}, skipping...")
            return
        parsivel_combined_path = parsivel_combined_files[0]

    # Check if parsivel file exists
    if not os.path.exists(parsivel_combined_path):
        utils.log(f"Parsivel file not found: {parsivel_combined_path}, skipping {PIPS_name}...")
        return

    utils.log(f"Loading card parsivel file: {parsivel_combined_path}")
    parsivel_combined_card_ds = xr.open_dataset(parsivel_combined_path, decode_timedelta=False)
    # Load into memory and close file to prevent permission issues when saving
    parsivel_combined_card_ds = parsivel_combined_card_ds.load()
    parsivel_combined_card_ds.close()

    utils.log(f"Finding real-time files for {PIPS_name}...")
    file_path_list_onePIPS = glob.glob(os.path.join(realtime_dir, f'*{PIPS_name}*nc'))
    file_list_onePIPS = [os.path.basename(file_path) for file_path in file_path_list_onePIPS]

    # Get different real-time file types
    file_list_onesec = [f for f in file_list_onePIPS if 'onesec' in f and 'current' not in f]
    file_list_ND = [f for f in file_list_onePIPS if 'ND' in f and 'current' not in f]
    file_list_spectrum = [f for f in file_list_onePIPS if 'spectrum' in f and 'current' not in f]
    file_list_telegram = [f for f in file_list_onePIPS if 'telegram' in f and 'current' not in f]

    # Filter by time range
    file_list_onesec = get_files(file_list_onesec, starttime_dt, endtime_dt, 'onesec')
    file_list_ND = get_files(file_list_ND, starttime_dt, endtime_dt, 'ND')
    file_list_spectrum = get_files(file_list_spectrum, starttime_dt, endtime_dt, 'spectrum')
    file_list_telegram = get_files(file_list_telegram, starttime_dt, endtime_dt, 'telegram')

    # Create full file paths for real-time data
    file_path_list_onesec = [os.path.join(realtime_dir, f) for f in file_list_onesec]
    file_path_list_ND = [os.path.join(realtime_dir, f) for f in file_list_ND]
    file_path_list_spectrum = [os.path.join(realtime_dir, f) for f in file_list_spectrum]
    file_path_list_telegram = [os.path.join(realtime_dir, f) for f in file_list_telegram]

    utils.log(f"Found {len(file_path_list_onesec)} real-time onesec files")
    utils.log(f"Found {len(file_path_list_ND)} real-time ND files")
    utils.log(f"Found {len(file_path_list_spectrum)} real-time spectrum files")
    utils.log(f"Found {len(file_path_list_telegram)} real-time telegram files")

    # Load real-time datasets
    utils.log("Loading real-time datasets...")
    onesec_rt_ds = xr.open_mfdataset(file_path_list_onesec, combine='nested', concat_dim='logger_datetime',
                                     decode_timedelta=False) if file_path_list_onesec else None
    if onesec_rt_ds is not None:
        onesec_rt_ds = onesec_rt_ds.rename({'logger_datetime': 'time'})
        onesec_rt_ds = onesec_rt_ds.drop_duplicates('time')
        onesec_rt_ds = onesec_rt_ds.load()

    ND_rt_ds = xr.open_mfdataset(file_path_list_ND, combine='nested', concat_dim='time',
                                 decode_timedelta=False) if file_path_list_ND else None
    if ND_rt_ds is not None:
        ND_rt_ds = ND_rt_ds.drop_duplicates('time')
        ND_rt_ds = ND_rt_ds.rename_dims({'diameter': 'diameter_bin'})
        ND_rt_ds = ND_rt_ds.load()

    spectrum_rt_ds = xr.open_mfdataset(file_path_list_spectrum, combine='nested', concat_dim='time',
                                       decode_timedelta=False) if file_path_list_spectrum else None
    if spectrum_rt_ds is not None:
        spectrum_rt_ds = spectrum_rt_ds.drop_duplicates('time')
        spectrum_rt_ds = spectrum_rt_ds.rename_dims({'diameter': 'diameter_bin', 'velocity': 'fallspeed_bin'})
        spectrum_rt_ds = spectrum_rt_ds.rename({'velocity': 'fallspeed'})
        spectrum_rt_ds = spectrum_rt_ds.load()

    telegram_rt_ds = xr.open_mfdataset(file_path_list_telegram, combine='nested', concat_dim='index',
                                       decode_timedelta=False) if file_path_list_telegram else None
    if telegram_rt_ds is not None:
        telegram_rt_ds = telegram_rt_ds.rename({'index': 'time'})
        telegram_rt_ds = telegram_rt_ds.drop_duplicates('time')

    # Conventional file already loaded earlier to get time range

    # Correct times from real-time files to match GPS times (already done for card files)
    if onesec_rt_ds is not None:
        utils.log("Correcting real-time file times using GPS data...")

        # Get first good GPS time in file
        first_good_GPS_time = onesec_rt_ds.where(onesec_rt_ds['GPS_status'].compute() == 'A', drop=True).isel(time=0)

        logger_datetime = first_good_GPS_time['time'].values
        logger_datetime = pd.to_datetime(logger_datetime).to_pydatetime()

        GPS_date = str(first_good_GPS_time['GPS_date'].values)
        GPS_time = str(first_good_GPS_time['GPS_time'].values)

        # Construct datetime object from GPS info
        gyear = int('20' + GPS_date[4:])
        gmonth = int(GPS_date[2:4])
        gday = int(GPS_date[:2])
        ghour = int(GPS_time[:2])
        gmin = int(GPS_time[2:4])
        gsec = int(GPS_time[4:6])

        GPS_datetime = datetime(gyear, gmonth, gday, ghour, gmin, gsec)
        GPS_offset = GPS_datetime - logger_datetime

        utils.log(f'GPS time: {GPS_datetime.strftime("%Y-%m-%d %H:%M:%S")}, Logger time: {logger_datetime.strftime("%Y-%m-%d %H:%M:%S")}')
        utils.log(f'GPS Offset: {str(GPS_offset)}')

        # Apply GPS correction to all real-time datasets
        old_times = pd.to_datetime(onesec_rt_ds['time']).to_pydatetime()
        new_times = old_times + GPS_offset
        onesec_rt_ds = onesec_rt_ds.assign_coords({'time': new_times})

        if ND_rt_ds is not None:
            old_ND_times = pd.to_datetime(ND_rt_ds['time']).to_pydatetime()
            new_times = old_ND_times + GPS_offset
            ND_rt_ds = ND_rt_ds.assign_coords({'time': new_times})

        if spectrum_rt_ds is not None:
            old_spectrum_times = pd.to_datetime(spectrum_rt_ds['time']).to_pydatetime()
            new_times = old_spectrum_times + GPS_offset
            spectrum_rt_ds = spectrum_rt_ds.assign_coords({'time': new_times})

        if telegram_rt_ds is not None:
            old_telegram_times = pd.to_datetime(telegram_rt_ds['time']).to_pydatetime()
            new_times = old_telegram_times + GPS_offset
            telegram_rt_ds = telegram_rt_ds.assign_coords({'time': new_times})

    # Correct real-time dewpoint if real-time data exists
    if onesec_rt_ds is not None:
        utils.log("Correcting real-time dewpoint calculation...")
        onesec_rt_ds = correct_realtime_dewpoint(onesec_rt_ds)

    # Create time ranges for merging
    utils.log("Creating merged time ranges...")

    # Use card data times as the authoritative time range (as in original notebook)
    onesec_card_starttime = onesec_card_ds['time'][0].values
    onesec_card_endtime = onesec_card_ds['time'][-1].values

    all_onesec_times = xr.date_range(
        start=onesec_card_starttime,
        end=onesec_card_endtime,
        freq='1S'
    )

    all_tensec_times = xr.date_range(
        start=parsivel_combined_card_ds['time'][0].values,
        end=parsivel_combined_card_ds['time'][-1].values,
        freq=f'{int(requested_interval)}S'
    )

    # Merge onesec data (notebook approach)
    utils.log("Merging one-second data...")
    if onesec_rt_ds is not None:
        # Direct merge using combine_first (real-time priority)
        onesec_merged_ds = onesec_rt_ds.combine_first(onesec_card_ds)
        # Reindex to full time range
        onesec_merged_full_ds = onesec_merged_ds.reindex({'time': all_onesec_times})
    else:
        # Only card data available
        onesec_merged_full_ds = onesec_card_ds.reindex({'time': all_onesec_times})

    # Recompute thermodynamic parameters for merged data (as in notebook)
    utils.log("Recomputing thermodynamic parameters...")
    onesec_merged_full_ds = pips.calc_thermo(onesec_merged_full_ds)

    # Merge parsivel data (notebook approach)
    utils.log("Merging parsivel data...")
    parsivel_combined_card_ds_full = parsivel_combined_card_ds.reindex({'time': all_tensec_times})
    parsivel_combined_merged_full_ds = parsivel_combined_card_ds_full.copy()

    # Merge ND data
    if ND_rt_ds is not None and 'ND' in parsivel_combined_card_ds.data_vars:
        utils.log("Merging ND data...")
        ND_merged_da = ND_rt_ds['ND'].combine_first(parsivel_combined_card_ds['ND'])
        ND_merged_da_full = ND_merged_da.reindex({'time': all_tensec_times})
        parsivel_combined_merged_full_ds['ND'] = ND_merged_da_full

    # Merge spectrum data
    if spectrum_rt_ds is not None and 'VD_matrix' in parsivel_combined_card_ds.data_vars:
        utils.log("Merging spectrum data...")
        spectrum_merged_da = spectrum_rt_ds['spectrum'].combine_first(parsivel_combined_card_ds['VD_matrix'])
        spectrum_merged_da_full = spectrum_merged_da.reindex({'time': all_tensec_times})
        parsivel_combined_merged_full_ds['VD_matrix'] = spectrum_merged_da_full

    # Merge telegram data
    # Merge telegram data (notebook approach)
    if telegram_rt_ds is not None:
        utils.log("Merging telegram data...")

        telegram_mapping = {
            'rain rate (mm per hr)': 'precipintensity',
            'rain accumulation (mm)': 'precipaccum',
            'radar reflectivity (dBZ)': 'parsivel_dBZ',
            'sample interval': 'sample_interval',
            'signal amplitude': 'signal_amplitude',
            'particle count': 'pcount',
            'sensor temp': 'sensor_temp',
            'power supply voltage': 'pvoltage'
        }

        for rt_varname, card_varname in telegram_mapping.items():
            if rt_varname in telegram_rt_ds.data_vars and card_varname in parsivel_combined_card_ds.data_vars:
                merged_da = telegram_rt_ds[rt_varname].combine_first(parsivel_combined_card_ds[card_varname])
                merged_da_full = merged_da.reindex({'time': all_tensec_times})
                parsivel_combined_merged_full_ds[card_varname] = merged_da_full

    # Resample merged onesec data to requested interval
    utils.log(f"Resampling merged conventional data to {requested_interval}-second intervals...")
    PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_merged_full_ds['VD_matrix'])
    sec_offset = PSD_datetimes[0].second
    resample_interval = int(requested_interval)
    conv_resampled_ds = pips.resample_conv_da('PIPS', resample_interval, sec_offset,
                                              onesec_merged_full_ds, gusts=True, gustintvstr='3S')

    # Update parsivel dataset with resampled conventional data
    for varname in conv_resampled_ds.data_vars:
        if varname in parsivel_combined_merged_full_ds.data_vars:
            parsivel_combined_merged_full_ds[varname] = conv_resampled_ds[varname]

    # Copy attributes (notebook approach)
    utils.log("Copying dataset attributes...")

    # Global attributes
    onesec_merged_full_ds.attrs = onesec_card_ds.attrs
    parsivel_combined_merged_full_ds.attrs = parsivel_combined_card_ds.attrs

    # Variable attributes
    for varname in onesec_merged_full_ds.data_vars:
        if varname in onesec_card_ds.data_vars:
            onesec_merged_full_ds[varname].attrs = onesec_card_ds[varname].attrs

    for varname in parsivel_combined_merged_full_ds.data_vars:
        if varname in parsivel_combined_card_ds.data_vars:
            parsivel_combined_merged_full_ds[varname].attrs = parsivel_combined_card_ds[varname].attrs

    # Create diagnostic plots if requested
    if diagnostic_plots and plot_dir is not None:
        create_diagnostic_plots(
            all_onesec_times, onesec_rt_ds, onesec_card_ds,
            all_tensec_times, parsivel_combined_merged_full_ds,
            PIPS_name, deployment_name, plot_dir, ND_rt_ds,
            parsivel_combined_card_ds
        )

    # Save merged datasets
    utils.log("Saving merged datasets...")

    # Output filenames with configurable tag
    tag_suffix = f'_{output_tag}' if output_tag else ''
    parsivel_combined_output_filename = f'parsivel_combined_{deployment_name}_{PIPS_name}_{int(requested_interval)}s{tag_suffix}.nc'
    parsivel_combined_output_path = os.path.join(output_dir, parsivel_combined_output_filename)

    onesec_output_filename = f'conventional_raw_{deployment_name}_{PIPS_name}{tag_suffix}.nc'
    onesec_output_path = os.path.join(output_dir, onesec_output_filename)

    utils.log(f"Saving {onesec_output_path}")
    try:
        onesec_merged_full_ds.to_netcdf(onesec_output_path)
    except PermissionError as e:
        utils.fatal(f"Permission denied when saving {onesec_output_path}: {e}. "
                   f"File may still be open or locked. Aborting.")
    except Exception as e:
        utils.fatal(f"Error saving {onesec_output_path}: {e}. Aborting.")

    utils.log(f"Saving {parsivel_combined_output_path}")
    try:
        parsivel_combined_merged_full_ds.to_netcdf(parsivel_combined_output_path)
    except PermissionError as e:
        utils.fatal(f"Permission denied when saving {parsivel_combined_output_path}: {e}. "
                   f"File may still be open or locked. Aborting.")
    except Exception as e:
        utils.fatal(f"Error saving {parsivel_combined_output_path}: {e}. Aborting.")

    utils.log(f"Merge completed successfully for {PIPS_name}!")


if __name__ == '__main__':
    main()