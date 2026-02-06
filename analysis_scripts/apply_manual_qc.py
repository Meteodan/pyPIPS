#!/usr/bin/env python
"""
Apply manual QC decisions to PIPS datasets.

Reads QC time ranges from JSON config and sets specified variables to NaN
in the identified time periods. Updates both conventional and parsivel datasets.

Usage:
    python apply_manual_qc.py <case_config.py> <manual_qc_decisions.json> [options]

Example:
    python apply_manual_qc.py configs/ICECHIP_IOP1_2025_10s.py \\
                              configs/manual_qc_decisions.json \\
                              --output-tag _manual_qc
"""
import os
import argparse
import json
import numpy as np
import pandas as pd
import xarray as xr
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import pyPIPS.utils as utils


def apply_manual_qc_to_dataset(ds, qc_entries, verbose=True):
    """
    Apply manual QC entries to a dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to apply QC to
    qc_entries : list of dict
        List of QC entries with 'variable', 'start_time', 'end_time', 'reason'
    verbose : bool
        Print QC actions

    Returns
    -------
    ds : xarray.Dataset
        Dataset with QC applied
    n_flagged_total : int
        Total number of points flagged across all variables
    """
    n_flagged_total = 0

    for entry in qc_entries:
        variable = entry['variable']
        start_time = pd.to_datetime(entry['start_time'])
        if start_time.tz is not None:
            start_time = start_time.tz_localize(None)
        end_time = pd.to_datetime(entry['end_time'])
        if end_time.tz is not None:
            end_time = end_time.tz_localize(None)
        reason = entry.get('reason', 'Manual QC')

        if variable not in ds:
            if verbose:
                utils.log(f"  Warning: Variable '{variable}' not in dataset, skipping")
            continue

        # Convert dataset times to timezone-naive pandas datetime for comparison
        ds_times = pd.to_datetime(ds.time.values)
        if hasattr(ds_times, 'tz') and ds_times.tz is not None:
            ds_times = ds_times.tz_localize(None)

        # Find time range to flag using pandas datetime comparison
        time_mask = (ds_times >= start_time) & (ds_times <= end_time)
        n_flagged = int(time_mask.sum())

        if n_flagged > 0:
            # Get the time indices where mask is True and set to NaN
            time_indices = np.where(time_mask)[0]

            # Create modified data array
            current_data = ds[variable].values.copy()
            current_data[time_indices] = np.nan

            # Rebuild the DataArray with new values
            new_var = xr.DataArray(
                current_data,
                coords=ds[variable].coords,
                dims=ds[variable].dims,
                attrs=ds[variable].attrs.copy()
            )
            new_var.encoding = {}  # Clear encoding to avoid fill value conflicts

            # Replace in dataset
            ds[variable] = new_var

            # Update attributes to track manual QC
            if 'manual_qc_applied' not in ds[variable].attrs:
                ds[variable].attrs['manual_qc_applied'] = []

            qc_info = f"{start_time.isoformat()}_{end_time.isoformat()}: {reason}"

            # Handle both list and string attribute types
            if isinstance(ds[variable].attrs.get('manual_qc_applied'), list):
                ds[variable].attrs['manual_qc_applied'].append(qc_info)
            else:
                ds[variable].attrs['manual_qc_applied'] = [qc_info]

            # Also set a simple flag for easy checking (as integer for NetCDF compatibility)
            ds[variable].attrs['has_manual_qc'] = 1

            n_flagged_total += n_flagged

            if verbose:
                utils.log(f"  {variable}: Flagged {n_flagged} points from {start_time} to {end_time}")
                if reason:
                    utils.log(f"    Reason: {reason}")
        else:
            if verbose:
                utils.log(f"  Warning: No data found in specified time range for {variable}")

    return ds, n_flagged_total


def align_times_with_parsivel(parsivel_ds, requested_start, requested_end, verbose=True):
    """
    Align requested trim times to parsivel dataset's time grid.
    Start time rounds UP to next parsivel time, end time rounds DOWN to previous parsivel time.

    Parameters
    ----------
    parsivel_ds : xarray.Dataset
        Parsivel dataset with coarser time resolution
    requested_start : pd.Timestamp or None
        Requested start time (None to keep dataset start)
    requested_end : pd.Timestamp or None
        Requested end time (None to keep dataset end)
    verbose : bool
        Print alignment information

    Returns
    -------
    aligned_start : pd.Timestamp
        Start time aligned to parsivel grid
    aligned_end : pd.Timestamp
        End time aligned to parsivel grid
    """
    parsivel_times = pd.to_datetime(parsivel_ds.time.values)
    # Ensure parsivel times are timezone-naive
    if hasattr(parsivel_times, 'tz') and parsivel_times.tz is not None:
        parsivel_times = parsivel_times.tz_localize(None)

    # Align start time (round UP to next parsivel time)
    if requested_start is None:
        aligned_start = parsivel_times[0]
    else:
        # Find first parsivel time >= requested start
        valid_starts = parsivel_times[parsivel_times >= requested_start]
        if len(valid_starts) > 0:
            aligned_start = valid_starts[0]
        else:
            # Requested start is after all parsivel times
            aligned_start = parsivel_times[-1]
            if verbose:
                utils.log(f"    Warning: Requested start time after parsivel end, using last parsivel time")

    # Align end time (round DOWN to previous parsivel time)
    if requested_end is None:
        aligned_end = parsivel_times[-1]
    else:
        # Find last parsivel time <= requested end
        valid_ends = parsivel_times[parsivel_times <= requested_end]
        if len(valid_ends) > 0:
            aligned_end = valid_ends[-1]
        else:
            # Requested end is before all parsivel times
            aligned_end = parsivel_times[0]
            if verbose:
                utils.log(f"    Warning: Requested end time before parsivel start, using first parsivel time")

    if verbose and (requested_start is not None or requested_end is not None):
        if requested_start is not None and aligned_start != requested_start:
            utils.log(f"    Aligned start: {requested_start} → {aligned_start} (next parsivel time)")
        if requested_end is not None and aligned_end != requested_end:
            utils.log(f"    Aligned end:   {requested_end} → {aligned_end} (previous parsivel time)")

    return aligned_start, aligned_end


def apply_manual_trim_to_dataset(ds, trim_entry, aligned_times=None, verbose=True):
    """
    Apply manual trimming to a dataset - trim to new time bounds.
    Similar to trim_fixed_deployment_periods in apply_QC.py but with user-specified times.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to trim
    trim_entry : dict
        Trim entry with 'new_start_time' and/or 'new_end_time', 'reason'
    aligned_times : tuple of pd.Timestamp, optional
        Pre-computed (start, end) times aligned to parsivel grid.
        If provided, these override times from trim_entry.
    verbose : bool
        Print trimming actions

    Returns
    -------
    ds : xarray.Dataset
        Trimmed dataset with updated global attributes
    n_trimmed : int
        Number of time points removed
    """
    original_length = len(ds.time)
    original_start = pd.to_datetime(ds.time.values[0])
    if hasattr(original_start, 'tz') and original_start.tz is not None:
        original_start = original_start.tz_localize(None)
    original_end = pd.to_datetime(ds.time.values[-1])
    if hasattr(original_end, 'tz') and original_end.tz is not None:
        original_end = original_end.tz_localize(None)

    reason = trim_entry.get('reason', 'Manual trimming')

    # Use aligned times if provided, otherwise use requested times from trim_entry
    if aligned_times is not None:
        new_start, new_end = aligned_times
    else:
        new_start_time = trim_entry.get('new_start_time')
        new_end_time = trim_entry.get('new_end_time')

        # Use original times as defaults
        if new_start_time is None:
            new_start = original_start
        else:
            new_start = pd.to_datetime(new_start_time)
            if new_start.tz is not None:
                new_start = new_start.tz_localize(None)

        if new_end_time is None:
            new_end = original_end
        else:
            new_end = pd.to_datetime(new_end_time)
            if new_end.tz is not None:
                new_end = new_end.tz_localize(None)

    if verbose:
        utils.log(f"  Trimming dataset:")
        utils.log(f"    Original: {original_start} to {original_end} ({original_length} points)")
        utils.log(f"    New:      {new_start} to {new_end}")
        if reason:
            utils.log(f"    Reason: {reason}")

    # Trim dataset to new time range
    ds = ds.sel(time=slice(new_start, new_end))
    new_length = len(ds.time)
    n_trimmed = original_length - new_length

    if n_trimmed > 0:
        # Update global attributes with new time bounds
        ds.attrs['starting_time'] = new_start.strftime('%Y%m%d%H%M%S')
        ds.attrs['ending_time'] = new_end.strftime('%Y%m%d%H%M%S')

        # Update time encoding to reference new start time
        ds['time'].encoding['units'] = f"seconds since {new_start.strftime('%Y-%m-%d %H:%M:%S')}"
        ds['time'].encoding['calendar'] = 'proleptic_gregorian'

        # Add manual trim provenance to global attributes
        if 'manual_trim_applied' not in ds.attrs:
            ds.attrs['manual_trim_applied'] = f"{new_start.isoformat()}_{new_end.isoformat()}: {reason}"
        else:
            ds.attrs['manual_trim_applied'] += f"; {new_start.isoformat()}_{new_end.isoformat()}: {reason}"

        if verbose:
            utils.log(f"    ✓ Trimmed {n_trimmed} points ({n_trimmed} seconds)")
            utils.log(f"    New length: {new_length} points")

    else:
        if verbose:
            utils.log(f"    Warning: No points trimmed (times may be outside dataset range)")

    return ds, n_trimmed


def main():
    parser = argparse.ArgumentParser(
        description='Apply manual QC decisions to PIPS datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Apply manual QC with default output tag
  python apply_manual_qc.py configs/ICECHIP_IOP1_2025_10s.py \\
                            configs/manual_qc_decisions.json

  # Apply with custom output tag
  python apply_manual_qc.py configs/ICECHIP_IOP1_2025_10s.py \\
                            configs/manual_qc_decisions.json \\
                            --output-tag _mqc

  # Overwrite original files (use with caution!)
  python apply_manual_qc.py configs/ICECHIP_IOP1_2025_10s.py \\
                            configs/manual_qc_decisions.json \\
                            --output-tag ""
        """
    )
    parser.add_argument('case_config_path',
                        help='Path to case configuration file')
    parser.add_argument('manual_qc_file',
                        help='Path to manual QC decisions JSON file')
    parser.add_argument('--output-tag', dest='output_tag', default='_manual_qc',
                        help='Tag to add to output filenames (default: _manual_qc). '
                             'Use empty string "" to overwrite original files.')
    parser.add_argument('--conventional-only', dest='conv_only', action='store_true',
                        help='Only process conventional datasets (skip parsivel files)')
    parser.add_argument('--parsivel-only', dest='parsivel_only', action='store_true',
                        help='Only process parsivel datasets (skip conventional files)')

    args = parser.parse_args()

    # Load case configuration
    utils.log(f"Loading case configuration: {args.case_config_path}")
    try:
        config = utils.import_all_from(args.case_config_path)
        utils.log("✓ Successfully imported case configuration")
    except Exception as e:
        utils.fatal(f"Unable to import case configuration: {e}")

    PIPS_dir = config.PIPS_IO_dict['PIPS_dir']
    deployment_names = config.PIPS_IO_dict['deployment_names']
    PIPS_names = config.PIPS_IO_dict['PIPS_names']
    conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', [])
    PIPS_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', [])

    # Load manual QC decisions
    utils.log(f"Loading manual QC decisions: {args.manual_qc_file}")
    try:
        with open(args.manual_qc_file, 'r') as f:
            all_data = json.load(f)

        # Separate QC entries from trim entries
        manual_trim_dict = all_data.pop('_manual_trim', {})
        manual_qc_dict = all_data  # Everything else is QC entries

        utils.log(f"✓ Loaded {len(manual_qc_dict)} dataset QC configurations")
        total_qc_entries = sum(len(v) for v in manual_qc_dict.values() if isinstance(v, list))
        total_trim_entries = len(manual_trim_dict)
        utils.log(f"  Total QC entries: {total_qc_entries}")
        utils.log(f"  Total trim entries: {total_trim_entries}")
    except FileNotFoundError:
        utils.fatal(f"Manual QC file not found: {args.manual_qc_file}")
    except json.JSONDecodeError as e:
        utils.fatal(f"Invalid JSON in manual QC file: {e}")

    if total_qc_entries == 0 and total_trim_entries == 0:
        utils.log("Warning: No manual QC or trim entries found in file. Nothing to do.")
        return

    # Track statistics
    datasets_processed = 0
    total_points_flagged = 0
    total_points_trimmed = 0

    # Process each PIPS/deployment pair
    # deployment_names is parallel to PIPS_names (each PIPS has its deployment label)
    for idx, (pips_name, deployment_name) in enumerate(zip(PIPS_names, deployment_names)):
        key = f"{pips_name}_{deployment_name}"

        has_qc_entries = key in manual_qc_dict and len(manual_qc_dict[key]) > 0
        has_trim_entry = key in manual_trim_dict

        if not has_qc_entries and not has_trim_entry:
            continue

        utils.log(f"\n{'='*70}")
        utils.log(f"Processing {key}...")
        if has_qc_entries:
            utils.log(f"  {len(manual_qc_dict[key])} manual QC entries to apply")
        if has_trim_entry:
            utils.log(f"  Manual trimming to apply")

        # Get trim_entry for easier access
        trim_entry = manual_trim_dict.get(key) if has_trim_entry else None

        # If trimming, compute aligned times using parsivel dataset first
        aligned_trim_times = None
        if has_trim_entry:
            # Get parsivel filename from config
            parsivel_file = None
            if idx < len(PIPS_filenames_nc) and PIPS_filenames_nc[idx]:
                parsivel_file = os.path.join(PIPS_dir, PIPS_filenames_nc[idx])
                if not os.path.exists(parsivel_file):
                    parsivel_file = None

            if parsivel_file:
                utils.log(f"  Aligning trim times with parsivel grid: {os.path.basename(parsivel_file)}")
                try:
                    parsivel_ds_temp = xr.open_dataset(parsivel_file)

                    # Get requested times
                    requested_start = trim_entry.get('new_start_time')
                    requested_end = trim_entry.get('new_end_time')

                    if requested_start is not None:
                        requested_start = pd.to_datetime(requested_start)
                        if requested_start.tz is not None:
                            requested_start = requested_start.tz_localize(None)
                    if requested_end is not None:
                        requested_end = pd.to_datetime(requested_end)
                        if requested_end.tz is not None:
                            requested_end = requested_end.tz_localize(None)

                    # Align to parsivel grid
                    aligned_trim_times = align_times_with_parsivel(
                        parsivel_ds_temp, requested_start, requested_end, verbose=True)

                    parsivel_ds_temp.close()

                    # Update trim_entry with aligned times and save back to file
                    if aligned_trim_times is not None:
                        aligned_start, aligned_end = aligned_trim_times
                        trim_entry['new_start_time'] = aligned_start.isoformat()
                        trim_entry['new_end_time'] = aligned_end.isoformat()

                        # Add note about alignment
                        if 'aligned_to_parsivel' not in trim_entry:
                            original_reason = trim_entry.get('reason', '')
                            if original_reason and not original_reason.endswith(')'):
                                trim_entry['reason'] = f"{original_reason} (aligned to parsivel grid)"
                            trim_entry['aligned_to_parsivel'] = True

                        # Save updated manual_qc/trim dictionaries back to file
                        try:
                            output_data = manual_qc_dict.copy()
                            output_data['_manual_trim'] = manual_trim_dict
                            with open(args.manual_qc_file, 'w') as f:
                                json.dump(output_data, f, indent=2)
                            utils.log(f"    ✓ Updated trim times saved to {args.manual_qc_file}")
                        except Exception as save_err:
                            utils.log(f"    Warning: Could not save updated times: {save_err}")

                except Exception as e:
                    utils.log(f"  Warning: Could not align times with parsivel: {e}")
                    utils.log(f"  Proceeding with requested times...")
            else:
                utils.log(f"  Note: Parsivel file not found, using requested times without alignment")

        # Process conventional dataset
        if not args.parsivel_only:
            # Get conventional filename from config
            conv_file = None
            if idx < len(conv_filenames_nc) and conv_filenames_nc[idx]:
                conv_file = os.path.join(PIPS_dir, conv_filenames_nc[idx])

            if conv_file and os.path.exists(conv_file):
                utils.log(f"  Loading conventional dataset: {os.path.basename(conv_file)}")
                try:
                    conv_ds = xr.open_dataset(conv_file)
                    conv_ds.load()  # Load into memory

                    # Apply manual trimming first (if specified)
                    n_trimmed = 0
                    if has_trim_entry:
                        conv_ds, n_trimmed = apply_manual_trim_to_dataset(
                            conv_ds, trim_entry, aligned_times=aligned_trim_times, verbose=True)
                        if n_trimmed > 0:
                            total_points_trimmed += n_trimmed

                    # Apply manual QC (if specified)
                    n_flagged = 0
                    if has_qc_entries:
                        conv_ds, n_flagged = apply_manual_qc_to_dataset(
                            conv_ds, manual_qc_dict[key], verbose=True)

                    if n_trimmed > 0 or n_flagged > 0:
                        # Determine output filename
                        if args.output_tag:
                            output_file = conv_file.replace('.nc', f'{args.output_tag}.nc')
                        else:
                            output_file = conv_file

                        utils.log(f"  Saving conventional data: {os.path.basename(output_file)}")
                        conv_ds.to_netcdf(output_file)

                        summary_parts = []
                        if n_trimmed > 0:
                            summary_parts.append(f"{n_trimmed} points trimmed")
                        if n_flagged > 0:
                            summary_parts.append(f"{n_flagged} points flagged")
                            total_points_flagged += n_flagged

                        utils.log(f"  ✓ Conventional dataset: {', '.join(summary_parts)}")
                        datasets_processed += 1
                    else:
                        utils.log(f"  No changes to conventional dataset")

                    conv_ds.close()

                except Exception as e:
                    utils.log(f"  Error processing conventional file: {e}")
            else:
                utils.log(f"  Warning: Conventional file not found")

        # Process parsivel dataset
        if not args.conv_only:
            # Get parsivel filename from config
            parsivel_file = None
            if idx < len(PIPS_filenames_nc) and PIPS_filenames_nc[idx]:
                parsivel_file = os.path.join(PIPS_dir, PIPS_filenames_nc[idx])

            if parsivel_file and os.path.exists(parsivel_file):
                utils.log(f"  Loading parsivel dataset: {os.path.basename(parsivel_file)}")
                try:
                    parsivel_ds = xr.open_dataset(parsivel_file)
                    parsivel_ds.load()  # Load into memory

                    # Apply manual trimming first (if specified)
                    n_trimmed = 0
                    if has_trim_entry:
                        parsivel_ds, n_trimmed = apply_manual_trim_to_dataset(
                            parsivel_ds, trim_entry, aligned_times=aligned_trim_times, verbose=True)
                        if n_trimmed > 0:
                            total_points_trimmed += n_trimmed

                    # Apply manual QC (if specified)
                    n_flagged = 0
                    if has_qc_entries:
                        parsivel_ds, n_flagged = apply_manual_qc_to_dataset(
                            parsivel_ds, manual_qc_dict[key], verbose=True)

                    if n_trimmed > 0 or n_flagged > 0:
                        # Determine output filename
                        if args.output_tag:
                            output_file = parsivel_file.replace('.nc', f'{args.output_tag}.nc')
                        else:
                            output_file = parsivel_file

                        utils.log(f"  Saving parsivel data: {os.path.basename(output_file)}")
                        parsivel_ds.to_netcdf(output_file)

                        summary_parts = []
                        if n_trimmed > 0:
                            summary_parts.append(f"{n_trimmed} points trimmed")
                        if n_flagged > 0:
                            summary_parts.append(f"{n_flagged} points flagged")
                            total_points_flagged += n_flagged

                        utils.log(f"  ✓ Parsivel dataset: {', '.join(summary_parts)}")
                        datasets_processed += 1
                    else:
                        utils.log(f"  No changes to parsivel dataset")

                    parsivel_ds.close()

                except Exception as e:
                    utils.log(f"  Error processing parsivel file: {e}")
            else:
                utils.log(f"  Warning: Parsivel file not found")

    # Summary
    utils.log(f"\n{'='*70}")
    utils.log("Manual QC application complete!")
    utils.log(f"  Datasets processed: {datasets_processed}")
    utils.log(f"  Total points trimmed: {total_points_trimmed}")
    utils.log(f"  Total points flagged: {total_points_flagged}")
    utils.log(f"{'='*70}")


if __name__ == '__main__':
    main()
