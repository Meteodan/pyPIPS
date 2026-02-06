# apply_QC.py
#
# This script applies quality control to the PIPS DSDs (netCDF version)
import os
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import pyPIPS.parsivel_params as pp
import pyPIPS.parsivel_qc as pqc
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.DSDlib as dsd
import pyPIPS.pips_io as pipsio
import pyPIPS.thermolib as thermo
import matplotlib.pyplot as plt
from scipy.stats.mstats import zscore

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']


def resample_compass_da(compass_dir, offset, intervalstr):
    """
    Resample compass directions using vector averaging to avoid discontinuity at 360/0 degrees.

    Parameters
    ----------
    compass_dir : xarray.DataArray
        1-Hz compass direction data
    offset : int
        Time offset in seconds for resampling alignment
    intervalstr : str
        Resampling interval string (e.g., '10S')

    Returns
    -------
    compass_dir_avg : xarray.DataArray
        Resampled compass direction data
    """
    offset_str = pips.get_interval_str(offset)
    # Compute x- and y-components of compass direction
    x = np.cos(np.deg2rad(-compass_dir + 270.))
    y = np.sin(np.deg2rad(-compass_dir + 270.))

    x_avg = x.resample(time=intervalstr, label='right', closed='right', offset=offset_str).mean()
    y_avg = y.resample(time=intervalstr, label='right', closed='right', offset=offset_str).mean()

    # Need to use %360 to keep direction between 0 and 360 degrees
    compass_dir_avg = (270.0 - (180. / np.pi) * np.arctan2(y_avg, x_avg)) % 360.
    return compass_dir_avg


def circular_mean_deg(angles_deg):
    """
    Calculate circular mean for compass/wind directions in degrees.

    Uses vector averaging to properly handle the circular nature of directional data.

    Parameters
    ----------
    angles_deg : xarray.DataArray or numpy.ndarray
        Angles in degrees (0-360)

    Returns
    -------
    circular_mean : float
        Circular mean in degrees (0-360)
    """
    # Remove NaN values
    angles_clean = angles_deg[~np.isnan(angles_deg)]

    if len(angles_clean) == 0:
        return np.nan

    # Convert to radians
    angles_rad = np.deg2rad(angles_clean)

    # Calculate mean resultant vector components
    C = np.mean(np.cos(angles_rad))
    S = np.mean(np.sin(angles_rad))

    # Calculate mean angle
    mean_rad = np.arctan2(S, C)
    mean_deg = np.rad2deg(mean_rad)

    # Ensure result is in 0-360 range
    if mean_deg < 0:
        mean_deg += 360.0

    return mean_deg


def circular_std_deg(angles_deg):
    """
    Calculate circular standard deviation for compass/wind directions in degrees.

    Uses the resultant vector length method to properly handle the circular nature
    of directional data (e.g., 359° and 1° are close, not far apart).

    Parameters
    ----------
    angles_deg : xarray.DataArray or numpy.ndarray
        Angles in degrees (0-360)

    Returns
    -------
    circular_std : float
        Circular standard deviation in degrees
    """
    # Remove NaN values
    angles_clean = angles_deg[~np.isnan(angles_deg)]

    if len(angles_clean) == 0:
        return np.nan

    # Convert to radians
    angles_rad = np.deg2rad(angles_clean)

    # Calculate mean resultant vector components
    C = np.mean(np.cos(angles_rad))
    S = np.mean(np.sin(angles_rad))

    # Calculate resultant vector length (R)
    R = np.sqrt(C**2 + S**2)

    # Circular standard deviation
    # Formula: sqrt(-2 * ln(R)) for R > 0
    # Convert back to degrees
    if R > 0:
        circular_std_rad = np.sqrt(-2 * np.log(R))
        circular_std_deg = np.rad2deg(circular_std_rad)
    else:
        # R = 0 means uniform distribution, std = infinite
        circular_std_deg = 180.0  # Maximum possible for circular data

    return circular_std_deg


def detect_deployment_periods(compass_dir, window_size=30, std_threshold=2.0, max_duration=120):
    """
    Detect rapid compass fluctuations at the beginning and end of timeseries that indicate
    probe deployment or retrieval.

    Parameters
    ----------
    compass_dir : xarray.DataArray
        1-Hz compass direction data
    window_size : int
        Window size in seconds for calculating standard deviation (default: 30)
    std_threshold : float
        Standard deviation threshold in degrees for flagging fluctuations (default: 2.0)
    max_duration : int
        Maximum duration in seconds to check at start/end (default: 120 = 2 minutes)

    Returns
    -------
    trim_start : pd.Timestamp or None
        First timestamp to keep (None if no trimming needed at start)
    trim_end : pd.Timestamp or None
        Last timestamp to keep (None if no trimming needed at end)
    """
    # Skip if too few data points
    if len(compass_dir) < window_size * 2:
        return None, None

    trim_start = None
    trim_end = None

    # Check beginning of timeseries
    check_duration_start = min(max_duration, len(compass_dir) // 2)
    for i in range(0, check_duration_start - window_size, window_size // 2):
        window_data = compass_dir.isel(time=slice(i, i + window_size))
        # Calculate circular standard deviation for compass data
        window_std = circular_std_deg(window_data.values)

        if window_std < std_threshold:
            # Found stable period - trim everything before this window
            if i > 0:
                trim_start = compass_dir.time.isel(time=i).values
                utils.log(f"    Detected deployment period: trimming first {i} seconds")
            break

    # Check end of timeseries
    check_duration_end = min(max_duration, len(compass_dir) // 2)
    for i in range(len(compass_dir) - window_size, len(compass_dir) - check_duration_end, -window_size // 2):
        if i < 0:
            break
        window_data = compass_dir.isel(time=slice(i, i + window_size))
        window_std = circular_std_deg(window_data.values)

        if window_std < std_threshold:
            # Found stable period - trim everything after this window
            if i + window_size < len(compass_dir):
                trim_end = compass_dir.time.isel(time=i + window_size - 1).values
                utils.log(f"    Detected retrieval period: trimming last {len(compass_dir) - (i + window_size)} seconds")
            break

    return trim_start, trim_end


def trim_fixed_deployment_periods(conv_ds, trim_duration=30):
    """
    Trim a fixed duration from the beginning and end of the conventional dataset timeseries.
    This removes deployment/retrieval periods without relying on compass variability detection.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset (1-Hz data)
    trim_duration : int or float
        Duration in seconds to trim from start and end (default: 30)

    Returns
    -------
    trim_start : pd.Timestamp or None
        First timestamp to keep (None if dataset too short for trimming)
    trim_end : pd.Timestamp or None
        Last timestamp to keep (None if dataset too short for trimming)
    """
    # Skip if dataset is too short
    if len(conv_ds.time) < trim_duration * 2:
        utils.log(f"    Dataset too short ({len(conv_ds.time)} records) for {trim_duration}s trimming")
        return None, None

    # Calculate trim timestamps
    trim_start = conv_ds.time.values[int(trim_duration)]
    trim_end = conv_ds.time.values[-int(trim_duration)]

    utils.log(f"    Trimming fixed duration: {trim_duration}s from start and end")
    utils.log(f"    Start: {pd.to_datetime(conv_ds.time.values[0])} -> {pd.to_datetime(trim_start)}")
    utils.log(f"    End: {pd.to_datetime(trim_end)} -> {pd.to_datetime(conv_ds.time.values[-1])}")

    return trim_start, trim_end


def trim_datasets_by_time(conv_ds, parsivel_ds, trim_start=None, trim_end=None):
    """
    Trim conventional and parsivel datasets to specified time range and ensure alignment.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset (1-Hz data)
    parsivel_ds : xarray.Dataset or None
        Parsivel dataset (resampled data)
    trim_start : pd.Timestamp or None
        First timestamp to keep (None = keep from beginning)
    trim_end : pd.Timestamp or None
        Last timestamp to keep (None = keep to end)

    Returns
    -------
    conv_ds_trimmed : xarray.Dataset
        Trimmed conventional dataset
    parsivel_ds_trimmed : xarray.Dataset or None
        Trimmed parsivel dataset (None if input was None)
    """
    # Start with full time range
    conv_time_slice = slice(None, None)

    # Apply trimming to conventional dataset
    if trim_start is not None or trim_end is not None:
        conv_time_slice = slice(trim_start, trim_end)
        conv_ds_trimmed = conv_ds.sel(time=conv_time_slice)
        utils.log(f"    Conventional dataset trimmed: {len(conv_ds.time)} -> {len(conv_ds_trimmed.time)} records")
    else:
        conv_ds_trimmed = conv_ds

    # If parsivel dataset exists, trim it to match and then align conventional to parsivel times
    if parsivel_ds is not None:
        # Trim parsivel to same time range
        parsivel_ds_trimmed = parsivel_ds.sel(time=conv_time_slice)
        utils.log(f"    Parsivel dataset trimmed: {len(parsivel_ds.time)} -> {len(parsivel_ds_trimmed.time)} records")

        # Further trim conventional dataset to align with parsivel times
        # (parsivel times are resampled, so conventional should cover the same range)
        if len(parsivel_ds_trimmed) > 0:
            parsivel_start = parsivel_ds_trimmed.time.values[0]
            parsivel_end = parsivel_ds_trimmed.time.values[-1]

            # Trim conventional to match parsivel time range
            conv_ds_trimmed = conv_ds_trimmed.sel(time=slice(parsivel_start, parsivel_end))
            utils.log(f"    Conventional dataset aligned to parsivel times: {len(conv_ds_trimmed)} records")
    else:
        parsivel_ds_trimmed = None

    return conv_ds_trimmed, parsivel_ds_trimmed


def update_time_attributes(dataset, time_format='%Y%m%d%H%M%S'):
    """
    Update starting_time and ending_time global attributes based on actual data times.
    Also updates the time coordinate encoding to use the new start time as reference.

    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset to update
    time_format : str
        Format string for time attributes (default: '%Y%m%d%H%M%S')

    Returns
    -------
    dataset : xarray.Dataset
        Dataset with updated time attributes and time coordinate encoding
    """
    if len(dataset.time) > 0:
        # Save original time attributes
        old_start = dataset.attrs.get('starting_time', 'N/A')
        old_end = dataset.attrs.get('ending_time', 'N/A')

        start_time = pd.to_datetime(dataset.time.values[0])
        end_time = pd.to_datetime(dataset.time.values[-1])

        dataset.attrs['starting_time'] = start_time.strftime(time_format)
        dataset.attrs['ending_time'] = end_time.strftime(time_format)

        # Update time coordinate encoding to use new start time as reference
        # This ensures netCDF serialization uses "seconds since [new_start_time]"
        time_units = f"seconds since {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
        dataset['time'].encoding['units'] = time_units
        dataset['time'].encoding['calendar'] = 'proleptic_gregorian'

        utils.log(f"    Updated time attributes:")
        utils.log(f"      Original: {old_start} to {old_end}")
        utils.log(f"      New:      {dataset.attrs['starting_time']} to {dataset.attrs['ending_time']}")
        utils.log(f"    Updated time coordinate encoding: {time_units}")

    return dataset


def apply_compass_qc(conv_ds, zscore_threshold=3.0, abs_threshold=20.0):
    """
    Apply quality control to compass headings using circular statistics and dual threshold.

    Uses both z-score and absolute deviation criteria to avoid over-cleaning stable data.
    A point is flagged as an outlier only if it exceeds BOTH thresholds.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset with compass_dir variable
    zscore_threshold : float
        Z-score threshold for outlier detection (default: 3.0)
    abs_threshold : float
        Minimum absolute deviation in degrees to be considered an outlier (default: 20.0)

    Returns
    -------
    cleaned_compass_dir : xarray.DataArray
        Compass direction with outliers replaced by NaN
    avg_compass_dir : float
        Circular mean compass direction after removing outliers
    circ_std : float
        Circular standard deviation of compass data
    n_outliers : int
        Number of outliers removed
    """
    compass_data = conv_ds['compass_dir'].values

    # Calculate circular mean and std
    circ_mean = circular_mean_deg(compass_data)
    circ_std = circular_std_deg(compass_data)

    # Calculate circular deviations (handling 360/0 wrap)
    deviations = compass_data - circ_mean
    # Adjust for circular wrap (e.g., 350° - 10° should be -20°, not 340°)
    deviations = np.where(deviations > 180, deviations - 360, deviations)
    deviations = np.where(deviations < -180, deviations + 360, deviations)

    # Calculate absolute deviations
    abs_deviations = np.abs(deviations)

    # Calculate z-scores using circular std
    if circ_std > 0:
        zscores = deviations / circ_std
    else:
        # If std is zero, data is perfectly stable - no outliers
        zscores = np.zeros_like(deviations)

    # Dual threshold: flag as outlier only if BOTH conditions are met
    is_outlier = (np.abs(zscores) > zscore_threshold) & (abs_deviations > abs_threshold)
    n_outliers = np.sum(is_outlier)

    # Create cleaned version
    cleaned_compass_dir = conv_ds['compass_dir'].copy()
    cleaned_compass_dir.values[is_outlier] = np.nan

    # Recalculate circular mean after removing outliers
    avg_compass_dir = circular_mean_deg(cleaned_compass_dir.values)

    return cleaned_compass_dir, avg_compass_dir, circ_std, n_outliers


def apply_slowtemp_qc(conv_ds, diff_threshold=2.0):
    """
    Apply quality control to slow temperature observations by removing data points
    where the difference between slowtemp and fasttemp exceeds a threshold.

    Large differences between the two temperature sensors may indicate sensor
    malfunctions, contamination, or rapid temperature changes where the sensors'
    different response times cause systematic bias.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset with slowtemp and fasttemp variables
    diff_threshold : float
        Maximum allowed absolute difference in °C (default: 2.0)

    Returns
    -------
    cleaned_slowtemp : xarray.DataArray
        Slow temperature with outliers replaced by NaN
    n_flagged : int
        Number of observations flagged
    mean_diff : float
        Mean difference (slowtemp - fasttemp) before filtering
    """
    slowtemp_data = conv_ds['slowtemp'].values
    fasttemp_data = conv_ds['fasttemp'].values

    # Calculate temperature differences
    temp_diff = np.abs(slowtemp_data - fasttemp_data)

    # Calculate mean difference before filtering (for diagnostics)
    mean_diff = np.nanmean(slowtemp_data - fasttemp_data)

    # Flag observations exceeding threshold
    is_outlier = temp_diff > diff_threshold
    n_flagged = np.sum(is_outlier)

    # Create cleaned version with outliers set to NaN
    cleaned_slowtemp = conv_ds['slowtemp'].copy()
    cleaned_slowtemp.values[is_outlier] = np.nan

    return cleaned_slowtemp, n_flagged, mean_diff


def recompute_dewpoint_RH_with_fallback(conv_ds, reason=''):
    """
    Recompute dewpoint and RH_derived using slowtemp with fasttemp fallback.

    For dewpoint calculation: uses slowtemp where available, fasttemp where slowtemp is NaN.
    For RH_derived calculation: uses fasttemp and the newly computed dewpoint.

    This handles cases where slowtemp has been set to NaN by QC or manual editing,
    ensuring thermodynamic consistency by falling back to fasttemp.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset with slowtemp, fasttemp, pressure, and RH variables
    reason : str, optional
        Reason for recomputation (added to variable attributes)

    Returns
    -------
    conv_ds : xarray.Dataset
        Dataset with recomputed dewpoint and RH_derived
    n_fallback : int
        Number of times fasttemp was used as fallback
    """
    # Check if required variables exist
    required_vars = ['slowtemp', 'fasttemp', 'pressure', 'RH', 'dewpoint', 'RH_derived']
    missing_vars = [var for var in required_vars if var not in conv_ds]
    if missing_vars:
        utils.log(f"  Warning: Cannot recompute dewpoint/RH, missing variables: {missing_vars}")
        return conv_ds, 0

    pressure = conv_ds['pressure']  # hPa
    slowtemp = conv_ds['slowtemp']  # °C
    fasttemp = conv_ds['fasttemp']  # °C
    RH = conv_ds['RH']  # %

    # Create temperature array for dewpoint calculation:
    # Use slowtemp where available, fasttemp where slowtemp is NaN
    temp_for_dewpoint = slowtemp.where(~np.isnan(slowtemp), fasttemp)

    # Count how many times fasttemp is used as fallback
    n_fallback = int(np.sum(np.isnan(slowtemp.values) & ~np.isnan(fasttemp.values)))

    # Calculate dewpoint from RH and temperature
    # thermo.calTdfromRH expects: pressure (Pa), temperature (K), RH (fraction)
    dewpoint = thermo.calTdfromRH(
        pressure * 100.,  # Convert hPa to Pa
        temp_for_dewpoint + 273.15,  # Convert °C to K
        RH / 100.  # Convert % to fraction
    ) - 273.15  # Convert K back to °C

    # Recompute RH_derived using fasttemp and newly calculated dewpoint
    # thermo.calRH expects: pressure (Pa), temperature (K), dewpoint (K)
    RH_derived = thermo.calRH(
        pressure * 100.,  # Convert hPa to Pa
        fasttemp + 273.15,  # Convert °C to K
        dewpoint + 273.15  # Convert °C to K
    ) * 100.  # Convert fraction to %

    # Update dataset with recomputed thermodynamic variables
    conv_ds['dewpoint'].data = dewpoint.data
    if reason:
        conv_ds['dewpoint'].attrs['recomputation_reason'] = reason

    conv_ds['RH_derived'].data = RH_derived.data
    if reason:
        conv_ds['RH_derived'].attrs['recomputation_reason'] = reason

    return conv_ds, n_fallback


def apply_slowtemp_bias_correction(conv_ds, slope, intercept):
    """
    Apply bias correction to slow temperature observations using pre-computed
    linear regression coefficients. Also recomputes dewpoint and RH_derived
    based on the corrected slowtemp.

    The correction inverts the linear relationship: slowtemp = slope * fasttemp + intercept
    to obtain: slowtemp_corrected = (slowtemp - intercept) / slope

    For dewpoint calculation, if slowtemp is NaN but fasttemp is available,
    fasttemp is used as a fallback.

    Parameters
    ----------
    conv_ds : xarray.Dataset
        Conventional dataset with slowtemp, fasttemp, pressure, and RH variables
    slope : float
        Linear regression slope from slowtemp vs fasttemp analysis
    intercept : float
        Linear regression intercept from slowtemp vs fasttemp analysis

    Returns
    -------
    conv_ds : xarray.Dataset
        Dataset with corrected slowtemp (overwritten), dewpoint, and RH_derived
    """
    # Apply bias correction to slowtemp: invert linear relationship
    # Original relationship: slowtemp = slope * fasttemp + intercept
    # Corrected: slowtemp_corrected = (slowtemp - intercept) / slope
    slowtemp_corrected = (conv_ds['slowtemp'] - intercept) / slope

    # Add slope and intercept as attributes
    slowtemp_corrected.attrs.update(conv_ds['slowtemp'].attrs)
    slowtemp_corrected.attrs['bias_correction_slope'] = slope
    slowtemp_corrected.attrs['bias_correction_intercept'] = intercept
    slowtemp_corrected.attrs['bias_corrected'] = True

    # Overwrite slowtemp with corrected values
    conv_ds['slowtemp'] = slowtemp_corrected

    # Recompute dewpoint and RH_derived using the reusable function
    conv_ds, n_fallback = recompute_dewpoint_RH_with_fallback(
        conv_ds, reason='recomputed_from_bias_corrected_slowtemp'
    )

    if n_fallback > 0:
        utils.log(f"  Used fasttemp fallback for {n_fallback} points where slowtemp was NaN")

    return conv_ds


def plot_compass_qc_diagnostics(compass_dir_original, winddirabs_original,
                                cleaned_compass_dir, avg_compass_dir,
                                winddirabs_new, PIPS_name, deployment_name, output_dir):
    """
    Generate diagnostic plots showing before/after compass and wind direction QC.

    Parameters
    ----------
    compass_dir_original : xarray.DataArray
        Original (unmodified) compass directions
    winddirabs_original : xarray.DataArray
        Original (unmodified) wind directions
    cleaned_compass_dir : xarray.DataArray
        QC'd compass directions
    avg_compass_dir : float
        Average compass direction after QC
    winddirabs_new : xarray.DataArray
        Recomputed wind directions
    PIPS_name : str
        PIPS station name
    deployment_name : str
        Deployment name
    output_dir : str
        Directory to save plots
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f'{PIPS_name} - Compass and Wind Direction QC', fontsize=14)

    # Plot 1: Original compass direction
    ax = axes[0, 0]
    compass_dir_original.plot(ax=ax, ls='None', marker='o', ms=1., label='Original')
    ax.set_title('Original Compass Direction')
    ax.set_ylabel('Compass Direction (deg)')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 2: Cleaned compass direction
    ax = axes[0, 1]
    cleaned_compass_dir.plot(ax=ax, ls='None', marker='o', ms=1., color='orange', label='QC\'d')
    ax.axhline(avg_compass_dir, color='red', ls='--', lw=1.5, label=f'Mean: {avg_compass_dir:.1f}°')
    ax.set_title('QC\'d Compass Direction')
    ax.set_ylabel('Compass Direction (deg)')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 3: Original wind direction
    ax = axes[1, 0]
    winddirabs_original.plot(ax=ax, ls='None', marker='o', ms=1., label='Original')
    ax.set_ylim(0., 360.)
    ax.set_title('Original Wind Direction')
    ax.set_ylabel('Wind Direction (deg)')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Plot 4: Corrected wind direction
    ax = axes[1, 1]
    winddirabs_new.plot(ax=ax, ls='None', marker='o', ms=1., color='green', label='QC\'d')
    ax.set_ylim(0., 360.)
    ax.set_title('QC\'d Wind Direction')
    ax.set_ylabel('Wind Direction (deg)')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plot_filename = os.path.join(output_dir, f'compass_wind_QC_{deployment_name}_{PIPS_name}.png')
    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
    utils.log(f"Saved compass QC diagnostic plot: {plot_filename}")
    plt.close(fig)


def plot_slowtemp_qc_diagnostics(slowtemp_original, fasttemp_original,
                                 slowtemp_cleaned, fasttemp_current,
                                 PIPS_name, deployment_name, output_dir):
    """
    Generate diagnostic plots showing before/after slowtemp QC.

    Parameters
    ----------
    slowtemp_original : xarray.DataArray
        Original (unmodified) slow temperature
    fasttemp_original : xarray.DataArray
        Original (unmodified) fast temperature
    slowtemp_cleaned : xarray.DataArray
        QC'd slow temperature
    fasttemp_current : xarray.DataArray
        Current fast temperature (for reference)
    PIPS_name : str
        PIPS station name
    deployment_name : str
        Deployment name
    output_dir : str
        Directory to save plots
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f'{PIPS_name} - Slow Temperature QC', fontsize=14)

    # Plot 1: Original slowtemp timeseries
    ax = axes[0, 0]
    slowtemp_original.plot(ax=ax, label='Slowtemp (original)', color='blue', lw=0.5)
    fasttemp_original.plot(ax=ax, label='Fasttemp', color='red', lw=0.5, alpha=0.7)
    ax.set_title('Original Temperature Data')
    ax.set_ylabel('Temperature (°C)')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')

    # Plot 2: QC'd slowtemp timeseries
    ax = axes[0, 1]
    slowtemp_cleaned.plot(ax=ax, label='Slowtemp (QC\'d)', color='green', lw=0.5)
    fasttemp_current.plot(ax=ax, label='Fasttemp', color='red', lw=0.5, alpha=0.7)
    ax.set_title('QC\'d Temperature Data')
    ax.set_ylabel('Temperature (°C)')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')

    # Plot 3: Original temperature difference
    ax = axes[1, 0]
    temp_diff_original = slowtemp_original - fasttemp_original
    temp_diff_original.plot(ax=ax, color='gray', lw=0.5)
    ax.axhline(0, color='black', ls='--', lw=1)
    ax.set_title('Original Temperature Difference (Slow - Fast)')
    ax.set_ylabel('Temperature Difference (°C)')
    ax.grid(True, alpha=0.3)

    # Plot 4: QC'd temperature difference
    ax = axes[1, 1]
    temp_diff_cleaned = slowtemp_cleaned - fasttemp_current
    temp_diff_cleaned.plot(ax=ax, color='green', lw=0.5)
    ax.axhline(0, color='black', ls='--', lw=1)
    ax.set_title('QC\'d Temperature Difference (Slow - Fast)')
    ax.set_ylabel('Temperature Difference (°C)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plot_filename = os.path.join(output_dir, f'slowtemp_QC_{deployment_name}_{PIPS_name}.png')
    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
    utils.log(f"Saved slowtemp QC diagnostic plot: {plot_filename}")
    plt.close(fig)


# Parse the command line options
description = "Calculates various derived parameters from PIPS DSDs (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--input-QC-tag', dest='input_QC_tag', default=None,
                    help='Tag for input ND variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--output-QC-tags', dest='output_QC_tags', nargs='*', default=['qc'],
                    help='list of QC groups to apply')
parser.add_argument('--output-file-tag', dest='output_file_tag', default='',
                    help='tag for output nc files to distinguish from original if desired')
parser.add_argument('--compass-qc', dest='compass_qc', action='store_true',
                    help='Apply compass quality control to remove outliers and recompute wind directions')
parser.add_argument('--compass-zscore-threshold', dest='compass_zscore_threshold', type=float, default=3.0,
                    help='Z-score threshold for compass outlier detection (default: 3.0)')
parser.add_argument('--compass-abs-threshold', dest='compass_abs_threshold', type=float, default=20.0,
                    help='Minimum absolute deviation in degrees for outlier detection (default: 20.0)')
parser.add_argument('--plot-compass-qc', dest='plot_compass_qc', action='store_true',
                    help='Generate diagnostic plots showing before/after compass and wind QC')
parser.add_argument('--trim-deployment-periods', dest='trim_deployment_periods', action='store_true',
                    help='Remove fixed duration from start/end of deployment to exclude deployment/retrieval periods')
parser.add_argument('--trim-duration', dest='trim_duration', type=float, default=30.0,
                    help='Duration in seconds to trim from start and end of deployment (default: 30.0)')
# Commented out: Old automatic detection-based trimming arguments (kept for potential future use)
# parser.add_argument('--trim-window-size', dest='trim_window_size', type=int, default=30,
#                     help='Window size in seconds for detecting compass fluctuations (default: 30)')
# parser.add_argument('--trim-std-threshold', dest='trim_std_threshold', type=float, default=2.0,
#                     help='Std dev threshold in degrees for flagging compass fluctuations (default: 2.0)')
# parser.add_argument('--trim-max-duration', dest='trim_max_duration', type=int, default=120,
#                     help='Maximum duration in seconds to check at start/end (default: 120)')
parser.add_argument('--slowtemp-qc', dest='slowtemp_qc', action='store_true',
                    help='Apply quality control to slowtemp by removing observations with large temp differences')
parser.add_argument('--slowtemp-diff-threshold', dest='slowtemp_diff_threshold', type=float, default=2.0,
                    help='Maximum allowed temperature difference between slowtemp and fasttemp in °C (default: 2.0)')
parser.add_argument('--plot-slowtemp-qc', dest='plot_slowtemp_qc', action='store_true',
                    help='Generate diagnostic plots showing before/after slowtemp QC')
parser.add_argument('--slowtemp-bias-correction', dest='slowtemp_bias_correction', action='store_true',
                    help='Apply bias correction to slowtemp using pre-computed slopes and intercepts')
parser.add_argument('--slowtemp-bias-config', dest='slowtemp_bias_config', type=str, default=None,
                    help='Path to config file with slowtemp bias correction parameters (default: None)')
parser.add_argument('--recompute-dewpoint-with-fallback', dest='recompute_dewpoint', action='store_true',
                    help='Recompute dewpoint and RH_derived for any points where slowtemp is NaN, '
                         'using fasttemp as fallback (handles manual QC done outside script)')

args = parser.parse_args()
if args.input_QC_tag:
    input_QC_tag = '_{}'.format(args.input_QC_tag)
    if 'RB15' in args.input_QC_tag and 'vshift' not in args.input_QC_tag:
        VD_tag = args.input_QC_tag.replace('RB15', 'RB15_vshift')
        VD_tag = '_{}'.format(VD_tag)
    else:
        VD_tag = '_{}'.format(args.input_QC_tag)
else:
    input_QC_tag = ''
    VD_tag = ''

# if not args.output_QC_tags:
#     output_QC_tags = [input_QC_tag]
# else:
#     output_QC_tags = ['_{}'.format(output_QC_tag) for output_QC_tag in args.output_QC_tags]

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
parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf in parsivel_combined_filenames]

# Get conventional data filenames if they exist in config
try:
    conv_filenames = config.PIPS_IO_dict.get('conv_filenames_nc', None)
    if conv_filenames is not None:
        conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames]
    else:
        conv_filelist = None
except Exception:
    conv_filelist = None

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print("Reading {}".format(parsivel_combined_file))
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)

    DSD_interval = parsivel_combined_ds.DSD_interval
    PIPS_name = parsivel_combined_ds.probe_name
    deployment_name = parsivel_combined_ds.deployment_name
    ND_name = 'ND{}'.format(input_QC_tag)
    ND_da = parsivel_combined_ds[ND_name]
    VD_name = 'VD_matrix{}'.format(VD_tag)
    vd_matrix_da = parsivel_combined_ds[VD_name]
    coord_to_combine = 'time'
    flagged_times = None

    # Load conventional dataset if compass/wind QC or trimming is requested
    conv_ds = None
    if (args.compass_qc or args.trim_deployment_periods or args.slowtemp_qc or args.recompute_dewpoint) and conv_filelist is not None:
        conv_file = conv_filelist[index]
        if os.path.exists(conv_file):
            conv_ds = xr.load_dataset(conv_file)
            utils.log(f"Loaded conventional data for {PIPS_name}")

            # Save copies of original compass/wind data for diagnostic plotting
            # (before any QC or trimming) if plotting is requested
            if args.plot_compass_qc and 'compass_dir' in conv_ds and 'winddirabs' in conv_ds:
                compass_dir_original = conv_ds['compass_dir'].copy(deep=True)
                winddirabs_original = conv_ds['winddirabs'].copy(deep=True)

            # Save copies of original temperature data for diagnostic plotting
            # (before any QC or trimming) if plotting is requested
            if args.plot_slowtemp_qc and 'slowtemp' in conv_ds and 'fasttemp' in conv_ds:
                slowtemp_original = conv_ds['slowtemp'].copy(deep=True)
                fasttemp_original = conv_ds['fasttemp'].copy(deep=True)
        else:
            utils.log(f"Warning: Conventional file not found: {conv_file}")

    # Do some QC on the V-D matrix. This will make a copy of the raw matrix. The netCDF file
    # will contain the original raw VD matrix and the new QC'ed matrix with the appropriate
    # QC tag applied
    for output_QC_tag in args.output_QC_tags:
        if output_QC_tag not in pqc.PIPS_qc_dict:
            continue

        strongwindQC = pqc.PIPS_qc_dict[output_QC_tag]['strongwindQC']
        splashingQC = pqc.PIPS_qc_dict[output_QC_tag]['splashingQC']
        marginQC = pqc.PIPS_qc_dict[output_QC_tag]['marginQC']
        rainfallQC = pqc.PIPS_qc_dict[output_QC_tag]['rainfallQC']
        rainonlyQC = pqc.PIPS_qc_dict[output_QC_tag]['rainonlyQC']
        hailonlyQC = pqc.PIPS_qc_dict[output_QC_tag]['hailonlyQC']
        graupelonlyQC = pqc.PIPS_qc_dict[output_QC_tag]['graupelonlyQC']

        vd_matrix_qc_da = vd_matrix_da.copy(deep=True)
        if strongwindQC:
            vd_matrix_qc_da, flagged_times = pqc.strongwindQC(vd_matrix_qc_da)
        if splashingQC:
            vd_matrix_qc_da = pqc.splashingQC(vd_matrix_qc_da)
        if marginQC:
            vd_matrix_qc_da = pqc.marginQC(vd_matrix_qc_da)
        if rainfallQC:
            fallspeedmask = pqc.get_fallspeed_mask(avg_diameter, avg_fall_bins)
            vd_matrix_qc_da = pqc.rainfallspeedQC(vd_matrix_qc_da, fallspeedmask)
        if rainonlyQC:
            vd_matrix_qc_da = pqc.rainonlyQC(vd_matrix_qc_da)
        if hailonlyQC:
            vd_matrix_qc_da = pqc.hailonlyQC(vd_matrix_qc_da)

        fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins,
                                                          correct_rho=True,
                                                          rho=parsivel_combined_ds['rho'])

        vd_matrix_qc_da = vd_matrix_qc_da.where(vd_matrix_qc_da > 0.0)
        ND_qc_da = pips.calc_ND(vd_matrix_qc_da, fallspeed_spectrum, DSD_interval)

        if input_QC_tag == output_QC_tag:
            new_VD_name = VD_name
            new_ND_name = ND_name
        else:
            new_VD_name = '{}_{}'.format(VD_name, output_QC_tag)
            new_ND_name = '{}_{}'.format(ND_name, output_QC_tag)

        parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, vd_matrix_qc_da,
                                                            name=new_VD_name)
        parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, ND_qc_da,
                                                            name=new_ND_name)
        # Update metadata
        for varname in [new_VD_name, new_ND_name]:
            parsivel_combined_ds[varname].attrs['strongwindQC'] = int(strongwindQC)
            parsivel_combined_ds[varname].attrs['splashingQC'] = int(splashingQC)
            parsivel_combined_ds[varname].attrs['marginQC'] = int(marginQC)
            parsivel_combined_ds[varname].attrs['rainfallQC'] = int(rainfallQC)
            parsivel_combined_ds[varname].attrs['rainonlyQC'] = int(rainonlyQC)
            parsivel_combined_ds[varname].attrs['hailonlyQC'] = int(hailonlyQC)
            parsivel_combined_ds[varname].attrs['graupelonlyQC'] = int(graupelonlyQC)

        # Add flagged times variable if available
        if flagged_times is not None:
            parsivel_combined_ds['flagged_times_{}'.format(output_QC_tag)] = (
                ('time',), flagged_times)
            parsivel_combined_ds['flagged_times_{}'.format(output_QC_tag)].attrs['description'] = (
                'Flagged times from QC: 0=good, 2=severe wind contamination')

    # =============================================================================================
    # COMPASS AND WIND DIRECTION QUALITY CONTROL
    # =============================================================================================
    if args.compass_qc and conv_ds is not None:
        utils.log(f"Applying compass quality control for {PIPS_name}...")

        # Apply compass QC to remove outliers using circular statistics and dual threshold
        cleaned_compass_dir, avg_compass_dir, circ_std, n_outliers = apply_compass_qc(
            conv_ds, zscore_threshold=args.compass_zscore_threshold,
            abs_threshold=args.compass_abs_threshold)

        original_mean = circular_mean_deg(conv_ds['compass_dir'].values)
        utils.log(f"  Original circular mean: {original_mean:.2f}°, circular std: {circ_std:.2f}°")
        utils.log(f"  QC'd circular mean: {avg_compass_dir:.2f}°")
        utils.log(f"  Outliers removed: {n_outliers} ({100*n_outliers/len(conv_ds['compass_dir']):.2f}%)")

        # Recompute absolute wind directions using cleaned compass
        winddirabs_new = np.mod(avg_compass_dir + conv_ds['winddirrel'], 360.)

        # Update conventional dataset with cleaned compass and recomputed winds
        conv_ds['compass_dir'] = cleaned_compass_dir
        conv_ds['compass_dir'].attrs['average'] = avg_compass_dir
        conv_ds['compass_dir'].attrs['circular_std'] = circ_std
        conv_ds['compass_dir'].attrs['qc_zscore_threshold'] = args.compass_zscore_threshold
        conv_ds['compass_dir'].attrs['qc_abs_threshold'] = args.compass_abs_threshold
        conv_ds['compass_dir'].attrs['n_outliers_removed'] = n_outliers
        conv_ds['winddirabs'] = winddirabs_new

        # If parsivel data exists, resample winds and compass to parsivel times
        if parsivel_combined_ds is not None and 'VD_matrix' in parsivel_combined_ds:
            try:
                PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
                sec_offset = PSD_datetimes[0].second
                DSD_interval = parsivel_combined_ds.DSD_interval
                interval_str = pips.get_interval_str(DSD_interval)

                # Resample winds with new corrected compass
                new_wind_ds = pips.resample_wind_da(winddirabs_new, conv_ds['windspd'],
                                                    interval_str, sec_offset,
                                                    gusts=True, gustintvstr='3S')

                # Resample cleaned compass directions
                cleaned_compass_dir_avg = resample_compass_da(cleaned_compass_dir,
                                                              sec_offset, interval_str)

                # Update parsivel dataset with new wind and compass data
                for key in new_wind_ds.data_vars:
                    parsivel_combined_ds[key] = new_wind_ds[key]

                parsivel_combined_ds['compass_dir'] = cleaned_compass_dir_avg
                avg_compass_dir_resampled = circular_mean_deg(cleaned_compass_dir_avg.values)
                parsivel_combined_ds['compass_dir'].attrs['average'] = avg_compass_dir_resampled
                parsivel_combined_ds['compass_dir'].attrs['qc_zscore_threshold'] = args.compass_zscore_threshold
                parsivel_combined_ds['compass_dir'].attrs['qc_abs_threshold'] = args.compass_abs_threshold

                utils.log(f"  Updated parsivel winds and compass (resampled mean: {avg_compass_dir_resampled:.2f}°)")

            except Exception as e:
                utils.log(f"Warning: Could not resample winds for {PIPS_name}: {e}")

    # =============================================================================================
    # SLOW TEMPERATURE QUALITY CONTROL
    # =============================================================================================
    if args.slowtemp_qc and conv_ds is not None:
        utils.log(f"Applying slow temperature quality control for {PIPS_name}...")

        # Check if required temperature variables exist
        if 'slowtemp' in conv_ds and 'fasttemp' in conv_ds:
            # Apply slowtemp QC to remove observations with large differences
            cleaned_slowtemp, n_flagged, mean_diff = apply_slowtemp_qc(
                conv_ds, diff_threshold=args.slowtemp_diff_threshold)

            n_total = np.sum(~np.isnan(conv_ds['slowtemp'].values))
            utils.log(f"  Mean difference (slowtemp - fasttemp): {mean_diff:.3f}°C")
            utils.log(f"  Observations flagged: {n_flagged} ({100*n_flagged/n_total:.2f}%)")
            utils.log(f"  Threshold used: {args.slowtemp_diff_threshold}°C")

            # Update conventional dataset with cleaned slowtemp
            conv_ds['slowtemp'] = cleaned_slowtemp
            conv_ds['slowtemp'].attrs['qc_diff_threshold'] = args.slowtemp_diff_threshold
            conv_ds['slowtemp'].attrs['n_flagged'] = n_flagged
            conv_ds['slowtemp'].attrs['mean_diff_before_qc'] = mean_diff

            # Recompute dewpoint and RH_derived for points where slowtemp was set to NaN
            if n_flagged > 0:
                utils.log(f"  Recomputing dewpoint/RH_derived for flagged points...")
                conv_ds, n_fallback = recompute_dewpoint_RH_with_fallback(
                    conv_ds, reason='recomputed_after_slowtemp_qc'
                )
                if n_fallback > 0:
                    utils.log(f"  Used fasttemp fallback for {n_fallback} points")

            # If parsivel data exists and has slowtemp/dewpoint/RH_derived, update them by resampling
            if parsivel_combined_ds is not None and 'slowtemp' in parsivel_combined_ds:
                try:
                    PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
                    sec_offset = PSD_datetimes[0].second
                    DSD_interval = parsivel_combined_ds.DSD_interval
                    interval_str = pips.get_interval_str(DSD_interval)
                    offset_str = pips.get_interval_str(sec_offset)

                    # Resample cleaned slowtemp to parsivel times
                    slowtemp_resampled = cleaned_slowtemp.resample(
                        time=interval_str, label='right', closed='right', offset=offset_str
                    ).mean()

                    # Update parsivel dataset
                    parsivel_combined_ds['slowtemp'] = slowtemp_resampled
                    parsivel_combined_ds['slowtemp'].attrs['qc_diff_threshold'] = args.slowtemp_diff_threshold

                    utils.log(f"  Updated parsivel slowtemp (resampled from cleaned data)")

                    # Also resample dewpoint and RH_derived if they were recomputed
                    if n_flagged > 0:
                        if 'dewpoint' in parsivel_combined_ds:
                            dewpoint_resampled = conv_ds['dewpoint'].resample(
                                time=interval_str, label='right', closed='right', offset=offset_str
                            ).mean()
                            parsivel_combined_ds['dewpoint'] = dewpoint_resampled
                            utils.log(f"  Updated parsivel dewpoint (resampled from recomputed data)")

                        if 'RH_derived' in parsivel_combined_ds:
                            RH_derived_resampled = conv_ds['RH_derived'].resample(
                                time=interval_str, label='right', closed='right', offset=offset_str
                            ).mean()
                            parsivel_combined_ds['RH_derived'] = RH_derived_resampled
                            utils.log(f"  Updated parsivel RH_derived (resampled from recomputed data)")

                except Exception as e:
                    utils.log(f"Warning: Could not resample slowtemp for {PIPS_name}: {e}")
        else:
            utils.log(f"Warning: slowtemp or fasttemp not found in {PIPS_name} dataset")

    # =============================================================================================
    # RECOMPUTE DEWPOINT/RH FOR ANY NaN SLOWTEMP VALUES (GENERAL CLEANUP)
    # =============================================================================================
    if args.recompute_dewpoint and conv_ds is not None:
        utils.log(f"Checking for NaN slowtemp values in {PIPS_name}...")

        # Check if required temperature variables exist
        if 'slowtemp' in conv_ds and 'fasttemp' in conv_ds:
            # Count how many points have NaN slowtemp but valid fasttemp
            slowtemp_data = conv_ds['slowtemp'].values
            fasttemp_data = conv_ds['fasttemp'].values
            n_nan_slowtemp = np.sum(np.isnan(slowtemp_data) & ~np.isnan(fasttemp_data))

            if n_nan_slowtemp > 0:
                utils.log(f"  Found {n_nan_slowtemp} points with NaN slowtemp and valid fasttemp")
                utils.log(f"  Recomputing dewpoint and RH_derived with fasttemp fallback...")

                conv_ds, n_fallback = recompute_dewpoint_RH_with_fallback(
                    conv_ds, reason='recomputed_with_fasttemp_fallback_for_nan_slowtemp'
                )

                utils.log(f"  ✓ Recomputed for {n_fallback} points")

                # If parsivel data exists, resample updated dewpoint/RH_derived
                if parsivel_combined_ds is not None:
                    try:
                        PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
                        sec_offset = PSD_datetimes[0].second
                        DSD_interval = parsivel_combined_ds.DSD_interval
                        interval_str = pips.get_interval_str(DSD_interval)
                        offset_str = pips.get_interval_str(sec_offset)

                        if 'dewpoint' in parsivel_combined_ds:
                            dewpoint_resampled = conv_ds['dewpoint'].resample(
                                time=interval_str, label='right', closed='right', offset=offset_str
                            ).mean()
                            parsivel_combined_ds['dewpoint'] = dewpoint_resampled
                            parsivel_combined_ds['dewpoint'].attrs.update(conv_ds['dewpoint'].attrs)

                        if 'RH_derived' in parsivel_combined_ds:
                            RH_derived_resampled = conv_ds['RH_derived'].resample(
                                time=interval_str, label='right', closed='right', offset=offset_str
                            ).mean()
                            parsivel_combined_ds['RH_derived'] = RH_derived_resampled
                            parsivel_combined_ds['RH_derived'].attrs.update(conv_ds['RH_derived'].attrs)

                        utils.log(f"  Updated parsivel dataset with recomputed values")

                    except Exception as e:
                        utils.log(f"Warning: Could not resample recomputed values for {PIPS_name}: {e}")
            else:
                utils.log(f"  No NaN slowtemp values found (all OK)")
        else:
            utils.log(f"Warning: slowtemp or fasttemp not found in {PIPS_name} dataset")

    # =============================================================================================
    # DEPLOYMENT/RETRIEVAL PERIOD TRIMMING
    # =============================================================================================
    if args.trim_deployment_periods and conv_ds is not None:
        utils.log(f"Trimming deployment/retrieval periods for {PIPS_name}...")

        # Trim fixed duration from start and end
        trim_start, trim_end = trim_fixed_deployment_periods(
            conv_ds,
            trim_duration=args.trim_duration
        )

        # Commented out: Old automatic detection-based trimming (kept for potential future use)
        # # Detect rapid fluctuations at start and end
        # trim_start, trim_end = detect_deployment_periods(
        #     conv_ds['compass_dir'],
        #     window_size=args.trim_window_size,
        #     std_threshold=args.trim_std_threshold,
        #     max_duration=args.trim_max_duration
        # )

        # If trimming is needed, apply it to both datasets
        if trim_start is not None or trim_end is not None:
            utils.log(f"  Applying trimming to datasets for {PIPS_name}...")

            # Trim both conventional and parsivel datasets
            conv_ds, parsivel_combined_ds = trim_datasets_by_time(
                conv_ds, parsivel_combined_ds, trim_start, trim_end
            )

            # Update time attributes for both datasets
            conv_ds = update_time_attributes(conv_ds)
            parsivel_combined_ds = update_time_attributes(parsivel_combined_ds)
        else:
            utils.log(f"  Dataset too short for trimming {PIPS_name}")
            # Old message for detection-based approach:
            # utils.log(f"  No deployment/retrieval periods detected for {PIPS_name}")

    # =============================================================================================
    # SLOW TEMPERATURE BIAS CORRECTION
    # =============================================================================================
    if args.slowtemp_bias_correction and conv_ds is not None:
        utils.log(f"Applying slow temperature bias correction for {PIPS_name}...")

        # Check if bias correction config was provided
        if args.slowtemp_bias_config is None:
            utils.log(f"Warning: --slowtemp-bias-correction specified but no config file provided")
            utils.log(f"         Use --slowtemp-bias-config to specify path to bias correction config")
        else:
            # Import bias correction parameters from config file
            try:
                bias_config = utils.import_all_from(args.slowtemp_bias_config)
                slowtemp_bias_dict = bias_config.slowtemp_bias

                # Get slope and intercept for this PIPS
                if PIPS_name in slowtemp_bias_dict:
                    slope, intercept = slowtemp_bias_dict[PIPS_name]

                    # Only apply correction if slope != 1.0 or intercept != 0.0
                    if slope != 1.0 or intercept != 0.0:
                        # Check if required variables exist
                        required_vars = ['slowtemp', 'fasttemp', 'pressure', 'RH', 'dewpoint', 'RH_derived']
                        if all(var in conv_ds for var in required_vars):
                            # Calculate pre-correction statistics
                            slowtemp_orig = conv_ds['slowtemp'].values.copy()
                            fasttemp_data = conv_ds['fasttemp'].values
                            valid_mask = ~np.isnan(slowtemp_orig) & ~np.isnan(fasttemp_data)
                            if np.sum(valid_mask) > 0:
                                bias_before = np.mean(slowtemp_orig[valid_mask] - fasttemp_data[valid_mask])
                            else:
                                bias_before = np.nan

                            # Apply bias correction
                            conv_ds = apply_slowtemp_bias_correction(conv_ds, slope, intercept)

                            # Calculate post-correction statistics
                            slowtemp_corrected = conv_ds['slowtemp'].values
                            valid_mask = ~np.isnan(slowtemp_corrected) & ~np.isnan(fasttemp_data)
                            if np.sum(valid_mask) > 0:
                                bias_after = np.mean(slowtemp_corrected[valid_mask] - fasttemp_data[valid_mask])
                            else:
                                bias_after = np.nan

                            utils.log(f"  Applied correction: slope={slope:.6f}, intercept={intercept:.6f}")
                            utils.log(f"  Bias before correction: {bias_before:.4f}°C")
                            utils.log(f"  Bias after correction: {bias_after:.4f}°C")
                            utils.log(f"  Recomputed dewpoint and RH_derived from corrected slowtemp")

                            # If parsivel data exists, resample corrected data
                            if parsivel_combined_ds is not None:
                                try:
                                    PSD_datetimes = pips.get_PSD_datetimes(parsivel_combined_ds['VD_matrix'])
                                    sec_offset = PSD_datetimes[0].second
                                    DSD_interval = parsivel_combined_ds.DSD_interval
                                    interval_str = pips.get_interval_str(DSD_interval)
                                    offset_str = pips.get_interval_str(sec_offset)

                                    # Resample corrected variables to parsivel times
                                    if 'slowtemp' in parsivel_combined_ds:
                                        slowtemp_resampled = conv_ds['slowtemp'].resample(
                                            time=interval_str, label='right', closed='right', offset=offset_str
                                        ).mean()
                                        parsivel_combined_ds['slowtemp'] = slowtemp_resampled
                                        parsivel_combined_ds['slowtemp'].attrs.update(conv_ds['slowtemp'].attrs)

                                    if 'dewpoint' in parsivel_combined_ds:
                                        dewpoint_resampled = conv_ds['dewpoint'].resample(
                                            time=interval_str, label='right', closed='right', offset=offset_str
                                        ).mean()
                                        parsivel_combined_ds['dewpoint'] = dewpoint_resampled
                                        parsivel_combined_ds['dewpoint'].attrs.update(conv_ds['dewpoint'].attrs)

                                    if 'RH_derived' in parsivel_combined_ds:
                                        RH_derived_resampled = conv_ds['RH_derived'].resample(
                                            time=interval_str, label='right', closed='right', offset=offset_str
                                        ).mean()
                                        parsivel_combined_ds['RH_derived'] = RH_derived_resampled
                                        parsivel_combined_ds['RH_derived'].attrs.update(conv_ds['RH_derived'].attrs)

                                    utils.log(f"  Updated parsivel dataset with corrected thermodynamic variables")

                                except Exception as e:
                                    utils.log(f"Warning: Could not resample corrected variables for {PIPS_name}: {e}")
                        else:
                            missing_vars = [var for var in required_vars if var not in conv_ds]
                            utils.log(f"Warning: Cannot apply bias correction, missing variables: {missing_vars}")
                    else:
                        utils.log(f"  No correction needed for {PIPS_name} (slope=1.0, intercept=0.0)")
                else:
                    utils.log(f"Warning: {PIPS_name} not found in bias correction config")

            except Exception as e:
                utils.log(f"Error loading bias correction config: {e}")

    # =============================================================================================
    # GENERATE DIAGNOSTIC PLOTS (after all QC and trimming)
    # =============================================================================================
    # Generate compass/wind diagnostic plots if plotting was requested and original data was saved
    # This works for:
    #   1. Compass QC applied (with or without trimming)
    #   2. Only trimming applied (shows effect of removing fixed duration from start/end)
    #   3. Both compass QC and trimming applied
    if args.plot_compass_qc and 'compass_dir_original' in locals() and 'winddirabs_original' in locals():
        utils.log(f"Generating compass QC diagnostic plots for {PIPS_name}...")

        # Use the current state of conv_ds, which reflects all operations applied
        # (compass QC, trimming, or both)
        final_compass_dir = conv_ds['compass_dir']
        final_avg_compass_dir = circular_mean_deg(final_compass_dir.values)
        final_winddirabs = conv_ds['winddirabs']

        plot_compass_qc_diagnostics(compass_dir_original, winddirabs_original,
                                   final_compass_dir, final_avg_compass_dir,
                                   final_winddirabs, PIPS_name, deployment_name, plot_dir)

    # Generate slowtemp diagnostic plots if plotting was requested and original data was saved
    # This works for slowtemp QC applied (with or without trimming)
    if args.plot_slowtemp_qc and 'slowtemp_original' in locals() and 'fasttemp_original' in locals():
        utils.log(f"Generating slowtemp QC diagnostic plots for {PIPS_name}...")

        # Use the current state of conv_ds, which reflects all operations applied
        final_slowtemp = conv_ds['slowtemp']
        final_fasttemp = conv_ds['fasttemp']

        plot_slowtemp_qc_diagnostics(slowtemp_original, fasttemp_original,
                                    final_slowtemp, final_fasttemp,
                                    PIPS_name, deployment_name, plot_dir)

    # =============================================================================================
    # SAVE UPDATED DATASETS
    # =============================================================================================
    # Save conventional dataset if it was loaded and modified
    if conv_ds is not None:
        conv_file = conv_filelist[index]
        conv_output_file = conv_file + args.output_file_tag
        utils.log(f"Saving updated conventional data: {conv_output_file}")
        conv_ds.to_netcdf(conv_output_file)

    # Save parsivel dataset
    parsivel_combined_output_file = parsivel_combined_file + args.output_file_tag
    print("Dumping {}".format(parsivel_combined_output_file))
    parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
