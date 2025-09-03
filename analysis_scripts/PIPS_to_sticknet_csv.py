# PIPS_to_sticknet_csv.py
#
# This script reads in conventional PIPS netCDF data for a given case and outputs corresponding
# CSV files in the same format as a Texas Tech (TTU) Sticknet file.

from __future__ import annotations

import argparse
import os

import pandas as pd
import xarray as xr
from pyPIPS import utils


def convert_pips_to_sticknet_csv(conv_ds, output_filepath, parsivel_ds=None,
                                 parsivel_interp='none', drop_nan=False):
    """
    Convert PIPS conventional dataset to TTU Sticknet CSV format.

    Parameters:
    -----------
    conv_ds : xarray.Dataset
        The conventional PIPS dataset containing the meteorological data.
    output_filepath : str
        The path where the CSV file will be saved.
    parsivel_ds : xarray.Dataset, optional
        The parsivel_combined PIPS dataset containing additional DSD-derived variables.
    parsivel_interp : str, optional
        Method for handling parsivel data at different time intervals:
        'interpolate' (linear interpolation), 'forward_fill' (forward fill),
        'none' (leave as NaN, default)
    drop_nan : bool, optional
        Whether to drop rows with any NaN values (default: False)

    Returns:
    --------
    None
    """
    # Extract the required variables from the dataset
    time = conv_ds['time']

    # Find temperature variable
    temperature_vars = ['fasttemp', 'slowtemp']
    temperature = None
    for var in temperature_vars:
        if var in conv_ds.variables:
            temperature = conv_ds[var]
            break

    if temperature is None:
        msg = "No temperature variable found in dataset"
        raise ValueError(msg)

    # Convert to Kelvin if needed
    if hasattr(temperature, 'units') and 'celsius' in temperature.units.lower():
        temperature += 273.15

    # Relative humidity
    rh_vars = ['RH_derived', 'RH']
    relative_humidity = None
    for var in rh_vars:
        if var in conv_ds.variables:
            relative_humidity = conv_ds[var]
            break

    if relative_humidity is None:
        msg = "No relative humidity variable found in dataset"
        raise ValueError(msg)

    # Pressure (convert from Pa to hPa if needed)
    pressure_vars = ['pressure']
    pressure = None
    for var in pressure_vars:
        if var in conv_ds.variables:
            pressure = conv_ds[var]
            break

    if pressure is None:
        msg = "No pressure variable found in dataset"
        raise ValueError(msg)

    # Convert from Pa to hPa if needed
    if (hasattr(pressure, 'units') and 'pa' in pressure.units.lower()
        and 'hpa' not in pressure.units.lower()):
        pressure /= 100.0

    # Wind speed
    ws_vars = ['windspd']
    wind_speed = None
    for var in ws_vars:
        if var in conv_ds.variables:
            wind_speed = conv_ds[var]
            break

    if wind_speed is None:
        msg = "No wind speed variable found in dataset"
        raise ValueError(msg)

    # Wind direction
    wd_vars = ['winddirabs']
    wind_direction = None
    for var in wd_vars:
        if var in conv_ds.variables:
            wind_direction = conv_ds[var]
            break

    if wind_direction is None:
        msg = "No wind direction variable found in dataset"
        raise ValueError(msg)

    # Create the base DataFrame with conventional variables
    df_dict = {
        'Time': pd.to_datetime(time.to_numpy()).strftime('%Y-%m-%d %H:%M:%S'),
        'T': temperature.to_numpy(),
        'RH': relative_humidity.to_numpy(),
        'P': pressure.to_numpy(),
        'WS': wind_speed.to_numpy(),
        'WD': wind_direction.to_numpy()
    }

    # Add parsivel variables if parsivel dataset is provided
    if parsivel_ds is not None:
        # Interpolate parsivel data to conventional time grid if needed
        conv_time = pd.to_datetime(time.to_numpy())
        parsivel_time = pd.to_datetime(parsivel_ds['time'].to_numpy())

        # Mean drop diameter (Dm43_roqc)
        if 'Dm43_roqc' in parsivel_ds.variables:
            parsivel_df = pd.DataFrame({
                'time': parsivel_time,
                'Dm43': parsivel_ds['Dm43_roqc'].to_numpy() * 1000.  # Convert to mm
            }).set_index('time')

            if parsivel_interp == 'interpolate':
                # Linear interpolation
                interp_df = parsivel_df.reindex(conv_time).interpolate(method='linear')
                df_dict['Dm43'] = interp_df['Dm43'].to_numpy()
            elif parsivel_interp == 'forward_fill':
                # Forward fill
                interp_df = parsivel_df.reindex(conv_time).ffill()
                df_dict['Dm43'] = interp_df['Dm43'].to_numpy()
            else:  # 'none'
                # Leave as NaN for missing times
                interp_df = parsivel_df.reindex(conv_time)
                df_dict['Dm43'] = interp_df['Dm43'].to_numpy()

        # Rain rate (rainrate_derived_roqc)
        if 'rainrate_derived_roqc' in parsivel_ds.variables:
            parsivel_df = pd.DataFrame({
                'time': parsivel_time,
                'RR': parsivel_ds['rainrate_derived_roqc'].to_numpy()
            }).set_index('time')

            if parsivel_interp == 'interpolate':
                interp_df = parsivel_df.reindex(conv_time).interpolate(method='linear')
                df_dict['RR'] = interp_df['RR'].to_numpy()
            elif parsivel_interp == 'forward_fill':
                interp_df = parsivel_df.reindex(conv_time).ffill()
                df_dict['RR'] = interp_df['RR'].to_numpy()
            else:  # 'none'
                interp_df = parsivel_df.reindex(conv_time)
                df_dict['RR'] = interp_df['RR'].to_numpy()

        # Reflectivity (reflectivity_derived_roqc)
        if 'reflectivity_derived_roqc' in parsivel_ds.variables:
            parsivel_df = pd.DataFrame({
                'time': parsivel_time,
                'Z': parsivel_ds['reflectivity_derived_roqc'].to_numpy()
            }).set_index('time')

            if parsivel_interp == 'interpolate':
                interp_df = parsivel_df.reindex(conv_time).interpolate(method='linear')
                df_dict['Z'] = interp_df['Z'].to_numpy()
            elif parsivel_interp == 'forward_fill':
                interp_df = parsivel_df.reindex(conv_time).ffill()
                df_dict['Z'] = interp_df['Z'].to_numpy()
            else:  # 'none'
                interp_df = parsivel_df.reindex(conv_time)
                df_dict['Z'] = interp_df['Z'].to_numpy()

    # Create a pandas DataFrame
    df = pd.DataFrame(df_dict)

    # Remove any rows with NaN values only if drop_nan is True
    if drop_nan:
        df = df.dropna()

    # Save to CSV with the specified header
    df.to_csv(output_filepath, index=False)
    utils.log(f"Saved Sticknet CSV file: {output_filepath}")


def convert_parsivel_to_sticknet_csv(parsivel_ds, output_filepath):
    """
    Convert PIPS parsivel_combined dataset to TTU Sticknet CSV format.

    Parameters:
    -----------
    parsivel_ds : xarray.Dataset
        The parsivel_combined PIPS dataset containing DSD-derived variables.
    output_filepath : str
        The path where the CSV file will be saved.

    Returns:
    --------
    None
    """
    # Extract the required variables from the dataset
    time = parsivel_ds['time']

    # Create the DataFrame with parsivel-only variables
    df_dict = {
        'Time': pd.to_datetime(time.to_numpy()).strftime('%Y-%m-%d %H:%M:%S'),
    }

    # Add parsivel variables
    # Mean drop diameter (Dm43_roqc)
    if 'Dm43_roqc' in parsivel_ds.variables:
        df_dict['Dm43'] = parsivel_ds['Dm43_roqc'].to_numpy()
    else:
        utils.log("Warning: Dm43_roqc not found in parsivel dataset")

    # Rain rate (rainrate_derived_roqc)
    if 'rainrate_derived_roqc' in parsivel_ds.variables:
        df_dict['RR'] = parsivel_ds['rainrate_derived_roqc'].to_numpy()
    else:
        utils.log("Warning: rainrate_derived_roqc not found in parsivel dataset")

    # Reflectivity (reflectivity_derived_roqc)
    if 'reflectivity_derived_roqc' in parsivel_ds.variables:
        df_dict['Z'] = parsivel_ds['reflectivity_derived_roqc'].to_numpy()
    else:
        utils.log("Warning: reflectivity_derived_roqc not found in parsivel dataset")

    # Create a pandas DataFrame
    df = pd.DataFrame(df_dict)

    # Remove any rows with NaN values
    df = df.dropna()

    # Save to CSV with the specified header
    df.to_csv(output_filepath, index=False)
    utils.log(f"Saved Parsivel Sticknet CSV file: {output_filepath}")


def main():
    """Main function to process PIPS data and convert to Sticknet CSV format."""

    # Parse the command line options
    description = ("Reads in conventional PIPS netCDF data for a given case and outputs "
                  "corresponding CSV files in TTU Sticknet format")
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                        help='The path to the case configuration file')
    parser.add_argument('--output-dir', dest='output_dir', type=str, default=None,
                        help='Output directory for CSV files (defaults to PIPS_dir)')
    parser.add_argument('--parsivel-interp', dest='parsivel_interp', type=str,
                        default='none', choices=['interpolate', 'forward_fill', 'none'],
                        help='Method for handling parsivel data at different time intervals: '
                             'interpolate (linear interpolation), forward_fill (forward fill), '
                             'none (leave as NaN, default)')
    parser.add_argument('--drop-nan', dest='drop_nan', action='store_true',
                        help='Drop rows with any NaN values (default: False)')

    args = parser.parse_args()

    # Dynamically import the case configuration file
    utils.log(f"Case config file is {args.case_config_path}")
    try:
        config = utils.import_all_from(args.case_config_path)
        utils.log("Successfully imported case configuration parameters!")
    except Exception:
        utils.fatal("Unable to import case configuration parameters! Aborting!")

    # Extract needed lists and variables from PIPS_IO_dict configuration dictionary
    dataset_name = config.PIPS_IO_dict.get('dataset_name', None)
    PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
    PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
    conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', None)
    parsivel_filenames_nc = config.PIPS_IO_dict.get('PIPS_filenames_nc', None)

    # Set output directory
    output_dir = args.output_dir if args.output_dir is not None else PIPS_dir

    # Get a list of the conventional and parsivel_combined netCDF data files
    conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames_nc]
    parsivel_filelist = ([os.path.join(PIPS_dir, pf) for pf in parsivel_filenames_nc]
                        if parsivel_filenames_nc else None)

    # Process each PIPS conventional file
    PIPS_locs = []
    for i, PIPS_name in enumerate(PIPS_names):
        conv_filepath = conv_filelist[i]
        parsivel_filepath = parsivel_filelist[i] if parsivel_filelist else None
        utils.log(f"Processing {PIPS_name}: {conv_filepath}")

        # Load the conventional dataset
        try:
            conv_ds = xr.load_dataset(conv_filepath)
        except Exception as e:
            utils.log(f"Error loading {conv_filepath}: {e}")
            continue

        # Load the parsivel dataset if available
        parsivel_ds = None
        if parsivel_filepath:
            try:
                parsivel_ds = xr.load_dataset(parsivel_filepath)
                utils.log(f"Also loading parsivel data: {parsivel_filepath}")
            except Exception as e:
                utils.log(f"Error loading {parsivel_filepath}: {e}")
                # Continue without parsivel data

        # Get PIPS location information

        if 'location' in conv_ds.attrs:
            PIPS_loc = eval(conv_ds.attrs['location'])
            utils.log(f"PIPS location for {PIPS_name}: {PIPS_loc}")
        else:
            PIPS_loc = None
            utils.log(f"No location information found for {PIPS_name}")
        PIPS_locs.append(PIPS_loc)

        # Create output filename
        output_filename = f'{PIPS_name}_{dataset_name}_sticknet.txt'
        output_filepath = os.path.join(output_dir, output_filename)

        # Convert to Sticknet CSV format
        try:
            convert_pips_to_sticknet_csv(conv_ds, output_filepath, parsivel_ds,
                                       args.parsivel_interp, args.drop_nan)
            utils.log(f"Successfully converted {PIPS_name} to Sticknet format")
        except Exception as e:
            utils.log(f"Error converting {PIPS_name}: {e}")
            continue

        # Also create a parsivel-only CSV file if parsivel data is available
        if parsivel_ds is not None:
            parsivel_output_filename = f'{PIPS_name}_{dataset_name}_parsivel_sticknet.txt'
            parsivel_output_filepath = os.path.join(output_dir, parsivel_output_filename)

            try:
                convert_parsivel_to_sticknet_csv(parsivel_ds, parsivel_output_filepath)
                utils.log(f"Successfully converted {PIPS_name} parsivel data to Sticknet format")
            except Exception as e:
                utils.log(f"Error converting {PIPS_name} parsivel data: {e}")
                continue

    # Save PIPS locations to a file
    locs_filepath = os.path.join(output_dir, f'{dataset_name}_PIPS_locations.csv')
    with open(locs_filepath, 'w', encoding='utf-8') as locs_file:
        # Write the header. Should look like: ID,Latitude,Longitude,Elevation,Array_Type
        locs_file.write("ID,Latitude,Longitude,Elevation,Array_Type\n")
        # Write each PIPS location, matching the expected format. Array_Type is assumed to be "Fine"
        for i, loc in enumerate(PIPS_locs):
            if loc is not None:
                locs_file.write(f"{PIPS_names[i]},{loc[0]:.6f},{loc[1]:.6f},{loc[2]:.2f},Fine\n")
            else:
                locs_file.write(f"{PIPS_names[i]},No location information\n")
    utils.log(f"Saved PIPS locations to {locs_filepath}")
    utils.log("Conversion complete!")


if __name__ == '__main__':
    main()
