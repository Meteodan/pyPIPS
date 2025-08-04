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


def convert_pips_to_sticknet_csv(conv_ds, output_filepath):
    """
    Convert PIPS conventional dataset to TTU Sticknet CSV format.

    Parameters:
    -----------
    conv_ds : xarray.Dataset
        The conventional PIPS dataset containing the meteorological data.
    output_filepath : str
        The path where the CSV file will be saved.

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

    # Create a pandas DataFrame
    df = pd.DataFrame({
        'Time': pd.to_datetime(time.to_numpy()).strftime('%Y-%m-%d %H:%M:%S'),
        'T': temperature.to_numpy(),
        'RH': relative_humidity.to_numpy(),
        'P': pressure.to_numpy(),
        'WS': wind_speed.to_numpy(),
        'WD': wind_direction.to_numpy()
    })

    # Remove any rows with NaN values
    df = df.dropna()

    # Save to CSV with the specified header
    df.to_csv(output_filepath, index=False)
    utils.log(f"Saved Sticknet CSV file: {output_filepath}")


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

    # Set output directory
    output_dir = args.output_dir if args.output_dir is not None else PIPS_dir

    # Get a list of the conventional netCDF data files
    conv_filelist = [os.path.join(PIPS_dir, cf) for cf in conv_filenames_nc]

    # Process each PIPS conventional file
    for PIPS_name, conv_filepath in zip(PIPS_names, conv_filelist):
        utils.log(f"Processing {PIPS_name}: {conv_filepath}")

        # Load the conventional dataset
        try:
            conv_ds = xr.load_dataset(conv_filepath)
        except Exception as e:
            utils.log(f"Error loading {conv_filepath}: {e}")
            continue

        # Create output filename
        output_filename = f'{dataset_name}_{PIPS_name}_sticknet.csv'
        output_filepath = os.path.join(output_dir, output_filename)

        # Convert to Sticknet CSV format
        try:
            convert_pips_to_sticknet_csv(conv_ds, output_filepath)
            utils.log(f"Successfully converted {PIPS_name} to Sticknet format")
        except Exception as e:
            utils.log(f"Error converting {PIPS_name}: {e}")
            continue

    utils.log("Conversion complete!")


if __name__ == '__main__':
    main()
