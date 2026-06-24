"""
This script checks the time attributes of PIPS netCDF files to ensure they are correctly set.
It reads the 'starting_time' and 'ending_time' attributes and compares them to the actual time range
of the data in the netCDF files, correcting any discrepancies if necessary.
"""

import requests
import pandas as pd
import glob
import os
import xarray as xr
import pyPIPS.PIPS as pips


if __name__ == "__main__":
    # Read and check time attributes of PIPS netCDF files in a series of subdirectories.
    dry_run = True  # Set to True to only print discrepancies without modifying files
    basedir = "/Users/dawson29/Projects/ICECHIP/obsdata/PIPS_data"
    IOP_dirs = glob.glob(os.path.join(basedir, "IOP*"))
    for IOP_dir in IOP_dirs:
        netcdf_dir = os.path.join(IOP_dir, 'netcdf')
        conv_filepaths = glob.glob(netcdf_dir + '/conv*nc')
        parsivel_filepaths = glob.glob(netcdf_dir + '/parsivel*nc')
        for conv_filepath in conv_filepaths:
            print(f"Processing {conv_filepath}...")
            ds = xr.load_dataset(conv_filepath)
            ds.close()  # Close the dataset to avoid file locking issues
            # Check and correct time attributes
            starting_time_attr = ds.attrs.get('starting_time')
            ending_time_attr = ds.attrs.get('ending_time')
            print(f"  Original starting_time: {starting_time_attr}, ending_time: {ending_time_attr}")
            data_times = pips.get_datetimes(ds)
            if data_times is not None:
                if starting_time_attr != data_times[0].strftime('%Y%m%d%H%M%S'):
                    ds.attrs['starting_time'] = data_times[0].strftime('%Y%m%d%H%M%S')
                    print("Starting time attribute was incorrect and has been updated.")
                if ending_time_attr != data_times[-1].strftime('%Y%m%d%H%M%S'):
                    ds.attrs['ending_time'] = data_times[-1].strftime('%Y%m%d%H%M%S')
                    print("Ending time attribute was incorrect and has been updated.")
                print(f"  Updated starting_time: {ds.attrs['starting_time']}, ending_time: {ds.attrs['ending_time']}")
                if not dry_run:
                    ds.to_netcdf(conv_filepath, mode='w')
        for parsivel_filepath in parsivel_filepaths:
            print(f"Processing {parsivel_filepath}...")
            ds = xr.load_dataset(parsivel_filepath)
            ds.close()  # Close the dataset to avoid file locking issues
            # Check and correct time attributes
            starting_time_attr = ds.attrs.get('starting_time')
            ending_time_attr = ds.attrs.get('ending_time')
            print(f"  Original starting_time: {starting_time_attr}, ending_time: {ending_time_attr}")
            data_times = pips.get_datetimes(ds)
            if data_times is not None:
                if starting_time_attr != data_times[0].strftime('%Y%m%d%H%M%S'):
                    ds.attrs['starting_time'] = data_times[0].strftime('%Y%m%d%H%M%S')
                    print("Starting time attribute was incorrect and has been updated.")
                if ending_time_attr != data_times[-1].strftime('%Y%m%d%H%M%S'):
                    ds.attrs['ending_time'] = data_times[-1].strftime('%Y%m%d%H%M%S')
                    print("Ending time attribute was incorrect and has been updated.")
                print(f"  Updated starting_time: {ds.attrs['starting_time']}, ending_time: {ds.attrs['ending_time']}")
                if not dry_run:
                    ds.to_netcdf(parsivel_filepath, mode='w')