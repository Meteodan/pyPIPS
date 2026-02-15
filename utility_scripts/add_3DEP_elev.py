"""
This script adds 3DEP elevation data to PIPS netCDF files as an additional global attribute.
It uses the latitude and longitude of the PIPS site to query the 3DEP API and retrieve the elevation
data.
"""

import requests
import pandas as pd
import glob
import os
import xarray as xr

def epqs_3dep_elev_m(lat, lon):
    """
    USGS EPQS point elevation query (3DEP). Returns elevation in meters.

    Notes:
      - EPQS expects lon/lat as x/y.
      - EPQS elevations are interpolated from 3DEP elevation service.
    """
    url = "https://epqs.nationalmap.gov/v1/json"
    params = {
        "x": float(lon),
        "y": float(lat),
        "units": "Meters",
        "output": "json"
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    j = r.json()

    # EPQS response schema nests things; be defensive:
    # Commonly: j["value"] exists, but some responses use a nested structure.
    if "value" in j:
        return float(j["value"])

    # Fallback for documented structure (varies by EPQS version/output schema):
    try:
        return float(j["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"])
    except Exception as e:
        raise RuntimeError(f"Unexpected EPQS response structure: {j}") from e

if __name__ == "__main__":
    # Add 3DEP elevation data to PIPS netCDF files in a series of subdirectories.
    basedir = "/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data"
    IOP_dirs = glob.glob(os.path.join(basedir, "IOP*"))
    for IOP_dir in IOP_dirs:
        netcdf_dir = os.path.join(IOP_dir, 'netcdf')
        conv_filepaths = glob.glob(netcdf_dir + '/conv*nc')
        parsivel_filepaths = glob.glob(netcdf_dir + '/parsivel*nc')
        for conv_filepath in conv_filepaths:
            print(f"Processing {conv_filepath}...")
            ds = xr.load_dataset(conv_filepath)
            ds.close()  # Close the dataset to avoid file locking issues
            lat, lon, elev = eval(ds.location)
            elev_3dep = epqs_3dep_elev_m(lat, lon)
            ds.attrs['elevation_3dep'] = elev_3dep
            ds.to_netcdf(conv_filepath, mode='w')
        for parsivel_filepath in parsivel_filepaths:
            print(f"Processing {parsivel_filepath}...")
            ds = xr.load_dataset(parsivel_filepath)
            ds.close()  # Close the dataset to avoid file locking issues
            lat, lon, elev = eval(ds.location)
            elev_3dep = epqs_3dep_elev_m(lat, lon)
            ds.attrs['elevation_3dep'] = elev_3dep
            ds.to_netcdf(parsivel_filepath, mode='w')