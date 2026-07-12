"""
This script adds 3DEP elevation data to PIPS netCDF files as an additional global attribute.
It uses the latitude and longitude of the PIPS site to query the 3DEP API and retrieve the elevation
data.
"""

from __future__ import annotations

import argparse
import ast
import glob
import importlib.util
import os
from collections import OrderedDict

import requests
import xarray as xr

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIGS_DIR = os.path.join(PROJECT_ROOT, 'configs')


def parse_args():
    """Parse command-line arguments."""
    description = (
        'Add 3DEP elevation metadata to PIPS netCDF files for one or more case configurations'
    )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'case_config_path',
        nargs='+',
        metavar='<path/to/case/config/file.py>',
        help=(
            'Configuration file path or glob pattern. Bare filenames are resolved in configs/, '
            'for example ICECHIP_IOP1_2025_10s.py or configs/ICECHIP_IOP*_2025_10s.py.'
        ),
    )
    return parser.parse_args()


def resolve_config_patterns(config_patterns):
    """Resolve one or more config file paths or glob patterns to absolute config paths."""
    resolved_paths = OrderedDict()

    for pattern in config_patterns:
        if os.path.isabs(pattern):
            search_pattern = pattern
        elif pattern.startswith(('configs/', f'configs{os.sep}')) or os.path.dirname(pattern):
            search_pattern = os.path.join(PROJECT_ROOT, pattern)
        else:
            search_pattern = os.path.join(CONFIGS_DIR, pattern)

        matches = sorted(glob.glob(search_pattern))
        if not matches:
            msg = f'No configuration files matched pattern: {pattern}'
            raise FileNotFoundError(msg)

        for match in matches:
            resolved_paths[os.path.abspath(match)] = None

    return list(resolved_paths)


def import_all_from(module_path):
    """Load a Python config file from disk as a module."""
    spec = importlib.util.spec_from_file_location('mod', module_path)
    if spec is None or spec.loader is None:
        msg = f'Unable to load module from {module_path}'
        raise ImportError(msg)

    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def get_case_filepaths(config_path):
    """Load a case config and return the configured Parsivel and conventional netCDF files."""
    print(f'Case config file is {config_path}')  # noqa: T201
    try:
        config = import_all_from(config_path)
        print('Successfully imported case configuration parameters!')  # noqa: T201
    except Exception:
        msg = 'Unable to import case configuration parameters! Aborting!'
        raise RuntimeError(msg) from None

    pips_io_dict = config.PIPS_IO_dict
    pips_dir = pips_io_dict.get('PIPS_dir')
    parsivel_filenames = pips_io_dict.get('PIPS_filenames_nc', [])
    conv_filenames = pips_io_dict.get('conv_filenames_nc', [])

    if pips_dir is None:
        msg = f'PIPS_dir is not defined in {config_path}'
        raise ValueError(msg)

    return [os.path.join(pips_dir, filename) for filename in parsivel_filenames + conv_filenames]


def add_3dep_attribute(netcdf_path):
    """Update a netCDF file with the 3DEP elevation attribute derived from its location."""
    print(f'Processing {netcdf_path}...')  # noqa: T201
    ds = xr.load_dataset(netcdf_path)
    try:
        location = ds.attrs.get('location')
        if location is None:
            msg = f'location attribute is missing from {netcdf_path}'
            raise ValueError(msg)

        lat, lon, _ = ast.literal_eval(location)
        elev_3dep = epqs_3dep_elev_m(lat, lon)
        ds.attrs['elevation_3dep'] = elev_3dep
    finally:
        ds.close()

    ds.to_netcdf(netcdf_path, mode='w')

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
        msg = f"Unexpected EPQS response structure: {j}"
        raise RuntimeError(msg) from e


if __name__ == "__main__":
    args = parse_args()
    config_paths = resolve_config_patterns(args.case_config_path)

    for config_path in config_paths:
        netcdf_paths = get_case_filepaths(config_path)
        for netcdf_path in netcdf_paths:
            add_3dep_attribute(netcdf_path)