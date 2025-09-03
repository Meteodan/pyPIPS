# pyPIPS_meteograms.py
#
# This script calculates radar retrievals from PIPS DSD data (netCDF version)
from __future__ import annotations

import argparse
import os

# from datetime import datetime, timedelta
# import numpy as np
# import pandas as pd
import numpy as np
import xarray as xr

# import pyPIPS.PIPS as pips
# import pyPIPS.parsivel_qc as pqc
# import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd

# import matplotlib.ticker as ticker
# import matplotlib.dates as dates
# import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp

# import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.polarimetric as dp
import pyPIPS.radarmodule as radar
from pyPIPS import utils
import contextlib

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Calculates radar retrievals from PIPS and collocated radar data (netCDF version)"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--QC-tag', dest='QC_tag', default=None,
                    help='QC tag for DSD variables in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--calc-for-SATP', action='store_true', dest='calc_for_SATP',
                    help='calculate for the SATP-filtered dataset')
parser.add_argument('--coefficients', nargs=3, metavar=('c1', 'c2', 'c3'), default=None, type=float,
                    dest='coefficients',
                    help='coefficients of mu-lambda polynomial (in increasing order of exponent')
parser.add_argument('--retrieval-tag', dest='retrieval_tag', default='',
                    help='nametag for the name of the mu-lambda relation (e.g. SATP, C08, Z01)')
parser.add_argument('--output-tag', dest='output_tag', default='',
                    help='tag for output nc files to distinguish from original if desired')
parser.add_argument('--update-CG-coeff-attrs', dest='update_CG_coeff_attrs', action='store_true',
                    default=False,
                    help='Write CG coefficients as attributes back to PIPS nc files?')
parser.add_argument('--use-filtered-fields', dest='use_filtered_fields', default=False,
                    action='store_true',
                    help='Whether to use previously filtered dBZ and ZDR fields for the retrieval')
parser.add_argument('--time-dim', dest='time_dim', default='time',
                    help='Name of the time dimension in the PIPS netCDF file')

args = parser.parse_args()
QC_tag = '' if not args.QC_tag else f'_{args.QC_tag}'

# Dynamically import the case configuration file
utils.log(f"Case config file is {args.case_config_path}")
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
start_times = config.PIPS_IO_dict.get('start_times', [None] * len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None] * len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None] * len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

# Extract needed lists and variables from the radar_dict configuration dictionary
comp_radar = config.radar_config_dict.get('comp_radar', False)
radar_name = config.radar_config_dict.get('radar_name', None)
radar_type = config.radar_config_dict.get('radar_type', 'NEXRAD')
radar_dir = config.radar_config_dict.get('radar_dir', None)
scatt_dir = config.radar_config_dict.get('scatt_dir', None)
wavelength = config.radar_config_dict.get('wavelength', 10.7)

scatt_file = os.path.join(scatt_dir, 'SCTT_RAIN_fw100.dat')

if not args.calc_for_SATP:
    # Get a list of the combined parsivel netCDF data files that are present in the PIPS directory
    parsivel_combined_filelist = [os.path.join(PIPS_dir, pcf) for pcf
                                  in parsivel_combined_filenames]
else:
    # Just read in the single combined SATP dataset
    parsivel_combined_filename = f'ND_avg_{dataset_name}{QC_tag}_{int(requested_interval):d}s.nc'
    parsivel_combined_filelist = [os.path.join(PIPS_dir, parsivel_combined_filename)]

for index, parsivel_combined_file in enumerate(parsivel_combined_filelist):
    print(f"Reading {parsivel_combined_file}")  # noqa: T201
    parsivel_combined_ds = xr.load_dataset(parsivel_combined_file)

    if index == 0:
        if args.coefficients:
            mu_lambda_coeff = args.coefficients
            retrieval_tag = args.retrieval_tag
        else:
            if args.retrieval_tag in {'C08', 'Z01'}:
                retrieval_tag = args.retrieval_tag
            else:
                retrieval_tag = f'{args.retrieval_tag}{QC_tag}'
            mu_lambda_coeff = parsivel_combined_ds.attrs[f'CG_coeff_{retrieval_tag}']

    if not args.calc_for_SATP:
        DSD_interval = parsivel_combined_ds.DSD_interval
        PIPS_name = parsivel_combined_ds.probe_name
        deployment_name = parsivel_combined_ds.deployment_name
        ND = parsivel_combined_ds[f'ND{QC_tag}']
        coord_to_combine = 'time'
    else:
        ND = parsivel_combined_ds[f'SATP_ND{dataset_name}']
        ND = ND.rename({'D0_RR_level_0': 'D0_idx', 'D0_RR_level_1': 'RR_idx'})
        DSD_interval = ND.DSD_interval
        ND = pipsio.reconstruct_MultiIndex(ND, ['D0_idx', 'RR_idx'], 'D0_RR')
        coord_to_combine = 'D0_RR'

    ND = ND.where(ND > 0.)
    D = parsivel_combined_ds.coords['diameter']
    dD = (parsivel_combined_ds.coords['max_diameter'] - parsivel_combined_ds.coords['min_diameter'])
    dualpol_dict = dp.calpolrain(wavelength, scatt_file, ND, dD)

    ZH = dualpol_dict['REF']
    ZDR = dualpol_dict['ZDR']
    fa2 = dualpol_dict['fa2']
    fb2 = dualpol_dict['fb2']

    retr_dict = dsd.retrieval_Cao_xr(ZH, ZDR, ND, D, dD, fa2, fb2, wavelength, mu_lambda_coeff,
                                     retrieval_tag=f'{args.retrieval_tag}{QC_tag}')

    retr_ds = xr.Dataset(retr_dict)
    parsivel_combined_ds.update(retr_ds)
    if args.calc_for_SATP:
        # Save attrs (may not need to do this)
        attrs = parsivel_combined_ds.attrs
        parsivel_combined_ds = parsivel_combined_ds.reset_index('D0_RR')
        parsivel_combined_ds.attrs = attrs
    # parsivel_combined_ds.attrs['CG_coeff_{}'.format(args.retrieval_tag)] = mu_lambda_coeff
    parsivel_combined_ds.attrs['retrieval_wavelength'] = wavelength

    # Save the CG coefficients as attributes if desired
    if args.update_CG_coeff_attrs:
        CG_attr_name = f'CG_coeff_{retrieval_tag}'
        parsivel_combined_ds.attrs[CG_attr_name] = mu_lambda_coeff

    # Now do the same for the radar data interpolated to the PIPS. This will overwrite the
    # existing interpolated retrieved fields in the PIPS file. Note that this is a good thing
    # because the interpolated retrieved fields are in some cases (especially N0) not very good.
    # So, it is better to first interpolate the radar ZH and ZDR to the PIPS location and *then*
    # do the retrieval. TODO: disable interpolation of any retrieval fields on the radar sweeps
    # to the PIPS location in radar_to_PIPS.py to avoid confusion.

    if comp_radar:
        radar_fields_at_PIPS_da = parsivel_combined_ds[f'{radar_name}_at_PIPS']
        fname_tag = QC_tag + f'_{radar_name}'
        # print(radar_fields_at_PIPS_da)
        # At least add the reflectivity
        # First get list of radar fields in DataArray
        dim_name = f'fields_{radar_name}'
        radar_fields = radar_fields_at_PIPS_da.coords[dim_name].to_numpy()
        num_times = radar_fields_at_PIPS_da.sizes[args.time_dim]
        # Then find which one corresponds to reflectivity
        # ref_name = next((fname for fname in radar_fields if fname in radar.REF_aliases))
        ref_name = radar.find_radar_field_name(radar_fields, radar.REF_aliases)
        zdr_name = radar.find_radar_field_name(radar_fields, radar.ZDR_aliases)
        if args.use_filtered_fields:
            ref_name += '_filtered'

        ZH_rad = radar_fields_at_PIPS_da.loc[{dim_name: ref_name}]
        ZDR_rad = radar_fields_at_PIPS_da.loc[{dim_name: zdr_name}]

        retr_dict = dsd.retrieval_Cao_xr(ZH_rad, ZDR_rad, ND, D, dD, fa2, fb2, wavelength,
                                         mu_lambda_coeff, retrieval_tag=f'{args.retrieval_tag}')

        # This is tricky, because the new retrieval fields will change the dimension size and the
        # coordinate names for the radar fields, we have to create a new DataArray from scratch
        # and then replace the old one in the dataset. This is a bit of a PITA TBH.

        retr_names = ['RR', 'D0', 'mu', 'lamda', 'N0', 'Nt', 'W', 'sigma', 'Dm43']
        retr_keys = [f'{retr_name}_{retrieval_tag}' for retr_name in retr_names]
        radar_fields = [field_name for field_name in radar_fields if field_name not in retr_keys]
        new_keys = list(radar_fields) + retr_keys
        new_len_fields = len(new_keys)

        new_radar_fields_at_PIPS_da = \
            xr.DataArray(np.empty((num_times, new_len_fields)),
                         coords={
                             args.time_dim: radar_fields_at_PIPS_da[args.time_dim],
                             dim_name: new_keys,
                             },
                         dims=[args.time_dim, dim_name],
                         attrs={
                             'radar_name': radar_name,
                             'PIPS_name': PIPS_name,
                             })
        # Copy the PIPS x and y coordinates (relative to the radar) to the new DataArray if they
        # exist
        with contextlib.suppress(AttributeError):
            new_radar_fields_at_PIPS_da.attrs['PIPS_x'] = radar_fields_at_PIPS_da.PIPS_x
            new_radar_fields_at_PIPS_da.attrs['PIPS_y'] = radar_fields_at_PIPS_da.PIPS_y
        # Copy the old radar fields to the new DataArray
        for field_name in radar_fields:
            new_radar_fields_at_PIPS_da.loc[{dim_name: field_name}] = \
                radar_fields_at_PIPS_da.loc[{dim_name: field_name}]
        # Add the new retrieval fields
        for retr_name, retr_key in zip(retr_names, retr_keys):
            new_radar_fields_at_PIPS_da.loc[{dim_name: retr_key}] = \
                retr_dict[f'{retr_name}_retr_{retrieval_tag}']

        # For some reason I chose to name the Dm field differently in the radar and PIPS files
        # radar_fields_at_PIPS_da.loc[{dim_name: f'Dm_{retrieval_tag}'}] = \
        #     retr_dict[f'Dm43_retr_{retrieval_tag}']

        # Get rid of the old radar fields at PIPS
        if f'{radar_name}_at_PIPS' in parsivel_combined_ds:
            parsivel_combined_ds = parsivel_combined_ds.drop_vars(f'{radar_name}_at_PIPS',
                                                                  errors='ignore')
            parsivel_combined_ds = parsivel_combined_ds.drop_vars(f'fields_{radar_name}',
                                                                  errors='ignore')
        if f'fields_{radar_name}' in parsivel_combined_ds.dims:
            parsivel_combined_ds = parsivel_combined_ds.drop_dims([f'fields_{radar_name}'],
                                                                  errors='ignore')

        # Save the old and new retrieved radar fields back to the PIPS file
        parsivel_combined_ds[f'{radar_name}_at_PIPS'] = new_radar_fields_at_PIPS_da

    if args.output_tag:
        parsivel_combined_output_file = \
            parsivel_combined_file.replace(".nc", f"_{args.output_tag}.nc")
        print(f"Dumping {parsivel_combined_output_file}")  # noqa: T201
        parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
    else:
        parsivel_combined_output_file = parsivel_combined_file
        print(f"Dumping {parsivel_combined_output_file}")  # noqa: T201
        parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
