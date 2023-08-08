# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import argparse
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp
import logging

# TODO: Figure out how to turn off the netCDF logging messages! Below doesn't work....
# logging.getLogger("xarray").setLevel(logging.ERROR)

qc_attr_names = ['strongwindQC', 'splashingQC', 'marginQC', 'rainfallQC', 'rainonlyQC',
                 'hailonlyQC', 'graupelonlyQC', 'basicQC']

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Parse the command line options
description = "Reads PIPS comma-delimited text data files and converts to netCDF"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('case_config_path', metavar='<path/to/case/config/file.py>',
                    help='The path to the case configuration file')
parser.add_argument('--output-dir', dest='output_dir', default=None,
                    help='Directory to put netCDF files. Defaults to input dir in config file')
parser.add_argument('--check-order', dest='check_order', action='store_true',
                    help='Check if there are any times out of order')
parser.add_argument('--sort-times', dest='sort_times', action='store_true',
                    help='Sort dataset by time')
# parser.add_argument('--apply-qc', dest='apply_qc', action='store_true', help='apply QC to DSDs?')
# parser.add_argument('--qc-tag', dest='qc_tag', default=None,
#                     help='name of QC tag for DSD variables')
args = parser.parse_args()
# apply_qc = args.apply_qc
# qc_tag = args.qc_tag
# if qc_tag is None:
#     apply_qc = False
# if not apply_qc:
#     qc_tag = None

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
deployment_names = config.PIPS_IO_dict.get('deployment_names', None)
PIPS_dir = config.PIPS_IO_dict.get('PIPS_dir', None)
input_txt_dir = config.PIPS_IO_dict.get('input_txt_dir', PIPS_dir)
plot_dir = config.PIPS_IO_dict.get('plot_dir', None)
PIPS_types = config.PIPS_IO_dict.get('PIPS_types', None)
PIPS_names = config.PIPS_IO_dict.get('PIPS_names', None)
PIPS_filenames = config.PIPS_IO_dict.get('PIPS_filenames', None)
start_times = config.PIPS_IO_dict.get('start_times', [None]*len(PIPS_names))
end_times = config.PIPS_IO_dict.get('end_times', [None]*len(PIPS_names))
geo_locs = config.PIPS_IO_dict.get('geo_locs', [None]*len(PIPS_names))
requested_interval = config.PIPS_IO_dict.get('requested_interval', 10.)

if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = PIPS_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
# Read in the PIPS data for the deployment
conv_df_list = []
conv_resampled_df_list = []
parsivel_df_list = []
vd_matrix_da_list = []

for index, PIPS_filename, PIPS_name, start_time, end_time, geo_loc, ptype, deployment_name in \
        zip(range(0, len(PIPS_filenames)), PIPS_filenames, PIPS_names, start_times, end_times,
            geo_locs, PIPS_types, deployment_names):

    tripips = (ptype == 'TriPIPS')
    PIPS_data_file_path = os.path.join(input_txt_dir, PIPS_filename)
    conv_df, parsivel_df, vd_matrix_da = pipsio.read_PIPS(PIPS_data_file_path,
                                                          start_timestamp=start_time,
                                                          end_timestamp=end_time, tripips=tripips,
                                                          sort=args.sort_times,
                                                          check_order=args.check_order)

    # We need the disdrometer locations. If they aren't supplied in the input control file, find
    # them from the GPS data
    if not geo_loc:
        geo_locs[index] = pipsio.get_PIPS_loc(conv_df['GPS_status'], conv_df['GPS_lat'],
                                              conv_df['GPS_lon'], conv_df['GPS_alt'])
    print("Lat/Lon/alt of {}: {}".format(PIPS_name, str(geo_locs[index])))

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Only do the following if we have parsivel data. Need to clean this up to make it more
    # graceful.

    if parsivel_df is not None and vd_matrix_da is not None:

        # if apply_qc:
        #     # Do some QC on the V-D matrix. This will make a copy of the raw matrix. The netCDF 
        #     # file
        #     # will contain both
        #     qc_attr_dict = {}
        #     for qc_attr_name in qc_attr_names:
        #         qc_attr_dict[qc_attr_name] = config.PIPS_qc_dict[qc_attr_name]

        #     if qc_attr_dict['basicQC']:
        #         qc_attr_dict['strongwindQC'] = True
        #         qc_attr_dict['splashingQC'] = True
        #         qc_attr_dict['marginQC'] = True

        #     vd_matrix_qc_da = vd_matrix_da.copy()
        #     if qc_attr_dict['strongwindQC']:
        #         vd_matrix_qc_da = pqc.strongwindQC(vd_matrix_qc_da)
        #     if qc_attr_dict['splashingQC']:
        #         vd_matrix_qc_da = pqc.splashingQC(vd_matrix_qc_da)
        #     if qc_attr_dict['marginQC']:
        #         vd_matrix_qc_da = pqc.marginQC(vd_matrix_qc_da)
        #     if qc_attr_dict['rainfallQC']:
        #         fallspeedmask = pqc.get_fallspeed_mask(avg_diameter, avg_fall_bins)
        #         vd_matrix_qc_da = pqc.rainfallspeedQC(vd_matrix_qc_da, fallspeedmask)
        #     if qc_attr_dict['rainonlyQC']:
        #         vd_matrix_qc_da = pqc.rainonlyQC(vd_matrix_qc_da)
        #     if qc_attr_dict['hailonlyQC']:
        #         vd_matrix_qc_da = pqc.hailonlyQC(vd_matrix_qc_da)
        #     # if graupelonlyQC:
        #     #     vd_matrix_da = pqc.graupelonlyQC(vd_matrix_da)

        # Compute the number density ND from the V-D matrix
        # Use measured fallspeed by default

        # Resample the parsivel data to a longer interval if desired
        if requested_interval > 10.:
            DSD_interval = pips.check_requested_resampling_interval(requested_interval, 10.)
            vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
            # if apply_qc:
            #     vd_matrix_qc_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_qc_da)
            parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
        else:
            DSD_interval = 10.

        PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
        if start_time is None:
            start_time = PSD_datetimes[0].strftime('%Y%m%d%H%M%S')
        if end_time is None:
            end_time = PSD_datetimes[-1].strftime('%Y%m%d%H%M%S')
        # Resample conventional data to the parsivel times
        sec_offset = PSD_datetimes[0].second
        conv_resampled_df = pips.resample_conv(ptype, DSD_interval, sec_offset, conv_df, gusts=True)
        conv_resampled_df_index = conv_resampled_df.index.intersection(parsivel_df.index)
        conv_resampled_df = conv_resampled_df.loc[conv_resampled_df_index]

        fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, 
                                                          correct_rho=True, 
                                                          rho=conv_resampled_df['rho'])
        vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
        ND_raw = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
        # if apply_qc:
        #     vd_matrix_qc_da = vd_matrix_qc_da.where(vd_matrix_qc_da > 0.0)
        #     ND = pips.calc_ND(vd_matrix_qc_da, fallspeed_spectrum, DSD_interval)

        # Convert resampled conventional data pd.DataFrame to xr.Dataset and add some metadata
        conv_resampled_ds = pipsio.conv_df_to_ds(conv_resampled_df)

        # Convert parsivel telegram data pd.DataFrame to xr.Dataset and add some metadata
        parsivel_ds = pipsio.parsivel_df_to_ds(parsivel_df)

        # Now merge the parsivel telegram Dataset and the resampled conventional Dataset together
        # Then add the vd_matrix_da (raw and QC'ed) and ND (raw and QC'ed) DataArrays.

        parsivel_combined_ds = xr.merge([parsivel_ds, conv_resampled_ds])
        parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, vd_matrix_da,
                                                            name='VD_matrix')
        parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, ND_raw, name='ND')
        # if qc_tag is not None:
        #     new_VD_name = 'VD_matrix_{}'.format(qc_tag)
        #     new_ND_name = 'ND_{}'.format(qc_tag)
        #     parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, 
        #                                                         vd_matrix_qc_da,
        #                                                         name=new_VD_name)
        #     parsivel_combined_ds = pipsio.combine_parsivel_data(parsivel_combined_ds, ND,
        #                                                         name=new_ND_name)
        #     for qc_attr_name in qc_attr_names:
        #         parsivel_combined_ds['VD_matrix_{}'.format(qc_tag)].attrs[qc_attr_name] = \
        #             int(qc_attr_dict[qc_attr_name])
        #         parsivel_combined_ds['ND_{}'.format(qc_tag)].attrs[qc_attr_name] = \
        #             int(qc_attr_dict[qc_attr_name])

        # Add some metadata
        parsivel_combined_ds.attrs['probe_name'] = PIPS_name
        parsivel_combined_ds.attrs['parsivel_angle'] = pp.probe_info[PIPS_name]['parsivel_angle']
        parsivel_combined_ds.attrs['deployment_name'] = deployment_name
        parsivel_combined_ds.attrs['location'] = str(geo_locs[index])
        parsivel_combined_ds.attrs['starting_time'] = start_time
        parsivel_combined_ds.attrs['ending_time'] = end_time
        parsivel_combined_ds.attrs['DSD_interval'] = DSD_interval

        # Dump to netCDF files

        ncfile_name = 'parsivel_combined_{}_{}_{:d}s.nc'.format(deployment_name, PIPS_name,
                                                                int(DSD_interval))
        ncfile_path = os.path.join(output_dir, ncfile_name)
        print("Dumping {}".format(ncfile_path))
        parsivel_combined_ds.to_netcdf(ncfile_path)

    # Convert conventional data pd.DataFrame to xr.Dataset and add some metadata
    conv_ds = pipsio.conv_df_to_ds(conv_df)

    # Now do the same for the 1-s conventional data (put in a separate Dataset and netCDF file
    # due to different time resolution)
    conv_ds.attrs['probe_name'] = PIPS_name
    conv_ds.attrs['parsivel_angle'] = pp.probe_info[PIPS_name]['parsivel_angle']
    conv_ds.attrs['deployment_name'] = deployment_name
    conv_ds.attrs['location'] = str(geo_locs[index])
    conv_ds.attrs['starting_time'] = start_time
    conv_ds.attrs['ending_time'] = end_time


    ncfile_name = 'conventional_raw_{}_{}.nc'.format(deployment_name, PIPS_name)
    ncfile_path = os.path.join(output_dir, ncfile_name)
    print("Dumping {}".format(ncfile_path))
    conv_ds.to_netcdf(ncfile_path)
