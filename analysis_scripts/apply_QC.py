# apply_QC.py
#
# This script applies quality control to the PIPS DSDs (netCDF version)
import os
import argparse
import numpy as np
import xarray as xr
import pyPIPS.parsivel_params as pp
import pyPIPS.parsivel_qc as pqc
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.DSDlib as dsd
import pyPIPS.pips_io as pipsio

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

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

        vd_matrix_qc_da = vd_matrix_da.copy()
        if strongwindQC:
            vd_matrix_qc_da = pqc.strongwindQC(vd_matrix_qc_da)
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

    parsivel_combined_output_file = parsivel_combined_file + args.output_file_tag
    print("Dumping {}".format(parsivel_combined_output_file))
    parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
