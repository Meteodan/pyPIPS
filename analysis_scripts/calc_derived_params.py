# pyPIPS_meteograms.py
#
# This script calculates dereived parameters from the PIPS DSDs (netCDF version)
import os
import argparse
import numpy as np
import xarray as xr
import pyPIPS.parsivel_params as pp
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.DSDlib as dsd

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
parser.add_argument('--ND-tag', dest='ND_tag', default=None,
                    help='Tag for ND variable in file (i.e., qc, RB15_vshift_qc, RB15_qc).')
parser.add_argument('--output-tag', dest='output_tag', default='',
                    help='tag for output nc files to distinguish from original if desired')

args = parser.parse_args()
if args.ND_tag:
    ND_tag = '_{}'.format(args.ND_tag)
    if 'RB15_qc' in args.ND_tag:
        VD_tag = '_RB15_vshift_qc'
    else:
        VD_tag = '_{}'.format(args.ND_tag)
else:
    ND_tag = ''
    VD_tag = ''

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
    ND = parsivel_combined_ds['ND{}'.format(ND_tag)]
    vd_matrix = parsivel_combined_ds['VD_matrix{}'.format(VD_tag)]
    coord_to_combine = 'time'

    vd_matrix = vd_matrix.where(vd_matrix > 0.)
    ND = ND.where(ND > 0.)

    # Compute rainrate using empirical fallspeed curve
    # TODO: allow for the use of the measured fallspeeds in the rainrate calculation.
    fallspeeds_emp = pips.calc_empirical_fallspeed(avg_diameter, correct_rho=True,
                                                   rho=parsivel_combined_ds['rho'])
    rainrate_bin = (6. * 10.**-4.) * np.pi * fallspeeds_emp * avg_diameter**3. * ND * bin_width
    rainrate = rainrate_bin.sum(dim='diameter_bin')
    parsivel_combined_ds['rainrate_derived{}'.format(ND_tag)] = rainrate
    parsivel_combined_ds['rainrate_derived{}'.format(ND_tag)].attrs['units'] = \
        'millimeters per hour'

    # Compute particle counts from raw or QC'ed VD matrix
    pcount = vd_matrix.sum(dim=['fallspeed_bin', 'diameter_bin'])
    parsivel_combined_ds['pcount_derived{}'.format(ND_tag)] = pcount

    # Compute (Rayleigh) radar reflectivity from raw or QC'ed ND
    reflectivity = dsd.calc_dBZ_from_bins(ND)
    parsivel_combined_ds['reflectivity_derived{}'.format(ND_tag)] = reflectivity
    parsivel_combined_ds['reflectivity_derived{}'.format(ND_tag)].attrs['units'] = 'dBZ'

    # Compute median volume diameter D0 (m)
    D0 = dsd.calc_D0_bin(ND)
    parsivel_combined_ds['D0{}'.format(ND_tag)] = D0
    parsivel_combined_ds['D0{}'.format(ND_tag)].attrs['units'] = 'meters'

    # Compute mass weighted mean diameter and spectral width
    D = ND['diameter']
    dD = ND['max_diameter'] - ND['min_diameter']

    Dm = dsd.calc_Dmpq_binned(4, 3, ND)
    sigma = dsd.calc_sigma(D, dD, ND)

    parsivel_combined_ds['Dm43{}'.format(ND_tag)] = Dm
    parsivel_combined_ds['Dm43{}'.format(ND_tag)].attrs['units'] = 'meters'

    parsivel_combined_ds['sigma{}'.format(ND_tag)] = sigma
    parsivel_combined_ds['sigma{}'.format(ND_tag)].attrs['units'] = 'meters'

    parsivel_combined_output_file = parsivel_combined_file + args.output_tag
    print("Dumping {}".format(parsivel_combined_output_file))
    parsivel_combined_ds.to_netcdf(parsivel_combined_output_file)
