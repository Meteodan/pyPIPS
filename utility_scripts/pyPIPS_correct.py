# pyPIPS_merge
#
# This script is based on Sean Waugh's MATLAB script "POD_dataMerger.m"
# Original description follows:
#
# This program is written to combine and process the retrieved SD card data
# from the mobile PIPS. It will merge the two data file types together to
# create a single, line matched record of both information streams. It also
# derives the ambient wind direction, dewpoint, and RH (using the method
# outlined by Richardson et al. (1998).

# Richardson, S. J., S. E. Fredrickson, F. V. Brock, and J. A. Brotzge, 1998:
# Combination temperature and relative humidity probes: avoiding large air
# temperature errors and associated relative humidity errors. Preprints,
# 10th Symposium on Meteorological Observations and Instrumentation,
# Phoenix, AZ, USA, American Meteorological Society, 278-283.

# Original MATLAB script written by Sean Waugh: 1/25/2016
# Contact: sean.waugh@noaa.gov or 405-312-7585
#
# Python version written by Dan Dawson: 09/13/2017
# Latest update 08/30/2020
# Contact: dandawson@purdue.edu

# Import required modules

import os
import argparse
from pyPIPS.pips_io import correct_PIPS

# Parse the command line options
description = "Attempts to correct corrupted Parsivel records in a PIPS comma-delimited file"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--PIPS-dir', dest='PIPS_dir', default=None,
                    help='Directory containing PIPS data file')
parser.add_argument('--input-file-name', dest='input_file_name', default=None,
                    help='Name of PIPS input file')
parser.add_argument('--output-file-name', dest='output_file_name', default=None,
                    help='Name of corrected PIPS output file')
parser.add_argument('--serial-number', dest='serialnum', default=None,
                    help='PIPS serial number')
parser.add_argument('--probe-type', dest='probe_type', default='PIPS',
                    help='Probe type (PIPS or TriPIPS), default PIPS')
args = parser.parse_args()

tripips = (args.probe_type == 'TriPIPS')
print("Attempting to correct corrupted Parsivel records!")
input_path = os.path.join(args.PIPS_dir, args.input_file_name)
output_path = os.path.join(args.PIPS_dir, args.output_file_name)
serialnum = args.serialnum
correct_PIPS(serialnum, input_path, output_path, tripips=tripips)
