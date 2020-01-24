# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import sys
import numpy as np
import matplotlib.dates as dates
import matplotlib.pyplot as plt
import pyPIPS.parsivel_params as pp
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.timemodule as tm
import pyPIPS.DSDlib as dsd
import xarray as xr
import pandas as pd

min_diameter = pp.parsivel_parameters['min_diameter_bins_mm']
max_diameter = pp.parsivel_parameters['max_diameter_bins_mm']
bin_width = max_diameter - min_diameter
avg_diameter = pp.parsivel_parameters['avg_diameter_bins_mm']
min_fall_bins = pp.parsivel_parameters['min_fallspeed_bins_mps']
max_fall_bins = pp.parsivel_parameters['max_fallspeed_bins_mps']
avg_fall_bins = pp.parsivel_parameters['avg_fallspeed_bins_mps']

# Create V-D relationship for rain based on Terry Schuur's relationship
# rainvd = dis.assignfallspeed(avg_diameter)

# -----------------------------------------------------------------------
#
#   Dynamically import pyPIPScontrol.py or user version of it.
#
# -----------------------------------------------------------------------

if len(sys.argv) > 1:   # Try to import user-defined plotcontrol module
    controlmodpath = sys.argv[1]
    utils.log("Input file is " + controlmodpath)
    pc = utils.import_all_from(controlmodpath)
    try:
        pc = utils.import_all_from(controlmodpath)
        utils.log("Successfully imported pyPIPS control parameters!")
    except Exception:
        utils.warning(
            "Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
        import pyPIPScontrol as pc
else:   # Read in default plotcontrol.py
    import pyPIPS.pyPIPScontrol as pc

# Parse command line argument and read in disdrometer/radar information from text file
# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
# version first, and if not, reading the textfile and pickling the read-in
# variables for future read-ins. Think about this.
# TODO: Change this input file to a python file

if len(sys.argv) == 1:
    argindex = 1
elif len(sys.argv) > 1:
    argindex = 2
else:
    sys.exit("No text input file defined! Quitting!")

ib = utils.readpyPIPSinput(sys.argv[argindex])

if not os.path.exists(ib.image_dir):
    os.makedirs(ib.image_dir)

# Create the directory for the meteogram plots if it doesn't exist
meteogram_image_dir = os.path.join(ib.image_dir, 'meteograms')
if not os.path.exists(meteogram_image_dir):
    os.makedirs(meteogram_image_dir)

# Read in the PIPS data for the deployment
conv_df_list = []
conv_resampled_df_list = []
parsivel_df_list = []
vd_matrix_da_list = []

# Outer disdrometer (and deployment) loop
for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc, ptype in \
        zip(range(0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list, ib.starttimes, ib.stoptimes,
            ib.centertimes, ib.dlocs, ib.type):

    if starttime == '-1':
        starttime = None
    if stoptime == '-1':
        stoptime = None

    tripips = (ptype == 'TriPIPS')
    PIPS_data_file_path = os.path.join(ib.dis_dir, dis_filename)
    conv_df, parsivel_df, vd_matrix_da = pipsio.read_PIPS(PIPS_data_file_path,
                                                          starttimestamp=starttime,
                                                          stoptimestamp=stoptime, tripips=tripips)

    # We need the disdrometer locations. If they aren't supplied in the input control file, find
    # them from the GPS data
    if np.int(dloc[0]) == -1:
        ib.dlocs[index] = pipsio.get_PIPS_loc(conv_df['GPS_status'], conv_df['GPS_lat'],
                                              conv_df['GPS_lon'], conv_df['GPS_alt'])
    print("Lat/Lon/alt of {}: {}".format(dis_name, str(ib.dlocs[index])))

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Do some QC on the V-D matrix
    strongwindQC = pc.strongwindQC
    splashingQC = pc.splashingQC
    marginQC = pc.marginQC
    rainfallQC = pc.rainfallQC
    rainonlyQC = pc.rainonlyQC
    hailonlyQC = pc.hailonlyQC
    graupelonlyQC = pc.graupelonlyQC

    if pc.basicQC:
        strongwindQC = True
        splashingQC = True
        marginQC = True

    if strongwindQC:
        vd_matrix_da = pqc.strongwindQC(vd_matrix_da)
    if splashingQC:
        vd_matrix_da = pqc.splashingQC(vd_matrix_da)
    if marginQC:
        vd_matrix_da = pqc.marginQC(vd_matrix_da)
    if rainfallQC:
        fallspeedmask = pqc.get_fallspeed_mask(avg_diameter, avg_fall_bins)
        vd_matrix_da = pqc.rainfallspeedQC(vd_matrix_da, fallspeedmask)
    if rainonlyQC:
        vd_matrix_da = pqc.rainonlyQC(vd_matrix_da)
    if hailonlyQC:
        vd_matrix_da = pqc.hailonlyQC(vd_matrix_da)
    # if graupelonlyQC:
    #     vd_matrix_da = pqc.graupelonlyQC(vd_matrix_da)

    # Compute the number density ND from the V-D matrix
    # Use measured fallspeed by default

    # Resample the parsivel data to a longer interval if desired
    if pc.DSD_interval > 10.:
        DSD_interval = pips.check_requested_resampling_interval(pc.DSD_interval, 10.)
        vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
        parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
    else:
        DSD_interval = 10.

    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    logND = np.log10(ND)

    # Compute various fits using the MM and TMM

    M2, _ = dsd.calc_moment_bin(ND, moment=2)
    M3, _ = dsd.calc_moment_bin(ND, moment=3)
    M4, _ = dsd.calc_moment_bin(ND, moment=4)
    M6, _ = dsd.calc_moment_bin(ND, moment=6)

    DSD_MM24 = dsd.fit_DSD_MM24(M2, M4)
    DSD_MM36 = dsd.fit_DSD_MM36(M3, M6)
    DSD_MM346 = dsd.fit_DSD_MM346(M3, M4, M6)
    DSD_MM246 = dsd.fit_DSD_MM246(M2, M4, M6)
    DSD_MM234 = dsd.fit_DSD_MM234(M2, M3, M4)

    D_min, D_max = dsd.get_max_min_diameters(ND)
    DSD_TMM246 = dsd.fit_DSD_TMM_xr(M2, M4, M6, D_min, D_max)

    # Wrap fits into a DataSet and dump to netCDF file
    data_arrays = []
    for da_tuple in [DSD_MM24, DSD_MM36, DSD_MM346, DSD_MM246, DSD_MM234, DSD_TMM246]:
        da_concat = xr.concat(da_tuple, pd.Index(['N0', 'lamda', 'alpha'], name='parameter'))
        data_arrays.append(da_concat)
    names = ['DSD_MM24', 'DSD_MM36', 'DSD_MM346', 'DSD_MM246', 'DSD_MM234', 'DSD_TMM246']
    fits_ds = xr.Dataset({name: da for name, da in zip(names, data_arrays)})

    ncfile_name = 'DSD_fits_{}_{}_{}.nc'.format(dis_name, starttime, stoptime)
    ncfile_path = os.path.join(ib.dis_dir, ncfile_name)
    print("Dumping {}".format(ncfile_path))
    fits_ds.to_netcdf(ncfile_path)
