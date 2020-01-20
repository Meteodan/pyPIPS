# plot_conv_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import sys
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.dates as dates
import pyPIPS.parsivel_params as pp
import pyPIPS.radarmodule as radar
import pyPIPS.plotmodule as pm
import pyPIPS.pips_io as pipsio
import pyPIPS.utils as utils
import pyPIPS.PIPS as pips
import pyPIPS.parsivel_qc as pqc
import pyPIPS.DSDlib as dsd
import pyPIPS.polarimetric as dp
import pyPIPS.timemodule as tm

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
parsivel_df_list = []
vd_matrix_da_list = []

for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc, ptype in zip(range(
        0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list, ib.starttimes, ib.stoptimes,
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

    # Resample the parsivel data to a longer interval if desired
    if pc.DSD_interval > 10.:
        DSD_interval = pips.check_requested_resampling_interval(pc.DSD_interval, 10.)
        vd_matrix_da = pips.resample_vd_matrix(DSD_interval, vd_matrix_da)
        parsivel_df = pips.resample_parsivel(DSD_interval, parsivel_df)
    else:
        DSD_interval = 10.

    conv_df_list.append(conv_df)
    parsivel_df_list.append(parsivel_df)
    vd_matrix_da_list.append(vd_matrix_da)

# Outer disdrometer (and deployment) loop
for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc, ptype, conv_df, \
    parsivel_df, vd_matrix_da in zip(range(0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list,
                                     ib.starttimes, ib.stoptimes, ib.centertimes, ib.dlocs,
                                     ib.type, conv_df_list, parsivel_df_list, vd_matrix_da_list):

    tripips = (ptype == 'TriPIPS')

    # Calculate some additional thermodynamic parameters from the conventional data
    conv_df = pips.calc_thermo(conv_df)

    # Get times for PIPS meteogram plotting
    conv_datetimes = pips.get_conv_datetimes(conv_df)
    conv_datetimes_nums = dates.date2num(conv_datetimes)
    timelimits = [dates.date2num(datetime.strptime(starttime, tm.timefmt3)),
                  dates.date2num(datetime.strptime(stoptime, tm.timefmt3))]
    # Interval for plotting in seconds.  Setting to larger intervals is useful to avoid
    # plot clutter for long datasets, but information will be lost.
    # TODO: Changing the plot interval from 1 is currently not implemented!
    plotinterval = 1
    plotintervalstr = '{:d}S'.format(int(plotinterval))
    sec_offset = conv_datetimes[0].second

    # Plot wind meteogram
    windavgintv = 60
    windgustintv = 3

    # Resample wind diagnostics array to flag as bad any time in the interval given by
    # plotinterval. Note, have to use numpy max() as a lambda function because the
    # pandas resample.max() does not propagate NaN's!
    # if pc.plot_diagnostics and (ptype == 'PIPS' or ptype == 'TriPIPS'):
    #     winddiag_resampled = conv_df['winddiag'].resample(plotintervalstr, label='right',
    #                                                         closed='right', base=sec_offset,
    #                                                         how=lambda x: utils.trymax(x.values))
    #     conv_plot_df['winddiag'] = winddiag_resampled.loc[
    #         winddiag_resampled.index.intersection(plottimeindex)]

    # Organize a bunch of stuff in a dictionary
    convmeteodict = {
        'plotinterval': plotinterval,
        'plottimes': conv_datetimes_nums,
        'windavgintv': windavgintv,
        'windgustintv': windgustintv,
        'conv_plot_df': conv_df,
        'xaxislimits': timelimits
    }

    pm.plotconvmeteograms(index, pc, ib, convmeteodict)
