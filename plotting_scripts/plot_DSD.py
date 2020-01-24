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

    # Get times for plotting
    PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])

    # Resample conventional data to the parsivel times
    sec_offset = PSD_datetimes[0].second
    conv_resampled_df = pips.resample_conv(ptype, DSD_interval, sec_offset, conv_df)
    conv_resampled_df_index = conv_resampled_df.index.intersection(parsivel_df.index)
    conv_resampled_df = conv_resampled_df.loc[conv_resampled_df_index]

    vel_D_image_dir = os.path.join(ib.image_dir, 'vel_D/{}'.format(dis_name))
    if not os.path.exists(vel_D_image_dir):
        os.makedirs(vel_D_image_dir)

    axdict = {
        'min_diameter': min_diameter,
        'avg_diameter': avg_diameter,
        'min_fall_bins': min_fall_bins,
        'xlim': pc.D_range,
        'ylim': pc.vel_range,
        'dis_name': dis_name
        }

    for t, time in enumerate(vd_matrix_da['time_10s'].to_index()):
        print(parsivel_df['pcount'].loc[time])
        if parsivel_df['pcount'].loc[time] > 0:
            axdict['time'] = time
            PSDdict = {
                'vd_matrix_da': vd_matrix_da[t, :],
                # 'flaggedtime': parsivel_df['flaggedtimes'].values[t]
            }
            fig, ax = pm.plot_vel_D(ib, axdict, PSDdict, conv_resampled_df['rho'].loc[time])
            image_name = dis_name + '_vel_D_{}_t{:04d}.png'.format(time.strftime(tm.timefmt3), t)
            image_path = os.path.join(ib.image_dir, 'vel_D/', dis_name, image_name)
            fig.savefig(image_path, dpi=200, bbox_inches='tight')
            plt.close(fig)
