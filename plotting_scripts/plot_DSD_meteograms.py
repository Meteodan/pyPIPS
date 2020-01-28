# pyPIPS_meteograms.py
#
# This script plots meteograms from the Portable Integrated Precipitation Stations (PIPS)
import os
import sys
from datetime import datetime, timedelta
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

# ------
# Grab radar data for comparison if desired

if pc.comp_radar:
    if pc.calc_dualpol:
        fieldnames = ['dBZ', 'ZDR', 'RHV', 'Vr']  # Removed KDP for now
    else:
        fieldnames = ['dBZ', 'Vr']

    # sb = radar.readsweeps2PIPS(fieldnames, pc, ib)
    if not pc.loadradopt:
        radar_dict = radar.read_sweeps(ib.radar_name, ib.radar_dir,
                                       ib.starttimerad.strftime(tm.timefmt3),
                                       ib.stoptimerad.strftime(tm.timefmt3),
                                       field_names=fieldnames, el_req=ib.el_req)
        dradlocs = []
        for dloc in ib.dlocs:
            dradloc = radar.get_PIPS_loc_relative_to_radar(dloc, radar_dict['radarsweeplist'][0])
            dradlocs.append(dradloc)
        radar_fields_at_PIPS_da = radar.interp_sweeps_to_PIPS(ib.radar_name,
                                                              radar_dict['radarsweeplist'],
                                                              ib.dis_name_list, dradlocs)
    else:
        print("Reading radar data at PIPS from netCDF file not yet implemented.")

# Outer disdrometer (and deployment) loop
for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc, ptype, conv_df, \
    parsivel_df, vd_matrix_da in zip(range(0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list,
                                     ib.starttimes, ib.stoptimes, ib.centertimes, ib.dlocs,
                                     ib.type, conv_df_list, parsivel_df_list, vd_matrix_da_list):

    tripips = (ptype == 'TriPIPS')

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
    empirical_fallspeed = pips.calc_empirical_fallspeed(avg_diameter)
    fallspeed_spectrum = pips.calc_fallspeed_spectrum(avg_diameter, avg_fall_bins, correct_rho=True,
                                                      rho=conv_df['rho'])
    vd_matrix_da = vd_matrix_da.where(vd_matrix_da > 0.0)
    ND = pips.calc_ND(vd_matrix_da, fallspeed_spectrum, DSD_interval)
    logND = np.log10(ND)

    # Get times for PIPS meteogram plotting
    PSD_datetimes = pips.get_PSD_datetimes(vd_matrix_da)
    PSD_datetimes_dict = pips.get_PSD_time_bins(PSD_datetimes)

    PSD_edgetimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_edges'])
    PSD_centertimes = dates.date2num(PSD_datetimes_dict['PSD_datetimes_centers'])

    # If we are comparing with radar, get the radar timestamps as well
    if pc.comp_radar:
        # plotx_rad = dates.date2num(sb.radtimes)
        plotx_rad = pd.to_datetime(radar_fields_at_PIPS_da['time'].values).to_pydatetime()
    # Pack plotting variables into dictionary
    disvars = {
        'min_diameter': min_diameter,
        'PSDstarttimes': PSD_edgetimes,
        'PSDmidtimes': PSD_centertimes,
        'logND': logND.T
    }

    # Compute additional derived parameters
    disvars['D_0'] = dsd.calc_D0_bin(ND) * 1000.  # Get to mm
    if pc.calc_dualpol:
        # Calculate polarimetric variables using the T-matrix technique
        # Note, may try to use pyDSD for this purpose.
        scattfile = ib.scattdir + 'SCTT_RAIN_fw100.dat'
        # Observed DSD
        dualpol_dis = dp.calpolrain(ib.wavelength, scattfile, ND, bin_width)
        for varname in ['dBZ', 'ZDR', 'KDP', 'RHV']:
            var = dualpol_dis.get(varname, np.empty((0)))
            if var.size:
                # If desired, perform centered running averages
                if pc.avgwindow and False:
                    window = int(pc.avgwindow / DSD_interval)
                    var = pd.Series(var).rolling(
                        window=window, center=True, win_type='triang',
                        min_periods=1).mean().values
                disvars[varname] = var
    # Set up axis parameters
    try:
        timelimits = [dates.date2num(datetime.strptime(starttime, tm.timefmt3)),
                      dates.date2num(datetime.strptime(stoptime, tm.timefmt3))]
    except ValueError:
        timelimits = [dates.date2num(PSD_edgetimes[0]), dates.date2num(PSD_edgetimes[-1])]

    try:
        diamlimits = pc.DSD_D_range
        diamytick = pc.DSD_D_ytick
    except Exception:
        diamlimits = [0.0, 9.0]
        diamytick = 1.0

    DSDtype = 'observed'
    locator = dates.MinuteLocator(byminute=[0, 15, 30, 45])
    minorlocator = dates.MinuteLocator(byminute=range(0, 60, 5))
    dateformat = '%H:%M'
    formatter = dates.DateFormatter(dateformat)

    axparams = {
        'majorxlocator': locator,
        'majorxformatter': formatter,
        'minorxlocator': minorlocator,
        'axeslimits': [timelimits, diamlimits],
        'majorylocator': ticker.MultipleLocator(base=diamytick),
        'axeslabels': [None, 'D (mm)']
    }

    # If we are comparing with radar, grab associated radar variables for plotting
    if pc.comp_radar:
        # At least add the reflectivity
        dBZ_D_plt = radar_fields_at_PIPS_da.sel(fields='REF', PIPS=dis_name)
        # indexrad = sb.outfieldnames.index('dBZ')
        # dBZ_D_plt = sb.fields_D_tarr[:, index, indexrad]
        radvars = {'radmidtimes': plotx_rad, 'REF': dBZ_D_plt}
        # Add other polarimetric fields
        if pc.calc_dualpol:
            for radvarname in ['ZDR', 'KDP', 'RHO']:
                if radvarname in radar_fields_at_PIPS_da.fields:
                    # indexrad = sb.outfieldnames.index(radvarname)
                    # dualpol_rad_var = sb.fields_D_tarr[:, index, indexrad]
                    dualpol_rad_var = radar_fields_at_PIPS_da.sel(fields=radvarname, PIPS=dis_name)
                    if dualpol_rad_var.size:
                        radvars[radvarname] = dualpol_rad_var
            if pc.clean_radar:
                # remove non-precipitation echoes from radar data
                gc_mask = np.where((radvars['RHO'] < 0.90), True, False)
                for radvarname in ['ZDR', 'REF', 'RHO']:
                    radvars[radvarname] = np.ma.masked_array(radvars[radvarname],
                                                             mask=gc_mask)

    # Make the plot
    dis_plot_name = dis_name + '_' + DSDtype
    pm.plotDSDmeteograms(dis_plot_name, meteogram_image_dir, axparams, disvars, radvars.copy())
