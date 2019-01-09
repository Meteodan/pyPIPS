# Plot_Disdrometer
# Runs through all of the FMCW days and does meteograms and scattergrams and calculates overall lam-mu relation
# This script plots several quantities based on disdrometer data from the Parsivel laser disdrometer
# command line: python pyPIPS_CG.py pyPIPScontrol.py

import numpy as N
from numpy import ma as ma
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
from datetime import datetime, timedelta
import modules.thermolib as thermo
import os
import modules.disdrometer_module as dis
import sys
import modules.radarmodule as radar
import pandas as pd
import modules.plotmodule as pm
import modules.utils as utils
import modules.DSDretrieval as DR
import modules.empirical_module as em
import modules.timemodule as tm

mu=[]
lamda=[]
muTMF = []
lamTMF=[]
sigm_obs=[]
sigm_dis=[]
sigm_rad=[]
Dm_obs=[]
Dm_dis=[]
Dm_rad=[]
R_list = []
R_mm = []
R_retr = []
D0_list = []
D0_mm = []
D0_retr = []
Nt_obs = []
Nt_mm = []
Nt_retr = []
W_obs = []
W_mm = []
W_retr = []
Zh_dis = []
Zh_retr = []
ND_list = N.empty((0,32))

fieldnames = ['dBZ','ZDR','KDP','RHV','Vr']
deg2rad = N.pi/180.

# Disdrometer bin and diameter information
min_diameter = dis.min_diameter
max_diameter = dis.max_diameter
bin_width = max_diameter - min_diameter
avg_diameter = dis.avg_diameter
fall_bins = dis.fall_bins
min_fall_bins = dis.min_fall_bins

# Should make this an input parameter
outer_image_dir = '/Volumes/depot/dawson29/data/VORTEXSE/obsdata/dan_images/'

if (not os.path.exists(outer_image_dir)):
    os.makedirs(outer_image_dir)

# Z, ZDR relation from Cao et al. (2008)
Zh_Cao = N.arange(20, 61, 1)
Zdr_Cao = 10**((-2.6857 * 10**-4 * Zh_Cao**2) + 0.04892 * Zh_Cao - 1.4287)


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
    import pyPIPScontrol as pc

# Adjust certain plotting flags for consistency
if(not pc.calc_DSD):
    pc.calc_dualpol = False

if(not pc.plot_DSD_meteo):
    pc.plot_DSDderived = False

if(not pc.plot_conv_meteo and not pc.plot_DSD_meteo):
    pc.plot_diagnostics = False

#-----------------------------------------------------------------------
#
#   Loop to run through each FMCW day.
#
#-----------------------------------------------------------------------

directories = ['/Users/bozell/pyPIPS_work/input/NEXRAD/','/Volumes/depot/dawson29/data/VORTEXSE/obsdata/2017/PIPS/PIPS2A_FMCW/']
for directory in directories:
    for root, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith(("040517.txt","040417.txt","041717.txt","042017.txt","042417.txt","FMCW.txt")):
                continue
            elif f.endswith(".txt"):
                print directory
                print f

                filepath = os.path.join(root,f)

                if(directory == '/Volumes/depot/dawson29/data/VORTEXSE/obsdata/2017/PIPS/PIPS2A_FMCW/'):
                    inputdict = {}
                    inputdict['dis_list'] = [f]
                    dis_name = 'PIPS_2A'
                    inputdict['dis_name_list'] = [dis_name]
                    GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc = dis.readPIPSloc(filepath)
                    inputdict['dlocs'] = [dloc]
                    print "Lat/Lon/alt of "+dis_name+": "+str(dloc)

                    # Read in the disdrometer data file
                    PIPS_dict = dis.readPIPS(filepath, basicqc=pc.basicQC, rainfallqc=pc.rainfallQC,
                                         rainonlyqc=pc.rainonlyQC, strongwindqc=pc.strongwindQC,
                                         DSD_interval=pc.DSD_interval)

                     # Unpack some stuff from the PIPS_dict
                    convtimestamps = PIPS_dict['convtimestamps']
                    PSDtimestamps = PIPS_dict['PSDtimestamps']
                    conv_df = PIPS_dict['conv_df']
                    PSD_df = PIPS_dict['PSD_df']
                    ND = PIPS_dict['ND']
                    ND_onedrop = PIPS_dict['ND_onedrop']

                    ND = ND.T
                    ND_onedrop = ND_onedrop.T


                    logND = N.ma.log10(ND)
                    logND = N.ma.masked_where(ND < ND_onedrop, logND)

                    DSD_index = PIPS_dict['DSD_index']
                    DSD_interval = PIPS_dict['DSD_interval']
                    DSD_interval_td = timedelta(seconds=DSD_interval)
                    DSD_halfinterval_td = timedelta(seconds=DSD_interval / 2.)

                    # Compute potential temperature, water vapor mixing ratio, and density

                    conv_df['pt'] = thermo.caltheta(conv_df['pressure']*100., conv_df['fasttemp']+273.15)
                    conv_df['qv'] = thermo.calqv(conv_df['RH_derived']/100., conv_df['pressure']*100.,
                                                   conv_df['fasttemp']+273.15)
                    conv_df['rho'] = thermo.calrho(conv_df['pressure']*100., conv_df['pt'], conv_df['qv'])

                    max_ref = N.max(PSD_df['reflectivity'].values)

                    print max_ref
                    print N.sum(PSD_df['pcount2'].values)
                    if(max_ref <= 20.):
                        print "SKIPPING NON-RAIN DAY"
                        continue

                    # Determine start and end times/indices for analysis

                    convtimestampsnums = dates.date2num(convtimestamps)
                    PSDtimestampsnums = dates.date2num(PSDtimestamps)

                    startindex = 0
                    pstartindex = 0
                    starttime = convtimestampsnums[startindex]
                    inputdict['starttimes'] = [starttime]
                    pstarttime = PSDtimestampsnums[startindex]
                    plotstarttime = starttime
                    radar_date = convtimestamps[0].strftime(tm.timefmt4).strip()

                    stopindex = N.size(convtimestampsnums)-1
                    pstopindex = N.size(PSDtimestampsnums)-1
                    stoptime = convtimestampsnums[stopindex]
                    inputdict['stoptimes'] = [stoptime]
                    pstoptime = PSDtimestampsnums[pstopindex]
                    plotstoptime = stoptime

                    if(pc.timetospace):
                        centertime = -1
                    inputdict['centertimes']=[-1]

                    PSDtimestamps_edge = [x - DSD_interval_td for x in PSDtimestamps]
                    # Add an extra 10 sec for the last time bin boundary
                    PSDtimestamps_edge.append(PSDtimestamps_edge[-1] + DSD_interval_td)
                    PSDtimestamps_avg = [x - DSD_halfinterval_td for x in PSDtimestamps]
                    PSDstarttimes = dates.date2num(PSDtimestamps_edge[pstartindex:pstopindex + 1])
                    PSDmidtimes = dates.date2num(PSDtimestamps_avg[pstartindex:pstopindex + 1])

                    inputdict['platform'] = 'NEXRAD'
                    inputdict['wavelength'] = 10.7

                    # Read in radar name,lat,lon
                    inputdict['radar_name'] = 'KHTX'
                    inputdict['rlat'] = None
                    inputdict['rlon'] = None
                    inputdict['ralt'] = None
                    inputdict['el_req'] = 0.5
                    inputdict['heading'] = None

                    inputdict['image_dir'] = '/Users/dawson29/pyPIPS_work/VORTEXSE/PIPS2A_FMCW/NEXRAD/'+radar_date+'/'
                    inputdict['radar_dir'] = '/Volumes/depot/dawson29/data/VORTEXSE/obsdata/2017/NEXRAD/PIPS2A_FMCW/'+radar_date+'/CFRadial/'
                    inputdict['scattdir'] = '/Users/dawson29/pyPIPS/tmatrix/S-band/'

                    raddate = convtimestamps[0].strftime(tm.timefmt5).strip().split(',')
                    raddate_int = map(int,raddate)
                    # Read in start and end times for the radar data analysis (int)
                    inputdict['starttimerad'] = datetime(raddate_int[0],raddate_int[1],raddate_int[2],int(00),int(00),int(00))
                    inputdict['stoptimerad'] = datetime(raddate_int[0],raddate_int[1],raddate_int[2],int(23),int(59),int(59))
                    print "radar starttime:", inputdict['starttimerad']
                    print "radar stoptime:", inputdict['stoptimerad']
                    inputdict['radlims'] = [0.0,250000.0,0.,360.]
                    plotxmin = -1
                    inputdict['plotlims'] = [-1,-1,-1,-1]
                    ib = utils.Bunch(inputdict)
                else:
                    ib = utils.readpyPIPSinput(filepath)

                    # We need the disdrometer locations. If they aren't supplied in the input control file, find them
                    # from the GPS data

                    for index, dis_name, dis_filename, dloc in \
                            zip(xrange(0, ib.numdis), ib.dis_name_list, ib.dis_list, ib.dlocs):

                        if(N.int(dloc[0]) == -1):
                            filepath = os.path.join(ib.dis_dir, dis_filename)
                            GPS_lats, GPS_lons, GPS_stats, GPS_alts, dloc = dis.readPIPSloc(
                                filepath)
                            ib.dlocs[index] = dloc

                        print "Lat/Lon/alt of " + dis_name + ": " + str(dloc)

                if (not os.path.exists(ib.image_dir)):
                    os.makedirs(ib.image_dir)

                if(pc.plot_conv_meteo or pc.plot_DSD_meteo):
                    # Create the directory for the plots if it doesn't exist
                    meteogram_image_dir = ib.image_dir + 'meteograms/'
                    if (not os.path.exists(meteogram_image_dir)):
                        os.makedirs(meteogram_image_dir)
                # ------
                # Grab radar data for comparison if desired

                if(pc.comp_radar):
                    sb = radar.readsweeps2PIPS(fieldnames, pc, ib)

                    # Outer disdrometer (and deployment) loop
                for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc in zip(xrange(
                        0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list, ib.starttimes, ib.stoptimes,
                        ib.centertimes, ib.dlocs):

                    if(directory == '/Users/bozell/pyPIPS_work/input/NEXRAD/'):
                        # Read in the disdrometer data file
                        filepath = os.path.join(ib.dis_dir, dis_filename)

                        PIPS_dict = dis.readPIPS(filepath, basicqc=pc.basicQC,
                                                 rainfallqc=pc.rainfallQC,
                                                 rainonlyqc=pc.rainonlyQC,
                                                 strongwindqc=pc.strongwindQC,
                                                 DSD_interval=pc.DSD_interval)
                        # Unpack some stuff from the PIPS_dict
                        convtimestamps = PIPS_dict['convtimestamps']
                        PSDtimestamps = PIPS_dict['PSDtimestamps']
                        conv_df = PIPS_dict['conv_df']
                        PSD_df = PIPS_dict['PSD_df']
                        ND = PIPS_dict['ND']
                        ND_onedrop = PIPS_dict['ND_onedrop']

                        ND = ND.T
                        ND_onedrop = ND_onedrop.T

                        logND = N.ma.log10(ND)
                        logND = N.ma.masked_where(ND < ND_onedrop, logND)

                        DSD_index = PIPS_dict['DSD_index']
                        DSD_interval = PIPS_dict['DSD_interval']
                        DSD_interval_td = timedelta(seconds=DSD_interval)
                        DSD_halfinterval_td = timedelta(seconds=DSD_interval / 2.)

                        # Compute potential temperature, water vapor mixing ratio, and density

                        conv_df['pt'] = thermo.caltheta(conv_df['pressure']*100., conv_df['fasttemp']+273.15)
                        conv_df['qv'] = thermo.calqv(conv_df['RH_derived']/100., conv_df['pressure']*100.,
                                                       conv_df['fasttemp']+273.15)
                        conv_df['rho'] = thermo.calrho(conv_df['pressure']*100., conv_df['pt'], conv_df['qv'])

                        # Determine start and end times/indices for analysis

                        convtimestampsnums = dates.date2num(convtimestamps)
                        PSDtimestampsnums = dates.date2num(PSDtimestamps)

                        startindex, stopindex = utils.getTimeWindow(starttime, stoptime, convtimestampsnums)
                        pstartindex, pstopindex = utils.getTimeWindow(starttime, stoptime, PSDtimestampsnums)

                        starttime = convtimestampsnums[startindex]
                        stoptime = convtimestampsnums[stopindex]
                        pstarttime = PSDtimestampsnums[pstartindex]
                        pstoptime = PSDtimestampsnums[pstopindex]

                        plotstarttime = starttime
                        plotstoptime = stoptime

                        PSDtimestamps_edge = [x - DSD_interval_td for x in PSDtimestamps]
                        # Add an extra 10 sec for the last time bin boundary
                        PSDtimestamps_edge.append(PSDtimestamps_edge[-1] + DSD_interval_td)
                        PSDtimestamps_avg = [x - DSD_halfinterval_td for x in PSDtimestamps]
                        PSDstarttimes = dates.date2num(PSDtimestamps_edge[pstartindex:pstopindex + 1])
                        PSDmidtimes = dates.date2num(PSDtimestamps_avg[pstartindex:pstopindex + 1])

                        # Store all time-related parameters in a dictionary (not actually used right now)
                        timedict = {'startindex': startindex, 'stopindex': stopindex, 'pstartindex': pstartindex,
                                    'starttime': starttime, 'stoptime': stoptime, 'pstarttime': pstarttime,
                                    'pstoptime': pstoptime, 'convtimestamps': convtimestamps,
                                    'convtimestampnums': convtimestampsnums, 'PSDtimestamps': PSDtimestamps,
                                    'PSDtimestampsnums': PSDtimestampsnums, 'PSDtimestamps_edge': PSDtimestamps_edge,
                                    'PSDtimestamps_avg': PSDtimestamps_avg}

                    if(pc.filter_bimodal):
                        for t in range(N.size(ND, axis=1)):
                            for d in xrange(16,30):
                                if (ND[d,t] > ND[d+1,t]*0.5+ND[d+1,t] and ND[d,t] > ND[d-1,t]):
                                    ND[:,t] = N.nan
                                    break

                    # TODO: Address this section
                    if(pc.timetospace and pc.comp_radar):
                        delta_t_rad = [(x - centertime).total_seconds() for x in sb.radtimes]
                        delta_t_dis = [(x - centertime).total_seconds() for x in PSDtimestamps]
                        delta_t_dis_start = [(x - centertime).total_seconds()
                                             for x in PSDtimestamps_edge]
                        delta_t_dis_avg = [(x - centertime).total_seconds()
                                           for x in PSDtimestamps_avg]

                        xloc_rad = [(x * pc.stormmotion) / 1000.0 for x in delta_t_rad]
                        xloc_dis = [(x * pc.stormmotion) / 1000.0 for x in delta_t_dis]
                        xloc_dis_start = [(x * pc.stormmotion) /
                                          1000.0 for x in delta_t_dis_start]
                        xloc_dis_avg = [(x * pc.stormmotion) / 1000.0 for x in delta_t_dis_avg]

                        plotx_rad = N.array(xloc_rad)
                        plotx_dis = N.array(xloc_dis)
                        plotx_dis_start = N.array(xloc_dis_start)
                        plotx_dis_avg = N.array(xloc_dis_avg)
                    else:
                        if(pc.comp_radar):
                            plotx_rad = dates.date2num(sb.radtimes)
                        plotx_dis = dates.date2num(PSDtimestamps)
                        plotx_dis_start = dates.date2num(PSDtimestamps_edge)
                        plotx_dis_avg = dates.date2num(PSDtimestamps_avg)

                    # Resample the 1-s data corresponding to each integrated DSD (for the
                    # given DSD interval). Note, most of these aren't used for right now.
                    sec_offset = PSDtimestamps[0].second

                    conv_resampled_df = dis.resampleconv(DSD_interval, sec_offset, conv_df)
                    conv_resampled_df = conv_resampled_df.loc[DSD_index]
                    rho_tDSD = conv_resampled_df['rho']

                    # Now for the fun part.  Calculate exponential and gamma distribution size distribution
                    # parameters using the method of moments, after Zhang et al. 2008 and Tokay and
                    # Short 1996
                    if(pc.calc_DSD):
                        synthbins, exp_DSD, gam_DSD, tmf_DSD, dis_DSD = dis.calc_DSD(pc,
                            min_diameter, avg_diameter, max_diameter, bin_width, ND, logND, rho_tDSD.values, pc.qrQC,
                            pc.qr_thresh, PSD_df['pcount2'].values, PSD_df['intensity'].values)

                        # Unpack needed values from returned tuples

                        ND_expDSD, N0_exp, lamda_exp, mu_exp, qr_exp, Ntr_exp, refl_DSD_exp, D_med_exp, D_m_exp = \
                            exp_DSD
                        ND_gamDSD, N0_gam, lamda_gam, mu_gam, qr_gam, Ntr_gam, refl_DSD_gam, D_med_gam, D_m_gam, \
                            LWC_gam, rainrate_gam = gam_DSD
                        ND_tmfDSD, N0_tmf, lamda_tmf, mu_tmf, qr_tmf, Ntr_tmf, refl_DSD_tmf, \
                            LWC_tmf, rainrate_tmf = tmf_DSD
                        ND, logND, D_med_disd, D_m_disd, D_mv_disd, D_ref_disd, QR_disd, refl_disd, LWC_disd, M0, rainrate = \
                            dis_DSD

                        ND_expDSD = ND_expDSD.T
                        logND_expDSD = N.ma.log10(ND_expDSD / 1000.)  # Get to log(m^-3 mm^-1)
                        logND_expDSD = N.ma.masked_where(ND_expDSD < ND_onedrop, logND_expDSD)

                        ND_gamDSD = ND_gamDSD.T
                        logND_gamDSD = N.ma.log10(ND_gamDSD / 1000.)  # Get to log(m^-3 mm^-1)
                        logND_gamDSD = N.ma.masked_where(ND_gamDSD < ND_onedrop, logND_gamDSD)

                        ND_tmfDSD = ND_tmfDSD.T
                        logND_tmfDSD = N.ma.log10(ND_tmfDSD / 1000.)  # Get to log(m^-3 mm^-1)
                        logND_tmfDSD = N.ma.masked_where(ND_tmfDSD < ND_onedrop, logND_tmfDSD)

                        if(pc.calc_dualpol):
                            # Calculate polarimetric variables using the T-matrix technique
                            # Note, may try to use pyDSD for this purpose.
                            scattfile = ib.scattdir + 'SCTT_RAIN_fw100.dat'
                            # Observed DSD: TODO, need to make ND, ND_expDSD, and ND_gamDSD have consistent units!
                            dualpol_dis = dis.calpolrain(ib.wavelength, scattfile, ND, bin_width)
                            # Exponential DSD
                            dualpol_exp = dis.calpolrain(ib.wavelength, scattfile, (ND_expDSD / 1000.), bin_width)
                            # Gamma DSD
                            dualpol_gam = dis.calpolrain(ib.wavelength, scattfile, (ND_gamDSD / 1000.), bin_width)

                    if(pc.plot_DSD_meteo):
                        PSDplottimes = PSDtimestampsnums[pstartindex: pstopindex + 1]
                        plotPSDindex = pd.DatetimeIndex(PSDtimestamps[pstartindex: pstopindex + 1])
                        PSD_plot_df = PSD_df.loc[plotPSDindex]
                        # Loop through the observed, exponential, and untruncated and truncated gamma fits
                        for DSDtype in ['observed', 'exponential', 'gamma', 'tmf']:
                            if(DSDtype == 'observed'):
                                logND_plot = logND[:, pstartindex:pstopindex + 1]
                                if(pc.calc_DSD):
                                    D_0_plot = D_med_disd
                                    refl_ray_plot = refl_disd
                                    dp = dualpol_dis
                                    #Zh, Zv, Zhv, dBZ, ZDR, KDP, RHV, intv, d, fa2, fb2 = dualpol_dis
                            elif(DSDtype == 'exponential'):
                                logND_plot = logND_expDSD[:, pstartindex:pstopindex + 1]
                                if(pc.calc_DSD):
                                    D_0_plot = D_med_exp
                                    refl_ray_plot = refl_DSD_exp
                                    dp = dualpol_exp
                                    #Zh, Zv, Zhv, dBZ, ZDR, KDP, RHV, intv, d, fa2, fb2 = dualpol_exp
                            elif(DSDtype == 'gamma'):
                                logND_plot = logND_gamDSD[:, pstartindex:pstopindex + 1]
                                if(pc.calc_DSD):
                                    D_0_plot = D_med_gam
                                    refl_ray_plot = refl_DSD_gam
                                    dp = dualpol_gam
                                    #Zh, Zv, Zhv, dBZ, ZDR, KDP, RHV, intv, d, fa2, fb2 = dualpol_gam

                            disvars = {'min_diameter': min_diameter, 'PSDstarttimes': PSDstarttimes,
                                       'PSDmidtimes': PSDmidtimes, 'logND': logND_plot}

                            # Plot one meteogram for each dualpol variable, otherwise just plot a meteogram for
                            # reflectivity based on 6th moment of disdrometer DSD (Rayleigh approximation).
                            # Important! Currently the dualpol calculation for the disdrometer data assumes rain
                            # only and otherwise ignores data in all bins larger than 9 mm.  For this reason it is
                            # recommended to first set the "rain_only_QC" flag to True in disdrometer_module.py

                            if(pc.calc_DSD):
                                # If desired, perform centered running averages
                                if(pc.avgwindow):
                                    window = int(pc.avgwindow / DSD_interval)

                                    D_0_plot = pd.Series(D_0_plot).rolling(
                                        window=window,
                                        center=True,
                                        win_type='triang',
                                        min_periods=1).mean().values  # use for gamma fit
                                    refl_ray_plot = pd.Series(refl_ray_plot).rolling(
                                        window=window,
                                        center=True,
                                        win_type='triang',
                                        min_periods=1).mean().values  # use for gamma fit
                                disvars['D_0'] = D_0_plot[pstartindex:pstopindex + 1]
                                disvars['dBZ_ray'] = refl_ray_plot[pstartindex:pstopindex + 1]

                                if(pc.calc_dualpol):
                                    # Computed centered running averages of dualpol variables
                                    for varname in ['dBZ', 'ZDR', 'KDP', 'RHV']:
                                        var = dp.get(varname, N.empty((0)))
                                        if(var.size):
                                            # If desired, perform centered running averages
                                            if(pc.avgwindow):
                                                var = pd.Series(var).rolling(
                                                    window=window, center=True, win_type='triang',
                                                    min_periods=1).mean().values
                                            disvars[varname] = var[pstartindex:pstopindex + 1]

                            # Mark flagged times with vertical lines, if desired
                            if(pc.plot_diagnostics):
                                if(DSDtype == 'observed'):
                                    # Extract indices for flagged times
                                    flaggedtimes_index = N.where(
                                        PSD_plot_df['flaggedtimes'].values > 0)[0]
                                    # these are the times with wind contamination
                                    flaggedtimes_plot = PSDmidtimes[flaggedtimes_index]
                                    hailflag_index = N.where(PSD_plot_df['hailflag'].values)[0]
                                    # these are the times with hail detected
                                    hailflag_plot = PSDmidtimes[hailflag_index]

                                disvars['flaggedtimes'] = flaggedtimes_plot
                                disvars['hailflag'] = hailflag_plot

                            # Get radar variables together for comparison, if desired
                            if(pc.comp_radar and DSDtype == 'observed'):
                                radvars = {}
                                # At least add the reflectivity
                                indexrad = sb.outfieldnames.index('dBZ')
                                dBZ_D_plt = sb.fields_D_tarr[:, index, indexrad]
                                radvars = {'radmidtimes': plotx_rad, 'dBZ': dBZ_D_plt}
                                # Add other polarimetric fields
                                if(pc.calc_dualpol):
                                    for radvarname in ['ZDR', 'KDP', 'RHV']:
                                        if radvarname in sb.outfieldnames:
                                            indexrad = sb.outfieldnames.index(radvarname)
                                            dualpol_rad_var = sb.fields_D_tarr[:, index, indexrad]
                                            if(dualpol_rad_var.size):
                                                radvars[radvarname] = dualpol_rad_var
                                    if(pc.clean_radar):
                                        # remove non-precipitation echoes from radar data
                                        gc_mask = N.where(radvars['RHV'] < 0.90, True, False)
                                        for radvarname in ['ZDR', 'dBZ', 'RHV']:
                                                radvars[radvarname] = ma.masked_array(radvars[radvarname],
                                                                                      mask=gc_mask)
                            if(DSDtype == 'observed'):
                                # Prepare axis parameters
                                timelimits = [plotstarttime, plotstoptime]
                                diamlimits = [0.0, 9.0]

                                axparams = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                                            'minorxlocator': pc.minorlocator,
                                            'axeslimits': [timelimits, diamlimits],
                                            'majorylocator': ticker.MultipleLocator(base=1.0),
                                            'axeslabels': [None, 'D (mm)']}
                                if(pc.plot_only_precip and pc.calc_dualpol):
                                    # set plot start and end time to the start and end of precipitation
                                    rainindex = N.where(disvars['RHV'] > 0.6)
                                    raintimes = PSDmidtimes[rainindex]
                                    if (len(raintimes) == 0.0):
                                        break
                                    plotstarttime = raintimes[0]
                                    plotstoptime = raintimes[-1]
                                    # set radar PPI plot start and end time for FMCW days
                                    sweeptimes = PSDtimestamps[rainindex]
                                    sweepstart = sweeptimes[0]
                                    sweepstop = sweeptimes[-1]
                                else:
                                    sweepstart = -1
                                    sweepstop = -1

                            # Ok, now we should have everything ready to go to plot the meteograms.
                            # Let'er rip!
                            dis_plot_name = dis_name+'_'+DSDtype
                            pm.plotDSDmeteograms(dis_plot_name, meteogram_image_dir,
                                                 axparams, disvars, radvars.copy())
                        # Plot some derived quantities from the Parsivel
                        if(pc.plot_DSDderived):
                            PSDderiveddict = {'PSDmidtimes': PSDmidtimes, 'PSD_plot_df': PSD_plot_df}

                            pm.plotDSDderivedmeteograms(index, pc, ib, **PSDderiveddict)

                    
                    # add index == 0 because for IOP days, only want to plot for one, not for each disdrometer?
                    if(pc.plot_radar and index==0):
                        radar.plotsweeps(pc, ib, sb, sweepstart=sweepstart, sweepstop=sweepstop)

                    # Extract some needed stuff from dictionary containing disdrometer-derived
                    # polarimetric variables (note, this can be changed if you are interested
                    # instead to compare values from the gamma or exponential fits). Just change
                    # "dualpol_dis" to "dualpol_gam" or "dualpol_exp" below.
                    Zh = dualpol_dis['ZH']
                    ZDR = dualpol_dis['ZDR']
                    dBZ = dualpol_dis['dBZ']
                    RHV = dualpol_dis['RHV']
                    d = dualpol_dis['d']
                    fa2 = dualpol_dis['fa2']
                    fb2 = dualpol_dis['fb2']
                    intv = dualpol_dis['intv']

                    lamda_gam = lamda_gam/1000.
                    lamda_tmf = lamda_tmf/1000.

                ### Added for shape-slope relation plots
                    mu.extend(mu_gam)
                    lamda.extend(lamda_gam)
                    muTMF.extend(mu_tmf)
                    lamTMF.extend(lamda_tmf)
    #
                ### Interpolate radar values to match disdrometer times before starting retrieval
                    dBZ_rad = pd.Series(N.array(radvars['dBZ']), index=radvars['radmidtimes'])
                    ZDR_rad = pd.Series(N.array(radvars['ZDR']), index=radvars['radmidtimes'])
                    idx = dBZ_rad.index.union(PSDmidtimes)
                    dBZ_rad = N.array(dBZ_rad.reindex(idx).interpolate(method='index',limit=8).reindex(PSDmidtimes))
                    idx2 = ZDR_rad.index.union(PSDmidtimes)
                    ZDR_rad = N.array(ZDR_rad.reindex(idx2).interpolate(method='index',limit=8).reindex(PSDmidtimes))

                    precip_mask = (ZDR_rad > 4.) | (RHV_rad < 0.9) | (ZDR > 4.) | (RHV < 0.9)
                    dBZ_rad = ma.masked_where(precip_mask, dBZ_rad)
                    ZDR_rad = ma.masked_where(precip_mask, ZDR_rad)
                    dBZ_dis = ma.masked_where(precip_mask, dBZ)
                    ZDR_dis = ma.masked_where(precip_mask, ZDR)

                    radrettuple = DR.retrieve_DSD(dBZ_rad, ZDR_rad, d, fa2, fb2, intv, ib.wavelength)
                    (R_rad_retr, D0_rad_retr, mu_rad_retr, lam_rad_retr, N0_rad_retr, Nt_rad_retr,
                     W_rad_retr, sigm_rad_retr, Dm_rad_retr, N_rad_retr) = radrettuple
                    disrettuple = DR.retrieve_DSD(dBZ_dis, ZDR_dis, d, fa2, fb2, intv, ib.wavelength)
                    (R_dis_retr, D0_dis_retr, mu_dis_retr, lam_dis_retr, N0_dis_retr, Nt_dis_retr,
                     W_dis_retr, sigm_dis_retr, Dm_dis_retr, N_dis_retr) = disrettuple

                    lam_dis_retr = N.array(lam_dis_retr)
                    R_dis_retr = ma.masked_where(lam_dis_retr > 20., R_dis_retr)
                    D0_dis_retr = ma.masked_where(lam_dis_retr > 20., D0_dis_retr)
                    mu_dis_retr = ma.masked_where(lam_dis_retr > 20., mu_dis_retr)
                    N0_dis_retr = ma.masked_where(lam_dis_retr > 20., N0_dis_retr)
                    Nt_dis_retr = ma.masked_where(lam_dis_retr > 20., Nt_dis_retr)
                    W_dis_retr = ma.masked_where(lam_dis_retr > 20., W_dis_retr)
                    sigm_dis_retr = ma.masked_where(lam_dis_retr > 20., sigm_dis_retr)
                    Dm_dis_retr = ma.masked_where(lam_dis_retr > 20., Dm_dis_retr)
                    N_dis_retr = ma.masked_where(lam_dis_retr > 20., N_dis_retr)
                    lam_dis_retr = ma.masked_where(lam_dis_retr > 20., lam_dis_retr)

                    lam_rad_retr = N.array(lam_rad_retr)
                    R_rad_retr = ma.masked_where(lam_rad_retr > 20., R_rad_retr)
                    D0_rad_retr = ma.masked_where(lam_rad_retr > 20., D0_rad_retr)
                    mu_rad_retr = ma.masked_where(lam_rad_retr > 20., mu_rad_retr)
                    N0_rad_retr = ma.masked_where(lam_rad_retr > 20., N0_rad_retr)
                    Nt_rad_retr = ma.masked_where(lam_rad_retr > 20., Nt_rad_retr)
                    W_rad_retr = ma.masked_where(lam_rad_retr > 20., W_rad_retr)
                    sigm_rad_retr = ma.masked_where(lam_rad_retr > 20., sigm_rad_retr)
                    Dm_rad_retr = ma.masked_where(lam_rad_retr > 20., Dm_rad_retr)
                    N_rad_retr = ma.masked_where(lam_rad_retr > 20., N_rad_retr)
                    lam_rad_retr = ma.masked_where(lam_rad_retr > 20., lam_rad_retr)


                    sigm_dis_obs = N.sqrt((mu_gam + 4.)/(lamda_gam**2.)) ### equation for spectrum width from Vivekanandan et al. 2004
                    sigm_obs.extend(sigm_dis_obs)
                    sigm_dis.extend(sigm_dis_retr)
                    sigm_rad.extend(sigm_rad_retr)
                    Dm_obs.extend(D_m_disd) # Shouldn't this be Dm_dis?
                    Dm_dis.extend(Dm_dis_retr)
                    Dm_rad.extend(Dm_rad_retr)
                    N_retr=N.array(N_dis_retr)
                    N_retr=N_retr.T

                    N_retr2=N.array(N_rad_retr)
                    N_retr2=N_retr2.T

    #

                    if(pc.plot_DSDs):
                        if (not os.path.exists(ib.image_dir + 'DSDs/' + dis_name)):
                            os.makedirs(ib.image_dir + 'DSDs/' + dis_name)

                        maxd_idx = d.size

                        axdict = {'times': PSDtimestamps, 'xbin_left': min_diameter[:maxd_idx],
                                  'xbin_mid': avg_diameter[:maxd_idx], 'xbin_right': max_diameter[:maxd_idx],
                                  'xlim': (0.0, 9.0), 'ylim': (10.**2., 10.**8.5), 'interval': int(DSD_interval),
                                  'dis_name': dis_name}

                        for t in range(N.size(ND, axis=1)):

                            axdict['time'] = t

                            PSDdict = {'ND': ND[:maxd_idx, t]}

                            PSDfitdict = {'Exponential': (ND_expDSD[:maxd_idx, t], 'Exp'),
                                          'Gamma': (ND_gamDSD[:maxd_idx, t], 'Gamma'),
                                          'TMF': (ND_tmfDSD[:maxd_idx, t], 'TMF'),
                                          'Dis Retr': (N_retr[:maxd_idx, t], 'Dis Retr'),
                                          'Rad Retr': (N_retr2[:maxd_idx, t], 'Rad Retr')}

                            PSDparamdict = {'Shape_gam': (mu_gam[t], r'$\mu$ (gamma)'),
                                            'Shape_rad': (mu_rad_retr[t], r'$\mu$ (radar)'),
                                            'Slope_gam': (lamda_gam[t], r'$\lambda$ (gamma)'),
                                            'Slope_rad': (lam_rad_retr[t], r'$\lambda$ (radar)'),
                                            'D0_dis': (D_med_disd[t], r'$D_0$ (obs, mm)'),
                                            'D0_gam': (D_med_gam[t], r'$D_0$ (gamma, mm)'),
                                            'D0_rad': (D0_rad_retr[t], r'$D_0$ (radar, mm)'),
                                            'dBZ_dis': (refl_disd[t], r'Z (dis, dBZ)'),
                                            'dBZ_rad': (dBZ_rad[t], r'Z (radar, dBZ)'),
                                            'ZDR_dis': (ZDR_dis[t], r'$Z_{DR}$ (dis, dB)'),
                                            'ZDR_rad': (ZDR_rad[t], r'$Z_{DR}$ (radar, dB)'),
                                            'RR_disd': (rainrate[t],
                                                        r'Rain rate (obs, mm hr$^{-1}$)'),
                                            'Particle count': (PSD_df['pcount2'].values[t],
                                                               r'Particle count (QC)')}

                            sum = N.ma.sum(ND[:maxd_idx, t])
                            if(sum > 0.0):
                                pm.plot_DSD(ib, axdict, PSDdict, PSDfitdict, PSDparamdict)

                    if (not os.path.exists(ib.image_dir+'scattergrams/')):
                        os.makedirs(ib.image_dir+'scattergrams/')

        ###     RELATIONS FROM CAO ET AL. 2008


                    Zh_rad = pow(10.,dBZ_rad/10)

                    ###     Radar measured
                    Nt_rad_emp,D0_rad_emp,W_rad_emp,R_rad_emp = em.empirical(Zh_rad,ZDR_rad)

                    ###     Disdrometer measured
                    Nt_dis_emp,D0_dis_emp,W_dis_emp,R_dis_emp = em.empirical(Zh,ZDR_dis)


                    pcounts = PSD_df['pcount2'].values

        ####     Timeseries and Figure 9a-c from Cao et al. and Figure 8a-c from Cao et al.
                    if(pc.plot_retrieved and len(raintimes) > 0.0):
                        name = 'mu_retr'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'mu']}
                        em.dis_retr_timeseries(mu_tmf,mu_rad_retr,mu_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        name = 'lam_retr'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'lambda']}
                        em.dis_retr_timeseries(lamda_tmf,lam_rad_retr,lam_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        name = 'dBZ'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'dBZ']}
                        em.zh_zdr_timeseries(dBZ_rad,dBZ_dis,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        name = 'ZDR'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'ZDR']}
                        em.zh_zdr_timeseries(ZDR_rad,ZDR_dis,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        name = 'R'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[10**0.,10**2.25]],'axeslabels':[pc.timelabel,r'Rain Rate']}
                        em.retr_timeseries(rainrate,rainrate_tmf,R_rad_retr,R_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        em.one2one(rainrate/Zh,rainrate_tmf/Zh,R_dis_retr/Zh,R_rad_retr/Zh_rad,ib.image_dir,dis_name,name,rainrate)
                        em.scatters(N.log10(rainrate/Zh),N.log10(rainrate_tmf/Zh),N.log10(R_dis_retr/Zh),N.log10(R_rad_retr/Zh_rad),ZDR_rad,ZDR,PSDmidtimes,ib.image_dir,dis_name,name)

                        name = 'D0'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,5.0]],'axeslabels':[pc.timelabel,r'D0']}
                        em.retr_timeseries(D_med_disd,D_med_gam,D0_rad_retr,D0_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        em.one2one(D_med_disd,D_med_gam,D0_dis_retr,D0_rad_retr,ib.image_dir,dis_name,name,rainrate)
                        em.scatters(D_med_disd,D_med_gam,N.array(D0_dis_retr),D0_rad_retr,ZDR_rad,ZDR,PSDmidtimes,ib.image_dir,dis_name,name)

                        name = 'Nt'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[10**1.,10**5.]],'axeslabels':[pc.timelabel,r'Nt']}
                        em.retr_timeseries(M0,Ntr_tmf,Nt_rad_retr,Nt_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        em.one2one(M0/Zh,Ntr_tmf/Zh,Nt_dis_retr/Zh,Nt_rad_retr/Zh_rad,ib.image_dir,dis_name,name,rainrate)
                        em.scatters(N.log10(M0/Zh),N.log10(Ntr_tmf/Zh),N.log10(Nt_dis_retr/Zh),N.log10(Nt_rad_retr/Zh_rad),ZDR_rad,ZDR,PSDmidtimes,ib.image_dir,dis_name,name)

                        name = 'W'
                        axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,8.0]],'axeslabels':[pc.timelabel,r'LWC']}
                        em.retr_timeseries(LWC_disd,LWC_tmf,W_rad_retr,W_dis_retr,pstartindex,pstopindex,PSDmidtimes,axparamdict1,ib.image_dir,dis_name,name)

                        em.one2one(LWC_disd/Zh,LWC_tmf/Zh,W_dis_retr/Zh,W_rad_retr/Zh_rad,ib.image_dir,dis_name,name,rainrate)
                        em.scatters(N.log10(LWC_disd/Zh),N.log10(LWC_tmf/Zh),N.log10(W_dis_retr/Zh),N.log10(W_rad_retr/Zh_rad),ZDR_rad,ZDR,PSDmidtimes,ib.image_dir,dis_name,name)

                    if (pc.plot_outer):
                        #collect variables in order to compute statistics of entire data set 
                        R_mm.extend(rainrate_tmf)
                        R_retr.extend(R_rad_retr)
                        D0_mm.extend(D_med_gam)
                        D0_retr.extend(D0_rad_retr)
                        Nt_obs.extend(M0)
                        Nt_mm.extend(Ntr_tmf)
                        Nt_retr.extend(Nt_rad_retr)
                        W_obs.extend(LWC_disd)
                        W_mm.extend(LWC_tmf)
                        W_retr.extend(W_rad_retr)
                        Zh_dis.extend(Zh)
                        Zh_retr.extend(Zh_rad)
                        R_list.extend(rainrate)
                        D0_list.extend(D_med_disd)
                        ND = ND.T
                        for t in range(N.size(ND,axis=0)):
                            ND_list = N.append(ND_list, [ND[t,:]], axis=0)

#
#
# ### not using right now...ignore...need to do with all cases
# #               fig1=plt.figure(figsize=(8,8))
# #               ax1=fig1.add_subplot(111)
# #               ax1.scatter(mu_gam, rainrate, color='m', marker='.', label='Method Moments')
# #               ax1.scatter(mu_dis_retr, R_dis_retr, color='c', marker='.', label='Dis Retrieval')
# #               ax1.scatter(mu_rad_retr, R_rad_retr, color='k', marker='.', label='Rad Retrieval')
# #               ax1.set_xlim(-2.0,15.0)
# #               ax1.set_yscale('log')
# #               ax1.set_ylim(10**-1,10**2)
# #               ax1.set_xlabel('Mu')
# #               ax1.set_ylabel('RainRate')
# #               plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
# #               plt.savefig(ib.image_dir+'scattergrams/'+dis_name+'_mu_R.png',dpi=200,bbox_inches='tight')
# #               plt.close(fig1)
# #
# #               fig1=plt.figure(figsize=(8,8))
# #               ax1=fig1.add_subplot(111)
# #               ax1.scatter(lamda_gam, rainrate, color='m', marker='.', label='Method Moments')
# #               ax1.scatter(lam_dis_retr, R_dis_retr, color='c', marker='.', label='Dis Retrieval')
# #               ax1.scatter(lam_rad_retr, R_rad_retr, color='k', marker='.', label='Rad Retrieval')
# #               ax1.set_xlim(0.0,15.0)
# #               ax1.set_yscale('log')
# #               ax1.set_ylim(10**-1,10**2)
# #               ax1.set_xlabel('Lamda')
# #               ax1.set_ylabel('RainRate')
# #               plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
# #               plt.savefig(ib.image_dir+'scattergrams/'+dis_name+'_lam_R.png',dpi=200,bbox_inches='tight')
# #               plt.close(fig1)
# #
# #               fig1=plt.figure(figsize=(8,8))
# #               ax1=fig1.add_subplot(111)
# #               ax1.scatter(D_med_gam, rainrate, color='m', marker='.', label='Method Moments')
# #               ax1.scatter(D0_dis_retr, R_dis_retr, color='c', marker='.', label='Dis Retrieval')
# #               ax1.scatter(D0_rad_retr, R_rad_retr, color='k', marker='.', label='Rad Retrieval')
# #               ax1.set_xlim(0.5,3.)
# #               ax1.set_yscale('log')
# #               ax1.set_ylim(10**-1,10**2)
# #               ax1.set_xlabel('D0')
# #               ax1.set_ylabel('RainRate')
# #               plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
# #               plt.savefig(ib.image_dir+'scattergrams/'+dis_name+'_D0_R.png',dpi=200,bbox_inches='tight')
# #               plt.close(fig1)
# #
# #               D0dict[dis_name+'_obs'] = D_med_disd
# #               D0dict[dis_name+'_mom'] = D_med_gam
# #               D0dict[dis_name+'_dis_retr'] = D0_dis_retr
# #               D0dict[dis_name+'_rad_retr'] = D0_rad_retr
# #               ZDRdict[dis_name+'_dis'] = ZDR
# #               ZDRdict[dis_name+'_rad'] = ZDR_rad
# #               Zhdict[dis_name+'_dis'] = Zh
# #               Zhdict[dis_name+'_rad'] = Zh_rad
# #               dBZdict[dis_name+'_] = dBZ
# #               Wdict[dis_name+'_obs'] = LWC_disd
# #               Rdict[dis_name+'_obs'] = intensities
# #               Ntdict[dis_name+'_obs'] = M0
# #
# #               name = 'D0'
# #               ymin = 0.0
# #               ymax = 5.0
# #               ylabel = 'D0'
# #               em.PIPS(D0dict['PIPS_1A_obs'],D0dict['PIPS_1B_obs'],D0dict['PIPS_2A_obs'],D0dict['PIPS_2B_obs'],ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,ib.image_dir,dis_name,name,ylabel)
# #
# #               name = 'W'
# #               ymin = -6.0
# #               ymax = -1.0
# #               ylabel = 'log(W/Zh)'
# #               em.PIPS(N.log10(Wdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Wdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Wdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Wdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,ib.image_dir,dis_name,name,ylabel)
# #
# #               name = 'R'
# #               ymin = -5.0
# #               ymax = 0.0
# #               ylabel = 'log(R/Zh)'
# #               em.PIPS(N.log10(Rdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Rdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Rdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Rdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,ib.image_dir,dis_name,name,ylabel)
# #
# #               name = 'Nt'
# #               ymin = -4.0
# #               ymax = 2.0
# #               ylabel = 'log(Nt/Zh)'
# #               em.PIPS(N.log10(Ntdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Ntdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Ntdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Ntdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,ib.image_dir,dis_name,name,ylabel)
# #
#
# ###    not using right now Figure 2 from Brandes et al. 2004
# #
# #           one_x = N.linspace(0.0,60.0)
# #           upper_y = N.exp(1.01*10**-4*one_x**3 - 7.09*10**-3*one_x**2 + 2.38*10**-1*one_x - 3.44)
# #           lower_y = N.exp(2.12*10**-4*one_x**2 + 6.48*10**-2*one_x - 3.87)
# #
# #           fig1=plt.figure(figsize=(8,8))
# #           ax1=fig1.add_subplot(111)
# #           ax1.scatter(ZDRdict['PIPS_1B'], ZDRdict['PIPS_1B_obs'], marker='.', label='PIPS 1B')
# #           ax1.scatter(ZDRdict['PIPS_2A'], ZDRdict['PIPS_2A_obs'], marker='.', label='PIPS 2A')
# #           ax1.scatter(ZDRdict['PIPS_2B'], ZDRdict['PIPS_2B_obs'], marker='.', label='PIPS 2B')
# #           ax1.set_xlim(0.0,60.0)
# #           ax1.set_ylim(-1.0,4.0)
# #           ax1.set_xlabel('Reflectivity (dBZ)')
# #           ax1.set_ylabel('ZDR (dB)')
# #           ax1.plot(one_x,upper_y,color='k')
# #           ax1.plot(one_x,lower_y,color='k')
# #           plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
# #           plt.savefig(ib.image_dir+'scattergrams/brandes.png',dpi=200,bbox_inches='tight')
# #           plt.close(fig1)
#

if (pc.plot_outer):
    name = 'R'
    em.outer_one2one(N.array(R_list)/N.array(Zh_dis),N.array(R_mm)/N.array(Zh_dis),N.array(R_retr)/N.array(Zh_retr),outer_image_dir,name,R_list)

    name = 'D0'
    em.outer_one2one(N.array(D0_list), N.array(D0_mm), N.array(D0_retr), outer_image_dir,name,R_list)

    name = 'Nt'
    em.outer_one2one(N.array(Nt_obs)/N.array(Zh_dis),N.array(Nt_mm)/N.array(Zh_dis),N.array(Nt_retr)/N.array(Zh_retr), outer_image_dir,name,R_list)

    name = 'W'
    em.outer_one2one(N.array(W_obs)/N.array(Zh_dis),N.array(W_mm)/N.array(Zh_dis),N.array(W_retr)/N.array(Zh_retr),outer_image_dir,name,R_list)

    ## Plot the untruncated lambda-mu relation and fit with 2nd order polynomial
    lamda = N.array(lamda)
    mu = N.array(mu)
    lamda = lamda[~N.isnan(lamda)]
    mu = mu[~N.isnan(mu)]
    poly=N.polynomial.polynomial.polyfit(lamda,mu,2)
    polynomial=N.polynomial.polynomial.Polynomial(poly)

    xx = N.linspace(0.0, 30.0)
    yy = polynomial(xx)
    y2 = -0.0201*xx**2. + 0.902*xx - 1.718
    y3 = -0.016*xx**2. + 1.213*xx - 1.957

    fig=plt.figure(figsize=(8,8))
    ax1=fig.add_subplot(111)
    plt.title('Untruncated Shape-Slope Relation')
    ax1.scatter(lamda,mu, color='k', marker='.')
    ax1.plot(xx,yy,label='Our Relation')
    ax1.plot(xx,y2,label='Cao Relation')
    ax1.plot(xx,y3,label='Zhang Relation')
    ax1.set_xlim(0.0,20.0)
    ax1.set_ylim(-5.0,20.0)
    ax1.set_xlabel('Slope parameter')
    ax1.set_ylabel('Shape parameter')
    ax1.text(0.05,0.90,'# of Points: %2.1f'%len(lamda), transform=ax1.transAxes, fontsize=12.)
    ax1.text(0.05,0.85,'%2.4f'%poly[2]+'*lam^2 + %2.4f'%poly[1]+'*lam + %2.4f'%poly[0], transform=ax1.transAxes, fontsize=12.)
    ax1.text(0.05,0.80,'Average mu: %2.2f'%N.mean(mu), transform=ax1.transAxes, fontsize = 12.)
    ax1.text(0.05,0.75,'Average lambda: %2.2f'%N.mean(lamda), transform=ax1.transAxes, fontsize = 12.)
    plt.legend(loc='upper left',numpoints=1,ncol=3,fontsize=12.)
    plt.savefig(outer_image_dir+'shape_slope.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

    ## Plot the truncated lambda-mu relation and fit with 2nd order polynomial
    print len(lamTMF)
    print len(muTMF)
    lamTMF = N.array(lamTMF)
    muTMF = N.array(muTMF)
    lamTMF = lamTMF[~N.isnan(lamTMF)]
    muTMF = muTMF[~N.isnan(muTMF)]
    print len(lamTMF)
    print len(muTMF)
    poly2=N.polynomial.polynomial.polyfit(lamTMF,muTMF,2)
    polynomial2=N.polynomial.polynomial.Polynomial(poly2)

    yyTMF = polynomial2(xx)

    fig=plt.figure(figsize=(8,8))
    ax1=fig.add_subplot(111)
    plt.title('Truncated Shape-Slope Relation')
    ax1.scatter(lamTMF,muTMF, color='k', marker='.')
    ax1.plot(xx,yyTMF,label='Our Relation')
    ax1.plot(xx,y2,label='Cao Relation')
    ax1.plot(xx,y3,label='Zhang Relation')
    ax1.set_xlim(0.0,20.0)
    ax1.set_ylim(-5.0,20.0)
    ax1.set_xlabel('Slope parameter')
    ax1.set_ylabel('Shape parameter')
    ax1.text(0.05,0.90,'# of Points: %2.1f'%len(lamTMF), transform=ax1.transAxes, fontsize=12.)
    ax1.text(0.05,0.85,'%2.4f'%poly2[2]+'*lam^2 + %2.4f'%poly2[1]+'*lam + %2.4f'%poly2[0], transform=ax1.transAxes, fontsize=12.)
    ax1.text(0.05,0.80,'Average mu: %2.2f'%N.mean(muTMF), transform=ax1.transAxes, fontsize = 12.)
    ax1.text(0.05,0.75,'Average lambda: %2.2f'%N.mean(lamTMF), transform=ax1.transAxes, fontsize = 12.)
    plt.legend(loc='upper left',numpoints=1,ncol=3,fontsize=12.)
    plt.savefig(outer_image_dir+'TMF_shape_slope.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

    #######################
    Dm_obs = N.array(Dm_obs)
    Dm_rad = N.array(Dm_rad)
    Dm_dis = N.array(Dm_dis)
    sigm_obs = N.array(sigm_obs)
    sigm_rad = N.array(sigm_rad)
    sigm_dis = N.array(sigm_dis)

    ### additional figures from Cao et al 2008
    fig=plt.figure(figsize=(8,8))
    ax1=fig.add_subplot(111)
    plt.title('D_m-Sigm_m Relation with Slope < 20')
    ax1.scatter(Dm_obs,sigm_obs, color='k', marker='.',label='Dis Obs')
    ax1.scatter(Dm_rad,sigm_rad, color='b', marker = '.',label='Radar Retr')
    #ax1.set_xlim(0.0,5.0)
    #ax1.set_ylim(0.0,5.0)
    ax1.set_xlabel('Dm , mm')
    ax1.set_ylabel('sigma m , mm')
    plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=10.)
    plt.savefig(outer_image_dir+'dm_sigma.png',dpi=200,bbox_inches='tight')
    plt.close(fig)

    fig1=plt.figure(figsize=(8,8))
    ax1=fig1.add_subplot(111)
    one_x = N.linspace(0.0,6.0)
    one_y = one_x
    bias_dis = 100 * ((N.nansum(sigm_dis-sigm_obs))/N.nansum(sigm_obs))
    bias_rad = 100 * ((N.nansum(sigm_rad-sigm_obs))/N.nansum(sigm_obs))
    cc_dis = pd.DataFrame({'dis': sigm_dis, 'obs': sigm_obs}).corr()
    cc_rad = pd.DataFrame({'rad': sigm_rad, 'obs': sigm_obs}).corr()
    ax1.scatter(sigm_obs, sigm_rad, color='g', marker='.', label='Rad Retrieval')
    ax1.scatter(sigm_obs, sigm_dis, color='c', marker='.', label='Dis Retreival')
    #ax1.set_xlim(0.0,4.0)
    #ax1.set_ylim(0.0,4.0)
    ax1.set_xlabel('Observed Sigma')
    ax1.set_ylabel('Calculated Sigma')
    ax1.plot(one_x,one_y,lw=2,color='k')
    ax1.text(0.6,0.20,'Dis Retr. Bias =%2.2f'%bias_dis+'%',transform=ax1.transAxes)
    ax1.text(0.6,0.15,'Rad Retr. Bias =%2.2f'%bias_rad+'%',transform=ax1.transAxes)
    ax1.text(0.6,0.10,'Dis Retr. Corr Coeff =%2.3f'%cc_dis.iloc[0,1],transform=ax1.transAxes)
    ax1.text(0.6,0.05,'Rad Retr. Corr Coeff =%2.3f'%cc_rad.iloc[0,1],transform=ax1.transAxes)
    plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=12.)
    plt.savefig(outer_image_dir+'sigma_one2one.png',dpi=200,bbox_inches='tight')
    plt.close(fig1)

    fig1=plt.figure(figsize=(8,8))
    ax1=fig1.add_subplot(111)
    one_x = N.linspace(0.0,6.0)
    one_y = one_x
    bias_dis = 100 * ((N.nansum(Dm_dis-Dm_obs))/N.nansum(Dm_obs))
    bias_rad = 100 * ((N.nansum(Dm_rad-Dm_obs))/N.nansum(Dm_obs))
    cc_dis = pd.DataFrame({'dis': Dm_dis, 'obs': Dm_obs}).corr()
    cc_rad = pd.DataFrame({'rad': Dm_rad, 'obs': Dm_obs}).corr()
    ax1.scatter(Dm_obs, Dm_rad, color='g', marker='.', label='Rad Retrieval')
    ax1.scatter(Dm_obs, Dm_dis, color='c', marker='.', label='Dis Retreival')
    #ax1.set_xlim(0.0,6.0)
    #ax1.set_ylim(0.0,6.0)
    ax1.set_xlabel('Observed Dm')
    ax1.set_ylabel('Calculated Dm')
    ax1.plot(one_x,one_y,lw=2,color='k')
    ax1.text(0.6,0.20,'Dis Retr. Bias =%2.2f'%bias_dis+'%',transform=ax1.transAxes)
    ax1.text(0.6,0.15,'Rad Retr. Bias =%2.2f'%bias_rad+'%',transform=ax1.transAxes)
    ax1.text(0.6,0.10,'Dis Retr. Corr Coeff =%2.3f'%cc_dis.iloc[0,1],transform=ax1.transAxes)
    ax1.text(0.6,0.05,'Rad Retr. Corr Coeff =%2.3f'%cc_rad.iloc[0,1],transform=ax1.transAxes)
    plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=12.)
    plt.savefig(outer_image_dir+'Dm_one2one.png',dpi=200,bbox_inches='tight')
    plt.close(fig1)

    ### Packages variables to be used in SATP.py
    R_tarr = N.array(R_list)
    D0_tarr = N.array(D0_list)

    radnpz_filename = 'SATP_all.npz'
    savevars={}
    savevars['Nc_bin'] = ND_list
    savevars['R'] = R_tarr
    savevars['D0'] = D0_tarr
    N.savez(radnpz_filename,**savevars)

#     #### Experimental plots 
#     fig1=plt.figure(figsize=(8,8))
#     ax1=fig1.add_subplot(111)
#     one_x = N.linspace(-100.0,100.0)
#     one_y = one_x
#     ax1.scatter(R_list, D0_list-D0_retr, marker='.')
#     ax1.set_xlabel('Rain Rate (mm/hr)')
#     ax1.set_ylabel('Observed D0 - Retrieved D0')
#     ax1.plot(one_x,one_y,lw=2,color='k')
#     plt.savefig(outer_image_dir+'D0_vs_R.png',dpi=200,bbox_inches='tight')
#     plt.close(fig1)
#     
#     fig1=plt.figure(figsize=(8,8))
#     ax1=fig1.add_subplot(111)
#     one_x = N.linspace(-100.0,100.0)
#     one_y = one_x
#     ax1.scatter(Zh, D0_list-D0_retr, marker='.')
#     ax1.set_xlabel('Zh (dBZ)')
#     ax1.set_ylabel('Observed D0 - Retrieved D0')
#     ax1.plot(one_x,one_y,lw=2,color='k')
#     plt.savefig(outer_image_dir+'D0_vs_R.png',dpi=200,bbox_inches='tight')
#     plt.close(fig1)

