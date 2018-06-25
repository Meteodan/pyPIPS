# pyPIPS.py
#
# This script analyzes and plots data from the Portable Integrated Precipitation Stations (PIPS)

import numpy as N
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as dates
from datetime import timedelta
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

fieldnames = ['dBZ', 'ZDR', 'KDP', 'RHV', 'Vr']

min_diameter = dis.min_diameter
max_diameter = dis.max_diameter
bin_width = max_diameter - min_diameter
avg_diameter = dis.avg_diameter
fall_bins = dis.fall_bins
min_fall_bins = dis.min_fall_bins

# Create V-D relationship for rain based on Terry Schuur's relationship
rainvd = dis.assignfallspeed(avg_diameter)

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

# Parse command line argument and read in disdrometer/radar information from text file
# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
# version first, and if not, reading the textfile and pickling the read-in
# variables for future read-ins.  Think about this.

if(len(sys.argv) == 1):
    argindex = 1
elif(len(sys.argv) > 1):
    argindex = 2
else:
    sys.exit("No text input file defined! Quitting!")

ib = utils.readpyPIPSinput(sys.argv[argindex])

if (not os.path.exists(ib.image_dir)):
    os.makedirs(ib.image_dir)

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

# ------
# Grab radar data for comparison if desired

if(pc.comp_radar):
    sb = radar.readsweeps2PIPS(fieldnames, pc, ib)

    if(pc.plot_radar):
        radar.plotsweeps(pc, ib, sb)

# Outer disdrometer (and deployment) loop
mu = []
lamda = []
Mu_retr = []
Lam_retr = []
D0dict = {}
ZDRdict = {}
Wdict = {}
Rdict = {}
Ntdict = {}

for index, dis_filename, dis_name, starttime, stoptime, centertime, dloc in zip(xrange(
        0, len(ib.dis_list)), ib.dis_list, ib.dis_name_list, ib.starttimes, ib.stoptimes,
        ib.centertimes, ib.dlocs):

    # Is this necessary?
    if(pc.timetospace):
        centertime = ib.centertimes[index]

    if(pc.plot_conv_meteo or pc.plot_DSD_meteo):
        # Create the directory for the meteogram plots if it doesn't exist
        meteogram_image_dir = ib.image_dir + 'meteograms/'
        if (not os.path.exists(meteogram_image_dir)):
            os.makedirs(meteogram_image_dir)

    # Read in the disdrometer data file
    filepath = os.path.join(ib.dis_dir, dis_filename)

    PIPS_dict = dis.readPIPS(filepath, basicqc=pc.basicQC, rainfallqc=pc.rainfallQC,
                             rainonlyqc=pc.rainonlyQC, hailonlyqc=pc.hailonlyQC,
                             strongwindqc=pc.strongwindQC, DSD_interval=pc.DSD_interval,
                             starttime=starttime, stoptime=stoptime)

    # Unpack some stuff from the PIPS_dict
    onesectimestamps = PIPS_dict['onesectimestamps']
    PSDtimestamps = PIPS_dict['PSDtimestamps']
    onesec_df = PIPS_dict['onesec_df']
    PSD_df = PIPS_dict['PSD_df']
    ND = PIPS_dict['ND']
    ND_onedrop = PIPS_dict['ND_onedrop']
    countsMatrix = PIPS_dict['countsMatrix']

    ND = ND.T
    ND_onedrop = ND_onedrop.T

    logND = N.ma.log10(ND)
    # Commented out below: we want to do the "one-drop" masking for the gamma and
    # exponential fits, *not* for the actual data. (You can't observe a fraction of a drop! and the
    # way it's computed (using the theoretical fall speed instead of measured) may actually end up
    # zeroing out some bins.)
    #logND = N.ma.masked_where(ND < ND_onedrop, logND)

    DSD_index = PIPS_dict['DSD_index']
    DSD_interval = PIPS_dict['DSD_interval']
    DSD_interval_td = timedelta(seconds=DSD_interval)
    DSD_halfinterval_td = timedelta(seconds=DSD_interval / 2.)

    # Compute potential temperature, water vapor mixing ratio, and density

    onesec_df['pt'] = thermo.caltheta(onesec_df['pressure']*100., onesec_df['fasttemp']+273.15)
    onesec_df['qv'] = thermo.calqv(onesec_df['RH_derived']/100., onesec_df['pressure']*100.,
                                   onesec_df['fasttemp']+273.15)
    onesec_df['rho'] = thermo.calrho(onesec_df['pressure']*100., onesec_df['pt'], onesec_df['qv'])

    # Determine start and end times/indices for analysis

    onesectimestampsnums = dates.date2num(onesectimestamps)
    PSDtimestampsnums = dates.date2num(PSDtimestamps)

    startindex, stopindex = utils.getTimeWindow(starttime, stoptime, onesectimestampsnums)
    pstartindex, pstopindex = utils.getTimeWindow(starttime, stoptime, PSDtimestampsnums)

    starttime = onesectimestampsnums[startindex]
    stoptime = onesectimestampsnums[stopindex]
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
                'pstoptime': pstoptime, 'onesectimestamps': onesectimestamps,
                'onesectimestampnums': onesectimestampsnums, 'PSDtimestamps': PSDtimestamps,
                'PSDtimestampsnums': PSDtimestampsnums, 'PSDtimestamps_edge': PSDtimestamps_edge,
                'PSDtimestamps_avg': PSDtimestamps_avg}

    # TODO: Address this section
    if(pc.timetospace and pc.comp_radar):
        delta_t_rad = [(x - centertime).total_seconds() for x in sb.radtimes]
        delta_t_dis = [(x - centertime).total_seconds() for x in PSDtimestamps]
        delta_t_dis_start = [(x - centertime).total_seconds()
                             for x in PSDtimestamps_edge]
        delta_t_dis_avg = [(x - centertime).total_seconds()
                           for x in PSDtimestamps_avg]
#         if(plot_thermo):
#             delta_t_thermo = [(x-centertime).total_seconds() for x in thermodates]
#             delta_t_thermo_1min = [(x-centertime).total_seconds() for x in thermodates_1min]
#             delta_t_wind = [(x-centertime).total_seconds() for x in winddates]
#             delta_t_wind_1min = [(x-centertime).total_seconds() for x in winddates_1min]

        xloc_rad = [(x * pc.stormmotion) / 1000.0 for x in delta_t_rad]
        xloc_dis = [(x * pc.stormmotion) / 1000.0 for x in delta_t_dis]
        xloc_dis_start = [(x * pc.stormmotion) /
                          1000.0 for x in delta_t_dis_start]
        xloc_dis_avg = [(x * pc.stormmotion) / 1000.0 for x in delta_t_dis_avg]
#         if(plot_thermo):
#             xloc_thermo = [(x*pc.stormmotion)/1000.0 for x in delta_t_thermo]
#             xloc_thermo_1min = [(x*pc.stormmotion)/1000.0 for x in delta_t_thermo_1min]
#             xloc_wind = [(x*pc.stormmotion)/1000.0 for x in delta_t_wind]
#             xloc_wind_1min = [(x*pc.stormmotion)/1000.0 for x in delta_t_wind_1min]

        plotx_rad = N.array(xloc_rad)
        plotx_dis = N.array(xloc_dis)
        plotx_dis_start = N.array(xloc_dis_start)
        plotx_dis_avg = N.array(xloc_dis_avg)
#         if(plot_thermo):
#             plotx_thermo = N.array(xloc_thermo)
#             plotx_thermo_1min = N.array(xloc_thermo_1min)
#             plotx_wind = N.array(xloc_wind)
#             plotx_wind_1min = N.array(xloc_wind_1min)
    else:
        if(pc.comp_radar):
            plotx_rad = dates.date2num(sb.radtimes)
        plotx_dis = dates.date2num(PSDtimestamps)
        plotx_dis_start = dates.date2num(PSDtimestamps_edge)
        plotx_dis_avg = dates.date2num(PSDtimestamps_avg)
#         if(plot_thermo):
#             plotx_thermo = dates.date2num(thermodates)
#             plotx_thermo_1min = dates.date2num(thermodates_1min)
#             plotx_wind = dates.date2num(winddates)
#             plotx_wind_1min = dates.date2num(winddates_1min)

    # Resample the 1-s data corresponding to each integrated DSD (for the
    # given DSD interval). Note, most of these aren't used for right now.
    sec_offset = PSDtimestamps[0].second

    onesec_resampled_df = dis.resampleonesec(DSD_interval, sec_offset, onesec_df)
    onesec_resampled_df = onesec_resampled_df.loc[DSD_index]
    rho_tDSD = onesec_resampled_df['rho']

    # Conventional data meteogram plotting
    if(pc.plot_conv_meteo):
        # Interval for plotting in seconds.  Setting to larger intervals is useful to avoid
        # plot clutter for long datasets, but information will be lost.
        plotinterval = 1
        plotintervalstr = '{:d}S'.format(int(plotinterval))
        sec_offset = onesectimestamps[0].second

        # Plot wind meteogram
        windavgintv = 60
        windgustintv = 3

        plottimes = onesectimestampsnums[startindex:stopindex + 1:plotinterval]
        plottimeindex = pd.DatetimeIndex(
            onesectimestamps[startindex:stopindex + 1:plotinterval])

        onesec_plot_df = onesec_df.loc[plottimeindex]

        # Resample wind diagnostics array to flag as bad any time in the interval given by
        # plotinterval. Note, have to use numpy max() as a lambda function because the
        # pandas resample.max() does not propagate NaN's!
        if(pc.plot_diagnostics):
            onesec_plot_df['winddiag'] = onesec_df['winddiag'].resample(
                plotintervalstr, label='right', closed='right', base=sec_offset,
                how=lambda x: utils.trymax(x.values)).loc[plottimeindex]

        # Organize a bunch of stuff in a dictionary
        convmeteodict = {'plotinterval': plotinterval, 'plottimes': plottimes,
                         'windavgintv': windavgintv, 'windgustintv': windgustintv,
                         'onesec_plot_df': onesec_plot_df}

        pm.plotonesecmeteograms(index, pc, ib, convmeteodict)

    # N.savetxt('ND1minavg.txt', ND1minavg) # Commented out for now

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
        ND, logND, D_med_disd, D_m_disd, D_mv_disd, D_ref_disd, QR_disd, refl_disd, LWC_disd, M0, rainrate = \
            dis_DSD

        ND_expDSD = ND_expDSD.T
        logND_expDSD = N.ma.log10(ND_expDSD / 1000.)  # Get to log(m^-3 mm^-1)
        logND_expDSD = N.ma.masked_where(ND_expDSD < ND_onedrop, logND_expDSD)

        ND_gamDSD = ND_gamDSD.T
        logND_gamDSD = N.ma.log10(ND_gamDSD / 1000.)  # Get to log(m^-3 mm^-1)
        logND_gamDSD = N.ma.masked_where(ND_gamDSD < ND_onedrop, logND_gamDSD)

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
        # Loop through the observed, exponential, and gamma fits
        for DSDtype in ['observed', 'exponential', 'gamma']:
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

            radvars = {}
            # Get radar variables together for comparison, if desired
            if(pc.comp_radar and DSDtype == 'observed'):
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
                        gc_mask = N.where((radvars['RHV'] < 0.90), True, False)
                        for radvarname in ['ZDR','dBZ','RHV']:
                                radvars[radvarname] = N.ma.masked_array(radvars[radvarname],
                                                                      mask=gc_mask)
            if(pc.plot_only_precip and pc.calc_dualpol):
                # set plot start and end time to the start and end of precipitation
                rainindex = N.where(disvars['RHV'] > 0.1)
                raintimes = DSDmidtimes[rainindex]
                plotstarttime = raintimes[0]
                plotstoptime = raintimes[len(raintimes)-1]
            if(DSDtype == 'observed'):
                # Prepare axis parameters
                timelimits = [plotstarttime, plotstoptime]
                try:
                    diamlimits = pc.DSD_D_range
                    diamytick = pc.DSD_D_ytick
                except:
                    diamlimits = [0.0, 9.0]
                    diamytick = 1.0

                axparams = {'majorxlocator': pc.locator, 'majorxformatter': pc.formatter,
                            'minorxlocator': pc.minorlocator,
                            'axeslimits': [timelimits, diamlimits],
                            'majorylocator': ticker.MultipleLocator(base=diamytick),
                            'axeslabels': [None, 'D (mm)']}

            # Ok, now we should have everything ready to go to plot the meteograms.
            # Let'er rip!
            dis_plot_name = dis_name+'_'+DSDtype
            pm.plotDSDmeteograms(dis_plot_name, meteogram_image_dir,
                                 axparams, disvars, radvars.copy())

        # Plot some derived quantities from the Parsivel
        if(pc.plot_DSDderived):
            PSDderiveddict = {'PSDmidtimes': PSDmidtimes, 'PSD_plot_df': PSD_plot_df}

            pm.plotDSDderivedmeteograms(index, pc, ib, **PSDderiveddict)

    if(pc.plot_DSDs):
        if (not os.path.exists(ib.image_dir + 'DSDs/' + dis_name)):
            os.makedirs(ib.image_dir + 'DSDs/' + dis_name)

        axdict = {'times': PSDtimestamps, 'xbin_left': min_diameter,
                  'xbin_mid': avg_diameter, 'xbin_right': max_diameter,
                  'xlim': (0.0, 9.0), 'ylim': (10.**2., 10.**8.5), 'interval': int(DSD_interval),
                  'dis_name': dis_name}

        for t in range(N.size(ND, axis=1)):

            axdict['time'] = t

            PSDdict = {'ND': ND[:, t]}

            PSDfitdict = {'Exponential': (ND_expDSD[:, t], 'Exp'),
                          'Gamma': (ND_gamDSD[:, t], 'Gamma')}

            PSDparamdict = {'Shape_gam': (mu_gam[t], r'$\mu$ (gamma)'),
                            'Slope_gam': (lamda_gam[t], r'$\lambda$ (gamma)'),
                            'D0_gam': (D_med_gam[t], r'$D_0$ (gamma, mm)'),
                            'Slope_exp': (lamda_exp[t], r'$\lambda$ (exp)'),
                            'N0_exp': (N0_exp[t], r'$N_0$ (exp, m$^{-4}$)'),
                            'D0_dis': (D_med_disd[t], r'$D_0$ (obs, mm)'),
                            'dBZ_dis': (refl_disd[t], r'Z (obs, dBZ)'),
                            'RR_disd': (PSD_df['intensity'].values[t],
                                        r'Rain rate (obs, mm hr$^{-1}$)'),
                            'Particle count': (PSD_df['pcount2'].values[t],
                                               r'Particle count (QC)')}

            if(pc.calc_dualpol):
                PSDparamdict['ZDR_dis'] = (dualpol_dis['ZDR'][t], r'$Z_{DR}$ (obs, dB)')

            sum = N.ma.sum(ND[:, t])
            if(sum > 0.0):
                pm.plot_DSD(ib, axdict, PSDdict, PSDfitdict, PSDparamdict)

    if(pc.plot_vel_D):
        if (not os.path.exists(ib.image_dir + 'vel_D/' + dis_name)):
            os.makedirs(ib.image_dir + 'vel_D/' + dis_name)

        axdict = {'times': PSDtimestamps, 'min_diameter': min_diameter,
                  'avg_diameter': avg_diameter, 'min_fall_bins': min_fall_bins,
                  'xlim': pc.D_range, 'ylim': pc.vel_range,
                  'dis_name': dis_name}

        for t in range(N.size(countsMatrix, axis=0)):
            if(PSD_df['pcount2'].values[t] > 0):
                axdict['time'] = t
                PSDdict = {'countsMatrix': countsMatrix[t, :],
                           'flaggedtime': PSD_df['flaggedtimes'].values[t]}
                pm.plot_vel_D(ib, axdict, PSDdict, rho_tDSD[t])
