# Plot_Disdrometer
# Runs through all of the IOPs in the input folder and does meteograms and scattergrams and calculates overall lam-mu relation
# This script plots several quantities based on disdrometer data from the Parsivel laser disdrometer
# command line: python pyPIPS_CG.py pyPIPScontrol.py

import netCDF4 as netcdf
import Nio
import matplotlib
#matplotlib.use('TkAgg')
import numpy as N
import numpy.ma as MA
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.dates import *
from matplotlib.ticker import *
from mpl_toolkits.basemap import Basemap
from datetime import datetime,timedelta
import pytz as pytz
from scipy import special
from scipy import signal
from scipy import ndimage
import modules.thermolib as thermo
import math as M
import imp
import glob
import os
import modules.obanmodule as oban
from mpl_toolkits.axes_grid1 import make_axes_locatable,host_subplot
from matplotlib import colorbar
import modules.disdrometer_module as dis
from matplotlib.font_manager import FontProperties
import modules.raymond_lowpass as raymond
import modules.DSDlib as dsd
import pickle
#import cressman as cress
import pdb
import sys
import modules.ctablesfrompyesviewer as ctables
import modules.radarmodule as radar
import pandas as pd
import modules.plotmodule as pm
from modules.utils import log,warning,fatal
import modules.DSDretrieval as DR
import modules.empirical_module as em

k=0
mu=[]
lamda=[]

indir = '/Users/bozell/pyPIPS_work/input/'
for root, dirs, filenames in os.walk(indir):
	for f in filenames:
		if f.endswith(".txt"):
			print f 
			k = k+1
			fieldnames = ['dBZ','ZDR','KDP','RHV','Vr']
		
			# Function definition
			def trymax(a,default=0):
				try:
					return a.max()
				except:
					return default


			deg2rad = N.pi/180.


			# Set global font size for axes and colorbar labels, etc.
			font = {'size':10}
			matplotlib.rc('font',**font)
			fontP = FontProperties()
			fontP.set_size('small')


			# Disdrometer bin and diameter information 
			min_diameter = dis.min_diameter
			max_diameter = dis.max_diameter
			bin_width = max_diameter-min_diameter
			avg_diameter = dis.avg_diameter
			fall_bins = dis.fall_bins
			min_fall_bins = dis.min_fall_bins


			#Time zone stuff
			Central = pytz.timezone('US/Central')
			fmt = '%Y-%m-%d %H:%M:%S %Z%z'
			fmt2 = '%Y-%m-%d %H:%M:%S'
			fmt3 = '%Y%m%d%H%M%S'



			#-----------------------------------------------------------------------
			#
			#   Dynamically import pyPIPScontrol.py or user version of it.
			#
			#-----------------------------------------------------------------------

			def import_all_from(module_path):
				"""Modified from http://grokbase.com/t/python/python-list/1172ahxp0s/from-module-import-using-import
				   Loads python file at "module_path" as module and adds contents to global namespace."""
				#mod = __import__(module_name)
				mod = imp.load_source('mod',module_path)
				return mod

			if len(sys.argv) > 1:   # Try to import user-defined plotcontrol module
				controlmodpath = sys.argv[1]
				log("Input file is "+controlmodpath)
				pc = import_all_from(controlmodpath)
				try:
					pc = import_all_from(controlmodpath)
					log("Successfully imported pyPIPS control parameters!")
				except:
					warning("Unable to import user-defined pyPIPS control parameters! Reverting to defaults.")
					import pyPIPScontrol as pc
			else:   # Read in default plotcontrol.py
				import pyPIPScontrol as pc

			# Parse command line argument and read in disdrometer/radar information from text file
			# Maybe in the future can find a more efficient way, such as checking to see if there is a pickled
			# version first, and if not, reading the textfile and pickling the read-in variables for future read-ins.  Think about this.

			if(len(sys.argv) == 1):
				argindex = 1
			elif(len(sys.argv) > 1):
				argindex = 2
			else:
				sys.exit("No text input file defined! Quitting!")

			with open(os.path.join(root, f),'r') as inputfile:
				# Read and discard header
				inputfile.readline()
				inputfile.readline()
				inputfile.readline()
				dis_dir = inputfile.readline().strip()
				image_dir = inputfile.readline().strip()
				inputfile.readline()
				numdis = N.int(inputfile.readline().strip().split()[0]) # Read in number of disdrometers
	
				# Read in disdrometer locations and names
				dis_list = []
				dis_ftype_list = []
				tprh_filenames = []
				wind_filenames = []
				dlocs=[]
				dis_name_list = []
				starttimes = []
				stoptimes = []
				centertimes = []
				for l in xrange(numdis):
					line = inputfile.readline().strip().split(',')
					dname = line[0] # Disdrometer name
					dfile = line[1] # Disdrometer filename
					starttime = line[2] # N.int(line[2])
					stoptime = line[3] # N.int(line[3])
					centertime = line[4] # N.int(line[4])
					lat = N.float(line[5])
					lon = N.float(line[6])
					alt = N.float(line[7])
					# Append to appropriate lists
					dis_list.append(dfile)
					dis_name_list.append(dname)
					starttimes.append(starttime)
					stoptimes.append(stoptime)
					centertimes.append(centertime)
					dlocs.append((lat,lon,alt))
		
				inputfile.readline()
				inputfile.readline()
				inputfile.readline()
				line = inputfile.readline().strip().split(',')
				platform = line[0]
				wavelength = N.float(line[1])
				inputfile.readline()
				radar_dir = inputfile.readline().strip()
				inputfile.readline()
				scattdir = inputfile.readline().strip()
	
				 # Read in start and end times for the radar data analysis
				inputfile.readline()
				line = inputfile.readline().strip().split(',')
				line_int = map(int,line)
				starttimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
				line = inputfile.readline().strip().split(',')
				line_int = map(int,line)
				stoptimerad = datetime(line_int[0],line_int[1],line_int[2],line_int[3],line_int[4],line_int[5])
	
				# Read in plot window bounds
				inputfile.readline()
				line = inputfile.readline().strip().split(',')
				line_float = map(float,line)
				plotxmin = line_float[0]
				plotxmax = line_float[1]
				plotymin = line_float[2]
				plotymax = line_float[3]
	
				# Read in radar name,lat,lon
				inputfile.readline()
				line = inputfile.readline().strip().split(',')
				radar_name = line[0]
				try:
					rlat = N.float(line[1])
				except:
					rlat = None
				try:
					rlon = N.float(line[2])
				except:
					rlon = None
				try:
					ralt = N.float(line[3])
				except:
					ralt = None
	
				#print line[4]
				#print N.float(line[4])
				try:
					el_req = N.float(line[4])
					print "requested elevation angle",el_req
				except:
					el_req = 0.5    # Default to 0.5 degrees
		
				# Read in min range,max range, min azimuth, and max azimuth for radar plotting (may deprecate this)
				inputfile.readline()
				line = inputfile.readline().strip().split(',')
				line_float = map(float,line)
				minrange = line_float[0]
				maxrange = line_float[1]
				minazim = line_float[2]
				maxazim = line_float[3]
	
				radlims = [minrange,maxrange,minazim,maxazim]
				plotlims = [plotxmin,plotxmax,plotymin,plotymax]
	
			# We need the disdrometer locations. If they aren't supplied in the input control file, find them 
			# from the GPS data

			for index,dis_name,dis_filename,dloc in \
				zip(xrange(0,len(dis_list)),dis_name_list,dis_list,dlocs):
	
				if(N.int(dloc[0]) == -1):
					filepath = dis_dir+dis_filename
					GPS_lats,GPS_lons,GPS_stats,GPS_alts,dloc = dis.readPIPSloc(filepath)
					dlocs[index] = dloc
	
				print "Lat/Lon/alt of "+dis_name+": "+str(dloc)

			# ------
			# Grab radar data for comparison if desired

			if(pc.comp_radar):
				if(not pc.loadradopt):

					radar_filelist,radtimes = radar.getradarfilelist(radar_dir,starttime=starttimerad,
											  stoptime=stoptimerad,platform=platform,el_req=el_req)        
					# Read in the radar data, sweep by sweep.
		
					fields_tlist = []
					range_start_tlist = []
					range_tlist = []
					azimuth_start_tlist = []
					azimuth_tlist = []
					rlat_tlist = []
					rlon_tlist = []
					ralt_tlist = []
					el_tlist = []
					fields_D_tlist = []
					dxy_tlist = []
		
					for index,path,sweeptime in zip(xrange(len(radar_filelist)),radar_filelist,radtimes):
						print "Processing file: "+path
			
						if(platform == 'NEXRAD'):
							outfieldnames,fields,range_start,radar_range,azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad = \
							radar.readCFRadial(True,el_req,rlat,rlon,ralt,path,sweeptime,fieldnames)
						else:
							outfieldnames,fields,range_start,radar_range,azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad = \
							radar.readUMXPnc(path,sweeptime,fieldnames)

		
						fields_tlist.append(fields)
						range_start_tlist.append(range_start)
						range_tlist.append(radar_range)
						azimuth_start_tlist.append(azimuth_start_rad)
						azimuth_tlist.append(azimuth_rad)
						rlat_tlist.append(rlat_rad)
						rlon_tlist.append(rlon_rad)
						ralt_tlist.append(ralt)
						el_tlist.append(el_rad)
			
						dxy_list,fields_D = dis.rad2DD2(fields,range_start,radar_range,
						azimuth_start_rad,azimuth_rad,rlat_rad,rlon_rad,ralt,el_rad,dlocs)
			
						fields_D_tlist.append(fields_D)
						dxy_tlist.append(dxy_list)
		
					fields_tarr = N.array(fields_tlist)
					range_start_tarr = N.array(range_start_tlist)
					range_tarr = N.array(range_tlist)
					azimuth_start_tarr = N.array(azimuth_start_tlist)
					azimuth_tarr = N.array(azimuth_tlist)
					rlat_tarr = N.array(rlat_tlist)
					rlon_tarr = N.array(rlon_tlist)
					ralt_tarr = N.array(ralt_tlist)
					el_tarr = N.array(el_tlist)
					fields_D_tarr = N.array(fields_D_tlist)
					dxy_tarr = N.array(dxy_tlist)
		
					if(pc.saveradopt):
						raddate_file=open(radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
										  stoptimerad.strftime(fmt3).strip()+'.txt','w')
						pickle.dump(radtimes,raddate_file)
						pickle.dump(radar_filelist,raddate_file)
						pickle.dump(outfieldnames,raddate_file)
						raddate_file.close()
			
						radnpz_filename = radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
										  stoptimerad.strftime(fmt3).strip()+str(el_req)+'.npz'           
						savevars={}
						savevars['fields_tarr'] = fields_tarr
						savevars['range_start_tarr'] = range_start_tarr
						savevars['range_tarr'] = range_tarr
						savevars['azimuth_start_tarr'] = azimuth_start_tarr
						savevars['azimuth_tarr'] = azimuth_tarr
						savevars['rlat_tarr'] = rlat_tarr
						savevars['rlon_tarr'] = rlon_tarr
						savevars['ralt_tarr'] = ralt_tarr
						savevars['el_tarr'] = el_tarr
						savevars['fields_D_tarr'] = fields_D_tarr
						savevars['dxy_tarr'] = dxy_tarr
			
						N.savez(radnpz_filename,**savevars)
		
				if(pc.loadradopt):
					raddate_file=open(radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+        \
										  stoptimerad.strftime(fmt3).strip()+'.txt','r')
					radtimes = pickle.load(raddate_file)
					radar_filelist = pickle.load(raddate_file)
					outfieldnames = pickle.load(raddate_file)
					raddate_file.close()
		
					radnpz_filename = radar_name+'_'+starttimerad.strftime(fmt3).strip()+'_'+    \
										  stoptimerad.strftime(fmt3).strip()+str(el_req)+'.npz'
					radnpz_file = N.load(radnpz_filename)
		
					fields_tarr = radnpz_file['fields_tarr']
					range_start_tarr = radnpz_file['range_start_tarr']
					range_tarr = radnpz_file['range_tarr']
					azimuth_start_tarr = radnpz_file['azimuth_start_tarr']
					azimuth_tarr = radnpz_file['azimuth_tarr']
					rlat_tarr = radnpz_file['rlat_tarr']
					rlon_tarr = radnpz_file['rlon_tarr']
					ralt_tarr = radnpz_file['ralt_tarr']
					el_tarr = radnpz_file['el_tarr']
					fields_D_tarr = radnpz_file['fields_D_tarr']
					dxy_tarr = radnpz_file['dxy_tarr']        
		


			# Outer disdrometer (and deployment) loop
			Mu_retr=[]
			Lam_retr=[]
			D0dict={}
			ZDRdict={}
			Wdict={}
			Rdict={}
			Ntdict={}
			MuDict={}
			LamDict={}

			for index,dis_filename,dis_name,starttime,stoptime,centertime,dloc in \
				zip(xrange(0,len(dis_list)),dis_list,dis_name_list,starttimes,stoptimes,centertimes,dlocs):

				if(pc.timetospace):
					centertime = centertimes[index]

				if (not os.path.exists(image_dir)):
					os.makedirs(image_dir)
	
				if(pc.plot_conv_meteo or pc.plot_DSD_meteo):
					# Create the directory for the plots if it doesn't exist
					meteogram_image_dir = image_dir+'meteograms/'
					if (not os.path.exists(meteogram_image_dir)):
						os.makedirs(meteogram_image_dir)

	
				# Read in the disdrometer data file
				filepath = dis_dir+dis_filename
				datetimesUTC,pdatetimesUTC,flaggedtimes,hailflag,intensities,preciptots,reflectivities,pcounts,pcounts2, \
						   sensortemps,concentrations,onedrop_concentrations,countsMatrix,amplitudes,windspds,winddirrels,winddirabss, \
						winddiags,fasttemps,slowtemps,dewpoints,RHs_derived,RHs,pressures,compass_dirs,    \
						GPS_lats,GPS_lons,GPS_stats,GPS_alts,voltages,DSD_interval,DSD_intervalstr,DSD_index = \
						dis.readPIPS(filepath,basicqc=pc.basicQC,rainfallqc=pc.rainfallQC,rainonlyqc=pc.rainonlyQC,
						strongwindqc=pc.strongwindQC,DSD_interval=pc.DSD_interval)

				DSD_interval_td = timedelta(seconds=DSD_interval)
				DSD_halfinterval_td = timedelta(seconds=DSD_interval/2.)

				# Convert to numpy arrays
				windspds = N.array(windspds)
				winddirrels = N.array(winddirrels)
				winddirabss = N.array(winddirabss)
				winddiags = N.array(winddiags)
				fasttemps = N.array(fasttemps)
				slowtemps = N.array(slowtemps)
				dewpoints = N.array(dewpoints)
				RHs_derived = N.array(RHs_derived)
				RHs = N.array(RHs)
				pressures = N.array(pressures)
				pcounts = N.array(pcounts)
				pcounts2 = N.array(pcounts2)
				intensities = N.array(intensities)
				voltages = N.array(voltages)
				amplitudes = N.array(amplitudes)

				Nc_bin = concentrations.T
				dropperbin = onedrop_concentrations.T
		
				logNc_bin = N.ma.log10(Nc_bin)
				logNc_bin = N.ma.masked_where(Nc_bin < dropperbin, logNc_bin)

				# Compute potential temperature, water vapor mixing ratio, and density
				pt = thermo.caltheta(pressures*100.,slowtemps+273.15)
				qv = thermo.calqv(RHs/100.,pressures*100.,slowtemps+273.15)
				rho=thermo.calrho(pressures*100.,pt,qv)
		
				# Determine start and end times/indices for analysis
	
				datetimesUTCnums = date2num(datetimesUTC)
				pdatetimesUTCnums = date2num(pdatetimesUTC)
	
				if(N.int(starttime) == -1):
					startindex = 0
					pstartindex = 0
					starttime = datetimesUTCnums[startindex]
					pstarttime = pdatetimesUTCnums[startindex]
					plotstarttime = starttime
				else:
					starttime = date2num(datetime(N.int(starttime[:4]),N.int(starttime[4:6]),
									N.int(starttime[6:8]),N.int(starttime[8:10]),N.int(starttime[10:12])))
					pstarttime = starttime
					plotstarttime = starttime
					try:
						startindex = next(i for i,t in enumerate(datetimesUTCnums) if t >= starttime)
						pstartindex = next(i for i,t in enumerate(pdatetimesUTCnums) if t >= starttime)
					except:
						startindex = 0
						pstartindex = 0
						starttime = datetimesUTCnums[startindex]
						pstarttime = pdatetimesUTCnums[startindex]
	
				if(N.int(stoptime) == -1):
					stopindex = N.size(datetimesUTCnums)-1
					pstopindex = N.size(pdatetimesUTCnums)-1
					stoptime = datetimesUTCnums[stopindex]
					pstoptime = pdatetimesUTCnums[pstopindex]
					plotstoptime = stoptime
				else:
					stoptime = date2num(datetime(N.int(stoptime[:4]),N.int(stoptime[4:6]),
									N.int(stoptime[6:8]),N.int(stoptime[8:10]),N.int(stoptime[10:12])))
					pstoptime = stoptime
					plotstoptime = stoptime
					try:
						stopindex = next(i for i,t in enumerate(datetimesUTCnums) if t >= stoptime)
						pstopindex = next(i for i,t in enumerate(pdatetimesUTCnums) if t >= stoptime)
					except:
						stopindex = N.size(datetimesUTCnums)-1
						pstopindex = N.size(pdatetimesUTCnums)-1
						stoptime = datetimesUTCnums[stopindex]
						pstoptime = pdatetimesUTCnums[pstopindex]
	
				disdates = pdatetimesUTC
				disdates_start = [x-DSD_interval_td for x in pdatetimesUTC]
				disdates_start.append(disdates_start[-1]+DSD_interval_td)    # Add an extra 10 sec for the last time bin boundary
				disdates_avg = [x-DSD_halfinterval_td for x in pdatetimesUTC]
	
				if(pc.timetospace and pc.comp_radar):
					delta_t_rad = [(x-centertime).total_seconds() for x in radtimes]
					delta_t_dis = [(x-centertime).total_seconds() for x in disdates]
					delta_t_dis_start = [(x-centertime).total_seconds() for x in disdates_start]
					delta_t_dis_avg = [(x-centertime).total_seconds() for x in disdates_avg]
		
					xloc_rad = [(x*pc.stormmotion)/1000.0 for x in delta_t_rad]
					xloc_dis = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis]
					xloc_dis_start = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis_start]
					xloc_dis_avg = [(x*pc.stormmotion)/1000.0 for x in delta_t_dis_avg]

					plotx_rad = N.array(xloc_rad)
					plotx_dis = N.array(xloc_dis)
					plotx_dis_start = N.array(xloc_dis_start)
					plotx_dis_avg = N.array(xloc_dis_avg)

				else:
					if(pc.comp_radar):
						plotx_rad = date2num(radtimes)
					plotx_dis = date2num(disdates)
					plotx_dis_start = date2num(disdates_start)
					plotx_dis_avg = date2num(disdates_avg)

	
				# Resample the 1-s data corresponding to each integrated DSD (for the given DSD interval)
				sec_offset = pdatetimesUTC[0].second

				dummy,windspds_tDSD,winddirs_tDSD,dummy,dummy,dummy,dummy,dummy,dummy,dummy = dis.resamplewind(datetimesUTC,sec_offset,winddirabss,
								windspds,DSD_intervalstr,gusts=False,gustintvstr='3S',center=False)
	
				windspds_tDSD = windspds_tDSD.loc[DSD_index].values
				winddirs_tDSD = winddirs_tDSD.loc[DSD_index].values
	
				fasttemps_tDSD = pd.Series(data=fasttemps,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				slowtemps_tDSD = pd.Series(data=slowtemps,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				dewpoints_tDSD = pd.Series(data=dewpoints,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				RHs_derived_tDSD = pd.Series(data=RHs_derived,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				RHs_tDSD = pd.Series(data=RHs,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				pressures_tDSD = pd.Series(data=pressures,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
				rho_tDSD = pd.Series(data=rho,index=datetimesUTC).resample(DSD_intervalstr,label='right',closed='right',base=sec_offset).mean().loc[DSD_index].values
	
	
				if(pc.calc_DSD):
					# Now for the fun part.  Calculate exponential and gamma distribution size distribution parameters
					# using the method of moments, after Zhang et al. 2008 and Tokay and Short 1996
		
	
					synthbins,exp_DSD,gam_DSD,dis_DSD = dis.calc_DSD(min_diameter,avg_diameter,
									max_diameter,bin_width,Nc_bin,logNc_bin,rho_tDSD,pc.qrQC,pc.qr_thresh,pcounts2)
	
					# Unpack needed values from returned tuples
	
					N_expDSD,N0_exp,lamda_exp,mu_exp,qr_exp,Ntr_exp,refl_DSD_exp,D_med_exp,D_m_exp = exp_DSD
					N_gamDSD,N0_gam,lamda_gam,mu_gam,qr_gam,Ntr_gam,refl_DSD_gam,D_med_gam,D_m_gam,LWC_gam,rainrate = gam_DSD
					Nc_bin,logNc_bin,D_med_disd,D_m_disd,D_mv_disd,D_ref_disd,QR_disd,refl_disd,LWC_disd,M0 = dis_DSD
					N_gamDSD = N_gamDSD.T
					logN_gamDSD = N.ma.log10(N_gamDSD/1000.) # Get to log(m^-3 mm^-1)
					logN_gamDSD = N.ma.masked_where(N_gamDSD < dropperbin, logN_gamDSD)


					if(pc.calc_dualpol):
						# Calculate polarimetric variables using the T-matrix technique
						scattfile = scattdir+'SCTT_RAIN_fw100.dat'

			#            dualpol_dis = dis.calpolrain(wavelength,scattfile,Nc_bin,bin_width) # use for raw distribution
						dualpol_dis = dis.calpolrain(wavelength,scattfile,(N_gamDSD/1000.),bin_width) # use for gamma fit

						# Unpack needed values from returned tuple
						Zh,Zv,Zhv,dBZ,ZDR,KDP,RHV,intv,d,fa2,fb2 = dualpol_dis
		
				if(pc.plot_DSD_meteo):
					# Prepare arrays for plotting
					DSDstarttimes = date2num(disdates_start[pstartindex:pstopindex+1])
					DSDmidtimes = date2num(disdates_avg[pstartindex:pstopindex+1])

			#        logNc_bin_plot = logNc_bin[:,pstartindex:pstopindex+1] # use for raw distribution
					logNc_bin_plot = logN_gamDSD[:,pstartindex:pstopindex+1] # use for gamma fit

					disvars = {'min_diameter':min_diameter,'DSDstarttimes':DSDstarttimes,
							   'DSDmidtimes':DSDmidtimes,'logNc_bin':logNc_bin_plot}

					# Plot one meteogram for each dualpol variable, otherwise just plot a meteogram for reflectivity
					# based on 6th moment of disdrometer DSD (Rayleigh approximation).
					# Important! Currently the dualpol calculation for the disdrometer data assumes rain only
					# and otherwise ignores data in all bins larger than 9 mm.  For this reason it is recommended
					# to first set the "rain_only_QC" flag to True in disdrometer_module.py

					if(pc.calc_DSD):
						# 3-min average of median volume diameter and radar reflectivity (Rayleigh approx.) from disdrometer

						window = int(180./DSD_interval)
			#            D_0_dis_avg = pd.Series(D_med_disd).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for raw distribution
						D_0_dis_avg = pd.Series(D_med_gam).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for gamma fit
			#            dBZ_ray_dis_avg = pd.Series(refl_disd).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for raw distribution
						dBZ_ray_dis_avg = pd.Series(refl_DSD_gam).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values # use for gamma fit
						disvars['D_0_dis'] = D_0_dis_avg[pstartindex:pstopindex+1]
						disvars['dBZ_ray'] = dBZ_ray_dis_avg[pstartindex:pstopindex+1]


						if(pc.calc_dualpol):
							# Computed centered running averages of dualpol variables
							for varname,var in zip(['dBZ','ZDR','KDP','RHV'],[dBZ,ZDR,KDP,RHV]):
								if(var.size):
									var_avg = pd.Series(var).rolling(window=window,center=True,win_type='triang',min_periods=1).mean().values        
									disvars[varname] = var_avg[pstartindex:pstopindex+1]
		
					# Mark flagged times with vertical lines, if desired
					if(pc.plot_diagnostics):
						flaggedtimes_index = N.where(flaggedtimes)[0] # Extract indices for flagged times
						flaggedtimes_plot = DSDmidtimes[flaggedtimes_index] # these are the times with wind contamination
						disvars['flaggedtimes']=flaggedtimes_plot[pstartindex:pstopindex+1]

						hailflag_index = N.where(hailflag)[0]
						hailflag_plot = DSDmidtimes[hailflag_index] # these are the times with hail detected
						disvars['hailflag']=hailflag_plot[pstartindex:pstopindex+1]

					# Get radar variables together for comparison, if desired
					radvars = {}
					if(pc.comp_radar):
						# At least add the reflectivity
						indexrad = outfieldnames.index('dBZ')
						dBZ_D_plt = fields_D_tarr[:,index,indexrad]
						radvars = {'radmidtimes':plotx_rad,'dBZ':dBZ_D_plt}
						# Add other polarimetric fields
						if(pc.calc_dualpol):
							for radvarname in ['ZDR','KDP','RHV']:
								if radvarname in outfieldnames:
									indexrad = outfieldnames.index(radvarname)
									dualpol_rad_var = fields_D_tarr[:,index,indexrad]
									if(dualpol_rad_var.size):
										radvars[radvarname] = dualpol_rad_var


					# Prepare axis parameters
					timelimits = [plotstarttime,plotstoptime]
					diamlimits = [0.0,9.0]
		
					axparams = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,
								'minorxlocator':pc.minorlocator,'axeslimits':[timelimits,diamlimits],
								'majorylocator':MultipleLocator(base=1.0),'axeslabels':[None,'D (mm)']}
		
					# Ok, now we should have everything ready to go to plot the meteograms.  Let'er rip!
	
					pm.plotDSDmeteograms(dis_name,meteogram_image_dir,axparams,disvars,radvars.copy())
		
		
				lamda_gam = lamda_gam/1000.
				
			###     Added for shape-slope relation plots       
				mu.extend(mu_gam)
				lamda.extend(lamda_gam)

# 			###     Interpolate radar values to match disdrometer times before starting retrieval    
				dBZ_rad = pd.Series(N.array(radvars['dBZ']), index= radvars['radmidtimes'])
				ZDR_rad = pd.Series(N.array(radvars['ZDR']), index= radvars['radmidtimes'])
				idx = dBZ_rad.index.union(DSDmidtimes)
				dBZ_rad = N.array(dBZ_rad.reindex(idx).interpolate(method='index').reindex(DSDmidtimes))
				idx2 = ZDR_rad.index.union(DSDmidtimes)
				ZDR_rad = N.array(ZDR_rad.reindex(idx2).interpolate(method='index').reindex(DSDmidtimes))
	
				R_rad_retr,D0_rad_retr,mu_rad_retr,lam_rad_retr,N0_rad_retr,Nt_rad_retr,W_rad_retr,sigm_rad_retr,Dm_rad_retr,N_rad_retr = DR.retrieve_DSD(dBZ_rad,ZDR_rad,d,fa2,fb2,intv)
				R_dis_retr,D0_dis_retr,mu_dis_retr,lam_dis_retr,N0_dis_retr,Nt_dis_retr,W_dis_retr,sigm_dis_retr,Dm_dis_retr,N_dis_retr = DR.retrieve_DSD(dBZ,ZDR,d,fa2,fb2,intv)
	
				Mu_retr.extend(mu_dis_retr)
				Lam_retr.extend(lam_dis_retr)

				name = 'mu_retr'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'mu']}
				em.dis_retr_timeseries(mu_gam,mu_rad_retr,mu_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)
	
				name = 'lam_retr'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'lambda']}
				em.dis_retr_timeseries(lamda_gam,lam_rad_retr,lam_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)

				name = 'dBZ'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'dBZ']}
				em.dis_retr_timeseries(dBZ,dBZ_rad,dBZ,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)

				name = 'ZDR'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],None],'axeslabels':[pc.timelabel,r'ZDR']}
				em.dis_retr_timeseries(ZDR,ZDR_rad,ZDR,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)

	
				N_retr=N.array(N_dis_retr)
				N_retr=N_retr.T
	
				N_retr2=N.array(N_rad_retr)
				N_retr2=N_retr2.T

		
				if(pc.plot_DSDs):
					if (not os.path.exists(image_dir+'DSDs/'+dis_name)):
						os.makedirs(image_dir+'DSDs/'+dis_name)
					for t in range(N.size(Nc_bin,axis=1)):
						sum = N.ma.sum(Nc_bin[:,t])
						if(sum != 0.0):
							fig1=plt.figure(figsize=(8,6))
							ax1=fig1.add_subplot(111)
							plt.title('Disdrometer 10-sec DSD Fits for Time '+disdates[t].strftime(fmt2)+' EST')
							ax1.bar(min_diameter,Nc_bin[:,t]*1000.0,max_diameter-min_diameter,log=True,color='tan')
							ax1.plot(synthbins,N_expDSD[t],lw=2,label='Exp')
							ax1.plot(avg_diameter,N_gamDSD[:,t],lw=2,label='Gamma')
							ax1.plot(d,N_retr[:,t]*1000.,lw=2,label='Dis Retr')
							ax1.plot(d,N_retr2[:,t]*1000.,lw=2,label='Rad Retr')
							ax1.set_yscale('log')
							ax1.set_ylim(10.**2.0,10.**8.5)
							ax1.set_xlim(0.0,9.0)
							ax1.set_xlabel('Drop Diameter (mm)')
							ax1.set_ylabel(r'N(D) $(m^{-4})$')
							ax1.text(0.50,0.95,'Shape parameter (gamma) = %2.2f'%mu_gam[t],transform=ax1.transAxes)
							ax1.text(0.50,0.90,'Mu = %2.2f'%mu_dis_retr[t],transform=ax1.transAxes)
							ax1.text(0.50,0.85,'Slope parameter (gamma) = %2.2f'%lamda_gam[t],transform=ax1.transAxes)
							ax1.text(0.50,0.80,'Lambda = %2.2f'%lam_dis_retr[t],transform=ax1.transAxes)
							ax1.text(0.50,0.75,'Median Volume Diameter (gamma) =%2.2f'%D_med_gam[t],transform=ax1.transAxes)
							ax1.text(0.50,0.70,'D0 = %2.2f'%D0_dis_retr[t],transform=ax1.transAxes)
							ax1.text(0.50,0.65,'Reflectivity =%2.2f'%refl_disd[t]+'dBZ',transform=ax1.transAxes)
							ax1.text(0.50,0.6,'ZDR = %2.2f'%ZDR[t]+'dB',transform=ax1.transAxes)
							ax1.text(0.50,0.55,'Intensity =%2.2f'%intensities[t]+'mm/hr',transform=ax1.transAxes)
							ax1.text(0.50,0.5,'Particle count (QC) = '+str(pcounts2[t]),transform=ax1.transAxes)
							plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
							plt.savefig(image_dir+'DSDs/'+dis_name+'/'+dis_name+'_t'+str(t)+'DSD_plot.png',dpi=200,bbox_inches='tight')
							plt.close(fig1)


				Zh_Cao = N.arange(20,61,1)
				Zdr_Cao = 10**((-2.6857*10**-4*Zh_Cao**2)+0.04892*Zh_Cao-1.4287)
	
				if (not os.path.exists(image_dir+'scattergrams/')):
					os.makedirs(image_dir+'scattergrams/')
			#     
			###     RELATIONS FROM CAO ET AL. 2008


				Zh_rad = pow(10.,dBZ_rad/10)
# 	
# 				###     Radar measured  
				Nt_rad_emp,D0_rad_emp,W_rad_emp,R_rad_emp = em.empirical(Zh_rad,ZDR_rad)  
# 	
# 				###     Disdrometer measured     
				Nt_dis_emp,D0_dis_emp,W_dis_emp,R_dis_emp = em.empirical(Zh,ZDR)
# 	

			###     add for disdrometer timeseries   

				name = 'D0_dis'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,5.0]],'axeslabels':[pc.timelabel,r'D0']}
				em.retr_timeseries(D_med_disd,D_med_gam,D0_rad_retr,D0_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)
	
				name = 'Nt_dis'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,20000.0]],'axeslabels':[pc.timelabel,r'Nt']}
				em.retr_timeseries(M0,Ntr_gam,Nt_rad_retr,Nt_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)

				name = 'W_dis'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,10.0]],'axeslabels':[pc.timelabel,r'LWC']}
				em.retr_timeseries(LWC_disd,LWC_gam,W_rad_retr,W_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)
	
				name = 'R_dis'
				axparamdict1 = {'majorxlocator':pc.locator,'majorxformatter':pc.formatter,'minorxlocator':pc.minorlocator,'axeslimits':[[plotstarttime,plotstoptime],[0.0,250.0]],'axeslabels':[pc.timelabel,r'Rain Rate']}
				em.retr_timeseries(intensities,rainrate,R_rad_retr,R_dis_retr,pstartindex,pstopindex,DSDmidtimes,axparamdict1,image_dir,dis_name,name)

# 			###     Figure 9a-c from Cao et al. 
# 	
				name = 'R'
				maxlim = 10**-2.
				minlim = 10**-4.
				yscale = 'log'
				label = 'R/Zh'
				em.one2one(intensities/Zh,rainrate/Zh,R_dis_retr/Zh,R_rad_retr/Zh_rad,maxlim,minlim,image_dir,dis_name,name,yscale,label)

				name = 'D0'
				maxlim = 0.0
				minlim = 6.0
				yscale = 'linear'
				label = 'D0'
				em.one2one(D_med_disd,D_med_gam,D0_dis_retr,D0_rad_retr,maxlim,minlim,image_dir,dis_name,name,yscale,label)

				name = 'Nt'
				maxlim = 10.0
				minlim = 10**-3.
				yscale = 'log'
				label = 'Nt/Zh'
				em.one2one(M0/Zh,Ntr_gam/Zh,Nt_dis_retr/Zh,Nt_rad_retr/Zh_rad,maxlim,minlim,image_dir,dis_name,name,yscale,label)

				name = 'W'
				maxlim = 10**-3.
				minlim = 10**-6.
				yscale = 'log'
				label = 'W/Zh'
				em.one2one(LWC_disd/Zh,LWC_gam/Zh,W_dis_retr/Zh,W_rad_retr/Zh_rad,maxlim,minlim,image_dir,dis_name,name,yscale,label)

				###     Figure 8a-c from Cao et al.

				name = 'D0'
				ymin = 0.0
				ymax = 5.0
				ylabel = 'D0'
				em.scatters(D_med_disd,D_med_gam,D0_dis_retr,D0_rad_retr,ZDR_rad,ZDR,DSDmidtimes,ymin,ymax,image_dir,dis_name,name,ylabel)

				name = 'W'
				ymin = -6.0
				ymax = -1.0
				ylabel = 'log(W/Zh)'
				em.scatters(N.log10(LWC_disd/Zh),N.log10(LWC_gam/Zh),N.log10(W_dis_retr/Zh),N.log10(W_rad_retr/Zh_rad),ZDR_rad,ZDR,DSDmidtimes,ymin,ymax,image_dir,dis_name,name,ylabel)

				name = 'R'
				ymin = -5.0
				ymax = 0.0
				ylabel = 'log(R/Zh)'
				em.scatters(N.log10(intensities/Zh),N.log10(rainrate/Zh),N.log10(R_dis_retr/Zh),N.log10(R_rad_retr/Zh_rad),ZDR_rad,ZDR,DSDmidtimes,ymin,ymax,image_dir,dis_name,name,ylabel)

				name = 'Nt'
				ymin = -4.0
				ymax = 2.0
				ylabel = 'log(Nt/Zh)'
				em.scatters(N.log10(M0/Zh),N.log10(Ntr_gam/Zh),N.log10(Nt_dis_retr/Zh),N.log10(Nt_rad_retr/Zh_rad),ZDR_rad,ZDR,DSDmidtimes,ymin,ymax,image_dir,dis_name,name,ylabel)
				


			#     D0dict[dis_name+'_obs'] = D_med_disd
			#     D0dict[dis_name+'_mom'] = D_med_gam
			#     D0dict[dis_name+'_dis_retr'] = D0_dis_retr
			#     D0dict[dis_name+'_rad_retr'] = D0_rad_retr
			#     ZDRdict[dis_name+'_dis'] = ZDR
			#     ZDRdict[dis_name+'_rad'] = ZDR_rad
			#     Zhdict[dis_name+'_dis'] = Zh
			#     Zhdict[dis_name+'_rad'] = Zh_rad
			#     dBZdict[dis_name+'_] = dBZ
			#     Wdict[dis_name+'_obs'] = LWC_disd
			#     Rdict[dis_name+'_obs'] = intensities
			#     Ntdict[dis_name+'_obs'] = M0
	
			# name = 'D0'
			# ymin = 0.0
			# ymax = 5.0
			# ylabel = 'D0'
			# em.PIPS(D0dict['PIPS_1A_obs'],D0dict['PIPS_1B_obs'],D0dict['PIPS_2A_obs'],D0dict['PIPS_2B_obs'],ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
			# 
			# name = 'W'
			# ymin = -6.0
			# ymax = -1.0
			# ylabel = 'log(W/Zh)'
			# em.PIPS(N.log10(Wdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Wdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Wdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Wdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
			# 
			# name = 'R'
			# ymin = -5.0
			# ymax = 0.0
			# ylabel = 'log(R/Zh)'
			# em.PIPS(N.log10(Rdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Rdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Rdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Rdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)
			# 
			# name = 'Nt'
			# ymin = -4.0
			# ymax = 2.0
			# ylabel = 'log(Nt/Zh)'
			# em.PIPS(N.log10(Ntdict['PIPS_1A_obs']/ZDRdict['PIPS_1A_rad']),N.log10(Ntdict['PIPS_1B_obs']/ZDRdict['PIPS_1B_rad']),N.log10(Ntdict['PIPS_2A_obs']/ZDRdict['PIPS_2A_rad']),N.log10(Ntdict['PIPS_2B_obs']/ZDRdict['PIPS_2B_rad']),ZDRdict['PIPS_1A_obs'],ZDRdict['PIPS_1B_obs'],ZDRdict['PIPS_2A_obs'],ZDRdict['PIPS_2B_obs'],ymin,ymax,image_dir,dis_name,name,ylabel)


			###     Figure 2 from Brandes et al. 2004
# 
# 			one_x = N.linspace(0.0,60.0)
# 			upper_y = N.exp(1.01*10**-4*one_x**3 - 7.09*10**-3*one_x**2 + 2.38*10**-1*one_x - 3.44)
# 			lower_y = N.exp(2.12*10**-4*one_x**2 + 6.48*10**-2*one_x - 3.87)
# 
# 			fig1=plt.figure(figsize=(8,8))
# 			ax1=fig1.add_subplot(111)
# 			ax1.scatter(ZDRdict['PIPS_1B'], ZDRdict['PIPS_1B_obs'], marker='.', label='PIPS 1B')
# 			ax1.scatter(ZDRdict['PIPS_2A'], ZDRdict['PIPS_2A_obs'], marker='.', label='PIPS 2A')
# 			ax1.scatter(ZDRdict['PIPS_2B'], ZDRdict['PIPS_2B_obs'], marker='.', label='PIPS 2B')
# 			ax1.set_xlim(0.0,60.0)
# 			ax1.set_ylim(-1.0,4.0)
# 			ax1.set_xlabel('Reflectivity (dBZ)')
# 			ax1.set_ylabel('ZDR (dB)')
# 			ax1.plot(one_x,upper_y,color='k')
# 			ax1.plot(one_x,lower_y,color='k')
# 			plt.legend(loc='upper left',numpoints=1,ncol=1,fontsize=8)
# 			plt.savefig(image_dir+'scattergrams/brandes.png',dpi=200,bbox_inches='tight')
# 			plt.close(fig1)


# Plot the lambda-mu relation and fit with 2nd order polynomial 
print len(lamda)
# Lam = []
# Mu = []
# for n1 in xrange(0,len(lamda)):
# 	lam = lamda[n1]
# 	if(lam > 0.0):
# 		Lam.append(lamda[n1])
# 		Mu.append(mu[n1])
lamda = N.array(lamda)
mu = N.array(mu)
lamda = lamda[~N.isnan(lamda)]
mu = mu[~N.isnan(mu)]	
N.savetxt('LAMDA.txt',lamda)
N.savetxt('MU.txt',mu)
poly=N.polynomial.polynomial.polyfit(lamda,mu,2)
polynomial=N.polynomial.polynomial.Polynomial(poly)

xx = N.linspace(0.0, 30.0)
yy = polynomial(xx)
y2 = -0.0201*xx**2. + 0.902*xx - 1.718
y3 = -0.016*xx**2. + 1.213*xx - 1.957
#yy2 = polynomial2(xx)
		
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(111)
plt.title('Shape-Slope Relation')
ax1.scatter(lamda,mu, color='k', marker='.')
#ax1.scatter(Lam_retr,Mu_retr,color='g')
ax1.plot(xx,yy,label='Our Relation')
ax1.plot(xx,y2,label='Cao Relation')
ax1.plot(xx,y3,label='Zhang Relation')
#ax1.plot(xx,yy2,color='b')
#ax1.set_xlim(0.0,20.0)
#ax1.set_ylim(-5.0,20.0)
ax1.set_xlabel('Slope parameter')
ax1.set_ylabel('Shape parameter')
ax1.text(0.05,0.93,'# of Points: %2.1f'%len(lamda), transform=ax1.transAxes, fontsize=8.)
ax1.text(0.05,0.90,'%2.4f'%poly[2]+'*lam^2 + %2.4f'%poly[1]+'*lam + %2.4f'%poly[0], transform=ax1.transAxes, fontsize=8.)
plt.legend(loc='upper left',numpoints=1,ncol=3,fontsize=8)
plt.savefig('shape_slope.png',dpi=200,bbox_inches='tight')
plt.close(fig)

print(poly)
print len(lamda)

Lam = []
Mu = []

for n1 in xrange(0,len(lamda)):
	lam = lamda[n1]
	if(lam < 15.):
		Lam.append(lamda[n1])
		Mu.append(mu[n1])

poly2=N.polynomial.polynomial.polyfit(Lam,Mu,2)
polynomial2=N.polynomial.polynomial.Polynomial(poly2)

xx = N.linspace(0.0, 30.0)
yy = polynomial2(xx)
y2 = -0.0201*xx**2. + 0.902*xx - 1.718
y3 = -0.016*xx**2. + 1.213*xx - 1.957
#yy2 = polynomial2(xx)
		
fig=plt.figure(figsize=(8,8))
ax1=fig.add_subplot(111)
plt.title('Shape-Slope Relation with Slope < 15')
ax1.scatter(Lam,Mu, color='k', marker='.')
#ax1.scatter(Lam_retr,Mu_retr,color='g')
ax1.plot(xx,yy,label='Our Relation')
ax1.plot(xx,y2,label='Cao Relation')
ax1.plot(xx,y3,label='Zhang Relation')
#ax1.plot(xx,yy2,color='b')
#ax1.set_xlim(0.0,20.0)
#ax1.set_ylim(-5.0,20.0)
ax1.set_xlabel('Slope parameter')
ax1.set_ylabel('Shape parameter')
ax1.text(0.05,0.93,'# of Points: %2.1f'%len(Lam), transform=ax1.transAxes,fontsize=8.)
ax1.text(0.05,0.90,'%2.4f'%poly2[2]+'*lam^2 + %2.4f'%poly2[1]+'*lam + %2.4f'%poly2[0], transform=ax1.transAxes, fontsize = 8.)
plt.legend(loc='upper left',numpoints=1,ncol=3,fontsize=8)
plt.savefig('shape_slope_15under.png',dpi=200,bbox_inches='tight')
plt.close(fig)
		
print(poly2)
print len(Lam)







