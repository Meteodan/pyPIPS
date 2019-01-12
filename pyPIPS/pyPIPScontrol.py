# pyPIPScontrol.py
# This module contains controlling parameters for the pyPIPS.py script
#
# This module is imported into pyPIPS.py by default with the line inside pyPIPS.py
# that reads "import pyPIPScontrol as pc".  However, it is highly recommended that the
# user copy this module and make their own version.  The module would then be dynamically
# imported by passing the user's version as a command line argument to pyPIPS.py

import matplotlib.dates as dates
import matplotlib.ticker as ticker

saveradopt = True      # Option to save reflectivity at DSD locations to file to avoid unneccessary
# recalculation
loadradopt = False       # Option to load the above

# Plot/analysis control flags

calc_DSD = True         # Calculate DSD fits (exp and gamma)?
calc_evap = False       # Calculate evaporation rates?
calc_dualpol = True     # Calculate polarimetric radar variables (currently assumes rain only)?
comp_radar = False       # Compare with reflectivity from 88D?
radar_save_dir = 'radar_files'


plot_opt = True        # Plot or not?
plot_DSD_meteo = True   # Plot meteograms of DSDs?
avgwindow = 180.        # Averaging window in seconds for plotting of DSD parameters on meteograms
plot_only_precip = False  # Restrict plotting window to start and end of precip in requested period?
plot_DSDs = False      # Plot individual DSDs and fits?
DSD_interval = 60.0     # Interval of DSD (Minimum 10 s, multiple of 10 s)
plot_DSDderived = True  # Plot some internal derived quantities from the Parsivel?
plot_conv_meteo = True  # Plot meteograms of conventional fields?
meteo_T_Td_range = [-10., 25.]
meteo_ws_range = [0.0, 15.0]
plot_radar = False       # Plot radar base scan with disdrometer locations?
clean_radar = True      # Remove non-precipitation values from radar timeseries?
plot_scat = False        # Plot scattergrams of Zdr vs Z?
plot_retrieved = True   # Plot meteograms and scattergrams of retrieved variables?
filter_bimodal = False  # Filter out times that have bimodal distributions
# Add diagnostic information to various meteograms (e.g. wind quality flag)
# and create additional diagnostic plots (e.g. battery voltage, laser signal amplitude, etc.)
plot_diagnostics = True
# Plot scattergrams/mu-lambda relation from entire data set and save variables for SATP?
plot_outer = True

# dateformat = '%H:%M' # '%d/%H'
dateformat = '%d/%H'
formatter = dates.DateFormatter(dateformat)
# locator = MinuteLocator(interval=15) # HourLocator(interval=3)
# minorlocator = MinuteLocator(interval=1)
# locator = MinuteLocator(byminute=[0,15,30,45]) # HourLocator(interval=1)
# minorlocator = MinuteLocator(byminute=range(0,60,5)) # MinuteLocator(byminute=[0,15,30,45])
locator = ticker.HourLocator(interval=3)
minorlocator = ticker.HourLocator(interval=1)
# timelabel = 'Time (HH:MM)' # Time (day/hour)'
timelabel = 'Time (day/hour) UTC'

strongwindQC = True      # Remove time records that are contaminated by strong wind?
splashingQC = False      # Remove drops that result from splashing?
marginQC = False         # Remove drops that result from margin falls?
basicQC = False          # Performs strongwind, splashing, and margin QC only
# Perform QC on rain fall speed? (rejects particles falling too fast or too slow to be rain)
rainfallQC = False
# Note that parameters controlling tolerance are in disdrometer_module.py
rainonlyQC = False       # Remove all particles that are probably not rain?
hailonlyQC = False      # Remove all particles that are probably not hail?
graupelonlyQC = False   # Remove all particles that are probably not graupel?

masklowcounts = False   # Mask out DSDs below a certain particle count/rainrate? (currently
# hardcoded in calc_DSD in disdrometer_module.py)

# Following may not be working yet
wind_thresh = 5.0       # wind speed threshold (1-min average) above which to throw out DSDs
qrQC = False             # Perform additional QC on disdrometer DSDs based on unrealistic qr?
qr_thresh = 15.0        # Threshold of qr above which to throw out DSD (g/kg)

# Note, the following are for plotting a V-D diagram with the QC masks, *not* for actually
# performing the QC. They aren't working right now anyway.
plot_splashingQC = False
plot_marginQC = False
plot_strongwindQC = False
plot_rainonlyQC = False
plot_rainfallspeedQC = False

# Not currently used, may be deprecated.
reverse_times = False     # Reverse time axis?
timetospace = False      # Perform a simple time-to-space conversion on the data?
# 10.0, (5 June) 12.55 (9 June), 13.0 (7 June)     # Approximate storm motion in m/s
stormmotion = 12.55
xstart = -2.0
xstop = 13.0
xtickintv = 1.0
enhancesmalldrops = False        # For testing, artificially boost concentration of small drops?
enhancefactor = 2.0
enhancethresh = 0.562           # Drop size at and below which to enhance the small drops
