# pyPIPScontrol.py
# This module contains controlling parameters for the pyPIPS.py script
#
# This module is imported into pyPIPS.py by default with the line inside pyPIPS.py
# that reads "import pyPIPScontrol as pc".  However, it is highly recommended that the
# user copy this module and make their own version.  The module would then be dynamically
# imported by passing the user's version as a command line argument to pyPIPS.py

from matplotlib.dates import *
from matplotlib.ticker import *

enhancesmalldrops = False        # For testing, artificially boost concentration of small drops?
enhancefactor = 2.0            
enhancethresh = 0.562           # Drop size at and below which to enhance the small drops

saveradopt = True      # Option to save reflectivity at DSD locations to file to avoid unneccessary
                        # recalculation
loadradopt = False       # Option to load the above

# Plot/analysis control flags

calc_DSD = True         # Calculate DSD fits (exp and gamma)?
calc_evap = False       # Calculate evaporation rates?
calc_dualpol = True       # Calculate polarimetric radar variables (currently assumes rain only)?
comp_radar = True       # Compare with reflectivity from 88D?

plot_opt = True        # Plot or not?
plot_DSD_meteo = True   # Plot meteograms of DSDs?
plot_DSDs = True      # Plot individual DSDs and fits?
DSD_interval = 10.0     # Interval of DSD (Minimum 10 s, multiple of 10 s)
plot_conv_meteo = True     # Plot meteograms of conventional fields?
plot_radar = True       # Plot radar base scan with disdrometer locations?
plot_scat = True        # Plot scattergrams of Zdr vs Z
reverse_times = False     # Reverse time axis?
timetospace = False      # Perform a simple time-to-space conversion on the data?
stormmotion = 12.55 # 10.0, (5 June) 12.55 (9 June), 13.0 (7 June)     # Approximate storm motion in m/s
xstart = -2.0
xstop = 13.0
xtickintv = 1.0
avg_data = True        # average disdrometer data in time and calculate new DSD fits?
                        # For now just sum over each deployment using start and end times.

dateformat = '%H:%M' # '%d/%H'
formatter = DateFormatter(dateformat)
#dateformat = '%d/%H'
#locator = MinuteLocator(interval=15) # HourLocator(interval=3)
#minorlocator = MinuteLocator(interval=1)
locator = MinuteLocator(byminute=[0,15,30,45]) # HourLocator(interval=1) 
minorlocator = MinuteLocator(byminute=range(0,60,5)) # MinuteLocator(byminute=[0,15,30,45])
timelabel = 'Time (HH:MM)' # Time (day/hour)'
#timelabel = 'Time (day/hour) UTC'

windQC = False          # Perform additional QC on disdrometer DSDs based on wind speed?
wind_thresh = 5.0       # wind speed threshold (1-min average) above which to throw out DSDs
qrQC = False             # Perform additional QC on disdrometer DSDs based on unrealistic qr?
qr_thresh = 15.0        # Threshold of qr above which to throw out DSD (g/kg)
plot_splashingQC = False
plot_marginQC = False
plot_strongwindQC = False
plot_rainonlyQC = False
plot_rainfallspeedQC = False