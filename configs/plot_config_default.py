"""Configuration file for the various plotting scripts"""
import matplotlib.dates as dates
import matplotlib.ticker as ticker

PIPS_plotting_dict = {
    # For DSD meteograms
    'DSD_D_range': (0.0, 9.0),
    'DSD_D_ytick': 1.0,
    'DSD_param_avg_window': 180.,
    'plot_only_precip': False,
    # For velocity-diameter plots
    'velocity_range': (0.0, 24.0),
    'diameter_range': (0.0, 26.0),
    # For meteograms
    'plot_diagnostics': False,
    'T_Td_range': (-10., 25.),
    'ws_range': (0.0, 15.0),
    # Axis locating and formatting parameters
    'majorxformatter': dates.DateFormatter('%H:%M'),
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45]),
    'minorxlocator': dates.MinuteLocator(byminute=range(0, 60, 5)),
    'xlabel': 'Time (HH:MM)'
}