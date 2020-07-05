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
    'velocity_range': (0.0, 15.0),
    'diameter_range': (0.0, 9.0),
    # For meteograms
    'avgwindow': 180.,
    'plot_diagnostics': True,
    'T_Td_range': (-10., 25.),
    'ws_range': (0.0, 15.0),
    # Axis locating and formatting parameters
    'majorxformatter': dates.DateFormatter('%d:%H'),
    'majorxlocator': dates.HourLocator(byhour=range(0, 24, 6)),
    'minorxlocator': dates.HourLocator(byhour=range(0, 24, 1)),
    'xlabel': 'Time (DD:HH)',
    'fontsize': 8,  # TODO: implement this to override default settings in plotmodule.py
    'fontP': 'x-small'
}