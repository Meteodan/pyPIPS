"""Configuration file for the various plotting scripts"""
import matplotlib.dates as dates
import matplotlib.ticker as ticker
import matplotlib.colors as colors

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
    'plot_diagnostics': False,
    'T_Td_range': (10., 25.),
    # For wind plotting
    'avgwind': True,
    'maskbadwind': True,
    'windavgintv': 60,
    'windgustintv': 3,
    'ws_range': (0.0, 25.0),
    # Axis locating and formatting parameters
    'majorxformatter': dates.DateFormatter('%H:%M'),
    'majorxlocator': dates.MinuteLocator(byminute=[0, 15, 30, 45]),
    'minorxlocator': dates.MinuteLocator(byminute=range(0, 60, 5)),
    'xlabel': 'Time (HH:MM)'
}

axparams_Dm = {
    'var_lims': [[0.0, 4.5], [0.0, 4.5]],
    'col_field': 'RR_obs',
    'col_field_lims': [0.1, 200.],
    'norm': colors.LogNorm(vmin=0.1, vmax=200.),
    'alpha': 0.75,
    'markersize': 5,
    'label_x': r'$D_m$ (obs; mm)',
    'label_y': r'$D_m$ (retrieved; mm)',
    'label_cb': r'$RR$ (mm h$^{-1}$)'
}

