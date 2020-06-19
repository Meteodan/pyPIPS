"""Configuration script for TriPIPS dataset"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'TriPIPS_082019',
    'deployment_names': ['082019'],
    'PIPS_dir': '/Volumes/scr_fast/Projects/TriPIPS/2019',
    'plot_dir': '/Volumes/scr_fast/Projects/TriPIPS/2019/plots',
    'PIPS_types': ['TriPIPS'],
    'PIPS_names': ['TriPIPS'],
    'PIPS_filenames': ['TriPIPS_082019_merged.txt'],
    'start_times': ['20190820173000'],
    'end_times': ['20190820194500'],
    'requested_interval': 60.,
    'geo_locs': [(40.47489, -86.99165, 214.)]}

PIPS_qc_dict = {
    'strongwindQC': True,
    'splashingQC': True,
    'marginQC': True,
    'rainfallQC': True,
    'rainonlyQC': True,
    'hailonlyQC': False,
    'graupelonlyQC': False,
    'basicQC': False,
}

# TODO: Update this
# NOTE: probably don't need this for this config file anyway, since this is just for loading
# the whole dataset for those calculations that work on all the DSDs (like SATP, fitting mu-lambda,
# etc.)
radar_config_dict = {
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': False,
    'clean_radar': False,
    'calc_dualpol': False,
    'radar_name': 'xtrra',
    'radar_type': 'XTRRA',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/XTRRA/2019',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 5.1,
    'radar_start_timestamp': '20160331220000',
    'radar_end_timestamp': '20160331234000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/X-Band/',
    'wavelength': 3.2
}