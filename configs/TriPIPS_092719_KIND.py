"""Configuration script for TriPIPS 09/27/2019"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'TriPIPS_092719',
    'deployment_names': ['TriPIPS_092719'],
    'PIPS_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/2019',
    'plot_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/2019/plots',
    'PIPS_types': ['TriPIPS'],
    'PIPS_names': ['TriPIPS'],
    'PIPS_filenames': ['TriPIPS_092719_merged.txt'],
    'start_times': ['20190927191500'],
    'end_times': ['20190927204500'],
    'requested_interval': 10.,
    'geo_locs': [(40.47489, -86.99165, 214.)]
}

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

radar_config_dict = {
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': True,
    'clean_radar': False,
    'calc_dualpol': True,
    'radar_name': 'KIND',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/2019/NEXRAD/09/27/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20190927190000',
    'radar_end_timestamp': '20190927204500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}