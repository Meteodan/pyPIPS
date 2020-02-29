"""Configuration script for VORTEX-SE IOP3 2016"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'IOP3_2016',
    'deployment_names': ['IOP3_D1_2016'] * 4,
    'PIPS_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset',
    'plot_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset/plots',
    'PIPS_types': ['PIPS', 'PIPS', 'PIPS', 'PIPS'],
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS_1A_IOP_3_D1.txt', 'PIPS_1B_IOP_3_D1.txt', 'PIPS_2A_IOP_3_D1.txt',
                       'PIPS_2B_IOP_3_D1.txt'],
    'start_times': [None] * 4, # ['20160331220000', '20160331220000', '20160331220000', '20160331220000'],
    'end_times': [None] * 4, # ['20160331234000', '20160331234000', '20160331234000', '20160331234000'],
    'requested_interval': 60.
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
    'radar_name': 'KGWX',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_3/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160331220000',
    'radar_end_timestamp': '20160331221000', # '20160331234000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}