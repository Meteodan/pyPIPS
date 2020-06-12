"""Configuration script for VORTEX-SE full PIPS dataset"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'TriPIPS_full_dataset',
    'PIPS_dir': '/Volumes/scr_fast/Projects/TriPIPS/2019',
    'plot_dir': '/Volumes/scr_fast/Projects/TriPIPS/2019/plots',
    'requested_interval': 60.
}

PIPS_file_path_list = glob(PIPS_IO_dict['PIPS_dir'] + '/TriPIPS*txt')
numfiles = len(PIPS_file_path_list)
PIPS_filenames = [os.path.basename(PIPS_file_path) for PIPS_file_path in PIPS_file_path_list]
PIPS_names = ['TriPIPS'] * numfiles
deployment_names = [PIPS_filename[8:14] for PIPS_filename in PIPS_filenames]
PIPS_IO_dict['deployment_names'] = deployment_names

PIPS_IO_dict['PIPS_types'] = ['TriPIPS'] * numfiles
PIPS_IO_dict['PIPS_names'] = PIPS_names
PIPS_IO_dict['PIPS_filenames'] = PIPS_filenames
PIPS_IO_dict['start_times'] = [None] * numfiles
PIPS_IO_dict['end_times'] = [None] * numfiles
PIPS_IO_dict['geo_locs'] = [(40.47489, -86.99165, 214.)] * numfiles

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
    'calc_dualpol': True,
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