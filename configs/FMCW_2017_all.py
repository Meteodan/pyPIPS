"""Configuration script for VORTEX-SE FMCW deployment"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'FMCW_2017_all',
    'input_txt_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_links',
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/plots',
    'requested_interval': 60.
}

PIPS_file_path_list = glob(PIPS_IO_dict['input_txt_dir'] + '/PIPS2A_FMCW*txt')
numfiles = len(PIPS_file_path_list)
PIPS_filenames = [os.path.basename(PIPS_file_path) for PIPS_file_path in PIPS_file_path_list]

# Figure out the deployment names
deployment_names = []
for PIPS_filename in PIPS_filenames:
    deployment_name = 'FMCW_2017_{}'.format(PIPS_filename[-10:-4])
    deployment_names.append(deployment_name)
PIPS_IO_dict['deployment_names'] = deployment_names
PIPS_filenames_nc = ['parsivel_combined_{}_PIPS2A_60s.nc'.format(deployment_name)
                     for deployment_name in deployment_names]
print(deployment_names)
print(PIPS_filenames_nc)
PIPS_IO_dict['PIPS_types'] = ['PIPS'] * numfiles
PIPS_IO_dict['PIPS_names'] = ['PIPS2A'] * numfiles
PIPS_IO_dict['PIPS_filenames'] = PIPS_filenames
PIPS_IO_dict['PIPS_filenames_nc'] = PIPS_filenames_nc
PIPS_IO_dict['start_times'] = [None] * numfiles
PIPS_IO_dict['end_times'] = [None] * numfiles

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
    'comp_radar': False,
    'clean_radar': False,
    'calc_dualpol': False,
    'plot_retrieval': False,
    'radar_name': 'KHTX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2017/NEXRAD/PIPS2A_FMCW/0430/HTX/CFRadial/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170430000000',
    'radar_end_timestamp': '20170430235959',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}