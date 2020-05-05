"""Configuration script for VORTEX-SE full PIPS dataset"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'full_dataset',
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/SATP_retr',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/SATP_retr/plots',
    'requested_interval': 60.
}

PIPS_file_path_list = glob(PIPS_IO_dict['PIPS_dir'] + '/parsivel_combined*nc')
numfiles = len(PIPS_file_path_list)
PIPS_filenames = [os.path.basename(PIPS_file_path) for PIPS_file_path in PIPS_file_path_list]
PIPS_names = [PIPS_filename[-13:-7] for PIPS_filename in PIPS_filenames]
print(PIPS_names)

# Figure out the deployment names
deployment_names = []
for PIPS_filename in PIPS_filenames:
    if 'FMCW' in PIPS_filename:
        deployment_name = 'FMCW_2017_{}'.format(PIPS_filename[-20:-14])
    elif 'IOP' in PIPS_filename:
        IOP_name = PIPS_filename[PIPS_filename.index('IOP'):PIPS_filename.index('_D')+3]
        if '2016' in PIPS_filename:
            IOP_year = '2016'
        elif '2017' in PIPS_filename:
            IOP_year = '2017'
        else:
            IOP_year = 'unknown'
        deployment_name = IOP_name + '_' + IOP_year
    else:
        deployment_name = 'unknown'
    deployment_names.append(deployment_name)
PIPS_IO_dict['deployment_names'] = deployment_names

PIPS_IO_dict['PIPS_types'] = ['PIPS'] * numfiles
PIPS_IO_dict['PIPS_names'] = PIPS_names
PIPS_IO_dict['PIPS_filenames'] = PIPS_filenames
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

# TODO: Update this
# NOTE: probably don't need this for this config file anyway, since this is just for loading
# the whole dataset for those calculations that work on all the DSDs (like SATP, fitting mu-lambda,
# etc.)
radar_config_dict = {
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': True,
    'clean_radar': False,
    'calc_dualpol': True,
    'radar_name': 'KGWX',
    'radar_dir': '/Users/ddawson/Dropbox/Projects/VORTEXSE/obs_data/NEXRAD/IOP3_2016/KGWX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160331220000',
    'radar_end_timestamp': '20160331234000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}