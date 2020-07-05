"""Configuration script for VORTEX-SE full PIPS dataset"""
import os
from glob import glob

PIPS_IO_dict = {
    'dataset_name': 'all_but_FMCW',
    'input_txt_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_links',
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/plots',
    'requested_interval': 60.
}

full_PIPS_file_path_list = glob(PIPS_IO_dict['input_txt_dir'] + '/PIPS*txt')
PIPS_file_path_list = [PIPS_path for PIPS_path in full_PIPS_file_path_list
                       if 'FMCW' not in PIPS_path]
numfiles = len(PIPS_file_path_list)
PIPS_filenames = [os.path.basename(PIPS_file_path) for PIPS_file_path in PIPS_file_path_list]
PIPS_names = [PIPS_filename[:6] for PIPS_filename in PIPS_filenames]

# Figure out the deployment names
deployment_names = []
for PIPS_filename in PIPS_filenames:
    if 'IOP' in PIPS_filename:
        IOP_name = PIPS_filename[PIPS_filename.index('IOP'):PIPS_filename.index('.txt')]
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
PIPS_filenames_nc = ['parsivel_combined_{}_{}_60s.nc'.format(deployment_name, PIPS_name)
                     for deployment_name, PIPS_name in zip(deployment_names, PIPS_names)]
PIPS_IO_dict['PIPS_types'] = ['PIPS'] * numfiles
PIPS_IO_dict['PIPS_names'] = PIPS_names
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
    'radar_name': 'KGWX',
    'radar_dir': '/Users/ddawson/Dropbox/Projects/VORTEXSE/obs_data/NEXRAD/IOP3_2016/KGWX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160331220000',
    'radar_end_timestamp': '20160331234000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}