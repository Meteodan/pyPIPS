"""Configuration script for VORTEX-SE IOP3B 2017 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP3B_2017',
    'deployment_names': ['IOP3B_D1_2017'],
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15/plots',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS1B'],
    'PIPS_filenames': ['PIPS1B_2016_IOP3B_D1.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP3B_D1_2017_PIPS1B_60s.nc'],
    'start_times': [None],
    'end_times': [None],
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
    'plot_retrieval': True,
    'radar_name': 'KHTX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2017/NEXRAD/IOP_3B/HTX/CFRadial/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170405230000',
    'radar_end_timestamp': '20170406010000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}

model_config_dict = {
}