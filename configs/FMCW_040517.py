"""Configuration script for VORTEX-SE FMCW deployment"""

PIPS_IO_dict = {
    'dataset_name': 'FMCW_2017_040517',
    'deployment_names': ['FMCW_2017_040517'],
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset/plots',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS2A'],
    'PIPS_filenames': ['PIPS2A_FMCW_040517.txt'],
    'start_times': ['20170405150000'],
    'end_times': ['20170405235959'],
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
    'comp_radar': False,
    'clean_radar': False,
    'calc_dualpol': False,
    'plot_retrieval': False,
    'radar_name': 'KHTX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2017/NEXRAD/PIPS2A_FMCW/0430/CFRadial/modified/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170405150000',
    'radar_end_timestamp': '20170405235959',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}