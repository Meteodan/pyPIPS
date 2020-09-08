"""Configuration script for VORTEX-SE FMCW deployment"""

PIPS_IO_dict = {
    'dataset_name': 'FMCW_2017_032717',
    'deployment_names': ['FMCW_2017_032717'],
    'input_txt_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_links',
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15/',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15/plots',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS2A'],
    'PIPS_filenames': ['PIPS2A_FMCW_032717.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_FMCW_2017_032717_PIPS2A_60s.nc'],
    'start_times': ['20170327190000'],
    'end_times': ['20170327235959'],
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
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2017/NEXRAD/PIPS2A_FMCW/0327/HTX/CFRadial/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170327190000',
    'radar_end_timestamp': '20170327235959',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}