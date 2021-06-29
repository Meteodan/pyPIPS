"""Configuration script for Hurricane Laura deployment"""

PIPS_IO_dict = {
    'dataset_name': 'Laura_2020',
    'deployment_names': ['Laura_2020'] * 3,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/Laura2020/merged/',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/Laura2020/netcdf/',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/Laura2020/plots',
    'PIPS_types': ['PIPS'] * 3,
    'PIPS_names': ['PIPS1A', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1A_merged.txt',
                       'PIPS1B_merged.txt', 'PIPS2B_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_Laura_2020_PIPS1A_60s.nc',
                          'parsivel_combined_Laura_2020_PIPS2A_60s.nc',
                          'parsivel_combined_Laura_2020_PIPS2B_60s.nc'],
    'start_times': [None] * 3,
    'end_times': [None] * 3,
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
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2017/NEXRAD/IOP_1A/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO', 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170325173000',
    'radar_end_timestamp': '20170325184500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}