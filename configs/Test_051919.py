"""Configuration script for 05/19/19 test deployment"""

PIPS_IO_dict = {
    'dataset_name': 'Test_051919',
    'deployment_names': ['Test_051919'] * 4,
    'input_txt_dir': '/Users/dawson29/Dropbox/PIPS_data/2019/SPOTTR2019/051919_test/',
    'PIPS_dir': '/Users/dawson29/Dropbox/PIPS_data/2019/SPOTTR2019/051919_test/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/PIPS_data/2019/SPOTTR2019/051919_test/plots',
    'PIPS_types': ['PIPS'] * 4,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1A_Test_051919_merged.txt', 'PIPS1B_Test_051919_merged.txt',
                       'PIPS2A_Test_051919_merged.txt', 'PIPS2B_Test_051919_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_Test_051919_PIPS1A_10s.nc',
                          'parsivel_combined_Test_051919_PIPS1B_10s.nc',
                          'parsivel_combined_Test_051919_PIPS2A_10s.nc',
                          'parsivel_combined_Test_051919_PIPS2B_10s.nc'],
    'start_times': ['201905192000'] * 4,
    'end_times': ['201905192200'] * 4,
    'requested_interval': 10.
}

PIPS_qc_dict = {
    'strongwindQC': True,
    'splashingQC': True,
    'marginQC': True,
    'rainfallQC': False,
    'rainonlyQC': False,
    'hailonlyQC': False,
    'graupelonlyQC': False,
    'basicQC': False,
}

radar_config_dict = {
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': False,
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