"""Configuration script for 08/11/22 test"""

PIPS_IO_dict = {
    'dataset_name': 'Test_081122',
    'deployment_names': ['Test_081122'] * 2,
    'input_txt_dir': '/Users/dawson29/Dropbox/PIPS_data/2022/Test_081122/csv',
    'PIPS_dir': '/Users/dawson29/Dropbox/PIPS_data/2022/Test_081122/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/PIPS_data/2022/Test_081122/plots',
    'PIPS_types': ['PIPS'] * 2,
    'PIPS_names': ['PIPS1A', 'PIPS1B'],
    'PIPS_filenames': ['PIPS1A_Test_081122_merged.txt', 'PIPS1B_Test_081122_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_Test_081122_PIPS1A_10s.nc',
                          'parsivel_combined_Test_081122_PIPS1B_10s.nc'],
    'start_times': ['20220811164000'] * 2,
    'end_times': ['20220811164500'] * 2,
    'requested_interval': 10.
}

PIPS_qc_dict = {
    'strongwindQC': True,
    'splashingQC': True,
    'marginQC': True,
    'rainfallQC': False,
    'rainonlyQC': True,
    'hailonlyQC': False,
    'graupelonlyQC': False,
    'basicQC': False,
}

radar_config_dict = {
    'comp_radar': False,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/PERiLS/obsdata/2022/NEXRAD/IOP1/KGWX/',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO'],  # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220322192000',
    'radar_end_timestamp': '20220322214500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}