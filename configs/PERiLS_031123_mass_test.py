"""Configuration script for PERiLS 2023 03/11/23 mass test"""

PIPS_IO_dict = {
    'dataset_name': '031123_mass_test',
    'deployment_names': ['031123_mass_test'] * 6,
    'input_txt_dir': '/Users/dawson29/PIPS_data/2023/031123_mass_test/csv',
    'PIPS_dir': '/Users/dawson29/PIPS_data/2023/031123_mass_test/netcdf',
    'plot_dir': '/Users/dawson29/PIPS_data/2023/031123_mass_test/plots/10s',
    'PIPS_types': ['PIPS'] * 6,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_031123_mass_test_merged.txt',
                       'PIPS1B_031123_mass_test_merged.txt',
                       'PIPS2A_031123_mass_test_merged.txt',
                       'PIPS2B_031123_mass_test_merged.txt',
                       'PIPS3A_031123_mass_test_merged.txt',
                       'PIPS3B_031123_mass_test_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_031123_mass_test_PIPS1A_10s.nc',
                          'parsivel_combined_031123_mass_test_PIPS1B_10s.nc',
                          'parsivel_combined_031123_mass_test_PIPS2A_10s.nc',
                          'parsivel_combined_031123_mass_test_PIPS2B_10s.nc',
                          'parsivel_combined_031123_mass_test_PIPS3A_10s.nc',
                          'parsivel_combined_031123_mass_test_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_031123_mass_test_PIPS1A.nc',
                          'conventional_raw_031123_mass_test_PIPS1B.nc',
                          'conventional_raw_031123_mass_test_PIPS2A.nc',
                          'conventional_raw_031123_mass_test_PIPS2B.nc',
                          'conventional_raw_031123_mass_test_PIPS3A.nc',
                          'conventional_raw_031123_mass_test_PIPS3B.nc'],
    'start_times': ['20230311234500'] * 6,
    'end_times': ['20230312160000'] * 6,
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
    'radar_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/NEXRAD/IOP2/KGWX',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220330234000',
    'radar_end_timestamp': '20220331014500',
    'scatt_dir': '/depot/dawson29/apps/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
