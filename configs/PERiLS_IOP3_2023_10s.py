"""Configuration script for PERiLS 2023 IOP3 (03/24/23)"""

PIPS_IO_dict = {
    'dataset_name': 'IOP3_032423',
    'deployment_names': ['IOP3_032423'] * 2,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP3_032423/csv',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP3_032423/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP3_032423/plots/10s',
    'PIPS_types': ['PIPS'] * 2,
    'PIPS_names': ['PIPS2A', 'PIPS3A'],
    'PIPS_filenames': ['PIPS2A_IOP3_032423_merged.txt',
                       'PIPS3A_IOP3_032423_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP3_032423_PIPS2A_10s.nc',
                          'parsivel_combined_IOP3_032423_PIPS3A_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP3_032423_PIPS2A.nc',
                          'conventional_raw_IOP3_032423_PIPS3A.nc'],
    'start_times': ['20230324234049', '20230324233519'],
    'end_times': ['20230325025059', '20230325024159'],
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
