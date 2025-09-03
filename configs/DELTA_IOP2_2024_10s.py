"""Configuration script for PERiLS 2023 IOP3 (03/24/23)"""

PIPS_IO_dict = {
    'dataset_name': 'DELTA_IOP2_032424',
    'deployment_names': ['DELTA_IOP2_032424'] * 2,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/DELTA/temp/incoming_after_DELTA_IOP2/csv',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/DELTA/temp/incoming_after_DELTA_IOP2/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/DELTA/temp/incoming_after_DELTA_IOP2/plots/10s',
    'PIPS_types': ['PIPS'] * 2,
    'PIPS_names': ['PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS3A_DELTA_IOP2_032424_merged.txt',
                       'PIPS3B_DELTA_IOP2_032424_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_DELTA_IOP2_032424_PIPS3A_10s.nc',
                          'parsivel_combined_DELTA_IOP2_032424_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_DELTA_IOP2_032424_PIPS3A.nc',
                          'conventional_raw_DELTA_IOP2_032424_PIPS3B.nc'],
    'start_times': ['20240324225000', '20240324225000'],
    'end_times': ['20240324232000', '20240324232000'],
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
