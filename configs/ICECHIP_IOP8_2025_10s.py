"""Configuration script for ICECHIP 2025 IOP8 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP8_060225',
    'deployment_names': ['IOP8_060225'] * 3,
    'input_txt_dir': '/Users/dawson29/Projects/ICECHIP/obsdata/PIPS_data/IOP8_060225/csv/',
    'PIPS_dir': '/Users/dawson29/Projects/ICECHIP/obsdata/PIPS_data/IOP8_060225/netcdf/',
    'plot_dir': '/Users/dawson29/Projects/ICECHIP/obsdata/PIPS_data/IOP8_060225/plots/10s/',
    'PIPS_types': ['PIPS'] * 3,
    'probe_set': 'ICECHIP_2025_B',
    'PIPS_names': ['PIPS1A', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP8_060225_merged.txt', 'PIPS3A_IOP8_060225_merged.txt',
                       'PIPS3B_IOP8_060225_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP8_060225_PIPS1A_10s.nc',
                          'parsivel_combined_IOP8_060225_PIPS3A_10s.nc',
                          'parsivel_combined_IOP8_060225_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP8_060225_PIPS1A.nc',
                          'conventional_raw_IOP8_060225_PIPS3A.nc',
                          'conventional_raw_IOP8_060225_PIPS3B.nc'],
    'start_times': ['20250602211500'] * 3,
    'end_times': ['20250602223500'] * 3,
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
    'radar_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2022/NEXRAD/IOP2/KGWX',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220330234000',
    'radar_end_timestamp': '20220331014500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
