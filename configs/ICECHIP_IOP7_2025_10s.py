"""Configuration script for ICECHIP 2025 IOP7 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP7_052925',
    'deployment_names': ['IOP7_052925'] * 3,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP7_052925/csv/',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP7_052925/netcdf/',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP7_052925/plots/10s/',
    'PIPS_types': ['PIPS'] * 3,
    'probe_set': 'ICECHIP_2025_B',
    'PIPS_names': ['PIPS1A', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP7_052925_merged.txt', 'PIPS3A_IOP7_052925_merged.txt',
                       'PIPS3B_IOP7_052925_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP7_052925_PIPS1A_10s.nc',
                          'parsivel_combined_IOP7_052925_PIPS3A_10s.nc',
                          'parsivel_combined_IOP7_052925_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP7_052925_PIPS1A.nc',
                          'conventional_raw_IOP7_052925_PIPS3A.nc',
                          'conventional_raw_IOP7_052925_PIPS3B.nc'],
    'start_times': ['20250529213000'] * 3,
    'end_times': ['20250529233500'] * 3,
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
