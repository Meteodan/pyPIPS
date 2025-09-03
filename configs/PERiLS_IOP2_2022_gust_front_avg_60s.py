"""Configuration script for PERiLS 2022 IOP2 deployment (gust front average version for 60-s data)"""

PIPS_IO_dict = {
    'dataset_name': 'IOP2_033022',
    'deployment_names': ['IOP2_033022'],
    'input_txt_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2022/PIPS_data/IOP2_033022/csv',
    'PIPS_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2022/PIPS_data/IOP2_033022/netcdf_60s',
    'plot_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2022/PIPS_data/IOP2_033022/plots/60s',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS'],
    'PIPS_filenames': [''],
    'PIPS_filenames_nc': ['parsivel_combined_IOP2_033022_PIPS_avg_60s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP2_033022_PIPS_avg.nc'],
    'start_times': ['20220330234000'],
    'end_times': ['20220331014500'],
    'requested_interval': 60.
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
    'comp_radar': True,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2022/NEXRAD/IOP2/KGWX',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['VEL'],  # ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220330234000',
    'radar_end_timestamp': '20220331014500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
