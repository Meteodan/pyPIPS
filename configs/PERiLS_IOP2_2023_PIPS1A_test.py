"""Configuration script for PERiLS 2022 IOP2 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP2_030323',
    'deployment_names': ['IOP2_030323'] * 1,
    'input_txt_dir': '/Users/dawson29/Dropbox/PIPS_data/2023/incoming/PIPS1A/',
    'PIPS_dir': '/Users/dawson29/Dropbox/PIPS_data/2023/incoming/PIPS1A/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/PIPS_data/2023/incoming/PIPS1A/plots/10s',
    'PIPS_types': ['PIPS'] * 1,
    'PIPS_names': ['PIPS1A'],
    'PIPS_filenames': ['PIPS1A_IOP2_030323_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP2_030323_PIPS1A_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP2_030323_PIPS1A.nc'],
    'start_times': ['20230303081700'] * 1,
    'end_times': ['20230303115100'] * 1,
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
