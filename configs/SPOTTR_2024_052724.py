"""Configuration script for SPOTTR 2024 05/27/24 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'SPOTTR_2024_052724',
    'deployment_names': ['SPOTTR_2024_052724'] * 1,
    'input_txt_dir': '/Users/dawson29/Projects/SPOTTR_2024/PIPS_data/052724/csv',
    'PIPS_dir': '/Users/dawson29/Projects/SPOTTR_2024/PIPS_data/052724/netcdf',
    'plot_dir': '/Users/dawson29/Projects/SPOTTR_2024/PIPS_data/052724/plots/10s',
    'PIPS_types': ['PIPS'] * 1,
    'PIPS_names': ['PIPS3A'],
    'PIPS_filenames': ['PIPS3A_SPOTTR_2024_052724.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_SPOTTR_2024_052724_PIPS3A_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_SPOTTR_2024_052724_PIPS3A.nc'],
    'start_times': ['20240527234000'] * 1,
    'end_times': ['20240528004500'] * 1,
    'requested_interval': 10.
}

# TODO: this is set up in the parsivel_qc.py file now. Maybe keep this here to override that?
# Or remove from here and stick with command-line argument approach?
PIPS_qc_dict = {
    'qc': {
        'strongwindQC': True,
        'splashingQC': True,
        'marginQC': True,
        'rainfallQC': False,
        'rainonlyQC': False,
        'hailonlyQC': False,
        'graupelonlyQC': False,
        'basicQC': False,
    },
    'roqc': {
        'strongwindQC': True,
        'splashingQC': True,
        'marginQC': True,
        'rainfallQC': False,
        'rainonlyQC': True,
        'hailonlyQC': False,
        'graupelonlyQC': False,
        'basicQC': False,
    }
}

radar_config_dict = {
    'comp_radar': True,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KFWS',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/Projects/SPOTTR_2024/NEXRAD_data/052724/KFWS/CFRadial/',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO', 'VEL'], # ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20240527234000',
    'radar_end_timestamp': '20240528004500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
