"""Configuration script for PERiLS 2022 IOP1 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP1_032222',
    'deployment_names': ['IOP1_032222'] * 4,
    'input_txt_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/PIPS_data/IOP1_032222/csv/',
    'PIPS_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/PIPS_data/IOP1_032222/netcdf',
    'plot_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/PIPS_data/IOP1_032222/plots',
    'PIPS_types': ['PIPS'] * 4,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1A_IOP1_032222_merged.txt', 'PIPS1B_IOP1_032222_merged.txt',
                       'PIPS2A_IOP1_032222_merged.txt', 'PIPS2B_IOP1_032222_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP1_032222_PIPS1A_10s.nc',
                          'parsivel_combined_IOP1_032222_PIPS1B_10s.nc',
                          'parsivel_combined_IOP1_032222_PIPS2A_10s.nc',
                          'parsivel_combined_IOP1_032222_PIPS2B_10s.nc'],
    'start_times': ['20220322192000'] * 4,
    'end_times': ['20220322214500'] * 4,
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
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': False,
    'clean_radar': False,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/NEXRAD/IOP1/KGWX/',
    'field_names': ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220322190000',
    'radar_end_timestamp': '20220322191500',
    'scatt_dir': '/depot/dawson29/apps/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}