"""Configuration script for PERiLS 2022 IOP2 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP2_030323',
    'deployment_names': ['IOP2_030323'] * 6,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP2_030323/csv',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP2_030323/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP2_030323/plots/60s',
    'PIPS_types': ['PIPS'] * 6,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP2_030323_merged.txt', 'PIPS1B_IOP2_030323_merged.txt',
                       'PIPS2A_IOP2_030323_merged.txt', 'PIPS2B_IOP2_030323_merged.txt',
                       'PIPS3A_IOP2_030323_merged.txt', 'PIPS3B_IOP2_030323_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP2_030323_PIPS1A_60s.nc',
                          'parsivel_combined_IOP2_030323_PIPS1B_60s.nc',
                          'parsivel_combined_IOP2_030323_PIPS2A_60s.nc',
                          'parsivel_combined_IOP2_030323_PIPS2B_60s.nc',
                          'parsivel_combined_IOP2_030323_PIPS3A_60s.nc',
                          'parsivel_combined_IOP2_030323_PIPS3B_60s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP2_030323_PIPS1A.nc',
                          'conventional_raw_IOP2_030323_PIPS1B.nc',
                          'conventional_raw_IOP2_030323_PIPS2A.nc',
                          'conventional_raw_IOP2_030323_PIPS2B.nc',
                          'conventional_raw_IOP2_030323_PIPS3A.nc',
                          'conventional_raw_IOP2_030323_PIPS3B.nc'],
    'start_times': ['20230303080000'] * 6,
    'end_times': ['20230303120000'] * 6,
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
    'radar_name': 'KLZK',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2023/NEXRAD/IOP2/LZK/CFRadial/',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['VEL'],  # ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20230303080000',
    'radar_end_timestamp': '20230303120000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
