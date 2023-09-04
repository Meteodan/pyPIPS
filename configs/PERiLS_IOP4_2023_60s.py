"""Configuration script for PERiLS 2023 IOP4 (03/31/23)"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4_033123',
    'deployment_names': ['IOP4_033123'] * 6,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP4_033123/csv',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP4_033123/netcdf',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/PERiLS/obsdata/2023/IOP4_033123/plots/60s',
    'PIPS_types': ['PIPS'] * 6,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP4_033123_merged.txt',
                       'PIPS1B_IOP4_033123_merged.txt',
                       'PIPS2A_IOP4_033123_merged.txt',
                       'PIPS2B_IOP4_033123_merged.txt',
                       'PIPS3A_IOP4_033123_merged.txt',
                       'PIPS3B_IOP4_033123_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP4_033123_PIPS1A_60s.nc',
                          'parsivel_combined_IOP4_033123_PIPS1B_60s.nc',
                          'parsivel_combined_IOP4_033123_PIPS2A_60s.nc',
                          'parsivel_combined_IOP4_033123_PIPS2B_60s.nc',
                          'parsivel_combined_IOP4_033123_PIPS3A_60s.nc',
                          'parsivel_combined_IOP4_033123_PIPS3B_60s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP4_033123_PIPS1A.nc',
                          'conventional_raw_IOP4_033123_PIPS1B.nc',
                          'conventional_raw_IOP4_033123_PIPS2A.nc',
                          'conventional_raw_IOP4_033123_PIPS2B.nc',
                          'conventional_raw_IOP4_033123_PIPS3A.nc',
                          'conventional_raw_IOP4_033123_PIPS3B.nc'],
    'start_times': ['20230401070000'] * 6,
    'end_times': ['20230401090000'] * 6,
    'requested_interval': 60.
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
    'comp_radar': True,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/Projects/PERiLS/obsdata/2023/NEXRAD/IOP4/GWX/CFRadial',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['VEL'], # ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20230401070000',
    'radar_end_timestamp': '20230401090000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
