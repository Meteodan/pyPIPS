"""Configuration script for PERiLS 2022 IOP2 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP2_033022',
    'deployment_names': ['IOP2_033022'] * 4,
    'input_txt_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/PIPS_data/IOP2_033022/csv/',
    'PIPS_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/PIPS_data/IOP2_033022/netcdf',
    'plot_dir': '/home/fvendl/Projects/PERiLS/IOP2_033022/PIPS_data/plots',
    'PIPS_types': ['PIPS'] * 4,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP2_033022_merged.txt', 'PIPS1B_IOP2_033022_merged.txt',
                       'PIPS2A_IOP2_033022_merged.txt', 'PIPS3B_IOP2_033022_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP2_033022_PIPS1A_60s.nc',
                          'parsivel_combined_IOP2_033022_PIPS1B_60s.nc',
                          'parsivel_combined_IOP2_033022_PIPS2A_60s.nc',
                          'parsivel_combined_IOP2_033022_PIPS3B_60s.nc'],
    'start_times': ['20220330234000'] * 4,
    'end_times': ['20220331014500'] * 4,
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
    'radar_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/NEXRAD/IOP2/KGWX',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220330234000',
    'radar_end_timestamp': '20220331014500',
    'scatt_dir': '/depot/dawson29/apps/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
