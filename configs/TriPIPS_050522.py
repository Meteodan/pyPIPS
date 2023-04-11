"""Configuration script for TriPIPS 05/05/22 precipitation event"""

PIPS_IO_dict = {
    'dataset_name': 'IOP3_040522',
    'deployment_names': ['IOP3_040522'] * 6,
    'input_txt_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/',
    'PIPS_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/netcdf_web',
    'plot_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/TriPIPS/plots',
    'PIPS_types': ['TriPIPS'],
    'PIPS_names': ['TriPIPS'],
    'PIPS_filenames': [''],
    'PIPS_filenames_nc': ['parsivel_combined_IOP3_040522_PIPS1A_10s.nc',
                          'parsivel_combined_IOP3_040522_PIPS1B_10s.nc',
                          'parsivel_combined_IOP3_040522_PIPS2A_10s.nc',
                          'parsivel_combined_IOP3_040522_PIPS2B_10s.nc',
                          'parsivel_combined_IOP3_040522_PIPS3A_10s.nc',
                          'parsivel_combined_IOP3_040522_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP3_040522_PIPS1A.nc',
                          'conventional_raw_IOP3_040522_PIPS1B.nc',
                          'conventional_raw_IOP3_040522_PIPS2A.nc',
                          'conventional_raw_IOP3_040522_PIPS2B.nc',
                          'conventional_raw_IOP3_040522_PIPS3A.nc',
                          'conventional_raw_IOP3_040522_PIPS3B.nc'],
    'start_times': ['20220405143500'] * 6,
    'end_times': ['20220405172000'] * 6,
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
    'comp_radar': True,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/depot/dawson29/data/Projects/PERiLS/obsdata/2022/NEXRAD/IOP3/KGWX/',
    'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'field_names': ['REF', 'ZDR', 'RHO'],  # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20220405143500',
    'radar_end_timestamp': '20220405172000',
    'scatt_dir': '/depot/dawson29/apps/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
