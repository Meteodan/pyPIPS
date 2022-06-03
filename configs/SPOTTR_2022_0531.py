"""Configuration script for SPOTTR 2022 05/31 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'SPOTTR_053122',
    'deployment_names': ['SPOTTR_053122'] * 2,
    'input_txt_dir': '/Users/dandawson/Dropbox/Teaching/2022/EAPS_591_SSFW/PIPS_data/053122/csv/',
    'PIPS_dir': '/Users/dandawson/Dropbox/Teaching/2022/EAPS_591_SSFW/PIPS_data/053122/netcdf',
    'plot_dir': '/Users/dandawson/Dropbox/Teaching/2022/EAPS_591_SSFW/PIPS_data/053122/plots',
    'PIPS_types': ['PIPS'] * 2,
    'PIPS_names': ['PIPS1B', 'PIPS2A'],
    'PIPS_filenames': ['PIPS1B_SPOTTR_053122_merged.txt', 'PIPS2A_SPOTTR_053122_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_SPOTTR_053122_PIPS1B_10s.nc',
                          'parsivel_combined_SPOTTR_053122_PIPS2A_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_SPOTTR_053122_PIPS1B.nc',
                          'conventional_raw_SPOTTR_053122_PIPS2A.nc'],
    'start_times': ['20220531230000'] * 2,
    'end_times': ['20220601000500'] * 2,
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
    'plot_retrieval': True,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2017/NEXRAD/IOP_1A/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO', 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170325173000',
    'radar_end_timestamp': '20170325184500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}