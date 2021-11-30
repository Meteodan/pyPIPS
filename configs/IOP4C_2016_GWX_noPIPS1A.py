"""Configuration script for VORTEX-SE IOP4C 2016 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4C_2016',
    'deployment_names': ['IOP4C_D1_2016'] * 3,
    'input_txt_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_links',
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15/plots',
    'PIPS_types': ['PIPS'] * 3,
    'PIPS_names': ['PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1B_2016_IOP4C_D1.txt',
                       'PIPS2A_2016_IOP4C_D1.txt', 'PIPS2B_2016_IOP4C_D1.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP4C_D1_2016_PIPS1B_60s.nc',
                          'parsivel_combined_IOP4C_D1_2016_PIPS2A_60s.nc',
                          'parsivel_combined_IOP4C_D1_2016_PIPS2B_60s.nc'],
    'start_times': ['20160430204500', '20160430193000', '20160430204500'] * 3,
    'end_times': ['20160430214500', '20160430200000', '20160430214500'] * 3,
    'requested_interval': 60.
}

PIPS_qc_dict = {
    'strongwindQC': True,
    'splashingQC': True,
    'marginQC': True,
    'rainfallQC': True,
    'rainonlyQC': True,
    'hailonlyQC': False,
    'graupelonlyQC': False,
    'basicQC': False,
}

radar_config_dict = {
    'load_radar_at_PIPS': True,
    'save_radar_at_PIPS': False,
    'comp_radar': True,
    'clean_radar': False,
    'calc_dualpol': True,
    'plot_retrieval': True,
    'radar_name': 'KGWX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_4C/GWX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO', 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160430193000',
    'radar_end_timestamp': '20160430214500',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}