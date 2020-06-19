"""Configuration script for VORTEX-SE IOP4A-4D 2016 PIPS1A deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4A_4D_2016_PIPS1A',
    'deployment_names': ['IOP4A_4D_D1_2016'],
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_new',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_new/plots',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS1A'],
    'PIPS_filenames': ['PIPS1A_2016_IOP4A_4D_D1.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP4A_4D_D1_2016_PIPS1A_60s.nc'],
    'start_times': [None],
    'end_times': [None],
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
    'comp_radar': False,
    'clean_radar': False,
    'calc_dualpol': False,
    'plot_retrieval': False,
    'radar_name': 'KHTX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_4A_4C/HTX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO', 'Dm', 'sigma', 'RR', 'mu', 'lamda'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160430120000',
    'radar_end_timestamp': '20160430235959',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}