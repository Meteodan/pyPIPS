"""Configuration script for VORTEX-SE IOP4B 2016 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4B_2016',
    'deployment_names': ['IOP4B_D1_2016'],
    'PIPS_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_new',
    'plot_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/full_PIPS_dataset_new/plots',
    'PIPS_types': ['PIPS'],
    'PIPS_names': ['PIPS1A'],
    'PIPS_filenames': ['PIPS1A_2016_IOP4A_4D_D1.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP4A_4D_D1_2016_PIPS1A_60s.nc'],
    'start_times': ['20160429211500'],
    'end_times': ['20160430030000'],
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
    'radar_name': 'KHTX',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Volumes/scr_fast/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_4A_4C/HTX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO', 'Dm', 'sigma', 'RR', 'mu', 'lamda'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160429211500',
    'radar_end_timestamp': '20160430030000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}