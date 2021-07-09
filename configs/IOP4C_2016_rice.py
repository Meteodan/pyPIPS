"""Configuration script for VORTEX-SE IOP4C 2016 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4C_2016',
    'deployment_names': ['IOP4C_D1_2016'] * 4,
    'PIPS_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset/SATP_retr',
    'plot_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset/SATP_retr/plots',
    'PIPS_types': ['PIPS'] * 4,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1A_2016_IOP4C_D1.txt', 'PIPS1B_2016_IOP4C_D1.txt',
                       'PIPS2A_2016_IOP4C_D1.txt', 'PIPS2B_2016_IOP4C_D1.txt'],
    'start_times': ['20160430193000'] * 4,
    'end_times': ['20160430224500'] * 4,
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
    'radar_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_4C/GWX/CFRadial/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160430190000',
    'radar_end_timestamp': '20160430221000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}