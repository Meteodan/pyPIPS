"""Configuration script for VORTEX-SE IOP3 2016"""

PIPS_IO_dict = {
    'deployment_name': 'IOP3_2016',
    'PIPS_dir': '/Users/ddawson/Dropbox/Projects/VORTEXSE/obs_data/PIPS/2016/IOP3',
    'plot_dir': '/Users/dawson29/pyPIPS_work/VORTEXSE/2016/IOP3/plots',
    'PIPS_types': ['PIPS', 'PIPS', 'PIPS', 'PIPS'],
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
    'PIPS_filenames': ['PIPS_1A_IOP_3_D1.txt', 'PIPS_1B_IOP_3_D1.txt', 'PIPS_2A_IOP_3_D1.txt',
                       'PIPS_2B_IOP_3_D1.txt'],
    'start_times': ['20160331220000', '20160331220000', '20160331220000', '20160331220000'],
    'end_times': ['20160331234000', '20160331234000', '20160331234000', '20160331234000'],
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
    'radar_name': 'KGWX',
    'radar_dir': '/Users/ddawson/Dropbox/Projects/VORTEXSE/obs_data/NEXRAD/IOP3_2016/KGWX/CFRadial',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20160331220000',
    'radar_end_timestamp': '20160331234000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}