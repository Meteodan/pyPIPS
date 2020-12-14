"""Configuration script for VORTEX-SE IOP4C 2017 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4C_2017',
    'deployment_names': ['IOP4C_D1_2017'] * 1,
    'PIPS_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset_RB15/',
    'plot_dir': '/depot/dawson29/data/Projects/VORTEXSE/simulations/ARPS/2017_IOP4C/EnKF/PIPS/New/plots/full_time/CCN300/',
    'PIPS_types': ['PIPS'] * 1,
    'PIPS_names': ['PIPS2A'],
    'PIPS_filenames': ['PIPS2A_FMCW_043017.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_FMCW_2017_043017_PIPS2A_60s.nc'],
    'start_times': [None] * 1,
    'end_times': [None] * 1,
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
    'radar_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/2017/NEXRAD/IOP_4C/level2/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170430180000',
    'radar_end_timestamp': '20170430220000',
    'scatt_dir': '/home/cbelak/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}

model_config_dict = {
    'runname': 'CCN_750_1km243x243_3km153x153_043017_NAM',
    'nens': 40,
    'fileformat': 'hdf',
    'microphys': 'ZVD',
    'model_dt': 60.,
    'model_dt_mean': 300.,
    'cycle': 'prior',
    'basedirname': '/scratch/rice/c/cbelak/Projects/VORTEXSE/simulations/ARPS/2017_IOP4C/EnKF/1km243x243_3km153x153_043017_CCN750_New/',
    'timestamp_model_init': '20170430060000',
    'timestamp_model_start': '20170430200000',
    'timestamp_model_stop': '20170430220000',
    'nproc_x': 6,
    'nproc_y': 6,
}
