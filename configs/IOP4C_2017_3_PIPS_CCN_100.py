"""Configuration script for VORTEX-SE IOP4C 2017 deployment 1-km CCN 100"""

PIPS_IO_dict = {
    'dataset_name': 'IOP4C_2017',
    'deployment_names': ['IOP4C_D1_2017'] * 3,
    'PIPS_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/full_PIPS_dataset',
    'plot_dir': '/depot/dawson29/data/Projects/VORTEXSE/simulations/ARPS/2017_IOP4C/EnKF/PIPS/plots/CCN100',
    'PIPS_types': ['PIPS'] * 3,
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2B'],
    'PIPS_filenames': ['PIPS1A_2017_IOP4C_D1.txt', 'PIPS1B_2017_IOP4C_D1.txt',
                       'PIPS2B_2017_IOP4C_D1.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP4C_D1_2017_PIPS1A_60s.nc',
                          'parsivel_combined_IOP4C_D1_2017_PIPS1B_60s.nc',
                          'parsivel_combined_IOP4C_D1_2017_PIPS2B_60s.nc'],
    'start_times': [None] * 3,
    'end_times': [None] * 3,
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
    'radar_dir': '/depot/dawson29/data/Projects/VORTEXSE/obsdata/2017/NEXRAD/IOP_4C/CFRadial/modified/',
    'field_names': ['REF', 'ZDR', 'RHO'],
    'el_req': 0.5,
    'radar_start_timestamp': '20170430170000',
    'radar_end_timestamp': '20170430210500',
    'scatt_dir': '/home/cbelak/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}

model_config_dict = {
    'runname': '5_min_CCN_100_1km243x243_3km153x153_043017_NAM',
    'nens': 40,
    'fileformat': 'hdf',
    'microphys': 'ZVD',
    'model_dt': 60.,
    'model_dt_mean': 900.,
   'basedirname': '/scratch/rice/c/cbelak/Projects/VORTEXSE/simulations/ARPS/2017_IOP4C/EnKF/1km243x243_3km153x153_043017_CCN100',
    'timestamp_model_init': '20170430060000',
    'timestamp_model_start': '20170430183000',
    'timestamp_model_stop': '20170430193000',
    'nproc_x': 6,
    'nproc_y': 6,
}
