"""Configuration script for ICECHIP 2025 IOP19 deployment"""

PIPS_IO_dict = {
    'dataset_name': 'IOP19_061725',
    'deployment_names': ['IOP19_061725'] * 5,
    'input_txt_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP19_061725/csv/',
    'PIPS_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP19_061725/netcdf/',
    'plot_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/PIPS_data/IOP19_061725/plots/10s/',
    'PIPS_types': ['PIPS'] * 5,
    'probe_set': 'ICECHIP_2025_B',
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames': ['PIPS1A_IOP19_061725_merged.txt', 'PIPS1B_IOP19_061725_merged.txt',
                       'PIPS2A_IOP19_061725_merged.txt', 'PIPS3A_IOP19_061725_merged.txt',
                       'PIPS3B_IOP19_061725_merged.txt'],
    'PIPS_filenames_nc': ['parsivel_combined_IOP19_061725_PIPS1A_10s.nc',
                          'parsivel_combined_IOP19_061725_PIPS1B_10s.nc',
                          'parsivel_combined_IOP19_061725_PIPS2A_10s.nc',
                          'parsivel_combined_IOP19_061725_PIPS3A_10s.nc',
                          'parsivel_combined_IOP19_061725_PIPS3B_10s.nc'],
    'conv_filenames_nc': ['conventional_raw_IOP19_061725_PIPS1A.nc',
                          'conventional_raw_IOP19_061725_PIPS1B.nc',
                          'conventional_raw_IOP19_061725_PIPS2A.nc',
                          'conventional_raw_IOP19_061725_PIPS3A.nc',
                          'conventional_raw_IOP19_061725_PIPS3B.nc'],
    'start_times': ['20250618000500'] * 5,
    'end_times': ['20250618030000'] * 5,
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
    'comp_radar': False,
    'calc_dualpol': True,
    'plot_retrieval': False,
    'radar_name': 'KLBB',
    'radar_type': 'NEXRAD',
    'radar_dir': '/Users/dawson29/Dropbox/Projects/ICECHIP/obsdata/NEXRAD/IOP12_060625/KLBB/CFRadial/',
    # 'radar_fname_pattern': '{rad_name}{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}_V06.nc',
    'radar_fname_pattern': 'cfrad.{year:04d}{month:02d}{day:02d}_{hour:02d}{min:02d}{sec:02d}.{dum1:03d}_to_{year2:04d}{month2:02d}{day2:02d}_{hour2:02d}{min2:02d}{sec2:02d}.{dum2:03d}_{rad_name}_SUR.nc',
    'field_names': ['REF', 'ZDR', 'RHO'], # , 'Dm_Z01', 'sigma_Z01', 'RR_Z01', 'mu_Z01', 'lamda_Z01'],
    'el_req': 0.5,
    'radar_start_timestamp': '20250606230000',
    'radar_end_timestamp': '20250607020000',
    'scatt_dir': '/Users/dawson29/Projects/pyPIPS/tmatrix/S-Band/',
    'wavelength': 10.7
}
