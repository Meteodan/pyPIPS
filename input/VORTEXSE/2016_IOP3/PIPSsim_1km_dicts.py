# List of field names to read
fieldnames = ['dBZ']
# Requested elevation angle
el_req = 0.5

init_radar_dict = {
    '033116': {
        'radardir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2016/NEXRAD/IOP_3/KGWX/CFRadial',
        'fieldnames': fieldnames,
        'el_req': el_req,
        'radstarttimestamp': '20160331220000',
        'radstoptimestamp': '20160331234000'
    }
}

init_dis_dict = {
    '033116': {
        'dis_dir': '/Users/dawson29/sshfs_mounts/depot/data/Projects/VORTEXSE/obsdata/2016/PIPS/processed/IOP3',
        'dis_types': ['PIPS', 'PIPS', 'PIPS', 'PIPS'],
        'dis_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS2B'],
        'disfilenames': ['PIPS_1A_IOP_3_D1.txt', 'PIPS_1B_IOP_3_D1.txt', 'PIPS_2A_IOP_3_D1.txt',
                         'PIPS_2B_IOP_3_D1.txt'],
        'convfilenames': ['PIPS_1A_IOP_3_D1.txt', 'PIPS_1B_IOP_3_D1.txt', 'PIPS_2A_IOP_3_D1.txt',
                          'PIPS_2B_IOP_3_D1.txt'],
        'starttimes': ['20160331220000', '20160331220000', '20160331220000', '20160331220000'],
        'stoptimes': ['20160331234000', '20160331234000', '20160331234000', '20160331234000'],
        'interval': 60.
    }
}

init_model_dict = {
    '033116': {
        'runname': '1km453x453_newse',
        'nens': 36,
        'fileformat': 'hdf',
        'microphys': 'ZVD',
        'model_dt': 60.,
        'model_dt_mean': 900.,
        'basedirname': '/Users/dawson29/sshfs_mounts/rice_scratch/VORTEXSE/simulations/ARPS/2016_IOP3/EnKF/1km453x453_newse',
        'timestamp_model_init': '20160331180000',
        'timestamp_model_start': '20160331220000',
        'timestamp_model_stop': '20160331234000',
    }
}
