#!/bin/bash

python read_PIPS_real_time.py PIPS1A --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS1A_plot.log &
python read_PIPS_real_time.py PIPS1B --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS1B_plot.log &
python read_PIPS_real_time.py PIPS2A --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS2A_plot.log &
python read_PIPS_real_time.py PIPS2B --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS2B_plot.log &
python read_PIPS_real_time.py PIPS3A --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS3A_plot.log &
python read_PIPS_real_time.py PIPS3B --netcdf-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/netcdf &> PIPS3B_plot.log &
