#!/bin/bash

python read_PIPS_real_time.py PIPS1A --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS1A_read.log &
python read_PIPS_real_time.py PIPS1B --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS1B_read.log &
python read_PIPS_real_time.py PIPS2A --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS2A_read.log &
python read_PIPS_real_time.py PIPS2B --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS2B_read.log &
python read_PIPS_real_time.py PIPS3A --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS3A_read.log &
python read_PIPS_real_time.py PIPS3B --base-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime &> PIPS3B_read.log &