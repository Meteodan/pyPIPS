#!/bin/bash
export HDF5_USE_FILE_LOCKING=FALSE
source /Users/dawson29/mambaforge/bin/activate
conda activate pyPIPS
cd /Users/dawson29/Projects/pyPIPS/realtime_scripts/
python plot_TriPIPS_real_time.py
