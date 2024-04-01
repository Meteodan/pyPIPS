#!/bin/bash

# python plot_PIPS_real_time.py PIPS1A --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS1A_plot.log &
# python plot_PIPS_real_time.py PIPS1B --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS1B_plot.log &
# python plot_PIPS_real_time.py PIPS2A --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS2A_plot.log &
# python plot_PIPS_real_time.py PIPS2B --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS2B_plot.log &
# python plot_PIPS_real_time.py PIPS3A --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS3A_plot.log &
# python plot_PIPS_real_time.py PIPS3B --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images &> PIPS3B_plot.log &

# Path to your Python script
PYTHON_SCRIPT="plot_PIPS_real_time.py"

# Array of command-line arguments for each script instance
# Each set of arguments for an instance should be a single string
# For example, if your script takes two arguments, an instance might look like: "arg1 arg2"
ARGUMENTS=("PIPS2B --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images"
           "PIPS3A --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images"
           "PIPS3B --netcdf-dir-url stormlab.eaps.purdue.edu/realtime/perils_realtime/netcdf --image-output-dir /Users/dawson29/sshfs_mounts/stormlab_web/perils_realtime/images"
)

# Array of PID files for each script instance
PIDFILES=("/Users/dawson29/temp/pidfiles/PIPS2B_plot.pid"
          "/Users/dawson29/temp/pidfiles/PIPS3A_plot.pid"
          "/Users/dawson29/temp/pidfiles/PIPS3B_plot.pid"
)

# Function to check and restart the process if needed
check_process() {
    local args=$1
    local pidfile=$2

    if [ ! -f $pidfile ]; then
        echo "PID file for instance with args '$args' not found. Starting the process."
        restart_process "$args" $pidfile
    else
        PID=$(cat $pidfile)
        if ! ps -p $PID > /dev/null 2>&1; then
            echo "Process for instance with args '$args' and PID $PID not found. Restarting."
            restart_process "$args" $pidfile
        else
            echo "Process for instance with args '$args' and PID $PID is running."
        fi
    fi
}

# Function to restart the process
restart_process() {
    local args=$1
    local pidfile=$2

    python $PYTHON_SCRIPT $args &
    PID=$!
    echo $PID > $pidfile
}

# First-time start-up for each script instance
for i in "${!ARGUMENTS[@]}"; do
    check_process "${ARGUMENTS[$i]}" "${PIDFILES[$i]}"
done

# Main loop to check the process status periodically
while true; do
    for i in "${!ARGUMENTS[@]}"; do
        check_process "${ARGUMENTS[$i]}" "${PIDFILES[$i]}"
    done
    sleep 60  # Check every 60 seconds
done
