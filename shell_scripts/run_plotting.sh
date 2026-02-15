#!/bin/bash

# run_plotting.sh
#
# Shell script to generate standard plots from PIPS data.
# This script runs various plotting scripts on processed netCDF files.
#
# Usage: ./run_plotting.sh <path/to/case/config/file.py> [options]

set -e  # Exit on any error

# Default settings - toggle individual plotting scripts on/off
RUN_CONV_METEOGRAMS=1
RUN_DSD_METEOGRAMS=1
RUN_DSD_PLOTS=1
RUN_PIPS_DIAG=1
RUN_RADAR_PPI=0  # Optional - disabled by default (requires radar data)
RUN_VEL_D_PLOTS=1

# Default plotting options
QC_TAGS="qc"
ND_TAGS="qc"
PLOT_START_TIME=""
PLOT_END_TIME=""
PLOT_DIR=""
PLOT_CONFIG_PATH="plot_config.py"
IMAGE_FMT="png"
VERBOSE=0
ENABLE_LOGGING=0
LOG_DIR=""

# DSD plot options
PLOT_RAW=0
PLOT_FULL=1
PLOT_SERIES=0
PLOT_MM_FITS=0

# Vel-D plot options
PLOT_VD_RAW=0
PLOT_VD_QC=1
PLOT_VD_FULL=1
PLOT_VD_SERIES=0

# Radar PPI options
RADAR_EL_REQ=""
RADAR_FNAME_VARIANT="V06"
RADAR_INPUT_TAG=""

# Script paths
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
PLOT_SCRIPT_DIR="${SCRIPT_DIR}/../plotting_scripts"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Helper functions
print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_usage() {
    cat << EOF
Usage: $0 <case_config_path_or_glob> [options]

Required Arguments:
    case_config_path_or_glob    Path to case configuration file or glob pattern
                                (e.g., configs/my_deployment.py or configs/ICECHIP_IOP*.py)

Options - Script Selection:
    --skip-conv-meteograms    Skip conventional variable meteograms
    --skip-dsd-meteograms     Skip DSD meteograms
    --skip-dsd-plots          Skip DSD histograms
    --skip-diag               Skip diagnostic plots
    --skip-vel-d              Skip velocity-diameter plots
    --enable-radar-ppi        Enable radar PPI plots (disabled by default)

Options - General Plotting:
    --qc-tags TAGS            Space-separated QC tags for DSD meteograms (default: "qc")
    --nd-tags TAGS            Space-separated ND tags for DSD plots (default: "qc")
    --plot-start-time TIME    Start time for plots (YYYYmmDDHHMMSS format)
    --plot-end-time TIME      End time for plots (YYYYmmDDHHMMSS format)
    --plot-dir DIR            Directory to store plots (overrides config file)
    --plot-config PATH        Path to plot configuration file (default: plot_config.py)
    --image-fmt FMT           Image format: png, pdf, eps (default: png)

Options - DSD Plot Settings:
    --plot-raw                Include raw (non-QC) DSDs in DSD plots
    --no-plot-full            Don't plot full-deployment DSDs
    --plot-series             Plot time series of individual DSDs
    --plot-mm-fits            Include method of moments fits in DSD plots

Options - Velocity-Diameter Plot Settings:
    --plot-vd-raw             Plot raw velocity-diameter matrices
    --no-plot-vd-qc           Don't plot QC'd velocity-diameter matrices
    --no-plot-vd-full         Don't plot full-deployment v-d matrices
    --plot-vd-series          Plot time series of v-d matrices

Options - Radar PPI Settings:
    --radar-el-req ANGLE      Radar elevation angle (degrees)
    --radar-fname-variant VAR Radar filename variant (default: V06)
    --radar-input-tag TAG     Radar input tag for file selection

Options - Workflow:
    --log-dir DIR             Enable logging and specify directory for log files
    --verbose                 Print commands being executed

Examples:
    # Standard plotting for single config
    $0 configs/IOP1_2016.py

    # Plot all ICECHIP IOPs with custom QC tags
    $0 "configs/ICECHIP_IOP*.py" --qc-tags "qc roqc hoqc"

    # Only plot meteograms (skip DSD plots)
    $0 configs/IOP1_2016.py --skip-dsd-plots --skip-vel-d

    # Full DSD analysis with MM fits
    $0 configs/IOP1_2016.py --plot-mm-fits --plot-series

    # Include radar PPI plots
    $0 configs/IOP1_2016.py --enable-radar-ppi --radar-el-req 0.5

    # Custom time range for specific event
    $0 configs/IOP1_2016.py --plot-start-time 20160101120000 --plot-end-time 20160101180000

Note: When using glob patterns, enclose them in quotes to prevent premature shell expansion.

EOF
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    print_error "No case config file provided"
    print_usage
    exit 1
fi

CASE_CONFIG_PATTERN="$1"
shift

# Expand glob pattern to get list of config files
shopt -s nullglob  # Make glob return empty array if no matches
CONFIG_FILES=($CASE_CONFIG_PATTERN)
shopt -u nullglob

# Check if any config files were found
if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
    print_error "No config files found matching: $CASE_CONFIG_PATTERN"
    exit 1
fi

# Report number of configs found
NUM_CONFIGS=${#CONFIG_FILES[@]}
if [ $NUM_CONFIGS -eq 1 ]; then
    print_step "Found 1 config file: ${CONFIG_FILES[0]}"
else
    print_step "Found $NUM_CONFIGS config files matching pattern: $CASE_CONFIG_PATTERN"
fi

# Parse remaining options
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-conv-meteograms)
            RUN_CONV_METEOGRAMS=0
            shift
            ;;
        --skip-dsd-meteograms)
            RUN_DSD_METEOGRAMS=0
            shift
            ;;
        --skip-dsd-plots)
            RUN_DSD_PLOTS=0
            shift
            ;;
        --skip-diag)
            RUN_PIPS_DIAG=0
            shift
            ;;
        --skip-vel-d)
            RUN_VEL_D_PLOTS=0
            shift
            ;;
        --enable-radar-ppi)
            RUN_RADAR_PPI=1
            shift
            ;;
        --qc-tags)
            QC_TAGS="$2"
            shift 2
            ;;
        --nd-tags)
            ND_TAGS="$2"
            shift 2
            ;;
        --plot-start-time)
            PLOT_START_TIME="$2"
            shift 2
            ;;
        --plot-end-time)
            PLOT_END_TIME="$2"
            shift 2
            ;;
        --plot-dir)
            PLOT_DIR="$2"
            shift 2
            ;;
        --plot-config)
            PLOT_CONFIG_PATH="$2"
            shift 2
            ;;
        --image-fmt)
            IMAGE_FMT="$2"
            shift 2
            ;;
        --plot-raw)
            PLOT_RAW=1
            shift
            ;;
        --no-plot-full)
            PLOT_FULL=0
            shift
            ;;
        --plot-series)
            PLOT_SERIES=1
            shift
            ;;
        --plot-mm-fits)
            PLOT_MM_FITS=1
            shift
            ;;
        --plot-vd-raw)
            PLOT_VD_RAW=1
            shift
            ;;
        --no-plot-vd-qc)
            PLOT_VD_QC=0
            shift
            ;;
        --no-plot-vd-full)
            PLOT_VD_FULL=0
            shift
            ;;
        --plot-vd-series)
            PLOT_VD_SERIES=1
            shift
            ;;
        --radar-el-req)
            RADAR_EL_REQ="$2"
            shift 2
            ;;
        --radar-fname-variant)
            RADAR_FNAME_VARIANT="$2"
            shift 2
            ;;
        --radar-input-tag)
            RADAR_INPUT_TAG="$2"
            shift 2
            ;;
        --log-dir)
            ENABLE_LOGGING=1
            LOG_DIR="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=1
            shift
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Build common plotting arguments
COMMON_PLOT_ARGS=""
[ -n "$PLOT_START_TIME" ] && COMMON_PLOT_ARGS="$COMMON_PLOT_ARGS --plot-start-time $PLOT_START_TIME"
[ -n "$PLOT_END_TIME" ] && COMMON_PLOT_ARGS="$COMMON_PLOT_ARGS --plot-end-time $PLOT_END_TIME"
[ -n "$PLOT_DIR" ] && COMMON_PLOT_ARGS="$COMMON_PLOT_ARGS --plot-dir $PLOT_DIR"
[ -n "$PLOT_CONFIG_PATH" ] && COMMON_PLOT_ARGS="$COMMON_PLOT_ARGS --plot-config-path $PLOT_CONFIG_PATH"

# Arrays to track results across all configs
declare -a SUCCESSFUL_CONFIGS
declare -a FAILED_CONFIGS

# Main processing loop - iterate through all matching config files
for CASE_CONFIG_PATH in "${CONFIG_FILES[@]}"; do

# Validate that the config file exists (should always be true after glob expansion)
if [ ! -f "$CASE_CONFIG_PATH" ]; then
    print_error "Config file not found: $CASE_CONFIG_PATH"
    FAILED_CONFIGS+=("$CASE_CONFIG_PATH (file not found)")
    continue
fi

# Extract config name for display and logging
CONFIG_NAME=$(basename "$CASE_CONFIG_PATH" .py)

# Print header for this config (especially important when processing multiple)
if [ $NUM_CONFIGS -gt 1 ]; then
    echo
    echo "###############################################"
    echo "# Processing config [$((${#SUCCESSFUL_CONFIGS[@]} + ${#FAILED_CONFIGS[@]} + 1))/$NUM_CONFIGS]: $CONFIG_NAME"
    echo "###############################################"
    echo
fi

# Set up logging if enabled
if [ $ENABLE_LOGGING -eq 1 ]; then
    # Create log directory if it doesn't exist
    if [ ! -d "$LOG_DIR" ]; then
        mkdir -p "$LOG_DIR"
        if [ $? -ne 0 ]; then
            print_error "Failed to create log directory: $LOG_DIR"
            FAILED_CONFIGS+=("$CASE_CONFIG_PATH (log directory creation failed)")
            continue
        fi
    fi

    # Generate timestamp for log files
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

    if [ $NUM_CONFIGS -eq 1 ]; then
        print_step "Logging enabled. Log files will be saved to: $LOG_DIR"
        print_step "Log file prefix: ${CONFIG_NAME}_${TIMESTAMP}"
    fi
fi

# Print plotting plan
echo "==============================================="
echo "pyPIPS Plotting Workflow"
echo "==============================================="
echo "Config file: $CASE_CONFIG_PATH"
echo "Plots to generate:"
[ $RUN_CONV_METEOGRAMS -eq 1 ] && echo "  ✓ Conventional variable meteograms"
[ $RUN_DSD_METEOGRAMS -eq 1 ] && echo "  ✓ DSD meteograms"
[ $RUN_DSD_PLOTS -eq 1 ] && echo "  ✓ DSD histograms"
[ $RUN_PIPS_DIAG -eq 1 ] && echo "  ✓ Diagnostic plots"
[ $RUN_VEL_D_PLOTS -eq 1 ] && echo "  ✓ Velocity-diameter plots"
[ $RUN_RADAR_PPI -eq 1 ] && echo "  ✓ Radar PPI plots"
[ -n "$PLOT_START_TIME" ] && echo "Plot start time: $PLOT_START_TIME"
[ -n "$PLOT_END_TIME" ] && echo "Plot end time: $PLOT_END_TIME"
[ -n "$PLOT_DIR" ] && echo "Plot directory: $PLOT_DIR"
echo "Image format: $IMAGE_FMT"
echo "==============================================="
echo

# Flag to track if this config succeeds
CONFIG_SUCCESS=1

# Plot 1: Conventional variable meteograms
if [ $RUN_CONV_METEOGRAMS -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating conventional variable meteograms..."

    cmd="python ${PLOT_SCRIPT_DIR}/plot_conv_meteograms_nc.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_conv_meteograms.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_conv_meteograms.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Conventional meteograms completed"
        else
            print_error "Conventional meteograms failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Conventional meteograms completed"
        else
            print_error "Conventional meteograms failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Plot 2: DSD meteograms
if [ $RUN_DSD_METEOGRAMS -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating DSD meteograms..."

    # Loop through QC tags
    for QC_TAG in $QC_TAGS; do
        print_step "  Processing QC tag: $QC_TAG"

        cmd="python ${PLOT_SCRIPT_DIR}/plot_DSD_meteograms_nc.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS --QC-tag $QC_TAG"
        [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

        if [ $ENABLE_LOGGING -eq 1 ]; then
            LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_dsd_meteograms_${QC_TAG}.stdout"
            LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_dsd_meteograms_${QC_TAG}.stderr"

            if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
                print_success "DSD meteograms ($QC_TAG) completed"
            else
                print_error "DSD meteograms ($QC_TAG) failed (see logs)"
                CONFIG_SUCCESS=0
                break
            fi
        else
            if eval $cmd; then
                print_success "DSD meteograms ($QC_TAG) completed"
            else
                print_error "DSD meteograms ($QC_TAG) failed"
                CONFIG_SUCCESS=0
                break
            fi
        fi
    done
    echo
fi

# Plot 3: DSD histograms
if [ $RUN_DSD_PLOTS -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating DSD histograms..."

    cmd="python ${PLOT_SCRIPT_DIR}/plot_DSD_nc.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS --ND-tags $ND_TAGS --image-fmt $IMAGE_FMT"
    [ $PLOT_RAW -eq 1 ] && cmd="$cmd --plot-raw"
    [ $PLOT_FULL -eq 1 ] && cmd="$cmd --plot-full"
    [ $PLOT_SERIES -eq 1 ] && cmd="$cmd --plot-series"
    [ $PLOT_MM_FITS -eq 1 ] && cmd="$cmd --plot-MM-fits"

    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_dsd_plots.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_dsd_plots.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "DSD histograms completed"
        else
            print_error "DSD histograms failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "DSD histograms completed"
        else
            print_error "DSD histograms failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Plot 4: Diagnostic plots
if [ $RUN_PIPS_DIAG -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating diagnostic plots..."

    cmd="python ${PLOT_SCRIPT_DIR}/plot_PIPS_diag.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_diag.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_diag.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Diagnostic plots completed"
        else
            print_error "Diagnostic plots failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Diagnostic plots completed"
        else
            print_error "Diagnostic plots failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Plot 5: Velocity-diameter plots
if [ $RUN_VEL_D_PLOTS -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating velocity-diameter plots..."

    cmd="python ${PLOT_SCRIPT_DIR}/plot_vel_D_nc.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS"
    [ $PLOT_VD_RAW -eq 1 ] && cmd="$cmd --plot-raw"
    [ $PLOT_VD_QC -eq 1 ] && cmd="$cmd --plot-qc"
    [ $PLOT_VD_FULL -eq 1 ] && cmd="$cmd --plot-full"
    [ $PLOT_VD_SERIES -eq 1 ] && cmd="$cmd --plot-series"
    [ -n "$QC_TAGS" ] && cmd="$cmd --QC-tag $(echo $QC_TAGS | awk '{print $1}')"  # Use first QC tag

    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_vel_d.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_vel_d.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Velocity-diameter plots completed"
        else
            print_error "Velocity-diameter plots failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Velocity-diameter plots completed"
        else
            print_error "Velocity-diameter plots failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Plot 6: Radar PPI plots (optional)
if [ $RUN_RADAR_PPI -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Generating radar PPI plots..."

    cmd="python ${PLOT_SCRIPT_DIR}/plot_radar_ppi.py $CASE_CONFIG_PATH $COMMON_PLOT_ARGS --image-fmt $IMAGE_FMT --fname-variant $RADAR_FNAME_VARIANT"
    [ -n "$RADAR_EL_REQ" ] && cmd="$cmd --el-req $RADAR_EL_REQ"
    [ -n "$RADAR_INPUT_TAG" ] && cmd="$cmd --input-tag $RADAR_INPUT_TAG"

    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_radar_ppi.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_radar_ppi.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Radar PPI plots completed"
        else
            print_warning "Radar PPI plots failed (this is optional, see logs)"
        fi
    else
        if eval $cmd; then
            print_success "Radar PPI plots completed"
        else
            print_warning "Radar PPI plots failed (this is optional)"
        fi
    fi
    echo
fi

# Track result for this config
if [ $CONFIG_SUCCESS -eq 1 ]; then
    SUCCESSFUL_CONFIGS+=("$CONFIG_NAME")
else
    FAILED_CONFIGS+=("$CONFIG_NAME")
fi

done  # End of config loop

# Final summary
echo
echo "==============================================="
if [ $NUM_CONFIGS -eq 1 ]; then
    if [ $CONFIG_SUCCESS -eq 1 ]; then
        print_success "pyPIPS plotting workflow completed!"
    else
        print_error "pyPIPS plotting workflow failed!"
    fi
    echo "==============================================="
    echo "Config file: $CASE_CONFIG_PATH"
else
    echo "pyPIPS Batch Plotting Summary"
    echo "==============================================="
    echo "Total configs processed: $NUM_CONFIGS"
    echo "Successful: ${#SUCCESSFUL_CONFIGS[@]}"
    echo "Failed: ${#FAILED_CONFIGS[@]}"
    echo

    if [ ${#SUCCESSFUL_CONFIGS[@]} -gt 0 ]; then
        print_success "Successfully processed configs:"
        for config in "${SUCCESSFUL_CONFIGS[@]}"; do
            echo "  ✓ $config"
        done
    fi

    if [ ${#FAILED_CONFIGS[@]} -gt 0 ]; then
        echo
        print_error "Failed configs:"
        for config in "${FAILED_CONFIGS[@]}"; do
            echo "  ✗ $config"
        done
    fi
    echo "==============================================="
fi

[ -n "$PLOT_DIR" ] && echo "Plots saved to: $PLOT_DIR"
echo "Plot format: $IMAGE_FMT"
echo

# Exit with error if any configs failed
if [ ${#FAILED_CONFIGS[@]} -gt 0 ]; then
    exit 1
fi
