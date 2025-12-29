#!/bin/bash

# run_initial_analysis.sh
#
# Shell script to perform initial PIPS deployment analysis starting from CSV files.
# This script orchestrates the standard pyPIPS analysis workflow:
# CSV → netCDF → QC → derived parameters → MM fits → (optional) radar interpolation
#
# Usage: ./run_initial_analysis.sh <path/to/case/config/file.py> [options]

set -e  # Exit on any error

# Default settings - toggle individual analysis steps on/off
RUN_CSV_TO_NC=1
RUN_APPLY_QC=1
RUN_CALC_DERIVED=1
RUN_MM_FITS=1
RUN_RADAR_INTERP=0  # Optional - disabled by default

# Default analysis options
QC_TAGS="qc roqc hoqc"
MM_MOMENT_COMBOS="24 234 246 346"
OUTPUT_TAG=""
VERBOSE=0
ENABLE_LOGGING=0
LOG_DIR=""

# Script paths
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
ANALYSIS_DIR="${SCRIPT_DIR}/../analysis_scripts"

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
Usage: $0 <case_config_path> [options]

Required Arguments:
    case_config_path    Path to case configuration file (e.g., configs/my_deployment.py)

Options:
    --skip-csv2nc       Skip CSV to netCDF conversion
    --skip-qc           Skip quality control application
    --skip-derived      Skip derived parameter calculation
    --skip-mm-fits      Skip method of moments fitting
    --enable-radar      Enable radar interpolation (disabled by default)

    --qc-tags TAGS      Space-separated QC tags to process (default: "qc")
    --mm-combos COMBOS  Space-separated moment combinations (default: "23 34 246")
    --output-tag TAG    Output file tag to distinguish from originals

    --log-dir DIR       Enable logging and specify directory for log files
    By default, runs: CSV→netCDF → QC → derived params → MM fits
    Radar interpolation is optional and must be enabled with --enable-radar

Examples:
    # Standard full analysis
    $0 configs/IOP1_2016.py

    # Skip CSV conversion (already have netCDF files)
    $0 configs/IOP1_2016.py --skip-csv2nc

    # Only run QC and derived parameters
    $0 configs/IOP1_2016.py --skip-csv2nc --skip-mm-fits

    # Full analysis including radar interpolation
    $0 configs/IOP1_2016.py --enable-radar --output-tag "with_radar"

EOF
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    print_error "No case config file provided"
    print_usage
    exit 1
fi

CASE_CONFIG_PATH="$1"
shift

# Check if config file exists
if [ ! -f "$CASE_CONFIG_PATH" ]; then
    print_error "Config file not found: $CASE_CONFIG_PATH"
    exit 1
fi

# Parse remaining options
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-csv2nc)
            RUN_CSV_TO_NC=0
            shift
            ;;
        --skip-qc)
            RUN_APPLY_QC=0
            shift
            ;;
        --skip-derived)
            RUN_CALC_DERIVED=0
            shift
            ;;
        --skip-mm-fits)
            RUN_MM_FITS=0
            shift
            ;;
        --enable-radar)
            RUN_RADAR_INTERP=1
            shift
            ;;
        --qc-tags)
            QC_TAGS="$2"
            shift 2
            ;;
        --mm-combos)
            MM_MOMENT_COMBOS="$2"
            shift 2
            ;;
        --output-tag)
            OUTPUT_TAG="$2"
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

# Construct common arguments
COMMON_ARGS="$CASE_CONFIG_PATH"
if [ -n "$OUTPUT_TAG" ]; then
    OUTPUT_ARG="--output-tag $OUTPUT_TAG"
else
    OUTPUT_ARG=""
fi

# Set up logging if enabled
if [ $ENABLE_LOGGING -eq 1 ]; then
    # Create log directory if it doesn't exist
    if [ ! -d "$LOG_DIR" ]; then
        mkdir -p "$LOG_DIR"
        if [ $? -ne 0 ]; then
            print_error "Failed to create log directory: $LOG_DIR"
            exit 1
        fi
    fi

    # Generate timestamp for log files
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    CONFIG_NAME=$(basename "$CASE_CONFIG_PATH" .py)

    print_step "Logging enabled. Log files will be saved to: $LOG_DIR"
    print_step "Log file prefix: ${CONFIG_NAME}_${TIMESTAMP}"
fi

# Print analysis plan
echo "==============================================="
echo "pyPIPS Initial Analysis Workflow"
echo "==============================================="
echo "Config file: $CASE_CONFIG_PATH"
echo "Steps to run:"
[ $RUN_CSV_TO_NC -eq 1 ] && echo "  ✓ CSV to netCDF conversion"
[ $RUN_APPLY_QC -eq 1 ] && echo "  ✓ Quality control application"
[ $RUN_CALC_DERIVED -eq 1 ] && echo "  ✓ Derived parameter calculation"
[ $RUN_MM_FITS -eq 1 ] && echo "  ✓ Method of moments fitting"
[ $RUN_RADAR_INTERP -eq 1 ] && echo "  ✓ Radar interpolation"
echo "QC tags: $QC_TAGS"
echo "MM moment combinations: $MM_MOMENT_COMBOS"
[ -n "$OUTPUT_TAG" ] && echo "Output tag: $OUTPUT_TAG"
[ $ENABLE_LOGGING -eq 1 ] && echo "Logging directory: $LOG_DIR"
echo "==============================================="
echo

# Step 1: Convert CSV to netCDF
if [ $RUN_CSV_TO_NC -eq 1 ]; then
    print_step "Converting CSV files to netCDF format..."

    cmd="python ${ANALYSIS_DIR}/PIPS_csv_to_nc.py $COMMON_ARGS --check-order --sort-times --fill-gaps $OUTPUT_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_csv2nc.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_csv2nc.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "CSV to netCDF conversion completed"
        else
            print_error "CSV to netCDF conversion failed (see logs)"
            exit 1
        fi
    else
        if eval $cmd; then
            print_success "CSV to netCDF conversion completed"
        else
            print_error "CSV to netCDF conversion failed"
            exit 1
        fi
    fi
    echo
fi

# Step 2: Apply Quality Control
if [ $RUN_APPLY_QC -eq 1 ]; then
    print_step "Applying quality control filters..."

    cmd="python ${ANALYSIS_DIR}/apply_QC.py $COMMON_ARGS --output-QC-tags $QC_TAGS $OUTPUT_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_qc.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_qc.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Quality control application completed"
        else
            print_error "Quality control application failed (see logs)"
            exit 1
        fi
    else
        if eval $cmd; then
            print_success "Quality control application completed"
        else
            print_error "Quality control application failed"
            exit 1
        fi
    fi
    echo
fi

# Step 3: Calculate derived parameters
if [ $RUN_CALC_DERIVED -eq 1 ]; then
    print_step "Calculating derived parameters..."

    cmd="python ${ANALYSIS_DIR}/calc_derived_params.py $COMMON_ARGS --QC-tags $QC_TAGS $OUTPUT_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_derived.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_derived.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Derived parameter calculation completed"
        else
            print_error "Derived parameter calculation failed (see logs)"
            exit 1
        fi
    else
        if eval $cmd; then
            print_success "Derived parameter calculation completed"
        else
            print_error "Derived parameter calculation failed"
            exit 1
        fi
    fi
    echo
fi

# Step 4: Calculate method of moments fits
if [ $RUN_MM_FITS -eq 1 ]; then
    print_step "Calculating method of moments DSD fits..."

    cmd="python ${ANALYSIS_DIR}/calc_MM_fits.py $COMMON_ARGS --QC-tags $QC_TAGS --moment-combos $MM_MOMENT_COMBOS"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_mmfits.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_mmfits.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Method of moments fitting completed"
        else
            print_error "Method of moments fitting failed (see logs)"
            exit 1
        fi
    else
        if eval $cmd; then
            print_success "Method of moments fitting completed"
        else
            print_error "Method of moments fitting failed"
            exit 1
        fi
    fi
    echo
fi

# Step 5: Radar interpolation (optional)
if [ $RUN_RADAR_INTERP -eq 1 ]; then
    print_step "Interpolating radar observations to PIPS locations..."

    cmd="python ${ANALYSIS_DIR}/radar_to_PIPS.py $COMMON_ARGS $OUTPUT_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_radar.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_radar.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Radar interpolation completed"
        else
            print_warning "Radar interpolation failed (this is optional, see logs)"
        fi
    else
        if eval $cmd; then
            print_success "Radar interpolation completed"
        else
            print_warning "Radar interpolation failed (this is optional)"
        fi
    fi
    echo
fi

echo "==============================================="
print_success "pyPIPS initial analysis workflow completed!"
echo "==============================================="
echo "Config file: $CASE_CONFIG_PATH"
echo "Output files should be in the directory specified in your config file"
[ -n "$OUTPUT_TAG" ] && echo "Look for files with tag: $OUTPUT_TAG"
echo