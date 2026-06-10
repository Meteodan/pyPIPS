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
RUN_MERGE_REALTIME=1
RUN_MANUAL_QC=0
RUN_APPLY_QC=1
RUN_CALC_DERIVED=1
RUN_MM_FITS=1
RUN_RADAR_INTERP=0  # Optional - disabled by default

# Default analysis options
QC_TAGS="qc roqc hoqc"
MM_MOMENT_COMBOS="24 234 246 346"
OUTPUT_TAG=""
REALTIME_DIR=""
MANUAL_QC_FILE=""
WRITE_ALIGNED_TRIM_TIMES=0
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
Usage: $0 <case_config_path_or_glob> [options]

Required Arguments:
    case_config_path_or_glob    Path to case configuration file or glob pattern
                                (e.g., configs/my_deployment.py or configs/ICECHIP_IOP*.py)

Options:
    --skip-csv2nc       Skip CSV to netCDF conversion
    --skip-merge        Skip real-time/card data merging
    --skip-manual-qc    Skip manual QC application
    --skip-qc           Skip quality control application
    --skip-derived      Skip derived parameter calculation
    --skip-mm-fits      Skip method of moments fitting
    --enable-radar      Enable radar interpolation (disabled by default)
    --enable-manual-qc  Enable manual QC application (disabled by default)
    --write-aligned-trim-times  Write parsivel-aligned trim times back to manual QC JSON (optional)

    --qc-tags TAGS      Space-separated QC tags to process (default: "qc")
    --mm-combos COMBOS  Space-separated moment combinations (default: "23 34 246")
    --output-tag TAG    Output file tag to distinguish from originals
    --realtime-dir DIR  Directory containing real-time netCDF files (required for merge step)
    --manual-qc-file FILE        Manual QC JSON decisions file (required if --enable-manual-qc)

    --log-dir DIR       Enable logging and specify directory for log files
    By default, runs: CSV→netCDF → merge realtime → QC → derived params → MM fits
    Manual QC is optional and runs within the QC step when enabled with --enable-manual-qc
    Radar interpolation is optional and must be enabled with --enable-radar

Examples:
    # Standard full analysis (single config)
    $0 configs/IOP1_2016.py

    # Process all ICECHIP IOPs sequentially
    $0 "configs/ICECHIP_IOP*.py"

    # Process specific IOP range
    $0 "configs/ICECHIP_IOP[1-5]*.py"

    # Skip CSV conversion (already have netCDF files)
    $0 configs/IOP1_2016.py --skip-csv2nc

    # Only run QC and derived parameters on multiple configs
    $0 "configs/ICECHIP_*.py" --skip-csv2nc --skip-mm-fits

    # Full analysis including radar interpolation
    $0 configs/IOP1_2016.py --enable-radar --output-tag "with_radar"

    # Analysis with real-time data merging on multiple configs
    $0 "configs/ICECHIP_IOP*.py" --realtime-dir /path/to/realtime/data

    # Enable manual QC (overwrites original files by default)
    $0 configs/IOP1_2016.py --enable-manual-qc --manual-qc-file configs/manual_qc_decisions.json

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
        --skip-csv2nc)
            RUN_CSV_TO_NC=0
            shift
            ;;
        --skip-merge)
            RUN_MERGE_REALTIME=0
            shift
            ;;
        --skip-manual-qc)
            RUN_MANUAL_QC=0
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
        --enable-manual-qc)
            RUN_MANUAL_QC=1
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
        --realtime-dir)
            REALTIME_DIR="$2"
            shift 2
            ;;
        --manual-qc-file)
            MANUAL_QC_FILE="$2"
            shift 2
            ;;
        --write-aligned-trim-times)
            WRITE_ALIGNED_TRIM_TIMES=1
            shift
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

# Construct common arguments
COMMON_ARGS="$CASE_CONFIG_PATH"
if [ -n "$OUTPUT_TAG" ]; then
    OUTPUT_ARG="--output-tag $OUTPUT_TAG"
else
    OUTPUT_ARG=""
fi

# Validate realtime-dir if merge step is enabled
if [ $RUN_MERGE_REALTIME -eq 1 ] && [ -z "$REALTIME_DIR" ]; then
    print_error "--realtime-dir is required when merge step is enabled"
    print_error "Use --skip-merge to disable merging or provide --realtime-dir path"
    exit 1
fi

if [ $RUN_MERGE_REALTIME -eq 1 ] && [ ! -d "$REALTIME_DIR" ]; then
    print_error "Real-time directory not found: $REALTIME_DIR"
    exit 1
fi

# Validate manual QC inputs if enabled
if [ $RUN_MANUAL_QC -eq 1 ] && [ -z "$MANUAL_QC_FILE" ]; then
    print_error "--manual-qc-file is required when manual QC is enabled"
    print_error "Use --skip-manual-qc to disable manual QC or provide --manual-qc-file path"
    exit 1
fi

if [ $RUN_MANUAL_QC -eq 1 ] && [ ! -f "$MANUAL_QC_FILE" ]; then
    print_error "Manual QC file not found: $MANUAL_QC_FILE"
    exit 1
fi

if [ $RUN_MANUAL_QC -eq 1 ] && [ $RUN_APPLY_QC -ne 1 ]; then
    print_error "Manual QC is integrated into apply_QC.py"
    print_error "Enable QC step or disable manual QC with --skip-manual-qc"
    exit 1
fi

if [ $WRITE_ALIGNED_TRIM_TIMES -eq 1 ] && [ $RUN_MANUAL_QC -ne 1 ]; then
    print_error "--write-aligned-trim-times requires manual QC to be enabled"
    print_error "Use --enable-manual-qc or remove --write-aligned-trim-times"
    exit 1
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

# Print analysis plan
echo "==============================================="
echo "pyPIPS Initial Analysis Workflow"
echo "==============================================="
echo "Config file: $CASE_CONFIG_PATH"
echo "Steps to run:"
[ $RUN_CSV_TO_NC -eq 1 ] && echo "  ✓ CSV to netCDF conversion"
[ $RUN_MERGE_REALTIME -eq 1 ] && echo "  ✓ Real-time/card data merging"
[ $RUN_APPLY_QC -eq 1 ] && echo "  ✓ Quality control application"
[ $RUN_MANUAL_QC -eq 1 ] && echo "  ✓ Manual QC (within QC step)"
[ $RUN_CALC_DERIVED -eq 1 ] && echo "  ✓ Derived parameter calculation"
[ $RUN_MM_FITS -eq 1 ] && echo "  ✓ Method of moments fitting"
[ $RUN_RADAR_INTERP -eq 1 ] && echo "  ✓ Radar interpolation"
echo "QC tags: $QC_TAGS"
echo "MM moment combinations: $MM_MOMENT_COMBOS"
[ -n "$OUTPUT_TAG" ] && echo "Output tag: $OUTPUT_TAG"
[ $RUN_MERGE_REALTIME -eq 1 ] && [ -n "$REALTIME_DIR" ] && echo "Real-time directory: $REALTIME_DIR"
[ $ENABLE_LOGGING -eq 1 ] && echo "Logging directory: $LOG_DIR"
echo "==============================================="
echo

# Flag to track if this config succeeds
CONFIG_SUCCESS=1

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
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "CSV to netCDF conversion completed"
        else
            print_error "CSV to netCDF conversion failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Step 2: Merge real-time and card data
if [ $RUN_MERGE_REALTIME -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Merging real-time and card-based data..."

    cmd="python ${ANALYSIS_DIR}/merge_PIPS_realtime.py $COMMON_ARGS --realtime-dir $REALTIME_DIR --diagnostic-plots $OUTPUT_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_merge.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_merge.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Real-time/card data merging completed"
        else
            print_error "Real-time/card data merging failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Real-time/card data merging completed"
        else
            print_error "Real-time/card data merging failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Step 3: Apply Quality Control
if [ $RUN_APPLY_QC -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
    print_step "Applying quality control filters..."

    MANUAL_QC_ARG=""
    if [ $RUN_MANUAL_QC -eq 1 ]; then
        MANUAL_QC_ARG="--manual-qc-file $MANUAL_QC_FILE"
    fi
    MANUAL_QC_ALIGN_ARG=""
    if [ $WRITE_ALIGNED_TRIM_TIMES -eq 1 ]; then
        MANUAL_QC_ALIGN_ARG="--write-aligned-trim-times"
    fi

    cmd="python ${ANALYSIS_DIR}/apply_QC.py $COMMON_ARGS --slowtemp-bias-correction --slowtemp-bias-config ../configs/ICECHIP_slowtemp_bias.py --compass-qc --compass-abs-threshold 10.0 --plot-compass-qc --slowtemp-qc --slowtemp-diff-threshold 5.0 --plot-slowtemp-qc --recompute-dewpoint-with-fallback --output-QC-tags $QC_TAGS $OUTPUT_ARG $MANUAL_QC_ARG $MANUAL_QC_ALIGN_ARG"
    [ $VERBOSE -eq 1 ] && echo "Command: $cmd"

    if [ $ENABLE_LOGGING -eq 1 ]; then
        LOG_STDOUT="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_qc.stdout"
        LOG_STDERR="${LOG_DIR}/${CONFIG_NAME}_${TIMESTAMP}_qc.stderr"
        print_step "Logs: $LOG_STDOUT, $LOG_STDERR"

        if eval "$cmd" > "$LOG_STDOUT" 2> "$LOG_STDERR"; then
            print_success "Quality control application completed"
        else
            print_error "Quality control application failed (see logs)"
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Quality control application completed"
        else
            print_error "Quality control application failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Step 4: Calculate derived parameters
if [ $RUN_CALC_DERIVED -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
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
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Derived parameter calculation completed"
        else
            print_error "Derived parameter calculation failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Step 5: Calculate method of moments fits
if [ $RUN_MM_FITS -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
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
            CONFIG_SUCCESS=0
        fi
    else
        if eval $cmd; then
            print_success "Method of moments fitting completed"
        else
            print_error "Method of moments fitting failed"
            CONFIG_SUCCESS=0
        fi
    fi
    echo
fi

# Step 6: Radar interpolation (optional)
if [ $RUN_RADAR_INTERP -eq 1 ] && [ $CONFIG_SUCCESS -eq 1 ]; then
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
        print_success "pyPIPS initial analysis workflow completed!"
    else
        print_error "pyPIPS initial analysis workflow failed!"
    fi
    echo "==============================================="
    echo "Config file: $CASE_CONFIG_PATH"
else
    echo "pyPIPS Batch Analysis Summary"
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

echo "Output files should be in the directory specified in your config file"
[ -n "$OUTPUT_TAG" ] && echo "Look for files with tag: $OUTPUT_TAG"
echo

# Exit with error if any configs failed
if [ ${#FAILED_CONFIGS[@]} -gt 0 ]; then
    exit 1
fi
