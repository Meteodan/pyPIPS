# pyPIPS Development Progress Log

This file tracks significant development progress, lessons learned, and agent activities for pyPIPS development sessions.

## Session: Early 2026 - AI Development Infrastructure and Shell Workflow Creation

### Key Code Updates

#### 1. Created Comprehensive AI Coding Agent Instructions
**Problem**: Need structured guidelines for AI assistant to effectively work with pyPIPS codebase
**Solution**:
- Created comprehensive `.github/copilot-instructions.md` with project overview, architecture details, and workflow patterns
- Documented core data flow: Raw data → processing → QC → output
- Established module hierarchy and critical dependencies
- Defined configuration-driven analysis patterns and data conventions

**Files Modified**: `.github/copilot-instructions.md`

#### 2. Built Initial Analysis Shell Workflow
**Problem**: Need automated pipeline for standard pyPIPS analysis workflow
**Solution**:
- Created `run_initial_analysis.sh` with complete workflow: CSV→netCDF→QC→derived→MM fits
- Implemented modular step control with `RUN_*` flags and `--skip-*` options
- Added logging support, parameter validation, and comprehensive help documentation
- Established color-coded output and error handling patterns

**Files Modified**: `shell_scripts/run_initial_analysis.sh`

### Key Lessons Learned

#### Documentation-Driven Development
**Lesson**: Comprehensive AI assistant instructions significantly improve development effectiveness by providing context about architecture, patterns, and conventions.
**Impact**: Enables consistent development practices and reduces time spent understanding project structure.

#### Modular Workflow Design
**Lesson**: Analysis workflows benefit from modular design where each step can be independently enabled/disabled with clear parameters and validation.
**Impact**: Provides flexibility for different analysis scenarios and easier debugging of individual steps.

### Technical Context for Future Development
- Established foundation for configuration-driven analysis using `PIPS_IO_dict` patterns
- Created template for shell script integration with Python analysis modules
- Defined standards for command-line interface design and help documentation

---

## Session: Mid-2026 - Jupyter Notebook to Command-Line Script Conversion

### Key Code Updates

#### 1. Converted merge_PIPS_realtime_nc.ipynb to Production Script
**Problem**: Critical notebook functionality needed as command-line script for automation
**Solution**:
- Created `merge_PIPS_realtime.py` with complete notebook functionality
- Implemented argument parsing with required and optional parameters
- Added comprehensive logging and error handling
- Maintained notebook's data processing logic while adapting to command-line interface

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

#### 2. Implemented GPS Time Correction
**Problem**: Real-time logger timestamps need GPS synchronization
**Solution**:
- Added GPS time correction using coordinates from GPS_status='A' records
- Implemented fallback to "last good" GPS positions
- Created robust time synchronization handling for real-time data streams

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

#### 3. Added Configurable Output Options
**Problem**: Need flexible output file naming and directory control
**Solution**:
- Made output_tag default to empty string instead of 'merged'
- Added configurable output directories with fallback to PIPS_dir
- Implemented tag suffix logic for distinguishing output files

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

### Key Lessons Learned

#### Notebook to Script Conversion Patterns
**Lesson**: Converting notebooks requires careful preservation of data processing logic while restructuring for command-line interface and error handling.
**Impact**: Enables automation of previously interactive analyses while maintaining scientific accuracy.

#### GPS Time Synchronization Complexity
**Lesson**: Real-time data systems require robust time correction with GPS coordinates, including handling of GPS signal loss and coordinate validation.
**Impact**: Critical for accurate temporal alignment of meteorological measurements.

### Technical Context for Future Development
- Established patterns for real-time vs card data merging using xarray.combine_first()
- Created template for GPS time correction in real-time data processing
- Defined standards for configurable output naming and directory management

---

## Session: Late 2026 - Data Type Compatibility and Diagnostic Plotting

### Key Code Updates

#### 1. Resolved Critical dtype Promotion Error
**Problem**: TypeError during telegram data merging - Float64DType vs TimeDelta64DType incompatibility
**Solution**:
- Root cause: Inconsistent xarray loading parameters causing automatic timedelta conversion
- Fixed by using `decode_timedelta=False` consistently for all xarray.open_dataset() calls
- Prevented automatic conversion of sample_interval to timedelta64[ns] dtype
- Ensured dtype consistency between real-time and card datasets

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

#### 2. Implemented Comprehensive Diagnostic Plotting
**Problem**: Need visual verification of data merging quality and missing data patterns
**Solution**:
- Created `create_diagnostic_plots()` function with 6 different plot types
- Implemented missing data comparison plots
- Added DSD meteogram visualizations for merged, real-time, and card data
- Used pyPIPS plotmodule for consistent scientific visualization

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

#### 3. Enhanced Configuration Support
**Problem**: Need flexible resampling intervals and deployment-specific parameters
**Solution**:
- Added support for `requested_interval` from configuration files
- Implemented configurable resampling with proper time offset handling
- Added comprehensive configuration parameter validation

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

### Key Lessons Learned

#### xarray Loading Parameter Consistency
**Lesson**: When merging datasets from different sources, all xarray loading operations must use identical parameters to prevent artificial dtype mismatches.
**Impact**: Critical for avoiding data type promotion errors that can halt analysis workflows.

#### Diagnostic Plotting for Data Quality
**Lesson**: Visual diagnostics are essential for verifying complex data merging operations, especially when combining real-time and archived datasets.
**Impact**: Enables rapid identification of data quality issues and merge effectiveness.

#### Configuration-Driven Flexibility
**Lesson**: Analysis scripts should extract timing and resampling parameters from configuration files rather than hardcoding values.
**Impact**: Supports different deployment scenarios and instrument configurations without code changes.

### Technical Context for Future Development
- Established xarray loading standards for consistent dtype handling
- Created diagnostic plotting patterns for data quality verification
- Defined configuration parameter extraction and validation patterns

---

## Session: January 16-17, 2026 - Real-time Data Merging Integration

### Key Code Updates

#### 1. Fixed Diagnostic Plotting Issues in merge_PIPS_realtime.py
**Problem**: Diagnostic plotting code had dataset usage errors and redundant parameter setup
**Solution**:
- Fixed Plot 6 (card DSD meteogram) to use `parsivel_combined_card_ds` instead of incorrectly using merged data
- Created helper functions to eliminate code duplication:
  - `setup_dsd_plot_params()`: Returns diameter bins and log_ND plotting parameters
  - `format_time_axis()`: Handles consistent time axis formatting
- Updated `create_diagnostic_plots()` function interface to accept original card dataset parameter
- Reduced code duplication by ~50 lines while ensuring data source accuracy

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

#### 2. Integrated Merge Script into Analysis Workflow
**Problem**: merge_PIPS_realtime.py existed as standalone script but wasn't integrated into standard analysis pipeline
**Solution**:
- Added merge step as Step 2 in `run_initial_analysis.sh` workflow (between CSV conversion and QC)
- Updated workflow sequence: CSV→netCDF → **merge realtime** → QC → derived params → MM fits
- Added `RUN_MERGE_REALTIME=1` flag and `--skip-merge` command-line option
- Updated help documentation and workflow descriptions
- Added logging support for merge step

**Files Modified**: `shell_scripts/run_initial_analysis.sh`

#### 3. Added Real-time Directory Support
**Problem**: merge script required `--realtime-dir` parameter but workflow script didn't support it
**Solution**:
- Added `--realtime-dir DIR` command-line option to `run_initial_analysis.sh`
- Implemented validation logic: option is required when merge step is enabled
- Added directory existence checking with clear error messages
- Updated help text with usage examples
- Added real-time directory to workflow summary display

**Files Modified**: `shell_scripts/run_initial_analysis.sh`

### Key Lessons Learned

#### Data Source Management in Diagnostic Plotting
**Lesson**: When creating diagnostic plots that compare different data sources (merged vs original), it's critical to use the correct dataset for each plot. Card-only plots must use original card data, not merged data.
**Impact**: Ensures diagnostic plots accurately represent data quality and merge effectiveness.

#### Helper Functions for Code Maintainability
**Lesson**: Repetitive plotting parameter setup should be extracted into helper functions to reduce duplication and improve consistency.
**Impact**: Easier maintenance, consistent appearance across plots, reduced chance of parameter mismatches.

#### Workflow Integration Dependencies
**Lesson**: When integrating new scripts into existing workflows, careful consideration of command-line parameter requirements and validation is essential.
**Impact**: Better user experience, clearer error messages, more robust automation.

#### Required Parameter Validation
**Lesson**: Scripts with required parameters should validate them early with clear error messages and suggested solutions.
**Impact**: Reduces user confusion and provides actionable guidance when parameters are missing.

### Technical Context for Future Development

#### Current Merge Workflow Architecture
- **merge_PIPS_realtime.py**: Handles real-time/card data merging with GPS time correction
- **Positioned optimally**: After CSV conversion but before QC, ensuring merged data goes through quality filters
- **Diagnostic plotting**: Enabled by default with proper dataset separation
- **Configuration-driven**: Uses PIPS_IO_dict with deployment-specific parameters

#### Shell Script Integration Patterns
- **Modular flags**: Each major step has `RUN_*` flag and corresponding `--skip-*` option
- **Logging support**: All steps support optional logging to separate stdout/stderr files
- **Parameter validation**: Required parameters are validated early with helpful error messages
- **Help documentation**: Comprehensive usage examples and workflow descriptions

#### Development Workflow Best Practices Established
1. **Diagnostic plotting**: Always include diagnostic capabilities for data merging/processing steps
2. **Helper functions**: Extract common parameter setup into reusable functions
3. **Dataset accuracy**: Use correct data sources in comparative plots (original vs processed)
4. **Command-line interface**: Provide both required and optional parameters with clear validation
5. **Integration testing**: Consider parameter dependencies when adding scripts to workflows

### Next Development Priorities
- Consider adding configuration validation to merge script
- Evaluate whether other analysis scripts need similar diagnostic plotting capabilities
- Assess whether helper function patterns should be applied to other plotting modules

---
*Last Updated: January 17, 2026*
*Session: Real-time Data Merging Integration*