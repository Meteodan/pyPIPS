# pyPIPS Development Progress Log

This file tracks significant development progress, lessons learned, and agent activities for pyPIPS development sessions.

## Session: January 18, 2026 - Compass QC Enhancement and Deployment Trimming

### Key Code Updates

#### 1. Restructured Data Flow Pattern in apply_QC.py
**Problem**: Complex file loading logic with multiple reads/writes of conventional dataset scattered throughout processing loop
**Solution**:
- Load `conv_ds` once at start of main loop if compass QC or trimming requested
- Process all operations (trimming, compass QC) on in-memory dataset
- Save both conventional and parsivel datasets together at end of loop
- Eliminated redundant file existence checks and intermediate file I/O operations

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 2. Fixed Logging for Time Dimension Operations
**Problem**: Dataset trimming logged total dataset size instead of time dimension length
**Solution**:
- Changed from `len(conv_ds)` to `len(conv_ds.time)` for accurate record counts
- Updated `update_time_attributes()` to show both original and new time ranges
- Format: "Original: START to END" / "New: START to END" for clear comparison
- Logs time coordinate encoding updates (e.g., "seconds since YYYY-MM-DD HH:MM:SS")

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 3. Updated netCDF Time Coordinate Encoding
**Problem**: Trimmed datasets still referenced original start time in netCDF time coordinate units
**Solution**:
- Added time coordinate encoding update in `update_time_attributes()` function
- Sets `dataset['time'].encoding['units']` to "seconds since [new_start_time]"
- Sets `dataset['time'].encoding['calendar']` to 'proleptic_gregorian'
- Ensures netCDF files serialize with correct temporal reference after trimming

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 4. Implemented Circular Statistics for Compass QC
**Problem**: Z-score outlier detection with regular statistics over-cleaned stable compass data (very small std dev → large z-scores for tiny fluctuations)
**Solution**:
- Created `circular_mean_deg()` function using vector averaging (arctan2 of mean sin/cos components)
- Created `circular_std_deg()` function using resultant vector length: `sqrt(-2 * ln(R))`
- Compute circular deviations handling 360°/0° wrap (e.g., 350° - 10° = -20°, not 340°)
- Calculate z-scores using circular std instead of regular std

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 5. Implemented Dual Threshold Compass Outlier Detection
**Problem**: Single z-score threshold flagged legitimate data in stable compass readings as outliers
**Solution**:
- Require BOTH conditions for outlier flagging:
  - Z-score > threshold (default: 3.0, increased from 0.5)
  - AND absolute deviation > threshold (default: 20.0°, adjusted from initial 5.0°)
- Prevents over-cleaning when std is very small (stable data)
- Still catches genuine large deviations (deployment/retrieval periods)
- Added `--compass-abs-threshold` command-line argument

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 6. Reordered Compass QC Before Deployment Trimming
**Problem**: Deployment period detection used noisy compass data, potentially missing or misidentifying unstable periods
**Solution**:
- Moved compass QC section before deployment trimming section
- Cleaned compass data fed into `detect_deployment_periods()` for more reliable circular std dev calculations
- Processing order: Load → Compass QC → Deployment Trimming → Save
- Improved detection accuracy by removing outliers that could mask deployment patterns

**Files Modified**: `analysis_scripts/apply_QC.py`

#### 7. Enhanced Logging and Metadata Storage
**Problem**: Insufficient diagnostics about QC operations and parameter values
**Solution**:
- Log original circular mean, circular std, and QC'd mean
- Log number and percentage of outliers removed
- Store both `qc_zscore_threshold` and `qc_abs_threshold` in dataset attributes
- Store `circular_std` and `n_outliers_removed` in compass_dir attributes
- Apply same attributes to resampled parsivel compass data

**Files Modified**: `analysis_scripts/apply_QC.py`

### Key Lessons Learned

#### Dual Threshold Approach for Stable Compass Data
**Lesson**: Single z-score thresholding fails for very stable data because small std dev makes tiny fluctuations produce large z-scores. Require both statistical significance (z-score) AND practical significance (absolute deviation >20°) to flag true outliers.
**Impact**: Preserves legitimate stable compass readings while still catching genuine outliers. Critical for datasets like PIPS 1A, 1B, 2A, 2B with very stable compass behavior.

#### Circular Statistics Essential for Directional Data
**Lesson**: Regular statistics (mean, std) fail for circular data due to 360°/0° discontinuity. Must use vector averaging for mean (arctan2 of sin/cos components) and resultant vector length for std (sqrt(-2*ln(R))).
**Impact**: Proper handling of compass/wind direction data prevents incorrect outlier flagging and ensures accurate QC. Establishes pattern for all future directional data analysis.

#### Processing Order Matters for Dependent Operations
**Lesson**: Compass QC should precede deployment trimming because trimming relies on circular std dev calculations. Cleaning outliers first improves detection algorithm performance.
**Impact**: Sequential dependency must be considered when organizing QC operations. Cleaner input data leads to more reliable automated detection.

#### Time Coordinate Encoding Requires Explicit Update
**Lesson**: When trimming xarray datasets, the time coordinate encoding (units string) must be explicitly updated to reference the new start time. Otherwise netCDF files contain misleading temporal reference.
**Impact**: Critical for data provenance and correct time interpretation by downstream tools. Must update `dataset['time'].encoding['units']` whenever time range changes.

### Technical Context for Future Development

#### Circular Statistics Functions
- `circular_mean_deg(angles_deg)`: Computes vector-averaged mean for 0-360° data, handles NaN values
- `circular_std_deg(angles_deg)`: Uses resultant vector length R = sqrt(C² + S²), returns std in degrees
- Both functions handle wrap-around: ensure results in 0-360° range, properly compute deviations across 360°/0° boundary

#### Dual Threshold QC Pattern
```python
is_outlier = (np.abs(zscores) > zscore_threshold) & (abs_deviations > abs_threshold)
```
- Prevents false positives in stable data while maintaining sensitivity to real outliers
- Tunable via command-line: `--compass-zscore-threshold` (3.0) and `--compass-abs-threshold` (20.0)
- Apply same pattern to parsivel resampled data for consistency

#### Data Flow Pattern in apply_QC.py
1. Load conventional dataset once at loop start if needed
2. Apply compass QC (outlier removal, wind recomputation)
3. Apply deployment trimming (detect periods, trim both datasets, update time attributes)
4. Save both conventional and parsivel datasets at loop end
- In-memory processing avoids repeated file I/O
- Maintains synchronization between conventional and parsivel datasets

#### Time Encoding Update Pattern
```python
time_units = f"seconds since {start_time.strftime('%Y-%m-%d %H:%M:%S')}"
dataset['time'].encoding['units'] = time_units
dataset['time'].encoding['calendar'] = 'proleptic_gregorian'
```
- Apply whenever time range changes (trimming, subsetting)
- Ensures netCDF temporal reference matches actual data extent

### Next Development Priorities
- Test dual threshold parameters across diverse field campaigns to validate 20° threshold
- Consider adaptive thresholding based on deployment history (stable vs mobile platforms)
- Evaluate if MAD-based outliers (Median Absolute Deviation) would provide additional robustness
- Add compass QC diagnostic plots to show circular mean/std and outlier distribution
- Document circular statistics functions for use with other directional variables (wind direction)

---

## Session: January 17, 2026 - Diagnostic Plotting Infrastructure and File Handling

### Key Code Updates

#### 1. Created Comprehensive Diagnostic Plotting Script
**Problem**: Need systematic visualization of diagnostic variables from PIPS netCDF files
**Solution**:
- Created `plot_PIPS_diag.py` with support for both conventional and parsivel diagnostic variables
- Implemented conventional diagnostics: GPS variables (lat/lon/alt/speed), battery voltage, compass direction, wind diagnostic
- Implemented parsivel diagnostics: particle counts (log scale), reflectivity, rain rate, sensor temp/voltage, signal amplitude, accumulated precipitation, sample interval
- Added command-line flags for selective plotting (`--plot-conv`, `--plot-parsivel`)
- Organized output into separate subdirectories for conventional vs parsivel data

**Files Modified**: `plotting_scripts/plot_PIPS_diag.py`

#### 2. Refactored Plotting Functions to plotmodule.py
**Problem**: Need to centralize plotting logic and follow established patterns
**Solution**:
- Moved 9 diagnostic plotting functions from script to `plotmodule.py` module
- Implemented proper `axparams` dictionary pattern matching existing meteogram functions
- Added `set_meteogram_axes()` integration for consistent axis formatting
- Eliminated redundant wrapper functions by calling existing plotmodule functions directly
- Functions use `global_plot_config_dict` for locators, formatters, and labels

**Files Modified**: `pyPIPS/plotmodule.py`, `plotting_scripts/plot_PIPS_diag.py`

#### 3. Fixed plotmeteogram Array Passing Pattern
**Problem**: TypeError due to incorrect data structure - single arrays passed instead of lists
**Solution**:
- Corrected all plotting functions to pass `[plottimes]` (list) instead of `plottimes` (array)
- Changed data passing to `fields = [data.to_numpy()]` pattern
- Follows established pattern from `plot_voltage_meteogram` and other existing functions
- Resolved "vertices must be 2D with shape (N, 2)" error

**Files Modified**: `plotting_scripts/plot_PIPS_diag.py`, `pyPIPS/plotmodule.py`

#### 4. Fixed plot_date() Format Argument Conflicts
**Problem**: TypeError - "got multiple values for argument 'fmt'" in matplotlib plot_date calls
**Solution**:
- Converted positional format strings (e.g., `'b-'`) to explicit keyword arguments (`ls='-', color='b'`)
- Maintained `fmt=""` keyword for all plot_date calls
- Applied fix to 7 different plot_date invocations across GPS, temperature, voltage, and other plots
- Consistent with matplotlib API requirements

**Files Modified**: `pyPIPS/plotmodule.py`

#### 5. Resolved File Locking Issues in merge_PIPS_realtime.py
**Problem**: PermissionError when saving datasets - files still open in memory when attempting to overwrite
**Solution**:
- Added `.load()` and `.close()` calls after opening both conventional and parsivel datasets
- Loads data into memory and releases file handles before saving operations
- Implemented try-except blocks with specific `PermissionError` handling
- Uses `utils.fatal()` to exit script with descriptive error messages on save failures
- Handles both output files (conventional and parsivel combined) with proper error reporting

**Files Modified**: `analysis_scripts/merge_PIPS_realtime.py`

### Key Lessons Learned

#### Plotting Function Organization Patterns
**Lesson**: Diagnostic plotting functions should be centralized in `plotmodule.py` rather than duplicated across scripts, with consistent use of axparams dictionaries and `set_meteogram_axes()`.
**Impact**: Enables code reuse, maintains consistency across all plotting scripts, and simplifies future modifications to plotting behavior.

#### Data Structure Requirements for plotmeteogram
**Lesson**: The `plotmeteogram()` function expects lists of arrays for both `xvals` and `zvals` parameters, even when plotting single variables. This enables multi-variable overlays on same plot.
**Impact**: Critical pattern to follow when creating new plotting functions - always wrap time and data arrays in lists: `[plottimes]`, `[data.to_numpy()]`.

#### matplotlib plot_date() API Constraints
**Lesson**: Cannot pass format strings both positionally and via `fmt` keyword argument. Must use either positional format OR explicit keyword arguments (ls, color, marker) with `fmt=""`.
**Impact**: Establishes clear pattern for datetime plotting that avoids API conflicts while maintaining backward compatibility.

#### xarray File Handle Management
**Lesson**: `xr.open_dataset()` keeps file handles open by default. When overwriting input files, must explicitly `.load()` into memory and `.close()` handles before saving. Always wrap `.to_netcdf()` in try-except with PermissionError handling.
**Impact**: Prevents file locking issues and provides clear error messages when files cannot be saved. Essential for any script that might overwrite its input files.

### Technical Context for Future Development

#### Diagnostic Plotting Architecture
- New functions in plotmodule.py: `plot_GPS_variables`, `plot_winddiag_meteogram`, `plot_parsivel_counts_meteogram`, `plot_parsivel_reflectivity_meteogram`, `plot_parsivel_rainrate_meteogram`, `plot_parsivel_temp_voltage_meteogram`, `plot_parsivel_signal_amplitude_meteogram`, `plot_parsivel_accumulation_meteogram`, `plot_parsivel_sample_interval_meteogram`
- All follow pattern: accept `global_plot_config_dict`, create `axparamdict` with locators/formatters/limits/labels, call `set_meteogram_axes()`
- Variable names from pips_io.py: `conv_df_to_ds()` and `parsivel_df_to_ds()` provide authoritative variable naming

#### File Handling Best Practices
- Load datasets into memory with `.load()` before closing file handles
- Always close datasets with `.close()` after loading
- Wrap all `.to_netcdf()` calls in try-except blocks
- Catch `PermissionError` separately from general exceptions
- Use `utils.fatal()` for script termination with error messages

### Next Development Priorities
- Consider adding flagged_times visualization (quality control flags)
- Implement missing_times diagnostic plotting in plot_PIPS_diag.py
- Add optional derived variable overlays (e.g., computed vs observed values)
- Create summary statistics output for diagnostic variables

---

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