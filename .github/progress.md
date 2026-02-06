# pyPIPS Development Progress Log

This file tracks significant development progress, lessons learned, and agent activities for pyPIPS development sessions.

## Session: February 4-5, 2026 - Interactive Notebook Optimization for VS Code

### Key Code Updates

#### Matplotlib Widget Integration and Figure Size Optimization
**Problem**: Manual QC inspector notebook's interactive matplotlib widgets worked in Jupyter browser but had issues in VS Code - initially widgets didn't activate at all, and once activated, horizontal scrolling wasn't supported in output cells, making it difficult to view full time series plots.

**Solution**: Enhanced VS Code compatibility for interactive workflow:
1. **Figure Size Adjustment**: Changed default from `figsize=(14, 8)` to `figsize=(11, 7)` to fit VS Code output cells without horizontal scrolling
2. **Configurable Size Parameter**: Added `figsize` parameter to both `inspect_variable()` and `next_dataset()` functions for flexibility:
   - Default `(11, 7)` for VS Code viewport
   - Can override with `(14, 8)` for browser or large screens
   - Can reduce to `(9, 6)` for compact displays
3. **Figure Cleanup**: Added `plt.close('all')` after imports in `inspect_variable()` to prevent figure accumulation across multiple inspections

**Workflow Impact**: Enables complete interactive manual QC workflow within VS Code:
- Interactive time selection with mouse clicks (hover to see time, left click for start/end, right click to clear)
- Zoom and pan tools work natively
- No need to switch between browser and VS Code anymore
- Eliminates previous workflow of: inspect in browser → edit JSON in VS Code → reload notebook → verify

**Files Modified**: `notebooks/manual_QC_inspector.ipynb` (Cell 6: inspect_variable function, Cell 7: next_dataset function)

#### Notebook Editing Tool Selection and Sync Issues
**Problem**: When editing notebook files, using the `edit_notebook_file` tool sometimes failed to persist changes when notebook was open in both VS Code and browser. Tool reported success but changes weren't visible on reload.

**Root Cause**: Notebook JSON editing conflicts when file is actively loaded in multiple editors. VS Code and Jupyter may cache notebook state differently, causing synchronization issues.

**Solution Pattern**: Switched to direct JSON editing via `replace_string_in_file` for notebook modifications:
- More reliable when notebook is open in multiple environments
- Requires manual save-reload-verify cycle but ensures changes persist
- Trade-off: Less validation (must ensure valid JSON structure), but necessary for multi-editor workflow
- Best practice: Close browser notebook when making structural edits, or use save-reload cycle

**Technical Context**: VS Code notebook API (`edit_notebook_file`) is preferred for structural operations (cell insertion/deletion, metadata changes) when notebook is exclusively in VS Code. For simple code additions in dual-editor scenarios, direct file editing with `replace_string_in_file` is more reliable.

**Impact**: Future notebook edits should consider:
- Single editor (VS Code only): Use `edit_notebook_file` for full validation and structural awareness
- Multiple editors (VS Code + browser): Use `replace_string_in_file` with care, verify JSON syntax
- Alternative: Establish clear editing boundaries (browser for execution, VS Code for editing, explicit save/reload)

### Key Lessons Learned

#### VS Code vs Browser for Interactive Notebooks
**Lesson**: VS Code's native notebook support has matured enough to handle matplotlib interactive widgets (`%matplotlib widget`), but with viewport constraints requiring adaptation.

**Context**:
- Previously assumed interactive widgets required browser-based Jupyter
- VS Code now supports ipympl/widget backend, enabling pan/zoom/click interactions
- Key limitation: Output cells have fixed width, no horizontal scrolling support
- Workaround: Adjust figure dimensions to fit viewport rather than relying on scrolling

**Best Practice**:
```python
# Design functions with configurable figure size
def plot_function(..., figsize=(11, 7)):  # Default for VS Code
    fig, axes = plt.subplots(..., figsize=figsize)
```

**Impact**: Enables full-featured interactive analysis within VS Code IDE, eliminating context switching. Particularly valuable for workflows combining code editing, execution, and interactive data exploration.

#### Notebook Tool Selection Strategy
**Lesson**: Different notebook editing tools have different reliability characteristics depending on workflow context.

**Tool Comparison**:
- **`edit_notebook_file`** (VS Code notebook API):
  - ✓ Understands notebook cell structure and IDs
  - ✓ Validates changes against notebook schema
  - ✓ Handles structural operations (cell insertion, deletion, reordering)
  - ✗ May have cache/sync issues with multi-editor workflows
  - ✗ Requires notebook to be in clean state (saved, not executing)

- **`replace_string_in_file`** (direct JSON editing):
  - ✓ More reliable for dual-editor scenarios (VS Code + browser)
  - ✓ Direct file manipulation, no caching layer
  - ✓ Works even when notebook is executing or has unsaved changes
  - ✗ No validation of notebook JSON structure
  - ✗ Fragile for structural changes (must manually handle cell IDs, metadata)
  - ✗ Requires understanding notebook JSON format

**Decision Framework**:
- Use `edit_notebook_file` when: Single editor, structural changes, want validation
- Use `replace_string_in_file` when: Multi-editor workflow, simple code additions, experiencing sync issues

**Impact**: Informs tool selection for future notebook modifications. Documents tradeoffs for agent sessions and human developers.

#### Figure Lifecycle Management in Interactive Notebooks
**Lesson**: Matplotlib figures persist across cell executions in notebooks, requiring explicit cleanup to avoid memory bloat and visual clutter.

**Context**: Repeated calls to plotting functions create new figure windows without closing previous ones. In interactive backends (`%matplotlib widget`), this accumulates figures in notebook output cells.

**Best Practice**:
```python
import matplotlib.pyplot as plt
plt.close('all')  # Close existing figures before creating new ones

def plot_function():
    plt.close('all')  # Or close specific figures with plt.close(fig)
    fig, ax = plt.subplots()
    ...
```

**Impact**: Applies to any interactive notebook workflow with repeated plotting. Prevents notebook file size bloat and UI clutter from accumulated figure objects.

### Technical Context for Future Development

**Interactive Notebook Design Pattern**:
When building interactive notebooks for VS Code:
1. Use `%matplotlib widget` for interactivity
2. Design figures to fit ~11-12 inch width to avoid horizontal scrolling
3. Add `figsize` parameters for user customization
4. Include `plt.close('all')` at start of plotting functions
5. Test in VS Code environment, not just browser

**Manual QC Inspector Architecture**:
Current workflow now fully VS Code-compatible:
1. Load config with `PIPS_IO_dict` → identifies datasets to inspect
2. Iterator function `next_dataset()` cycles through PIPS/IOP combinations
3. Interactive selection: matplotlib event handlers capture mouse clicks → store timestamps globally
4. Decision recording: `add_manual_qc()` saves to JSON with ISO format timestamps
5. Batch application: `apply_manual_qc.py` reads JSON → modifies NetCDF → preserves provenance

**Notebook vs Script Division**:
- **Notebook**: Interactive exploration, visual identification of issues, decision recording
- **Script**: Batch processing, automated application of recorded decisions, reproducible workflow
- **Separation rationale**: Notebooks for human judgment, scripts for automation and version control

### Next Development Priorities

1. **Interactive Selection Enhancement**: Consider adding visual feedback for selected time ranges (shaded regions on plot) to confirm selection before saving to JSON

2. **Multi-variable QC Support**: Extend inspector to allow flagging multiple variables simultaneously when same time period affects multiple sensors

3. **QC Decision Validation**: Add diff/comparison view showing before/after when applying manual QC to verify correct time ranges being modified

4. **Progress Persistence**: Save iterator state to JSON so inspection session can be resumed after interruption

5. **Keyboard Shortcuts**: Add key bindings for common operations (next dataset, save QC, undo last entry) to speed up inspection workflow

## Session: February 3-4, 2026 - Manual QC Data Assignment Fix and xarray Encoding Issues

### Key Code Updates

#### Manual QC Variable Modification: xarray DataArray Assignment
**Problem**: The `apply_manual_qc_to_dataset()` function in `apply_manual_qc.py` was not successfully setting flagged time ranges to NaN. Multiple approaches tried (`.where()`, direct `.values` assignment, `.data` assignment, `.loc` indexer) all appeared to work (attributes were updated, debug showed NaN values in arrays), but the saved NetCDF files did not contain the expected NaN values. The issue was eventually traced to problems in the time selection logic in the notebook (unrelated to the script itself), combined with xarray's NetCDF encoding behavior.

**Solution**: Implemented robust variable modification approach:
1. **Time Selection**: Convert dataset times explicitly to timezone-naive pandas datetime for comparison
   ```python
   ds_times = pd.to_datetime(ds.time.values)
   if hasattr(ds_times, 'tz') and ds_times.tz is not None:
       ds_times = ds_times.tz_localize(None)
   time_mask = (ds_times >= start_time) & (ds_times <= end_time)
   ```

2. **Data Modification**: Use numpy indexing on copied data
   ```python
   time_indices = np.where(time_mask)[0]
   current_data = ds[variable].values.copy()
   current_data[time_indices] = np.nan
   ```

3. **DataArray Reconstruction**: Rebuild DataArray from scratch with cleared encoding
   ```python
   new_var = xr.DataArray(
       current_data,
       coords=ds[variable].coords,
       dims=ds[variable].dims,
       attrs=ds[variable].attrs.copy()
   )
   new_var.encoding = {}  # Clear to avoid fill value conflicts
   ds[variable] = new_var
   ```

**Key Insight**: NetCDF encoding attributes (`_FillValue`, `missing_value`, dtype specifications) can interfere with NaN handling during `to_netcdf()`. Clearing the encoding dictionary before assignment ensures xarray uses default NaN handling for modified variables.

**Files Modified**: `analysis_scripts/apply_manual_qc.py` (apply_manual_qc_to_dataset function)

**Debugging Process**: Extensive debugging revealed that the code logic was correct throughout - the actual issue was in how time ranges were selected in the notebook's interactive interface (matplotlib date handling and click coordinate conversion), not in the QC application script itself. Once the notebook issue was resolved, the script worked as designed.

### Key Lessons Learned

#### xarray DataArray Assignment Patterns
**Lesson**: When modifying xarray DataArrays that will be saved to NetCDF, rebuilding the DataArray from scratch with cleared encoding is more reliable than in-place modification.

**Context**: xarray's lazy evaluation and NetCDF encoding preservation can cause unexpected behavior when modifying variable values. Direct assignment to `.values` or `.data` may work in memory but fail to persist correctly to disk if encoding attributes conflict with the modified data type (e.g., NaN vs integer fill values).

**Best Practice**:
```python
# Copy data, modify, rebuild DataArray
modified_data = da.values.copy()
modified_data[mask] = np.nan
new_da = xr.DataArray(modified_data, coords=da.coords, dims=da.dims, attrs=da.attrs.copy())
new_da.encoding = {}  # Clear encoding
ds[var_name] = new_da
```

**Impact**: Applies to any workflow where xarray variables are modified after loading from NetCDF and must be saved back to disk. Particularly important when setting values to NaN in variables that originally had non-NaN fill values.

#### Debugging Complex Data Pipelines
**Lesson**: When data appears correct at every step but final output is wrong, the issue may be in an earlier stage (data selection/preparation) rather than the processing logic being debugged.

**Context**: Spent significant time debugging the QC application logic with various assignment approaches, all of which appeared to work when tested. The actual problem was that incorrect time ranges were being captured in the notebook's interactive selection, so the "wrong" times were being flagged (correctly according to the JSON, but not according to user intent).

**Debugging Strategy**:
1. Verify input data (time ranges in JSON) matches expectations
2. Test with known-good inputs to isolate whether issue is in current step or upstream
3. Add comprehensive debug output at each stage (inputs, intermediate results, outputs)
4. Don't assume earlier stages are correct just because they "worked before"

**Impact**: Reinforces importance of validating assumptions about input data and testing workflows end-to-end with known test cases.

### Technical Context for Future Development

**Manual QC Variable Modification Pattern**:
The final working implementation in `apply_manual_qc_to_dataset()`:
1. Extract and timezone-normalize dataset times for comparison
2. Create boolean mask with pandas datetime comparison
3. Find array indices with `np.where(time_mask)[0]`
4. Copy variable data, modify with numpy indexing
5. Rebuild DataArray with cleared encoding, replace in dataset
6. Update provenance attributes (manual_qc_applied, has_manual_qc)

**xarray Encoding Management**:
- Encoding dictionaries control how variables are written to NetCDF
- Include: dtype, fill value, compression, chunking, scale/offset
- Can cause issues when data type changes (e.g., introducing NaN to integer data)
- Clearing encoding (`new_var.encoding = {}`) lets xarray auto-detect appropriate settings
- Alternative: Explicitly set encoding with NaN-compatible fill value

**Interactive Notebook Time Selection**:
User resolved issues with matplotlib coordinate conversion and date handling in the inspection notebook. Details of that fix would be in the notebook's development history. The script assumes time ranges in JSON are correct and processes them faithfully.

### Next Development Priorities

1. **End-to-End Manual QC Testing**: Document complete workflow from notebook inspection through script application with test cases to verify correct time range capture and application

2. **Encoding Best Practices Documentation**: Create guide for when to preserve vs clear xarray encoding when modifying variables, especially for NaN handling scenarios

3. **Interactive Selection Improvements**: Consider adding visual confirmation in notebook when time ranges are selected (highlight selected regions, show exact timestamps before saving to JSON)

4. **Batch QC Validation**: Tool to compare before/after statistics for manual QC applications to verify expected number of points flagged per variable/time range

5. **Config Validation Enhancement**: Extend validation to check that filename lists match in length and that files exist before processing begins

## Session: February 3, 2026 - Manual QC Workflow Fixes and Timezone/NetCDF Compatibility

### Key Code Updates

#### 1. Manual QC Notebook: Date Axis Display and Dataset Iterator Fixes
**Problem**: Interactive manual QC inspector notebook displayed incorrect time ranges (showing 1970 epoch instead of actual 2025 data dates). Dataset iterator logic incorrectly treated `deployment_names` as a list of unique IOPs, creating nested loops that produced 16 iterations instead of expected 4 (for 1 IOP with 4 PIPS). Click event handlers initially broke during matplotlib date conversion fixes.

**Solution**: Fixed matplotlib date display and corrected fundamental data structure understanding:
- **Date Display Fix**:
  - Converted time coordinate explicitly: `time_pd = pd.to_datetime(ds.time.values)`
  - Used `mdates.date2num()` for plotting: `time_mpl = mdates.date2num(time_pd)`
  - Set explicit x-axis limits: `ax1.set_xlim(time_mpl[0], time_mpl[-1])` to prevent default epoch display
  - Preserved click handlers by using `mdates.num2date(event.xdata)` for coordinate conversion
  - Result: Plots now show correct 2025 date ranges with working interactive time selection

- **Iterator Structure Fix**:
  - Corrected understanding: `deployment_names` is **parallel** to `PIPS_names` (1:1 mapping), not a list of unique IOPs
  - Changed from nested loops `for iop in deployment_names: for pips in PIPS_names:`
  - To flat iterator: `dataset_list = list(zip(PIPS_names, deployment_names))`
  - Added unique IOP detection: `unique_iops = list(dict.fromkeys(deployment_names))` for progress display
  - Single loop through dataset_list with proper progress tracking (IOP number, PIPS within IOP, overall count)
  - Result: Correct iteration count (4 datasets for 1 IOP with 4 PIPS)

**Files Modified**: `notebooks/manual_QC_inspector.ipynb` (Cell 6: inspect_variable function, Cell 7: next_dataset iterator)

#### 2. Apply Manual QC Script: Data Structure and Filename Handling
**Problem**: Production script `apply_manual_qc.py` had incorrect nested loop structure (same misunderstanding as notebook). Used inefficient pattern-based filename searching (`parsivel_combined_{deployment}_{pips}_*.nc`) instead of utilizing existing config file structure with explicit filenames.

**Solution**: Aligned script with corrected data structure and config-based file access:
- **Loop Structure**:
  - Changed to `for idx, (pips_name, deployment_name) in enumerate(zip(PIPS_names, deployment_names)):`
  - Added index for accessing parallel filename lists
  - Removed nested loops entirely

- **Filename Access**:
  - Load filename lists from config: `conv_filenames_nc = config.PIPS_IO_dict.get('conv_filenames_nc', [])`
  - Direct lookup: `conv_file = os.path.join(PIPS_dir, conv_filenames_nc[idx])`
  - Removed all pattern matching logic (parsivel_patterns lists, os.path.exists() searching loops)
  - Removed filename construction from deployment/PIPS names
  - Result: More reliable, matches config structure exactly, clearer intent

**Files Modified**: `analysis_scripts/apply_manual_qc.py` (loop structure, config loading section, file access logic)

#### 3. Timezone-Naive Timestamp Handling
**Problem**: Script crashed with "Cannot compare tz-naive and tz-aware timestamps" when applying QC. JSON-parsed timestamps or dataset times had timezone information that didn't match, preventing pandas comparison operations.

**Solution**: Strip timezone information from all timestamps before comparison:
- **In `apply_manual_qc_to_dataset()`**: Convert start_time and end_time from JSON to timezone-naive using `tz_localize(None)` after `pd.to_datetime()`
- **In `align_times_with_parsivel()`**: Strip timezone from parsivel dataset times: `parsivel_times.tz_localize(None)`
- **In `apply_manual_trim_to_dataset()`**: Strip timezone from original dataset start/end times and parsed trim times
- **In main loop**: Strip timezone from requested_start/requested_end before alignment

**Pattern**: Check `if obj.tz is not None:` then apply `obj.tz_localize(None)` for consistency across all time comparisons

**Files Modified**: `analysis_scripts/apply_manual_qc.py` (5 locations handling time parsing/comparison)

#### 4. NetCDF Boolean Attribute Compatibility
**Problem**: Script failed when saving with "illegal data type for attribute 'has_manual_qc', got b1". NetCDF doesn't support Python boolean types directly - requires specific numeric types from `['S1', 'i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8', 'f4', 'f8']`.

**Solution**: Changed boolean flag to integer:
- `ds[variable].attrs['has_manual_qc'] = True` → `ds[variable].attrs['has_manual_qc'] = 1`
- Comment updated to note NetCDF compatibility requirement
- Result: Files save successfully, flag still provides easy QC status checking (1 = has manual QC, 0 or missing = no manual QC)

**Files Modified**: `analysis_scripts/apply_manual_qc.py` (apply_manual_qc_to_dataset function)

#### 5. Slowtemp QC: Dewpoint/RH Resampling to Parsivel Dataset
**Problem**: When slowtemp QC flagged bad data points, `recompute_dewpoint_RH_with_fallback()` updated dewpoint and RH_derived in the conventional (1-Hz) dataset. However, these recomputed values were not resampled to the parsivel dataset (10-second resolution), creating inconsistency between conventional and parsivel files.

**Solution**: Added dewpoint and RH_derived resampling after slowtemp QC:
- After resampling cleaned slowtemp to parsivel dataset, check if `n_flagged > 0` (QC was applied)
- If dewpoint exists in parsivel dataset, resample from conventional: `conv_ds['dewpoint'].resample(...).mean()`
- If RH_derived exists in parsivel dataset, resample from conventional: `conv_ds['RH_derived'].resample(...).mean()`
- Uses same resampling parameters as slowtemp: interval_str, label='right', closed='right', offset_str
- Result: Parsivel and conventional datasets remain thermodynamically consistent after slowtemp QC

**Files Modified**: `analysis_scripts/apply_QC.py` (slowtemp QC section, lines ~1010-1030)

### Key Lessons Learned

#### Data Structure Pattern: Parallel Lists vs Nested Structures
**Lesson**: Config file's `deployment_names` is a **label array** parallel to `PIPS_names`, not a list of unique deployment identifiers for nested iteration.

**Context**:
- Correct: `zip(PIPS_names, deployment_names)` where each PIPS has its deployment label
- Example: `PIPS_names = ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS3A']`, `deployment_names = ['IOP1', 'IOP1', 'IOP1', 'IOP1']`
- Incorrect: Iterating `for iop in deployment_names:` as if it's `['IOP1', 'IOP2', 'IOP3']` with nested `for pips in PIPS_names:`

**Impact**: Critical for understanding config structure. Applies to all scripts using these lists (analysis, plotting, processing). Config-driven workflows rely on parallel indexing: `PIPS_names[i]`, `deployment_names[i]`, `conv_filenames_nc[i]`, `PIPS_filenames_nc[i]`.

#### Timezone Handling in Scientific Data Pipelines
**Lesson**: Always normalize timestamps to timezone-naive before comparisons in scientific workflows. Mixed timezone awareness causes comparison failures in pandas/numpy operations.

**Pattern**:
```python
timestamp = pd.to_datetime(value)
if hasattr(timestamp, 'tz') and timestamp.tz is not None:
    timestamp = timestamp.tz_localize(None)
```

**Context**: JSON serialization, xarray coordinate parsing, and pandas datetime operations may produce timezone-aware objects. NetCDF datasets typically use timezone-naive timestamps referenced to a base time (e.g., "seconds since YYYY-MM-DD HH:MM:SS").

**Impact**: Apply this pattern whenever parsing timestamps from external sources (JSON, user input) or comparing against dataset time coordinates.

#### NetCDF Attribute Type Restrictions
**Lesson**: NetCDF attributes must use specific numeric types - booleans must be stored as integers (0/1), not Python bool type.

**Allowed types**: `['S1', 'i1', 'u1', 'i2', 'u2', 'i4', 'u4', 'i8', 'u8', 'f4', 'f8']` (strings, signed/unsigned integers, floats)

**Pattern**: Use `1` for True, `0` for False when storing boolean flags in NetCDF attributes.

**Impact**: Affects any code setting dataset or variable attributes that will be saved to NetCDF. Lists are also problematic - convert to strings or separate attributes when needed.

#### Config-Based File Access vs Pattern Matching
**Lesson**: When config files already specify exact filenames in parallel lists, use direct index-based lookup rather than constructing patterns and searching with `os.path.exists()`.

**Benefit**:
- More reliable (no ambiguity about which file to use if multiple matches)
- Clearer intent (explicit config specification)
- Faster (no filesystem searching)
- Matches data structure (parallel arrays indexed together)

**Pattern**: Load `conv_filenames_nc` and `PIPS_filenames_nc` from config, access via `enumerate()` index in processing loop.

#### Thermodynamic Consistency Across Resampled Datasets
**Lesson**: When recomputing derived thermodynamic variables (dewpoint, RH) in high-resolution datasets, **must resample** to lower-resolution datasets to maintain consistency.

**Context**: Conventional dataset (1-Hz) is primary source, parsivel dataset (10-second) contains resampled met variables. If QC invalidates temperature measurements requiring dewpoint/RH recalculation, both datasets need the updated values.

**Pattern**: After any thermodynamic recomputation at high resolution:
1. Check if variable exists in lower-resolution dataset
2. Resample updated high-resolution variable using same parameters as original processing
3. Update lower-resolution dataset variable

**Impact**: Prevents situations where conventional and parsivel datasets have mismatched thermodynamic states, ensuring analyses using either dataset produce consistent results.

### Technical Context for Future Development

**Manual QC Workflow Architecture**:
1. **Inspection Phase** (Jupyter notebook, browser-based):
   - Load config → initialize tracking → iterate with `next_dataset()`
   - Interactive matplotlib with mdates for proper datetime display
   - Hover/click to select time ranges → `add_manual_qc()` with printed commands
   - Incremental saves with `save_manual_qc()`, management functions (view, delete)
   - User must restart kernel before testing cell changes (browser sync issue)

2. **Application Phase** (command-line script):
   - Read JSON decisions → load config → iterate via `zip(PIPS_names, deployment_names)`
   - Use `enumerate()` index to access parallel filename lists from config
   - Apply trimming first (if specified), then variable QC
   - Handle timezone conversion consistently across all time comparisons
   - Use integer attributes for NetCDF compatibility (not booleans)
   - Add provenance attributes: `manual_qc_applied` (list), `has_manual_qc` (1/0)

**Key File Relationships**:
- Config files: Contain parallel lists `PIPS_names`, `deployment_names`, `conv_filenames_nc`, `PIPS_filenames_nc`
- Manual QC JSON: Keys are `f"{pips_name}_{deployment_name}"`, values are lists of QC entries
- Special key `_manual_trim` in JSON for trimming operations (separate from variable QC)
- Dataset keys and config iteration must use identical `(pips_name, deployment_name)` pairs

**Iterator Pattern for All Scripts**:
```python
for idx, (pips_name, deployment_name) in enumerate(zip(PIPS_names, deployment_names)):
    key = f"{pips_name}_{deployment_name}"
    conv_file = os.path.join(PIPS_dir, conv_filenames_nc[idx])
    parsivel_file = os.path.join(PIPS_dir, PIPS_filenames_nc[idx])
    # Process files...
```

### Next Development Priorities

1. **Automated QC Report Generation**: Create summary reports showing QC statistics (points flagged, variables affected) across all deployments/IOPs for campaign-level quality assessment

2. **Manual QC JSON Merging**: Tool to combine multiple manual QC JSON files (e.g., from different analysts or incremental inspection sessions) with conflict detection/resolution

3. **Bias Correction Validation Plots**: Automated plotting comparing slowtemp vs fasttemp before/after bias correction to visually verify correction effectiveness

4. **Thermodynamic Variable Recomputation**: Extract dewpoint/RH recalculation into standalone utility for use outside QC workflows (e.g., after merging datasets, correcting pressure readings)

5. **Config File Validation**: Tool to verify config file structure (parallel list lengths match, files exist, required keys present) before processing to catch errors early

## Session: February 2, 2026 - Comprehensive QC Framework with Bias Correction and Manual QC Tools

### Key Code Updates

#### 1. Slow Temperature Bias Correction System
**Problem**: PIPS slow temperature sensors exhibited systematic bias relative to fast temperature sensors across ICECHIP campaign. Initial QC notebook analysis (filtering rapid changes and large differences) revealed consistent slope/intercept relationships between slowtemp and fasttemp that varied by PIPS unit. No automated workflow existed to apply pre-computed bias corrections to operational datasets while maintaining proper thermodynamic consistency (dewpoint, RH_derived recomputation).

**Solution**: Implemented comprehensive bias correction system in `apply_QC.py`:
- Added `pyPIPS.thermolib` import for thermodynamic calculations (`calTdfromRH`, `calRH`)
- Created `apply_slowtemp_bias_correction()` function that:
  - Inverts linear relationship: `slowtemp_corrected = (slowtemp - intercept) / slope`
  - **Overwrites** existing `slowtemp` variable (no `_corrected` suffix in production)
  - Adds `bias_correction_slope`, `bias_correction_intercept`, `bias_corrected` attributes
  - Recomputes `dewpoint` using corrected slowtemp with fasttemp fallback when slowtemp is NaN
  - Recomputes `RH_derived` from fasttemp and new dewpoint
  - Updates both conventional (1-Hz) and parsivel (resampled) datasets
- Command-line arguments: `--slowtemp-bias-correction`, `--slowtemp-bias-config <path/to/config.py>`
- Reads slope/intercept from config file: `configs/ICECHIP_slowtemp_bias.py` with `slowtemp_bias` dict keyed by PIPS name
- Only applies correction if `slope != 1.0` or `intercept != 0.0` (allows neutral placeholders)
- Positioned AFTER trimming in processing pipeline (QC → trim → bias correct)
- Reports bias before/after correction for validation

**Files Modified**:
- `analysis_scripts/apply_QC.py` (added function, args, processing section)
- Config example: `configs/ICECHIP_slowtemp_bias.py`

**Workflow**: First pass with `--slowtemp-qc --trim-deployment-periods` → Analyze in notebook → Compute slopes/intercepts → Second pass with `--slowtemp-bias-correction --slowtemp-bias-config`

#### 2. Manual QC Inspection and Application Framework
**Problem**: Automated QC (threshold-based filtering, compass outlier detection, deployment trimming) catches systematic issues but misses localized anomalies requiring human judgment: sensor malfunctions mid-deployment, debris contamination, brief electronic glitches, rapid environmental changes. No interactive workflow existed for visually identifying and documenting time ranges needing manual QC, nor for applying those decisions reproducibly to datasets.

**Solution**: Created two-tool manual QC system:

**A. Interactive Notebook (`notebooks/manual_QC_inspector.ipynb`)**:
- Interactive matplotlib (`%matplotlib widget`) for zoom/pan identification of exact time ranges
- Sequential iterator (`next_dataset()`) loads each PIPS/IOP combination automatically
- Dual-panel plotting: temperatures + differences with threshold reference lines (±2°C)
- `add_manual_qc(pips, iop, variable, start, end, reason)` function records decisions to JSON
- `save_manual_qc()` for incremental saves (avoid losing work)
- `view_manual_qc()`, `delete_last_qc()` helper functions for management
- JSON structure: `{"PIPS1A_IOP1_2025": [{"variable": "slowtemp", "start_time": "...", "end_time": "...", "reason": "...", "flagged_on": "..."}]}`
- Progress tracking through deployments with status messages
- Configurable via `CONFIG_FILE` variable (import any case config)

**B. Application Script (`analysis_scripts/apply_manual_qc.py`)**:
- Command-line tool: `python apply_manual_qc.py <config.py> <manual_qc.json> [--output-tag _manual_qc]`
- Loads QC decisions from JSON, iterates through deployments/PIPS
- Sets data to NaN in specified time ranges for specified variables
- Updates both conventional and parsivel datasets
- Adds provenance attributes: `manual_qc_applied` (list of time ranges + reasons), `has_manual_qc` (bool flag)
- Flexible output: new files with tag (default `_manual_qc`) or overwrite with empty string
- Options: `--conventional-only`, `--parsivel-only` for selective processing
- Statistics reporting: datasets processed, total points flagged
- Handles missing files and variables gracefully with warnings

**Files Created**:
- `notebooks/manual_QC_inspector.ipynb` (10 cells, ~400 lines)
- `analysis_scripts/apply_manual_qc.py` (~350 lines)

**Workflow**: Inspect in notebook → Mark problematic ranges → Save JSON → Apply with script → Use QC'd files in downstream analysis

#### 3. Enhanced QC Processing Order and Documentation
**Problem**: Bias correction was initially placed after slowtemp QC but before trimming. This meant bias coefficients would be computed from data including deployment/retrieval periods (contaminated with handling artifacts). Optimal processing order unclear from code structure.

**Solution**: User repositioned bias correction section to AFTER trimming in `apply_QC.py` main loop. Processing order now:
1. Compass QC (outlier removal, wind direction recomputation)
2. Slow temperature QC (threshold-based filtering)
3. Deployment/retrieval period trimming (fixed duration or detection-based)
4. Slow temperature bias correction (using coefficients from already-cleaned data)
5. Diagnostic plotting (shows final state after all operations)
6. Dataset saving

**Rationale**: Bias analysis (notebook) should use data already cleaned of outliers and edge artifacts. This ensures slope/intercept calculations represent true operational behavior, not contamination from setup/teardown or sensor issues.

**Files Modified**: `analysis_scripts/apply_QC.py` (section reordering by user)

### Key Lessons Learned

#### Sequential QC Pipeline Design for Scientific Data
**Lesson**: Quality control operations must be ordered to ensure each step operates on appropriately cleaned data. Bias correction coefficients computed from data containing outliers or deployment artifacts will be corrupted. The correct sequence is: remove outliers → remove edge artifacts → calculate corrections from clean data → apply corrections. This differs from industrial QC where corrections might be applied first.

**Impact**: When designing multi-stage QC workflows for scientific instruments, document the intended processing order and rationale. Consider whether each operation should work on raw data, partially processed data, or fully cleaned data. For pyPIPS: automated filters first (catch systematic issues), then trimming (remove transient contamination), then bias correction (refine calibration). Manual QC can occur at any stage depending on what's being examined.

#### Thermodynamic Consistency After Sensor Corrections
**Lesson**: Correcting one thermodynamic variable (slowtemp) breaks relationships with derived quantities (dewpoint, RH_derived). Must recompute dependent variables using corrected inputs to maintain physical consistency. Additionally, when sensors fail (slowtemp → NaN), thermodynamic calculations need fallback strategies (use fasttemp for dewpoint calculation).

**Impact**: Any correction to a thermodynamic sensor requires propagating changes through the dependency chain:
- Slowtemp correction → recalculate dewpoint (uses slowtemp + RH) → recalculate RH_derived (uses fasttemp + dewpoint)
- Must handle missing data: `temp_for_dewpoint = slowtemp.where(~np.isnan(slowtemp), fasttemp)`
- Pattern applies to other sensor groups: if correcting pressure → recalculate air density, potential temperature, etc.
- Document these dependencies so future corrections trigger proper recalculations

#### Interactive + Batch Tool Pattern for Manual Curation
**Lesson**: Manual QC requires two distinct tool types: (1) interactive notebook for exploration and decision-making, (2) batch script for reproducible application. The notebook enables human pattern recognition (zoom, pan, inspect context), while the script ensures decisions are documented (JSON) and applied consistently (same time ranges, same logic, version controlled).

**Impact**: This pattern is applicable beyond temperature QC:
- Manual rain/no-rain event classification for DSD analysis
- Identifying storm cells vs stratiform precipitation periods
- Flagging instrument interference periods (e.g., observer blocking sensor)
- Marking periods requiring special handling (mixed phase, high winds)
- General pattern: notebook for annotation → JSON for storage → script for application
- Separating exploration from application enables review, sharing, and reapplication if processing changes

#### In-Place Variable Updates vs. Separate Variables
**Lesson**: Notebook analysis created `slowtemp_corrected` separate variables to preserve original data for comparison. Production workflow overwrites `slowtemp` directly to avoid proliferation of variable versions (`slowtemp`, `slowtemp_corrected`, `slowtemp_qc`, `slowtemp_qc_corrected`, etc.). Decision depends on use case: research → keep all versions for provenance; operations → use latest with metadata attributes.

**Impact**: For pyPIPS operations, attributes track modifications: `qc_diff_threshold`, `n_flagged`, `bias_correction_slope`, `manual_qc_applied`. This approach:
- Reduces file size (one temperature variable, not three)
- Simplifies downstream code (always use `slowtemp`, not "which version?")
- Maintains provenance through netCDF attributes
- Requires careful workflow design (can't recompute corrections without original data)
- Alternative: keep `slowtemp_raw` and update `slowtemp` working copy

### Technical Context for Future Development

**Bias Correction Architecture**:
- Configuration-driven: PIPS-specific coefficients in separate config file
- Modular application: function accepts slope/intercept, returns updated dataset
- Thermodynamic propagation: changes to slowtemp trigger dewpoint/RH_derived recalculations
- Dual dataset updates: conventional (1-Hz) gets corrected → parsivel (resampled) inherits via resampling
- Attribute documentation: `bias_correction_slope`, `bias_correction_intercept`, `bias_corrected=True`
- Fallback handling: `slowtemp.where(~np.isnan(slowtemp), fasttemp)` for dewpoint when slowtemp missing

**Manual QC JSON Schema**:
```json
{
  "PIPS_NAME_DEPLOYMENT_NAME": [
    {
      "variable": "slowtemp",
      "start_time": "2025-01-15T10:30:00",  // ISO 8601 format
      "end_time": "2025-01-15T11:45:00",
      "reason": "Sensor malfunction during storm",
      "flagged_on": "2026-02-02T15:23:45"  // Timestamp of QC decision
    }
  ]
}
```
- Key format: `{PIPS_name}_{deployment_name}` matches xarray dataset naming conventions
- Time format: ISO 8601 strings for universal parsing (`pd.to_datetime()` compatible)
- Reason field: free text for human documentation (not parsed)
- Flagged_on: audit trail for when decision was made (who: implicit from file ownership)
- Multiple entries per dataset: list accumulates all QC decisions
- Extension: could add `flagged_by` field for multi-user workflows

**QC Processing Pipeline**:
1. **Compass QC** (`--compass-qc`): Circular statistics (mean, std) → dual threshold (z-score AND absolute deviation) → flag outliers → recalculate mean → recompute wind directions
2. **Slowtemp QC** (`--slowtemp-qc`): Absolute difference |slowtemp - fasttemp| > threshold → flag → update attributes
3. **Trimming** (`--trim-deployment-periods`): Fixed duration (default 30s) from start/end → update time coordinates → update global attributes (starting_time, ending_time)
4. **Bias Correction** (`--slowtemp-bias-correction`): Load coefficients → apply to slowtemp → recompute dewpoint → recompute RH_derived → update attributes
5. **Manual QC**: Apply after automated pipeline via separate script
6. **Plotting**: Diagnostic plots show final state after all operations

Each step optional via command-line flags. Dataset attributes document which operations were performed.

**Error Handling Strategy**:
- Automated QC (`apply_QC.py`): Fail-fast with informative error messages. If config file missing or corrupt, abort. If required variables missing, log warning and skip that operation.
- Manual QC (`apply_manual_qc.py`): Graceful degradation. If variable not in dataset, log warning and continue. If time range has no data, log warning. Ensures maximum processing even with partial data availability.
- Notebook (`manual_QC_inspector.ipynb`): User-facing error messages. If file not found, print path and available files. If variable missing, list available variables.

### Next Development Priorities

1. **Extend Manual QC System to Other Variables**: Current implementation focuses on temperature (`slowtemp`), but framework supports any variable. Create helper functions or notebook sections for:
   - Rain rate anomalies (e.g., splashing, debris on sensor)
   - Wind speed spikes (e.g., observer interference)
   - Compass direction discontinuities (e.g., magnetic interference)
   - DSD matrix artifacts (e.g., specific diameter bins with persistent noise)

2. **Add Bias Correction Diagnostics**: Before/after scatter plots showing:
   - Slowtemp vs fasttemp (pre-correction): should show slope/intercept
   - Slowtemp vs fasttemp (post-correction): should show 1:1 line
   - Time series comparison: original slowtemp, corrected slowtemp, fasttemp
   - Similar to existing diagnostic plots for compass and slowtemp QC
   - Integration: add `--plot-bias-correction` flag to `apply_QC.py`

3. **Generalize Bias Correction Framework**: Current implementation is slowtemp-specific. Create generic `apply_sensor_bias_correction()` function:
   - Parameters: dataset, variable name, slope, intercept, dependent variables to recompute
   - Enables bias correction for any sensor with linear relationship to reference
   - Examples: pressure sensor drift, wind speed calibration, rain gauge catch efficiency
   - Config structure: nested dicts by variable type, PIPS name, then coefficients

4. **Automated Bias Coefficient Computation**: Currently manual workflow (notebook → visual inspection → record coefficients). Create script:
   - Load QC'd and trimmed data (after first pass)
   - Compute linear fits for each PIPS unit
   - Generate diagnostic plots (scatter, residuals, time series)
   - Output coefficients to config file format
   - Flag units with poor fits (low R², large residuals) for manual review
   - Integration: `python compute_bias_corrections.py <config.py> --output configs/bias_corrections.py`

5. **Manual QC Notebook Enhancements**:
   - Add multi-variable inspection: plot slowtemp + fasttemp + rain rate in stacked panels
   - Keyboard shortcuts: press 'a' to add QC for visible time range, 's' to save
   - Time range selection: click-drag on plot to select start/end times
   - QC summary visualization: timeline view showing all flagged periods across deployments
   - Export QC summary to report: markdown or HTML with plots + decisions table

6. **QC Audit Trail and Versioning**: Track full QC history:
   - Global attribute: `qc_history` with chronological list of operations applied
   - Example: `["2026-02-02T10:00:00: Compass QC applied (3.0σ, 20° threshold, 12 outliers)", ...]`
   - Version QC decision files: `manual_qc_decisions_v1.json`, `_v2.json`, etc.
   - Git integration: commit JSON files after each inspection session
   - Enables reconstruction: "What QC was applied to this dataset?" and "How have QC decisions evolved?"

## Session: January 29 - February 1, 2026 - Temperature and Humidity Meteogram Enhancements

### Key Code Updates

#### 1. Created Individual Temperature Meteogram Functions
**Problem**: The existing `plot_temperature_dewpoint_meteogram` function plotted both fasttemp and dewpoint together, but users sometimes need to visualize temperature sensors individually for diagnostics or sensor comparison analysis. No separate functions existed for plotting fasttemp-only or slowtemp-only meteograms.

**Solution**: Created two new plotting functions in `plotmodule.py` following the established meteogram pattern:
- `plot_fasttemp_meteogram()`: Plots only the fasttemp data with temperature-specific styling
- `plot_slowtemp_meteogram()`: Plots only the slowtemp data with support for corrected variables (`slowtemp_corrected`)
- Both functions accept same parameters as `plot_temperature_dewpoint_meteogram` for consistency
- Use existing `temp_params` styling (red fill) and temperature range from `global_plot_config_dict`
- Include 0°C reference line for phase transition context
- Return standard `(fig, ax1)` tuple for downstream manipulation

**Files Modified**: `pyPIPS/plotmodule.py`

#### 2. Created Original RH Meteogram Function
**Problem**: `plot_RH_meteogram` function automatically selected between `RH_derived` (for PIPS) and `RH` (for CU/NV2 probes), with additional logic for corrected versions. This made it impossible to directly plot the original `RH` variable for PIPS units, which is needed for comparing raw sensor output with derived calculations or validating thermodynamic consistency.

**Solution**: Created `plot_RH_original_meteogram()` function that always plots the `RH` variable:
- Removed ptype-dependent logic that switches between RH and RH_derived
- Removed corrected variable checking - always uses raw `RH` field
- Maintains same plotting parameters, axis configuration, and return signature as `plot_RH_meteogram`
- Enables direct sensor output visualization without derived calculations

**Files Modified**: `pyPIPS/plotmodule.py`

#### 3. Integrated New Functions into Meteogram Plotting Workflow
**Problem**: New plotting functions existed but weren't integrated into production workflow. Users would need to manually call functions or modify scripts to use them.

**Solution**: Added plotting steps to `plot_conv_meteograms_nc.py` for all three new functions:
- Generates `*_fasttemp.png` plots after temperature/dewpoint plots
- Generates `*_slowtemp.png` plots after fasttemp
- Generates `*_RH_original.png` plots after derived RH plots
- Follows existing naming conventions with descriptive suffixes
- Maintains same configuration (`pc.PIPS_plotting_dict`) and parameter passing as existing plots

**Files Modified**: `plotting_scripts/plot_conv_meteograms_nc.py`

#### 4. Added Comprehensive Error Handling to Plotting Workflow
**Problem**: If any plotting step failed (e.g., due to missing data variables, invalid time ranges, or file I/O issues), the entire script would abort with an exception. This meant partial dataset failures prevented generation of all other plots, losing valuable diagnostic information.

**Solution**: Wrapped all eight meteogram plotting steps in individual try-except blocks:
- Each plotting section catches generic `Exception` to handle any error type
- Logs descriptive warning messages using `utils.warning()` with PIPS name and error details
- Continues to next plotting step after logging failure
- Ensures maximum plot generation even with partial data availability
- Applies to: wind, temperature/dewpoint, RH, fasttemp, slowtemp, original RH, pressure, and compass direction meteograms

**Files Modified**: `plotting_scripts/plot_conv_meteograms_nc.py`

### Key Lessons Learned

#### Template-Based Function Creation for Consistency
**Lesson**: Creating new plotting functions by using existing functions as templates (e.g., `plot_fasttemp_meteogram` from `plot_temperature_dewpoint_meteogram`) ensures API consistency and reduces bugs. All meteogram functions now share the same parameter signature, return values, and internal structure patterns. This makes them interchangeable in workflows and predictable for users.

**Impact**: Users familiar with any meteogram function can immediately use new functions without consulting documentation. Future function additions can follow same template pattern. Suggests creating explicit "function template" documentation showing canonical structure for different function types (meteograms, validation functions, I/O functions, etc.).

#### Graceful Degradation in Batch Processing
**Lesson**: Adding error handling that logs failures but continues processing maximizes output from partially corrupt or incomplete datasets. In field campaigns, data dropouts are common - wind sensor failures, power interruptions, communication losses. Scripts that abort on first error are fragile; scripts that generate all possible plots from available data provide maximum value for quality control and event reconstruction.

**Impact**: Users can now run plotting scripts on entire campaigns knowing they'll get all possible plots even if some PIPS units have data gaps. Error messages help identify which specific plots failed and why, enabling targeted debugging. Pattern should be extended to all batch processing scripts - QC application, parameter calculations, etc. Trade-off: must ensure errors are visible (logged clearly) so users don't miss important failures.

#### Sensor-Specific Plotting for Validation
**Lesson**: Combining multiple variables in single plot (like temp + dewpoint) is excellent for operational monitoring, but validation and debugging often require isolating individual sensors. Creating sensor-specific functions (`plot_fasttemp_meteogram`, `plot_slowtemp_meteogram`, `plot_RH_original_meteogram`) enables direct comparison between sensors, identification of sensor drift, and validation of derived quantities against raw measurements.

**Impact**: Enables new diagnostic workflows - comparing fasttemp vs slowtemp to identify thermal lag, comparing RH vs RH_derived to validate dewpoint calculations, examining individual sensor time series for anomaly detection. Suggests creating similar functions for other sensor pairs/groups (wind speed vs GPS speed, compass heading vs wind direction, etc.). Pattern applicable to any multi-sensor system where sensor intercomparison matters.

### Technical Context for Future Development

**Meteogram Function Architecture**: All meteogram functions follow consistent pattern:
1. Accept parameters: `plottimes`, `conv_plot_ds`, `global_plot_config_dict`, `xlimits`, `ptype`, `use_corrected_vars` (optional), `use_plot_date`
2. Select variable(s) to plot with logic for corrected versions if applicable
3. Create figure and axis: `fig = plt.figure(figsize=(8, 3))`, `ax1 = fig.add_subplot(111)`
4. Call `plotmeteogram()` helper function with field data and parameter dictionaries
5. Configure axis with `set_meteogram_axes()` using `axparamdict` with locators, formatters, limits, labels
6. Return `(fig, ax1)` or `(fig, ax1, ax2)` for dual-axis plots

**Error Handling Strategy**: Try-except blocks implemented at plotting step level rather than individual function level. This design choice means:
- Functions can remain simple and assume valid inputs (fail-fast design)
- Caller controls error handling policy (abort vs continue)
- Logging happens at script level with context (PIPS name, plot type)
- Future enhancement could add finer-grained error handling with custom exception types for missing variables, invalid time ranges, etc.

**Variable Naming Conventions**: Temperature and humidity variables follow patterns:
- Raw sensors: `fasttemp`, `slowtemp`, `RH`
- Derived quantities: `dewpoint`, `RH_derived` (calculated from dewpoint/temperature)
- Corrected versions: `{variable}_corrected` (e.g., `slowtemp_corrected`, `dewpoint_corrected`, `RH_derived_corrected`)
- Functions check for corrected versions with `if use_corrected_vars and '{variable}_corrected' in conv_plot_ds.data_vars`

**Plot File Naming Pattern**: Output files use structured naming: `{PIPS_name}_{deployment_name}_{start_time}_{end_time}_{plot_type}.png`
- Example: `PIPS3A_IOP1_20160101000000_20160102000000_fasttemp.png`
- New plot types: `fasttemp`, `slowtemp`, `RH_original` (distinguishes from existing `RH` which is RH_derived)
- Pattern enables automated organization, glob pattern matching, and timestamp-based filtering

### Next Development Priorities

1. **Wind sensor comparison plots**: Create functions for comparing anemometer wind vs GPS-derived speed, and compass heading vs calculated wind direction
2. **Sensor drift detection**: Implement functions that plot sensor differences (fasttemp - slowtemp, RH - RH_derived) to identify calibration drift
3. **Time series anomaly detection**: Add statistical overlays to meteograms showing mean/std bounds or running percentiles to highlight outliers
4. **Multi-PIPS comparison plots**: Create functions that plot same variable from multiple PIPS units on single axis for spatial pattern analysis
5. **Corrected variable visualization**: Add flag to plotting script to generate separate plots showing original vs corrected variables side-by-side
6. **Extend error handling**: Apply same try-except pattern to DSD plotting scripts and radar plotting workflows

## Session: January 28-29, 2026 - Batch Plotting Shell Script Creation

### Key Code Updates

#### 1. Created run_plotting.sh for Automated Plot Generation
**Problem**: No unified workflow script for generating standard plots from processed PIPS data. Users had to manually run multiple plotting scripts individually with correct arguments, making it tedious to regenerate plots after reprocessing data or to create plots for multiple deployments.

**Solution**: Created new `run_plotting.sh` script modeled after `run_initial_analysis.sh` architecture:
- Orchestrates six plotting scripts: `plot_conv_meteograms_nc.py`, `plot_DSD_meteograms_nc.py`, `plot_DSD_nc.py`, `plot_PIPS_diag.py`, `plot_radar_ppi.py`, `plot_vel_D_nc.py`
- Each plotting script conditionally executed via flags (`--skip-conv-meteograms`, `--enable-radar-ppi`, etc.)
- Inherits glob pattern support for batch processing multiple configs
- Tracks success/failure per config with summary reporting
- DSD meteograms loop through multiple QC tags automatically

**Files Created**: `shell_scripts/run_plotting.sh`

**Key Features**:
- **General options**: `--plot-start-time`, `--plot-end-time`, `--plot-dir`, `--image-fmt`, `--plot-config`
- **QC tag control**: `--qc-tags` for DSD meteograms, `--nd-tags` for DSD histograms
- **DSD plot options**: `--plot-raw`, `--plot-series`, `--plot-mm-fits`, `--no-plot-full`
- **Velocity-diameter options**: `--plot-vd-raw`, `--plot-vd-series`, `--no-plot-vd-qc`, `--no-plot-vd-full`
- **Radar PPI options**: `--radar-el-req`, `--radar-fname-variant`, `--radar-input-tag`
- **Workflow**: `--log-dir`, `--verbose` for debugging

**Usage Examples**:
```bash
# Standard plotting for single config
./run_plotting.sh configs/IOP1_2016.py

# Batch plotting with multiple QC tags
./run_plotting.sh "configs/ICECHIP_IOP*.py" --qc-tags "qc roqc hoqc"

# Only meteograms (skip histogram plots)
./run_plotting.sh configs/IOP1_2016.py --skip-dsd-plots --skip-vel-d

# Full DSD analysis with time range
./run_plotting.sh configs/IOP1_2016.py --plot-start-time 20160101120000 \
  --plot-end-time 20160101180000 --plot-mm-fits --plot-series

# Include radar PPI plots
./run_plotting.sh configs/IOP1_2016.py --enable-radar-ppi --radar-el-req 0.5
```

### Key Lessons Learned

#### Script Template Reusability Accelerates Development
**Lesson**: Creating `run_plotting.sh` by adapting the `run_initial_analysis.sh` template was extremely efficient. The core structure (glob pattern handling, config looping, success/failure tracking, logging integration) transferred directly with only script-specific customizations needed (plotting script calls and their arguments). This demonstrates value of establishing robust patterns that can be replicated.

**Impact**: Reduces development time for new workflow scripts from hours to minutes. Establishes consistent user interface across pyPIPS automation tools - users familiar with one script immediately understand others. Suggests creating library of template scripts for common workflow patterns (sequential processing, parallel processing, validation, etc.).

#### Separating Analysis from Visualization Workflows
**Lesson**: Having independent scripts for data processing (`run_initial_analysis.sh`) and visualization (`run_plotting.sh`) provides flexibility users need. Analysis is computationally expensive and infrequent; plotting is fast and may be repeated with different parameters/formats. Separation allows iterating on plot aesthetics without reprocessing data, and enables generating multiple plot variants (different QC tags, time ranges, formats) from single processed dataset.

**Impact**: Users can now run full ICECHIP analysis once overnight, then quickly generate multiple plot sets the next day for papers, presentations, diagnostics. Script supports both "quick look" (default plots only) and "publication quality" (all QC tags, multiple formats) workflows. Pattern applicable to any scientific workflow with expensive computation and flexible visualization.

#### Optional vs Required Plotting Script Arguments
**Lesson**: Different plotting scripts have different required/optional arguments. DSD meteograms require `--QC-tag`, DSD histograms use `--ND-tags`, radar PPIs need elevation angles. Solution: script provides sensible defaults and allows per-script customization via command-line flags. DSD meteograms loop through space-separated QC tags; other scripts use first tag or specific flags.

**Impact**: Balances convenience (works out-of-box with defaults) with control (all script options accessible). Documentation crucial - help text provides examples showing which options affect which plots. Future enhancement: could add "presets" like `--preset publication` that configures all scripts for paper-quality output.

### Technical Context for Future Development

**Plotting Script Argument Patterns**: All plotting scripts share common arguments (`case_config_path`, `--plot-config-path`, `--plot-dir`, `--plot-start-time`, `--plot-end-time`) but diverge in specifics. DSD scripts use QC tags to select processed data versions; radar scripts need radar-specific parameters (elevation angle, filename conventions). Understanding these patterns was key to designing flexible common argument handling (`COMMON_PLOT_ARGS`) plus script-specific argument building.

**QC Tag Iteration Strategy**: DSD meteograms script called once per QC tag in loop, generating separate plots for each. Alternative approach would be passing all tags to script and having it loop internally. Current approach gives finer control (can log each tag separately, fail gracefully on missing data) but requires more shell logic. Trade-off appropriate for batch processing where monitoring individual steps matters.

**Radar PPI as Optional Component**: Radar plots disabled by default because (1) require external radar data not always available, (2) are slow due to cartopy/pyart overhead, (3) less frequently needed than core PIPS plots. Failure treated as warning rather than error. Pattern establishes precedent for integrating optional/experimental plotting scripts without cluttering default workflow.

**Plot Configuration File Handling**: All plotting scripts accept `--plot-config-path` for customizing aesthetics (color schemes, fonts, axis limits). Wrapper script passes this uniformly to all plotting scripts. Enables users to define campaign-specific plotting styles once in config and apply consistently across all plot types. Current implementation uses `plot_config.py` in working directory by default.

### Next Development Priorities

1. **Plot presets**: Add `--preset` option with predefined combinations (e.g., `quicklook`, `publication`, `diagnostic`) that configure multiple scripts at once
2. **Selective plotting by PIPS unit**: Add ability to plot only specific PIPS units from multi-unit deployments (e.g., `--pips-names "PIPS3A PIPS3B"`)
3. **Parallel plotting**: Since plotting scripts are independent, could parallelize across configs or across plot types for faster batch processing
4. **Plot inventory report**: Generate summary showing which plots were created for each config (useful for verifying complete plot generation)
5. **Animation generation**: Add option to create time-lapse animations from time series plots (requires ffmpeg integration)
6. **Web gallery generation**: Auto-generate HTML page organizing all plots by config/plot type for easy browsing

## Session: January 25-26, 2026 - Batch Analysis Shell Script Enhancement

### Key Code Updates

#### 1. Enhanced Shell Script for Glob Pattern Support
**Problem**: `run_initial_analysis.sh` only processed single config files, requiring manual iteration to analyze multiple deployments (e.g., all ICECHIP IOPs). This was inefficient for batch processing campaigns with many similar deployments.

**Solution**: Modified script to accept glob patterns and process multiple configs sequentially:
- Parse glob pattern using bash array expansion: `CONFIG_FILES=($CASE_CONFIG_PATTERN)`
- Wrap entire analysis workflow in for-loop over matched configs
- Track success/failure per config in arrays: `SUCCESSFUL_CONFIGS` and `FAILED_CONFIGS`
- Skip remaining steps for a config if earlier step fails (check `CONFIG_SUCCESS` flag)
- Continue to next config on failure instead of exiting entire script
- Display batch summary showing which configs succeeded/failed

**Files Modified**: `shell_scripts/run_initial_analysis.sh`

**Key Features Added**:
- Glob pattern matching with nullglob to handle no-match cases
- Progress indicators showing "Processing config [N/M]" for multi-config runs
- Per-config timestamp for log file differentiation
- Final summary table listing successful and failed configs
- Exit code reflects whether any configs failed (non-zero if failures occurred)

**Usage Examples**:
```bash
# Process all ICECHIP IOPs
./run_initial_analysis.sh "configs/ICECHIP_IOP*.py"

# Process specific IOP range
./run_initial_analysis.sh "configs/ICECHIP_IOP[1-5]*.py"

# Single config (backward compatible)
./run_initial_analysis.sh configs/ICECHIP_IOP5_2025_10s.py
```

### Key Lessons Learned

#### Shell Script Glob Pattern Handling Best Practices
**Lesson**: When designing CLI tools that accept file paths, supporting glob patterns provides major productivity gains for batch operations. Key implementation points:
1. Use bash arrays with nullglob option to expand patterns safely
2. Quote glob patterns in documentation examples to prevent premature expansion
3. Provide clear feedback about number of files matched before processing
4. Continue on individual failures rather than aborting entire batch

**Impact**: Enables efficient batch processing without wrapper scripts or manual loops. Pattern applicable to any workflow script that processes config files or data files. Particularly valuable for field campaign data where many similar deployments need identical processing.

#### Failure Handling in Sequential Processing Pipelines
**Lesson**: Multi-step analysis pipelines should distinguish between "fail entire batch" vs "fail one config and continue" strategies. For independent configs, continuing through failures and reporting summary is more useful than stopping at first error. Within a single config's pipeline, should skip remaining steps if dependency fails.

**Impact**: Establishes robust pattern for pyPIPS batch processing. User can run overnight batch jobs and get comprehensive success/failure report in morning rather than finding script stopped at first error. Logging becomes critical for post-mortem debugging of failures. Trade-off: may need explicit "fail-fast" mode for development/testing.

#### Maintaining Script Backward Compatibility
**Lesson**: When enhancing scripts with new features (glob support), preserve exact behavior for existing use cases. Single file path should work identically before and after modifications. This allows gradual adoption and prevents breaking existing workflows/documentation.

**Impact**: pyPIPS analysis workflows documented in papers, internal procedures, and shell history all continue working. Users can adopt glob patterns gradually as they learn about feature. Suggests general principle: additive enhancements over breaking changes unless absolutely necessary.

### Technical Context for Future Development

**Shell Script Architecture**: `run_initial_analysis.sh` orchestrates standard pyPIPS analysis chain: CSV→netCDF → merge realtime → QC → derived params → MM fits → (optional) radar interp. Each step conditionally executed based on flags (`--skip-csv2nc`, etc.). All steps share common config file argument and optional output tag.

**Config File Organization**: pyPIPS configs follow pattern `configs/<campaign>_<event>_<temporal_res>.py`. ICECHIP campaign has configs like `ICECHIP_IOP{1,2,3,...,23}_2025_10s.py` making them ideal candidates for glob pattern matching. Each config defines `PIPS_IO_dict` with deployment metadata and `PIPS_qc_dict` with QC parameters.

**Logging Integration**: Script supports optional logging via `--log-dir` flag. When enabled, each step's stdout/stderr captured to separate files with pattern `{config_name}_{timestamp}_{step}.{stdout|stderr}`. Critical for debugging batch failures. Timestamp generated once per config to group all step logs together.

**Error Propagation**: Uses `CONFIG_SUCCESS` flag (1=success, 0=failure) that propagates through pipeline. Each step checks flag before executing and sets to 0 on failure. This prevents cascade of failures from intermediate missing files while allowing independent configs to proceed.

### Next Development Priorities

1. **Add parallel processing option**: Consider `--parallel N` flag to run N configs simultaneously using GNU parallel or bash background jobs, significantly speeding up batch analysis
2. **Config validation step**: Add optional `--validate` mode that checks config files for required parameters before starting analysis
3. **Resume capability**: Track completed configs in state file to allow resuming interrupted batch runs
4. **Progress bars**: Integrate progress indicators for long-running steps (especially for multi-config batches)
5. **Email notifications**: Add optional notification on batch completion/failure for overnight runs
6. **Dry-run mode**: Add `--dry-run` to show what would be processed without executing

## Session: January 21-25, 2026 - Thermodynamic QC with Gust Front Filtering

### Key Code Updates

#### 1. Fixed xarray Concatenation Error with 'flagged_times' Coordinate
**Problem**: ValueError when concatenating temperature DataArrays from multiple ICECHIP deployments: `coordinate 'flagged_times' not present in all datasets`
**Solution**:
- Added conditional coordinate dropping before concatenation: `if 'flagged_times' in da.coords: da = da.drop_vars('flagged_times')`
- Applied to both slowtemp and fasttemp DataArrays in concatenation loop
- Allows clean merging of datasets with heterogeneous coordinates

**Files Modified**: `notebooks/pyPIPS_thermo_QC_slowtemp_only_all_ICECHIP.ipynb` (Cell 5)

#### 2. Implemented Rapid Temperature Change Detection Algorithm
**Problem**: Scatter plots of slowtemp vs fasttemp showed "filament-like" structures from gust fronts where sensors diverged due to different response times, contaminating systematic bias estimation
**Solution**:
- Calculate temperature derivative: `dT_dt = np.abs(fasttemp_diff / time_diff_seconds)` in °C/s
- Apply 60-timestep rolling window mean to identify sustained rapid changes (not just noise)
- Threshold at 0.01 °C/s (36 °C/hour) to flag gust front passages
- Filter both temperature arrays using mask: `fasttemp_filtered = all_fasttemp_da.where(~rapid_change_mask_full, drop=True)`
- Recompute regression statistics on filtered data showing improved R² and reduced RMSE

**Files Modified**: `notebooks/pyPIPS_thermo_QC_slowtemp_only_all_ICECHIP.ipynb` (Cells 7-8)

**Key Parameters**:
- `rapid_change_threshold = 0.01  # °C/s`
- `window_size = 60  # timesteps (~60 seconds at 1s intervals)`

#### 3. Added Deployment Boundary Detection for Concatenated Data
**Problem**: Concatenating data from multiple IOPs creates large time gaps that could cause rolling window to span across deployment boundaries, mixing unrelated data
**Solution**:
- Detect deployment boundaries: `deployment_boundary_mask = time_diff_seconds > max_normal_gap` (60s threshold)
- Set dT_dt to NaN at boundaries: `dT_dt = dT_dt.where(~deployment_boundary_mask)`
- Prevents rolling mean from calculating across deployment gaps
- Reports number and locations of detected boundaries with gap durations

**Files Modified**: `notebooks/pyPIPS_thermo_QC_slowtemp_only_all_ICECHIP.ipynb` (Cell 7)

#### 4. Fixed xarray/NumPy Type Mismatch in Boundary Detection
**Problem**: AttributeError: `'numpy.int64' object has no attribute 'values'` when treating numpy array as xarray DataArray
**Solution**:
- Keep `time_diff_seconds` as xarray DataArray: `time_diff_seconds = time_diff / np.timedelta64(1, 's')` (removed premature `.values`)
- `deployment_boundary_mask` inherits DataArray type from comparison operation
- Access underlying numpy arrays only when needed: `deployment_boundary_mask.values` for indexing
- Convert scalar results properly: `int(deployment_boundary_mask.sum())` instead of `.sum().values`

**Files Modified**: `notebooks/pyPIPS_thermo_QC_slowtemp_only_all_ICECHIP.ipynb` (Cell 7)

### Key Lessons Learned

#### Sensor Response Time Creates Transient Bias During Rapid Changes
**Lesson**: Different thermal response times between slow and fast temperature sensors create apparent temperature differences during rapid environmental changes (gust fronts) that are NOT representative of steady-state systematic bias. These transient effects appear as "filament" structures in scatter plots and must be filtered before regression analysis.
**Impact**: Critical distinction between systematic bias (correctable with linear regression) and transient response artifacts (must be excluded from calibration). Algorithm filters ~1-10% of data depending on weather conditions while preserving true bias signal. Applicable to any paired sensor calibration where response times differ.

#### Rolling Window Analysis Requires Boundary Protection
**Lesson**: When applying rolling window operations to concatenated timeseries from multiple deployments, must explicitly detect and handle temporal discontinuities. Setting values to NaN at boundaries prevents window from spanning gaps and mixing unrelated data.
**Impact**: Establishes pattern for all future multi-deployment analyses. Simple gap detection (time_diff > threshold) combined with NaN masking ensures rolling statistics respect temporal structure. Alternative approach would be to group by deployment first, but masking is computationally simpler.

#### xarray/NumPy Type Preservation in Calculation Chains
**Lesson**: When building calculation chains with xarray DataArrays, premature conversion to numpy arrays (`.values`) breaks xarray operations like `.where()`. Keep data as DataArrays through computation pipeline and only extract numpy arrays at final output/indexing steps.
**Impact**: Debugging type mismatches can be avoided by understanding xarray broadcasting rules and preserving DataArray types. Performance impact is minimal due to lazy evaluation. Pattern: compute in xarray, extract to numpy only for non-xarray libraries or final outputs.

#### Circular Statistics for Directional Data (Compass QC)
**Lesson**: Regular statistics fail for circular/directional data due to 360°/0° discontinuity. Must use vector averaging for mean (arctan2 of sin/cos components) and resultant vector length for standard deviation calculation.
**Impact**: Critical for compass and wind direction QC. Dual-threshold approach (z-score AND absolute deviation) prevents over-cleaning of stable data while catching genuine outliers.

### Technical Context for Future Development

**Thermodynamic QC Workflow**: Temperature correction uses inverse linear regression: `slowtemp_corrected = (slowtemp - intercept) / slope`. This corrected temperature then propagates through dewpoint calculation (`calTdfromRH`) and RH_derived recalculation (`calRH`), ensuring thermodynamic consistency.

**Filtering vs Correction Decision**: Gust front filtering is applied during CALIBRATION (regression fitting) but not during CORRECTION (applying regression to data). This preserves all temporal data while ensuring calibration parameters represent steady-state bias only.

**Data Concatenation Pattern**: Multi-deployment analysis concatenates individual deployment DataArrays using `xr.concat(das, dim='time')`. Must check for and drop non-universal coordinates (like 'flagged_times') before concatenation to avoid ValueError.

**Statistical Comparison Framework**: Side-by-side scatter plots with comprehensive statistics tables (N, bias, RMSE, r, R², slope, intercept) provide clear before/after validation. Filtering improved R² by ~0.01-0.05 percentage points and reduced RMSE by ~5-15% depending on weather conditions.

### Next Development Priorities

1. **Apply filtered regression parameters**: Consider using filtered regression results (slope_filt, intercept_filt) instead of original parameters for final temperature correction
2. **Complete thermodynamic calculations**: Compute remaining derived parameters (potential temperature, water vapor mixing ratio, density) using corrected slowtemp
3. **Temporal resampling**: Resample corrected 1s conventional data to 10s parsivel temporal resolution for merged analysis
4. **Save corrected datasets**: Write corrected thermodynamic variables back to netCDF files with proper attributes documenting QC methods
5. **Extend to other PIPS units**: Apply same QC workflow to PIPS1A, 1B, 2A, 2B, 3B for complete ICECHIP dataset quality control

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