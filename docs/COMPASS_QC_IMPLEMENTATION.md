# Compass Quality Control Implementation Summary

## Changes Made

Successfully integrated compass quality control functionality from `pyPIPS_wind_and_compass_QC.ipynb` into the `analysis_scripts/apply_QC.py` production script.

## Modified Files

### 1. analysis_scripts/apply_QC.py

**New Imports:**
- `matplotlib.pyplot as plt` - for diagnostic plotting
- `from scipy.stats.mstats import zscore` - for outlier detection

**New Command-Line Arguments:**
- `--compass-qc`: Flag to enable compass quality control
- `--compass-zscore-threshold FLOAT`: Threshold for outlier detection (default: 0.5)
- `--plot-compass-qc`: Flag to generate diagnostic plots

**New Functions:**

1. **`resample_compass_da(compass_dir, offset, intervalstr)`**
   - Resamples compass directions using vector averaging
   - Handles 360°/0° discontinuity properly
   - Returns resampled compass direction DataArray

2. **`apply_compass_qc(conv_ds, zscore_threshold=0.5)`**
   - Applies z-score based outlier detection to compass data
   - Replaces outliers with NaN
   - Returns cleaned compass direction and mean value

3. **`plot_compass_qc_diagnostics(...)`**
   - Generates 2x2 diagnostic plots (before/after compass and wind)
   - Saves PNG files to plot directory
   - Shows original vs QC'd timeseries

**New Processing Section:**

Added compass QC workflow after DSD QC loop (~70 lines):
- Loads conventional dataset for each PIPS
- Applies compass QC to remove outliers
- Recomputes absolute wind directions
- Resamples data to parsivel times if available
- Updates both conventional and parsivel datasets
- Generates diagnostic plots if requested
- Saves updated datasets

**Configuration Support:**
- Reads `conv_filenames_nc` from config dictionary
- Falls back gracefully if conventional files not available
- Compatible with existing config structure

## Created Files

### 2. docs/COMPASS_QC_USAGE.md

Comprehensive user documentation including:
- Overview and features
- Command-line usage examples
- Technical details of the algorithm
- Configuration requirements
- Troubleshooting guide
- Integration with existing workflow

## Key Features

### 1. Automatic Processing
- Processes all PIPS in a case configuration automatically
- No need to manually specify PIPS names or files
- Integrates seamlessly with existing apply_QC.py workflow

### 2. Configurable Thresholds
- Default z-score threshold: 0.5 (suitable for PIPS3A/3B)
- User can adjust based on data quality needs
- More conservative or aggressive filtering as needed

### 3. Comprehensive Updates
- Updates 1-Hz conventional datasets
- Resamples and updates parsivel datasets (10s, etc.)
- Preserves all metadata and adds QC attributes
- Recomputes all wind-related variables consistently

### 4. Diagnostic Visualization
- Optional before/after plots for quality assessment
- Shows both compass and wind direction changes
- One plot per PIPS for easy comparison
- Saved to plot directory specified in config

### 5. Robust Error Handling
- Graceful fallback if conventional files missing
- Handles cases with no parsivel data
- Logs informative messages for debugging
- Continues processing even if individual PIPS fails

## Usage Pattern

### Basic Usage:
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --compass-qc \
    --output-QC-tags qc
```

### With Diagnostic Plots:
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --compass-qc \
    --plot-compass-qc \
    --output-QC-tags qc \
    --output-file-tag _compass_qc
```

## Technical Approach

### Outlier Detection
Uses z-score method from scipy:
```python
compass_zscores = zscore(conv_ds['compass_dir'], nan_policy='omit')
cleaned_compass_dir = conv_ds['compass_dir'].where(np.abs(compass_zscores) < threshold)
```

### Wind Direction Correction
Recomputes using cleaned mean compass:
```python
avg_compass_dir = cleaned_compass_dir.mean(dim='time', skipna=True)
winddirabs_new = np.mod(avg_compass_dir + conv_ds['winddirrel'], 360.)
```

### Vector Averaging for Resampling
Avoids 360°/0° discontinuity:
```python
x = np.cos(np.deg2rad(-compass_dir + 270.))
y = np.sin(np.deg2rad(-compass_dir + 270.))
# ... resample x and y components ...
compass_dir_avg = (270.0 - (180. / np.pi) * np.arctan2(y_avg, x_avg)) % 360.
```

## Updated Dataset Variables

### Conventional Dataset:
- `compass_dir`: Outliers → NaN, attrs added (average, qc_zscore_threshold)
- `winddirabs`: Recomputed from cleaned compass

### Parsivel Dataset:
- `compass_dir`: Resampled cleaned data, attrs updated
- `winddirabs`: Resampled corrected wind direction
- `windspd`, `windspdavgvec`, `winddirunitavgvec`: Updated
- `windgust`, `uavg`, `vavg`, `unit_uavg`, `unit_vavg`: Updated

## Benefits

1. **Automation**: Single command processes entire case
2. **Reproducibility**: Consistent QC across all deployments
3. **Traceability**: QC parameters stored in file metadata
4. **Validation**: Diagnostic plots enable visual verification
5. **Integration**: Works within existing pyPIPS workflow
6. **Flexibility**: Adjustable thresholds for different conditions

## Testing Recommendations

1. Run with `--plot-compass-qc` first to verify QC effectiveness
2. Check diagnostic plots to ensure appropriate outliers removed
3. Adjust `--compass-zscore-threshold` if needed
4. Compare original vs QC'd wind directions
5. Verify resampled parsivel data updated correctly

## Future Enhancements (Optional)

Possible additions from notebook TODO:
- Interpolation across NaN gaps (currently just removes outliers)
- Handle compass crossing 360°/0° during interpolation
- Add option to use median instead of mean
- Support for multiple QC iterations

## Related Documentation

- Original notebook: `notebooks/pyPIPS_wind_and_compass_QC.ipynb`
- Usage guide: `docs/COMPASS_QC_USAGE.md`
- Example configs: `configs/ICECHIP_IOP*_2025_10s.py`
