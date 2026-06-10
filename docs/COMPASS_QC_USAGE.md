# Compass Quality Control in apply_QC.py

## Overview

The `apply_QC.py` script now includes functionality to perform quality control on PIPS compass headings to remove outliers and recompute wind directions. This is particularly important for PIPS3A and PIPS3B which can have intermittent compass issues.

## Features

1. **Outlier Detection**: Uses z-score method to identify and remove compass heading outliers
2. **Wind Direction Correction**: Recomputes absolute wind directions using cleaned compass data
3. **Data Resampling**: Updates both 1-Hz conventional data and resampled parsivel data
4. **Diagnostic Plotting**: Optional before/after plots showing compass and wind direction QC

## Command-Line Arguments

### New Compass QC Options

- `--compass-qc`: Enable compass quality control (default: disabled)
- `--compass-zscore-threshold FLOAT`: Z-score threshold for outlier detection (default: 0.5)
  - Lower values (e.g., 0.3) remove more outliers
  - Higher values (e.g., 1.0) are more conservative
- `--plot-compass-qc`: Generate diagnostic plots (default: no plots)

## Usage Examples

### Basic Compass QC

Apply compass QC with default settings (z-score threshold = 0.5):

```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --compass-qc \
    --output-QC-tags qc \
    --output-file-tag _compass_qc
```

### Compass QC with Custom Threshold

Use a more conservative threshold (1.0 standard deviations):

```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --compass-qc \
    --compass-zscore-threshold 1.0 \
    --output-QC-tags qc \
    --output-file-tag _compass_qc
```

### Compass QC with Diagnostic Plots

Generate before/after diagnostic plots:

```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --compass-qc \
    --plot-compass-qc \
    --output-QC-tags qc \
    --output-file-tag _compass_qc
```

### Combined DSD and Compass QC

Apply both standard DSD QC and compass QC:

```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --output-QC-tags qc \
    --compass-qc \
    --compass-zscore-threshold 0.5 \
    --plot-compass-qc \
    --output-file-tag _full_qc
```

## How It Works

### 1. Outlier Detection

The script uses z-score statistics to identify compass outliers:

```python
compass_zscores = zscore(conv_ds['compass_dir'], nan_policy='omit')
cleaned_compass_dir = conv_ds['compass_dir'].where(np.abs(compass_zscores) < threshold)
```

Any compass reading more than `threshold` standard deviations from the mean is replaced with NaN.

### 2. Mean Compass Calculation

The average compass direction is computed from the cleaned data:

```python
avg_compass_dir = cleaned_compass_dir.mean(dim='time', skipna=True)
```

### 3. Wind Direction Correction

Absolute wind directions are recomputed using the cleaned compass:

```python
winddirabs_new = np.mod(avg_compass_dir + conv_ds['winddirrel'], 360.)
```

### 4. Data Resampling (if parsivel data exists)

The cleaned compass and corrected winds are resampled to match parsivel times:

- Compass directions use vector averaging to handle 360°/0° discontinuity
- Winds are resampled using the standard `resample_wind_da()` function
- Both conventional and parsivel datasets are updated

## Output Files

### Modified Variables

**Conventional Dataset (`conventional_raw_*.nc`)**:
- `compass_dir`: Updated with outliers replaced by NaN
  - New attribute: `average` - mean compass direction after QC
  - New attribute: `qc_zscore_threshold` - threshold used for outlier detection
- `winddirabs`: Recomputed absolute wind direction

**Parsivel Dataset (`parsivel_combined_*_10s.nc`)**:
- `compass_dir`: Resampled cleaned compass directions
  - Updated attribute: `average` - mean resampled compass direction
  - New attribute: `qc_zscore_threshold` - threshold used
- `winddirabs`: Resampled corrected wind direction
- `windspd`, `windspdavgvec`, `winddirunitavgvec`: Updated wind variables
- `windgust`, `uavg`, `vavg`: Updated wind components

### Diagnostic Plots

When `--plot-compass-qc` is used, a 2x2 plot is generated for each PIPS showing:

1. **Top-left**: Original compass direction timeseries
2. **Top-right**: QC'd compass direction with mean marked
3. **Bottom-left**: Original wind direction timeseries
4. **Bottom-right**: QC'd wind direction timeseries

Plots are saved to the `plot_dir` specified in the config file:
- Filename format: `compass_wind_QC_{deployment_name}_{PIPS_name}.png`

## Configuration Requirements

Your config file must include `conv_filenames_nc` for compass QC to work:

```python
PIPS_IO_dict = {
    'dataset_name': 'IOP2_051925',
    'PIPS_dir': '/path/to/netcdf/',
    'plot_dir': '/path/to/plots/',
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS2A', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames_nc': [
        'parsivel_combined_IOP2_051925_PIPS1A_10s.nc',
        'parsivel_combined_IOP2_051925_PIPS1B_10s.nc',
        # ... more files
    ],
    'conv_filenames_nc': [
        'conventional_raw_IOP2_051925_PIPS1A.nc',
        'conventional_raw_IOP2_051925_PIPS1B.nc',
        # ... more files
    ],
    # ... other settings
}
```

## Workflow Integration

Typical processing workflow with compass QC:

1. **Combine raw data** → netCDF files
2. **Apply QC** (this script with `--compass-qc`)
3. **Calculate derived parameters** → `calc_derived_params.py`
4. **Generate visualizations** → plotting scripts

## Troubleshooting

### Warning: "Conventional file not found"

If you see this warning, the script cannot find the conventional data file:
- Check that `conv_filenames_nc` is defined in your config
- Verify the files exist in `PIPS_dir`
- Compass QC will be skipped for that PIPS

### No Parsivel Data

If parsivel data doesn't exist for a PIPS:
- Conventional data will still be updated with compass QC
- Only the 1-Hz conventional dataset will be saved
- No resampled data will be generated

### Choosing the Right Threshold

Start with the default (0.5) and check diagnostic plots:
- If too many good values are removed, increase threshold (e.g., 0.7, 1.0)
- If obvious outliers remain, decrease threshold (e.g., 0.3)
- For PIPS3A/3B with severe compass issues, 0.5 typically works well

## Technical Details

### Vector Averaging

Compass directions are resampled using vector averaging to properly handle the 360°/0° discontinuity:

```python
x = np.cos(np.deg2rad(-compass_dir + 270.))
y = np.sin(np.deg2rad(-compass_dir + 270.))
x_avg = x.resample(...).mean()
y_avg = y.resample(...).mean()
compass_dir_avg = (270.0 - (180. / np.pi) * np.arctan2(y_avg, x_avg)) % 360.
```

This ensures that averaging compass directions like 350° and 10° correctly yields ~0° rather than 180°.

### Metadata Preservation

The script preserves all original metadata and adds QC-specific attributes:
- Z-score threshold used
- Mean compass direction before and after QC
- QC flags remain with all variables

## Related Scripts

- **Notebook**: `notebooks/pyPIPS_wind_and_compass_QC.ipynb` - Original development/testing
- **Config examples**: `configs/ICECHIP_IOP*_2025_10s.py` - Sample configuration files
- **Plotting**: `plotting_scripts/plot_PIPS_diag.py` - View diagnostic variables including compass

## References

Original implementation developed in notebook to address compass issues observed during ICECHIP and PERiLS field campaigns, particularly affecting PIPS3A and PIPS3B units.
