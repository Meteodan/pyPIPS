# Automatic Deployment Period Detection and Trimming

## Overview

The `apply_QC.py` script can automatically detect and remove periods at the beginning and end of PIPS deployments where the probe is being placed or retrieved. These periods often show rapid fluctuations in compass headings that contaminate the dataset.

## Key Features

1. **Automatic Detection**: Identifies rapid compass fluctuations at start/end of timeseries
2. **Time Alignment**: Ensures conventional and parsivel datasets remain synchronized
3. **Metadata Updates**: Automatically updates `starting_time` and `ending_time` attributes
4. **Configurable**: Adjustable detection parameters for different deployment conditions

## How It Works

### Detection Algorithm

1. **Sliding Window Analysis**: Examines compass data in overlapping windows
2. **Standard Deviation Check**: Calculates std dev of compass direction in each window
3. **Threshold Comparison**: Flags windows exceeding the stability threshold
4. **Trim Point Identification**: Finds first/last stable periods

### Trimming Process

1. **Conventional Dataset**: Remove unstable periods from 1-Hz data
2. **Parsivel Dataset**: Remove corresponding resampled periods
3. **Time Alignment**: Trim conventional to match parsivel time bounds
4. **Attribute Update**: Update time metadata in both datasets

## Command-Line Arguments

### Enable Trimming
```bash
--trim-deployment-periods
```
Enable automatic detection and removal of deployment/retrieval periods.

### Detection Parameters

**`--trim-window-size SIZE`** (default: 60)
- Window size in seconds for calculating compass standard deviation
- Larger windows = more robust to short-term fluctuations
- Smaller windows = more sensitive to changes

**`--trim-std-threshold THRESHOLD`** (default: 10.0)
- Standard deviation threshold in degrees for flagging fluctuations
- Higher values = more data kept (less aggressive trimming)
- Lower values = more data removed (more aggressive trimming)

**`--trim-max-duration DURATION`** (default: 600)
- Maximum duration in seconds to check at start and end
- Default checks first/last 10 minutes
- Prevents excessive trimming if compass is always unstable

## Usage Examples

### Basic Usage

Trim deployment periods with default settings:
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --trim-deployment-periods \
    --output-QC-tags qc \
    --output-file-tag _trimmed
```

### Sensitive Detection

More aggressive trimming for problematic deployments:
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --trim-deployment-periods \
    --trim-std-threshold 5.0 \
    --trim-window-size 30 \
    --output-QC-tags qc
```

### Conservative Detection

Less aggressive trimming, only remove obvious problems:
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --trim-deployment-periods \
    --trim-std-threshold 15.0 \
    --trim-window-size 120 \
    --output-QC-tags qc
```

### Combined with Compass QC

Apply both trimming and compass QC (recommended workflow):
```bash
python analysis_scripts/apply_QC.py \
    configs/ICECHIP_IOP2_2025_10s.py \
    --trim-deployment-periods \
    --compass-qc \
    --plot-compass-qc \
    --output-QC-tags qc \
    --output-file-tag _full_qc
```

## Output

### Modified Files

**Conventional Dataset** (`conventional_raw_*{output_file_tag}.nc`):
- Trimmed to stable period
- Updated `starting_time` and `ending_time` attributes

**Parsivel Dataset** (`parsivel_combined_*{output_file_tag}.nc`):
- Trimmed to match conventional dataset
- Updated time attributes
- All DSD variables remain aligned

### Log Messages

Example output:
```
Detecting deployment/retrieval periods for PIPS3A...
    Detected deployment period: trimming first 120 seconds
    Detected retrieval period: trimming last 85 seconds
    Trimming datasets for PIPS3A...
    Conventional dataset trimmed: 8640 -> 8435 records
    Parsivel dataset trimmed: 864 -> 844 records
    Conventional dataset aligned to parsivel times: 8435 records
    Updated time attributes: 20250519120200 to 20250519235959
    Saving trimmed conventional data: /path/to/conventional_raw_IOP2_051925_PIPS3A_trimmed.nc
```

## Configuration Requirements

Your config file must include `conv_filenames_nc`:

```python
PIPS_IO_dict = {
    'PIPS_dir': '/path/to/netcdf/',
    'PIPS_names': ['PIPS1A', 'PIPS1B', 'PIPS3A', 'PIPS3B'],
    'PIPS_filenames_nc': ['parsivel_combined_*.nc', ...],
    'conv_filenames_nc': ['conventional_raw_*.nc', ...],  # Required!
}
```

## Workflow Integration

### Recommended Processing Order

1. **Trim deployment periods** (this feature)
2. **Apply compass QC** (if needed)
3. **Apply DSD QC** (standard QC)
4. **Calculate derived parameters**
5. **Generate visualizations**

### Combined Command

```bash
# Full QC workflow
python analysis_scripts/apply_QC.py \
    configs/YOUR_CONFIG.py \
    --trim-deployment-periods \
    --compass-qc \
    --compass-zscore-threshold 0.5 \
    --output-QC-tags qc \
    --output-file-tag _full_qc

# Then calculate derived parameters on QC'd data
python analysis_scripts/calc_derived_params.py \
    configs/YOUR_CONFIG.py \
    --input-QC-tag qc \
    --input-file-tag _full_qc
```

## Choosing Parameters

### Window Size (`--trim-window-size`)

| Value | Use Case |
|-------|----------|
| 30s | Fast deployments, need sensitivity to quick changes |
| 60s | **Default** - good balance for typical deployments |
| 120s | Slow deployments, more robust averaging |

### Std Dev Threshold (`--trim-std-threshold`)

| Value | Description |
|-------|-------------|
| 5.0° | Aggressive - only keep very stable periods |
| 10.0° | **Default** - good for typical PIPS3A/3B issues |
| 15.0° | Conservative - only remove obvious problems |

### Max Duration (`--trim-max-duration`)

| Value | Use Case |
|-------|----------|
| 300s | Check only first/last 5 minutes |
| 600s | **Default** - first/last 10 minutes |
| 1200s | Extended check for problematic deployments |

## Technical Details

### Time Synchronization

The trimming algorithm ensures proper time alignment:

1. Detects trim points using 1-Hz conventional compass data
2. Applies same time trim to resampled parsivel data
3. Further trims conventional to match parsivel time bounds
4. Result: Both datasets cover identical time period

### Circular Statistics

Standard deviation is calculated directly on compass angles. For datasets crossing 360°/0°, this may produce artificially high std dev values. If this is an issue, the algorithm can be enhanced to use circular statistics.

### Edge Cases

- **No unstable periods detected**: Datasets unchanged, normal processing continues
- **Very short deployments**: Skipped if less than 2× window size
- **Entire deployment unstable**: Only trims up to `max_duration` from each end
- **Missing conventional data**: Trimming skipped for that PIPS

## Troubleshooting

### Too Much Data Removed

If excessive trimming occurs:
- Increase `--trim-std-threshold` (e.g., from 10.0 to 15.0)
- Increase `--trim-window-size` (e.g., from 60 to 120)
- Reduce `--trim-max-duration` (e.g., from 600 to 300)

### Not Enough Data Removed

If deployment/retrieval periods remain:
- Decrease `--trim-std-threshold` (e.g., from 10.0 to 5.0)
- Decrease `--trim-window-size` (e.g., from 60 to 30)
- Increase `--trim-max-duration` (e.g., from 600 to 1200)

### Verify Results

Check the log output for:
- Number of seconds trimmed
- Before/after record counts
- Updated time attributes

Visually inspect:
```bash
# Plot compass timeseries before and after
python plotting_scripts/plot_PIPS_diag.py \
    configs/YOUR_CONFIG.py \
    --plot-conv
```

## Benefits

1. **Cleaner Data**: Removes contaminated deployment/retrieval periods
2. **Automatic**: No manual inspection of each deployment needed
3. **Consistent**: Same criteria applied across all PIPS
4. **Traceable**: Log messages document what was removed
5. **Synchronized**: Maintains alignment between conventional and parsivel data

## Limitations

- Relies on compass stability as proxy for deployment issues
- May not detect all types of deployment problems
- Cannot recover data from trimmed periods
- Requires conventional dataset with compass data

## Related Features

- **Compass QC** (`--compass-qc`): Removes outliers from stable periods
- **DSD QC** (`--output-QC-tags`): Standard precipitation data quality control
- **Diagnostic Plots** (`--plot-compass-qc`): Visualize compass QC results

## References

Developed to address deployment/retrieval artifacts observed in ICECHIP and PERiLS field campaigns, particularly affecting PIPS3A and PIPS3B units with intermittent compass issues during handling.
