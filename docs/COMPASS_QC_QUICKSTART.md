# Quick Start: Compass QC in apply_QC.py

## What it does
Removes compass heading outliers (common in PIPS3A/3B) and recomputes wind directions.

## Quick Commands

### 1. Basic compass QC (recommended first run)
```bash
python analysis_scripts/apply_QC.py \
    configs/YOUR_CONFIG.py \
    --compass-qc \
    --plot-compass-qc \
    --output-QC-tags qc
```

Check the plots in your `plot_dir` to verify the QC is working properly!

### 2. Adjust threshold if needed
If too many good points removed, increase threshold:
```bash
python analysis_scripts/apply_QC.py \
    configs/YOUR_CONFIG.py \
    --compass-qc \
    --compass-zscore-threshold 1.0 \
    --plot-compass-qc \
    --output-QC-tags qc
```

If outliers remain, decrease threshold:
```bash
python analysis_scripts/apply_QC.py \
    configs/YOUR_CONFIG.py \
    --compass-qc \
    --compass-zscore-threshold 0.3 \
    --plot-compass-qc \
    --output-QC-tags qc
```

### 3. Production run (no plots, save with tag)
```bash
python analysis_scripts/apply_QC.py \
    configs/YOUR_CONFIG.py \
    --compass-qc \
    --output-QC-tags qc \
    --output-file-tag _compass_qc
```

## Config Requirements

Your config file must have `conv_filenames_nc`:

```python
PIPS_IO_dict = {
    'PIPS_dir': '/path/to/netcdf/',
    'plot_dir': '/path/to/plots/',
    'PIPS_filenames_nc': ['parsivel_combined_*.nc', ...],
    'conv_filenames_nc': ['conventional_raw_*.nc', ...],  # ← Required!
    # ... other settings
}
```

## What gets updated

- **Conventional files**: compass_dir (outliers→NaN), winddirabs (recomputed)
- **Parsivel files**: All resampled wind variables updated with corrected compass
- **Plots** (if requested): `compass_wind_QC_{deployment}_{PIPS}.png`

## Recommended Workflow

1. Run with `--plot-compass-qc` first
2. Check plots to verify QC quality
3. Adjust threshold if needed
4. Run production version with `--output-file-tag`
5. Use QC'd files for downstream analysis

## Typical threshold values

- **0.3**: Aggressive (removes more points) - for severe outlier issues
- **0.5**: Default - good for PIPS3A/3B typical issues
- **1.0**: Conservative (removes fewer points) - for relatively clean data

## More info

See `docs/COMPASS_QC_USAGE.md` for detailed documentation.
