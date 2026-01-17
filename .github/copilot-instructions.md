# pyPIPS AI Coding Assistant Instructions

## Project Overview
pyPIPS analyzes precipitation data from Portable In-situ Precipitation Stations (PIPS) - disdrometer networks collecting Drop Size Distribution (DSD) measurements. This is a scientific computing project combining Python, Fortran, and real-time data processing.

## Architecture & Key Components

### Core Data Flow
1. **Raw Data**: Parsivel disdrometer text files → `pips_io.combine_parsivel_data()`
2. **Processing**: DSD matrices → `DSDlib.py` calculations → derived parameters
3. **Quality Control**: `parsivel_qc.py` applies multiple QC filters
4. **Output**: NetCDF files with xarray DataArrays for analysis/visualization

### Critical Module Hierarchy
- **`pyPIPS/PIPS.py`**: Core data structures, thermodynamic calculations, wind processing
- **`pyPIPS/DSDlib.py`**: Drop size distribution computations (1900+ lines) - gamma fitting, radar moments, microphysics
- **`pyPIPS/parsivel_params.py`**: Instrument specifications (diameter/velocity bins, sensor area)
- **`pyPIPS/pips_io.py`**: Data ingestion from raw Parsivel files to pandas/xarray
- **`pyPIPS/parsivel_qc.py`**: Quality control filters (wind, splashing, margin effects)

### Build System (Hybrid Python/Fortran)
- Uses **scikit-build-core** + CMake for Fortran compilation via f2py
- HPC-aware compiler detection: environment variables (FC, CC) → CMake → common paths → HPC wrappers
- Key Fortran modules: `dualpara.f90`, `global_module.f90` (compiled into Python extensions)

## Essential Workflows

### Configuration-Driven Analysis
- **Configs pattern**: `configs/*.py` files define dataset-specific parameters via `PIPS_IO_dict`
- Example: `PIPS_IO_dict['PIPS_filenames']` lists input files, `PIPS_IO_dict['deployment_names']` organizes by field campaign
- QC settings in `PIPS_qc_dict` (strongwindQC, splashingQC, marginQC, etc.)

### Analysis Script Pattern
```python
# Standard imports for analysis scripts
import pyPIPS.PIPS as pips
import pyPIPS.DSDlib as dsd
import pyPIPS.parsivel_qc as pqc
import pyPIPS.pips_io as pipsio
```
- Scripts use `argparse` with `case_config_path` pointing to configs
- Common workflow: load config → apply QC → calculate derived parameters → save NetCDF

### Data Conventions
- **Diameter bins**: `avg_diameter` (mm), `min_diameter`/`max_diameter` arrays from parsivel_parameters
- **Velocity bins**: `avg_fall_bins` (m/s), `min_fall_bins`/`max_fall_bins`
- **Variable naming**: Use suffixes like `_qc`, `_RB15_qc` for QC'd versions; tags like `input_QC_tag = '_{}'.format(args.input_QC_tag)`

## Development Practices

### Environment Management
- **Conda (recommended)**: `conda env create -f environment.yml` + `pip install -e . --no-deps` (avoids conflicts)
- **Pure pip**: `pip install -e .[scientific]` includes all scientific dependencies
- **Key dependencies**: numpy, xarray, pandas, numba (JIT compilation), cartopy, metpy

### Data Processing Patterns
- **Time resampling**: Uses pandas `resample()` with `offset` strings, `label='right'`, `closed='right'`
- **Wind calculations**: `wind_dir_and_speed_from_u_and_v()` handles meteorological conventions (270° - atan2())
- **Gamma distribution fitting**: `DSDlib.py` uses `scipy.special.gamma` extensively for DSD moment calculations

### File Organization
- **analysis_scripts/**: Production analysis workflows (apply_QC.py, calc_derived_params.py, etc.)
- **plotting_scripts/**: Visualization utilities using `plotmodule.py`
- **realtime_scripts/**: Live data ingestion (`pips_realtime.py` with urllib3/BeautifulSoup)
- **shell_scripts/**: Automation workflows (e.g., `tripips.sh` for data scraping)

### Plotting Conventions
- Uses `plotmodule.py` with MetPy colortables: `ctables.registry.get_with_steps('NWSReflectivity')`
- Contour levels: `clevels_ref = np.arange(0.0, 85.0, 5.0)` for reflectivity, custom ZDR/velocity ranges
- Font settings: `mpl.rc('font', size=9)` globally set

## Integration Points

### Real-time Data Ingestion
- HTTP scraping from PIPS servers: `scrape_onesec_data()` uses urllib3 → BeautifulSoup → pandas
- GPS coordinate handling with fallback to "last good" positions
- Data streams processed in 1-second intervals, resampled to analysis intervals (typically 60s)

### Cross-component Communication
- Config files bridge deployment metadata with processing scripts
- NetCDF serves as interchange format between processing stages
- QC tags propagate through analysis chain to maintain data provenance

When implementing new features, follow the config-driven pattern, maintain xarray compatibility for scientific datasets, and preserve QC tag conventions for traceability.

## Progress Tracking

### Update Progress Command
When the user says "Update progress", update the `.github/progress.md` file with a new session entry documenting:

1. **Session metadata**: Date range and session theme/focus
2. **Key code updates**: Specific files modified, problems solved, and solutions implemented
3. **Lessons learned**: Technical insights, best practices discovered, and patterns established
4. **Technical context**: Architecture decisions, integration patterns, and development workflows
5. **Next priorities**: Suggested future development directions

**Progress Entry Format**:
```markdown
## Session: [Date Range] - [Session Theme]

### Key Code Updates
#### [Update Category]
**Problem**: [Brief problem description]
**Solution**: [Key solution points as bullets]
**Files Modified**: [List of files]

### Key Lessons Learned
#### [Lesson Category]
**Lesson**: [Specific lesson or insight]
**Impact**: [Why this matters for future development]

### Technical Context for Future Development
[Architecture notes, patterns, workflow context]

### Next Development Priorities
[Suggested next steps or improvements]
```

**Session Tracking Guidelines**:
- Document significant code changes, not minor edits
- Focus on lessons that inform future development decisions
- Include context that helps new agent sessions understand current state
- Update existing entries if work continues on same topics
- Maintain chronological order with most recent sessions at top