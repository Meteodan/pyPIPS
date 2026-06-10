# pyCRMtools Installation Guide

This guide explains how to install pyCRMtools in different environments while avoiding dependency conflicts.

## Installation Options

### 1. Conda Environment (Recommended for most users)

If you're using conda/mamba with an environment that already has scientific packages:

```bash
# Install only pyCRMtools without dependencies (avoids conda conflicts)
pip install -e . --no-deps

# OR install with minimal dependencies
pip install -e .
```

**Benefits:**
- Preserves conda-managed packages
- No version conflicts
- Uses conda's optimized scientific stack

**Prerequisites:** Ensure your conda environment has:
```bash
conda install numpy scipy matplotlib pandas xarray netcdf4 cartopy
conda install -c conda-forge metpy pytz
```

### 2. Pure pip Environment

If you're in a pure pip environment or want pip to manage all dependencies:

```bash
# Install with scientific dependencies
pip install -e .[scientific]
```

### 3. Development Installation

For developers who need additional tools:

```bash
# Conda environment (recommended)
pip install -e .[dev] --no-deps

# Pure pip environment
pip install -e .[scientific,dev]
```

### 4. Minimal Installation

For users who only need the core functionality:

```bash
pip install -e .
```

This only installs numpy (required for Fortran compilation) and pyCRMtools itself.

## Dependency Management

### What's Included by Default
- `numpy>=1.21.0` (required for f2py compilation)

### Optional Dependencies

#### `[scientific]` - Scientific Computing Stack
- `scipy` - Scientific computing library
- `matplotlib` - Plotting library
- `pandas` - Data analysis library
- `xarray` - N-dimensional labeled arrays
- `netCDF4` - NetCDF file I/O
- `cartopy` - Cartographic projections
- `metpy` - Meteorological calculations
- `pytz` - Timezone handling

#### `[dev]` - Development Tools
- `autopep8`, `flake8`, `pylint` - Code formatting and linting
- `pytest`, `pytest-cov` - Testing framework
- `ipython`, `jupyter` - Interactive development
- `joblib` - Parallel computing

#### `[docs]` - Documentation Tools
- `sphinx`, `sphinx-rtd-theme`, `numpydoc` - Documentation generation

#### `[complete]` - Additional Optional Tools
- `sharppy` - Atmospheric sounding analysis

## Environment Setup Examples

### Scientific Conda Environment

```bash
# Create environment with scientific stack
conda create -n pycrmtools python=3.10
conda activate pycrmtools

# Install scientific dependencies via conda
conda install numpy scipy matplotlib pandas xarray netcdf4
conda install -c conda-forge cartopy metpy pytz

# Install pyCRMtools without dependency conflicts
cd /path/to/pyCRMtools
pip install -e . --no-deps
```

### Development Environment

```bash
# Create development environment
conda create -n pycrmtools-dev python=3.10
conda activate pycrmtools-dev

# Install scientific stack via conda
conda install numpy scipy matplotlib pandas xarray netcdf4
conda install -c conda-forge cartopy metpy pytz

# Install development tools
conda install jupyter ipython

# Install pyCRMtools with dev tools (excluding scientific deps)
cd /path/to/pyCRMtools
pip install -e .[dev] --no-deps
```

### Pure Pip Environment

```bash
# Create virtual environment
python -m venv pycrmtools-env
source pycrmtools-env/bin/activate  # Linux/macOS
# pycrmtools-env\Scripts\activate     # Windows

# Install everything via pip
cd /path/to/pyCRMtools
pip install -e .[scientific,dev]
```

## Avoiding Dependency Conflicts

### The Problem
Mixing conda and pip can cause:
- Version conflicts
- Conda packages overwritten by pip
- Broken dependency tracking

### The Solution
1. **Conda environments**: Use `--no-deps` with pip install
2. **Check existing packages**: `conda list` before installation
3. **Separate environments**: Keep conda and pip environments separate

### Checking Your Installation

```bash
# Check what's installed
conda list | grep -E "(numpy|scipy|matplotlib|pandas|xarray|netcdf4|cartopy|metpy)"

# Verify pyCRMtools installation
python -c "import pyCRMtools; print(pyCRMtools.__version__)"

# Test Fortran module
python -c "from pyCRMtools.modules import dualpara; print('Success!')"
```

## Troubleshooting

### Dependency Conflicts
If you see pip dependency conflicts:
```bash
# Uninstall and reinstall without dependencies
pip uninstall pyCRMtools
pip install -e . --no-deps
```

### Missing Dependencies
If you get import errors:
```bash
# Install missing packages via conda (preferred)
conda install package_name

# Or via pip if conda unavailable
pip install package_name
```

### Compilation Issues
If Fortran compilation fails:
```bash
# Check compiler detection
FC=gfortran pip install -e . -v --no-deps

# Force f2py autodetection
pip install -e . -C cmake.args="-DFORCE_F2PY_AUTODETECT=ON" --no-deps
```

## Best Practices

1. **Use conda for scientific packages** - Better optimized and tested
2. **Use pip for pure Python packages** - Easier to manage
3. **Keep environments separate** - Avoid mixing package managers
4. **Pin critical versions** - Use environment.yml for reproducibility
5. **Test after installation** - Verify imports and basic functionality

## Quick Reference

| Scenario | Command |
|----------|---------|
| Conda env with scientific packages | `pip install -e . --no-deps` |
| Conda env missing packages | `conda install <packages>` then `pip install -e . --no-deps` |
| Pure pip environment | `pip install -e .[scientific]` |
| Development setup | `pip install -e .[dev] --no-deps` |
| Minimal installation | `pip install -e .` |
| All optional features | `pip install -e .[scientific,dev,docs,complete]` |