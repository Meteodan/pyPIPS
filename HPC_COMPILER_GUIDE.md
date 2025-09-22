# HPC Compiler Configuration Guide

This guide explains how the modern build system supports flexible compiler detection for HPC environments.

## Flexible Compiler Detection

The build system now supports multiple ways to specify compilers, making it compatible with various HPC systems:

### 1. Environment Variables (Highest Priority)
Set these environment variables before running `pip install`:

```bash
# Fortran compiler
export FC=gfortran
export F90=gfortran

# C compiler
export CC=gcc

# Examples for different HPC systems:
export FC=ifort          # Intel Fortran
export FC=ftn            # Cray wrapper
export FC=mpif90         # MPI wrapper
export FC=pgfortran      # PGI/NVIDIA
```

### 2. HPC Module Systems
The build system automatically detects common HPC compiler wrappers:

```bash
# After loading modules on HPC systems:
module load intel/2023.1
module load gcc/12.2.0
module load cray-ftn

# The build system will automatically find:
# - ftn (Cray wrapper)
# - mpif90, mpifort (MPI wrappers)
# - ifort (Intel)
# - pgfortran (PGI/NVIDIA)
```

### 3. CMake Options
You can also specify compilers via CMake options:

```bash
pip install -e . -C cmake.args="-DCMAKE_Fortran_COMPILER=ifort"
pip install -e . -C cmake.args="-DFORCE_F2PY_AUTODETECT=ON"
```

### 4. Available CMake Options

- `FORCE_F2PY_AUTODETECT=ON`: Let f2py automatically detect compilers
- `VERBOSE_F2PY=ON`: Enable verbose compilation output
- `CMAKE_Fortran_COMPILER`: Specify Fortran compiler directly
- `CMAKE_C_COMPILER`: Specify C compiler directly

## Installation Examples

### Local Development (macOS/Linux)
```bash
pip install -e .
```

### HPC with Intel Compilers
```bash
module load intel
export FC=ifort
export CC=icc
pip install -e .
```

### HPC with Cray Wrappers
```bash
module load PrgEnv-gnu
export FC=ftn
export CC=cc
pip install -e .
```

### HPC with MPI Wrappers
```bash
module load openmpi
export FC=mpif90
export CC=mpicc
pip install -e .
```

### Force Auto-Detection
```bash
pip install -e . -C cmake.args="-DFORCE_F2PY_AUTODETECT=ON"
```

## Compiler Detection Priority

1. **Environment Variables** (`$FC`, `$F90`, `$CC`)
2. **CMake Detection** (find_program for specified compiler)
3. **Common Locations** (/usr/bin, /opt/homebrew/bin, etc.)
4. **HPC Wrappers** (ftn, mpif90, mpifort, pgfortran)
5. **f2py Auto-Detection** (as fallback)

## Troubleshooting

### Verbose Output
Enable verbose compilation to see detected compilers:
```bash
pip install -e . -v -C cmake.args="-DVERBOSE_F2PY=ON"
```

### Check Detected Compilers
The CMake output will show which compilers were detected:
```
-- Using Fortran compiler: /path/to/gfortran
-- Using C compiler: /path/to/gcc
```

### Force Specific Compiler
If automatic detection fails:
```bash
export FC=/full/path/to/fortran/compiler
pip install -e .
```

## Compatibility

This build system is tested and compatible with:

- **Local Development**: gfortran, clang, gcc
- **Intel OneAPI**: ifort, icx
- **Cray Systems**: ftn, cc (wrappers)
- **MPI Systems**: mpif90, mpifort, mpicc
- **NVIDIA HPC SDK**: pgfortran, pgcc
- **GNU Compilers**: gfortran, gcc on various systems

The flexible detection ensures your package builds correctly across different computing environments without manual intervention.