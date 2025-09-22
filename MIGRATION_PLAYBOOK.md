# pyCRMtools Modernization Playbook
## A Step-by-Step Guide for Applying to Other Projects

This document captures the complete modernization process we applied to pyCRMtools, so it can be systematically applied to other similar projects.

## Overview of What We Accomplished

### 1. **Replaced Deprecated numpy.distutils with scikit-build-core**
- **Problem**: numpy.distutils deprecated, setup.py overly complex
- **Solution**: Modern pyproject.toml + scikit-build-core + CMake
- **Result**: Single `pip install -e .` command handles everything

### 2. **Implemented Flexible Compiler Detection**
- **Problem**: Hardcoded compiler paths, HPC incompatibility
- **Solution**: Environment variable detection + HPC wrapper support
- **Result**: Works on local dev, HPC clusters, different compilers

### 3. **Solved Conda/Pip Dependency Conflicts**
- **Problem**: Pip overwriting conda-managed packages
- **Solution**: Minimal dependencies + `--no-deps` installation strategy
- **Result**: Clean installations respecting existing conda environments

## Step-by-Step Migration Process

### Phase 1: Analysis and Preparation

#### 1.1 Assess Current Project Structure
```bash
# Document current state
ls -la /path/to/project/
find . -name "*.f90" -o -name "*.f" | head -10
cat setup.py  # if exists
cat requirements.txt  # if exists
```

**Key Questions:**
- What Fortran files need compilation?
- What's the current build system?
- What are the dependencies?
- Any legacy build artifacts?

#### 1.2 Backup Current State
```bash
git checkout -b modernization-backup
git add . && git commit -m "Backup before modernization"
git checkout -b modernize-packaging
```

### Phase 2: Create Modern Build System

#### 2.1 Create pyproject.toml
**Template based on pyCRMtools:**

```toml
[build-system]
requires = [
    "scikit-build-core",
    "numpy>=1.21.0",
]
build-backend = "scikit_build_core.build"

[project]
name = "YOUR_PROJECT_NAME"
version = "0.1.0"  # or current version
description = "YOUR_DESCRIPTION"
readme = "README.rst"  # or README.md
license = {text = "YOUR_LICENSE"}
authors = [
    {name = "YOUR_NAME", email = "your.email@domain.com"},
]
requires-python = ">=3.8"
dependencies = [
    # Only essential deps to avoid conda conflicts
    "numpy>=1.21.0",
]

[project.optional-dependencies]
scientific = [
    # Add your scientific dependencies here
    "scipy",
    "matplotlib", 
    "pandas",
    # etc.
]

[tool.scikit-build]
cmake.version = ">=3.17.2"
cmake.args = ["-DCMAKE_BUILD_TYPE=Release"]
build.requires = ["numpy"]
wheel.packages = ["YOUR_PACKAGE_NAME"]
```

**Customization needed:**
- Replace ALL_CAPS placeholders
- Adjust dependencies based on your project needs
- Update package name and paths

#### 2.2 Create CMakeLists.txt
**Template for Fortran compilation:**

```cmake
cmake_minimum_required(VERSION 3.17.2)

project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

# Find required packages
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(PkgConfig REQUIRED)

# Compiler detection with HPC support
option(FORCE_F2PY_AUTODETECT "Force f2py to auto-detect compilers" OFF)
option(VERBOSE_F2PY "Enable verbose f2py compilation" OFF)

# Flexible compiler detection
if(DEFINED ENV{FC})
    set(CMAKE_Fortran_COMPILER $ENV{FC})
    message(STATUS "Using Fortran compiler from FC environment variable: ${CMAKE_Fortran_COMPILER}")
elseif(DEFINED ENV{F90})
    set(CMAKE_Fortran_COMPILER $ENV{F90})
    message(STATUS "Using Fortran compiler from F90 environment variable: ${CMAKE_Fortran_COMPILER}")
endif()

# HPC wrapper detection
if(NOT CMAKE_Fortran_COMPILER)
    foreach(compiler IN ITEMS ftn mpif90 mpifort pgfortran gfortran ifort)
        find_program(${compiler}_EXECUTABLE ${compiler})
        if(${compiler}_EXECUTABLE)
            set(CMAKE_Fortran_COMPILER ${${compiler}_EXECUTABLE})
            message(STATUS "Found HPC/system compiler: ${CMAKE_Fortran_COMPILER}")
            break()
        endif()
    endforeach()
endif()

# Enable languages after compiler detection
enable_language(Fortran)

# Get Python extension suffix dynamically
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Function to compile Fortran modules
function(compile_fortran_module module_name source_files)
    set(module_target "${module_name}${PYTHON_EXT_SUFFIX}")
    
    # Build f2py command
    set(f2py_cmd ${Python_EXECUTABLE} -m numpy.f2py)
    
    # Add compiler flags if not auto-detecting
    if(NOT FORCE_F2PY_AUTODETECT)
        if(CMAKE_Fortran_COMPILER)
            list(APPEND f2py_cmd --f90exec=${CMAKE_Fortran_COMPILER})
        endif()
        if(CMAKE_C_COMPILER)
            list(APPEND f2py_cmd --c_compiler=${CMAKE_C_COMPILER})
        endif()
    endif()
    
    # Add verbose flag if requested
    if(VERBOSE_F2PY)
        list(APPEND f2py_cmd --verbose)
    endif()
    
    # Add source files and compilation flags
    list(APPEND f2py_cmd -c ${source_files} -m ${module_name})
    
    add_custom_command(
        OUTPUT ${module_target}
        COMMAND ${f2py_cmd}
        DEPENDS ${source_files}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Compiling Fortran module ${module_name}"
    )
    
    add_custom_target(${module_name}_target ALL DEPENDS ${module_target})
    
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${module_target}
            DESTINATION ${SKBUILD_PROJECT_NAME}/modules)
endfunction()

# Compile your Fortran modules
# Example: compile_fortran_module(your_module "path/to/source.f90;path/to/other.f90")
compile_fortran_module(your_module_name "${CMAKE_CURRENT_SOURCE_DIR}/path/to/your/fortran/files.f90")
```

**Customization needed:**
- Replace `your_module_name` with actual module name
- Update paths to your Fortran source files
- Add multiple `compile_fortran_module()` calls if you have multiple modules

#### 2.3 Create Documentation Files
Copy and customize these files:
- `INSTALLATION_GUIDE.md`
- `HPC_COMPILER_GUIDE.md`
- Update `README.rst` with installation instructions

### Phase 3: Testing and Validation

#### 3.1 Clean Up Old Build System
```bash
# Remove old build artifacts (adapt to your project)
rm -rf build/ dist/ *.egg-info/
rm -f setup.py setup.cfg MANIFEST.in  # if they exist
# Remove any old .so files
find . -name "*.so" -delete
```

#### 3.2 Test New Build System
```bash
# Test compilation
pip install -e . --no-deps -v

# Test imports
python -c "import your_package; print('Success!')"
python -c "from your_package.modules import your_module; print('Fortran module works!')"
```

#### 3.3 Test Different Compiler Configurations
```bash
# Test environment variable detection
FC=gfortran pip install -e . --force-reinstall --no-deps

# Test HPC compatibility (if available)
export FC=ifort  # or ftn, mpif90, etc.
pip install -e . --force-reinstall --no-deps
```

### Phase 4: Environment Compatibility

#### 4.1 Test Conda Environment Installation
```bash
# In conda environment with existing scientific packages
conda list | grep -E "(numpy|scipy|matplotlib)"
pip install -e . --no-deps  # Should not disturb conda packages
```

#### 4.2 Test Pure Pip Environment
```bash
# In fresh virtual environment
python -m venv test_env
source test_env/bin/activate
pip install -e .[scientific]
```

## Key Decision Points

### 1. **Dependency Strategy**
- **Minimal approach**: Only numpy in main dependencies
- **Full approach**: Include all scientific dependencies
- **Recommendation**: Use minimal + optional dependencies for conda compatibility

### 2. **Compiler Detection Level**
- **Basic**: Just environment variables
- **Advanced**: HPC wrapper detection + fallbacks
- **Recommendation**: Include HPC support for maximum flexibility

### 3. **Documentation Completeness**
- **Minimal**: Just update README
- **Complete**: Full installation guides + HPC instructions
- **Recommendation**: Complete documentation prevents future confusion

## Common Gotchas and Solutions

### 1. **Module Path Issues**
- **Problem**: Python can't find compiled modules
- **Solution**: Ensure CMakeLists.txt installs to correct package directory
- **Check**: Verify `install(FILES ... DESTINATION ${SKBUILD_PROJECT_NAME}/modules)`

### 2. **Compiler Detection Failures**
- **Problem**: CMake can't find compilers
- **Solution**: Add more fallback options in CMakeLists.txt
- **Debug**: Use `pip install -e . -v` for verbose output

### 3. **Conda/Pip Conflicts**
- **Problem**: Pip overwrites conda packages
- **Solution**: Always use `--no-deps` in conda environments
- **Prevention**: Keep dependencies minimal in pyproject.toml

### 4. **Legacy Artifacts**
- **Problem**: Old build files interfere
- **Solution**: Thorough cleanup before testing
- **Check**: Look for `.so` files, `build/` directories, `.egg-info/`

## Success Criteria

✅ **Build System Works**
- Single `pip install -e .` command succeeds
- No manual compilation steps required
- Works with `--no-deps` flag

✅ **Cross-Platform Compatibility**
- Works on local development machine
- Compatible with HPC systems
- Handles different compiler configurations

✅ **Environment Respect**
- Doesn't disturb conda-managed packages
- Works in both conda and pip environments
- No version conflicts

✅ **Maintainability**
- Clear documentation for future developers
- Modern, standard packaging approach
- Easy to modify and extend

## Next Steps After Migration

1. **Update CI/CD** (if you have it) to use new build system
2. **Test on target HPC systems** to verify compatibility
3. **Update contributor documentation** with new build process
4. **Consider creating release workflow** with modern packaging

## Files to Reference from pyCRMtools

When migrating, refer to these files from the pyCRMtools project:
- `pyproject.toml` - Modern packaging configuration
- `CMakeLists.txt` - Flexible Fortran compilation
- `INSTALLATION_GUIDE.md` - Comprehensive user instructions
- `HPC_COMPILER_GUIDE.md` - HPC-specific guidance
- `README.rst` - Updated installation section