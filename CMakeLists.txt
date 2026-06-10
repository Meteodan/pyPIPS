cmake_minimum_required(VERSION 3.17.2)

project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

# Find required packages
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(PkgConfig REQUIRED)

# Compiler detection and HPC support options
option(FORCE_F2PY_AUTODETECT "Force f2py to auto-detect compilers" OFF)
option(VERBOSE_F2PY "Enable verbose f2py compilation" OFF)

message(STATUS "Configuring ${SKBUILD_PROJECT_NAME} Fortran compilation")

# Flexible compiler detection with HPC support
# Priority: Environment variables > CMake detection > Common locations > HPC wrappers

# 1. Check environment variables first (highest priority)
if(DEFINED ENV{FC})
    set(CMAKE_Fortran_COMPILER $ENV{FC})
    message(STATUS "Using Fortran compiler from FC environment variable: ${CMAKE_Fortran_COMPILER}")
elseif(DEFINED ENV{F90})
    set(CMAKE_Fortran_COMPILER $ENV{F90})
    message(STATUS "Using Fortran compiler from F90 environment variable: ${CMAKE_Fortran_COMPILER}")
endif()

if(DEFINED ENV{CC})
    set(CMAKE_C_COMPILER $ENV{CC})
    message(STATUS "Using C compiler from CC environment variable: ${CMAKE_C_COMPILER}")
endif()

# 2. CMake compiler detection (if not set by environment)
if(NOT CMAKE_Fortran_COMPILER)
    # Let CMake find compilers in standard locations
    find_program(GFORTRAN_EXECUTABLE gfortran PATHS /usr/bin /opt/homebrew/bin /usr/local/bin)
    find_program(IFORT_EXECUTABLE ifort)

    if(GFORTRAN_EXECUTABLE)
        set(CMAKE_Fortran_COMPILER ${GFORTRAN_EXECUTABLE})
        message(STATUS "Found gfortran: ${CMAKE_Fortran_COMPILER}")
    elseif(IFORT_EXECUTABLE)
        set(CMAKE_Fortran_COMPILER ${IFORT_EXECUTABLE})
        message(STATUS "Found ifort: ${CMAKE_Fortran_COMPILER}")
    endif()
endif()

# 3. HPC wrapper detection (lowest priority, but important for HPC systems)
if(NOT CMAKE_Fortran_COMPILER)
    message(STATUS "Checking for HPC compiler wrappers...")
    foreach(compiler IN ITEMS ftn mpifort mpif90 pgfortran pgf90)
        find_program(${compiler}_EXECUTABLE ${compiler})
        if(${compiler}_EXECUTABLE)
            set(CMAKE_Fortran_COMPILER ${${compiler}_EXECUTABLE})
            message(STATUS "Found HPC wrapper: ${CMAKE_Fortran_COMPILER}")
            break()
        endif()
    endforeach()
endif()

# Similar detection for C compiler
if(NOT CMAKE_C_COMPILER)
    find_program(GCC_EXECUTABLE gcc PATHS /usr/bin /opt/homebrew/bin /usr/local/bin)
    find_program(CLANG_EXECUTABLE clang PATHS /usr/bin /opt/homebrew/bin /usr/local/bin)
    find_program(ICC_EXECUTABLE icc)

    if(GCC_EXECUTABLE)
        set(CMAKE_C_COMPILER ${GCC_EXECUTABLE})
        message(STATUS "Found gcc: ${CMAKE_C_COMPILER}")
    elseif(CLANG_EXECUTABLE)
        set(CMAKE_C_COMPILER ${CLANG_EXECUTABLE})
        message(STATUS "Found clang: ${CMAKE_C_COMPILER}")
    elseif(ICC_EXECUTABLE)
        set(CMAKE_C_COMPILER ${ICC_EXECUTABLE})
        message(STATUS "Found icc: ${CMAKE_C_COMPILER}")
    endif()
endif()

# HPC C compiler wrappers
if(NOT CMAKE_C_COMPILER)
    foreach(compiler IN ITEMS cc mpicc pgcc)
        find_program(${compiler}_EXECUTABLE ${compiler})
        if(${compiler}_EXECUTABLE)
            set(CMAKE_C_COMPILER ${${compiler}_EXECUTABLE})
            message(STATUS "Found HPC C wrapper: ${CMAKE_C_COMPILER}")
            break()
        endif()
    endforeach()
endif()

# Enable languages after compiler detection
enable_language(Fortran)

# Get Python extension suffix dynamically for cross-platform compatibility
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Python extension suffix: ${PYTHON_EXT_SUFFIX}")

# Function to compile Fortran modules with f2py
function(compile_fortran_module module_name source_files)
    set(module_target "${module_name}${PYTHON_EXT_SUFFIX}")

    message(STATUS "Configuring compilation of Fortran module: ${module_name}")
    message(STATUS "Source files: ${source_files}")

    # Build f2py command
    set(f2py_cmd ${Python_EXECUTABLE} -m numpy.f2py)

    # Add source files and module compilation flags - let f2py auto-detect compilers
    list(APPEND f2py_cmd --backend=distutils -c ${source_files} -m ${module_name})

    # Create custom command to compile the Fortran module
    add_custom_command(
        OUTPUT ${module_target}
        COMMAND ${f2py_cmd}
        DEPENDS ${source_files}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Compiling Fortran module ${module_name} with f2py"
        VERBATIM
    )

    # Create target that depends on the compiled module
    add_custom_target(${module_name}_target ALL DEPENDS ${module_target})

    # Install the compiled module to the correct package location
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${module_target}
            DESTINATION ${SKBUILD_PROJECT_NAME})

    message(STATUS "Module ${module_name} will be installed to: ${SKBUILD_PROJECT_NAME}/")
endfunction()

# Compile the dualpara Fortran module
# This is the same module as in pyCRMtools, containing dual-polarization radar calculations
message(STATUS "Setting up dualpara module compilation...")

set(FORTRAN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/pyPIPS/global_module.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/pyPIPS/dualpara.f90"
)

# Verify source files exist
foreach(source_file ${FORTRAN_SOURCES})
    if(NOT EXISTS ${source_file})
        message(FATAL_ERROR "Fortran source file not found: ${source_file}")
    endif()
    message(STATUS "Found Fortran source: ${source_file}")
endforeach()

# Compile the dualpara module
compile_fortran_module(dualpara "${FORTRAN_SOURCES}")

message(STATUS "pyPIPS Fortran module configuration complete")