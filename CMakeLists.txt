# Molecular Orbital PACkage (MOPAC)
# Copyright (C) 2021, Virginia Polytechnic Institute and State University
#
# MOPAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MOPAC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Using Cmake version 3.11 for better CUDA & Fortran support
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

# Code coverage option
option(ENABLE_COVERAGE "Code coverage flag" OFF)

# Set a safe default build type
if(ENABLE_COVERAGE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build.")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build.")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

# Specify project name & programming language
project(openmopac LANGUAGES Fortran)

# location for CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Define a common library and executables for MOPAC & PARAM
add_library(mopac_core SHARED)
set_target_properties(mopac_core PROPERTIES OUTPUT_NAME "mopac")
add_executable(mopac)
add_executable(param)
target_link_libraries(mopac mopac_core)
target_link_libraries(param mopac_core)

# WinMOPAC build on Windows w/ the Intel Fortran compiler
option(BUILD_WINMOPAC "Build MOPAC with QuickWin GUI (ifort required)" OFF)
if(BUILD_WINMOPAC)
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"))
    message(FATAL_ERROR "WinMOPAC compilation requires the Intel Fortran compiler (ifort)")
  endif()
  if(NOT WIN32)
    message(FATAL_ERROR "WinMOPAC is only available on Windows")
  endif()
  add_executable(winmopac)
  target_link_libraries(winmopac mopac_core)
  target_compile_options(winmopac PRIVATE /libs:qwin)
  target_link_options(winmopac PRIVATE /subsystem:windows)
endif()

# MOPAC shared library ABI compatibility version
set_target_properties(mopac_core PROPERTIES SOVERSION 1)

# Install the executables and library
option(LIB_SUFFIX "Suffix for library directory (lib)" OFF)
install(TARGETS mopac DESTINATION bin)
install(TARGETS param DESTINATION bin)
if(LIB_SUFFIX)
  install(TARGETS mopac_core DESTINATION lib${LIB_SUFFIX})
else()
  install(TARGETS mopac_core DESTINATION lib)
endif()

# Pass OS information to MOPAC's Fortran source using preprocessor variables
if(APPLE)
  add_definitions(-DMOPAC_OS='Mac')
elseif(WIN32)
  add_definitions(-DMOPAC_OS='Windows')
elseif(UNIX)
  add_definitions(-DMOPAC_OS='Linux')
endif()

# Enable C language support for MKL?  (C is only needed for auto-detection of MKL by find_package)
option(MKL "Turn on C language support for MKL" ON)
if(MKL)
  enable_language(C)
endif()

# Try to use CMake's system of finding BLAS & LAPACK libraries
option(AUTO_BLAS "Use find_package to detect BLAS & LAPACK libraries" ON)
if(AUTO_BLAS)
  find_package(BLAS)
  find_package(LAPACK)

  if((NOT BLAS_FOUND) OR (NOT LAPACK_FOUND))
    message(FATAL_ERROR "find_package failed to find BLAS & LAPACK libraries, adjust environment variables & system paths or set AUTO_BLAS=OFF & manually define a link line with MOPAC_LINK")
  endif()

  # Use BLAS/LAPACK information provided by CMake's find_package system
  target_link_libraries(mopac_core PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  target_compile_options(mopac PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
  target_compile_options(param PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
else()
  message(STATUS "AUTO_BLAS=OFF, the user is responsible for linking BLAS & LAPACK libraries (e.g. define compiler link line with MOPAC_LINK)")
endif()

# Additional link line
set(MOPAC_LINK "" CACHE STRING "Additional link line for MOPAC executable")
separate_arguments(MOPAC_LINK)
target_link_libraries(mopac_core PUBLIC ${MOPAC_LINK})

# Add list of source files from src directory
add_subdirectory(src)

# Set up code coverage
set(GCOV_PROGRAM_NAME gcov CACHE STRING "Name of GNU code coverage analysis tool.")
set(LCOV_PROGRAM_NAME lcov CACHE STRING "Name of GNU code coverage visualization tool.")
if(ENABLE_COVERAGE)
  if(NOT ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU"))
    message(FATAL_ERROR "Code coverage requires the GNU Fortran compiler! Aborting...")
  endif()
  target_compile_options(mopac_core PUBLIC -fprofile-arcs -ftest-coverage -O0)
  target_link_options(mopac_core PUBLIC -fprofile-arcs -ftest-coverage)
  find_program(GCOV_PATH ${GCOV_PROGRAM_NAME})
  find_program(LCOV_PATH ${LCOV_PROGRAM_NAME})
  if(NOT GCOV_PATH)
    message(FATAL_ERROR "${GCOV_PROGRAM_NAME} not found! Aborting...")
  endif()
  if(NOT LCOV_PATH)
    message(FATAL_ERROR "${LCOV_PROGRAM_NAME} not found! Aborting...")
  endif()
  configure_file("cmake/CTestCustom.cmake" "${CMAKE_BINARY_DIR}/CTestCustom.cmake")
endif()

# MDI functionality (in development)
option(MDI "MDI build flag" OFF)
if(MDI)
  include(MDI)
endif()

# GPU functionality (not functional at the moment)
option(GPU "GPU build flag" OFF)
if(GPU)
  add_definitions(-DGPU)
endif()

# Add tests
enable_testing()
add_subdirectory(tests)

# some old stuff that should probably be deleted or cleaned up:

# compiler-specific debug flags
#if(${CMAKE_Fortran_COMPILER} MATCHES "pgfortran")
#  target_compile_options(mopac PRIVATE -g -Mbackslash)
#elseif(${CMAKE_Fortran_COMPILER} MATCHES "gfortran")
#    target_compile_options(mopac PRIVATE -g -O -Wuninitialized)
#    target_compile_options(mopac PRIVATE --coverage)
#  target_compile_options(mopac-lib PUBLIC "-fprofile-arcs" "-ftest-coverage")
#  target_link_libraries(mopac-lib PUBLIC gcov)
#     target_compile_options(mopac-all PUBLIC -g -Og -static -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -std=gnu -Wall -Wextra)
#     target_compile_options(mopac PUBLIC -g -Og -static -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -std=gnu -Wall -Wextra)
#     target_compile_options(param PUBLIC -g -Og -static -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -std=gnu -Wall -Wextra)
#    target_compile_options(mopac-lib PUBLIC -g -static -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -std=gnu -Wall -Wextra -Wunused-parameter -Waliasing -Wampersand -Wcharacter-truncation -Wreal-q-constant)
#    target_compile_options(mopac-lib PUBLIC -g -static -finit-local-zero -ffpe-trap=invalid,zero,overflow -fbacktrace -std=gnu -Wall -Wunused-parameter -fcheck=all -Waliasing -Wampersand -Wcharacter-truncation -Wreal-q-constant)
#    target_compile_options(mopac-lib PRIVATE -g -static -finit-local-zero -ffpe-trap=invalid,zero,overflow -fbacktrace -std=gnu -Wall -Wunused-parameter -fcheck=all -Waliasing -Wampersand -Wcharacter-truncation -Wreal-q-constant)
#elseif(${CMAKE_Fortran_COMPILER} MATCHES "ifort")
#     target_compile_options(mopac-all PUBLIC -g -traceback -check all -debug all -fpe0 -warn all)
#    target_compile_options(mopac PRIVATE -g -traceback -check all -debug all -fpe0 -warn all)
#    target_compile_options(mopac PRIVATE -g -traceback -check all -debug all -fpe0 -warn all -init=zero,arrays)
#else()
#  message("compiler not recognized, using default compiler options")
#endif()

#target_link_libraries(mopac PUBLIC -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)

#target_link_libraries(mopac PUBLIC $ENV{MKLROOT}/lib/libmkl_intel_ilp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a 
#                      $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
# ^^^ with -liomp5, which is causing linking problems
# THIS IS TOO IFORT SPECIFIC: target_link_libraries(mopac PUBLIC -mkl=parallel -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.3.185/mac/compiler/lib)
# ^^^ embeds path to dynamic OpenMP library

