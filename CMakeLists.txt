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
cmake_minimum_required(VERSION 3.11)

# location for CMake modules (MDI)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Specify project name & programming language (C is only needed for auto-detection of MKL by find_package)
project(openmopac LANGUAGES Fortran C)

# Define a combined library and executables for MOPAC & PARAM
add_library(mopac-all)
add_executable(mopac)
add_executable(param)
target_link_libraries(mopac mopac-all)
target_link_libraries(param mopac-all)

# Pass OS information to MOPAC's Fortran source using preprocessor variables
if(APPLE)
  add_definitions(-DMOPAC_OS='M')
elseif(UNIX)
  add_definitions(-DMOPAC_OS='L')
elseif(WIN32 OR MINGW OR MSYS OR CYGWIN)
  add_definitions(-DMOPAC_OS='W')
else()
  add_definitions(-DMOPAC_OS='X')
endif()

# Try to use CMake's system of finding BLAS & LAPACK libraries
find_package(BLAS)
find_package(LAPACK)

  set(BLAS_LIBRARIES "")
  set(LAPACK_LIBRARIES "")
  set(BLAS_LINKER_FLAGS "")
  set(LAPACK_LINKER_FLAGS "")

# Attempt to manually set libraries if find_package fails
if((NOT BLAS_FOUND) OR (NOT LAPACK_FOUND) AND 0)
  message(WARNING "find_package could not find BLAS or LAPACK libraries, attempting to find Intel MKL with find_library (setting MKLROOT might help)")
  if(DEFINED ENV{MKLROOT})
    message(STATUS "MKLROOT = " $ENV{MKLROOT})
    file(TO_CMAKE_PATH "$ENV{MKLROOT}/lib" MKL_PATH1)
    file(TO_CMAKE_PATH "$ENV{MKLROOT}/bin" MKL_PATH2)
    file(TO_CMAKE_PATH "$ENV{MKLROOT}/../compiler/lib" MKL_PATH3)
  endif()
  if(BLA_STATIC)
    find_library(INTEL_LP64 NAMES libmkl_intel_lp64.a mkl_intel_lp64.lib PATHS ${MKL_PATH1})
    find_library(INTEL_THREAD NAMES libmkl_intel_thread.a mkl_intel_thread.lib PATHS ${MKL_PATH1})
    find_library(INTEL_CORE NAMES libmkl_core.a mkl_core.lib PATHS ${MKL_PATH1})
  else()
    find_library(INTEL_LP64 mkl_intel_lp64 PATHS ${MKL_PATH1} ${MKL_PATH2})
    find_library(INTEL_THREAD mkl_intel_thread PATHS ${MKL_PATH1} ${MKL_PATH2})
    find_library(INTEL_CORE mkl_core PATHS ${MKL_PATH1} ${MKL_PATH2})
  endif()
  message("path check:" ${MKL_PATH1} ${MKL_PATH2} ${MKL_PATH3})
  find_library(INTEL_OPENMP NAMES iomp5 libiomp5md.dll PATHS ${MKL_PATH1} ${MKL_PATH2} ${MKL_PATH3})
  if(NOT INTEL_OPENMP)
    message("attempting to hard-code iomp for Windows")
    set(INTEL_OPENMP "${MKL_PATH2}/libiomp5md.dll")
  endif()
  set(BLAS_LIBRARIES
      ${INTEL_LP64}
      ${INTEL_THREAD}
      ${INTEL_CORE}
      ${INTEL_OPENMP})
  set(LAPACK_LIBRARIES "")
  set(BLAS_LINKER_FLAGS "")
  set(LAPACK_LINKER_FLAGS "")
  if(APPLE OR UNIX)
    find_library(LIB_PTHREAD pthread)
    find_library(LIB_MATH m)
    find_library(LIB_DL dl)
    list(APPEND BLAS_LIBRARIES ${LIB_PTHREAD} ${LIB_MATH} ${LIB_DL})
  endif()
endif()

# Use BLAS/LAPACK information provided by CMake's find_package system
target_link_libraries(mopac ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(param ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_compile_options(mopac PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
target_compile_options(param PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})

# Add list of source files from src directory
add_subdirectory(src)

# MDI functionality
option(MDI "MDI build flag" OFF)
if(MDI)
  include(MDI)
endif(MDI)

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

