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
target_link_libraries(mopac PUBLIC mopac-all)
target_link_libraries(param PUBLIC mopac-all)

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

# Link w/ MKL if it is available
if(DEFINED ENV{MKLROOT})
  set(CMAKE_BUILD_TYPE Release)
  message("environment variable MKLROOT is defined, assuming that Intel MKL is available")
  if(DEFINED ENV{OMPROOT})
    message("environment variable OMPROOT is defined, building threaded MKL w/ OMPROOT as the path to libiomp5")
    target_link_libraries(mopac PUBLIC -Wl,-rpath,$ENV{OMPROOT}
                          $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a 
                          $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{OMPROOT} -liomp5 -lpthread -lm -ldl)
    target_link_libraries(param PUBLIC -Wl,-rpath,$ENV{OMPROOT}
                          $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a 
                          $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{OMPROOT} -liomp5 -lpthread -lm -ldl)
#                          -dynamiclib -lgfortran -static-libgfortran
#                          /usr/local/opt/gcc/lib/gcc/8/libquadmath.a -static-libgcc)
  else()
    message("environment variable OMPROOT is not defined, building serial MKL")
    target_link_libraries(mopac PUBLIC $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_sequential.a
                          $ENV{MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
  endif()
# Otherwise try generic BLAS/LAPACK libraries
else()
  message("environment variable MKLROOT is not defined, relying on CMake to find working BLAS & LAPACK libraries")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(mopac PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  target_link_libraries(param PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

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

# Force static linking of GNU libraries by request
option(STATIC_GNU "STATIC_GNU" OFF)
if(STATIC_GNU)
  find_library(LIBQUADMATH libquadmath.a /usr/local/Cellar/gcc/10.2.0/lib/gcc/10 REQUIRED)
  target_link_libraries(mopac PUBLIC -lgfortran ${LIBQUADMATH} -static-libgfortran -static-libgcc)
  message("LIBQUADMATH = " ${LIBQUADMATH})
endif()

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
