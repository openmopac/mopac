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

# Specify project name & programming language
project(math_wrapper LANGUAGES C)

# Define a common library and executables for MOPAC & PARAM
add_library(mopac_math SHARED)

# MOPAC shared library ABI compatibility version
set_target_properties(mopac_math PROPERTIES SOVERSION 1)

# Install the executables and library
option(LIB_SUFFIX "Suffix for library directory (lib)" "")
if(LIB_SUFFIX)
  install(TARGETS mopac_math DESTINATION lib${LIB_SUFFIX})
else()
  install(TARGETS mopac_math DESTINATION lib)
endif()

# Check for MKLROOT environment variable
if(NOT DEFINED ENV{MKLROOT})
  message(FATAL_ERROR "Define the MKLROOT environment variable as the path to the Intel MKL library")
endif()

# Use hard-coded system-specific MKL link lines from the Link Line Advisor (& adjustments for alternate Anaconda paths)
if(APPLE)
  target_link_libraries(mopac_math PUBLIC
    $ENV{MKLROOT}/lib/libmkl_intel_lp64.a
    $ENV{MKLROOT}/lib/libmkl_intel_thread.a
    $ENV{MKLROOT}/lib/libmkl_core.a
    -liomp5 -lpthread -lm -ldl)
elseif(WIN32)
  target_link_libraries(mopac_math PUBLIC
    mkl_intel_lp64.lib
    mkl_intel_thread.lib
    mkl_core.lib
    libiomp5md.lib)
elseif(UNIX)
  find_library(MKL_LP64 libmkl_intel_lp64.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
  find_library(MKL_THREAD libmkl_intel_thread.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
  find_library(MKL_CORE libmkl_core.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
  target_link_libraries(mopac_math PUBLIC
    -Wl,--start-group
    ${MKL_LP64} ${MKL_THREAD} ${MKL_CORE}
    -Wl,--end-group
    -liomp5 -lpthread -lm -ldl)
endif()

target_sources(mopac_math PRIVATE math_wrap.c)
