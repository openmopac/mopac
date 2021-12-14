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

# Look for static libraries of Intel MKL
if(UNIX)
  find_library(MKL_LP64 libmkl_intel_lp64.a)
  find_library(MKL_THREAD libmkl_intel_thread.a)
  find_library(MKL_CORE libmkl_core.a)
  find_library(INTEL_OMP iomp5)
elseif(WIN32)
  find_library(MKL_LP64 mkl_intel_lp64.lib)
  find_library(MKL_THREAD mkl_intel_thread.lib)
  find_library(MKL_CORE mkl_core.lib)
  find_file(INTEL_OMP libomp-fallback-cstring.obj)
else()
  message(FATAL_ERROR "Unknown OS, failed to find static libraries for Intel MKL")
endif()
if(APPLE)
  target_link_libraries(mopac_math PUBLIC
    ${MKL_LP64} ${MKL_THREAD} ${MKL_CORE}
    ${INTEL_OMP}
    -lpthread -lm -ldl)
elseif(WIN32)
  target_link_libraries(mopac_math PUBLIC ${MKL_LP64} ${MKL_THREAD} ${MKL_CORE})
  add_library(openmp_obj OBJECT IMPORTED GLOBAL)
  set_property(TARGET openmp_obj PROPERTY IMPORTED_OBJECTS ${INTEL_OMP})
  target_sources(mopac_math PRIVATE $<TARGET_OBJECTS:openmp_obj>)
elseif(UNIX)
  target_link_libraries(mopac_math PUBLIC
    -Wl,--start-group
    ${MKL_LP64} ${MKL_THREAD} ${MKL_CORE}
    -Wl,--end-group
    ${INTEL_OMP}
    -lpthread -lm -ldl)
else()
  message(FATAL_ERROR "Unknown OS, failed to link static libraries for Intel MKL")
endif()

target_sources(mopac_math PRIVATE math_wrap.c)
