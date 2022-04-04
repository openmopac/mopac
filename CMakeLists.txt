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
project(openmopac VERSION 22.0.0 LANGUAGES Fortran)

# location for CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Define a common library and executables for MOPAC & PARAM
add_library(mopac_core SHARED)
set_target_properties(mopac_core PROPERTIES OUTPUT_NAME "mopac")
add_executable(mopac)
add_executable(param)
add_executable(makpol)
target_link_libraries(mopac mopac_core)
target_link_libraries(param mopac_core)

# WinMOPAC build on Windows w/ the Intel Fortran compiler
option(BUILD_WINMOPAC "Build MOPAC & BZ with QuickWin GUI (ifort required)" OFF)
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
  add_executable(winbz)
  target_compile_options(winbz PRIVATE /libs:qwin)
  target_link_options(winbz PRIVATE /subsystem:windows)
endif()

# MOPAC shared library ABI compatibility version
set_target_properties(mopac_core PROPERTIES SOVERSION 1)

# Pass OS information to MOPAC's Fortran source
if(APPLE)
  target_compile_definitions(mopac_core PRIVATE MOPAC_OS="Mac")
elseif(WIN32)
  target_compile_definitions(mopac_core PRIVATE MOPAC_OS="Windows")
elseif(UNIX)
  target_compile_definitions(mopac_core PRIVATE MOPAC_OS="Linux")
endif()

# Pass version & commit information to MOPAC's Fortran source
target_compile_definitions(mopac_core PRIVATE MOPAC_VERSION_FULL="${PROJECT_VERSION}")
execute_process(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_HASH)
if(GIT_HASH)
  string(STRIP ${GIT_HASH} GIT_HASH)
  target_compile_definitions(mopac_core PRIVATE MOPAC_GIT_HASH="${GIT_HASH}")
endif()

# Enable C language support for MKL?  (C is only needed for auto-detection of MKL by find_package)
option(ENABLE_MKL "Turn on C language support for MKL" ON)
if(ENABLE_MKL)
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
set(MOPAC_LINK "" CACHE STRING "Additional link line for MOPAC executable.")
set(MOPAC_LINK_PATH "" CACHE STRING "Path of external dependencies to bundle.")
separate_arguments(MOPAC_LINK)
separate_arguments(MOPAC_LINK_PATH)
target_link_libraries(mopac_core PUBLIC ${MOPAC_LINK})
target_link_directories(mopac_core PUBLIC ${MOPAC_LINK_PATH})

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

# CPack options
set(CPACK_PACKAGE_NAME "mopac")
set(CPACK_PACKAGE_VENDOR "openmopac GitHub organization")
set(CPACK_PACKAGE_CONTACT "openmopac@gmail.com")
file(COPY "${CMAKE_SOURCE_DIR}/README.md" DESTINATION "${CMAKE_BINARY_DIR}")
file(RENAME "${CMAKE_BINARY_DIR}/README.md" "${CMAKE_BINARY_DIR}/README.txt")
file(COPY "${CMAKE_SOURCE_DIR}/COPYING.lesser" DESTINATION "${CMAKE_BINARY_DIR}")
file(RENAME "${CMAKE_BINARY_DIR}/COPYING.lesser" "${CMAKE_BINARY_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_BINARY_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")

# Install the executables and library
set(LIB_SUFFIX "" CACHE STRING "Suffix for library install directory.")
install(TARGETS mopac RUNTIME)
install(TARGETS param RUNTIME)
install(TARGETS makpol RUNTIME)
if(WIN32)
  install(TARGETS mopac_core RUNTIME)
  install(TARGETS winmopac RUNTIME)
  install(TARGETS winbz RUNTIME)
else()
  set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR}${LIB_SUFFIX})
  install(TARGETS mopac_core LIBRARY)
endif()

# Add packaging
include(CPack)
cpack_add_component(Unspecified)
