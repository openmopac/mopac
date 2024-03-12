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

# Using Cmake version 3.14 for NumPy support
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# Ensure we are building out-of-source so that the tests work (issue 120)
get_filename_component(srcdir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(bindir "${CMAKE_BINARY_DIR}" REALPATH)
if("${srcdir}" STREQUAL "${bindir}")
  message(FATAL_ERROR "MOPAC should not be configured & built in the source directory
Please create a separate build directory and rerun cmake from the build directory
For example: mkdir build; cd build; cmake ..")
endif()

# Read version number
file(READ "CITATION.cff" CITATION)
string(REGEX MATCH "[^-]version: ([0-9\.]*)" _ ${CITATION})

# Specify project name & programming languages
project(openmopac VERSION ${CMAKE_MATCH_1} LANGUAGES Fortran)

# Enable C language support for MKL?  (C is only needed for auto-detection of MKL by find_package)
option(ENABLE_MKL "Turn on C language support for MKL" ON)
if(ENABLE_MKL)
  enable_language(C)
endif()

# Follow GNU conventions for installing directories (languages need to be enabled at this stage)
include(GNUInstallDirs)

# Failsafe definition of install directories (not guaranteed behavior in some older CMake versions)
if(NOT DEFINED CMAKE_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "bin")
endif()
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR "lib")
endif()

# Static or dynamic build
option(STATIC_BUILD "Build a static executable" OFF)
if (STATIC_BUILD)
  message(STATUS "Static executables")
  message(WARNING "CMake cannot verify the ability to link static executables in the pre-build phase. "
    "Please verify that you are able to link static executables with a \"-static\" linker flag. "
    "OpenMP does not support static builds by default, so try THREADS_KEYWORD=OFF if the static build fails to turn off OpenMP.")
else()
  message(STATUS "Shared library and dynamic executables")
endif()

# Code coverage option
option(ENABLE_COVERAGE "Code coverage flag" OFF)

# Set a safe default build type if there is no user input
if(NOT CMAKE_BUILD_TYPE)
  if (ENABLE_COVERAGE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  else()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  endif()
else()
# otherwise, check conflict with ENABLE_COVERAGE
  if(ENABLE_COVERAGE AND NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    message(FATAL_ERROR "Option ENABLE_COVERAGE requires CMAKE_BUILD_TYPE=Debug")
  endif()
endif()

message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

# location for CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set relative RPATH for installation
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Install-time options
if(APPLE)
  set(CMAKE_INSTALL_RPATH "\@loader_path/../lib")
elseif(UNIX)
  set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib")
endif()

# Define a common library and executables for MOPAC & PARAM
if (STATIC_BUILD)
  # Static build
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  add_library(mopac-core STATIC)
  target_link_options(mopac-core PUBLIC "-static")
  set(BLA_STATIC ON)
else()
  add_library(mopac-core SHARED)
endif()
set_target_properties(mopac-core PROPERTIES OUTPUT_NAME "mopac")
add_executable(mopac)
add_executable(mopac-param)
add_executable(mopac-makpol)
target_link_libraries(mopac mopac-core)
target_link_libraries(mopac-param mopac-core)
if (STATIC_BUILD)
  target_link_options(mopac PUBLIC "-static")
  target_link_options(mopac-param PUBLIC "-static")
  target_link_options(mopac-makpol PUBLIC "-static")
endif()

# Add CITATION.cff to source files so that version updates are not ignored
target_sources(mopac-core PRIVATE CITATION.cff)

# BZ build on Windows w/ the Intel Fortran compiler
option(BUILD_BZ "Build BZ with QuickWin GUI (ifort required)" OFF)
if(BUILD_BZ)
  if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"))
    message(FATAL_ERROR "BZ compilation requires the Intel Fortran compiler (ifort)")
  endif()
  if(NOT WIN32)
    message(FATAL_ERROR "BZ is only available on Windows")
  endif()
  add_executable(mopac-bz)
  target_compile_options(mopac-bz PRIVATE /libs:qwin)
  target_link_options(mopac-bz PRIVATE /subsystem:windows)
  if (STATIC_BUILD)
    target_link_options(mopac-bz PUBLIC "-static")
  endif()
endif()

# MOPAC shared library ABI compatibility version
set_target_properties(mopac-core PROPERTIES SOVERSION 1)

# Adjust legacy PGI compiler behavior: disable C-style backslash characters & warn about CRLF handling
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|Flang")
  target_compile_options(mopac-core PUBLIC "-fno-backslash")
  target_compile_options(mopac-makpol PUBLIC "-fno-backslash")
  if(UNIX)
    message(NOTICE "Legacy PGI compiler detected. Set environment variable FORTRANOPT=crlf at compile time to parse Windows-style CRLF endline sequences in MOPAC input files.")
  endif()
endif()

# Pass OS information to MOPAC's Fortran source
if(APPLE)
  target_compile_definitions(mopac-core PRIVATE MOPAC_OS="MacOS")
elseif(WIN32)
  target_compile_definitions(mopac-core PRIVATE MOPAC_OS="Windows")
elseif(UNIX)
  target_compile_definitions(mopac-core PRIVATE MOPAC_OS="Linux")
endif()

# Pass compiler information to MOPAC's Fortran source
if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  option(F2003_INTRINSICS "Replace non-standard intrinsics w/ Fortran 2003 standards" ON)
else()
  option(F2003_INTRINSICS "Replace non-standard intrinsics w/ Fortran 2003 standards" OFF)
endif()
if(F2003_INTRINSICS)
  target_compile_definitions(mopac-core PRIVATE MOPAC_F2003)
  target_compile_definitions(mopac-makpol PRIVATE MOPAC_F2003)
  if(BUILD_BZ)
    target_compile_definitions(mopac-bz PRIVATE MOPAC_F2003)
  endif()
endif()

# Pass version & commit information to MOPAC's Fortran source
target_compile_definitions(mopac-core PRIVATE MOPAC_VERSION_FULL="${PROJECT_VERSION}")
if(NOT GIT_HASH)
  execute_process(COMMAND git rev-parse HEAD OUTPUT_VARIABLE GIT_HASH)
endif()
if(GIT_HASH)
  string(STRIP ${GIT_HASH} GIT_HASH)
  target_compile_definitions(mopac-core PRIVATE MOPAC_GIT_HASH="${GIT_HASH}")
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
  target_link_libraries(mopac-core PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  target_link_options(mopac PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
  target_link_options(mopac-param PUBLIC ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
else()
  message(STATUS "AUTO_BLAS=OFF, the user is responsible for linking BLAS & LAPACK libraries (e.g. define compiler link line with MOPAC_LINK)")
endif()

# Use OpenMP API in MOPAC for "THREADS" keyword
option(THREADS_KEYWORD "Enable THREADS keyword using OpenMP API" ON)
if(THREADS_KEYWORD)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    target_compile_options(mopac-core PUBLIC ${OpenMP_Fortran_FLAGS})
    target_link_options(mopac-core PUBLIC ${OpenMP_Fortran_FLAGS})
  endif()
endif()

# Additional link line
set(MOPAC_LINK "" CACHE STRING "Additional link line for MOPAC executable.")
set(MOPAC_LINK_PATH "" CACHE STRING "Path of external dependencies to bundle.")
separate_arguments(MOPAC_LINK)
separate_arguments(MOPAC_LINK_PATH)
target_link_libraries(mopac-core PUBLIC ${MOPAC_LINK})
target_link_directories(mopac-core PUBLIC ${MOPAC_LINK_PATH})

# Add list of source files from src directory
add_subdirectory(src)

# Set up code coverage
set(GCOV_PROGRAM_NAME gcov CACHE STRING "Name of GNU code coverage analysis tool.")
set(LCOV_PROGRAM_NAME lcov CACHE STRING "Name of GNU code coverage visualization tool.")
if(ENABLE_COVERAGE)
  if(NOT ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU"))
    message(FATAL_ERROR "Code coverage requires the GNU Fortran compiler! Aborting...")
  endif()
  target_compile_options(mopac-core PUBLIC -fprofile-arcs -ftest-coverage -O0)
  target_link_options(mopac-core PUBLIC -fprofile-arcs -ftest-coverage)
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
option(TESTS "Enable tests (needs Python3 and Numpy)" ON)
if(TESTS)
  find_package (Python3 COMPONENTS Interpreter NumPy)
  if(NOT Python3_FOUND)
    if(Python3_Interpreter_FOUND AND (NOT Python3_NumPy_FOUND))
      message(NOTICE "Python3 found without NumPy at path: ${Python3_EXECUTABLE}\n"
              "(HINT: install NumPy for this Python3 instance with pip or give higher priority to an instance with NumPy installed)")
    endif()
    message(FATAL_ERROR "Python3 and Numpy are required for MOPAC testing (testing can be disabled with -DTESTS=OFF)")
  endif()
  enable_testing()
  add_subdirectory(tests)
endif()

# CPack options
set(CPACK_PACKAGE_NAME "mopac")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Molecular Orbital PACkage (MOPAC)")
set(CPACK_PACKAGE_VENDOR "openmopac GitHub organization")
set(CPACK_PACKAGE_CONTACT "openmopac@gmail.com")
if(APPLE)
  set(CPACK_SYSTEM_NAME "mac")
elseif(WIN32)
  set(CPACK_SYSTEM_NAME "win")
elseif(UNIX)
  set(CPACK_SYSTEM_NAME "linux")
endif()
set(CPACK_PROJECT_CONFIG_FILE "${CMAKE_SOURCE_DIR}/cmake/CPackOptions.cmake")

# Package the Intel OpenMP library for the QT installer (hard-coded for each OS)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if(WIN32)
    find_file(INTEL_OMP_DLL libiomp5md.dll HINT CMAKE_Fortran_COMPILER PATH_SUFFIXES ../../redist/intel64_win/compiler)
    if(INTEL_OMP_DLL STREQUAL "INTEL_OMP_DLL-NOTFOUND")
      message("WARNING: Cannot locate Intel's OpenMP dll for packaging")
    endif()
    install(PROGRAMS ${INTEL_OMP_DLL} DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT qtifw EXCLUDE_FROM_ALL)
  else()
    find_library(INTEL_OMP_LIB iomp5 HINT CMAKE_Fortran_COMPILER PATH_SUFFIXES ../../compiler/lib ../../compiler/lib/intel64_lin)
    if(INTEL_OMP_LIB STREQUAL "INTEL_OMP_LIB-NOTFOUND")
      message("WARNING: Cannot locate Intel's OpenMP shared library for packaging")
    endif()
    # resolve symbolic links
    get_filename_component(INTEL_OMP_LIB_PATH ${INTEL_OMP_LIB} REALPATH)
    message(STATUS "OMP packaging: ${INTEL_OMP_LIB} -> ${INTEL_OMP_LIB_PATH}")
    install(FILES ${INTEL_OMP_LIB_PATH} DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT redist EXCLUDE_FROM_ALL)
  endif()
endif()

# IFW-specific CPack options
set(CPACK_IFW_PACKAGE_NAME "MOPAC")
set(CPACK_IFW_PRODUCT_URL "http://openmopac.net")
set(CPACK_IFW_PACKAGE_CONTROL_SCRIPT "${CMAKE_SOURCE_DIR}/.github/qtifw_controller.qs")
set(CPACK_IFW_PACKAGE_MAINTENANCE_TOOL_NAME "uninstall-mopac")
set(CPACK_IFW_VERBOSE ON)
if(WIN32)
  set(CPACK_IFW_TARGET_DIRECTORY "C:/Program Files/MOPAC")
else()
  set(CPACK_IFW_TARGET_DIRECTORY "/opt/mopac")
endif()
install(FILES "${CMAKE_SOURCE_DIR}/.github/mopac.ico" DESTINATION "." COMPONENT qtifw EXCLUDE_FROM_ALL)
install(FILES "${CMAKE_SOURCE_DIR}/CITATION.cff"
              "${CMAKE_SOURCE_DIR}/COPYING"
              "${CMAKE_SOURCE_DIR}/COPYING.lesser" DESTINATION "." COMPONENT redist EXCLUDE_FROM_ALL)

# Install the executables and library
install(TARGETS mopac RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT main)
install(TARGETS mopac-param RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT main)
install(TARGETS mopac-makpol RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT main)
if(WIN32)
  install(TARGETS mopac-core RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT main)
  if(BUILD_BZ)
    install(TARGETS mopac-bz RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT main)
  endif()
else()
  install(TARGETS mopac-core LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT main)
endif()

# Add IFW packaging
include(CPack)
include(CPackIFW)
cpack_ifw_configure_component(main DISPLAY_NAME "main program" SORTING_PRIORITY 2 FORCED_INSTALLATION
  LICENSES "primary license (LGPL3)" "${CMAKE_SOURCE_DIR}/COPYING.lesser"
  "secondary license terms (GPL3)" "${CMAKE_SOURCE_DIR}/COPYING"
  SCRIPT "${CMAKE_SOURCE_DIR}/.github/qtifw_component.qs"
  USER_INTERFACES "${CMAKE_SOURCE_DIR}/.github/pathcheckboxform.ui" "${CMAKE_SOURCE_DIR}/.github/filecheckboxform.ui" "${CMAKE_SOURCE_DIR}/.github/iconcheckboxform.ui")
cpack_ifw_configure_component(qtifw VIRTUAL FORCED_INSTALLATION)
cpack_ifw_configure_component(redist VIRTUAL FORCED_INSTALLATION)
