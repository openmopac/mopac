# Minimum CMake version requirement [to use target_sources()]
cmake_minimum_required(VERSION 3.1)

# Specify project name & programming language
project(MOPAC LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MKL REQUIRED)

# Define a library and an executable for MOPAC
add_definitions(-DBITS64) # the other valid flag is BITS32
#add_library(MOPAC_lib "")
add_executable(MOPAC_exe "")

# Add list of source files from src directory
add_subdirectory(src)

# Build commands ...
