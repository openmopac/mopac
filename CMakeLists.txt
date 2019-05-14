# Minimum CMake version requirement [to use target_sources()]
cmake_minimum_required(VERSION 3.1)

# Specify project name & programming language
project(MOPAC LANGUAGES Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(MKL REQUIRED)

# Define a library and an executable for MOPAC
add_definitions(-DBITS64) # the other valid flag is BITS32
#add_library(MOPAC_lib "")
add_executable(MOPAC_exe "")

#target_link_libraries(MOPAC_exe PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_directories(MOPAC_exe PUBLIC /opt/intel/compilers_and_libraries_2018.3.185/mac/mkl/lib)
#target_link_libraries(MOPAC_exe PUBLIC -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
target_link_libraries(MOPAC_exe PUBLIC -mkl=parallel -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.3.185/mac/compiler/lib)
# ^^^ embeds path to dynamic OpenMP library

# Add list of source files from src directory
add_subdirectory(src)

# Add tests
enable_testing()
add_subdirectory(tests)
