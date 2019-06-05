# Minimum CMake version requirement [to use target_sources()]
cmake_minimum_required(VERSION 3.1)

# Specify project name & programming language
project(MOPAC LANGUAGES Fortran)

# Define a library and an executable for MOPAC
add_definitions(-DBITS64) # the other valid flag is BITS32
#add_library(MOPAC_lib "") encapsulate a MOPAC library at some point ...
add_executable(mopac.exe "")

message("compiler = " ${CMAKE_Fortran_COMPILER})

# Link w/ MKL if it is available
if(DEFINED ENV{MKLROOT})
  message("environment variable MKLROOT is defined, assuming that Intel MKL is available")
  add_definitions(-DMKL)
  if(DEFINED ENV{OMPROOT})
    message("environment variable OMPROOT is defined, building threaded MKL w/ OMPROOT as the path to lib/libiomp5")
    target_link_libraries(mopac.exe PUBLIC $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a 
                          $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{OMPROOT}/lib/ -liomp5 -lpthread -lm -ldl)
  else()
    message("environment variable OMPROOT is not defined, building serial MKL")
    target_link_libraries(mopac.exe PUBLIC $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_sequential.a
                          $ENV{MKLROOT}/lib/libmkl_core.a -lpthread -lm -ldl)
  endif()
# Otherwise try generic BLAS/LAPACK libraries
else()
  message("environment variable MKLROOT is not defined, relying on CMake to find working BLAS & LAPACK libraries")
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  target_link_libraries(mopac.exe PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif()

#target_link_libraries(mopac.exe PUBLIC -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)

#target_link_libraries(mopac.exe PUBLIC $ENV{MKLROOT}/lib/libmkl_intel_ilp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a 
#                      $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
# ^^^ with -liomp5, which is causing linking problems
# THIS IS TOO IFORT SPECIFIC: target_link_libraries(mopac.exe PUBLIC -mkl=parallel -Wl,-rpath,/opt/intel/compilers_and_libraries_2018.3.185/mac/compiler/lib)
# ^^^ embeds path to dynamic OpenMP library

# Add list of source files from src directory
add_subdirectory(src)

# Add tests
enable_testing()
add_subdirectory(tests)
