# download/ build MDI library
# always build static library with -fpic
# support cross-compilation and ninja-build
include(ExternalProject)
ExternalProject_Add(mdi_build
  URL     "https://github.com/MolSSI-MDI/MDI_Library/archive/v1.4.26.tar.gz"
  URL_MD5 "3124bb85259471e2a53a891f04bf697a"
  CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
  -Dlanguage=Fortran
  -Dlibtype=SHARED
  -Dmpi=OFF
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_BYPRODUCTS "<BINARY_DIR>/MDI_Library/libmdi${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )

# where is the compiled library?
ExternalProject_get_property(mdi_build BINARY_DIR)
set(MDI_BINARY_DIR "${BINARY_DIR}/MDI_Library")
# workaround for older CMake versions
file(MAKE_DIRECTORY ${MDI_BINARY_DIR})

# create imported target for the MDI library
add_library(mdi-lib UNKNOWN IMPORTED)
add_dependencies(mdi-lib mdi_build)
set_target_properties(mdi-lib PROPERTIES
  IMPORTED_LOCATION "${MDI_BINARY_DIR}/libmdi${CMAKE_SHARED_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES ${MDI_BINARY_DIR}
  )

link_directories(${MDI_BINARY_DIR})
include_directories(${MDI_BINARY_DIR})

add_dependencies(mopac-core mdi_build)
target_link_libraries(mopac-core PUBLIC mdi-lib)
