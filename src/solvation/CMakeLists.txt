target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/afmm_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/cosmo_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/cosmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/linear_cosmo.F90
)
