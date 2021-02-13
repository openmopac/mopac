target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/polar_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/CPE_Energy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/bonds.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chrge.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dimens.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dipole.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mpcsyb.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mullik.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/polar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/static_polarizability.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/thermo.F90
)
