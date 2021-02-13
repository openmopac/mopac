target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/ef_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/maps_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/big_swap.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/drc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/drcout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ef.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/grid.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pathk.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/paths.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/picopt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/prtdrc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/quadr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/react1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/set_up_RAPID.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/write_trajectory.F90
)
