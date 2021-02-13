target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/iter_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/capcor.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/cnvg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/fock1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/fock2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/iter.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/local.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pulay.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/resolv.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/swap.F90
)
