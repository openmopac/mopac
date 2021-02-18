target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/analyt_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mndod_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/analyt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/bfn.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ccrep.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/coe.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/diat.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/h1elec.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/hcore.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/helect.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/jab.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/kab.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mndod.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/repp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rotate.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/set.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/solrot.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/sp_two_electron.F90
)
