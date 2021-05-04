target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/to_screen_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/denrot.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/eigen.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/enpart.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geoutg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/matou1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/matout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/modchg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/modgra.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/output_rama.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pdbout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pinout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/prtpka.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/superd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/timout.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/to_screen.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/writmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/xxx.F90
)
