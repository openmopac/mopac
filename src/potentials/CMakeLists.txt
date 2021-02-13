target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/esp_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/esp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/genun.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/new_esp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/pmep.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rotate_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/setupg.F90
)
