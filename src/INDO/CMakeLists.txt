target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/reimers_C.F90
# interfaces
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/ci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/cosci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/excite.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/input.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/integral.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ovlp.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/scf.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/subs.F90
)
