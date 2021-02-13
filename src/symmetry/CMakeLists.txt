target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/symmetry_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/bldsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charst.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/charvi.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dtran2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/dtrans.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/maksym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/molsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mult33.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/rotmol.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symh.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symoir.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symopr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/sympop.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symr.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/symtrz.F90
)
