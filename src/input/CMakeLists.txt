target_sources(mopac-all PRIVATE
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/datin.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/digit.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geo_ref.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getdat.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getgeg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getgeo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getpdb.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getsym.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/gettxt.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/getval.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/init_filenames.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/myword.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/newflg.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/nuchar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/readmo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/wrt_diffs.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/wrtkey.F90
)
