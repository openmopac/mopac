target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/meci_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/aababc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/aabacd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/aabbcd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/babbbc.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/babbcd.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ciosci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/diagi.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/formxy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ijkl.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/meci.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecid.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecih.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/mecip.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/partxy.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/perm.F90
)
