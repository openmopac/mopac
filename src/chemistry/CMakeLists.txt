target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/mod_atomradii.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/add_hydrogen_atoms.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/atomrs.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chkion.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/chklew.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/findn1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk_bits_1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/geochk_bits_2.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/greek.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/lewis.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/ligand.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/lyse.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/names.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/nxtmer.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/orient_water.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/reseq.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/txtype.F90
)
