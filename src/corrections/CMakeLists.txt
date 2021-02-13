target_sources(mopac-all PRIVATE
# modules
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_EC_plus_ER.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_EH_plus.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_PM6_DH_Dispersion.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_PM6_DH_type.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bond_correction_bits.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/H_bonds4.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/Hydrogen_bond_corrections.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/disp_DnX.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/post_scf_corrections.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/set_up_dentate.F90
)
