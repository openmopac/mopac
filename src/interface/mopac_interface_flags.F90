module mopac_interface_flags
!dec$ attributes dllexport :: run_mopac
  implicit none

  ! flag that denotes that the API is used
  logical :: use_api = .false.

  ! flag that keeps track of if mopac
  logical :: reset_mopac_initialize=.false.

  ! flag that keeps track of if an SCF is available
  logical :: SCF_avail = .false.

  ! flags to reset all needed subroutines/functions
  logical :: reset_iter_L = .false.
  logical :: reset_fock2_L = .false.
  logical :: reset_fock2z_L = .false.
  logical :: reset_fock2zn_L = .false.
  logical :: reset_diagg2_L = .false.
  logical :: reset_diagg1_L = .false.
  logical :: reset_dipole_for_MOZYME_L = .false.
  logical :: reset_iter_for_MOZYME_L = .false.
  logical :: reset_outer2_L = .false.
  logical :: reset_tidy_L = .false.
  logical :: reset_scfcri_L = .false.
  logical :: reset_pulay_L = .false.
  logical :: reset_cnvg_L = .false.
  logical :: reset_geochk_L = .false.
  logical :: reset_check_h_L
  logical :: reset_check_CVS_L = .false.
  logical :: reset_compfg_L = .false.
  logical :: reset_Hbond_corr_PM6_DH_L = .false.
  logical :: reset_dftd3_L = .false.
  logical :: reset_prt_hbonds_L = .false.
  logical :: reset_Hydrogen_bond_corrections_L = .false.
  logical :: reset_H_bonds4_L = .false.
  logical :: reset_PM6_DH_Disp_L = .false.
  logical :: reset_deritr_L = .false.
  logical :: reset_dfock2_L = .false.
  logical :: reset_dcart_L = .false.
  logical :: reset_derp_L = .false.
  logical :: reset_deri1_L = .false.
  logical :: reset_dhc_L = .false.
  logical :: reset_deri2_L = .false.
  logical :: reset_deriv_L = .false.
  logical :: reset_dernvo_L = .false.
  logical :: reset_gmetry_L = .false.
  logical :: reset_axis_L = .false.
  logical :: reset_haddon_L = .false.
  logical :: reset_hcore_L = .false.
  logical :: reset_solrot_L = .false.
  logical :: reset_trunk_L = .false.
  logical :: reset_diat_L = .false.
  logical :: reset_rotatd_L = .false.
  logical :: reset_interp_L = .false.
  logical :: reset_meci_L = .false.
  logical :: reset_powsq_L = .false.
  logical :: reset_search_L = .false.
  logical :: reset_flepo_L = .false.
  logical :: reset_dfpsav_L = .false.
  logical :: reset_linmin_L = .false.
  logical :: reset_locmin_L = .false.
  logical :: reset_dipole_L = .false.
  logical :: reset_mullik_L = .false.
  logical :: reset_ef_L = .false.
  logical :: reset_overlp_L = .false.
  logical :: reset_formd_L = .false.
  logical :: reset_updhes_L = .false.
  logical :: reset_drcout_L = .false.
  logical :: reset_prtdrc_L = .false.
  logical :: reset_molsym_L = .false.

contains

  subroutine mopac_interface_flags_set_true()
    implicit none
    use_api = .true.
  
    reset_iter_L = .true.
    reset_fock2_L = .true.
    reset_fock2z_L = .true.
    reset_fock2zn_L = .true.
    reset_diagg2_L = .true.
    reset_diagg1_L = .true.
    reset_dipole_for_MOZYME_L = .true.
    reset_iter_for_MOZYME_L = .true.
    reset_outer2_L = .true.
    reset_tidy_L = .true.
    reset_scfcri_L = .true.
    reset_pulay_L = .true.
    reset_cnvg_L = .true.
    reset_geochk_L = .true.
    reset_check_h_L = .true.
    reset_check_CVS_L = .true.
    reset_compfg_L = .true.
    reset_Hbond_corr_PM6_DH_L = .true.
    reset_dftd3_L = .true.
    reset_prt_hbonds_L = .true.
    reset_Hydrogen_bond_corrections_L = .true.
    reset_H_bonds4_L = .true.
    reset_PM6_DH_Disp_L = .true.
    reset_deritr_L = .true.
    reset_dfock2_L = .true.
    reset_dcart_L = .true.
    reset_derp_L = .true.
    reset_deri1_L = .true.
    reset_dhc_L = .true.
    reset_deri2_L = .true.
    reset_deriv_L = .true.
    reset_dernvo_L = .true.
    reset_gmetry_L = .true.
    reset_axis_L = .true.
    reset_haddon_L = .true.
    reset_hcore_L = .true.
    reset_solrot_L = .true.
    reset_trunk_L = .true.
    reset_diat_L = .true.
    reset_rotatd_L = .true.
    reset_interp_L = .true.
    reset_meci_L = .true.
    reset_powsq_L = .true.
    reset_search_L = .true.
    reset_flepo_L = .true.
    reset_dfpsav_L = .true.
    reset_linmin_L = .true.
    reset_locmin_L = .true.
    reset_dipole_L = .true.
    reset_mullik_L = .true.
    reset_ef_L = .true.
    reset_overlp_L = .true.
    reset_formd_L = .true.
    reset_updhes_L = .true.
    reset_drcout_L = .true.
    reset_prtdrc_L = .true.
    reset_molsym_L = .true.
    return
  end subroutine mopac_interface_flags_set_true

  subroutine  mopac_interface_flags_set_false()
    implicit none
    use_api = .false.
  
    SCF_avail = .false.
    reset_iter_L = .false.
    reset_fock2_L = .false.
    reset_fock2z_L = .false.
    reset_fock2zn_L = .false.
    reset_diagg2_L = .false.
    reset_diagg1_L = .false.
    reset_dipole_for_MOZYME_L = .false.
    reset_iter_for_MOZYME_L = .false.
    reset_outer2_L = .false.
    reset_tidy_L = .false.
    reset_scfcri_L = .false.
    reset_pulay_L = .false.
    reset_cnvg_L = .false.
    reset_geochk_L = .false.
    reset_check_h_L = .false.
    reset_check_CVS_L = .false.
    reset_compfg_L = .false.
    reset_Hbond_corr_PM6_DH_L = .false.
    reset_dftd3_L = .false.
    reset_prt_hbonds_L = .false.
    reset_Hydrogen_bond_corrections_L = .false.
    reset_H_bonds4_L = .false.
    reset_PM6_DH_Disp_L = .false.
    reset_deritr_L = .false.
    reset_dfock2_L = .false.
    reset_dcart_L = .false.
    reset_derp_L = .false.
    reset_deri1_L = .false.
    reset_dhc_L = .false.
    reset_deri2_L = .false.
    reset_deriv_L = .false.
    reset_dernvo_L = .false.
    reset_gmetry_L = .false.
    reset_axis_L = .false.
    reset_haddon_L = .false.
    reset_hcore_L = .false.
    reset_solrot_L = .false.
    reset_trunk_L = .false.
    reset_diat_L = .false.
    reset_rotatd_L = .false.
    reset_interp_L = .false.
    reset_meci_L = .false.
    reset_powsq_L = .false.
    reset_search_L = .false.
    reset_flepo_L = .false.
    reset_dfpsav_L = .false.
    reset_linmin_L = .false.
    reset_locmin_L = .false.
    reset_dipole_L = .false.
    reset_mullik_L = .false.
    reset_ef_L = .false.
    reset_overlp_L = .false.
    reset_formd_L = .false.
    reset_updhes_L = .false.
    reset_drcout_L = .false.
    reset_prtdrc_L = .false.
    reset_molsym_L = .false.
    return
  end subroutine mopac_interface_flags_set_false


end module mopac_interface_flags
