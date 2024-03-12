! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.


subroutine post_scf_corrections(correction, l_grad)
!
!    Add dispersion, hydrogen bonding, extra H-H repulsion, and other energies to
!    improve intermolecular interaction geometries and energies.
!
  use molkst_C, only : keywrd, E_disp, E_hb, E_hh, method_pm7, P_Hbonds, &
    method_pm6_dh_plus, method_pm6_dh2, method_pm6_d3h4, method_pm6_dh2x, method_pm6_d3h4x, &
    method_pm6_d3, method_pm6_d3_not_h4, method_pm7_hh, method_pm7_minus, method_pm6_org, method_PM8
  use common_arrays_C, only: dxyz
  implicit none
  double precision, intent(out) ::  correction
  logical, intent (in) :: l_grad
!
! Local variables
!
  logical ::  prt
  double precision, external :: & !                  Original references
                          !
  energy_corr_hh_rep,   & ! Rezac J., Hobza P., "Advanced Corrections of Hydrogen Bonding and
                          ! Dispersion for Semiempirical Quantum Mechanical Methods", J. Chem.
                          ! Theory and Comp 8, 141-151 (2012).
                          !
                          !
  PM6_DH_Dispersion,    & ! S. Grimme, "Accurate description of van der Waals complexes by
                          ! density functional theory including empirical corrections."
                          ! J Comput Chem. 2004 Sep;25(12):1463-73.
                          !
  dftd3,                & ! Grimme S., Antony J., Ehrlich S., Krieg H., "A consistent and accurate
                          ! ab initio parametrization of density functional dispersion correction
                          ! (DFT-D) for the 94 elements H-Pu", J. Chem. Phys. 132, 154104:154104
                          ! (2010).
                          !
  disp_DnX,             & ! (DH2X) Rezac and Hobza's correction: "A halogen-bonding correction
                          ! for the semiempirical PM6 method" Chem. Phys. Lett. 506 286-289 (2011)
                          !
                          ! (D3H4X) Rezac J., Hobza P., "Advanced Corrections of Hydrogen Bonding
                          ! and Dispersion for Semiempirical Quantum Mechanical Methods", J. Chem.
                          ! Theory and Comp 8, 141-151 (2012)
                          !
  Hydrogen_bond_corrections !(DH2) Korth M., Pitonak M., Rezac J., Hobza P., "A Transferable
                          ! H-bonding Correction for Semiempirical Quantum-Chemical Methods",
                          ! J. Chem. Theory and Computation 6, 344-352 (2010).
                          !
                          ! (DH+) Korth M., "Third-Generation Hydrogen-Bonding Corrections for
                          ! Semiempirical QM Methods and Force Fields", J. Chem. Theory Comput.
                          ! 6, 3808-3816 (2010).
!
! PM6-D3H4X: Brahmkshatriya, P. S.; Dobes, P.; Fanfrlik, J.; Rezac, J.; Paruch, K.; Bronowska,
! A.; Lepsik, M.; Hobza, P. "Quantum Mechanical Scoring: Structural and Energetic Insights into
! Cyclin-Dependent Kinase 2 Inhibition by Pyrazolo[1,5-a]pyrimidines" Curr. Comput.-Aid. Drug.
! 2013 , 9 (1), 118�129.
!
  prt = (index(keywrd," 0SCF ") + index(keywrd," PRT ") /= 0 .and. index(keywrd," DISP") /= 0)
  correction = 0.d0
  E_hb       = 0.d0
  E_hh       = 0.d0
  E_disp     = 0.d0
  P_Hbonds   = 0
  if (.not. allocated(dxyz)) allocate (dxyz(1))
!
! All the hydrogen-bond corrections are in Hydrogen_bond_corrections
!
  if (method_pm6_d3h4x) then
    correction = correction + dftd3(l_grad, dxyz)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
    correction = correction + disp_DnX(l_grad)
  else if (method_pm6_d3h4) then
    correction = correction + dftd3(l_grad, dxyz)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
  else if (method_pm6_d3_not_h4) then
    correction = correction + dftd3(l_grad, dxyz)
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
  else if (method_pm6_d3) then
    correction = correction + dftd3(l_grad, dxyz)
  else if (method_pm6_dh_plus) then
    correction = correction + PM6_DH_Dispersion(l_grad)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
  else if (method_pm6_dh2) then
!
! Hydrogen_bond_corrections uses partial atomic charges if method_pm6_dh2 .or. method_pm6_dh2x
!
    correction = correction + PM6_DH_Dispersion(l_grad)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
  else if (method_pm6_dh2x) then
    correction = correction + PM6_DH_Dispersion(l_grad)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
    correction = correction + disp_DnX(l_grad)
  else if (method_pm7_hh) then
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
    correction = correction + PM6_DH_Dispersion(l_grad)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
  else if (method_pm7_minus) then
    return
  else if (method_pm6_org) then
    correction = correction + dftd3(l_grad, dxyz)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
  else if (method_pm8) then
    correction = correction + dftd3(l_grad, dxyz)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
    correction = correction + energy_corr_hh_rep(l_grad, dxyz)
  else if (method_pm7) then
    correction = correction + PM6_DH_Dispersion(l_grad)
    correction = correction + Hydrogen_bond_corrections(l_grad, prt)
  end if
  if (index(keywrd, " SILENT") == 0) then
    if (prt .and. P_Hbonds > 0) call print_post_scf_corrections
  end if
  return
end subroutine post_scf_corrections
