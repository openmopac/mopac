! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module Common_arrays_C

!
!  This module contains all the arrays relating to the system being calculated.
!  Every entry is unique, and has the same meaning in every subroutine.
!
!
!  This module can also be regarded as a dictionary of the arrays relating to the system.
!  Quantities which have only a one-line descriptionis are of no interest outside MOPAC.
!  The data-type for each quantity is given at the start of the statement
!
  implicit none
  integer, dimension (:), allocatable :: &
  &  labels,     & !  Term          Atomic numbers (all atoms, real, dummy, Tv)
                   !  Units         None
                   !  Min inclusive 0
                   !  Max inclusive 107
                   !
  &  nat,        & !  Term          Atomic numbers (real atoms)
                   !  Definition    Atomic numbers of real atoms.
                   !  Min inclusive 0
                   !  Max inclusive 85
                   !
  &  na, nb, nc, & !  Term          Connectivities used in Z-matrix
                   !  Definition    Connectivity for bonds, angles, and dihedrals
  &  na_store,   & !  Store of "na" for cases where na must be temporarily zero
                   !
  &  nfirst,     & !  Term          Starting index of atomic orbitals for each atom
                   !
  &  nlast,      & !  Term          Ending index of atomic orbitals for each atom
                   !
  &  nw,         & !
  &  ifact,      & !
  &  i1fact,     & !
  &  nbonds,     & ! Number of atoms bonded to a given atom
  &  set_a,      & ! Set of atoms in first fragment, Used by PM6-DH2
  &  set_b,      & ! Set of atoms in second fragment, Used by PM6-DH2
  &  acceptor_a, & ! List of O and N atoms in first fragment. Used by PM6-DH2
  &  acceptor_b, & ! List of O and N atoms in second fragment. Used by PM6-DH2
  &  bonding_a_h,& ! List of Hydrogen atoms attached to O or N in first fragment. Used by PM6-DH2
  &  bonding_b_h,& ! List of Hydrogen atoms attached to O or N in second fragment. Used by PM6-DH2
  &  ipKa_sorted,& ! List of ionizable hydrogen atoms, sorted by pKa
  &ipKa_unsorted,& ! List of ionizable hydrogen atoms, not yet sorted by pKa
  &  dummy
  integer, dimension (:,:), allocatable :: &
  &  loc,        & !  Indices of atoms and coordinates marked for optimization
  &  ibonds,     & !  Atom numbers of atoms bonded to a given atom
  &  lopt,       & !  optimization flags
  &  lopt_store, & !  Optimization flags in GEO_REF
  &  pibonds,    & !  Pi bonds supplied by the user (requires keyword SETPI)
  &  hblist,     & !  List of hydrogen bonds X-H --- Y and types used by PM6-DH2
  &  ijpars,     & !
  &  dummy2
  double precision, dimension(3,3) :: tvec = 0.d0 !  Translation vectors (for solid-state)
  integer ::         &
   cell_ijk(3),      &   !  Dispacement translation indices (for solid-state)
   breaks(400),      &   !  Locations of breaks in PDB structure
   time_start(8),    &   !
   time_end(8)           !
  Double Precision :: Vab(3)           !  Distance vector between two atoms in a solid.(for solid-state)
  character ::chains(100)*1            !  Names of chains in proteins
  double precision :: break_coords(3,400) = 0.d0   !  Locations of breaks in a PDB structure
!
  double precision, dimension (:), allocatable :: &
  &  xparam,     & !  Values of geometric parameters marked for optimization
  &  xparef,     & !  Reference copy of xparam
  &  uspd,       & !  Initial values of one-electron diagonal (Uss, Upp, and Udd)
  &  pdiag,      & !  initial values of diagonal of density matrix
                   !
  &  h,          & !  Term          One electron Hamiltonian matrix
                   !  Definition    NDDO one-electron matrix, packed
                   !  Units         eV
                   !
  &  hb,         & !  Term          One electron beta Hamiltonian matrix
                   !  Definition    NDDO one-electron beta matrix, packed
                   !  Units         eV
                   !
  &  f,          & !  Term          Fock matrix, or alpha Fock matrix, if UHF
                   !  Definition    one-electron plus two-electron Hamiltonian matrix, packed
                   !  Units         eV
                   !
  &  fb,         & !  Term          Beta Fock matrix, if UHF
                   !  Definition    one-electron plus two-electron beta Hamiltonian matrix, packed
                   !  Units         eV
                   !
  &  w,          & !  Term          Two-electron integrals
                   !  Definition    One and two center, "J" two-electron integrals
                   !  Units         eV
                   !
  &  wk,         & !  Term          "K"-type two-electron integrals (used in solid state only)
                   !  Definition    One and two center, "J" two-electron integrals
                   !  Units         eV
                   !
  &  p,          & !  Term          Total density matrix
                   !  Definition    Coulson density matrix
                   !  Description   Calculated in SCF from M.O.s
                   !  Units         Electrons
                   !  Min inclusive -1.0 (off diagonal) 0.0 (on diagonal)
                   !  Max inclusive  1.0 (off diagonal) 2.0 (on diagonal)
                   !
  &  pa, pb,     & !  Term          Alpha and beta density matrix
                   !  Definition    Coulson alpha and beta on-electron density matrix
                   !  Description   Calculated in SCF from M.O.s
                   !  Units         Electrons
                   !  Min inclusive -1.0 (off diagonal) 0.0 (on diagonal)
                   !  Max inclusive  1.0 (off diagonal) 1.0 (on diagonal)
                   !  units:        alpha and beta electrons
                   !
  &  q,          & !  Term          Net charge on atoms
                   !  Definition    Core charge on an atom minus number of valence electrons
                   !  Units         Electrons
                   !
  &  eigs,       & !  Term          M.O. eigenvalues, or Alpha M.O. eigenvalues, if UHF
                   !  Definition    Energies of molecular orbitals or alpha M.O.s if UHF
                   !  Units         eV
                   !
  &  eigb,       & !  Term          Beta M.O. eigenvalues (if UHF)
                   !  Definition    Energies of beta molecular orbitals
                   !  Units         eV
                   !
  &  dxyz,       & !  Cartesian gradient vector in kcal/mol/Angstrom
  &  grad,       & !  Gradient vector (kcal/mol/A, or kcal/mol/radian)
  &  chrg,       & !  Mulliken charges (used by to_screen)
  &  bondab,     & !  Bond order matrix
                   !
  &  atmass,     & !  Term          Atomic masses
                   !  Definition    Masses of individual atoms
                   !  Units         Atomic mass units (Hydrogen = 1.00794 on this scale)
                   !  Default       Atomic weights given in Periodic Table
  &  ch,         & !
  &  hesinv,     & !  inverse hessian
  &  fmatrx,     & !  Hessian Force Matrix
  &  gnext1,     & !
  &  gmin1,      & !
  &  errfn,      & !
  &  profil,     & !
  &  pKa_sorted, & !  pKa of ionizable hydrogen atoms, sorted by pKa
  &  parsij,     &
  &  aicorr,     & !
  &  H_energy,   & !  Energies of hydrogen bonds found in system
  &  T_range,    & !  Term          Temperatures in thermochemistry (in thermo.F90)
                   !  Units         Kelvin
  &  HOF_tot,    & !  Term          Heats of formation at temperatures in T_range
                   !  Units         kilo-calories per mole
  &  H_tot,      & !  Term          Enthalpy at temperatures in T_range
                   !  Units         cal/mol
  &  Cp_tot,     & !  Term          Heat capacity at temperatures in T_range
                   !  Units         cal/(K.Mol)
  &  S_tot         !  Term          Entropy at temperatures in T_range
                   !  Units         cal/(K.Mol)
  double precision, dimension (:,:), allocatable :: &
  &  geo,        & !  Term          Geometry
                   !  Definition    Coordinates of all atoms in units defined by input
                   !  Units         Angstroms and radians
                   !
  &  geoa,       & !  Term          Alternative Geometry
                   !  Definition    Alternative geometry coordinates of all atoms in units defined by input
                   !  Units         Angstroms and radians
                   !
  &  coord,      & !  Term          Cartesian coordinates
                   !  Definition    Cartesian coordinates of all atoms
                   !  Units         Angstroms
                   !
  &  coorda,     & !  Term          Alternative Cartesian coordinates
                   !  Definition    Alternative Cartesian coordinates of all atoms
                   !  Units         Angstroms
                   !
  &  c,          & !  Term          Eigenvectors, or alpha eigenvectors, if UHF
                   !  Definition    Molecular orbital coefficients of alpha M.O.s.  c**2 = 1.0
                   !  Units         dimensionless
                   !  Min inclusive -1.0
                   !  Max inclusive  1.0
                   !
  &  cb,         & !  Term          Eigenvectors, or beta eigenvectors, if UHF
                   !  Definition    Molecular orbital coefficients of beta M.O.s.  c**2 = 1.0
                   !  Units         dimensionless
                   !  Min inclusive -1.0
                   !  Max inclusive  1.0
  &  dxyz2,      & !  dxyz, but 2-D
  &  ptot2,      & !
     xmat,       &
  & rdummy2

  character, dimension (:), allocatable :: &
  &  simbol*10,      & !  Term          Symbol for geometric variable
                       !  Definition    Gaussian format symbolic
                       !  Unit          10 characters
  & all_comments*81, & !  Term          Comments
                       !  Definition    Comments in data set describing the data set. Not used in program working.
                       !                (Useful for holding PDB information)
  &  pibonds_txt*50, & !  pi bonds defined using SETPI, as text
  &  H_txt*120,      & !  Term          Description of hydrogen bonds
                       !  Definition    Atoms involved in H - bond, hydrogen bond-length and energy.
  &  txtatm*27,      & !  Term          Atom discriptor
                       !  Definition    A textual description of each atom.  Not used in program working.
  &  txtatm1*27        !  Term          Original atom discriptor

  logical, dimension (:), allocatable :: &
  &  l_atom            !  If true, then use atom in output, otherwise don't print

end module Common_arrays_C
