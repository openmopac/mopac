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

module MOZYME_C
!
!  This module holds all the data that is specific to the MOZYME Localized
!  Molecular Orbital (LMO) procedure.
!
    implicit none
    integer :: i
!
!                               Data for the system
!
  integer :: &
    noccupied,   & ! Number of occupied molecular orbitals
    nvirtual,    & ! Number of virtual molecular orbitals
    nres,        & ! Number of amino-acid residues in a polypeptide or protein
    morb,        & ! Maximum number of orbitals on any atom, usually 4 or 9
    ipad2,       & ! Estimate of the average number of atoms in a LMO
    ipad4,       & ! Estimate of the average number of atomic orbital coefficients in a LMO
    uni_res,     & ! Number of unique residues in a protein
    Lewis_tot,   & ! Total number of Lewis elements (occupied plus virtual)
    Lewis_max      ! Maximum number of Lewis elements (size of Lewis_elem buffer)
!
!                              Data on atoms
!
  integer, allocatable, dimension (:) :: &
    iorbs,       &
    jopt
!
!                              Data on diatomic interactions
!
  integer ::     &
    ij_dim         ! Number of interactions treated by NDDO approximations
                   ! = size of arrays iijj and ijall
  logical :: &
    lijbo          ! TRUE if array ijbo is to be used,
                   ! FALSE if the ijbo function is to be used
  integer, dimension (:), allocatable :: &
    iijj,        & !
    iii,         &
    iij,         &
    numij,       &
    ijall          !
  integer, allocatable, dimension (:,:) :: &
    nijbo          ! If it exists, nijbo holds the statring address of array elements
                   ! for atoms i and j in otomic orbital arrays such as P, H, F, etc.
!
!  Protein-specific data
!
  integer, parameter :: &
    maxres = 9999    !  Maximum number of residues allowed in a protein

  integer, dimension(:), allocatable :: &
                  !
                  ! The unique residue number of atom i is at_res(i)
    at_res,     & ! This number is independent of the residue number printed
                  ! and is only used inside MOPAC - it is never printed out
    iz,         & ! Number of valence electrons available for Lewis structure, on each atom
    ib,         & ! Number of valence orbitals available for Lewis structure, on each atom
    ions          ! Charge on ionized atoms, cations = 1 neutral = 0

                  !
  integer, dimension(:,:), allocatable :: &
    Lewis_elem       ! Atoms and atom pairs involved in making the Lewis structure
  integer ::     &
    iatom,       & ! Used in working out residue names
    jatom,       & ! Used in working out residue names
    mxeno,       & ! Used in working out residue names
    k,           &
    loop,        &
    icharges       ! Number of charged sites
   integer :: &
    nxeno(4,11), & ! Up to 10 xeno groups allowed. nxeno(1,*) = number of C,
                   ! 2 = No. on N, 3 = #O, 4 = #S
    nbackb(4),   & ! for each peptide line, nbackb(1) = atom number of C of CHR,
                   ! 2 = # of C of C=O, 3 = # of O of C=O, 4 = # of N
                        !
                        ! For each entry in in_res, res_start contains the atom
    res_start(-20:maxres),& ! number of the first atom in the residue. It is
                        ! only used inside MOPAC - it is never printed out
                        !
    start_res(-20:maxres), &   !
                        !
    bbone(3,-20:maxres) = 0! Atom numbers of the backbone atoms in a polypeptide

   double precision :: &
    angles(3, -20:maxres) !  Phi, Psi, and Omega for backbone angles (Ramachandran)


   character (len=4), dimension (-999:maxres) :: &
     allres        ! Names of all the residues
   data (allres(i), i = -999,0) /1000*"    "/
   data (allres(i), i = 1,maxres) /maxres*"    "/

   integer, parameter :: size_mres = 23
   character :: &
     tyres(size_mres)*3, tyr(size_mres), start_letter(500), allr(-999:maxres), ch*2

   character (len=40), dimension (11) :: &
     txeno         ! Names of xeno groups

!
!   Arrays that hold the localized molecular orbitals
!
!                               Occupied set
!
  integer ::     &
    cocc_dim,    & ! Size of the array cocc
    icocc_dim      ! Size of the array icocc
  integer, dimension (:), allocatable :: &
    ncf,         & ! Number of atoms involved in the LMO
    nncf,        & !  Starting address of the atom numbers of the atoms in the LMO
                   !  nncf(1) = 0
    icocc,       & ! Atom numbers of the atoms in the LMO's
    ncocc          ! Starting address of te atomic orbital coefficients in each LMO
  double precision, dimension (:), allocatable :: &
    cocc           ! Atomic orbital coefficients of the LMO's
!
!                                Virtual set
!
  integer ::     &
    cvir_dim,    & ! Size of the array cvir
    icvir_dim      ! Size of the array icvir
  integer, dimension (:), allocatable :: &
    nce,         & ! Number of atoms involved in the LMO
    nnce,        & ! Starting address of the atom numbers of the atoms in the LMO
                   ! nncf(1) = 0
    icvir,       & ! Atom numbers of the atoms in the LMO's
    ncvir          ! Starting address of te atomic orbital coefficients in each LMO
  double precision, dimension (:), allocatable :: &
    cvir           ! Atomic orbital coefficients of the LMO's
!
!                            Data for SCF and diagonalization
!
  integer, dimension (:,:), allocatable :: &
    ifmo           !
  integer, dimension (:), allocatable :: &
    idiag,       & !

    nfmo           ! Number of filled LMO's that interact with a virtual LMO
  double precision, dimension (:), allocatable :: &
    partf,       & !
    p1,          & !
    p2,          & !
    p3,          & !
    fmo,         & !
    partp,       & !
    parth,       & !
    ws             !

!
  integer ::     &
    nelred,      & ! Number of electrons to be used in the SCF during a geometry optimization.
                   ! In the first SCF, numred = nelecs
    norred,      & ! Number of atomic orbitals to be used in the SCF during a geometry optimization
                   ! In the first SCF, norred = norbs
    numred,      & ! Number of atoms to be used in the SCF during a geometry optimization
                   ! In the first SCF, numred = numat
    fmo_dim,     & ! Size of arrays holding the occupied M.O. - virtual M.O. interactions
    mode,        & !  0 for a simple MOZYME,
                   ! -1 if depleted arrays are to be constructed (RAPID only)
                   ! +1 if arrays are to be built using depleted arrays (RAPID only)
    ijc
  double precision, allocatable, dimension (:,:) ::  &
     part_dxyz
  double precision :: &
    cutofs,      & ! Overlap cut-off distance (normally 7 A)
    thresh,      & ! Criterion for deciding to do Euler rotation to annihilate
                   ! LMO occupied - virtual interaction
    energy_diff, &
    tiny,        &
    sumt,        &
    sumb,        &
    shift,       &
    pmax,        &
    ovmax,       &
    scfref,      &
    refnuc,      &
    refout,      &
    unused
  logical :: &
    rapid,       &      !  True if RAPID technique to be used
    odd_h,       &      !  Control print of banner for quentionable number of hydrogen atoms
    lstart_res(-20:maxres), &  !
    use_three_point_extrap
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!                    Delete or move everything below here
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    logical ::     &
      direct,      &
      lredop,      &
      semidr

    integer, allocatable, dimension (:) :: &
      isort,       &
      kopt,        &
      iopt
    data tyres / "GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", &
         & "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", &
         & "TYR", "MET", "PRO", "PRO", "PRO", "UNK" /
    data tyr / "G", "A", "V", "L", "I", "S", "T", "D", "N", "K", "E", "Q", &
         & "R", "H", "F", "C", "W", "Y", "M", "P", "P", "P", "?" /
    character*4 :: atomname(20,14), afn(20)
    integer :: n_add(20)
!
!  Names of all non-hydrogen atoms in each residue
!
  data (atomname(1,i),i=1,4)/" N  "," CA "," C  "," O  "/                                                                          !   GLY
  data (atomname(2,i),i=1,5)/" N  "," CA "," C  "," O  "," CB "/                                                                   !   ALA
  data (atomname(3,i),i=1,6)/" N  "," CA "," C  "," O  "," CB "," OG "/                                                            !   SER
  data (atomname(4,i),i=1,6)/" N  "," CA "," C  "," O  "," CB "," SG "/                                                            !   CYS
  data (atomname(5,i),i=1,7)/" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"/                                                     !   VAL
  data (atomname(6,i),i=1,7)/" N  "," CA "," C  "," O  "," CB "," OG1"," CG2"/                                                     !   THR
  data (atomname(7,i),i=1,8)/" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1"/                                              !   ILE
  data (atomname(8,i),i=1,7)/" N  "," CA "," C  "," O  "," CB "," CG "," CD "/                                                     !   PRO
  data (atomname(9,i),i=1,8)/" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE "/                                              !   MET
  data (atomname(10,i),i=1,8)/" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2"/                                             !   ASP
  data (atomname(11,i),i=1,8)/" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2"/                                             !   ASN
  data (atomname(12,i),i=1,8)/" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"/                                             !   LEU
  data (atomname(13,i),i=1,9)/" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ "/                                      !   LYS
  data (atomname(14,i),i=1,9)/" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2"/                                      !   GLU
  data (atomname(15,i),i=1,9)/" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2"/                                      !   GLN
  data (atomname(16,i),i=1,11)/" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2"/                       !   ARG
  data (atomname(17,i),i=1,10)/" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2"/                              !   HIS
  data (atomname(18,i),i=1,11)/" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "/                       !   PHE
  data (atomname(19,i),i=1,12)/" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH "/                !   TYR
  data (atomname(20,i),i=1,14)/" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"/  !   TRP
!
! "long" names of each residue
!
  data afn/'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO','MET','ASP','ASN','LEU','LYS','GLU','GLN','ARG','HIS','PHE','TYR','TRP'/
!
!  number of non-hydrogen atoms in each residue
!
  data n_add /4,    5,    6,    6,    7,    7,    8,    7,    8,    8,    8,    8,    9,    9,    9,   11,   10,   11,   12,   14/
    save
!
!  All arrays used in "big_swap" - these are the arrays that define the LMOs, the geometry, and some details of the SCF.
!
  integer, allocatable :: &
    nbonds_1(:), ibonds_1(:,:), icocc_1(:), icvir_1(:), ncocc_1(:), ncvir_1(:), &
    nncf_1(:), nncv_1(:), nnce_1(:), nce_1(:), ncf_1(:),  &
    nijbo_1(:,:), iijj_1(:), iij_1(:), ij_dim_1(:), ijall_1(:), numij_1(:), iorbs_1(:), &
    nbonds_2(:), ibonds_2(:,:), icocc_2(:), icvir_2(:), ncocc_2(:), ncvir_2(:), &
    nncf_2(:), nncv_2(:), nnce_2(:), nce_2(:), ncf_2(:),  &
    nijbo_2(:,:), iijj_2(:), iij_2(:), ij_dim_2(:), ijall_2(:), numij_2(:), iorbs_2(:)
    double precision, allocatable :: geo_1(:,:), geo_2(:,:), &
    cocc_1(:), cvir_1(:), cocc_2(:), cvir_2(:), xparam_1(:), xparam_2(:), &
    partf_1(:), partp_1(:), f_1(:), p_1(:), pold_1(:), dxyz_1(:), dxyz_2(:), &
    partf_2(:), partp_2(:), f_2(:), p_2(:), pold_2(:), parth_1(:), parth_2(:)
end module MOZYME_C
