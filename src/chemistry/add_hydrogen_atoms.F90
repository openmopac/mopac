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
  subroutine add_hydrogen_atoms()
!
!   Adds hydrogen atoms.
!
!   Intended for use in converting PDB files into files suitable for running in MOPAC
!
  use molkst_C,  only : numat, natoms, keywrd, maxatoms, moperr, refkey, line, id
  use common_arrays_C,  only : nat,  coord,  nbonds,  ibonds, txtatm, labels, atmass, geo, &
    lopt, na, l_atom, loc, xparam, nfirst, nlast, txtatm1, coorda, tvec
  use parameters_C, only : ams
  use MOZYME_C, only : atomname, afn, n_add
  use atomradii_C, only: atom_radius_covalent, is_metal
  use funcon_C,  only : pi
  use elemts_C, only: elemnt
  implicit none
  integer :: i,  j,  k, l, ii, nH, icc, store_numat, nb_icc, nc_icc, nd_icc, nmetals,  i_add(20,14), b_add(20,14), &
    an, store_nmetals
  character*4 :: labl
  double precision :: const, bond_length, angle, dihedral, internal_dihedral, store_atom_radius_covalent(107), &
  tetrahedral
  logical :: ionized = .false., l_reseq
  integer, allocatable :: store_nat(:), store_labels(:), store_lopt(:,:), hybrid(:), metals(:)
  logical, allocatable :: store_l_atom(:)
  double precision, allocatable :: store_coord(:,:), store_atmass(:)
  character, allocatable :: store_txtatm(:)*27, store_txtatm1(:)*27
  double precision, external :: distance
  integer, external :: nheavy
!
!  Number of hydrogen atoms on each atom in a residue
!
  data (i_add(1,i),i=1,4)/1,2,0,0/                        !   GLY
  data (i_add(2,i),i=1,5)/1,1,0,0,3/                      !   ALA
  data (i_add(3,i),i=1,6)/1,1,0,0,2,1/                    !   SER
  data (i_add(4,i),i=1,6)/1,1,0,0,2,1/                    !   CYS
  data (i_add(5,i),i=1,7)/1,1,0,0,1,3,3/                  !   VAL
  data (i_add(6,i),i=1,7)/1,1,0,0,1,1,3/                  !   THR
  data (i_add(7,i),i=1,8)/1,1,0,0,1,2,3,3/                !   ILE
  data (i_add(8,i),i=1,7)/0,1,0,0,2,2,2/                  !   PRO
  data (i_add(9,i),i=1,8)/1,1,0,0,2,2,0,3/                !   MET
  data (i_add(10,i),i=1,8)/1,1,0,0,2,0,0,1/               !   ASP
  data (i_add(11,i),i=1,8)/1,1,0,0,2,0,0,2/               !   ASN
  data (i_add(12,i),i=1,8)/1,1,0,0,2,1,3,3/               !   LEU
  data (i_add(13,i),i=1,9)/1,1,0,0,2,2,2,2,2/             !   LYS
  data (i_add(14,i),i=1,9)/1,1,0,0,2,2,0,0,1/             !   GLU
  data (i_add(15,i),i=1,9)/1,1,0,0,2,2,0,0,2/             !   GLN
  data (i_add(16,i),i=1,11)/1,1,0,0,2,2,2,1,0,2,1/        !   ARG
  data (i_add(17,i),i=1,10)/1,1,0,0,2,0,1,1,1,0/          !   HIS
  data (i_add(18,i),i=1,11)/1,1,0,0,2,0,1,1,1,1,1/        !   PHE
  data (i_add(19,i),i=1,12)/1,1,0,0,2,0,1,1,1,1,0,1/      !   TYR
  data (i_add(20,i),i=1,14)/1,1,0,0,2,0,1,0,1,0,1,1,1,1/  !   TRP
!
! Number of bonds (if less than minimum, this is a non-standard residue)
!
!    First four atoms: N CA C O
!
  data (b_add(1,i),i=1,4)/2,2,3,1/                        !   GLY
  data (b_add(2,i),i=1,5)/2,3,3,1,1/                      !   ALA
  data (b_add(3,i),i=1,6)/2,3,3,1,2,1/                    !   SER
  data (b_add(4,i),i=1,6)/2,3,3,1,2,1/                    !   CYS
  data (b_add(5,i),i=1,7)/2,3,3,1,3,1,1/                  !   VAL
  data (b_add(6,i),i=1,7)/2,3,3,1,3,1,1/                  !   THR
  data (b_add(7,i),i=1,8)/2,3,3,1,3,2,1,1/                !   ILE
  data (b_add(8,i),i=1,7)/2,3,3,1,2,2,2/                  !   PRO
  data (b_add(9,i),i=1,8)/2,3,3,1,2,2,2,1/                !   MET
  data (b_add(10,i),i=1,8)/2,3,3,1,2,3,1,1/               !   ASP
  data (b_add(11,i),i=1,8)/2,3,3,1,2,3,1,1/               !   ASN
  data (b_add(12,i),i=1,8)/2,3,3,1,2,3,1,1/               !   LEU
  data (b_add(13,i),i=1,9)/2,3,3,1,2,2,2,2,1/             !   LYS
  data (b_add(14,i),i=1,9)/2,3,3,1,2,2,2,1,1/             !   GLU
  data (b_add(15,i),i=1,9)/2,3,3,1,2,2,3,1,1/             !   GLN
  data (b_add(16,i),i=1,11)/2,3,3,1,2,2,2,2,3,1,1/        !   ARG
  data (b_add(17,i),i=1,10)/2,3,3,1,2,3,2,2,2,2/          !   HIS
  data (b_add(18,i),i=1,11)/2,3,3,1,2,3,2,2,2,2,2/        !   PHE
  data (b_add(19,i),i=1,12)/2,3,3,1,2,3,2,2,2,2,3,1/      !   TYR
  data (b_add(20,i),i=1,14)/2,3,3,1,2,3,2,3,2,3,2,2,2,2/  !   TRP
!
  tetrahedral = 109.4712206d0
!
!  Check for type of hydrogenation:  neutral or ionized
!
    if (index(keywrd, " IONI") /= 0) then
!
!  Set residues for ionized form
!
      i_add(13,9)  = 3                        ! Lysine nitrogen (+)
      i_add(16,11) = 2                        ! Arginine nitrogen (+)
      i_add(17,10) = 1                        ! Histidine nitrogen (+)
      i_add(10,8)  = 0                        ! Aspartate oxygen (-)
      i_add(14,9)  = 0                        ! Glutamate oxygen (-)
      ionized = .true.
    end if
    allocate (store_nat(natoms), store_coord(3,natoms))
    allocate (store_labels(natoms), store_atmass(natoms), store_lopt(3,natoms))
    allocate (store_l_atom(natoms), metals(natoms), store_txtatm(natoms), store_txtatm1(natoms))
    an = 0
    store_nmetals = 0
!
! Remove any hydrogen atoms present
!
    j = 0
    do i = 1, natoms
      if (nat(i) /= 1) then
        j = j + 1
        store_nat(j) = nat(i)
        store_coord(:,j) = coord(:,i)
        store_atmass(j)  = atmass(i)
        store_lopt(:,j)  = lopt(:,i)
        store_l_atom(j)  = l_atom(i)
        store_txtatm(j)  = txtatm(i)
        store_txtatm1(j) = txtatm1(i)
        txtatm1(j) = txtatm1(i)
        txtatm(j) = txtatm(i)
        coorda(:,j) = coorda(:,i)
      end if
    end do
    if (i - 1 /= j) call reset_breaks()
    if (moperr) return
    k = 0
    do i = 1, natoms
      if (labels(i) /= 99 .and. labels(i) /= 1) then
        k = k + 1
        store_labels(k) = labels(i)
      end if
    end do
    numat = j - id
    natoms = j
    nmetals = 0
    store_atom_radius_covalent = atom_radius_covalent
    do i = 1, numat
      do j = 1, 102
        if (is_metal(store_nat(i))) exit
      end do
      if (j < 103) then
        nmetals = nmetals + 1
        metals(nmetals) = i
        atom_radius_covalent(store_nat(i)) = -1.d0
      end if
    end do
!
!  Make extra storage room for hydrogen atoms
!
    deallocate (nat, txtatm, txtatm1, coord, nbonds, ibonds, labels, atmass, geo, lopt, na, l_atom, loc, xparam)
    deallocate (nfirst, nlast)
    maxatoms = numat*5
    allocate (nat(maxatoms), txtatm(maxatoms), txtatm1(maxatoms), coord(3,maxatoms), atmass(maxatoms))
    allocate (nbonds(maxatoms), ibonds(15,maxatoms), labels(maxatoms), geo(3,maxatoms))
    allocate (lopt(3,maxatoms), na(maxatoms), l_atom(maxatoms), loc(2,3*maxatoms), xparam(3*maxatoms))
    allocate (nfirst(maxatoms), nlast(maxatoms))
    na = 0
    nat(:numat) = store_nat(:numat)
    txtatm(:numat) = store_txtatm(:numat)
    txtatm1 = " "
    txtatm1(:numat) = store_txtatm1(:numat)
    coord(:,:numat) = store_coord(:,:numat)
    labels(:natoms) = store_labels(:natoms)
    atmass(:numat) = store_atmass(:numat)
    lopt(:,:numat) = store_lopt(:,:numat)
    l_atom(:numat) = store_l_atom(:numat)
    deallocate (store_nat, store_coord, store_labels, store_atmass, store_l_atom, &
    store_txtatm, store_txtatm1)
    store_numat = numat
!
!  Re-calculate connectivity.  This is simpler than trying to remove hydrogen atoms from nbonds and ibonds.
!
     l_reseq = (index(keywrd, " RESEQ") /= 0)
     if (l_reseq) call l_control("RESEQ", len("RESEQ"), -1)
     call geochk()                               ! This needs to be called, to work out residue names
     if (l_reseq) call l_control("RESEQ", len("RESEQ"), 1)
     if (moperr) return
     line = trim(refkey(1))
     call upcase(line, len_trim(line))
     i = index(line, " CVB")
     if (i > 0) then
       j = index(line(i:), ") ") + i
       line = line(i + 1: j - 1)
       ii = len_trim(line)
!
!    If everything is in quotation marks, use the CVB keyword
!
       k = 0
       do
         k = k + 1
         if (k >= ii) exit
         if (line(k:k) == '"') then
           do
             k = k + 1
             if (line(k:k) == '"') exit
           end do
         end if
         if (line(k:k) >= "0" .and. line(k:k) <= "9") exit
       end do
       if (k > ii - 3) &
       call l_control(line, ii, 1)
     end if
     call lewis(.false.)                         ! Needed again, because geochk breaks S-S bonds
     call check_CVS(.false.)                     ! and changes to bonding might be made by CVB
!
!   "hybrid" is intended for a two-pass method for assigning hydrogen atoms.
!   In the first pass, all atoms are assigned a hybridization value
!   In the second pass, atoms are assigned a hybrid value, based on their environment and their neighbor hybrid values,
!   and any hydrogen atoms necessary are assigned.
!
!
    allocate (hybrid(numat))
    hybrid = 0
    dihedral = 0.d0
    internal_dihedral = 0.d0
    do icc = 1, numat
      call h_type(22, icc, ionized, nH, nb_icc, nc_icc, nd_icc, bond_length, angle, dihedral, internal_dihedral, &
        hybrid, metals, nmetals)
    end do
    do icc = 1, numat
      bond_length = 0.d0
      angle = 0.d0
      do j = 1, 20
        if (afn(j) == txtatm(icc)(18:20)) exit
      end do
      if (j < 21) then
        do an = 1, n_add(j)
          labl = elemnt(nat(icc))//txtatm(icc)(15:16)
          if (labl(3:4) == "XT") labl = labl(1:2)
          if (labl == atomname(j,an)) exit
        end do
!
!   Check that the number of bonds is correct.  For peptide "N" and "C" allow for terminal group.
!
        if (nbonds(icc) /= b_add(j,an) .or. &
        (an == 1 .and. nbonds(icc) == 1) .or. &
        (an == 3 .and. nbonds(icc) == 2)) then
          j = 21
        end if
      end if
      if (j == 11 .and. txtatm(icc)(14:15) == "ND") j = 21 ! force planarity of primary amide in Asn
      if (j == 15 .and. txtatm(icc)(14:15) == "NE") j = 21 ! force planarity of primary amide in Gln
      if (j == 17 .and. (txtatm(icc)(14:15) == "ND" .or. txtatm(icc)(14:15) == "NE")) j = 21  ! His' imidazole N are too hard - solve directly
      if (nat(icc) == 8 .and. nbonds(icc) == 1) then    ! Carboxylic acid
    cooh: do i = 1, nbonds(icc)
          k = ibonds(i,icc)
          if (nat(k) == 6) then
            if (nbonds(k) == 3) then
              l = 0
              do ii = 1, nbonds(k)
                if (nat(ibonds(ii,k)) == 8) l = l + 1
              end do
              if (l == 2) then
                j = 21
                exit cooh
              end if
            end if
          end if
        end do cooh
      end if
      if (nat(icc) == 7 .and. txtatm(icc)(14:20) /= "NE  ARG") then                   !  Arginine
       arg: do i = 1, nbonds(icc)
          k = ibonds(i,icc)
          if (nat(k) == 6) then
            if (nbonds(k) == 3) then
              l = 0
              do ii = 1, nbonds(k)
                if (nat(ibonds(ii,k)) == 7) l = l + 1
              end do
              if (l == 3) then
                j = 21
                exit arg
              end if
            end if
          end if
        end do arg
      end if
!
! Atom icc is a hetero-atom, icc.e., not one of the 20 common residues
!
       call h_type(j, icc, ionized, nH, nb_icc, nc_icc, nd_icc, bond_length, angle, dihedral, internal_dihedral, &
        hybrid, metals, nmetals)
       if (j < 21) then
         bond_length = 0.d0
!
!  For atom icc, residue name is afn(j)
!
!
!  Check if it is a special case, e.g., a terminus or a hetero group
!
        if ( an > n_add(j)) then
          if (labl == " OXT") then
            if (ionized) then
              nH = 0  !  Terminus: -COO(-)
            else
              nH = 1  !  Terminus: -COOH
            end if
          else
            nH = 0
          end if
        else
!
!  For atom icc, atom name is atomname(an)
!
!  Number of hydrogen atoms to add: i_add(j,an)
!
!
!   Get number of hydrogen atoms to be added from amino-acid residue
!
          nH = i_add(j,an)
!
! if "an" is 4, then this is a peptide oxygen atom.
! check to see if it is a terminal atom, i.e., it is ...-(CR)-C-O
!
          if (an == 4) then
            if (nheavy(ibonds(1,icc)) == 2) then
              if (distance(icc, ibonds(1,icc)) > 1.35d0) nH = 1
            end if
          end if
          if (nat(icc) == 8 .and. nbonds(icc) > 1) nH = 0
          if (an == 1 .and. nbonds(icc) == 1) then
            if (ionized) then
              nH = 3  !  Terminus: make -NH2 -NH3(+)
            else
              nH = 2  !  Terminus: make -NH -NH2
            end if
          end if
!
!  Special cases where the default number of hydrogen atoms needs to be changed
!
          if (              an == 3 .and. nbonds(icc) == 2) nH = 1  !  Terminus: make -CO -CHO
          if (j == 3  .and. an == 6 .and. nbonds(icc)  > 1) nH = 0
          if (j == 4  .and. an == 6 .and. nbonds(icc)  > 1) nH = 0
          if (j == 8  .and. an == 1 .and. nbonds(icc) == 2) nH = 1  ! Pro terminus, make it >NH
          if (j == 11 .and. an == 8 .and. nbonds(icc) == 2) nH = 1
          if (nat(icc) == 6 .and. an >  4 .and. nbonds(icc) == 1) nH = 3  ! Convert truncated atom into methyl
        end if
        select case (nh)
          case (1)
            if (nd_icc == 0) then
              angle = 120.d0
              dihedral = 180.d0
            else
              bond_length = 1.09d0
            end if
          case (2)
            angle = tetrahedral
            dihedral = 120.d0
            internal_dihedral = 120.d0
          case (3)
            angle = tetrahedral
            dihedral = 120.d0
            internal_dihedral = 120.d0
        end select
        if (nat(icc) == 6) then
          bond_length = 1.09d0
        else if (nat(icc) == 7) then
          bond_length = 1.01d0
        else if (nat(icc) == 16) then
          bond_length = 1.34d0
        else
          bond_length = 0.99d0
        end if
      end if
      if (abs(bond_length) < 0.1d0) cycle
      angle = angle*pi/180.d0
      dihedral = dihedral*pi/180.d0
      internal_dihedral = internal_dihedral*pi/180.d0
      if (j < 21 .and. nat(icc) /= 8 .and. nat(icc) /= 16) then
        store_nmetals = nmetals
        nmetals = 0
      end if
      select case (nH)
        case (1)
          if (nc_icc /= 0 .and. nd_icc /= 0) then
            call add_a_sp3_hydrogen_atom(icc, nb_icc, nc_icc, nd_icc, bond_length, metals, nmetals)
          else
            call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals)
            if (moperr) return
!
!  Check to see if the system is water
!
            if (nat(icc) == 8 .and. nbonds(icc) == 2) then
              if (nat(ibonds(1,icc)) == 1 .and. nat(ibonds(2,icc)) == 1) then
!
! System is water.  Oxygen atom number is "icc", H1 is ibonds(1,icc), and H2 is ibonds(2,icc)
!
                i = ibonds(1,icc)
                j = ibonds(2,icc)
                call orient_water(icc, i, j, const)
                const = 1.d0 - 0.98d0/distance(icc, i)
                coord(:, i) = coord(:, i) - const*(coord(:, i) - coord(:, icc))
                const = 1.d0 - 0.98d0/distance(icc, j)
                coord(:, j) = coord(:, j) - const*(coord(:, j) - coord(:, icc))
              end if
            end if
          end if
        case (2)
          call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals)
          internal_dihedral = internal_dihedral + dihedral
          call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, internal_dihedral, metals, nmetals)
        case (3)
          call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, pi, metals, nmetals)
          call add_a_generic_hydrogen_atom(icc, nb_icc, numat, bond_length, angle, internal_dihedral, metals, nmetals)
          call add_a_generic_hydrogen_atom(icc, nb_icc, numat, bond_length, angle, internal_dihedral, metals, nmetals)
      end select
      if (j < 21 .and. nat(icc) /= 8 .and. nat(icc) /= 16) then
        nmetals = store_nmetals
      end if
    end do
!
!  Set all new atoms to hydrogen
!
    do i = store_numat + 1, numat
      nat(i) = 1
      labels(i) = 1
      atmass(i) = ams(1)
      j = ibonds(1,i)
      write(txtatm(i),'(a,i5,a)') txtatm(j)(:6),i,"  H  "//txtatm(j)(17:)
    end do
!
!  and fill arrays that might be used by PDBOUT
!
    natoms = max(numat + id, natoms)
    labels(numat + 1: numat + id) = 107
    geo(:,:numat) = coord(:,:numat)
    lopt(:,store_numat + 1: numat) = 1
    do i = 1, id
      lopt(:, numat + id) = store_lopt(:,store_numat + id)
    end do
    deallocate (store_lopt)
    l_atom(store_numat + 1: numat + id) = .true.
    na(:numat) = 0
    atom_radius_covalent = store_atom_radius_covalent
!
! At this point it's necessary to re-calculate nbonds and ibonds
!
    call set_up_dentate()
    if (index(keywrd, "SITE") /= 0) call geochk()
    call l_control("NOSITE", len("NOSITE"), 1)
    call reset_breaks()
    do i = numat + 1, numat + id
      labels(i) = 107
      nat(i) = 107
      geo(:,i)  = tvec(:, i - numat)
      coord(:,i)  = tvec(:, i - numat)
    end do
    call delete_ref_key("ADD-H", len_trim("ADD-H"), ' ', 1)
    return
  end subroutine add_hydrogen_atoms
  subroutine h_type(jjj, icc, ionized, nH, nb_icc, nc_icc, nd_icc, bond_length, angle, dihedral, internal_dihedral, &
    hybrid, metals, nmetals)
!
!  Work out the type and number of hydrogen atom(s) attached to atom icc
!  Work out bond-lengths and angles.
!
! bond_length:       The distance, in Angstroms, from atom icc to the hydrogen atom
!
! angle:             The angle between the hydrogen atom to be added and atoms icc and nb_icc
!
! dihedral :         The torsion angle between the hydrogen atom to be added and atoms icc,
!                    nb_icc, and nc_icc
!
! internal_dihedral: The torsion angle between the hydrogen atom to be added and atoms
!                    icc, nb_icc, and the previous hydrogen atom added, when two or three hydrogen atoms are added
!
  use molkst_C,  only : numat
  use funcon_C,  only : pi
  use common_arrays_C,  only : nat,  coord,  nbonds,  ibonds, txtatm, labels
  implicit none
  integer, intent (in) :: jjj, icc, nmetals, metals(nmetals)
  integer, intent (out) :: nb_icc, nc_icc, nd_icc
  integer, intent (out) :: nH
  integer, intent (inout) :: hybrid(numat)
  logical, intent (in) :: ionized
  double precision, intent (out) :: bond_length, angle, dihedral, internal_dihedral
!
  integer :: i, j, k, l, kk, ll, ii, jj, search(3,6)
  double precision :: Rab, Rac, Rad, sum, bits, tetrahedral = 109.4712206d0
  double precision, external :: distance
!
!   CN: distance half way between C-N and C=N
!   CC: distance half way between C-C and C=C
!   CO: distance half way between C-O and C=O
  double precision, parameter :: cyanide = 1.3d0, CC = 1.46d0, CN = 1.405d0, CO = 1.33d0
  logical :: l_tmp
  logical, external :: aromatic, aromatic_5, near_a_metal, guanidine
  data search /1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1/
  save :: tetrahedral
!
!   If hybrid(i) = 0, then hybridization unknown
!   If hybrid(i) = 1, then hybridization is sp
!   If hybrid(i) = 2, then hybridization is sp2
!   If hybrid(i) = 3, then hybridization is sp3
!
    bond_length = 0.d0
    angle = 0.d0       !  Not necessary, but useful when debugging
    nb_icc = 0
    nc_icc = 0
    nd_icc = 0
    nH = 0
    Rab = 0.d0
    Rac = 0.d0
    Rad = 0.d0
    if (nbonds(icc) == 0) then
!
!                                         Atom is isolated
!
      if (hybrid(icc) == 0) then
        hybrid(icc) = 3
        return
      end if
      select case (nat(icc))
        case (6)
          nH = 3                        ! Methane
          bond_length = 1.09d0
          angle = tetrahedral
          dihedral = 120.d0
          internal_dihedral = 120.d0
        case (7)
          nH = 2                        ! Ammonia
          bond_length = 1.0d0
          angle = 107.d0
          dihedral = 114.5d0
          internal_dihedral = dihedral
        case (8)
          nH = 1                        ! Water
          bond_length = 0.957d0
          angle = 104.4d0
        case (14)
          nH = 3                        ! Silane
          bond_length = 1.48d0
          angle = tetrahedral
          dihedral = 120.d0
          internal_dihedral = 120.d0
        case (15)
          nH = 2                        ! Phosphine
          bond_length = 1.42d0
          angle = 93.5d0
          dihedral = 93.8d0
          internal_dihedral = dihedral
        case (16)
          nH = 1                        ! H2S
          bond_length = 1.336d0
          angle = 92.1d0
        case (9, 17, 35, 53)
          return                        ! Halogen ions
        case (33)
          nH = 2                        ! Arsine
          bond_length = 1.519d0
          angle = 91.8d0
          dihedral = 92d0
          internal_dihedral = dihedral
        case (34)
          nH = 1                        ! H2Se
          bond_length = 1.46d0
          angle = 91.d0
        case (51)
          nH = 2                        ! Stibine
          bond_length = 1.707d0
          angle = 91.8d0
          dihedral = 92d0
          internal_dihedral = dihedral
        case (52)
          nH = 1                        ! H2Se
          bond_length = 1.69d0
          angle = 90.d0
        case default
          nH = 0                        ! If the atom is not specified, don't assign hydrogen atoms to it.
          return
      end select
!
!                                         Common to all isolated atoms that have hydrogen atoms
!
      numat = numat + 1
      nat(numat) = 1
      labels(numat) = 1
      do i = 1, 6
        coord(:, numat) = coord(:,icc) + bond_length*search(:,i)
        if ( .not. near_a_metal(icc, numat, metals, nmetals)) exit
      end do
      if (i > 6) then
        nH = 0
        if (nat(icc) == 8) then
          if (nbonds(icc) < -2) then
            numat = numat -1
            return
          end if
!
!  Water near to a metal atom.  Carefully add two hydrogen atoms
!
          coord(:, numat) = coord(:,icc) + bond_length*search(:,1)
          ibonds(1,icc) = numat
          ibonds(1,numat) = icc
          numat = numat + 1
          coord(:, numat) = coord(:,icc) + bond_length*search(:,2)
          ibonds(2,icc) = numat
          ibonds(1,numat) = icc
          i = ibonds(1,icc)
          j = ibonds(2,icc)
          call orient_water(icc, i, j, sum)
          sum = 1.d0 - 0.98d0/distance(icc, i)
          coord(:, i) = coord(:, i) - sum*(coord(:, i) - coord(:, icc))
          sum = 1.d0 - 0.98d0/distance(icc, j)
          coord(:, j) = coord(:, j) - sum*(coord(:, j) - coord(:, icc))
          return
        end if
        numat = numat - 1
        return
      end if
      nb_icc = numat
      nbonds(icc) = 1
      nbonds(numat) = 1
      ibonds(1,icc) = numat
      ibonds(1,numat) = icc
    else
      ii = ibonds(1,icc)
      if (nbonds(icc) == 3) then
        kk = ibonds(3,icc)
        jj = ibonds(2,icc)
!
!  Sort into ascending order
!
        if (nat(ii) > nat(jj)) then
          i = jj
          jj = ii
          ii = i
        end if
        if (nat(ii) > nat(kk)) then
          i = kk
          kk = ii
          ii = i
        end if
          if (nat(jj) > nat(kk)) then
            i = kk
            kk = jj
            jj = i
        end if
!
! Atomic numbers increase ii <= jj <= kk
!
        Rab = distance(icc, ii)
        Rac = distance(icc, jj)
        Rad = distance(icc, kk)
        if (nat(ii) == nat(jj)) then
          if (Rab > Rac) then
            i = ii
            ii = jj
            jj =i
            sum = Rab
            Rab = Rac
            Rac = sum
          end if
        end if
      else if (nbonds(icc) == 2) then
        jj = ibonds(2,icc)
        if (nat(jj) < nat(ii)) then
          i  = ii
          ii = jj
          jj = i
        end if
!
! Atomic numbers increase ii <= jj
!
        Rab = distance(icc, ii)
        Rac = distance(icc, jj)
        if (nat(ii) == 6 .and. nat(jj) == 6) then
          if (Rab > Rac) then   ! reverse order so Rab is smaller than Rac
            Rad = Rab
            Rab = Rac
            Rac = Rad
            i = ii
            ii = jj
            jj = i
          end if
        end if
      else
        Rab = distance(icc, ii)
      end if
!
!
!                   Now start to assign the number of hydrogen atoms to each non-hydrogen atom.
!                   The order of the assignment is as follows:
!                   First, select by atomic number of atoc "icc"
!                   Next, select by number of atoms attached to the atom "icc"
!                   Next, select by atomic number of atoms attached to atom attached to atom "icc",
!                     in the order: lowest atomic number first, then next, then highest.
!                   Then use bond-lengths to decide the bond order, and thus the number of hydrogen atoms
!
!                   Do high-reliability tests first, e.g, is an atom part of an aromatic ring, or attached to
!                   an aromatic ring?
!
!                   When editing this section, preserve the order of assignment, in particular
!                   the increasing order of atomic number.  If this order is lost, then debugging becomes hard.
!
!             "Typical" bond-lengths  ("T" means "triple bond")
!
!                     C-C 1.54     C=C 1.34  CTC 1.20
!                     C-N 1.43     C=N 1.38  CTN 1.16
!                     C-O 1.43     C=O 1.23  CTO 1.13
!                     C-S 1.83     C=S
!                     N-N 1.47     N=N 1.24  NTN 1.10
!                     N-O 1.36     N=O 1.22
!                     O-O 1.48     O=O 1.21
!                     S-S 2.08     S=S
!
!             ii = atom attached to icc
!             jj = atom attached to icc, and its atomic number is greater than or equal to that of ii
!             kk = atom attached to icc, and its atomic number is greater than or equal to that of jj
!             Rab = interatomic distance icc to ii
!             Rac = interatomic distance icc to jj
!             Rad = interatomic distance icc to kk
!
!
      bond_length = 1.09d0
      if (nbonds(icc) == 1) then
        nb_icc = ii
        do i = 1, nbonds(ii)
          kk = ibonds(i,ii)
          if (kk /= icc .and. nat(kk) > 1) nc_icc = kk
        end do
        if (nc_icc == 0) then
          do i = 1, nbonds(ii)
            kk = ibonds(i,ii)
            if (kk /= icc ) nc_icc = kk
          end do
        end if
      else if (nbonds(icc) == 2) then
        nb_icc = ii
        nc_icc = jj
      else
        nb_icc = ibonds(1,icc)
        nc_icc = ibonds(2,icc)
        nd_icc = ibonds(3,icc)
      end if
      if (jjj < 21) return
      select case (nat(icc))
        case (5)
          if (nbonds(icc) == 1) then
                nH = 3                !  sp3 boron
                bond_length = 1.07d0
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 3
          else if (nbonds(icc) == 2) then
                nH = 1                !  sp3 boron
                bond_length = 1.07d0
                angle = 120
                dihedral = 180.d0
                hybrid(icc) = 1
          end if
                                        ! Carbon is very important
        case (6)                        ! so there are a lot of options for C-H bonds
                                        ! (Add more as necessary)
          select case (nbonds(icc))
          case (1)                      ! Carbon bonded to one other atom
            if (nat(ii) == 6) then
!
!  C-C bonds are common. So look at single, double, and triple bonds
!
              jj = 0
              do i = 1, nbonds(ii)
                if (nat(ibonds(i,ii)) == 6) jj = jj + 1
              end do
              if (jj == 3) then
!
! See if atom icc is attached to is already sp3
!
                 jj = 0
                 kk = 0
                 do i = 1, nbonds(ii)
                   ll = ibonds(i,ii)
                   if (nat(ll) == 6 .and. jj == 0 .and. ll /= icc .and. (kk /= ll .or. kk == 0)) jj = ibonds(i,ii)
                   if (nat(ll) == 6 .and. kk == 0 .and. ll /= icc .and. (jj /= ll .or. jj == 0)) kk = ibonds(i,ii)
                 end do
                call dihed(coord, icc, ii, jj, kk, sum)
                sum = abs(sum)
                if (hybrid(ii) == 3 .or. min(pi - sum, sum) > 0.4d0) then  ! 22.9 degrees
                  nH = 3
                  angle = tetrahedral                   !  If flat, nc
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  return
                end if
              end if
              jj = 0
              do i = 1, nbonds(ii)
                kk = ibonds(i,ii)
                if (kk /= icc) jj = kk
              end do
              if (Rab > CC .or. hybrid(ii) == 3) then
                nH = 3  ! Methyl
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 3
                if (nc_icc == 0) then
                  do i = 1, nbonds(ii)
                    j = ibonds(i,ii)
                    if (j /= icc) then
                      nc_icc = j
                      return
                    end if
                  end do
                end if
              else if (Rab > 11.35d0) then
!
!  Why was this test here?  The original test was Rab > 1.35d0
!
                nH = 3                !  sp3 carbon
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 3
              else if (Rab > 1.25d0) then
                nH = 2                 ! ethylene
                bond_length = 1.08d0
                angle = 123.d0
                dihedral = 180.d0
                internal_dihedral = 180.d0
                hybrid(icc) = 2
              else
                nH = 1                  ! acetylene
                bond_length = 1.07d0
                angle = 179.d0
                hybrid(icc) = 1
              end if
            else if (nat(ii) == 7) then
              if (hybrid(ii) == 3 .or. Rab > CN) then
                nH = 3                ! Methyl
                bond_length = 1.09d0
                angle = tetrahedral
                internal_dihedral = 120.d0
                return
              end if
              hybrid(icc) = 1
              if (ionized) then
                nH = 0                    ! Cyanide ion
              else
                nH = 1                    ! Hydrogen cyanide
                angle = 179.d0
              end if
            else if (nat(ii) == 8) then
              if (Rab > 1.37d0 .or. nbonds(ii) > 1) then
                nH = 3                ! Methanol
                bond_length = 1.09d0
                angle = tetrahedral
                internal_dihedral = 120.d0
                hybrid(icc) = 3
              else
                nH = 2                ! Formaldehyde
                angle = 120.d0
                internal_dihedral = 180.d0
                hybrid(icc) = 2
              end if
            else
                nH = 3                ! -CH3 (catch-all)
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 3
            end if
!
          case (2)                      ! Carbon bonded to two other atoms
!
            if(nat(ii) == 6) then
              if (nat(jj) == 6) then
                if (aromatic_5(icc, ii, jj, hybrid)) then
                  if (Rab < 1.44d0) then
                    nH = 1                  ! C=C-C
                    angle = 120.d0
                    dihedral = 180.d0
                    hybrid(icc) = 2
                    return
                  else
                    nH = 2                  ! C-C-C
                    angle = tetrahedral
                    dihedral = 120.d0
                    internal_dihedral = 120.d0
                    hybrid(icc) = 3
                    return
                  end if
                end if
                if (aromatic(icc, ii, jj)) then
                  nH = 1                  ! Double bond
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                end if
                if (Rab < 1.27d0) then
                  nH = 0
                  hybrid(icc) = 1
                  return
                end if
                nH = 2
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 120.d0
                hybrid(icc) = 3
                if (hybrid(ii) == 3 .and. hybrid(jj) == 3) return
!
!   Carbon bonded to two other carbon atoms, e.g. C-C-C or C-C=C, etc.
!   Check for double bond
!
                call bangle (coord, ii, icc, jj, sum)    !
                if (sum > 3.d0 .and. Rac < 1.35d0) then ! Acetylenic
                 nH = 0 ! Acetylenic
                 hybrid(icc) = 2
                 return
                end if
                if ( Rab < 1.44d0) then
                  nH = 1                  ! Double bond
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                end if
                i = 0
                if (sum > 2.03d0 .and. (Rab < 1.39d0) .or. sum > 2.2d0 .and. (Rab < 1.41d0)) then
                  if (nbonds(ii) /= 4) then
                    if (nbonds(ii) == 3) then
                      j = 0
                      do k = 1,3
                        if (nat(ibonds(k,ii)) == 8) j = j + 1
                      end do
                      if (j /= 2) i = 1
                    else if (nbonds(ii) == 2) then
                      call bangle (coord, ibonds(1,ii), ii, ibonds(2,ii), sum)
                      if (sum > 2.d0) i = 1
                    else
                      i = 1 ! over 116.3 degrees and short bond-lengths
                    end if
                  end if
                end if
                if (i == 1) then
                  nH = 1                               !  sp2 carbon
                  angle = 120.d0
                  internal_dihedral = 180.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                end if
              else if (nat(jj) == 7) then
!
!   Carbon bonded to one carbon and one nitrogen atoms, e.g. C-C-N or C-C=N
                if (Rac < 1.2d0) then
                  hybrid(icc) = 1
                  return
!   Check for double bond
!
                 else if (hybrid(ii) == 3 .and. hybrid(jj) == 3) then
                  nH = 2                  ! C-C-N
                  bond_length = 1.09d0
                  angle = tetrahedral
                  internal_dihedral = 120.d0
                  dihedral = 120.d0
                  hybrid(icc) = 3
                else if (Rab < 1.43d0 .or. Rab + Rac < 2.9d0 .and. aromatic_5(icc, ii, jj, hybrid) .and. nbonds(jj) == 2) then
                  nH = 1                  ! Double bond
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                else if (Rab + Rac < 2.8d0 .and. aromatic_5(icc, ii, jj, hybrid)) then
                  nH = 1                  ! Double bond
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                else if (aromatic(icc,ii,jj) .or. Rac < CN) then
                  nH = 1                  ! Double bond
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                else
                  nH = 2                  ! C-C-N
                  bond_length = 1.09d0
                  angle = tetrahedral
                  internal_dihedral = 120.d0
                  dihedral = 120.d0
                  hybrid(icc) = 3
                end if
              else if (nat(jj) == 8) then
!
!   Carbon bonded to one carbon and one oxygen atom, e.g. C-C-O or C-C=O
!   Check for double bond
!
                if (aromatic_5(icc,ii,jj, hybrid)) then
                  if (Rab < 1.44d0) then
                    nH = 1                  ! C=C-O
                    angle = 120.d0
                    dihedral = 180.d0
                    hybrid(icc) = 2
                    return
                  else
                    nH = 2                  ! C-C-O
                    angle = tetrahedral
                    dihedral = 120.d0
                    internal_dihedral = 120.d0
                    hybrid(icc) = 3
                    return
                  end if
                end if
                if (Rab < 1.4d0 .or. (Rac < 1.3d0 .and. nbonds(jj) == 1)) then
                  if (Rab < 1.4d0 .and. (Rac < 1.3d0 .and. nbonds(jj) == 1)) return ! ketene
                  nH = 1                  ! C=C-O or C-C=O
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                end if
                sum = 0.d0
                j = 0
                do i = 1, nbonds(ii)
                  if (nat(ibonds(i,ii)) /= 1) j = j + 1
                end do
                if (nat(ii) == 6 .and. j == 2) then
                  jj = ibonds(1,ii)
                  if (jj == icc) jj = ibonds(2,ii)
                  call bangle (coord, icc, ii, jj, sum)
                end if
                if (Rac < 1.3d0 .and. sum > 2.03d0) then
                  nH = 1                  ! C-C=O
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                else
                  nH = 2                  ! C-C-O
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  return
                end if
              else if (nat(jj) == 9 .or. nat(jj) == 17 .or. nat(jj) == 35 .or. nat(jj) == 53) then
                if (Rab < 1.25d0) then
                  nH = 0                  ! CTC-X
                  hybrid(icc) = 1
                  return
                elseif (Rab < 1.44d0) then
                  nH = 1                  ! C=CH-X
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                else
                  nH = 2                  ! C-CH2-X
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  return
                end if
              else if (nat(jj) == 15) then
!
!   Carbon bonded to one carbon and one phosphorus atom, e.g. C-C-P
!   Check for double bond
!
                nH = 2                    ! C-C-P
                bond_length = 1.09d0
                angle = tetrahedral
                internal_dihedral = 120.d0
                dihedral = 120.d0
                hybrid(icc) = 3
              else if (nat(jj) == 16 .or. nat(jj) == 34) then
!
!   Carbon bonded to one carbon and one sulfur or selenium atom, e.g. C-C-S
!   Check for double bond
!
                if (Rab > 1.35d0) then
                  nH = 2
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  hybrid(icc) = 2
                else if (Rab > 1.25d0) then
                  if (Rab < 1.4d0 .and. (Rac < 1.6d0 .and. nbonds(jj) == 1)) return ! thioketene
                  nH = 1                  ! C=CH-S
                  bond_length = 1.07d0
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 1
                else
                  nH = 0                    ! CTC-S
                  hybrid(icc) = 3
                end if
              end if
            else if (nat(ii) == 7) then
              if (nat(jj) == 8) then
                if (Rac < 1.3d0) then
                  nH = 1                    ! aromatic_5 - type
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                else
                  nH = 2                    ! C-C-S
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  return
                end if
              else if (nat(jj) == 7) then
                if ((hybrid(ii) == 3 .or. hybrid(jj) == 3) .and. (Rab > 1.4d0 .and. Rac > 1.4d0)) then
                  nH = 2                    ! N-C-N
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                  return
                else
                  nH = 1                    ! N=C-N
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                  return
                end if
              else if (nat(jj) == 14) then  ! N-C-S
                if (Rab > 1.4d0 .and. Rac > 1.8d0) then
                  nH = 2                    ! C-C-S
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                else
                  nH = 1
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                end if
                return
              end if
              if (nat(jj) == 16) then
!
!   Carbon bonded to one nitrogen and one sulfur atom, e.g. N-C-S
!
                nH = 2                    ! N-C-S
                angle = tetrahedral
                dihedral = 120.d0
                internal_dihedral = 120.d0
                hybrid(icc) = 3
                return
              end if
              if (nbonds(icc) == 2 .and. Rab < cyanide) then
                nH = 0                    ! Cyanide (nitrile)
              else
                nH = 1                    ! aromatic_5 - type
                angle = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 2
              end if
            else if (nat(ii) == 8) then
              if (nat(jj) == 9 .or. nat(jj) == 17 .or. nat(jj) == 35 .or. nat(jj) == 53) then
                if (Rab > 1.3d0) then
                  nH = 2
                  angle = tetrahedral
                  dihedral = 120.d0
                  internal_dihedral = 120.d0
                  hybrid(icc) = 3
                else
                  nH = 1
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 2
                end if
                return
              end if
              if (nbonds(ii) == 2 .and. nbonds(jj) == 2) then
                nH = 2                    ! X-O-CH2-O-X
                angle = tetrahedral
                dihedral = 120.d0
                internal_dihedral = 120.d0
                hybrid(icc) = 3
                return
              end if
              call bangle (coord, ibonds(1,icc), icc, ibonds(2,icc), sum)
              if (nbonds(ii) /= 1 .or. nbonds(jj) /= 1 .or. &
                Rab > 1.28d0 .or. (nat(jj) == 8 .and. sum < 2.618d0) ) then
                nH = 1                    ! HCO2
                angle = 120.d0
                dihedral = 180.d0
                hybrid(icc) = 2
              else
                nH = 0
                hybrid(icc) = 2
              end if
              hybrid(icc) = 2
            else if (nat(ii) == 15) then
              nH = 0
              hybrid(icc) = 2
            else
              hybrid(icc) = 3
            end if
!
          case (3)                      ! Carbon bonded to three other atoms
!
            if (nat(jj) == 7 .and. Rac < CN .and. nat(kk) == 8) then
              nH = 0 ! C-(CO)=N already three ligands
              return
            end if
            if (nat(jj) == 6 .and. Rad < 1.31d0 .and. nat(kk) == 7) then
              nH = 0 ! C-(C)=N already three ligands
              return
            end if
            if (nat(ii) == 6 .and. Rab < CC .and. nat(jj) == 6 .and. nat(kk) == 16) then
              nH = 0 ! S-C=C and already three ligands
              return
            end if
            call dihed(coord, icc, ii, jj, kk, sum)
!
!  A dihedral angle involving the current atom and the three joined to it
!  can be used in deciding whether or not it should be sp2 or sp3
!
!  The tetrahedral dihedral angle, that is the torsion angle from an atom
!  to the first ligand, then to the second, then to the third, is ~35.26439 degrees.
!  The sine of this number is: 0.57735
!
            if (min(2*pi - sum, sum) < 0.05d0) then
              nH = 0    ! Very flat, so it must be sp2
              hybrid(icc) = 2
              return
            end if
            if (min(2*pi - sum, sum) < 0.15d0) then
              if (aromatic(icc,ii,jj) .or. aromatic(icc,ii,kk) .or. aromatic(icc,jj,kk)) then
                nH = 0                               ! aromatic ring
                hybrid(icc) = 2
                return
              end if
              if (aromatic_5(icc,ii,jj, hybrid) .or. aromatic_5(icc,ii,kk, hybrid) &
                .or. aromatic_5(icc,jj,kk, hybrid)) then
                nH = 0                               ! aromatic five membered ring
                hybrid(icc) = 2
                return
              end if
            end if
            if (min(2*pi - sum, sum) < 0.40d0) then
              bits = 0.d0
              if (nat(ii) == 7 .or. nat(ii) == 8) bits = bits + 0.1d0
              if (nat(jj) == 7 .or. nat(jj) == 8) bits = bits + 0.1d0
              if (nat(kk) == 7 .or. nat(kk) == 8) bits = bits + 0.1d0
              if (bits + Rab + Rac + Rad < 4.5d0) then
                nH = 0                      !  If flat, and short bond lengths, make sp2
                hybrid(icc) = 2
                return
              end if
            end if
            if (min(2*pi - sum, sum) > 0.25d0) then
              nH = 1                      !  If pyrimidal, make sp3, and don't check anything else
              hybrid(icc) = 3
              return
            end if
            nH = 1
            hybrid(icc) = 3
            if (nat(ii) == 6) then
              if (nat(jj) == 6) then
                if (nat(kk) == 6) then              !  carbon bonded to three carbon atoms
                  sum = min(Rab, Rac, Rad)
                  if (sum < 1.44d0) then
                    nH = 0                            ! Result is obvious
                    hybrid(icc) = 2
!
!  Consider using "call dihed (coord, icc, ii, jj, kk, bits)" if there is an ambiguity.
!
                  end if
                else if (nat(kk) == 7) then         !  carbon bonded to two carbon atoms and one nitrogen
                  if (hybrid(ii) == 3 .and. hybrid(jj) == 3 .and. hybrid(kk) == 3) return ! Obvious
                  j = 0
                  do i = 1, nbonds(kk)
                    if (nat(ibonds(i,kk)) /= 1) j = j + 1  ! Number of heavy atoms on nitrogen
                  end do
                  if (min(Rab, Rac) < 1.41d0 .or. Rad < 1.39d0 .and. j < 3) then
                    nH = 0
                    hybrid(icc) = 2
                  end if
                else if (nat(kk) == 8) then         !  carbon bonded to two carbon atoms and one oxygen
                  if (min(Rab, Rac) < 1.45d0 .or. Rad < 1.30d0) then
                    nH = 0
                    hybrid(icc) = 2
                  end if
                else if (nat(kk) == 9 .or. nat(kk) == 17 .or. nat(kk) == 35 .or. nat(kk) == 53) then
                  if (min(Rab, Rac) < 1.45d0) then     !  carbon bonded to two carbon atoms and one halogen
                    nH = 0
                    hybrid(icc) = 2
                  end if
                else if (nat(kk) == 16) then
                  sum = min(Rab, Rac)
                  if (aromatic(icc,ii,jj) .or. sum < 1.44d0) then
                    nH = 0                               ! aromatic ring
                    hybrid(icc) = 2
                    return
                  end if
                end if
              else if (nat(jj) == 7) then
                if (nat(kk) == 7) then
                  sum = min(Rab, Rac)
                  if (sum < 1.405d0 .or. hybrid(ii) == 2 .or. hybrid(jj) == 2 .or. hybrid(kk) == 2) then
                    nH = 0
                    hybrid(icc) = 2
                  end if
                else if (nat(kk) == 8) then
                  if (Rac > 1.45d0 .and. Rad > 1.28d0) then       !  Carbon bonded to carbon, sulfur and nitrogen
                    return
                  else
                    nH = 0                                        !  Most likely, a peptide
                    hybrid(icc) = 2
                  end if
                else if (nat(kk) == 16) then
                  call dihed(coord, icc, ii, jj, kk, sum)
                  if (min(2*pi - sum, sum) < 0.15d0) then
                    nH = 0                               ! aromatic ring
                    hybrid(icc) = 2
                    return
                  end if
                else
                  nH = 0                          !  Most likely, a peptide
                  hybrid(icc) = 2
                end if
              else if (nat(jj) == 8) then
                if (nat(kk) == 8 .or. Rac < 1.25d0) then
                  if (nbonds(jj) == 2 .and. nbonds(kk) == 2) then
                    if (Rab < 1.44d0) nH = 0
                    return
                  else
                    nH = 0
                    hybrid(icc) = 2
                  end if
                end if
              else if (nat(jj) == 9 .or. nat(jj) == 17 .or. nat(jj) == 35 .or. nat(jj) == 53) then
                if (Rab < 1.44d0) nH = 0  ! X-C=C
              end if
            else if (nat(ii) == 7) then
              if (nat(jj) == 7) then
                select case (nat(kk))
                case (7)
                  nH = 0
                  hybrid(icc) = 2
                case(8)
                  if (Rad < 1.33d0) then
                    nH = 0
                    hybrid(icc) = 2
                  end if
                case(16)
                  if (Rad < 1.7d0) then
                    nH = 0
                    hybrid(icc) = 2
                  end if
                end select
              end if
            else if (nat(ii) == 8 .or. nat(ii) == 16) then
              nH = 0                                ! Carbon bonded to three oxygen or sulfur atoms !
            end if
          case default
            hybrid(icc) = 3
          end select
                                        ! Nitrogen is very important
        case (7)                        ! so there are a lot of options for N-H bonds
                                        ! (Add more as necessary)
          if (nbonds(icc) == 1) then
            nH = 2                                ! -NH2
            bond_length = 1.02d0
            angle = tetrahedral
            dihedral = 120.d0
            internal_dihedral = 120.d0
            ii = ibonds(1,icc)
            if (nat(ii) == 6) then
              if (Rab < cyanide .and. nbonds(ii) < 3) then
                nH = 0
                hybrid(icc) = 1                          ! Cyanide
                return
              end if
              sum = 1.29d0
              if (nat(ii) == 6 .and. nbonds(ii) == 3) then
                do i = 1,3
                  if (nat(ibonds(i,ii)) == 8) sum = 1.1d0
                end do
              end if
              if (Rab < sum) then
                nH = 1                                   !  C=NH
                angle = 120.d0
                internal_dihedral = 180.d0
                dihedral = 0.d0
                hybrid(icc) = 1
                return
              end if
              ! primary amide group, -C(=O)-NH2
              if (nbonds(ii) == 3) then
                l_tmp = .false.
                do i = 1,3
                  if (nat(ibonds(i,ii)) == 8 .and. nbonds(ibonds(i,ii)) == 1) then
                    l_tmp = .true.
                    nc_icc = ibonds(i,ii)
                  end if
                end do
                if (l_tmp) then
                  angle = 120.d0
                  internal_dihedral = 180.d0
                  dihedral = 180.d0
                  return
                end if
              end if
!
!  For Arg, make one -N -NH2, then decide on the other later on.
!
              if (guanidine(icc, ionized, nH, angle, internal_dihedral, dihedral, hybrid)) then
                continue
                return
                end if
               do i = 1, nbonds(ii)
                jj = ibonds(i,ii)
                if (nat(jj) /= 6 .and. nat(jj) /= 7) exit
                if (nat(jj) == 7 .and. jj /= icc) then  !  N-C-N structure
                  k = 0
                  l = 0
                  do j = 1, nbonds(jj)
                    if (nat(ibonds(j,jj)) == 6) k = k + 1
                    if (nat(ibonds(j,jj)) == 1) l = l + 1
                  end do
                  if (k == 2) cycle               !  Exclude C-N-C
                  if (l == 2) then
                    hybrid(icc) = 2
                    if (ionized) then
                      nH = 2                      !   ionized and an Arg
                      angle = 120.d0
                      internal_dihedral = 180.d0
                      dihedral = 0.d0
                    else
                      nH = 1
                      angle = 120.d0
                      internal_dihedral = 180.d0
                      dihedral = 0.d0
                    end if
                    return                         !  Unconditional
                  else
                    nH = 2
                    angle = 120.d0
                    internal_dihedral = 180.d0
                    dihedral = 0.d0
                    hybrid(icc) = 3
                    return
                  end if
                end if
               end do
            else if (nat(ii) == 7) then
              nH = 0  ! Azide
              return
            end if
            if (ionized) then
              nH = 3  !  Ionized  (what about cyanide?)
              if (nat(ii) == 6) then
                do i = 1, nbonds(ii)
                  if (nat(ibonds(i,ii)) == 8) exit
                end do
                if (i <= nbonds(ii)) nH = 2          !  Exclude -CO-NH3
              end if
            else
              nH = 2                                 ! Generic C-NH2
              angle = 120.d0
              internal_dihedral = 120.d0
              dihedral = 120.d0
              hybrid(icc) = 3
            end if
          else if (nbonds(icc) == 2) then
!
!  Nitrogen bonded to two atoms
!
            if (nat(ii) == 7) then
              if (nat(jj) == 7) then
                sum = min(Rab, Rac)
                if (sum < 1.35d0 .or. nbonds(ii) + nbonds(jj) > 4) then
                  nH = 0
                  hybrid(icc) = 2
                else
                  nH = 1
                  angle = 120.d0
                  dihedral = 180.d0
                  hybrid(icc) = 3
                end if
                return
              end if
            end if
            if (nat(ii) == 6) then
              if (nat(jj) == 6) then
                if (Rab < 1.3d0) then
                  nH = 0
                  continue
                  return
                end if
                if (Rab < CN) then
!
!  Test for double bond
!
                  do i = 1, nbonds(ii)
                    kk = ibonds(i,ii)
                    if (kk /= icc) then
                      if (nat(kk) == 6 .and. distance(kk, ii) < CC .or. &
                        nat(kk) == 7 .or. &
                        nat(kk) == 8 .and. distance(kk, ii) < CO) exit
                    end if
                  end do
                  if (i > nbonds(ii)) then
                    nH = 0
                    continue
                    return
                  end if
                end if
              end if
              if (guanidine(icc, ionized, nH, angle, internal_dihedral, dihedral, hybrid)) return
              if (nat(jj) == 7) then
                if (Rab < 1.39d0 .or. Rac < 1.40d0) then
                  nH = 0
                  return
                end if
              end if
            end if
            if (nat(ii) == 8) then
              if (nat(jj) == 8) then
                if (Rab < 1.30d0 .or. Rac < 1.30d0) then
                  nH = 0
                  return
                end if
              end if
            end if
            if (aromatic(icc,ii,jj)) then
              nH = 0                               ! aromatic ring
              hybrid(icc) = 2
              return
            end if
!
!
!  Locate the N-C-N structure, kk = distant N
!
            kk = 0
            do i = 1, nbonds(ii)
              j = ibonds(i,ii)
              if (nat(j) == 7 .and. j /= icc) then
                kk = j
                if (nbonds(j) == 2) then
                  k = ii   ! N-C-N structure found.
                  ii = jj  ! ii is the C in N-C-N
                  jj = k   ! C(jj) - N(icc) - C(ii) - N(kk)
                end if
              end if

            end do
            if (kk == 0) then
              do i = 1, nbonds(jj)
                j = ibonds(i,jj)
                if (nat(j) == 7 .and. j /= icc) then
                  if (nbonds(j) == 2) then
                    k = ii   ! N-C-N structure found.
                    ii = jj  ! ii is the C in N-C-N
                    jj = k   ! C(jj) - N(icc) - C(ii) - N(kk)
                  end if
                  kk = j   ! kk is the distant N
                end if
              end do
            end if
            if (nat(ii) == 6) then
              i = nbonds(ii)
              do j = 1,i
                k = ibonds(j,ii)
                if (nat(k) == 8) then
                  Rab = distance(ii,k)
                  if (Rab < 1.33d0) kk = 0  ! >C=O detected, therefore not conjugated
                end if
              end do
            end if
            if (kk /= 0) then
              l = 0
              if (nbonds(kk) >= 2) then
                ii = ibonds(1,kk)
                if (nat(ii) == 1 .or. nat(ii) == 6) l = 1
                ii = ibonds(2,kk)
                if (nat(ii) == 1 .or. nat(ii) == 6) l = l + 1
                if (nbonds(kk) == 3) then
                  ii = ibonds(3,kk)
                  if (nat(ii) == 1 .or. nat(ii) == 6) l = l + 1
                end if
                if (l == 3) then                     ! distant nitrogen in a histidine (imidazole) ring
                  hybrid(icc) = 3                    ! has a hydrogen or a carbon attached
                  if (ionized) then
                    nH = 1
                    angle = 125.5d0
                    dihedral = 180.d0
                    return
                  else
                    nH = 0
                    return
                  end if
                else
                  nH = 1
                  angle = 125.5d0
                  dihedral = 180.d0
                  hybrid(icc) = 3
                  return
                end if
              else
                if (distance(icc, ii) > 1.4d0 .or. distance(icc, jj) > 1.4d0) then
!
!  Not aromatic/conjugated
!
                  nH = 1
                  angle = 125.5d0
                  dihedral = 180.d0
                  hybrid(icc) = 3
                  return
                else
                  nH = 0
                  return
                end if
              end if
            end if
            hybrid(icc) = 3
            if (nat(ii) == 6 .and. nat(jj) == 6) then
              if (aromatic(icc,ii,jj)) then
                nH = 0                               ! aromatic ring
                hybrid(icc) = 2
                return
              end if
              if (aromatic_5(icc,ii,jj, hybrid)) then
                nH = 1                               ! five-membered aromatic ring
                angle = 125.5d0
                dihedral = 180.d0
                return
              end if
            end if
            if (nat(ii) == 6 .and. nat(jj) == 8) then
               if (nbonds(jj) == 1) then
                 nH = 0                              ! Nitroso
               else
                nH = 1                               ! -O-NH-C
                angle = 125.5d0
                dihedral = 180.d0
              end if
              return
            end if
            if (nat(jj) == 15) then
              if (Rac < 1.7d0) then
                nH = 0
              else
                nH = 1
                angle = 120.d0
                dihedral = 180.d0
              end if
              return
            end if
            nH = 1
            angle = 120.d0
            dihedral = 180.d0
          else
            hybrid(icc) = 3
          end if
          Rab = Rab
                                        ! Oxygen is very important
        case (8)                        ! so there are a lot of options for O-H bonds
                                        ! (Add more as necessary)
          if (nbonds(icc) == 1) then
            if (nat(ii) == 5) then
              nH = 1
              angle = 110.3d0
              return
            end if
            if (nat(ii) == 6) then
!
!  Check for single or double bond
!
              bond_length = 0.975d0                  ! Ethanol
              angle = 110.3d0                        ! Ethanol
              nH = 1                                 ! Default
              hybrid(icc) = 3
              if (nbonds(ii) == 4) then
                hybrid(icc) = 3
                return
              end if
!
!  See if oxygen is attached to an aromatic ring, if so, then it's a hydroxyl.
!
              if (nbonds(ii) == 3) then
                do i = 1,3
                  jj = ibonds(i,ii)
                  if (jj /= icc) exit
                end do
                do i = i + 1,3
                  kk = ibonds(i,ii)
                  if (kk /= icc .and. kk /= jj) exit
                end do

                if (aromatic(ii, jj, kk)) then
                  if (distance(ii,icc) > 1.3d0) then
                    nH = 1                             ! aromatic ring, therefore hydroxyl
                    hybrid(icc) = 3
                  else
                    nH = 0                           ! Looks like hydroxyl, but C-O bond is too short
                    hybrid(icc) = 2
                  end if
                  return
                end if
              end if
!
!  Check for -COOH
!
              jj = nbonds(ii)
              if (jj == 3) then
                k = 0
                i = 0
                do l = 1, jj
                  kk = ibonds(l,ii)
                  if (nat(kk) == 8) then             !  Found an oxygen atom
                    if (nbonds(kk) == 2) i = i + 1
                    k = k + 1
                    if (kk /= icc) j = kk
                  end if
                end do
                if (k == 2) then                     !  Two oxygen atoms attached to a carbon
                  if (txtatm(ii)(15:16) == " ") then
!
!  This carbon is unlabeled, therefore it is a terminal carbon
!
                    l_tmp = (icc > j)
                  else
                    call find_polar_atom(icc, ii, j, Rab)
                    Rac = distance(ii, j)              !  Pick the longer of the two C-O distances
                    if (Rab < 0.51d0) Rab = distance(ii,icc) ! No polar atom near to icc or j
                    l_tmp =  (Rab > Rac .and. i == 0)
                  end if
                  if (l_tmp) then   !  and neither of them has another atom attached
                    if (ionized) then
                      nH = 0
                    else
                      bond_length = 1.00d0
                      nH = 1
                      angle = 110.d0
                      dihedral = 0.d0
                      do l = 1, jj
                        kk = ibonds(l,ii)
                        if (nat(kk) == 8 .and. kk /= icc) nc_icc = kk
                      end do
                    end if
                    return
                  else
                    nH = 0
                    hybrid(icc) = 2
                    return
                  end if
                end if
              end if
              if (Rab > 1.37d0) then
                nH = 1                               ! Almost certainly hydroxyl
                bond_length = 0.97d0
                angle = 110.d0
                hybrid(icc) = 3
                return
              end if
              if (Rab < 1.3d0 .and. hybrid(ii) == 2) then
                nH = 0                               ! Almost certainly carbonyl
                hybrid(icc) = 2
                return
              end if
!
! start of low-reliability tests
!
              if (Rab < 1.23d0) then
                  nH = 0
                  return
                end if
                if (Rab > 1.17d0) then
                nH = 0
                j = 0
                do i = 1, nbonds(ii)
                  if (nat(ibonds(i,ii)) /= 1) j = j + 1
                end do
                if (nat(ii) == 6 .and. j == 2) then
                  jj = ibonds(1,ii)
                  if (jj == icc) jj = ibonds(2,ii)
                  call bangle (coord, icc, ii, jj, sum)
                  if (sum < 2.03d0) then
                    nH = 1
                    bond_length = 0.97d0
                    angle = 110.d0
                    hybrid(icc) = 3
                    return
                  end if
                end if
              end if
              if (Rab > 1.27d0) then
                if (nat(ii) == 6 .and. nbonds(ii) == 3) then  ! check for planarity
                  jj = ibonds(1,ii)
                  kk = ibonds(2,ii)
                  ll = ibonds(3,ii)
                  i = 0
                  if (nat(jj) == 6) i = 1
                  if (nat(kk) == 6) i = i + 1
                  if (nat(ll) == 6) i = i + 1
                  if (i == 2) then
                    call dihed(coord, ii, jj, kk, ll, sum)
                    if (min(2*pi - sum, sum) < 0.25d0) then
                      nH = 1                         !  If flat, make sp3, and don't check anything else
                      hybrid(icc) = 3
                      return
                    end if
                  end if                             !  Not sure
                end if
              else
                nH = 0                               ! >C=O
                hybrid(icc) = 2
              end if
            else if (nat(ii) == 7) then
!
! O-N  is it HNO3?
!
              if (.not. ionized) then
                if (nbonds(ii) == 3) then
                  j = 0
                  do i = 1, 3
                    if (nat(ibonds(i,ii)) == 8) j = j + 1
                  end do
                  if (j == 3) then  ! Yes, it's a NO3
  !
  ! First, break O-O bonds
  !
                    do i = 1, 3
                      k = ibonds(i,ii)
                      if (nbonds(k) /= 2) cycle
                      if (ibonds(1,k) /= ii) then
                        j = ibonds(1,k)
                        ibonds(1,k) = ibonds(2,k)
                        ibonds(2,k) = j
                      end if
                    end do
                    do ll = 1, 2
                      i = ibonds(ll,ii)
                      if (nbonds(i) == 2) then
                        jj = ibonds(2,i)
                        do k = ll + 1,3
                          j = ibonds(k,ii)
                          if (nbonds(j) == 2) then
                            kk = ibonds(2,j)
                            if (i == kk) then
                              nbonds(i) = 1
                              nbonds(j) = 1
                            end if
                          end if
                        end do
                      end if
                    end do
  !
  ! now check for oxygen bonded to other atoms
  !
                    j = 0
                    do i = 1,3
                      if (nbonds(ibonds(i,ii)) > 1) j = 1
                    end do
                    if (j == 0) then
                      nH = 1  ! NO3, therefore add a hydrogen
                      angle = 120.d0
                      dihedral = 180.d0
                      hybrid(icc) = 3
                      do i = 1,3
                        if (ibonds(i,ii) /= icc) exit
                      end do
                      jj = ibonds(i,ii)
                      do i = i + 1,3
                        if (ibonds(i,ii) /= icc) exit
                      end do
                      jj = ibonds(i,ii)
                    end if
                  else
                    if (Rab > 1.3d0) then
                      nH = 1  !  Long bond, therefore hydroxyl
                      angle = 120.d0
                      dihedral = 180.d0
                      hybrid(icc) = 3
                    else
                      nH = 0
                    end if
                  end if
                end if
              end if
            else if (nat(ii) == 8) then
              nH = 1
              angle = 104.d0
              dihedral = 90.d0
            else if (nat(ii) == 12 .or. nat(ii) == 20) then
              nH = 1                                 ! Mg or Ca
              angle = 150.d0
              hybrid(icc) = 3
            else if (nat(ii) == 15 .or. nat(ii) == 33) then
              hybrid(icc) = 3
               if (nbonds(ii) == 3) then
                 nH = 1
                 angle = 120.d0
                 dihedral = 180.d0
                 do i = 1, 3
                 nd_icc = ibonds(i,ii)
                   if (nd_icc /= nb_icc .and. nd_icc /= icc .and. nd_icc /= nc_icc) exit
                 end do
                 hybrid(icc) = 3
!
!  Test for phosphate and arsenate
!
              else if (nbonds(ii) == 4) then
                jj = 0
                j = 0
                do i = 1, nbonds(ii)
                  kk = ibonds(i,ii)
                  if (nbonds(kk) == 1 .and. nat(kk) == 8) jj = jj + 1
                  if (nat(kk) == 8) j = j + 1
                end do
!
!  jj is the number of oxygen atoms that are not bonded to another atom
!
                if (ionized .and. j == 4) then
                  nH =0
                  return
                end if
                if (jj > 1) then
                  nH = 1
                  angle = 150.d0
                  dihedral = 0.d0
                end if
              end if
            else if (nat(ii) == 16 .or. nat(ii) == 34) then
!
!  Test for sulfate, sulfite, selenate, and selenite
!
              hybrid(icc) = 3
              if (nbonds(ii) == 2) then
                if (nat(ibonds(1,ii)) == 6 .or. nat(ibonds(2,ii)) == 6) then
                  nH = 1
                  angle = 120.d0
                  dihedral = 90.d0
                end if
              else if (nbonds(ii) > 2) then
                jj = 0
                j = 0
                do i = 1, nbonds(ii)
                  kk = ibonds(i,ii)
                  if (nbonds(kk) == 1 .and. nat(kk) == 8) jj = jj + 1
                  if (nat(kk) == 8) j = j + 1
                end do
                if (ionized .and. j == 4) then
                  nH =0
                  return
                end if
!
!  jj is the number of oxygen atoms that are not bonded to another atom
!
                if (jj > nbonds(ii) - 2) then
                  nH = 1
                  angle = 150.d0
                  dihedral = 0.d0
                end if
              else
                hybrid(icc) = 3
              end if
            else
              hybrid(icc) = 3
            end if
          else
            hybrid(icc) = 2
          end if
!
!
!
           case (14, 32, 50)
             bond_length = 1.3d0
             select case (nbonds(icc))
             case (1)
               nH = 3
               angle = tetrahedral
               internal_dihedral = 120.d0
               dihedral = 180.d0
               hybrid(icc) = 3
             case (2)
               nH = 2
               angle = tetrahedral
               internal_dihedral = 180.d0
               dihedral = 90.d0
               hybrid(icc) = 3
             case (3)
               nH = 1
               angle = 125.5d0
               dihedral = 180.d0
               hybrid(icc) = 3
             case (4, 5, 6)
               nH = 0
             end select
             return
           case (15)
           hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            nH = 2                                   ! -PH2
            bond_length = 1.42d0
            angle = 93.5d0
            internal_dihedral = 93.8d0
          else if (nbonds(icc) == 2) then
            nH = 1                                   ! >PH
            bond_length = 1.42d0
            angle = 122.d0
            dihedral = 180.d0
          end if
        case (16)
          hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            if (Rab < 2.00d0 .and. nat(ii) == 16 .or. &
                Rab < 1.69d0 .and. nat(ii) ==  6) then
              nH = 0
            else
              if (Rab < 1.65d0) then
                nH = 0
              else
                nH = 1                                   ! -S-H
                bond_length = 1.3d0
                angle = 100.d0
                dihedral = 90.d0
              end if
            end if
          end if
        case (33)
          hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            nH = 2                                   ! -AsH2
            bond_length = 1.52d0
            angle = 92.d0
            internal_dihedral = 92.d0
          else if (nbonds(icc) == 2) then
            nH = 1                                   ! >AsH
            bond_length = 1.52d0
            angle = 95.d0
            internal_dihedral = 95.d0
          end if
        case (34)
          hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            nH = 1                                   ! -Se-H
            bond_length = 1.46d0
            angle = 91.d0
          end if
        case (51)
          hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            nH = 2                                   ! -SbH2
            bond_length = 1.71d0
            angle = 92.d0
            internal_dihedral = 92.d0
          else if (nbonds(icc) == 2) then
            nH = 1                                   ! >SbH
            bond_length = 1.71d0
            angle = 95.d0
            internal_dihedral = 95.d0
            hybrid(icc) = 3
          end if
        case (52)
          hybrid(icc) = 3
          if (nbonds(icc) == 1) then
            nH = 1                                   ! -Te-H
            bond_length = 1.69d0
            angle = 90.d0
          end if
        case default
          hybrid(icc) = 3
        end select
    end if
  end subroutine h_type

  subroutine add_a_generic_hydrogen_atom(na, nb, nc, bond_length, angle, dihedral, metals, nmetals)
!
!  Add a single hydrogen atom to a heavy atom.
!  (This subroutine is generic, and is derived from subroutine geout)
!
    use common_arrays_C, only : coord, nbonds, ibonds, nat, txtatm
    use molkst_C, only : numat, id, temp_1, temp_2, temp_3, line
    use chanel_C, only: iw
    implicit none
    double precision, intent (in) :: bond_length, angle, dihedral
    integer, intent (in) :: na, nb, nc, nmetals, metals(nmetals)
!
    double precision :: xb, yb, zb, rbc, xa, ya, za, xpa, xpb, costh, sinth, &
      sinph, cosph, zqa, yza, xyb, ypa, coskh, sinkh, sina, cosd, &
      xd, yd, zd, xpd, ypd, zqd, xqd, yqd, xrd, cosa, sind, zpd
    integer :: k
    double precision, external :: distance
    logical, external :: near_a_metal
      cosa = cos(angle)
      if (id > 0) then
        xb = distance(nb, na)
        xb = temp_1
        yb = temp_2
        zb = temp_3
      else
        xb = coord(1,nb) - coord(1,na)
        yb = coord(2,nb) - coord(2,na)
        zb = coord(3,nb) - coord(3,na)
      end if
      rbc = xb*xb + yb*yb + zb*zb
      if (rbc < 1.d-4) then
        write(line,'(a, i6, a, i6, a)')"Atoms",na, " and", nb, " have the same coordinates."
        call mopend(trim(line))
        write(iw,'(/10x, a, i6, a, a, 3f10.3)') &
          "Atom No.:", na, "  Label: """//trim(txtatm(na))//"""", " Coordinates:", coord(:,na)
        write(iw,'(10x, a, i6, a, a, 3f10.3)') &
          "Atom No.:", nb, "  Label: """//trim(txtatm(nb))//"""", " Coordinates:", coord(:,nb)
        write(iw,'(/10x,a)')"Correct error before continuing."
        return
      end if
      rbc = 1.0D00/sqrt(rbc)
      if (nc /= 0) then
        if (id > 0) then
          xa = distance(nc, na)
          xa = temp_1
          ya = temp_2
          za = temp_3
        else
          xa = coord(1,nc) - coord(1,na)
          ya = coord(2,nc) - coord(2,na)
          za = coord(3,nc) - coord(3,na)
        end if
      else
        xa = 1.d0
        ya = 2.d0
        za = 3.d0
      end if
!
!     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!
      xyb = sqrt(xb*xb + yb*yb)
      k = -1
      if (xyb <= 0.009d0) then
        xpa = za
        za = -xa
        xa = xpa
        xpb = zb
        zb = -xb
        xb = xpb
        xyb = sqrt(xb*xb + yb*yb)
        k = 1
      end if
!
!     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!
      costh = xb/xyb
      sinth = yb/xyb
      xpa = xa*costh + ya*sinth
      ypa = ya*costh - xa*sinth
      sinph = zb*rbc
      cosph = sqrt(abs(1.D00 - sinph*sinph))
      zqa = za*cosph - xpa*sinph
!
!     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!
      yza = sqrt(ypa**2 + zqa**2)
      if (yza >= 1.D-4) then
        coskh = ypa/yza
        sinkh = zqa/yza
      else
!
!   ANGLE TOO SMALL TO BE IMPORTANT
!
        coskh = 1.D0
        sinkh = 0.D0
      end if
!
!     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!     NONE ARE NEGATIVE.
!     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
!
      sina = sin(angle)
      sind = -sin(dihedral)
      cosd = cos(dihedral)
      xd = bond_length*cosa
      yd = bond_length*sina*cosd
      zd = bond_length*sina*sind
!
!     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!
      ypd = yd*coskh - zd*sinkh
      zpd = zd*coskh + yd*sinkh
      xpd = xd*cosph - zpd*sinph
      zqd = zpd*cosph + xd*sinph
      xqd = xpd*costh - ypd*sinth
      yqd = ypd*costh + xpd*sinth
      if (k >= 1) then
        xrd = -zqd
        zqd = xqd
        xqd = xrd
      end if
      numat = numat + 1
      coord(1,numat) = xqd + coord(1,na)
      coord(2,numat) = yqd + coord(2,na)
      coord(3,numat) = zqd + coord(3,na)
      nbonds(na) = nbonds(na) + 1
      nbonds(numat) = 1
      nat(numat) = 1
      ibonds(1,numat) = na
      ibonds(nbonds(na),na) = numat
      if (near_a_metal(na, numat, metals, nmetals)) numat = numat - 1
      return
  end subroutine add_a_generic_hydrogen_atom
!
!
  subroutine add_a_sp3_hydrogen_atom(icc, nb_icc, nc_icc, nd_icc, bond_length, metals, nmetals)
    use molkst_C,  only : numat
    use common_arrays_C,  only : coord, nbonds, ibonds, nat
    implicit none
    integer,  intent (in) :: icc, nmetals, metals(nmetals)
    integer ::  nb_icc, nc_icc, nd_icc
    double precision :: bond_length
    double precision :: x, y, z, ax, ay, az, r1, r2, sum, angle, dihedral
    logical, external :: near_a_metal
!
!  Adding a hydrogen atom to a heavy atom.  The heavy atom is atom number "icc"
!  In this subroutine, a hydrogen atom is added to a single heavy atom.
!  Four atoms are needed: atom icc, which will have the hydrogen atom attached,
!  and atoms nb_icc, nc_icc, and nd_icc which are attached to icc.
!  These atoms form a triangle, this triangle is used in determining the direction that the hydrogen
!  atom is pointed.
!
    ax = coord(1, icc)
    ay = coord(2, icc)
    az = coord(3, icc)
    x = (coord(1,nb_icc) + coord(1,nc_icc) + coord(1,nd_icc))/3.d0
    y = (coord(2,nb_icc) + coord(2,nc_icc) + coord(2,nd_icc))/3.d0
    z = (coord(3,nb_icc) + coord(3,nc_icc) + coord(3,nd_icc))/3.d0
    r1 = sqrt((x-ax)**2+(y-ay)**2+(z-az)**2)
    if (r1 < 0.4d0 .and. nat(icc) == 6 .or. r1 < 0.3d0) then
!
!  System is almost flat, so put hydrogen perpendicular to the plane
!
      angle = 1.570796d0
      dihedral = 1.57d0
      call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals)
      r1 =         (coord(1,nb_icc) - coord(1,numat))**2 + &
        (coord(2,nb_icc) - coord(2,numat))**2 + (coord(3,nb_icc) - coord(3,numat))**2
      r1 = min(r1, (coord(1,nc_icc) - coord(1,numat))**2 + &
        (coord(2,nc_icc) - coord(2,numat))**2 + (coord(3,nc_icc) - coord(3,numat))**2)
      r1 = min(r1, (coord(1,nd_icc) - coord(1,numat))**2 + &
        (coord(2,nd_icc) - coord(2,numat))**2 + (coord(3,nd_icc) - coord(3,numat))**2)
      nbonds(icc) = nbonds(icc) - 1
      dihedral = -dihedral
      call add_a_generic_hydrogen_atom(icc, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals)
      r2 =         (coord(1,nb_icc) - coord(1,numat))**2 + &
        (coord(2,nb_icc) - coord(2,numat))**2 + (coord(3,nb_icc) - coord(3,numat))**2
      r2 = min(r2, (coord(1,nc_icc) - coord(1,numat))**2 + &
        (coord(2,nc_icc) - coord(2,numat))**2 + (coord(3,nc_icc) - coord(3,numat))**2)
      r2 = min(r2, (coord(1,nd_icc) - coord(1,numat))**2 + &
        (coord(2,nd_icc) - coord(2,numat))**2 + (coord(3,nd_icc) - coord(3,numat))**2)
      if (r2 > r1) coord(:,numat - 1) = coord(:,numat)
      numat = numat - 1
      return
    end if
    sum=bond_length/r1
    numat = numat + 1
    coord(1, numat) = ax + sum*(ax-x)
    coord(2, numat) = ay + sum*(ay-y)
    coord(3, numat) = az + sum*(az-z)
    nbonds(icc) = nbonds(icc) + 1
    nbonds(numat) = 1
    nat(numat) = 1
    ibonds(1,numat) = icc
    ibonds(nbonds(icc),icc) = numat
    if (near_a_metal(icc, numat, metals, nmetals)) numat = numat - 1
  end subroutine add_a_sp3_hydrogen_atom
!
!
  logical function aromatic(atom_1, atom_2, atom_6)
!
!  Given three adjacent atoms, atom_2 - atom_1 - atom_6, determine if they belong
!  to an aromatic ring.  If so, "aromatic" return .true.
!
    use common_arrays_C,  only : coord,  nbonds,  ibonds, nat
    use funcon_C,  only : pi
    implicit none
    integer, intent (in) :: atom_1, atom_2, atom_6
    integer :: atom_3, atom_5, i, j, k, l, m, nii, njj, kka, atom_4, ii
    double precision :: bent = 0.4d0, sum, Rcc
    double precision, external :: distance
      aromatic = .false.
      nii = nbonds(atom_2)
      njj = nbonds(atom_6)
      do i = 1, nii
        atom_3 = ibonds(i,atom_2)
        if (atom_3 == atom_1) cycle
        do j = 1, njj
          atom_5 = ibonds(j,atom_6)
          if (atom_5 == atom_1) cycle
          do k = 1, nbonds(atom_5)
            kka = ibonds(k,atom_5)
            if (kka == atom_6) cycle
            do l = 1, nbonds(atom_3)
              atom_4 = ibonds(l,atom_3)
              if (atom_4 == atom_2) cycle
              if (atom_4 == atom_5) then
!
!  Five-membered ring, check for N, O or S
!
                if (nat(atom_4) == 7 .or. nat(atom_4) == 8 .or. nat(atom_4) == 16 .or. &
                    nat(atom_3) == 7 .or. nat(atom_3) == 8 .or. nat(atom_3) == 16) then
                  call dihed(coord, atom_1, atom_2, atom_4, atom_6, sum)
                  if (min(2*pi - sum, sum) < bent) then
                    call dihed(coord, atom_1, atom_2, atom_3, atom_4, sum)
                    if (nat(atom_1) /= 6) then
                      Rcc = 2.d0 ! Definitely not aromatic
                    else
                      Rcc = min(distance(atom_1, atom_2), distance(atom_1, atom_6))
                    end if
                    aromatic = (min(2*pi - sum, sum) < bent .and. Rcc < 1.43d0 )
                  end if
                  return
                end if
              end if
              if (kka == atom_4) then
                if (atom_1 == atom_3 .or. atom_1 == atom_4 .or. atom_1 == atom_5) exit
                if (atom_2 == atom_3 .or. atom_2 == atom_4 .or. atom_2 == atom_5) exit
                if (atom_3 == atom_4 .or. atom_3 == atom_5 .or. atom_3 == atom_6) exit
                if (atom_4 == atom_5 .or. atom_4 == atom_6 .or. atom_5 == atom_6) exit
                m = max(nat(atom_1), nat(atom_2), nat(atom_3), nat(atom_4), nat(atom_5), nat(atom_6))
                if (m == 8 .or. m == 14 .or. m == 16) exit
                m = 0
                if (nat(atom_1) == 7 .and. nbonds(atom_1) > 2) m = m + 1
                if (nat(atom_2) == 7 .and. nbonds(atom_2) > 2) m = m + 1
                if (nat(atom_3) == 7 .and. nbonds(atom_3) > 2) m = m + 1
                if (nat(atom_4) == 7 .and. nbonds(atom_4) > 2) m = m + 1
                if (nat(atom_5) == 7 .and. nbonds(atom_5) > 2) m = m + 1
                if (nat(atom_6) == 7 .and. nbonds(atom_6) > 2) m = m + 1
!
!  Check for six-membered rings that have one or more >C=O (exclude >C-O)
!
                if (nat(atom_1) == 7) then
                  if (nat(atom_2) == 6 .and. nbonds(atom_2) == 3) then
                    do ii = 1, 3
                      if (nat(ibonds(ii, atom_2)) == 8) exit
                    end do
                    if (ii < 4) then
                      if (distance(atom_2, ibonds(ii, atom_2)) < 1.29d0) m = m + 2
                    end if
                  end if
                  if (nat(atom_6) == 6 .and. nbonds(atom_6) == 3) then
                    do ii = 1, 3
                      if (nat(ibonds(ii, atom_6)) == 8) exit
                    end do
                    if (ii < 4) then
                      if (distance(atom_6, ibonds(ii, atom_6)) < 1.29d0) m = m + 2
                    end if
                  end if
                end if
                if (m > 1) exit
!
!  Six-membered ring, with atoms numbered atom_1, atom_2, atom_3, atom_4 = kka, atom_5, atom_6
!
                call dihed(coord, atom_1, atom_2, atom_3, atom_4, sum)
                if (min(2*pi - sum, sum) < bent) then
                  call dihed(coord, atom_2, atom_3, atom_4, atom_5, sum)
                  if (min(2*pi - sum, sum) < bent) then
                    call dihed(coord, atom_3, atom_4, atom_5, atom_6, sum)
                    if (min(2*pi - sum, sum) < bent) then
                      call dihed(coord, atom_4, atom_5, atom_6, atom_1, sum)
                      if (min(2*pi - sum, sum) < bent) then
                        call dihed(coord, atom_5, atom_6, atom_1, atom_2, sum)
                        if (min(2*pi - sum, sum) < bent) then
                          call dihed(coord, atom_6, atom_1, atom_2, atom_3, sum)
                          if (min(2*pi - sum, sum) < bent) then
!
! Flat six-membered ring
!
                          aromatic = .true.
                          return
                          end if
                        end if
                      end if
                    end if
                  end if
                end if
              end if
            end do
          end do
        end do
      end do
      return
  end function aromatic
  logical function aromatic_5(atom_1, atom_2, atom_5, hybrid)
!
!  Given three adjacent atoms, atom_2 - atom_1 - atom_5, determine if they belong
!  to an five-membered aromatic ring.  If so, "aromatic_5" return .true.
!
    use common_arrays_C,  only : coord,  nbonds,  ibonds
    use molkst_C, only : numat
    use funcon_C,  only : pi
    implicit none
    integer, intent (in) :: atom_1, atom_2, atom_5, hybrid(numat)
    integer :: atom_3, i, j, k, nii, njj, kka, atom_4, l
    double precision :: bent = 0.4d0, sum
      if (numat < 5) then
        aromatic_5 = .false.
        return
      end if
      nii = nbonds(atom_2)
      njj = nbonds(atom_5)
      do i = 1, nii
        atom_3 = ibonds(i,atom_2)
        if (atom_3 == atom_1) cycle
        do j = 1, njj
          atom_4 = ibonds(j,atom_5)
          if (atom_4 == atom_1) cycle
          do k = 1, nbonds(atom_4)
            kka = ibonds(k,atom_4)
            if (kka == atom_3) then
              l = max(hybrid(2), hybrid(3), hybrid(4), hybrid(5))
              if (l == 3) cycle  ! if any atom in the ring is sp3, the ring cannot be aromatic
!
!  five-membered ring, with atoms numbered atom_1, atom_2, atom_3 = kka, atom_4, atom_5
!
              call dihed(coord, atom_1, atom_2, atom_3, atom_4, sum)
              if (min(2*pi - sum, sum) < bent) then
                call dihed(coord, atom_2, atom_3, atom_4, atom_5, sum)
                if (min(2*pi - sum, sum) < bent) then
                  call dihed(coord, atom_3, atom_4, atom_5, atom_1, sum)
                  if (min(2*pi - sum, sum) < bent) then
                    call dihed(coord, atom_4, atom_5, atom_1, atom_2, sum)
                    if (min(2*pi - sum, sum) < bent) then
                      call dihed(coord, atom_5, atom_1, atom_2, atom_3, sum)
                      if (min(2*pi - sum, sum) < bent) then
!
! Flat five-membered ring
!
                        aromatic_5 = .true.
                        return
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end do
        end do
      end do
      aromatic_5 = .false.
      return
  end function aromatic_5
  logical function near_a_metal(icc, ii, metals, nmetals)
!
!  Work out whether the atom (icc) or the hydrogen atom (numat) is near
!  to a metal atom (metals(natoms)).  If so, set near_a_metal = .false.
!  and delete bonds formed.
!
  use common_arrays_C,  only : nbonds,  nat
  implicit none
  integer, intent (in) :: icc, ii, nmetals, metals(nmetals)
  integer :: i
  logical :: Sulfur, Oxygen, Nitrogen
  double precision :: Rab, Rac
  double precision, external :: distance
    do i = 1, nmetals
      Rab = distance(metals(i), icc)
      Rac = distance(metals(i), ii)
!
!  Check for Sulfur and Oxygen.  Maybe check for (not carbon)?
!
      Sulfur =   (nat(icc) == 16 .and. Rab < 2.6d0)
      Oxygen =   (nat(icc) == 8 .and. Rab < 2.1d0 .and. Rac <= Rab)
      Nitrogen = (nat(icc) == 7 .and. Rab < 2.5d0 .and. Rac <= Rab)
      if (Rac < 1.5d0) exit
      if ((Sulfur .or. Oxygen .or. Nitrogen)) exit
    end do
    near_a_metal = (i <= nmetals)
    if (.not. near_a_metal) return
    nbonds(icc) = nbonds(icc) - 1
    return
  end function near_a_metal



  subroutine bridge_H
!
!  Detect and print unusually short hydrogen bonds
!
  use molkst_C,  only : numat, line
  use common_arrays_C,  only : nat, nbonds,  ibonds, txtatm, &
    breaks
  use chanel_C, only : iw, log, ilog
  use elemts_C, only: elemnt
  implicit none
  integer :: i,ii, j, k, l, kk, jj, nchain2
  double precision :: bond_length, bond_length_2, sum
  integer, allocatable :: h1(:), h2(:)
  double precision, allocatable :: Rab1(:), Rab2(:)
  double precision, external :: distance
    allocate (h1(numat), h2(numat),Rab1(numat),Rab2(numat))
    Rab1 = 100.d0
    Rab2 = 100.d0
    k = 0
    ii = 0
    do i = 1, numat
      if (nat(i) /= 1) cycle
      bond_length = 10.d0
      bond_length_2 = 10.d0
      nchain2 = 1
      l = 0
      k = 0
      do j = 1, numat
        if (j == breaks(nchain2)) nchain2 = nchain2 + 1
        if (nat(j) == 1) cycle
        sum = distance(i,j)
!
!  Exclude water
!
        jj = 0
        if (nat(j) == 8) then
          jj = 0
          do kk = 1, nbonds(j)
            if (nat(ibonds(kk,j)) == 1) jj = jj + 1
          end do
        end if
        if (jj < 2) then
          if (sum < bond_length_2) then
            if (sum < bond_length) then
              bond_length_2 = bond_length
              l = k
              bond_length = sum
              k = j
            else
              bond_length_2 = sum
              l = j
            end if
          end if
        end if
      end do
      if (k > 0) then
!
! Rab1 = bond length from hydrogen to the nearer heavy atom
! Rab2 = bond length from hydrogen to the more distant heavy atom
! h1   = atom that hydrogen is nearer to.
! h2   = atom that hydrogen is more distant from
!
        Rab1(i) = bond_length
        h1(i) = k
        Rab2(i) = bond_length_2
        h2(i) = l
      end if
    end do
    do j = 1, 20
      bond_length_2 = 100.d0
      do i = 1, numat
        if (Rab2(i) < bond_length_2) then
          bond_length_2 = Rab2(i)
          k = i
        end if
      end do
      if (bond_length_2 > 1.5d0) exit
      if (ii == 0) then
        ii = 1
        write(line,"(16x,a)")"  UNUSUALLY SHORT HYDROGEN BOND DISTANCES"
        write(iw,"(/,a)")trim(line)
        if(log) write(ilog,"(/,a)")trim(line)
        write(iw,"(a)")
        if(log) write(ilog,"(a)")
        write(line,"(a)")"  Hydrogen   R(H-A)  Atom A    R(H-B)  Atom B     R(A-B)      Loc-A      Loc-B"
        write(iw,"(a)")trim(line)
        if(log) write(ilog,"(a)")trim(line)
      end if
      write(line,"(i8,f10.2,2x,a,i5,f9.2,2x,a,i5, f10.2,5x,a,2x,a)") k, Rab2(k), elemnt(nat(h2(k))), &
          h2(k), Rab1(k), elemnt(nat(h1(k))), h1(k), distance(h1(k), h2(k)), &
          txtatm(h2(k))(18:20)//txtatm(h2(k))(23:27)//" "//txtatm(h2(k))(22:22), &
          txtatm(h1(k))(18:20)//txtatm(h1(k))(23:27)//" "//txtatm(h1(k))(22:22)
      write(iw,"(a)")trim(line)
      if(log) write(ilog,"(a)")trim(line)
      Rab2(k) = 10.d0
    end do
    if (ii == 1) write(iw,*)
    if (ii == 1 .and. log) write(ilog,*)
    return
  end subroutine bridge_H
  subroutine reset_breaks
    use molkst_C,  only : numat, nbreaks, line, keywrd
    use common_arrays_C,  only : coord, txtatm, break_coords, breaks
    use chanel_C, only : iw
    implicit none
    integer :: i, j, k, l
    double precision :: sum
    logical :: first = .true.
    double precision, external :: reada
    save :: first
!
! Re-set "breaks" - they move when RESEQ or ADD_H is called
!
      do nbreaks = 1, 50
        numat_loop: do j = 1, numat
          sum = abs(coord(1,j) - break_coords(1,nbreaks)) + &
                abs(coord(2,j) - break_coords(2,nbreaks)) + &
                abs(coord(3,j) - break_coords(3,nbreaks))
          if (sum < 0.1d0) then
            exit numat_loop
          end if
        end do numat_loop
        if (j > numat) exit
  !
  !  Now move to the end of the residue
  !
        line = txtatm(j)(23:)
        k = nint(reada(line,1))
        do
          j = j + 1
          if (j > numat) exit
          line = txtatm(j)(23:)
          l = nint(reada(line,1))
          if (l /= k) exit
        end do
        breaks(nbreaks) = j - 1
      end do
!
!  Sequence breaks into increasing order
!
      j = 1
      nbreaks = nbreaks - 1
      do i = 1, nbreaks
        l = breaks(i)
        do k = i + 1, nbreaks
          if (breaks(k) < l) then
            breaks(i) = breaks(k)
            breaks(k) = l
            l = breaks(i)
          end if
        end do
        if (i > 1) then
          if (breaks(i) /= breaks(i - 1)) then
            j = j + 1
            breaks(j) = breaks(i)
          end if
        end if
      end do
      nbreaks = j
      breaks(nbreaks + 1) = 0
!
!  The following block is probably never used - delete it sometime (note date: 12/6/2015)
!
      if (index(keywrd, " NORES") == 0) then
        do i = 1, nbreaks
          do j = i + 1, nbreaks
            if (breaks(j) < breaks(i)) then
              if (first .and. index(keywrd, "GEO-OK") == 0) then
                first = .false.
!
!   An entire block has moved, so the breaks are no longer valid
!
!   Print a warning and quit
!
                call mopend("WARNING: After RESEQ, PDB ""TER"" locations are incorrect")
                write(iw,'(5x,a)')" To correct this error, edit the data set or PDB file"// &
                " or add keyword NORESEQ"
                write(iw,'(5x,a)')" To suppress this error, add ""GEO-OK"" to the keywords"
              end if
            end if
          end do
        end do
      end if
      j = 1
      do i = 2, nbreaks
        if (breaks(i) - breaks(i - 1) > 0) then
          j = j + 1
          breaks(j) = breaks(i)
        end if
      end do
      nbreaks = j
      breaks(nbreaks + 1) = 0
      j = 1
      do i = 1, numat
        write(line,'(a6,i5,a)')txtatm(i)(:6),i + j - 1,txtatm(i)(12:)
        txtatm(i) = trim(line)
        if (i == breaks(j)) j = j + 1
      end do
  end subroutine reset_breaks
  subroutine find_polar_atom(O1, C, O2, Rab)
!
! Given a -COO group, -(C)(O1)(O2), find the nearest polar oxygen.
! If it's O1, set Rab large
! If it's O2, set Rab small
! If there isn't one, set Rab very small
!
  use molkst_C,  only : numat
  use common_arrays_C,  only : nat,  coord,  nbonds,  ibonds
  implicit none
  integer, intent (in) :: O1, C, O2
  double precision, intent (inout) :: Rab
  integer :: i, j, k, n1, n2, n1_near(10), n2_near(10)
  double precision :: RO1, RO2, angle, RO1_near(10), RO2_near(10)
  double precision, external :: distance
    n1 = 0
    n2 = 0
    RO1 = 0.d0
    do i = 1, numat
      if (nat(i) == 8) then
        if (i /= O1) then
          do j = 1, nbonds(i)
            k = ibonds(j,i)
            if (nat(k) == 6) then
              RO1 = distance(O1, i)
              if (RO1 < 4.d0) then
                call bangle(coord, C, O1, i, angle)
                if (angle > 1.570796d0) then ! 90 degrees
                  call dihed (coord, i, O1, C, O2, angle)
                  if (angle < 1.570796d0 .or. angle > 4.71239d0) then
                    n1 = n1 + 1
                    n1_near(n1) = i
                    RO1_near(n1) = RO1
                  end if
                end if
              end if
            end if
          end do
        end if
        if (i /= O2) then
          do j = 1, nbonds(i)
            k = ibonds(j,i)
            if (nat(k) == 6) then
              RO2 = distance(O2, i)
              if (RO2 < 4.d0) then
                call bangle(coord, C, O2, i, angle)
                if (angle > 1.570796d0) then ! 90 degrees
                  call dihed (coord, i, O2, C, O1, angle)
                  if (angle < 1.570796d0 .or. angle > 4.71239d0) then
                    n2 = n2 + 1
                    n2_near(n2) = i
                    RO2_near(n2) = RO1
                  end if
                end if
              end if
            end if
          end do
        end if
      end if
    end do
    if (n1 == 0) then
      if (n2 == 0) return
!
!  There is an oxygen near to O2.  Set Rab small, so that the hydrogen is put on O2
!
      Rab = 0.52d0
    else
      if (n2 == 0) then
!
!  There is an oxygen near to O1.  Set Rab large, so that the hydrogen is put on O1
!
        Rab = 4.02d0
      else
!
!  There are no oxygen atoms near tp O1 or O2.  Set Rab large, so that the C-O1 and C-O2
!  distances will be used in deciding which oxygen atom has the hydrogen
!
        RO1 = 10.d0
        do i = 1, n1
          if (RO1 > RO1_near(i)) then
            RO1 = RO1_near(i)
            j = n1_near(i)
          end if
        end do
        RO2 = 10.d0
        do i = 1, n2
          if (RO2 > RO2_near(i)) then
            RO2 = RO2_near(i)
            k = n2_near(i)
          end if
        end do
        if (RO1 < RO2) then
!
!  There is an oxygen near to both O1 and O2.  Use the O1-O and O2-O distances
!  in deciding where to put the hydrogen atom.
!
           Rab = 4.1d0
        else
           Rab = 0.5d0
        end if
      end if
    end if
    return
  end subroutine find_polar_atom
  logical function guanidine (icc, ionized, nH, angle, internal_dihedral, dihedral, hybrid)
!
!   Test nitrogen atom "icc" to check that it's a "nu-1" or "nu-2" of a guanidine group in an Arg
!
  use common_arrays_C,  only : nat,  nbonds,  ibonds, txtatm
  use molkst_C,  only : numat
  implicit none
  integer :: nH, icc, hybrid(numat)
  double precision :: angle, internal_dihedral, dihedral
  logical :: ionized
!
!  Local
!
  integer :: i, j, k, ii, jj, loop
  guanidine = .false.
  if (.not. (txtatm(icc)(13:16) == " NH1" .or. txtatm(icc)(13:16) == " NH2")) return
  do loop = 1, nbonds(icc)
    ii = ibonds(loop,icc)
    do i = 1, nbonds(ii)
      jj = ibonds(i,ii)
      if (nat(jj) /= 6 .and. nat(jj) /= 7) exit
      if (nat(jj) == 7 .and. jj /= icc) then  !  N-C-N structure
        k = 0
        do j = 1, nbonds(jj)
          if (nat(ibonds(j,jj)) == 6) k = k + 1
        end do
        if (k == 2) cycle               !  Exclude C-N-C
        if (nbonds(jj) == 3) then
          hybrid(icc) = 2
          guanidine = .true.
          if (ionized) then
            nH = 2                      !   ionized and an Arg
            angle = 120.d0
            internal_dihedral = 180.d0
            dihedral = 0.d0
          else
            nH = 1
            angle = 120.d0
            internal_dihedral = 180.d0
            dihedral = 0.d0
          end if
          return                         !  Unconditional
        else
          nH = 2
          angle = 120.d0
          internal_dihedral = 180.d0
          dihedral = 0.d0
          hybrid(icc) = 3
          return
        end if
      end if
    end do
  end do
  return
  end function guanidine
