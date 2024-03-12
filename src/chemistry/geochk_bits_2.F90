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

  subroutine site(neutral, chain, res, alt, charge, nres, max_sites, allkey)
  use common_arrays_C, only : geo, coord, labels, nat, nbonds, ibonds, &
    atmass, txtatm, tvec, lopt
  use parameters_C, only : ams
  use chanel_C, only : iw
  use elemts_C, only : elemnt
  use molkst_C, only : numat, natoms, line, keywrd, id, moperr,maxtxt
  use funcon_C,  only : pi
  implicit none
  logical, intent(in) :: neutral(100)
  integer, intent (in) :: nres, res(nres), max_sites
  character, intent (in) :: chain(nres)*1, alt(nres)*1
  character, intent (inout) :: charge(max_sites,3)*1, allkey*3000
  integer :: i_atom, i_charge
  integer :: i, j, k, l, m, n, o(2), jres, metals(1), nmetals = 0, nadd, ndel, p, ncarbon, &
    nb_icc, nc_icc, nd_icc
  logical :: l_res, bug = .false., first, let, l_type, first_2
  character :: res_txt*16, txt*26, num*1, line_1*1200
  character, allocatable :: changes(:)*30
  double precision :: bond_length, angle, dihedral, internal_dihedral, tetrahedral, sum
  logical, allocatable :: l_used(:)
  integer, allocatable :: hybrid(:)
  double precision, external :: reada, distance
  i = max(20,numat)
  allocate (changes(i), l_used(natoms), hybrid(natoms))
  l_used = .false.
  labels(numat + 1:) = 0
  tetrahedral = 109.4712206d0
  let = (index(keywrd," LET ") /= 0)
  first = (.not. let)
  first_2 = .true.
  i = index(keywrd," SITE=")
  j = index(keywrd(i:), ") ") + i
  k = index(allkey, "SALT")
  if (k /= 0) then
    k = index(allkey, " SITE")
    l = index(allkey(k:), ") ") + k - 1
    allkey = " """//allkey(k + 1:l)//""" => "//keywrd(i:j)
  end if
  changes = " "
  nadd = 0
  ndel = 0
  ncarbon = 0
  lopt(:,natoms + 1:) = 1
  geo(:,:natoms) = coord(:,:natoms)
  i_charge = 100
  do i_atom = 1, numat
    if (txtatm(i_atom)(22:22) /= " " .and. txtatm(i_atom)(22:22) /= "A") exit
  end do
  if (i_atom > numat) then
!
! No chain letters, so insert letter "A" into all atom labels
!
    do i_atom = 1, numat
      txtatm(i_atom)(22:22) = "A"
    end do
  end if
  call check_cvs(.true.)
  jres = 0
  do i_atom = 1, numat
    l_res = .false.
    if (nres > 0) then
      do jres = 1, nres
        if (txtatm(i_atom)(22:22) /= chain(jres)) cycle
        if (txtatm(i_atom)(27:27) /= alt(jres)) then
          i = i + 1
          cycle
        end if
        i = nint(reada(txtatm(i_atom), 23))
        if (i == res(jres) .and. charge(jres,1) /= "*") exit
      end do
      l_res = (jres <= nres)
      if (l_res) then
        if (charge(jres,1) == "+") then
          i_charge = 1
        else if (charge(jres,1) == "-") then
          i_charge = -1
        else
          i_charge = 0
        end if
      else
        i_charge = 100 ! this means "not l_res"
      end if
    end if
!
!  Deliberately add or remove hydrogen atoms as required by keyword "SITE"
!  This request uses sites defined in the "C" style in the manual, e.g., "[HIS]239:A.NE2"(+)
!
    if (neutral(5) .or. neutral(6) .or. l_res) then !  Must be -Arg(+)- or -Arg(0)-
!
!  Check for Arg(0) or Arg(+)
!
      l_type = (txtatm(i_atom)(18:20) == "ARG")
      if (l_type .and. nat(i_atom) == 6 .and. nbonds(i_atom) == 3) then
        if (nat(ibonds(1, i_atom)) == 7 .and. &
            nat(ibonds(2, i_atom)) == 7 .and. &
            nat(ibonds(3, i_atom)) == 7 ) then! System is C bonded to three N, check for number of H on N
          i = 0
          do j = 1, 3
            i = i + nbonds(ibonds(j,i_atom))
          end do
          if (i_charge /= 100 .and. (i == 8 .or. i == 9)) charge(jres,1:2) = charge(jres,2:3)
          if (i /= 9 .and. (neutral(6) .or. i_charge == 0) .or. i /= 8 .and. (neutral(5) .or. i_charge == 1)) cycle
          if (neutral(6) .or. i_charge == 0) then
!
! Found the system -NH-C(NH2)2, now delete a hydrogen atom
!
            do j = 1, 3
              if (nbonds(ibonds(j, i_atom)) == 3) then
                l = ibonds(j, i_atom)
                k = 0
                if (nat(ibonds(1,l)) == 1) k = k + 1
                if (nat(ibonds(2,l)) == 1) k = k + 1
                if (nat(ibonds(3,l)) == 1) k = k + 1
                if (k == 2) then
                  if (nat(ibonds(1,l)) == 1) then
                    k = ibonds(1,l)
                  else
                    k = ibonds(2,l)
                  end if
                  nat(k) = 99
                  ibonds(nbonds(i_atom), i_atom) = 0
                  nbonds(i_atom) = nbonds(i_atom) - 1
                  exit
                end if
              end if
            end do
          else if (neutral(5) .or. i_charge == 1) then
!
! Found the system -NH-C(NH2)(NH), now add a hydrogen atom
!
            do j = 1, 3
              if (nbonds(ibonds(j, i_atom)) == 2) then ! Found a N with two bonds
                l = ibonds(j, i_atom)
                call add_sp2_H(ibonds(1,l), l, ibonds(2,l))
                numat = numat + 1
                if (i_charge /= 100) charge(jres,1:2) = charge(jres,2:3)
                nbonds(i_atom) = nbonds(i_atom) + 1
                ibonds(nbonds(i_atom), i_atom) = numat
                nadd = nadd + 1
                changes(nadd)(:15) = trim(txtatm(l)(12:))
                write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(l)(:6), natoms, txtatm(l)(15:)
                exit
              end if
            end do
          end if
        end if
      end if
    end if
    if (neutral(11) .or. neutral(12) .or. l_res) then !  Must be -Lys(+)- or -Lys(0)-
!
!  Check for Lys(0) or Lys(+)
!
      l_type = (txtatm(i_atom)(18:20) == "LYS")
      if (l_type .and. txtatm(i_atom)(14:15) == "NZ" ) then
        if (neutral(12) .or. i_charge == 0) then
!
! Found the system -NH3(+), now delete a hydrogen atom
!
          do j = 1, 3
            if (nbonds(ibonds(j, i_atom)) == 3) then
              l = ibonds(j, i_atom)
              k = 0
              if (nat(ibonds(1,l)) == 1) k = k + 1
              if (nat(ibonds(2,l)) == 1) k = k + 1
              if (nat(ibonds(3,l)) == 1) k = k + 1
              if (k == 2) then
                if (nat(ibonds(1,l)) == 1) then
                  k = ibonds(1,l)
                else
                  k = ibonds(2,l)
                end if
                nat(k) = 99
                ibonds(nbonds(i_atom), i_atom) = 0
                nbonds(i_atom) = nbonds(i_atom) - 1
                exit
              end if
            end if
          end do
        else if ((neutral(11) .or. i_charge == 1) .and. nbonds(i_atom) == 3) then
!
!  Nitrogen bonded to three atoms, add a hydrogen atom.
!
          if (nbonds(i_atom) == 3 .and. i_charge /= 100) charge(jres,1:2) = charge(jres,2:3)
          call add_sp3_H(ibonds(1,i_atom), i_atom, ibonds(2,i_atom), ibonds(3,i_atom))
          nbonds(i_atom) = nbonds(i_atom) + 1
          ibonds(nbonds(i_atom), i_atom) = numat
          numat = numat + 1
          nadd = nadd + 1
          changes(nadd)(:15) = trim(txtatm(i_atom)(12:))
          write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(i_atom)(:6), natoms, txtatm(i_atom)(15:)
        end if
      end if
    end if
    if (neutral(1)  .or. neutral(2) .or. l_res) then
      if (nat(i_atom) == 6 .and. (nbonds(i_atom) == 3 .or. (nbonds(i_atom) > 2 .and. l_res))) then
!
! check for -COOH or -COO(-) or -CH2OH or -CH2O(-)
!
        k = 0
        l = 0
        m = 0
        do i = 1,3
          j = ibonds(i,i_atom)
          if (nat(j) == 8) then
            if (nbonds(j) == 1) then
              k = k + 1
              m = m + 1
            else if (nbonds(j) == 2 .and. (nat(ibonds(1,j)) == 1 .or. nat(ibonds(2,j)) == 1)) then
              k = k + 1
            end if
            l = j
          end if
        end do
!
!  k is the number of oxygen atoms that are ionizable or are ionized
!
        if (k > 0) then
          if (k == 2 .or. k == 1 .and. (nbonds(i_atom) == 4 .or. &! A bit ugly, but it checks for ">C=O" of a peptide bond
            (l_res .and. nbonds(i_atom) == 3 .and. distance(i_atom, l) > 1.3d0))) then
            if (i_charge == 1) then
!
!  Panic: -COO found, but charge is impossible.  Check for another charge on the same residue
!
              do i = jres + 1, nres
                if (chain(i) == chain(jres) .and. res(i) == res(jres)) then
!
!  Found one.  Swap charges and hope for the best
!
                    do k = 1, 3
                      num = charge(i,k)
                      charge(i,k) = charge(jres,k)
                      charge(jres,k) = num
                    end do
                    if (charge(jres,1) == "+") then
                      i_charge = 1
                    else if (charge(jres,1) == "-") then
                      i_charge = -1
                    else
                      i_charge = 0
                    end if
                    exit
                  end if
              end do
            end if
            if (m == 2 .and. i_charge == -1) then
              if (.not. bug) write(iw,*)
              write(iw,'(10x,a)')"A charge of -1 requested for -COO in """//txtatm(i_atom)(18:)// &
              & """, but current charge is already -1"
              bug = .true.
            end if
            if(neutral(2) .or. i_charge == -1) then ! Make it -COO(-) or -CH2O(-)
              do i = 1, nbonds(i_atom)
                j = ibonds(i,i_atom)
                o = 0
                o(1) = nat(ibonds(1,j))
                if (nbonds(j) == 2) o(2) = nat(ibonds(2,j))
                if (nat(j) == 8 .and. nbonds(j) == 2 .and. (o(1) == 1 .or. o(2) == 1)) then
                  if(nat(ibonds(1,j)) == 1) nat(ibonds(1,j)) = 99
                  if(nat(ibonds(2,j)) == 1) nat(ibonds(2,j)) = 99
                  ibonds(nbonds(j), j) = 0
                  nbonds(j) = nbonds(j) - 1
                  if (jres > 0) charge(jres,1:2) = charge(jres,2:3)
                  exit
                end if
              end do
            else if(neutral(1) .or. i_charge == 0) then
              k = 0
              l = 0
              do i = 1, nbonds(i_atom)
                j = ibonds(i,i_atom)
                if (nat(j) == 8 .and. nbonds(j) == 2) then
                  if (nat(ibonds(1,j)) == 1 .or. nat(ibonds(2,j)) == 1) k = k + 1
                end if
                if (nat(j) == 8 .and. nbonds(j) == 1) then
                  l = l + 1
                  o(l) = j
                end if
              end do
              if (k == 0 .and. (l == 2 .and. nbonds(i_atom) == 3 .or. l == 1 .and. nbonds(i_atom) == 4)) then
!
!  -COO(-) or -CH2O(-) group identified.  Now add a hydrogen
!
                call add_sp_H(j, o(1), o(2))
                m = 1
                call h_type(22, j, .false., m, nb_icc, nc_icc, nd_icc, bond_length, angle, dihedral, internal_dihedral, &
        hybrid, metals, nmetals)
                if (abs(bond_length - 1.d0) < 1.d-3) then
                  angle = angle*pi/180.d0
                  call add_a_generic_hydrogen_atom(j, nb_icc, nc_icc, bond_length, angle, dihedral, metals, nmetals)
                  geo(:,numat) = coord(:,numat)
                  nbonds(j) = nbonds(j) + 1
                  ibonds(nbonds(j), j) = numat
                  nadd = nadd + 1
                  changes(nadd)(:15) = trim(txtatm(j)(12:))
                  charge(jres,1:2) = charge(jres,2:3)
                  write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(j)(:6), natoms, txtatm(j)(15:)
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    if (neutral(3) .or. neutral(4) .or. l_res) then !  Must be -NH3(+) or -NH2
      if (nat(i_atom) == 7) then
!
!  Exclude -NH-C-(NH2) type
!
        do j = 1, nbonds(i_atom)
          if (nat(ibonds(j,i_atom)) == 6) then
            i = ibonds(j,i_atom)
            l = 0
            do k = 1, nbonds(i)
              if (nat(ibonds(k,i)) == 7) l = l + 1
            end do
            if (l == 3) exit
          end if
        end do
        if (l == 3) cycle
        i = 0
        l = 0
        do j = 1, nbonds(i_atom)
          if (nat(ibonds(j,i_atom)) == 1) i = i + 1
          do k = 1, nbonds(ibonds(j,i_atom))
            m = ibonds(k, ibonds(j,i_atom))
            if (nat(m) == 7 .and. m /= i_atom) l = 1
          end do
        end do
        if (l == 1) cycle
        if (i == 2 .and. (neutral(3) .or. i_charge == 1) .and. nbonds(i_atom) == 3) then
!
!  Selectively exclude -CO-NH2 if all -NH2 are flagged for ionization
!
          if (neutral(3)) then
            do i = 1, nbonds(i_atom)
              j = ibonds(i,i_atom)
              if (nat(j) == 6) exit
            end do
            if (i <= nbonds(i_atom)) then
              do i = 1, nbonds(j)
                if (nat(ibonds(i,j)) == 8) exit
              end do
              if (i <= nbonds(j)) cycle
            end if
          end if
!
!  Nitrogen bonded to two hydrogen atoms, add a hydrogen atom.
!
          if ((i == 2 .and. nbonds(i_atom) == 3 .or. i == 3 .and. nbonds(i_atom) == 4) &
          .and. i_charge /= 100) charge(jres,1:2) = charge(jres,2:3)
          call add_sp3_H(ibonds(1,i_atom), i_atom, ibonds(2,i_atom), ibonds(3,i_atom))
          nbonds(i_atom) = nbonds(i_atom) + 1
          ibonds(nbonds(i_atom), i_atom) = numat
          numat = numat + 1
          nadd = nadd + 1
          changes(nadd)(:15) = trim(txtatm(i_atom)(12:))
          write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(i_atom)(:6), natoms, txtatm(i_atom)(15:)
          continue
        else if (i == 3 .and. (neutral(4) .or. i_charge == 0) .and. nbonds(i_atom) == 4) then
!
!  Nitrogen bonded to three hydrogen atoms, delete a hydrogen atom.
!
          do i = 1, 4
            if (nat(ibonds(i, i_atom)) == 1) then
              j = ibonds(i, i_atom)
              nat(j) = 99
              ibonds(nbonds(i_atom), i_atom) = 0
              nbonds(i_atom) = nbonds(i_atom) - 1
              charge(jres,1:2) = charge(jres,2:3)
              exit
            end if
          end do
        end if
      end if
    end if
    if (neutral(7) .or. neutral(8) .or. l_res) then!  Must be -His(+) or -His(0)-
!
!  Check for His(0) or His(+)
!
      if (nat(i_atom) == 6 .and. nbonds(i_atom) == 3) then
        j = 0
        k = 0
        do i = 1, 3
          if (nat(ibonds(i, i_atom)) == 7 ) j = j + 1
          if (nat(ibonds(i, i_atom)) == 1 ) k = k + 1
        end do
        if (j == 2 .and. k == 1) then! System is C bonded to two N and one H, therefore His
          i = 0
          do j = 1,3
            k = ibonds(j,i_atom)
            do l = 1, nbonds(k)
              m = ibonds(l,k)
              if (nat(m) == 1) i = i + 1
            end do
          end do
          if (i /= 2 .and. (neutral(8) .or. i_charge == 0) .or. i /= 1 .and. (neutral(7) .or. i_charge == 1)) cycle
          if (i_charge /= 100) charge(jres,1:2) = charge(jres,2:3)
          if (neutral(8) .or. i_charge == 0) then
!
! Found a His(+), now delete a hydrogen atom
!
            outer_loop: do j = 1,3
              k = ibonds(j,i_atom)
              do l = 1, nbonds(k)
                m = ibonds(l,k)
                if (nat(m) == 1) then
                  nat(m) = 99
                  ibonds(nbonds(i_atom), i_atom) = 0
                  nbonds(i_atom) = nbonds(i_atom) - 1
                  exit outer_loop
                end if
              end do
            end do outer_loop
          else if (neutral(7) .or. i_charge == 1) then
!
! Found the system -C-N-C-, now add a hydrogen atom
!
            do j = 1, 3
              if (nbonds(ibonds(j, i_atom)) == 2) then ! Found a N with two bonds
                l = ibonds(j, i_atom)
                call add_sp2_H(ibonds(1,l), l, ibonds(2,l))
                numat = numat + 1
                nbonds(l) = nbonds(l) + 1
                ibonds(nbonds(l), l) = numat
                nadd = nadd + 1
                changes(nadd)(:15) = trim(txtatm(l)(12:))
                write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(l)(:6), natoms, txtatm(l)(15:)
                exit
              end if
            end do
          end if
        end if
      end if
    end if
    if (neutral(9)) then  !  Check for sulfate.  If the structure S-O-H exists, delete the "H"
      if (nat(i_atom) == 8 .and. nbonds(i_atom) == 2) then
        do j = 1, 2
          if (nat(ibonds(j,i_atom)) == 16) exit  !  Look for a sulfur atom attached to the oxygen
        end do
        if (j < 3) then
          do k = 1, 2
            if (nat(ibonds(k,i_atom)) == 1) exit  !  Look for a hydrogen atom attached to the oxygen
          end do
          if (k < 3) then
            l = ibonds(j,i_atom)
            m = 0
            do n = 1, nbonds(l)
              if (nat(ibonds(n,l)) == 8) m = m + 1
            end do
            if (m == 4) then
              m = ibonds(k,i_atom)
              nat(m) = 99
              ibonds(k, i_atom) = 0
              nbonds(i_atom) = nbonds(i_atom) - 1
            end if
          end if
        end if
      end if
    end if
    if (neutral(10)) then
!  Check for phosphate.
!
!  Aim for [PO4](-), e.g., CH3-PO4-H
!
      if (nat(i_atom) == 8 .and. nbonds(i_atom) == 2) then
        do j = 1, 2
          if (nat(ibonds(j,i_atom)) == 15) exit  !  Look for a phosphorus atom attached to the oxygen
        end do
        if (j < 3) then
          if(nbonds(ibonds(j,i_atom)) == 4) then
!
!   oxygen atoms on phosphorus should be attached to only two ligands
!
            m = ibonds(j,i_atom)
            l = 0
            do k = 1, 4
              if (nat(ibonds(k,m)) == 8 .and. nbonds(ibonds(k,m)) == 2) then
                if (nat(ibonds(1, ibonds(k,m))) /= 99 .and. nat(ibonds(2, ibonds(k,m))) /= 99) l = l + 1
              end if
            end do
            if (l > 2) then
              do k = 1, 2
                if (nat(ibonds(k,i_atom)) == 1) exit  !  Look for a hydrogen atom attached to the oxygen
              end do
              if (k < 3) then
                l = ibonds(j,i_atom)
                m = 0
                do n = 1, nbonds(l)
                  if (nat(ibonds(n,l)) == 8) m = m + 1
                end do
                if (m == 4) then
                  nat(ibonds(k,i_atom)) = 99
                  ibonds(nbonds(i_atom), i_atom) = 0
                  nbonds(i_atom) = nbonds(i_atom) - 1
                end if
              end if
            end if
          end if
        end if
      end if
    end if
    if (l_res) then
      if (charge(jres,1) /= "*") then
!
!  Check for sulfate
!
        if (nat(i_atom) == 16 .and. nbonds(i_atom) == 4) then
          k = 2
          do i = 1, 4
            j = ibonds(i, i_atom)
            if (nat(j) /= 8) exit
            if (nbonds(j) == 1) k = k - 1
          end do
          if (i == 5) then

!
!  Found a sulfate.
!
            i_charge = min(0, i_charge)
            if (charge(jres,2) == "-") i_charge = -2
!
!  Work out the change in the number of ionized atoms
!  If i_charge becomes negative, delete hydrogen atoms
!
            i_charge = i_charge - k
            if (i_charge < 0) then
!
!   Now delete hydrogens, as necessary.
!
              do i = 1, -i_charge
                do l = 1, 4
                  j = ibonds(l, i_atom)
                  if (nbonds(j) == 2) then
                    if(nat(ibonds(2,j)) == 1) exit
                  end if
                end do
                if (l < 5) then
                  nat(ibonds(2,j)) = 99
                  ibonds(nbonds(i_atom), i_atom) = 0
                  nbonds(i_atom) = nbonds(i_atom) - 1
                end if
              end do
            else if (i_charge > 0) then
!
!   Now add hydrogens, as necessary.
!
              l = 1
              do i = 1, 4
                j = ibonds(i, i_atom)
                if (nbonds(j) == 1) then
                  do k = 1, 4
                    if (k /= i) exit
                  end do
                  call add_sp_H(j, i_atom, ibonds(k,i_atom))
                  numat = numat + 1
                  nbonds(i_atom) = nbonds(i_atom) + 1
                  ibonds(nbonds(i_atom), i_atom) = numat
                  nadd = nadd + 1
                  changes(nadd)(:15) = trim(txtatm(j)(12:))
                  l = l + 1
                  if (l > i_charge) exit
                end if
              end do
            end if
            charge(jres,1) = "*"
          end if
        end if
      end if
!
!  Check for phosphate
!
      if (nat(i_atom) == 15 .and. nbonds(i_atom) == 4) then
        k = 1
        do i = 1, 4
          j = ibonds(i, i_atom)
          if (nat(j) /= 8) exit
          if (nbonds(j) == 1) k = k - 1
        end do
        if (i == 5) then

!
!  Found a phosphate.
!
          i_charge = min(0, i_charge)
          if (charge(jres,2) == "-") i_charge = -2
!
!  Work out the change in the number of ionized atoms
!  If i_charge becomes negative, delete hydrogen atoms
!
          i_charge = i_charge - k
          if (i_charge < 0) then
!
!   Now delete hydrogens, as necessary.
!
            do i = 1, -i_charge
              do l = 1, 4
                j = ibonds(l, i_atom)
                if (nbonds(j) == 2) then
                  if(nat(ibonds(2,j)) == 1) exit
                end if
              end do
              if (l < 5) then
                nat(ibonds(2,j)) = 99
                charge(jres,1:2) = charge(jres,2:3)
                ibonds(nbonds(i_atom), i_atom) = 0
                nbonds(i_atom) = nbonds(i_atom) - 1
              end if
            end do
          else if (i_charge > 0) then
!
!   Now add hydrogens, as necessary.
!
            l = 1
            do i = 1, 4
              j = ibonds(i, i_atom)
              if (nbonds(j) == 1) then
                do k = 1, 4
                  if (k /= i) exit
                end do
                call add_sp_H(j, i_atom, ibonds(k,i_atom))
                numat = numat + 1
                nbonds(i_atom) = nbonds(i_atom) + 1
                ibonds(nbonds(i_atom), i_atom) = numat
                nadd = nadd + 1
                changes(nadd)(:15) = trim(txtatm(j)(12:))
                write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(j)(:6), natoms, txtatm(j)(15:)
                charge(jres,1:2) = charge(jres,2:3)
                l = l + 1
                if (l > i_charge) exit
              end if
            end do
          else
            charge(jres,1:2) = charge(jres,2:3)
          end if
        end if
      end if
    end if
  end do
!
!  Now check for explicit atoms
!
  i = index(keywrd," SITE=")
  j = index(keywrd(i:), ") ") + i
  if (j == i) then
    line = "ERROR IN KEYWORD ""SITE"" - no closing "")"""
    call mopend(trim(line))
    return
  end if
!
!  In the next block "txt_to_atom_no" cannot be used, because the PDB text is
!  needed for recognizing the atom.
!
  do
    k = Index (keywrd(i:), ") ") + i
    j = index(keywrd(i:k),"""") + i
    if (j == i) exit
    l = index(keywrd(j:),"""") + j
    line = keywrd(j:l - 2)
    if (line(1:1) == "[") then
!
! Atom defined using Jmol format
!
      n = index(keywrd(j:l - 2),".") + j
      k = index(keywrd(j:l - 2),":") + j
      m = index(keywrd(j:l - 2),"]") + j
      if (n == j .or. m == j) then
        if (n == j) write(iw,'(/10x,a)')"""."" Missing in JSmol atom label: """//keywrd(j:l - 2)//'"'
        if (m == j) write(iw,'(/10x,a)')"""]"" Missing in JSmol atom label: """//keywrd(j:l - 2)//'"'
        call mopend("ERROR IN JSMOL ATOM LABEL IN KEYWORD ""SITE""")
        return
      end if
      if (k == j) then
        line = keywrd(n:l - 2)//keywrd(j + 1:m - 2)//keywrd(j + 5:n - 2)
      else
        line = keywrd(n:l - 2)//keywrd(j + 1:m - 2)//keywrd(k:k)//keywrd(m:k - 2)
      end if
    end if
    keywrd(j:) = trim(line)//trim(keywrd(l - 1:))
    i = j + len_trim(line) + 2
  end do
  i = index(keywrd," SITE=")
  j = index(keywrd(i:), ") ") + i
  line = keywrd(i + 6:j)
  do
    i = index(line, """")
    if (i == 0) exit
    do j = i + 1, len_trim(line)
      if (line(j:j) == """") exit
    end do
!
! Convert one PDB label into an atom-number
! First, run a check to confirm that the PDB label exists in the data-set
! and that it is only used for one atom.
!
    line_1 = trim(line)
    call txt_to_atom_no(line_1, i, .false.)
    if (moperr) return
    res_txt = line(i + 1: j - 1)
    if (line(j + 2:j + 2) == "+") i_charge = 1
    if (line(j + 2:j + 2) == "0") i_charge = 0
    if (line(j + 2:j + 2) == "-") i_charge = -1
    if (line(j + 3:j + 3) == "2") i_charge = i_charge*2
    if (line(j + 3:j + 3) == "3") i_charge = i_charge*3
    line = line(j + 1:)
    m = 0
    do k = 1, len_trim(res_txt)
      if (res_txt(k:k) /= " ") then
        m = m + 1
        res_txt(m:m) = res_txt(k:k)
      end if
    end do
    res_txt(m + 1:) = " "
    do i = 1, numat
      txt = trim(txtatm(i))
      n = 0
      do k = 13, maxtxt
        if (txt(k:k) /= " ") then
          n = n + 1
          txt(n:n) = txt(k:k)
        end if
      end do
      txt(n + 1:) = " "
      n = max(n, m)
      do k = 1, n
        if (res_txt(k:k) /= txt(k:k) .and. res_txt(k:k) /= "*") exit
      end do
      if (k > n) exit
    end do
    if (i > numat) then
      if (first) then
        write(iw,'(/10x,a)')"An atom defined by keyword SITE does not match an atom in the system."
        write(iw,'(/10x,a)')"To ignore messages of this type, add keyword ""LET"""
        first = .false.
      end if
      k = index(allkey, " SITE")
      l = index(allkey(k:), ") ") + k - 1
      do
        k = k + 1
        if (k >= l) exit
        if (allkey(k:k) == '"') then
!
!  Start of an atom-label
!
          i = k
          do
            k = k + 1
            if (allkey(k:k) == '"') exit
          end do
!
!  End of atom-label
!
          if (allkey(i + 1:i + 1) == "[") then
            n = index(allkey(i:k),".") + i
            p = index(allkey(i:k),":") + i
            m = index(allkey(i:k),"]") + i
            if (n == i  .or. m == i) then
              if (n == i) write(iw,'(/10x,a)')"""."" Missing in JSmol atom label: "//allkey(i:k)
              if (m == i) write(iw,'(/10x,a)')"""]"" Missing in JSmol atom label: "//allkey(i:k)
              call mopend("ERROR IN JSMOL ATOM LABEL IN KEYWORD ""SITE""")
              return
            end if
            if (p == i) then
!
!  Chain-letter missing, so:
!
              txt = allkey(n:k - 1)//allkey(i + 2:i + 4)//allkey(i + 6:n - 2)
            else
              txt = allkey(n:k - 1)//allkey(i + 2:m - 2)//allkey(p:p)//allkey(m:p - 2)
            end if
          else
            txt = allkey(i + 1:k - 1)
          end if
          m = 0
          do n = 1, len_trim(txt)
            if (txt(n:n) /= " ") then
              m = m + 1
              txt(m:m) = txt(n:n)
            end if
          end do
          txt(m + 1:) = " "
          if (trim(txt) == trim(res_txt)) exit
        end if
      end do
      if (.not. let .and. first_2) call mopend("FAULT DETECTED IN KEYWORD ""SITE"" AT: "//allkey(i:k))
      first_2 = .false.
      cycle
    end if
    if (abs(i_charge) > 1 .and. nat(i) /= 6) then
      call mopend("IN ""SITE"" THE USE OF MULTIPLE CHARGES IS ONLY ALLOWED FOR CARBON ATOMS")
      return
    end if
    if (l_used(i)) then
      num = char(ichar("2") + int(log10(i*1.0001)))
      if (nat(i) == 6) then
        num = char(ichar("2") + int(log10(i*1.0001))) 
        write(iw,'(/5x, a, i'//num//', a)')"Two requests were made to modify the number of hydrogen atoms on atom", i, &
          ", a carbon atom, PDB label: """//txtatm(i)(:26)//"""."
        write(iw,'(5x, a)')"When more than one hydrogen atom is to be added or deleted use "// &
          """+3"" or ""+2"" or ""-2"" or ""-3"" instead of ""+"" or ""-""."
        call mopend("MULTIPLE REQUESTS WERE MADE USING ""SITE"" TO CHANGE THE NUMBER "// &
        "OF HYDROGEN ATOMS ON A CARBON ATOM.")
      else
        call mopend("MULTIPLE REQUESTS WERE MADE USING ""SITE"" TO CHANGE THE NUMBER "// &
        "OF HYDROGEN ATOMS ON A NON-CARBON ATOM.") 
        call mopend("THIS IS NOT ALLOWED.  USE SEPARATE JOBS.") 
      end if
      return
    else
      l_used(i) = .true.
    end if
    j = nat(i)
    select case (j)
    case (6) ! Add and delete hydrogen atoms on Carbon
      if (i_charge == 1) then
        nadd = nadd + 1
        changes(nadd)(:15) = trim(txtatm(i)(12:))
      end if
      ncarbon = ncarbon + 1
!
!  Delete existing hydrogen atoms
!
        m = 0
        do j = 1, nbonds(i)
          if (nat(ibonds(j,i)) == 1 .or. nat(ibonds(j,i)) == 99) then
             nat(ibonds(j,i)) = 99
             ibonds(j,i) = 0
             nbonds(i) = nbonds(i) - 1
             m = m + 1
          end if
        end do
!
! There were "m" hydrogen atoms, so new number of hydrogen atoms should be
!
        i_charge =  i_charge + m
        if (4 - nbonds(i) < i_charge) then
          num = char(ichar("2") + int(log10(i*1.0001)))
          write(iw,'(/5x, a, i'//num//', a)')"An attempt to add a hydrogen atom to atom", i, &
            ", a carbon atom, PDB label: """//txtatm(i)(:26)//""","
          write(iw,'(5x, a, i'//num//', a,/)')"failed because the atom was already saturated. "// &
            "Atoms currently bonded to atom", i, " are:"
          write(iw,'(10x, a, i'//num//')')" Atom No.    <------ PDB Label ------->      Distance to atom", i
          sum = 0.d0
          do j = 1, nbonds(i)
            bond_length = distance(i, ibonds(j,i))
            if (bond_length > sum) then
              sum = bond_length
              k = ibonds(j,i)
            end if
            write(iw,'(i5, i11, a, f12.4, a)')j, ibonds(j,i), "      """//txtatm(ibonds(j,i))(:26)//"""", bond_length, " Angstroms"
          end do
          write(iw,'(/2x,a)')"[ The simplest way to correct this fault would be to use the"// &
            & " CVB keyword, e.g., CVB=("""//txtatm(i)(13:26)//""":-"""//txtatm(k)(13:26)//""") ]"
          call mopend("ATTEMPT WAS MADE TO ADD A HYDROGEN ATOM TO A CARBON ATOM THAT WAS ALREADY SATURATED")
          return
        end if
        j = ibonds(1,i)
!
! Need to jump over any atoms that were hydrogen atoms.
! At this point, hydrogen atoms have ibonds(a,b) = 0, so they are easily identified.
!
        m = 0
        do l = 1, nbonds(i)
          m = m + 1
          do n = 1, 6
            k = ibonds(m,i)
            if (k /= 0) exit
            m = m + 1
          end do
          if (k /= i .and. k /= j) exit
        end do
        bond_length = 1.09d0
        select case (i_charge)
          case (1)
!
! Carbon with one hydrogen atom
!
            select case (nbonds(i))
              case (1)
!
! Carbon bonded to one other non-hydrogen atom
!
                angle = pi
                dihedral = pi
              case (2)
!
! Carbon bonded to two other non-hydrogen atoms
!
                angle = 120.d0*pi/180.d0
                dihedral = pi
              case (3)
!
! Carbon bonded to three other non-hydrogen atoms
!
                angle = 0.d0
            end select
            if (angle < 1.d0) then
              call add_a_sp3_hydrogen_atom(i, ibonds(1,i), ibonds(2,i), ibonds(3,i), bond_length, metals, nmetals)
            else
              call add_a_generic_hydrogen_atom(i, j, k, bond_length, angle, dihedral, metals, nmetals)
            end if
            geo(:,numat) = coord(:,numat)
            natoms = numat
          case (2)
!
! Carbon with two hydrogen atoms
!
            select case (nbonds(i))
              case (1)
!
! Carbon bonded to one other non-hydrogen atom
!
                angle = 120.d0*pi/180.d0
                internal_dihedral = pi
                dihedral = pi
              case (2)
!
! Carbon bonded to two other non-hydrogen atoms
!
                angle = tetrahedral*pi/180.d0
                internal_dihedral = 120.d0*pi/180.d0
                dihedral = 120.d0*pi/180.d0
            end select
            call add_a_generic_hydrogen_atom(i, j, k, bond_length, angle, dihedral, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
            internal_dihedral = internal_dihedral + dihedral
            call add_a_generic_hydrogen_atom(i, j, k, bond_length, angle, internal_dihedral, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
          case (3)
!
! Carbon with three hydrogen atoms,
! this can only bond to one other atom.
!
            j = ibonds(1,i)
            do l = 1, nbonds(j)
              k = ibonds(l,j)
              if (k > 0 .and. k /= i) exit
            end do
            angle = tetrahedral*pi/180.d0
            internal_dihedral = 120.d0*pi/180.d0
            dihedral = pi
            call add_a_generic_hydrogen_atom(i, j, k, bond_length, angle, pi, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
            call add_a_generic_hydrogen_atom(i, j, numat, bond_length, angle, internal_dihedral, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
            call add_a_generic_hydrogen_atom(i, j, numat, bond_length, angle, internal_dihedral, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
        end select
      case (7)! Select for Nitrogen in different coordination numbers
        select case (nbonds(i))
          case (1)
            j = ibonds(1,i)
!
!  Search for an atom attached to the atom that's attached to the nitrogen that's not that nitrogen
!
            do k = 1, nbonds(j)
              if (ibonds(k,j) /= i) exit
            end do
            k = ibonds(k,j)
            if (k == 0) then
              bond_length = 1.09d0
              angle = 2.6d0
              dihedral = pi
              call add_a_generic_hydrogen_atom(i, j, j, bond_length, angle, dihedral, metals, nmetals)
              geo(:,numat) = coord(:,numat)
            else
              call add_sp_H(i, j, k)
              numat = numat + 1
              nbonds(i) = nbonds(i) + 1
              ibonds(nbonds(i),i) = numat
            end if
            natoms = numat
            nadd = nadd + 1
            changes(nadd)(:15) = trim(txtatm(i)(12:))
            write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(i)(:6), natoms, txtatm(i)(15:)
          case (2)
          if(i_charge == 1) then
            call add_sp2_H(ibonds(1,i), i, ibonds(2,i))
            numat = numat + 1
            coord(:,numat) = geo(:,numat)
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            nadd = nadd + 1
            changes(nadd)(:15) = trim(txtatm(i)(12:))
            write(txtatm(natoms),'(a,i5," 3H",a)')txtatm(i)(:6), natoms, txtatm(i)(15:)
          else
!
!  Delete a hydrogen atom, and hope the user knows what's going on.
!
            do j = 1, nbonds(i)
              if (nat(ibonds(j,i)) == 1) exit
            end do
            if (j <= nbonds(i)) then
              nat(ibonds(j,i)) = 99
              ibonds(nbonds(i),i) = 0
              nbonds(i) = nbonds(i) - 1
            end if
          end if
        case (3)
          if (i_charge == 1) then
!
!  Three-coordinate nitrogen, therefore make into a quaternary nitrogen
!
            call add_a_sp3_hydrogen_atom(i, ibonds(1,i), ibonds(2,i), ibonds(3,i), 1.1d0, metals, nmetals)
            natoms = natoms + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i),i) = numat
            geo(:,numat) = coord(:,numat)
            nadd = nadd + 1
            changes(nadd)(:15) = trim(txtatm(i)(12:))
            write(txtatm(numat),'(a,i5," 3H",a)')txtatm(i)(:6), natoms, txtatm(i)(15:)
          else if (i_charge == 0) then
!
!  Three-coordinate nitrogen, therefore neutralize (assume pyridinium)
!
            k = 0
            do j = 1, nbonds(i)
              if (nat(ibonds(j,i)) == 1) then
                k = k + 1
                l = j
              end if
            end do
            if (k == 0) then
              write(iw,'(/1x,a,/)')" A nitrogen atom defined by keyword SITE is not bonded to any hydrogen atoms."
              write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), &
              "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
              call mopend(trim(line_1))
              bug = .true.
            else if (k == 2 .and. txtatm(i)(18:20) /= "ARG") then
              write(iw,'(3(/,a),/)')"  A nitrogen atom defined by keyword SITE is bonded to two hydrogen atoms, ", &
              "  but a charge of zero has been specified.  If a hydrogen atom were to be deleted", &
              "  the resulting geometry would not be chemically sensible."
              write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), &
              "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
              call mopend(trim(line_1))
              bug = .true.
            else
              nat(ibonds(l,i)) = 99
              ibonds(nbonds(i),i) = 0
              nbonds(i) = nbonds(i) - 1
            end if
          else
            write(iw,'(/1x,a)')"  A nitrogen atom defined by keyword SITE cannot have a negative charge."
            write(iw,'(a,/)')"  (If a hydrogen atom is to be deleted from a nitrogen, use ""(0)"" instead of ""(-)"")"
            write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), &
            "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
            call mopend(trim(line_1))
            bug = .true.
          end if
        case (4)
          if (i_charge == 0) then
!
!  Quaternary nitrogen, therefore neutralize
!
            do j = 1, 4
              if (nat(ibonds(j,i)) == 1) exit
            end do
            if (j > 4) then
              write(iw,'(/1x,a,/)')" A nitrogen atom defined by keyword SITE is not bonded to any hydrogen atoms."
              write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
              call mopend(trim(line_1))
              bug = .true.
            end if
            nat(ibonds(j,i)) = 99
            ibonds(nbonds(i),i) = 0
            nbonds(i) = nbonds(i) - 1
          else if (i_charge == 1) then
            write(iw,'(/1x,a,/)')" A nitrogen atom defined by keyword SITE is already bonded to four atoms."
            write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
            call mopend(trim(line_1))
          end if
        end select
      case (8, 16)! Select for Oxygen in different coordination numbers
  99    continue
        select case (nbonds(i))
        case (0)
          if (i_charge == 0 .or. i_charge == 1) then
            line = "An attempt has been made to change the number of hydrogen atoms attached to an isolated oxygen atom"
            call mopend(trim(line))
            write(iw, '(10x,a)')"(Either use CVB to make a bond to the oxygen, or make this change outside MOPAC)"
            return
          end if
        case (1)
          if (i_charge == 0 .or. i_charge == 1) then
            j = ibonds(1,i)
            do k = 1, nbonds(j)  !  First search for -COO structure
              if (ibonds(k,j) /= i .and. nat(ibonds(k,j)) == 8) exit
            end do
            if (k > nbonds(j)) then
              do k = 1, nbonds(j)   ! If no -COO, search for anything other that atom i
                if (ibonds(k,j) /= i) exit
              end do
            end if
            if (k > nbonds(j)) then
              do k = numat, 1, -1
                if (k /= j .and. k /= i) exit
              end do
              if (k == 0) then
                line = "A hydrogen atom cannot be added to a diatomic.  Do this outside MOPAC"
                call mopend(trim(line))
                return
              end if
            else
              k = ibonds(k,j)
            end if
            call add_sp_H(i, j, k)
            numat = numat + 1
            nadd = nadd + 1
            changes(nadd)(:15) = trim(txtatm(i)(12:))
            if (i_charge == 0) then
              nbonds(i) = 2
              ibonds(2,i) = natoms
            end if
          else if (i_charge == -1) then
            num = char(ichar("1") + int(log10(i*1.001)))
            write(iw,'(/,a,i'//num//',a,/,a,/)')"  Oxygen atom ",i, &
              " defined by keyword SITE is already bonded to only one atom,", &
              "  so a hydrogen atom cannot be deleted, as requested by the option ""(-)"""
            write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), &
              "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
            call mopend(trim(line_1))
            bug = .true.
          end if
        case (2)
          if (i_charge == -1) then
            do j = 1, nbonds(i)
              if (nat(ibonds(j,i)) == 1) exit
            end do
            if (j > nbonds(i)) then
              num = char(ichar("1") + int(log10(i*1.001)))
              write(iw,'(/1x,a,i'//num//',a,/)')"  Oxygen atom ",i, " defined by keyword SITE is not bonded to any hydrogen atoms."
              write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
              call mopend(trim(line_1))
              bug = .true.
            else
              nat(ibonds(j,i)) = 99
              ibonds(nbonds(i),i) = 0
              nbonds(i) = nbonds(i) - 1
            end if
          else if (i_charge == 0) then
            num = char(ichar("1") + int(log10(i*1.001)))
            write(iw,'(/1x,a,i'//num//',a,/)')"  Oxygen atom ",i, " defined by keyword SITE is already bonded to two atoms."
            write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
            call mopend(trim(line_1))
            bug = .true.
          else
            call add_sp2_H(ibonds(1,i), i, ibonds(2,i))
            numat = numat + 1
            nbonds(i) = nbonds(i) + 1
            ibonds(nbonds(i), i) = numat
            nadd = nadd + 1
            changes(nadd)(:15) = trim(txtatm(i)(12:))
          end if
        case (3)
          if (i_charge == -1 .or. i_charge == 0) then
            do j = 1, nbonds(i)
              if (nat(ibonds(j,i)) == 1) exit
            end do
            if (j > nbonds(i)) then
              write(iw,'(/1x,a,/)')" An oxygen atom defined by keyword SITE is not bonded to any hydrogen atoms."
              write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
              call mopend(trim(line_1))
              bug = .true.
            else
              nat(ibonds(j,i)) = 99
              nbonds(i) = 2
              if (j == 2) ibonds(2,i) = ibonds(3,i)
              if (i_charge == -1) goto 99
            end if
          else
            write(iw,'(/1x,a,/)')" An oxygen atom defined by keyword SITE is already bonded to three atoms."
            write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
            call mopend(trim(line_1))
            bug = .true.
          end if
        end select
      case default
        write(iw,'(/1x,a,/)')" An atom defined by keyword SITE is not a nitrogen, oxygen, or sulfur."
        write(line_1,'(a,4x,3f8.3,a)')'Faulty atom = "'//txtatm(i), coord(:,i), "  1.00  0.00      PROT"//elemnt(nat(i))//'"'
        call mopend(trim(line_1))
        bug = .true.
   end select
  end do
!
!  Check that all sites have been recognized, whether they are modified or not
!
  k = 0
  do j = 1 , nres
    if (charge(j,1) /= "*") then
      if (k == 0) then
        k = index(keywrd," SITE=")
        l = index(keywrd(k:), ") ") + k
        if (allkey /= " ") write(iw,'(/10x,a)')"Faulty SITE Keyword: '"//trim(allkey)//"'"
      end if
      k = index(keywrd," SITE=")
      if (j == 1) then
       do k = k, len_trim(keywrd)
         if (keywrd(k:k) == "(") exit
       end do
       k = k + 1
      else
        do i = 1, j - 1
          k = index(keywrd(k + 1:l), ",") + k + 1
        end do
      end if
      i = index(keywrd(k + 1:l), ",")
      if (i == 0) i = index(keywrd(k + 1:l + 1), ") ")
      i = i + k - 1
      write(line_1,'(a,i3,a)')"Site:", j, ", '"//keywrd(k:i)//"', was not properly recognized."
      call mopend(trim(line_1))
      write(iw,'(/10x,a)')"(Check that the residue is not duplicated and can be ionized or de-ionized)"
      if (index(keywrd, " 0SCF") /= 0) moperr = .false.
    end if
  end do
  if (bug) then
    if (allkey /= " ") then
      write(line_1,'(a)')"FAULTY SITE KEYWORD: '"//trim(allkey)//"'"
    else
      write(line_1,'(a)')"FAULTY SITE KEYWORD"
    end if
    call mopend(trim(line_1))
    write(iw,'(/10x,a,/10x,a)')"(Check that the atom is N or O and can be ionized or de-ionized)"
    if (index(keywrd, " 0SCF") /= 0) return
  end if
  do i = 1, natoms - id
    if (nat(i) == 99 .and. nbonds(i) > 0) then
      do j = 1, nadd
        if (changes(j)(:15) == trim(txtatm(ibonds(1,i))(12:))) exit
      end do
      if (j > nadd) then
        ndel = ndel + 1
        changes(ndel)(16:) = trim(txtatm(ibonds(1,i))(12:))
      end if
    end if
  end do
  j = 0
  do i = 1, ndel
    do k = 1, i - 1
      if (changes(k)(16:) == changes(i)(16:)) exit
    end do
    if ( k == i) then
      j = j + 1
      changes(j)(16:) = changes(i)(16:)
    end if
  end do
  ndel = j
  j = 0
  do i = 1, numat
    if (nat(i) /= 99) then
      j = j + 1
      geo(:,j) = geo(:,i)
      coord(:,j) = geo(:,i)
      nat(j) = nat(i)
      atmass(j) = ams(nat(i))
      txtatm(j) = txtatm(i)
      lopt(:,j) = lopt(:,i)
    end if
  end do
  numat = j
  natoms = numat
  labels(:numat) = nat(:numat)
  call set_up_dentate()
  call update_txtatm(.true., .false.)
  j = max(nadd, ndel)
  if (j > 0) then
    write(iw,'(/7x,a,/)')"Changes in Ionization caused by keyword SITE"
    n = 85
    k = len_trim(allkey)
    i = min(k,n - 8)
    l = index(allkey(i:), ",")
    if (l == 0) l = index(allkey(i:), "))")
    i = min(k, l + i) - 1
    if (i == k - 1) then
      write(iw,'(7x,a)')"Keyword:"//allkey(1:k)
    else
      write(iw,'(7x,a)')"Keyword:"//allkey(1:i)
      do
        l = i + 1
        i = min(l + n, k)
        i = min(k, index(allkey(i:), ",") + i)
        if (i == k) then
          if (allkey(i - 1:i - 1) == ",") allkey(i - 1:i) = ") "
          write(iw,'(7x, a)')allkey(l:i)
          exit
        else
          i = i - 1
          write(iw,'(7x, a)')allkey(l:i)
        end if
      end do
    end if
    write(iw,'(/5x,a,/)')" Hydrogen atoms added to     Hydrogen atoms deleted from"
    k = 0
    do i = 1, j
      if (index(changes(i)(17:18), "O") /= 0 .or. index(changes(i)(17:18), "N") /= 0) k = k - 1
      if (index(changes(i)( 2: 3), "O") /= 0 .or. index(changes(i)( 2: 3), "N") /= 0) k = k + 1
      if (changes(i)(:15) == " ") changes(i)(:15) = "       -"
      if (changes(i)(16:) == " ") changes(i)(16:) = "       -"
      write(iw,'(i3,5x,a)')i, changes(i)(:15)//"               "//trim(changes(i)(16:))
    end do
    if (mod(ncarbon,2) == 0) write(iw,'(/10x,a,SP,i4)')"Change in net ionization:", k
  end if
  do i = numat + 1, numat + id
    labels(i) = 107
    nat(i) = 107
    geo(:,i)  = tvec(:, i - numat)
  end do
  natoms = natoms + id
  return
  end subroutine site
