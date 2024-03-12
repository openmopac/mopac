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

subroutine PDB3
!
! Work out the atom label of a hydrogen atom, as defined by the PDB format definition 3.0 and higher
!
USE molkst_C, ONLY: numat, line
use MOZYME_C, only : tyres
use common_arrays_C, only : nat, txtatm, coord, nbonds, ibonds
use funcon_C, only : pi
implicit none
double precision :: torsion
integer :: i, ii, j, k, l, m, n
logical :: swap
logical, allocatable :: used(:)
!
! see: http://publications.iupac.org/pac/pdf/1998/pdf/7001x0117.pdf
!
! Greek letters are used as atom identifiers. The CA or the substituent closer to CA (in
! the order CA,CB,CG,. . . ) takes precedence over atoms in branches in defining stereochemical relationships. For example, if
! tetrahedral carbon C has four substituents X, Y, Z, and Z' (with priority X > Y > Z = Z'; i.e., Z and Z are diastereotopic sub-
! stituents designated provisionally as unprimed and primed), their numbering is derived as follows: if one sights down the
! X-C axis (with the X atom toward the viewer), the equivalent atoms, Z and Z', are designated Z2 and Z3, such that Y, Z2,and
! Z3 follow a clockwise orientation. The side-chain -NH2 nitrogens of Arg are designated as NHl and NH2 by their relationship
! (cis or trans, respectively) to CD. The hydrogen atoms of the side-chain -NH2 groups of Asn, Gln, and Arg are distinguished
! by numbers (1 or 2) on the basis of their relationship (cis or trans, respectively) to the heavy atom three bonds closer to the
! main chain (CB for Asn, CG for Gln, NE for Arg). Thus, each -NH2 hydrogen of Arg is distinguished by two numbers, the first
! indicating the nitrogen to which it is attached and the second indicating the stereochemistry of the hydrogen itself. Numbering
! of Phe and Tyr rings gives higher priority to the atom with the smaller absolute value of the Xi^2 torsion angle (ref. 20). For
! example, the ring carbons of Phe and Tyr lying in the plane with the smaller Xi^2 torsion angle are designated as CD1 and CE1
!
! Special action in the case of ARG
!
do i = 1, numat
  if (txtatm(i)(18:20) == "ARG") then
    if (txtatm(i)(14:16) == "NH1") then
      swap = .false.
!
! Check that the side-chain -NH1 is cis to CD
!
      do j = 1, nbonds(i)
        k = ibonds(j,i)
        if (txtatm(k)(14:15) == "CZ") then
          do l = 1, nbonds(k)
            m = ibonds(l,k)
            if (txtatm(m)(14:15) == "NE") then
              do n = 1, nbonds(m)
                ii = ibonds(n,m)
                if (txtatm(ii)(14:15) == "CD") then
                  call dihed (coord, i, k, m, ii, torsion)
                  if (torsion > pi) torsion = torsion - 2*pi
                  swap = (torsion > pi*0.5d0)
!
!  NH1 is trans, it should be cis, so swap labels
!
                  go to 10
                end if
              end do
            end if
          end do
        end if
10      if (swap) then
          do l = 1, nbonds(k)
            m = ibonds(l,k)
            if (txtatm(m)(14:16) == "NH2") exit
          end do
          txtatm(i)(16:16) = "2"
          do l = 1, nbonds(i)
            n = ibonds(l,i)
            if (txtatm(n)(14:15) == "HH") txtatm(n)(16:16) = "2"
          end do
          txtatm(m)(16:16) = "1"
          do l = 1, nbonds(m)
            n = ibonds(l,m)
            if (txtatm(n)(14:15) == "HH") txtatm(n)(16:16) = "1"
          end do
          exit
        end if
      end do
    end if
  end if
end do
allocate (used(numat))
used = .false.
do i = 1, numat
  if (used(i)) cycle
  if (nat(i) == 1) then
    do j = 1, 20
      if (txtatm(i)(18:20) == tyres(j)) exit
    end do
    if (j == 21) cycle
!
!  Convert hydrogen atom from old-style to Greek letter only
!
    do j = 12, 15
      if (txtatm(i)(j:j) == "H") exit
    end do
    do k = j, 17
      if (txtatm(i)(k:k) == " ") exit
    end do
    line = " H"//txtatm(i)(j + 1:k - 1)
    txtatm(i)(13:16) = trim(line)
  end if
end do
do i = 1, numat
  if (nat(i) == 1) then
    do j = 1, 20
      if (txtatm(i)(18:20) == tyres(j)) exit
    end do
    if (j == 21) then
!
! Not one of the 20 amino-acids, so treat it as a generic atom.
! If a number is present at the start of the atom symbol, move the number to the end of the atom symbol.
!
      do j = 13, 15
        if (txtatm(i)(j:j) /= " ") exit
      end do
      if (txtatm(i)(j:j) > "0" .and. txtatm(i)(j:j) < "9") then
        line = txtatm(i)(j:16)
        if (len_trim(line) == 4) then
          txtatm(i)(j:16) = trim(line(2:))
          j = j + 3
          txtatm(i)(j:j) = line(1:1)
        else
          txtatm(i)(j:15) = " "//trim(line(2:))
          j = j + len_trim(line) 
          txtatm(i)(j:j) = line(1:1)
        end if
      end if
      cycle
    end if
!
! Identify the atom that the hydrogen atom is attached to
!
    j = ibonds(1,i)
    l = 0
    do k = 1, nbonds(j)
      if (nat(ibonds(k, j)) == 1) l = l + 1
    end do
    if (i == 21) then
      continue
    end if
    if (nbonds(j) == 4) then
!
! "j" is an atom that is bonded to four atoms, one of which is hydrogen atom "i"
!
      if (l == 3) then
!
! Three hydrogen atoms attached to atom "j"
!
        l = 0
        do m = 1, 4
          n = ibonds(m,j)
          if (nat(n) == 1 .and. .not. used(n)) then
            l = l + 1
            do k = 15, 17
              if (txtatm(n)(k:k) == " ") exit
            end do
            if (k == 17) then
              line = txtatm(n)(14:16)
              txtatm(n)(13:16) = line(1:4)
              k = 16
            end if
            used(n) = .true.
            write(txtatm(n)(k:k), '(i1)')l
          end if
        end do
      else if (l == 2) then
        call two_atoms(i,j)
      end if
    else if (nbonds(j) == 3 .and. l == 2) then
!
! "j" is an atom that is bonded to three atoms, two of these are hydrogen, and one is hydrogen atom "i"
!
      call two_atoms(i,j)
    end if
  end if
end do
return
end subroutine PDB3
subroutine two_atoms(i, j)
!
! two_atoms assignes the PDB-3.n labels to two hydrogen atoms attached to atom "j"
!
! On input:
!          "i" = atom-number of a hydrogen atom
!          "j" = atom-number of the atom that "i" is attached to
! On output:
!   txtatm(i) is updated with the PDB-3.n label
!
use common_arrays_C, only : txtatm, coord, nbonds, ibonds
use funcon_C, only : pi
USE molkst_C, ONLY : line, numat
implicit none
integer, intent (in) :: i, j
integer :: jj, k, l, m, n, ii, priority(4), order(3)
double precision :: torsion
logical :: swap
character :: types(24)*1, nos(3)*1, no1*1
data types / "A", "B", "G", "D", "E", "Z", "H", "T", "I", "K", "L", "M", &
   & "N", "X", "O", "P", "R", "S", "T", "U", "F", "C", "Y", "W" /
data nos / "1", "2", "3" /
! The hydrogen atoms of the side-chain -NH2 groups of Asn, Gln, and Arg are distinguished
! by numbers (1 or 2) on the basis of their relationship (cis or trans, respectively) to the heavy atom three bonds closer to the
! main chain (CB for Asn, CG for Gln, NE for Arg). Thus, each -NH2 hydrogen of Arg is distinguished by two numbers, the first
! indicating the nitrogen to which it is attached and the second indicating the stereochemistry of the hydrogen itself.
!
! Two hydrogen atoms attached to atom "j"
  if (txtatm(j)(18:20) == "ARG" .and. txtatm(j)(14:15) == "NH") then
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
!
! If a hydrogen atom on atom j, i.e., NH, has already been defined, then
! automatically assign the current hydrogen atom, atom i.
!
      if (txtatm(k)(13:14) == "HH") then
        if (k < i) then
          line = txtatm(i)(14:16)
          if (txtatm(k)(16:16) == "1") then
            txtatm(i)(13:16) = line(1:3)//"2"
          else
            txtatm(i)(13:16) = line(1:3)//"1"
          end if
          return
        end if
      end if
    end do
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
      if (txtatm(k)(14:15) == "CZ") then
        do l = 1, nbonds(k)
          m = ibonds(l,k)
          if (txtatm(m)(14:15) == "NE") then
            call dihed (coord, i, j, k, m, torsion)
            if (torsion > pi) torsion = torsion - 2*pi
            swap = (abs(torsion) > pi*0.5d0)
            go to 10
          end if
        end do
      end if
    end do
10  continue
    if (swap) then
      no1 = "2"
    else
      no1 = "1"
    end if
    line = txtatm(i)(14:16)
    txtatm(i)(13:16) = line(1:3)//no1
    return
  end if
!
  if (txtatm(j)(18:20) == "ASN" .and. txtatm(j)(14:15) == "ND") then
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
!
! If a hydrogen atom on atom j, i.e., ND, has already been defined, then
! automatically assign the current hydrogen atom, atom i.
!
      if (txtatm(k)(13:14) == "HD") then
        if (k < i) then
          line = txtatm(i)(14:16)
          if (txtatm(k)(16:16) == "1") then
            txtatm(i)(13:16) = line(1:3)//"2"
          else
            txtatm(i)(13:16) = line(1:3)//"1"
          end if
          return
        end if
      end if
    end do
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
      if (txtatm(k)(14:15) == "CG") then
        do l = 1, nbonds(k)
          m = ibonds(l,k)
          if (txtatm(m)(14:15) == "CB") then
            call dihed (coord, i, j, k, m, torsion)
            if (torsion > pi) torsion = torsion - 2*pi
            swap = (abs(torsion) > pi*0.5d0)
            go to 20
          end if
        end do
      end if
    end do
20  continue
    if (swap) then
      no1 = "2"
    else
      no1 = "1"
    end if
    line = txtatm(i)(14:16)
    txtatm(i)(13:16) = line(1:3)//no1
    return
  end if
  !
  if (txtatm(j)(18:20) == "GLN" .and. txtatm(j)(14:15) == "NE") then
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
!
! If a hydrogen atom on atom j, i.e., NE, has already been defined, then
! automatically assign the current hydrogen atom, atom i.
!
      if (txtatm(k)(13:14) == "HE") then
        if (k < i) then
          line = txtatm(i)(14:16)
          if (txtatm(k)(16:16) == "1") then
            txtatm(i)(13:16) = line(1:3)//"2"
          else
            txtatm(i)(13:16) = line(1:3)//"1"
          end if
          return
        end if
      end if
    end do
    do jj = 1, nbonds(j)
      k = ibonds(jj,j)
      if (txtatm(k)(14:15) == "CD") then
        do l = 1, nbonds(k)
          m = ibonds(l,k)
          if (txtatm(m)(14:15) == "CG") then
            call dihed (coord, i, j, k, m, torsion)
            if (torsion > pi) torsion = torsion - 2*pi
            swap = (abs(torsion) > pi*0.5d0)
            go to 30
          end if
        end do
      end if
    end do
30  continue
    if (swap) then
      no1 = "2"
    else
      no1 = "1"
    end if
    line = txtatm(i)(14:16)
    txtatm(i)(13:16) = line(1:3)//no1
    return
  end if
!
  if (txtatm(j)(18:20) == "GLY" .and. txtatm(j)(14:15) == "CA") then
!
! Assumed convention: The first hydrogen on CA is HA2, the second is HA3
!
!  Find the first HA atom in this GLY residue. 
!
    do l = j + 1,  numat
      if (txtatm(l)(18:20) == "GLY" .and. txtatm(l)(14:15) == "HA") exit
    end do    
    if (l == i) then
      txtatm(i)(14:16) = "HA2"
    else
      txtatm(i)(14:16) = "HA3"
    end if
    return
  end if

!
  priority = 100
  do k = 1, nbonds(j)
    l = ibonds(k, j)
    if (l == i) then
      m = 25
    else
      jj = 15
      if (txtatm(l)(13:13) /= " ") jj = 14
      if (txtatm(l)(jj:jj) /= " ") then
        do m = 1, 24
          if (txtatm(l)(jj:jj) == types(m)) exit
        end do
      else
        m = 24
      end if
    end if
    priority(k) = m
  end do
!
!  Put the atoms in order.
!
  do k = 1, 2
    l = 110
    do m = 1, 4
      if (priority(m) < l .and. ibonds(m,j) /= i) then
        l = priority(m)
        n = ibonds(m,j)
        ii = m
      end if
    end do
    priority(ii) = 100
    order(k) = n
  end do
!
! Evaluate dihedral i - j - order(1) - order(2)
!
  call dihed (coord, i, j, order(1), order(2), torsion)
  if (torsion > pi) torsion = torsion - 2*pi
  do k = 15, 17
    if (txtatm(i)(k:k) == " ") exit
  end do
  if (k == 17) then
    line = txtatm(i)(14:16)
    txtatm(i)(13:16) = line(1:4)
    k = 16
  end if
  if (torsion < 0.d0) then
    txtatm(i)(k:k) = nos(nbonds(j) - 2)
  else
    txtatm(i)(k:k) = nos(nbonds(j) - 1)
  end if
  return
end subroutine two_atoms
