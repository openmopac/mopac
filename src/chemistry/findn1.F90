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

subroutine findn1 (n1, ioptl, io, delta_res)
!**********************************************************************
!
!  Locate the nitrogen atom of a peptide linkage.
!
!  On input, ioptl is a logical array of atoms which might contain a
!  peptide linkage.
!
!  On output:
!
!  n1 is the atom-number of the nitrogen atom of the peptide
!  or zero, if a nitrogen could not be found.
!
!  io is the oxygen of the "C=O" of the peptide linkage.
!
!  delta_res is the change in residue number caused by moving along
!  the chain to the N end.  This is normally zero.
!
!**********************************************************************

    use molkst_C, only: numat, isok
    use common_arrays_C, only : nat, nbonds, ibonds
    implicit none
    integer, intent (out) :: n1, io, delta_res
    logical, dimension (numat), intent (in) :: ioptl
    logical :: peptide, is_ring, ok = .true.
!
    integer :: i, iatom, ii, ires_loc, j, k, l, jj, kk, ll, mm, last_n1
    logical, external :: ring_6
    integer, external :: nheavy
    save :: ok
!
    n1 = 0
    delta_res = 0
    last_n1 = 0
!
!   FIND A N-C-C=O SYSTEM
!
    outer_loop: do iatom = 1, numat
      if ( .not. ioptl(iatom)) then
!
!   START WITH OXYGEN OF C=O
!
        if (nat(iatom) == 8 .and. nbonds(iatom) == 1) then
          do i = 1, nbonds(iatom)
            if (nat(ibonds(i, iatom)) == 6) then
!
!    THERE IS A CARBON ATTACHED TO THE OXYGEN
!
              l = ibonds(i, iatom)
      loop_k: do k = 1, nbonds(l)
                if (nat(ibonds(k, l)) == 6) then
!
!  THERE IS A CARBON ATTACHED TO THE CARBON ATTACHED
!  TO THE OXYGEN
!
                  j = ibonds(k, l)
                  if (nheavy(j) == 4) cycle loop_k
                  if (j == 4159) then
                   j = j
                   end if
!
! Make sure that C(alpha) is not part of a six-membered ring
!
                  do jj = 1, nbonds(j)
                    kk = ibonds(jj,j)
                    if (kk /= l) then
                      if (ring_6(j, l, kk)) cycle loop_k
                    end if
                    if (nat(kk) /= 6 .and. nat(kk) /= 1 .and. nat(kk) /=7) cycle loop_k
                  end do
                  do ii = 1, nbonds(j)
                    if (nheavy(ibonds(ii, j)) < 4 .and. nat(ibonds(ii, j)) == 7 .and. .not. ioptl(ibonds(ii, j))) then

!
!  THERE IS A NITROGEN ATTACHED TO THE CARBON ATTACHED TO THE
!  CARBON ATTACHED TO THE OXYGEN
!
!
!   YES, A PEPTIDE LINKAGE HAS BEEN FOUND.
!   THE NITROGEN IS N1
!
                      n1 = ibonds(ii, j)
!
!  Make sure it's not part of a ring
!
                      is_ring = .false.
                      do jj = 1, nbonds(n1)
                        kk = ibonds(jj,n1)
                        do ll = 1, nbonds(n1)
                          mm = ibonds(ll,n1)
                          if (nat(kk) > 1 .and. nat(mm) > 1 .and. kk > mm) then
                            if (ring_6(n1, mm, kk)) is_ring = .true.
                          end if
                        end do
                      end do
                      if ( is_ring) then
                        n1 = 0
                      else
                        n1 = n1
                        exit outer_loop
                      end if
                    end if
                  end do
                end if
              end do loop_k
            end if
          end do
        end if
      end if
    end do outer_loop
    io = iatom
    if (n1 == 0) return
!
!  N1 IS A NITROGEN IN AN O=C-C-N SYSTEMS.
!
    if (nat(n1) == 7) last_n1 = n1
    j = 0
    do i = 1, nbonds(n1)
      if (nat(ibonds(i, n1)) /= 1) then
        j = j + 1
      end if
    end do
    if (j == 1) return
!
!  AT THIS POINT, THE SEQUENCE O-C-C-N HAS BEEN IDENTIFIED.
!  NOW TO MOVE ALONG BACKBONE TO THE START OF THE PROTEIN
!
    peptide = .false.
    do ires_loc = 1, 1000000
      do i = 1, nbonds(n1)
        j = ibonds(i, n1)
        if (nat(j) == 6) then
          do k = 1, nbonds(j)
            if (nat(ibonds(k, j)) == 8) go to 1000
          end do
          cycle
!
!   J IS THE CARBON OF C=O
!
1000      do k = 1, nbonds(j)
            if (nat(ibonds(k, j)) == 6) then
              peptide = .true.
              go to 1010
            end if
          end do
        end if
      end do
      cycle
1010  j = ibonds(k, j)
      do k = 1, nbonds(j)
        if (nat(ibonds(k, j)) == 7) go to 1020
      end do
      if (peptide) then
!
!  The residue at the start does not contain a nitrogen atom.
!  So, temporarily, artificially convert a hydrogen into a nitrogen.
!  This will allow the second residue to be correctly identified.
!  (The first residue might be almost correctly identified)
!
        do k = 1, nbonds(j)
          if (nat(ibonds(k, j)) == 1) n1 = ibonds(k, j)
        end do
        ok = (ok .and. n1 == ibonds(k, j))
        isok = ok
      end if
      exit
1020  if (n1 == ibonds(k, j)) exit
      if (.not. ioptl(ibonds(k, j))) then
        n1 = ibonds(k, j)
        delta_res = delta_res + 1
      end if
      j = 0
      do ii = 1, nbonds(n1)
        if (nat(ibonds(ii, n1)) /= 1) then
          j = j + 1
        end if
      end do
      if (j == 1) exit
    end do
    if (nat(n1) /= 7) then
      n1 = last_n1
    end if
    return
end subroutine findn1
logical function ring_6(atom_1, atom_2, atom_6)
!
!  Given three adjacent atoms, atom_2 - atom_1 - atom_6, determine if they belong
!  to an aromatic ring.  If so, "aromatic" return .true.
!
    use common_arrays_C,  only : nbonds,  ibonds
    implicit none
    integer, intent (in) :: atom_1, atom_2, atom_6
    integer :: atom_3, atom_5, i, j, k, l, nii, njj, kka, atom_4
      ring_6 = .false.
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
              if (kka == atom_4) then
                if (atom_1 == atom_3 .or. atom_1 == atom_4 .or. atom_1 == atom_5) exit
                if (atom_2 == atom_3 .or. atom_2 == atom_4 .or. atom_2 == atom_5) exit
                if (atom_3 == atom_4 .or. atom_3 == atom_5 .or. atom_3 == atom_6) exit
                if (atom_4 == atom_5 .or. atom_4 == atom_6 .or. atom_5 == atom_6) exit
!
!  Six-membered ring, with atoms numbered atom_1, atom_2, atom_3, atom_4 = kka, atom_5, atom_6
!

                ring_6 = .true.
                return
              end if
            end do
          end do
        end do
      end do
      return
  end function ring_6
