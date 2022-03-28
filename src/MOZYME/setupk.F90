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

subroutine setupk (nocc1)
    use MOZYME_C, only : icocc, ncf, nncf, kopt
    use molkst_C, only : numat
    implicit none
    integer, intent (in) :: nocc1
    integer :: i, j, k, l
   !********************************************************************
   !
   !   SETUPK DETERMINES WHICH ATOMS NEED TO BE CONSIDERED IN
   !   RE-CONSTRUCTING THE FOCK MATRIX.
   !
   !   THE SET OF ATOMS TO BE USED IS PUT INTO KOPT.  ALL ATOMS USED IN
   !   THE LMOs IN THE SCF CALCULATION ARE PUT INTO KOPT.
   !
   !   (Note:  Perhaps a better test would be 'All atoms that make a
   !   significant contribution to the LMOs in the SCF calculation?'
   !
   !********************************************************************
    kopt = 0
    do i = 1, nocc1
      j = nncf(i)
      do k = 1, ncf(i)
        kopt(icocc(k+j)) = 1
      end do
    end do
!
!   If atom i is to be used in constructing the Fock matrix, then kopt(i) = 1
!   otherwise kopt(i) = 0
!
!   Now compress kopt, so that all the atoms used for the Fock matrix are in
!   order.
!
    l = 0
    do i = 1, numat
      if (kopt(i) == 1) then
        l = l + 1
        kopt(l) = i
      end if
    end do
!
!  Force a zero in after the last atom.  This can be used in finding the end of the atom list.
!
    if (l /= numat) kopt(l+1) = 0
end subroutine setupk
