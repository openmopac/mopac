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

      subroutine rotmol(numat, coord, sina, cosa, i, j, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat
      integer , intent(in) :: i
      integer , intent(in) :: j
      double precision , intent(in) :: sina
      double precision , intent(in) :: cosa
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k
      double precision :: buff
!-----------------------------------------------
      call symopr (numat, coord, -1, r)
      do k = 1, 3
        buff = (-sina*r(k,i)) + cosa*r(k,j)
        r(k,i) = cosa*r(k,i) + sina*r(k,j)
        r(k,j) = buff
      end do
      call symopr (numat, coord, 1, r)
      return
      end subroutine rotmol
