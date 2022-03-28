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

      subroutine symopr(numat, coord, jump, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: numat
      integer , intent(in) :: jump
      double precision , intent(inout) :: coord(3,numat)
      double precision , intent(in) :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      double precision, dimension(3) :: help
!-----------------------------------------------
!***********************************************************************
!
!   SYMOPR performs the symmetry operation stored in the 3 by 3 array
!          R on the coordinates.   If JUMP is +1, the operation is
!          X.R, if JUMP is -1, the operation is R.X
!
!   On output, COORD = COORD.R or R.COORD
!
!***********************************************************************
      if (jump >= 0) then
        do i = 1, numat
          help = coord(:,i)
          do j = 1, 3
            coord(j,i) = 0.D0
            coord(j,i) = coord(j,i) + sum(r(:,j)*help)
          end do
        end do
        return
      end if
      do i = 1, numat
        help = coord(:,i)
        do j = 1, 3
          coord(j,i) = 0.D0
          coord(j,i) = coord(j,i) + sum(r(j,:)*help)
        end do
      end do
      return
      end subroutine symopr
