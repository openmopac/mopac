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

      subroutine mult33(fmat, iplace)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE symmetry_C, only : elem
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iplace
      double precision , intent(in) :: fmat(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k
      double precision, dimension(3,3) :: help
!-----------------------------------------------
!***********************************************************************
!
!   MULT33 performs a eulerian rotation on a symmetry operation.
!
!    The symmetry operation IPLACE, stored in ELEM is subjected
!    to the operation stored in FMAT.
!
!***********************************************************************
      do i = 1, 3
        do j = 1, 3
          help(i,j) = 0.D0
          do k = 1, 3
            help(i,j) = help(i,j) + sum(fmat(i,:)*fmat(j,k)*elem(:,k,iplace))
          end do
        end do
      end do
      elem(:,:,iplace) = help
      return
      end subroutine mult33
