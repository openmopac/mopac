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

      subroutine supdot(s, h, g, n)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(inout) :: s(*)
      double precision , intent(in) :: h(*)
      double precision , intent(in) :: g(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j
      double precision :: sum, gi
!-----------------------------------------------
!     (S)=(H)*(G) WITH  H  IN PACKED FORM (CANONICAL ORDER).
      k = 0
      do i = 1, n
        sum = 0.D0
        do j = 1, i
          sum = sum + g(j)*h(k+j)
        end do
        s(i) = sum
        k = k + i
      end do
      if (n == 1) return
      k = 1
      do i = 2, n
        gi = g(i)
        s(:i-1) = s(:i-1) + h(k+1:i-1+k)*gi
        k = k + i
      end do
      return
      end subroutine supdot
