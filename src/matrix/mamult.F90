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

      subroutine mamult(a, b, c, n, one)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: one
      double precision , intent(in) :: a(*)
      double precision , intent(in) :: b(*)
      double precision , intent(inout) :: c(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, i, ii, j, jj, k, kk
      double precision :: sum
!-----------------------------------------------
!***********************************************************************
!
!   MAMULT MULTIPLIES A BY B AND PUTS THE RESULT IN C
!
!***********************************************************************
      l = 0
      do i = 1, n
        ii = ((i - 1)*i)/2
        do j = 1, i
          jj = ((j - 1)*j)/2
          l = l + 1
          sum = 0.D0
          do k = 1, j
            sum = sum + a(ii+k)*b(jj+k)
          end do
          do k = j + 1, i
            sum = sum + a(ii+k)*b(((k-1)*k)/2+j)
          end do
          do k = i + 1, n
            kk = (k*(k - 1))/2
            sum = sum + a(kk+i)*b(kk+j)
          end do
          c(l) = sum + one*c(l)
        end do
      end do
      return
      end subroutine mamult
