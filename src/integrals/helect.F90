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

      double precision function helect (n, p, h, f)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: p(*)
      double precision , intent(in) :: h(*)
      double precision , intent(in) :: f(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, nn, i, jj
      double precision :: ed, ee
!-----------------------------------------------
!***********************************************************************
!
!    SUBROUTINE CALCULATES THE ELECTRONIC ENERGY OF THE SYSTEM IN EV.
!
!    ON ENTRY N = NUMBER OF ATOMIC ORBITALS.
!             P = DENSITY MATRIX, PACKED, LOWER TRIANGLE.
!             H = ONE-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
!             F = TWO-ELECTRON MATRIX, PACKED, LOWER TRIANGLE.
!    ON EXIT
!        HELECT = ELECTRONIC ENERGY.
!
!    NO ARGUMENTS ARE CHANGED.
!
!***********************************************************************
      ed = 0.0D00
      ee = 0.0D00
      k = 0
      nn = n + 1
      do i = 2, nn
        k = k + 1
        jj = i - 1
        ed = ed + p(k)*(h(k)+f(k))
        if (i == nn) cycle
        if (jj > 0) then
          ee = ee + sum(p(k+1:jj+k)*(h(k+1:jj+k)+f(k+1:jj+k)))
          k = jj + k
        end if
      end do
      ee = ee + .5D00*ed
      helect = ee
      return
!
      end function helect
