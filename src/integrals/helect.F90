! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
