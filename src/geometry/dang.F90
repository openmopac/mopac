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

      subroutine dang(a1, a2, b1, b2, rcos)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(inout) :: a1
      double precision , intent(inout) :: a2
      double precision , intent(inout) :: b1
      double precision , intent(inout) :: b2
      double precision , intent(out) :: rcos
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: zero, anorm, bnorm, sinth, costh
!-----------------------------------------------
!*********************************************************************
!
!    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!
!*********************************************************************
      zero = 1.0D-6
      if (abs(a1) >= zero .or. abs(a2) >= zero) then
        if (abs(b1) >= zero .or. abs(b2) >= zero) then
          anorm = 1.0D0/sqrt(a1**2 + a2**2)
          bnorm = 1.0D0/sqrt(b1**2 + b2**2)
          a1 = a1*anorm
          a2 = a2*anorm
          b1 = b1*bnorm
          b2 = b2*bnorm
          sinth = a1*b2 - a2*b1
          costh = a1*b1 + a2*b2
          costh = min(1.0D0,costh)
          costh = dmax1(-1.0D0,costh)
          rcos = acos(costh)
          if (abs(rcos) >= 4.0D-5) then
            if (sinth > 0.D0) rcos = 6.28318530717959D0 - rcos
            rcos = -rcos
            return
          end if
        end if
      end if
      rcos = 0.0D0
      return
      end subroutine dang
