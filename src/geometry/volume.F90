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

      double precision function volume (vecs, ndim)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ndim
      double precision , intent(in) :: vecs(3,ndim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: a, b, gamma, sing
!-----------------------------------------------
!**********************************************************************
!
! VOLUME RETURNS (A) THE VOLUME OF A UNIT CELL IF NDIM=3
!                (B) THE AREA   OF A UNIT CELL IF NDIM=2
!                (C) THE LENGTH OF A UNIT CELL IF NDIM=1
!
!  ON INPUT VECS = ARRAY OF POINTS MARKING THE ENDS OF UNIT CELL VECTORS
!                  THE ORIGIN BEING (0.0, 0.0, 0.0)
!
!**********************************************************************
      a = sqrt(vecs(1,1)**2+vecs(2,1)**2+vecs(3,1)**2)
!
!    CASE 1: SYSTEM IS A POLYMER
!
      if (ndim == 1) then
        volume = a
        return
      end if
      b = sqrt(vecs(1,2)**2+vecs(2,2)**2+vecs(3,2)**2)
      gamma = sqrt((vecs(1,1)-vecs(1,2))**2+(vecs(2,1)-vecs(2,2))**2+(vecs(3,1)&
        -vecs(3,2))**2)
!
!     SING = SIN OF ANGLE BETWEEN FIRST AND SECOND VECTORS
!            OBTAINED FROM COSINE RULE C**2=A**2+B**2-2*A*B*COS(C)
!
      sing = sqrt(1.D0 - ((a*a + b*b - gamma*gamma)/(2.D0*a*b))**2)
!
!    CASE 2: SYSTEM IS A LAYER STRUCTURE
!
      if (ndim == 2) then
!
!   AREA OF A PARALLELOGRAM = BASE * HEIGHT
!
        volume = a*b*sing
        return
      end if
!
!    CASE 3: SYSTEM IS A SOLID
!
         !
         !      VOLUME = A.(B x C) IN VECTOR NOTATION
         !
        volume = Abs (((vecs(2, 1)*vecs(3, 2)-vecs(3, 1)*vecs(2, 2))*vecs(1, 3) + &
                       (vecs(3, 1)*vecs(1, 2)-vecs(1, 1)*vecs(3, 2))*vecs(2, 3) + &
                       (vecs(1, 1)*vecs(2, 2)-vecs(2, 1)*vecs(1, 2))*vecs(3, 3)))
      return
      end function volume
