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

      double precision function dot (x, y, n)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: x(n)
      double precision , intent(in) :: y(n)
!For MOPAC BLAS
      double precision, external :: ddot
!

!-----------------------------------------------
!***********************************************************************
!
!   DOT FORMS THE SCALAR PRODUCT OF TWO VECTORS.
!
!   ON INPUT     X   =    FIRST VECTOR, OF LENGTH N.
!                Y   =    SECOND VECTOR, OF LENGTH N.
!
!   ON RETURN    DOT =    DOT PRODUCT OF X AND Y.
!
!***********************************************************************
!For MOPAC BLAS
!     dot = dot_product(x(:n),y(:n))
      dot = ddot(n,x(1:n),1,y(1:n),1)
!
      return
      end function dot
