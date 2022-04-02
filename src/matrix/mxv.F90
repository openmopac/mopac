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

      subroutine mxv(a, nar, vecx, nbr, vecy)
      implicit none
      integer, parameter :: incy = 1, incx = 1
      integer  :: nar
      integer  :: nbr
      double precision  :: a(nar,nbr)
      double precision  :: vecx(nbr)
      double precision  :: vecy(nar)
!
!     RECTANGULAR MATRIX-VECTOR PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.

      call dgemv ('N', nar, nbr, 1.0d0, a, nar, vecx, incx, 0.0d0, vecy, incy)
      return
      end subroutine mxv
