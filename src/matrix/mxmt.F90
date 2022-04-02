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

      subroutine mxmt(a, nar, b, nbr, c, ncc)
      implicit none
      integer  :: nar
      integer  :: nbr
      integer  :: ncc
      double precision  :: a(nar,nbr)
      double precision  :: b(ncc,nbr)
      double precision  :: c(nar,ncc)
!
!     MATRIX PRODUCT C(NAR,NCC) = A(NAR,NBR) * (B(NCC,NBR))'
!     ALL MATRICES RECTANGULAR , PACKED.
!
      call dgemm ('N', 'T', nar, ncc, nbr, 1.0D0, a, nar, b, ncc, 0.0D0, c, nar)
      return
      end subroutine mxmt
