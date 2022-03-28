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

      subroutine mtxmc(a, nar, b, nbr, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nar
      integer  :: nbr
      double precision  :: a(nbr,nar)
      double precision  :: b(nbr,nar)
      double precision  :: c(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, i
!-----------------------------------------------
!     MATRIX PRODUCT C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR)
!     A AND B RECTANGULAR , PACKED,
!     C LOWER LEFT TRIANGLE ONLY, PACKED IN CANONICAL ORDER.
!  NOTE ... THIS IS THE BEST VERSION ON CRAY 1.
      l = 1
      do i = 1, nar
        call mxm (a(1,i), 1, b, nbr, c(l), i)
        l = l + i
      end do
      return
      end subroutine mtxmc
