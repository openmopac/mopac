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

      subroutine mat33(a, b, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!  A subroutine that will multiply two 3 by 3 matricies in the following
!     fashion:    C = A(transpose) B A
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(9)
      double precision , intent(in) :: b(9)
      double precision , intent(out) :: c(9)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision, dimension(9) :: t
!-----------------------------------------------
!
!
      t(1) = b(1)*a(1) + b(2)*a(4) + b(3)*a(7)
      t(2) = b(1)*a(2) + b(2)*a(5) + b(3)*a(8)
      t(3) = b(1)*a(3) + b(2)*a(6) + b(3)*a(9)
      t(4) = b(4)*a(1) + b(5)*a(4) + b(6)*a(7)
      t(5) = b(4)*a(2) + b(5)*a(5) + b(6)*a(8)
      t(6) = b(4)*a(3) + b(5)*a(6) + b(6)*a(9)
      t(7) = b(7)*a(1) + b(8)*a(4) + b(9)*a(7)
      t(8) = b(7)*a(2) + b(8)*a(5) + b(9)*a(8)
      t(9) = b(7)*a(3) + b(8)*a(6) + b(9)*a(9)
!
      c(1) = a(1)*t(1) + a(4)*t(4) + a(7)*t(7)
      c(2) = a(1)*t(2) + a(4)*t(5) + a(7)*t(8)
      c(3) = a(1)*t(3) + a(4)*t(6) + a(7)*t(9)
      c(4) = a(2)*t(1) + a(5)*t(4) + a(8)*t(7)
      c(5) = a(2)*t(2) + a(5)*t(5) + a(8)*t(8)
      c(6) = a(2)*t(3) + a(5)*t(6) + a(8)*t(9)
      c(7) = a(3)*t(1) + a(6)*t(4) + a(9)*t(7)
      c(8) = a(3)*t(2) + a(6)*t(5) + a(9)*t(8)
      c(9) = a(3)*t(3) + a(6)*t(6) + a(9)*t(9)
!
      return
      end subroutine mat33
