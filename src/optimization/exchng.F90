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

      subroutine exchng(a, b, c, d, x, y, t, q, n)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: a
      double precision , intent(out) :: b
      double precision , intent(in) :: c
      double precision , intent(out) :: d
      double precision , intent(in) :: t
      double precision , intent(out) :: q
      double precision , intent(in) :: x(*)
      double precision , intent(out) :: y(*)
!-----------------------------------------------
!********************************************************************
!
! THE CONTENTS OF A, C, T, AND X ARE STORED IN B, D, Q, AND Y!
!
!   THIS IS A DEDICATED ROUTINE, IT IS CALLED BY LINMIN AND LOCMIN ONLY.
!
!********************************************************************
      b = a
      d = c
      q = t
      y(:n) = x(:n)
      return
!
      end subroutine exchng
