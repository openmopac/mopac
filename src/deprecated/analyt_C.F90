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

      module analyt_C
!
!  Variables used in STO-6G expansion of Slater orbitals
      implicit none
      integer, dimension(107) :: nztype
      integer, dimension(30) :: mtype
      integer :: ltype
      double precision, dimension(16) :: ds
      double precision, dimension(22) :: dg
      double precision, dimension(100) :: dr
      double precision, dimension(3) :: tdx, tdy, tdz
      double precision, dimension(22) :: g
      double precision, dimension(3) :: tx, ty, tz
      end module analyt_C
