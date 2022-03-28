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

      module funcon_C
      double precision :: &
      &  fpc_1,     & !
      &  fpc_2,     & !
      &  a0,        & !
      &  ev,        & !
      &  fpc_5,     & !
      &  fpc_6,     & !
      &  fpc_7,     & !
      &  fpc_8,     & !
      &  fpc_9,     & !
      &  fpc_10
      double precision, dimension (10) :: fpc = 0.d0
      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: twopi = 2.0d0 * pi
      equivalence (fpc(1), fpc_1), (fpc(2), fpc_2), (fpc(3), a0), (fpc(4), ev), &
      & (fpc(5), fpc_5), (fpc(6), fpc_6), (fpc(7), fpc_7), (fpc(8), fpc_8), &
      & (fpc(9), fpc_9), (fpc(10), fpc_10)
      end module funcon_C
