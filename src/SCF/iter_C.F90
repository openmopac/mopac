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

      module iter_C
      double precision, dimension(:), allocatable :: pold, pold2, pold3, pbold, &
      pbold2, pbold3,  pgasa, pgasb, psona, psonb
!
! Arrays used in interp.F90
!
      double precision, dimension(:), allocatable :: h_ai, vecl_ai, h_bi, vecl_bi
      double precision, dimension(:,:), allocatable :: vec_ai, fock_ai, &
       & p_ai, vec_bi, fock_bi, p_bi

      end module iter_C
