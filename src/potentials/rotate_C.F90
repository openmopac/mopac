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

      module rotate_C
!
!   Set up an equivalence between the core-core array ccore and the individual
!   terms, for use in ROTATE
!
      double precision :: css1, csp1, cpps1, cppp1, css2, csp2, cpps2, cppp2
      double precision, dimension(4,2) :: ccore = 0.d0
      equivalence (ccore(1,1), css1), (ccore(2,1), csp1), (ccore(3,1), cpps1), &
      & (ccore(4,1), cppp1), (ccore(1,2), css2), (ccore(2,2), csp2), &
      & (ccore(3,2), cpps2), (ccore(4,2), cppp2)
      end module rotate_C
