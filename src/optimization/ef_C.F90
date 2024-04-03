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

      module ef_C
      integer :: ef_mode, nstep, negreq, iprnt, iloop
      double precision :: ddx, rmin, rmax, omin = 0.d0, xlamd, xlamd0, skal, x0, x1, x2
      double precision, dimension(:), allocatable :: oldf, d, vmode, pmat, uc, hessc
      double precision, dimension(:,:), allocatable :: hess, bmat, u, oldhss, oldu, alparm
      end module ef_C
      module drc_C
        double precision, dimension(:), allocatable :: vref, vref0, now
        double precision, dimension(:,:), allocatable :: allxyz, allvel, xyz3, vel3, allgeo, geo3, georef
        double precision, dimension(:), allocatable :: parref
        double precision :: time
      end module drc_C
      module derivs_C
        double precision, dimension(:), allocatable :: wmat, hmat, fmat
        double precision, dimension(:,:), allocatable  :: b, ab, fb
          double precision, dimension(:), allocatable :: aidref, work2
      end module derivs_C
