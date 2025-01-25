! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
