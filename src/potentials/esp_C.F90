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

      module esp_C
      integer  :: nesp, idip, iz, ipx, isc, is_esp, icd, ipe, npr, ic, ip, ncc
      double precision :: dens, scale, cf, &
      & rms, rrms, dx, dy, dz, den
      double precision, dimension(0:8,821) :: fv
      double precision, dimension(-1:96) :: dex
      double precision, dimension(0:2) :: tf
      integer, dimension (10) :: ixn, iyn, izn, jxn, jyn, jzn
      double precision, dimension (0:7) :: fac

      integer, dimension(:), allocatable :: ind, itemp, ird, indc

      integer, dimension(:,:), allocatable :: iam

      double precision, dimension(:), allocatable :: cc, ex, temp, rad, es, &
      & b_esp, esp_array, cesp, cespml, al, qesp, td, dx_array, dy_array, dz_array, &
      & ptd, pexs, pce, pexpn, pewcx, pewcy, pewcz, pf0, pf1, pf2, pmtd, fc, &
      & qsc, espx, espp

      double precision, dimension(:,:), allocatable :: cespm, cen, potpc, co, &
      & potpt, ovl, cespm2, a2, exs, ce, expn, ewcx, ewcy, ewcz, f0, f1, u_esp, &
      & rnai, rnai1, rnai2, espi, exsr
      logical, dimension(:,:), allocatable :: cequiv
      data tf/33.d0, 37.d0, 41.d0/
      data fac / 1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0 /
!                S    P             D
!                s  x y z   xx yy zz xy xz yz
      data ixn /  0, 4,0,0,  8, 0, 0, 4, 4, 0 /
      data iyn /  0, 0,4,0,  0, 8, 0, 4, 0, 4 /
      data izn /  0, 0,0,4,  0, 0, 8, 0, 4, 4 /
      data jxn /  0, 1,0,0,  2, 0, 0, 1, 1, 0 /
      data jyn /  0, 0,1,0,  0, 2, 0, 1, 0, 1 /
      data jzn /  0, 0,0,1,  0, 0, 2, 0, 1, 1 /
      end module esp_C
