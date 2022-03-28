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

      module overlaps_C
      double precision :: &
      cutof1,   &  !  For distances beyond cutof1 overlaps are set to zero
      cutof2       !  For distances beyond cutof2, NDDO is replaced by point-charges
      double precision, dimension(60,6) :: ccc, zzz
      double precision, dimension(6,6,2) :: allc, allz
      integer :: isp, ips
      double precision, dimension(7) :: a, b
      double precision :: sa, sb
      double precision, dimension(0:17) :: fact    ! Factorials:  fact(n) = n!
      data fact/ 1.d0, 1.D0, 2.D0, 6.D0, 24.D0, 120.D0, 720.D0, 5040.D0, 40320.D0, &
        362880.D0, 3628800.D0, 39916800.D0, 479001600.D0, 6227020800.D0, &
        8.71782912D10, 1.307674368D12, 2.092278989D13, 3.556874281D14/
      end module overlaps_C
