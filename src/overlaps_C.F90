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
