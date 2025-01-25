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

      subroutine mtxmc(a, nar, b, nbr, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nar
      integer  :: nbr
      double precision  :: a(nbr,nar)
      double precision  :: b(nbr,nar)
      double precision  :: c(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, i
!-----------------------------------------------
!     MATRIX PRODUCT C(NAR,NAR) = (A(NBR,NAR))' * B(NBR,NAR)
!     A AND B RECTANGULAR , PACKED,
!     C LOWER LEFT TRIANGLE ONLY, PACKED IN CANONICAL ORDER.
!  NOTE ... THIS IS THE BEST VERSION ON CRAY 1.
      l = 1
      do i = 1, nar
        call mxm (a(1,i), 1, b, nbr, c(l), i)
        l = l + i
      end do
      return
      end subroutine mtxmc
