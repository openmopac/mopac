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

      subroutine rotmol(numat, coord, sina, cosa, i, j, r)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat
      integer , intent(in) :: i
      integer , intent(in) :: j
      double precision , intent(in) :: sina
      double precision , intent(in) :: cosa
      double precision  :: coord(3,numat)
      double precision  :: r(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k
      double precision :: buff
!-----------------------------------------------
      call symopr (numat, coord, -1, r)
      do k = 1, 3
        buff = (-sina*r(k,i)) + cosa*r(k,j)
        r(k,i) = cosa*r(k,i) + sina*r(k,j)
        r(k,j) = buff
      end do
      call symopr (numat, coord, 1, r)
      return
      end subroutine rotmol
