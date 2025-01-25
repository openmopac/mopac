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

      subroutine sympop(h, i, iskip, deldip)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : nsym, ipo
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i
      integer , intent(out) :: iskip
      double precision  :: h(*)
      double precision  :: deldip(3,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j
!-----------------------------------------------
      do j = 1, nsym
        if (ipo(i,j) >= i) cycle
        call symh (h, deldip, i, j, ipo)
        iskip = 3
!  atom ipo(i,j) is suitable for transition dipole calc'n
!
!            K=I*3-2
!            WRITE(IW,*)' Transition dipole after symmetry operation'
!            WRITE(IW,'(3f12.5)')(deldip(l,k),l=1,3)
!            WRITE(IW,'(3f12.5)')(deldip(l,k+1),l=1,3)
!            WRITE(IW,'(3f12.5)')(deldip(l,k+2),l=1,3)
!
!   INSERT DELDIP ROTATION HERE
!
        go to 20
      end do
      iskip = 0
   20 continue
      return
      end subroutine sympop
