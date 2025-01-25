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

      subroutine mult33(fmat, iplace)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE symmetry_C, only : elem
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iplace
      double precision , intent(in) :: fmat(3,3)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k
      double precision, dimension(3,3) :: help
!-----------------------------------------------
!***********************************************************************
!
!   MULT33 performs a eulerian rotation on a symmetry operation.
!
!    The symmetry operation IPLACE, stored in ELEM is subjected
!    to the operation stored in FMAT.
!
!***********************************************************************
      do i = 1, 3
        do j = 1, 3
          help(i,j) = 0.D0
          do k = 1, 3
            help(i,j) = help(i,j) + sum(fmat(i,:)*fmat(j,k)*elem(:,k,iplace))
          end do
        end do
      end do
      elem(:,:,iplace) = help
      return
      end subroutine mult33
