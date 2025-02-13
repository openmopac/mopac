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

      double precision function aabacd (iocca1, ioccb1, iocca2, ioccb2, nmos, xy)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmos
      integer , intent(in) :: iocca1(nmos)
      integer , intent(in) :: ioccb1(nmos)
      integer , intent(in) :: iocca2(nmos)
      integer , intent(in) :: ioccb2(nmos)
      double precision , intent(in) :: xy(nmos,nmos,nmos,nmos)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij, i, j, k, l
      double precision :: sum
!-----------------------------------------------
!**********************************************************************
!
! AABACD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
!       BY TWO ALPHA MOS. ONE MICROSTATE HAS ALPHA ELECTRONS IN
!       M.O.S PSI(I) AND PSI(J) FOR WHICH THE OTHER MICROSTATE HAS
!       ELECTRONS IN PSI(K) AND PSI(L)
!
!**********************************************************************
      ij = 0
      do i = 1, nmos
        if (iocca1(i) >= iocca2(i)) cycle
        exit
      end do
      do j = i + 1, nmos
        if (iocca1(j) < iocca2(j)) exit
        ij = ij + iocca2(j) + ioccb2(j)
      end do
      do k = 1, nmos
        if (iocca1(k) <= iocca2(k)) cycle
        exit
      end do
      do l = k + 1, nmos
        if (iocca1(l) > iocca2(l)) exit
        ij = ij + iocca1(l) + ioccb1(l)
      end do
      ij = ij + ioccb2(i) + ioccb1(k)
      sum = xy(i,k,j,l) - xy(i,l,k,j)
      if (mod(ij,2) == 1) sum = -sum
      aabacd = sum
      return
      end function aabacd
