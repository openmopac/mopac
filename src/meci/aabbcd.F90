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

      double precision function aabbcd (iocca1, ioccb1, iocca2, ioccb2, nmos, &
        xy)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE meci_C, only : ispqr, is, iiloop, jloop
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
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
      integer :: i, j, k, l, m, ij
      double precision :: xr
!-----------------------------------------------
!**********************************************************************
!
! AABBCD EVALUATES THE C.I. MATRIX ELEMENT FOR TWO MICROSTATES DIFFERING
!       BY TWO SETS OF M.O.S. ONE MICROSTATE HAS AN ALPHA ELECTRON
!       IN PSI(I) AND A BETA ELECTRON IN PSI(K) FOR WHICH THE OTHER
!       MICROSTATE HAS AN ALPHA ELECTRON IN PSI(J) AND A BETA ELECTRON
!       IN PSI(L)
!
!**********************************************************************
      do i = 1, nmos
        if (iocca1(i) == iocca2(i)) cycle
        exit
      end do
      do j = i + 1, nmos
        if (iocca1(j) == iocca2(j)) cycle
        exit
      end do
      do k = 1, nmos
        if (ioccb1(k) == ioccb2(k)) cycle
        exit
      end do
      do l = k + 1, nmos
        if (ioccb1(l) == ioccb2(l)) cycle
        exit
      end do
      if (i==k .and. j==l .and. iocca1(i)/=ioccb1(i)) then
        ispqr(iiloop,is) = jloop
        is = is + 1
      end if
      if (iocca1(i) < iocca2(i)) then
        m = i
        i = j
        j = m
      end if
      if (ioccb1(k) < ioccb2(k)) then
        m = k
        k = l
        l = m
      end if
      xr = xy(i,j,k,l)
!#      WRITE(IW,'(4I5,F12.6)')I,J,K,L,XR
!
!   NOW UNTANGLE THE MICROSTATES
!
      ij = 1
      if (i>k .and. j>l .or. i<=k .and. j<=l) ij = 0
      if (i > k) ij = ij + iocca1(k) + ioccb1(i)
      if (j > l) ij = ij + iocca2(l) + ioccb2(j)
      if (i > k) then
        m = i
        i = k
        k = m
      end if
      ij = ij + sum(ioccb1(i:k)+iocca1(i:k))
      if (j > l) then
        m = j
        j = l
        l = m
      end if
      ij = ij + sum(ioccb2(j:l)+iocca2(j:l))
!
!   IJ IN THE PERMUTATION NUMBER, .EQUIV. -1 IF IJ IS ODD.
!
      if (mod(ij,2) == 1) xr = -xr
      aabbcd = xr
      return
      end function aabbcd
