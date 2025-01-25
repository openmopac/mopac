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

      subroutine swap(c, n, mdim, nocc, ifill)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE chanel_C, only : iw
      use iter_C, only : psi => pbold, stdpsi => pbold2
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: mdim
      integer , intent(in) :: nocc
      integer , intent(inout) :: ifill
      double precision , intent(inout) :: c(mdim,mdim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  i, jfill
      double precision :: sum, summax, x
      double precision, external :: ddot

!-----------------------------------------------
!******************************************************************
!
!        SWAP ENSURES THAT A NAMED MOLECULAR ORBITAL IFILL IS FILLED
! ON INPUT
!          C = EIGENVECTORS IN A MDIM*MDIM MATRIX
!          N = NUMBER OF ORBITALS
!          NOCC = NUMBER OF OCCUPIED ORBITALS
!          IFILL = FILLED ORBITAL
!******************************************************************
      if (ifill <= 0) then
!
!     WE NOW DEFINE THE FILLED ORBITAL
!
        if (allocated(psi)) deallocate(psi)
        if (allocated(stdpsi)) deallocate(stdpsi)
        allocate (psi(n), stdpsi(n))
        ifill = -ifill
        stdpsi(:n) = c(:n,ifill)
        psi(:n) = c(:n,ifill)
        return
      end if
      sum = 0.D0
      sum = ddot(n,psi(:n),1,c(:n,ifill),1)
!
      if (abs(sum) <= 0.707106781187D0) then
!
!     IFILL HAS MOVED!
!
        summax = 0.D0
        do ifill = 1, n
          sum = 0.D0
          do i = 1, n
            sum = sum + stdpsi(i)*c(i,ifill)
          end do
          sum = abs(sum)
          if (sum > summax) jfill = ifill
          summax = dmax1(sum,summax)
          if (sum > 0.707106781187D0) go to 90
        end do
        do ifill = 1, n
          sum = 0.D0
          do i = 1, n
            sum = sum + psi(i)*c(i,ifill)
          end do
          sum = abs(sum)
          if (sum > summax) jfill = ifill
          summax = dmax1(sum,summax)
          if (sum > 0.7071D0) go to 90
        end do
        write (iw, 80) summax, jfill
   80   format(/,' CAUTION !!! SUM IN SWAP VERY SMALL, SUMMAX =',f10.5,&
          ' JFILL=',i3)
        ifill = jfill
      end if
   90 continue
      if (ifill <= nocc) return
!
!    ITS EMPTY, SO SWAP IT WITH THE HIGHEST FILLED
!
      do i = 1, n
        x = c(i,nocc)
        c(i,nocc) = c(i,ifill)
        c(i,ifill) = x
      end do
      return
      end subroutine swap
