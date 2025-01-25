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

      subroutine chrge(p, q)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only:  nfirst, nlast
      use molkst_C, only: numat, norbs, mozyme
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      double precision , intent(out) :: q(numat)
      double precision , intent(in) :: p((norbs*(norbs+1))/2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, ia, ib, j
!-----------------------------------------------
!***********************************************************************
!
!      CHRGE STORES IN Q THE TOTAL ELECTRON DENSITIES ON THE ATOMS
!
!      ON INPUT P      = DENSITY MATRIX
!
!      ON OUTPUT Q     = ATOM ELECTRON DENSITIES
!
!***********************************************************************
      if (mozyme) then
        call chrge_for_MOZYME(p, q)
        return
      end if
      k = 0
      do i = 1, numat
        ia = nfirst(i)
        ib = nlast(i)
        q(i) = 0.D0
        do j = ia, ib
          k = k + j
          q(i) = q(i) + p(k)
        end do
      end do
      return
      end subroutine chrge

! TODO: to be parallelized
