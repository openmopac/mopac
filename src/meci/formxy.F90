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


      subroutine formxy(w, kr, wca, wcb, ca, na, cb, nb)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numcal
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(inout) :: kr
      integer , intent(in) :: na
      integer , intent(in) :: nb
      double precision , intent(in) :: w(*)
      double precision , intent(inout) :: wca(45)
      double precision , intent(inout) :: wcb(45)
      double precision , intent(in) :: ca(45)
      double precision , intent(in) :: cb(45)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, ij, i, j, kl, k, l, nna, nnb, n1, n2
      double precision :: aa, sum, bb
      integer, dimension (45) :: in

      save icalcn
!-----------------------------------------------
!***********************************************************************
!
!    EACH OF THE NA ELEMENTS OF WCA WILL ADD ON THE NB ELECTROSTATIC
!    TERMS FROM ATOM B IN CB
!
!    EACH OF THE NB ELEMENTS OF WCB WILL ADD ON THE NA ELECTROSTATIC
!    TERMS FROM ATOM A IN CA
!
!    BOTH SUMS WILL INVOLVE THE NA*NB TERMS IN ARRAY W.  ONCE USED,
!    W WILL BE INCREMENTED BY NA*NB.
!
! NA=NUMBER OF ATOMIC ORBITALS ON ATOM 'A'.
! NB=NUMBER OF ATOMIC ORBITALS ON ATOM 'B'.
!
!***********************************************************************
      data icalcn/ 0/
      in(1)  = 1
      in(10) = 4
      in(45) = 9
      if (icalcn /= numcal) icalcn = numcal
      ij = 0
      n1 = 0
      nna = in(na)
      nnb = in(nb)
      do i = 1, nna
        aa = 1.d0
        do j = 1, i
          n1 = n1 + 1
          if (i == j) then
            aa = 0.5d0
          end if
          ij = ij + 1
          sum = 0.d0
          kl = 0
          n2 = 0
          do k = 1, nnb
            bb = 1.d0
            do l = 1, k
              n2 = n2 + 1
              if (k == l) then
                bb = 0.5d0
              end if
              kl = kl + 1
              sum = sum + cb(kl) * w((n1-1)*nb+n2) * bb
            end do
          end do
          wca(ij) = wca(ij) + sum * aa
        end do
      end do
      ij = 0
      n1 = 0
      do i = 1, nnb
        aa = 1.d0
        do j = 1, i
          n1 = n1 + 1
          if (i == j) then
            aa = 0.5d0
          end if
          ij = ij + 1
          sum = 0.d0
          kl = 0
          n2 = 0
          do k = 1, nna
            bb = 1.d0
            do l = 1, k
              n2 = n2 + 1
              if (k == l) then
                bb = 0.5d0
              end if
              kl = kl + 1
              sum = sum + ca(kl) * w((n2-1)*nb+n1) * bb
            end do
          end do
          wcb(ij) = wcb(ij) + sum * aa
        end do
      end do
      kr = kr + na * nb
      return
      end subroutine formxy
