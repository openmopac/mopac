! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

      subroutine partxy(c34, pq34, w)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numcal, numat, lm61
      use common_arrays_C, only : nfirst, nlast
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision  :: c34(*)
      double precision  :: pq34(lm61)
      double precision  :: w(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(-1:8) :: nb
      integer :: icalcn, iband, kr, ls, ii, ia, ib, lsp, i, &
        j, lsw, k, l, jband, js, jj
      double precision :: aa, sum, bb
      save nb, icalcn
!-----------------------------------------------
!------------------------------------------------------------------
!
!    PARTXY WORKS OUT  IN MNDO FORMALISM THE FIRST 2-INDICES TRANSFO.
!          REQUIRED IN THE COMPUTATION OF 2-ELECTRONS REPULSION OVER M.O
!  INPUT
!     C34   : VECTOR OF THE CURRENT CHARGE DISTRIBUTION BETWEEN TWO M.O.
!  OUTPUT
!     PQ34(PQ) : <P(1),Q(1)|C3(2),C4(2)> WHERE P ,Q  ARE A.O.
!                                          AND C3,C4 ARE M.O.
!                P AND Q RUN IN CANONICAL ORDER OVER THE A.O BELONGING
!                TO AN ATOM 'A' ONLY (BASIC ASSUMPTION OF MNDO SCHEME)
!                AND 'A' RUNS OVER THE ATOMS OF THE SYSTEM.
!     D.L. (DEWAR GROUP) 1986
!----------------------------------------------------------------------
      data nb/ 0, 1, 0, 0, 10, 0, 0, 0, 0, 45/
      data icalcn/ 0/
      if (numcal /= icalcn) then
        icalcn = numcal
      end if
!     KK    : POINTER OF SUPPORTING ATOM, SPARKLES SKIPPED OUT.
 !     KK    : POINTER OF SUPPORTING ATOM, SPARKLES SKIPPED OUT.
   !
   !     LOOP OVER OUTER ATOM A, SPARKLES EXCLUDED.
   !     ------------------------------------------
    iband = 1
    kr = 1
    ls = 0
    pq34 = 0.d0
    do ii = 1, numat
      ia = nfirst(ii)
      ib = nlast(ii)
      if (ib >= ia) then
        ls = ls + iband
        iband = nb(ib-ia)
         !
         !     PQ34(IJ) = <IJ|KL> * C34(KL)  , 1-CENTRE CONTRIBUTIONS.
        lsp = ls - 1
         !
         !        LOOP OVER CHARGE DISTRIBUTION OF INNER ATOMS  B < A .
         !        -----------------------------------------------------
         !        PQ34(IJ)=<IJ|KL>*C34(KL) 2-CENTRES CONTRIBUTIONS.
         !
        jband = 1
        js = 0
        do jj = 1, ii - 1
          js = js + jband
          jband = nb(nlast(jj)-nfirst(jj))
          if (jband /= 0) &
          call formxy (w(kr), kr, pq34(ls), pq34(js), c34(ls), iband,  c34(js), jband)
        end do
        lsp = ls - 1
        do i = ia, ib
          aa = 1.d0
          do j = ia, i
            if (i == j) then
              aa = 0.5d0
            end if
            lsp = lsp + 1
               !
               !    'J' Address IJ in W
               !
            sum = 0.d0
            lsw = ls - 1
            do k = ia, ib
              bb = 1.d0
              do l = ia, k
                if (l == k) then
                  bb = 0.5d0
                end if
                lsw = lsw + 1
                     !
                     !   The term itself
                     !
                sum = sum + c34(lsw) * w(kr) * bb
                kr = kr + 1
              end do
            end do
            pq34(lsp) = pq34(lsp) + sum * aa
          end do
        end do
      end if
    end do
end subroutine partxy
