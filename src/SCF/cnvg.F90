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

      subroutine cnvg(pnew, p, p1, niter, pl)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : norbs, numcal, keywrd, mpack, method_indo
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: niter
      double precision , intent(out) :: pl
      double precision , intent(inout) :: pnew(mpack)
      double precision , intent(inout) :: p(mpack)
      double precision , intent(inout) :: p1(norbs)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, ii, k, i, ie, j
      double precision :: rhfuhf, faca, damp, facb, fac, sum1, a, sa, sum2, sum0, &
        sum
      logical :: extrap

      save rhfuhf, icalcn
!-----------------------------------------------
!***********************************************************************
!
!  CNVG IS A TWO-POINT INTERPOLATION ROUTINE FOR SPEEDING CONVERGENCE
!       OF THE DENSITY MATRIX.
!
! ON OUTPUT P      = NEW DENSITY MATRIX
!           P1     = DIAGONAL OF OLD DENSITY MATRIX
!           PL     = LARGEST DIFFERENCE BETWEEN OLD AND NEW DENSITY
!                    MATRIX DIAGONAL ELEMENTS
!***********************************************************************
      data icalcn/ 0/
      if (icalcn /= numcal) then
        icalcn = numcal
        if (index(keywrd,' UHF') /= 0) then
          rhfuhf = 1.d0
        else
          rhfuhf = 2.d0
        end if
      end if
      pl = 0.0d0
      faca = 0.0d0
      damp = 1.D10
      if (niter > 3) damp = 0.05d0
      if (method_indo) then
        if (niter > 40) damp = 0.01D0
        if (niter > 200) damp = 0.002D0
        if (niter > 350) damp = 0.001D0
      end if
      facb = 0.0d0
      fac = 0.0d0
      ii = mod(niter,3)
      extrap = ii /= 0
      sum1 = 0.d0
      k = 0
      do i = 1, norbs
        k = k + i
        a = pnew(k)
        sum1 = sum1 + a
        sa = Abs (a - p(k))
        if (sa > pl) then
          pl = sa
        end if
        if ( .not. extrap) then
          faca = faca + sa ** 2
          facb = facb + (a - 2.d0*p(k) + p1(i)) ** 2
        end if
        p1(i) = p(k)
        p(k) = a
      end do
      if (facb > 1.d-10 .and. faca < (100.d0*facb)) then
        fac = Sqrt (faca/facb)
      end if
      ie = 0
      sum2 = 0.d0
      do i = 1, norbs
        ii = i - 1
        do j = 1, ii
          ie = ie + 1
          a = pnew(ie)
          p(ie) = a + fac*(a-p(ie))
          pnew(ie) = p(ie)
        end do
        ie = ie + 1
        if (Abs (p(ie) - p1(i)) > damp) then
          p(ie) = p1(i) + Sign (damp, p(ie) - p1(i))
        else
          p(ie) = p(ie) + fac * (p(ie) - p1(i))
        end if
        p(ie) = Min (rhfuhf, Max (p(ie), 0.d0))
        sum2 = sum2 + p(ie)
        pnew(ie) = p(ie)
      end do
   !
   !   RE-NORMALIZE IF ANY DENSITY MATRIX ELEMENTS HAVE BEEN TRUNCATED
   !
      sum0 = sum1
      do
        if (sum2 > 1.d-3) then
          sum = sum1/sum2
        else
          sum = 0.d0
        end if
        sum1 = sum0
        if (sum2 < 1.d-3 .or. Abs (sum-1.d0) < 1.d-5) exit
        sum2 = 0.d0
        do i = 1, norbs
          j = (i*(i + 1))/2
           !
           !   ADD ON A SMALL NUMBER IN CASE AN OCCUPANCY IS EXACTLY ZERO
           !
          p(j) = p(j)*sum + 1.d-20
          p(j) = Max (p(j), 0.d0)
           !
           !  SET UP RENORMALIZATION OVER PARTLY OCCUPIED M.O.'S ONLY.  FULL M.O.'S
           !  CAN'T BE FILLED ANY MORE
           !
          if (p(j) > rhfuhf) then
            p(j) = rhfuhf
            sum1 = sum1 - rhfuhf
          else
            sum2 = sum2 + p(j)
          end if
          pnew(j) = p(j)
        end do
      end do
      end subroutine cnvg
