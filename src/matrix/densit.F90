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

      subroutine densit(c, mdim, norbs, nocc, occ, nfract, fract, p, mode)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mdim
      integer , intent(in) :: norbs
      integer , intent(in) :: nocc
      integer , intent(in) :: nfract
      integer , intent(in) :: mode
      double precision , intent(in) :: occ, fract
      double precision , intent(in) :: c(mdim,mdim)
      double precision , intent(inout) :: p((norbs*(norbs + 1))/2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: norbs2, nl2, nu2, nl1, nu1, l, i, j
      double precision :: sign, frac, const, sum2, sum1
!-----------------------------------------------
!***********************************************************************
!
!   DENSIT COMPUTES THE DENSITY MATRIX GIVEN THE EIGENVECTOR MATRIX, AND
!          INFORMATION ABOUT THE M.O. OCCUPANCY.
!
!  INPUT:  C     = SQUARE EIGENVECTOR MATRIX, C IS OF SIZE MDIM BY MDIM
!                  AND THE EIGENVECTORS ARE STORED IN THE TOP LEFT-HAND
!                  CORNER.
!          NORBS = NUMBER OF ORBITALS
!          NOCC  = NUMBER OF FULLY-OCCUPIED M.O.S
!          occ   = NUMBER OF FRACTIONALLY OCCUPIED M.O.S.
!          MODE  = 2 IF POSITRON EQUIVALENT IS NOT TO BE USED
!
!   ON EXIT: P   = DENSITY MATRIX
!
!***********************************************************************
!
! SET UP LIMITS FOR SUMS
!  NL1 = BEGINING OF ONE ELECTRON SUM
!  NU1 = END OF SAME
!  NL2 = BEGINING OF TWO ELECTRON SUM
!  NU2 = END OF SAME
!
      norbs2 = norbs/2
      if (nocc /= 0 .and. nfract > norbs2 .and. mode /= 2) then
!
!    TAKE POSITRON EQUIVALENT
!
        sign = -1.d0
        frac = occ - fract
        const = occ
        nl2 = nfract + 1
        nu2 = norbs
        nl1 = nocc + 1
        nu1 = nfract
      else
!
!    TAKE ELECTRON EQUIVALENT
!
        sign = 1.d0
        frac = fract
        const = 0.d0
        nl2 = 1
        nu2 = nocc
        nl1 = nocc + 1
        nu1 = nfract
      end if
      l = 0
      do i = 1, norbs
        do j = 1, i
          l = l + 1
          sum1 = 0.D0
          sum2 = sum(c(i,nl2:nu2)*c(j,nl2:nu2))
          sum2 = sum2*occ
          sum1 = sum(c(i,nl1:nu1)*c(j,nl1:nu1))
          p(l) = (sum2 + sum1*frac)*sign
        end do
        p(l) = const + p(l)
      end do
      return
      end subroutine densit
