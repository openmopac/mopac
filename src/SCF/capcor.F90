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

      double precision function capcor (nat, nfirst, nlast, p, h)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : numat
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nat(numat)
      integer , intent(in) :: nfirst(numat)
      integer , intent(in) :: nlast(numat)
      double precision , intent(in) :: p(*)
      double precision , intent(in) :: h(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ni, il, iu, j, ii, k, jl, kk
      double precision :: sum
!-----------------------------------------------
!*****************************************************************
!
!    CORRECTION TO ELECTRONIC ENERGY DUE TO CAPPED BONDS
!
!*****************************************************************
      sum = 0.D0
      do i = 1, numat
        ni = nat(i)
        il = nfirst(i)
        iu = nlast(i)
        if (ni == 102) then
!
!   DO ENTIRE ROW - NO NEED TO CHECK FURTHER.
!
          j = (nlast(i)*(nlast(i)+1))/2
          ii = iu - 1
          do k = 1, ii
            j = j - 1
            sum = sum + p(j)*h(j)
          end do
        else
          do j = 1, i
            jl = nfirst(j)
            if (nat(j) /= 102) cycle
            do k = il, iu
              kk = (k*(k - 1))/2 + jl
              sum = sum + p(kk)*h(kk)
            end do
          end do
        end if
      end do
!
!   DOUBLE SUM SINCE WE ONLY CALCULATED LOWER HALF, AND CAPCOR
!   WILL APPEAR IN 1/2*P(H+F).  ONLY H PART OF F WILL BE USED.
      capcor = -sum*2.D0
      return
      end function capcor
