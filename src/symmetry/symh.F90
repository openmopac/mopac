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

      subroutine symh(h, dip, i, n, ipo)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : r
      USE molkst_C, ONLY: numat
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i
      integer , intent(in) :: n
      integer, dimension (numat, 120), intent (in) :: ipo
      double precision , intent(inout) :: h(*)
      double precision , intent(inout) :: dip(3,*)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i3, i6, k3, k6, j, l, l3, l6, j3, j6, iel33, im1t3, istart
      double precision, dimension(9) :: temp, temp2
      double precision :: fact
!-----------------------------------------------
!****************************************************************
!
!  INPUT:   H()   A packed lower triangular hessian
!           DIP(,) A MATRIX OF DIPOLE TENSORS TO BE SYMM
!           R(,)  A matrix of symmetry operations
!           IPO(,) A matrix of atomic mapping according to R
!           I     The atom (row and column) to add to H()
!           N     The symmetry operation to use to generate I
!
!  OUTPUT:  H()   A packed lower triangular Hessian with information
!                   about atom I added
!           DIP(,) A MATRIX OF DIPOLE TENSORS THAT HAVE BEEN SYMM
!
!****************************************************************
!
!
!  This subroutine will add all necessary information to the Hessian
!    concerning atom I.  Since the Hessian is a packed lower half
!    triangle, the existing information for atom pair (K,L) where K,L
!    < I is fully known, (K > I and L < I) or (vice versa) is half
!    known, K,L > I is completely unknown.
!    Therefore, start in unknown region and make it half known.  Double
!    known values, and move in the diagonal element at full strength.
!
!
!
!
!
!  Variables used:  (n represents the number of atomic centers)
!     H(3n,3n):  Input/output matrix.  It is a packed lower half triangu
!        matrix.  Commonly, the Hessian.
!     TEMP(9), TEMP2(9):   Temporary matricies used to hold small parts
!          larger matricies for specific matrix operations.
!
!    For the next two items, the last indicy represents the symmetry
!        operation number.
!     R(14,*):   The first 9 elements of each record are a packed 3 by 3
!          array of a given symmetry operations.  Elements 10 - 14 are t
!          users input describing the symmetry operation.
!     IPO(n,*):  A vector that contains the symmetry mapping of atomic c
!
!
      k = ipo(i,n)
      i3 = 3*i
      i6 = 6*i
      k3 = 3*k
      k6 = 6*k
!
!  Now, to climb up the matrix
      do j = numat, i + 1, -1
        l = ipo(j,n)
        l3 = 3*l
        l6 = 6*l
        j3 = 3*j
        j6 = 6*j
!
!  Now, to actually perform R H R
!
!  Do this multiplication in a 3 by 3 block at a time.  Store H(i,j) in
!    H( IPO(I,N), IPO(J,N))
!
        if (k > l) then
          iel33 = (k3*(k3 - 1))/2 + l3
          temp(9) = 0.5D0*h(iel33)
          temp(8) = 0.5D0*h(iel33-1)
          temp(7) = 0.5D0*h(iel33-2)
          temp(6) = 0.5D0*h(iel33-k3+1)
          temp(5) = 0.5D0*h(iel33-k3)
          temp(4) = 0.5D0*h(iel33-k3-1)
          temp(3) = 0.5D0*h(iel33-k6+3)
          temp(2) = 0.5D0*h(iel33-k6+2)
          temp(1) = 0.5D0*h(iel33-k6+1)
        else
          iel33 = (l3*(l3 - 1))/2 + k3
          fact = 1.0D0
          if (l < i) fact = 0.5D0
          temp(9) = fact*h(iel33)
          temp(6) = fact*h(iel33-1)
          temp(3) = fact*h(iel33-2)
          temp(8) = fact*h(iel33-l3+1)
          temp(5) = fact*h(iel33-l3)
          temp(2) = fact*h(iel33-l3-1)
          temp(7) = fact*h(iel33-l6+3)
          temp(4) = fact*h(iel33-l6+2)
          temp(1) = fact*h(iel33-l6+1)
        end if
!
        call mat33 (r(1,n), temp, temp2)
!
        iel33 = (j3*(j3 - 1))/2 + i3
        h(iel33) = temp2(9)
        h(iel33-j3+1) = temp2(8)
        h(iel33-j6+3) = temp2(7)
        h(iel33-1) = temp2(6)
        h(iel33-j3) = temp2(5)
        h(iel33-j6+2) = temp2(4)
        h(iel33-2) = temp2(3)
        h(iel33-j3-1) = temp2(2)
        h(iel33-j6+1) = temp2(1)
      end do
!
!  Now, to do the diagonal term
!
      iel33 = (k3*(k3 + 1))/2
      temp(9) = 0.5D0*h(iel33)
      temp(8) = 0.5D0*h(iel33-1)
      temp(7) = 0.5D0*h(iel33-2)
      temp(6) = temp(8)
      temp(5) = 0.5D0*h(iel33-k3)
      temp(4) = 0.5D0*h(iel33-k3-1)
      temp(3) = temp(7)
      temp(2) = temp(4)
      temp(1) = 0.5D0*h(iel33-k6+1)
!
      call mat33 (r(1,n), temp, temp2)
!
      iel33 = (i3*(i3 + 1))/2
      h(iel33) = temp2(9)
      h(iel33-1) = temp2(8)
      h(iel33-2) = temp2(7)
      h(iel33-i3) = temp2(5)
      h(iel33-i3-1) = temp2(4)
      h(iel33-i6+1) = temp2(1)
!
!   NOW, TO ROTATE THE DIPOLE TENSOR TERM
!
      temp(9) = dip(3,k*3)
      temp(8) = dip(2,k*3)
      temp(7) = dip(1,k*3)
      temp(6) = dip(3,k*3-1)
      temp(5) = dip(2,k*3-1)
      temp(4) = dip(1,k*3-1)
      temp(3) = dip(3,k*3-2)
      temp(2) = dip(2,k*3-2)
      temp(1) = dip(1,k*3-2)
!
      call mat33 (r(1,n), temp, temp2)
!
      dip(3,i*3) = temp2(9)
      dip(2,i*3) = temp2(8)
      dip(1,i*3) = temp2(7)
      dip(3,i*3-1) = temp2(6)
      dip(2,i*3-1) = temp2(5)
      dip(1,i*3-1) = temp2(4)
      dip(3,i*3-2) = temp2(3)
      dip(2,i*3-2) = temp2(2)
      dip(1,i*3-2) = temp2(1)
!
!   Now, to double all existing values going across
!
      im1t3 = (i - 1)*3
      istart = (im1t3*(im1t3 + 1))/2 + 1
      h(istart:iel33) = h(istart:iel33) + h(istart:iel33)
!  Everything is now done for this symmetry element.
!
      return
      end subroutine symh
