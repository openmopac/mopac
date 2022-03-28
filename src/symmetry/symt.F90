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

      subroutine symt(h, deldip, ha)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : nsym, ipo, r
      use molkst_C, only : numat
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
      double precision , intent(inout) :: h(*)
      double precision , intent(inout) :: deldip(3,*)
      double precision , intent(inout) :: ha(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, n, i3, i6, j, k, l, j3, k3, k6, l3, l6, iel33
      double precision, dimension(9) :: temp, temp2
      double precision, dimension(3,3*numat) :: deltmp
!-----------------------------------------------
!****************************************************************
!
!   ON INPUT   H    = HESSIAN MATRIX, PACKED LOWER HALF TRIANGLE
!              R    = SYMMETRY OPERATIONS
!              IPO  = MAP OF ATOMS MOVED
!              NSYM = NUMBER OF SYMMETRY OPERATIONS
!
!   ON OUTPUT  H    = SYMMETRIZED HESSIAN MATRIX
!
!****************************************************************
!  A subroutine that will symmetrize the Hamiltonian, or other matrix
!     by successive application of group operations.  The method used
!     is R H R  added to HA then divided by the total number of symmetry
!     operations used.  This in effects averages all the values in a
!     symmetry correct fashion.
!
!
!
!  Variables used:  (n represents the number of atomic centers)
!     H(3n,3n):  Input/output matrix.  It is a packed lower half triangu
!        matrix.  Commonly, the Hessian.
!     HA(3n,3n): An internal matrix used to sum the symatrized Hessian
!     NSYM:      Input, the value of this symmetry operation.
!     TEMP(9), TEMP2(9):   Temporary matricies used to hold small parts
!          larger matricies for specific matrix operations.
!
!    For the next two items, the last indicy represents the symmetry
!        operation number.
!     IPO(n,*):  A vector that contains the symmetry mapping of atomic c
!
!   Skip this subroutine if NSYMM <= 0.  This implies that only E is pre
      if (nsym < 2) return
!
      ha(:numat*(9*numat+3)/2) = 0.D0
!
      deltmp(1,:numat*3) = 0.D0
      deltmp(2,:numat*3) = 0.D0
      deltmp(3,:numat*3) = 0.D0
!
      do n = 1, nsym
!
!  Now, to actually perform R H R
        do i = 1, numat
          i3 = 3*i
          i6 = 6*i
          do j = 1, i - 1
!
!  Do this multiplication in a 3 by 3 block at a time.  Store H(i,j) in
!    HA( IPO(I,N), IPO(J,N)) or HS( IPO(I,N), IPO(J,N))
!
            k = ipo(i,n)
            l = ipo(j,n)
            j3 = 3*j
            k3 = 3*k
            k6 = 6*k
            l3 = 3*l
            l6 = 6*l
            if (k > l) then
              iel33 = (k3*(k3 - 1))/2 + l3
              temp(9) = h(iel33)
              temp(8) = h(iel33-1)
              temp(7) = h(iel33-2)
              temp(6) = h(iel33-k3+1)
              temp(5) = h(iel33-k3)
              temp(4) = h(iel33-k3-1)
              temp(3) = h(iel33-k6+3)
              temp(2) = h(iel33-k6+2)
              temp(1) = h(iel33-k6+1)
            else
              iel33 = (l3*(l3 - 1))/2 + k3
              temp(9) = h(iel33)
              temp(6) = h(iel33-1)
              temp(3) = h(iel33-2)
              temp(8) = h(iel33-l3+1)
              temp(5) = h(iel33-l3)
              temp(2) = h(iel33-l3-1)
              temp(7) = h(iel33-l6+3)
              temp(4) = h(iel33-l6+2)
              temp(1) = h(iel33-l6+1)
            end if
!
            call mat33 (r(1,n), temp, temp2)
!
            iel33 = (i3*(i3 - 1))/2 + j3
            ha(iel33) = temp2(9) + ha(iel33)
            ha(iel33-1) = temp2(8) + ha(iel33-1)
            ha(iel33-2) = temp2(7) + ha(iel33-2)
            ha(iel33-i*3+1) = temp2(6) + ha(iel33-i3+1)
            ha(iel33-i*3) = temp2(5) + ha(iel33-i3)
            ha(iel33-i*3-1) = temp2(4) + ha(iel33-i3-1)
            ha(iel33-6*i+3) = temp2(3) + ha(iel33-i6+3)
            ha(iel33-6*i+2) = temp2(2) + ha(iel33-i6+2)
            ha(iel33-6*i+1) = temp2(1) + ha(iel33-i6+1)
          end do
          k = ipo(i,n)
          k3 = 3*k
          k6 = 6*k
          iel33 = (k3*(k3 + 1))/2
          temp(9) = h(iel33)
          temp(8) = h(iel33-1)
          temp(7) = h(iel33-2)
          temp(6) = temp(8)
          temp(5) = h(iel33-k3)
          temp(4) = h(iel33-k3-1)
          temp(3) = temp(7)
          temp(2) = temp(4)
          temp(1) = h(iel33-k6+1)
!
          call mat33 (r(1,n), temp, temp2)
!
          iel33 = (i3*(i3 + 1))/2
          ha(iel33) = temp2(9) + ha(iel33)
          ha(iel33-1) = temp2(8) + ha(iel33-1)
          ha(iel33-2) = temp2(7) + ha(iel33-2)
          ha(iel33-i*3) = temp2(5) + ha(iel33-i3)
          ha(iel33-i*3-1) = temp2(4) + ha(iel33-i3-1)
          ha(iel33-6*i+1) = temp2(1) + ha(iel33-i6+1)
!
!  APPLY SYMMETRY TO DIPOLE TERM AS WELL
!
          temp(9) = deldip(3,k*3)
          temp(8) = deldip(2,k*3)
          temp(7) = deldip(1,k*3)
          temp(6) = deldip(3,k*3-1)
          temp(5) = deldip(2,k*3-1)
          temp(4) = deldip(1,k*3-1)
          temp(3) = deldip(3,k*3-2)
          temp(2) = deldip(2,k*3-2)
          temp(1) = deldip(1,k*3-2)
!
          call mat33 (r(1,n), temp, temp2)
!
          deltmp(3,i*3) = temp2(9) + deltmp(3,i*3)
          deltmp(2,i*3) = temp2(8) + deltmp(2,i*3)
          deltmp(1,i*3) = temp2(7) + deltmp(1,i*3)
          deltmp(3,i*3-1) = temp2(6) + deltmp(3,i*3-1)
          deltmp(2,i*3-1) = temp2(5) + deltmp(2,i*3-1)
          deltmp(1,i*3-1) = temp2(4) + deltmp(1,i*3-1)
          deltmp(3,i*3-2) = temp2(3) + deltmp(3,i*3-2)
          deltmp(2,i*3-2) = temp2(2) + deltmp(2,i*3-2)
          deltmp(1,i*3-2) = temp2(1) + deltmp(1,i*3-2)
!
        end do
      end do
!
      h(:numat*(9*numat+3)/2) = ha(:numat*(9*numat+3)/2)/nsym
!
      deldip(1,:3*numat) = deltmp(1,:3*numat)/nsym
      deldip(2,:3*numat) = deltmp(2,:3*numat)/nsym
      deldip(3,:3*numat) = deltmp(3,:3*numat)/nsym
!
      return
      end subroutine symt
