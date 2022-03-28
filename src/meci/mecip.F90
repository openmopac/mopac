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

      subroutine mecip()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use meci_C, only : nmos, lab, occa, nalmat, microa, microb, nelec, &
      & nstate, vectci, deltap
      use common_arrays_C, only : c, p
      use molkst_C, only : norbs

!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, id, jd, ix, iy, ij, k
      double precision :: sum,  delta(norbs, nmos)
!-----------------------------------------------
!***********************************************************************
!
!   MECIP WILL CORRECT THE TOTAL DENSITY MATRIX FOR THE EFFECT OF THE
!   C.I.
!              ON INPUT
!
!  COEFFS       : ALL M.O.'S (NORBS M.O.S)
!  NORBS        : NUMBER OF MOLECULAR ORBITALS = NUMBER OF A.O.'S
!  P            : TOTAL DENSITY MATRIX
!  NMOS         : NUMBER OF M.O.'S IN ACTIVE SPACE
!  VECTCI       : STATE VECTOR OF LENGTH LAB
!  MICROA(I,J)  : ALPHA OCCUPANCY OF M.O. 'I' IN MICROSTATE 'J'
!  MICROB(I,J)  : BETA  OCCUPANCY OF M.O. 'I' IN MICROSTATE 'J'
!
!  NOTE: THIS IS A MODIFICATION OF CODE ORIGINALLY WRITTEN BY
!        PROF. DANIEL LIOTARD
!***********************************************************************
!     INITIALIZE WITH THE OPPOSITE OF THE 'SCF' DENSITY.
      do i = 1, nmos
        deltap(i,i) = -occa(i)*2.D0
        deltap(i,:i-1) = 0.D0
      end do
!
!     ADD THE C.I. CORRECTION
      do id = 1, lab
        do jd = 1, id
!     CHECK SPIN AGREEMENT
          if (nalmat(id) /= nalmat(jd)) cycle
          ix = 0
          iy = 0
          do j = 1, nmos
            ix = ix + abs(microa(j,id)-microa(j,jd))
            iy = iy + abs(microb(j,id)-microb(j,jd))
          end do
!     CHECK NUMBER OF DIFFERING M.O.
          if (ix + iy > 2) cycle
          if (ix == 2) then
!        DETERMINANTS ID AND JD DIFFER BY M.O I IN ID AND M.O J IN JD:
            do i = 1, nmos
              if (microa(i,id) == microa(i,jd)) cycle
              exit
            end do
            ij = microb(i,id)
            do j = i + 1, nmos
              if (microa(j,id) /= microa(j,jd)) exit
              ij = ij + microa(j,id) + microb(j,id)
            end do
!        IJ GIVES THE SIGN OF THE PERMUTATION
            sum = 0.D0
            do k = 1, nstate
              sum = sum + vectci(id+(k-1)*lab)*vectci(jd+(k-1)*lab)
            end do
            deltap(j,i) = deltap(j,i) + sum*dble(1 - 2*mod(ij,2))/nstate
          else if (iy == 2) then
!        DETERMINANTS ID AND JD DIFFER BY M.O J IN ID AND M.O I IN JD:
            do i = 1, nmos
              if (microb(i,id) == microb(i,jd)) cycle
              exit
            end do
            ij = 0
            do j = i + 1, nmos
              if (microb(j,id) /= microb(j,jd)) exit
              ij = ij + microa(j,id) + microb(j,id)
            end do
            ij = ij + microa(j,id)
            sum = 0.D0
            do k = 1, nstate
              sum = sum + vectci(id+(k-1)*lab)*vectci(jd+(k-1)*lab)
            end do
            deltap(j,i) = deltap(j,i) + sum*dble(1 - 2*mod(ij,2))/nstate
          else
!        DETERMINANTS ID AND JD ARE IDENTICAL:
            sum = 0.D0
            do k = 1, nstate
              sum = sum + vectci(id+(k-1)*lab)**2
            end do
            do i = 1, nmos
              deltap(i,i) = deltap(i,i) + (microa(i,id)+microb(i,id))*sum/nstate
            end do
          end if
        end do
      end do
!
!     BACK TRANSFORM INTO A.O. BASIS.
!     -------------------------------
!     P(C.I.) = P(SCF) + C * DELTAP * C'
      do i = 1, nmos
        deltap(:i-1,i) = deltap(i,:i-1)
      end do
!     STEP 1: DELTAP = C * DELTAP
      call mxm (c(1,nelec+1), norbs, deltap, nmos, delta, nmos)
!     STEP 2: P = P + DELTAP * C'
      ij = 0
      do i = 1, norbs
        do j = 1, i
          ij = ij + 1
          sum = 0.D0
          do k = 1, nmos
            sum = sum + delta(i,k)*c(j,nelec+k)
          end do
          p(ij) = p(ij) + sum
        end do
      end do
!     NOTE FROM D.L.: AT THIS POINT THE 'NATURAL ORBITALS' OF THIS STATE
!     CAN BE OBTAINED STRAIGHTWAY AS EIGENVECTORS OF THE DENSITY MATRIX.
      return
      end subroutine mecip
