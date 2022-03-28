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

      subroutine dijkl2(dc)
      use meci_C, only : dijkl, xy, nmos
      use molkst_C, only : norbs
      implicit none
      double precision  :: dc(norbs,nmos)
!
      integer :: ij, i, j, kl, k, ll, l
      double precision :: val, val2
      logical :: lij, lkl
      double precision, external :: ddot
!***********************************************************************
!     RELAXATION OF 2-ELECTRONS INTEGRALS IN M.O BASIS.
!
!   INPUT
!   DC(NORBS,NMOS) : C.I-ACTIVE M.O DERIVATIVES IN M.O BASIS, IN COLUMN.
!   NORBS          : TOTAL NUMBER OF M.O.
!   NMOS           : NUMBER OF C.I-ACTIVE M.O.
!   DIJKL(I,J,KL)  : <I(1),J(1)|K(2),L(2)> WITH
!                     I              OVER     ALL    M.O.
!                     J,KL CANONICAL OVER C.I-ACTIVE M.O.
!   OUTPUT
!     xy(I,J,K,L)= d< I(1),J(1) | K(2),L(2) >
!                   = <dI,J|K,L> + <I,dJ|K,L> + <I,J|dK,L> + <I,J|K,dL>
!                     WITH I,J,K,L OVER ALL C.I-ACTIVE M.O.
!     WRITTEN BY DANIEL LIOTARD
! (NOTE BY JJPS: AS THIS CODE IS HIGHLY EFFICIENT, NO CHANGES WERE MADE)
!***********************************************************************
!
      ij = 0
      do i = 1, nmos
        do j = 1, i
          ij = ij + 1
          lij = i == j
          kl = 0
          do k = 1, i
            if (k == i) then
              ll = j
            else
              ll = k
            end if
            do l = 1, ll
              kl = kl + 1
              lkl = k == l
              val = ddot(norbs,dc(1,i),1,dijkl(1,j,kl),1)
              if (lij .and. lkl .and. j==k) then
                val = val*4.D0
              else
                if (lij) then
                  val = val*2.D0
                else
                  val = val + ddot(norbs,dc(1,j),1,dijkl(1,i,kl),1)
                end if
                val2 = ddot(norbs,dc(1,k),1,dijkl(1,l,ij),1)
                if (lkl) then
                  val = val + val2*2.D0
                else
                  val = val + val2 + ddot(norbs,dc(1,l),1,dijkl(1,k,ij),1)
                end if
              end if
              xy(i,j,k,l) = val
              xy(i,j,l,k) = val
              xy(j,i,k,l) = val
              xy(j,i,l,k) = val
              xy(k,l,i,j) = val
              xy(k,l,j,i) = val
              xy(l,k,i,j) = val
              xy(l,k,j,i) = val
            end do
          end do
        end do
      end do
      return
      end subroutine dijkl2
