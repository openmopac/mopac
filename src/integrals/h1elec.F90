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

      subroutine h1elec(ni, nj, xi, xj, smat)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use parameters_C, only : natorb, betas, betap, betad
      use MOZYME_C, only : cutofs
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ni
      integer  :: nj
      double precision , intent(in) :: xi(3)
      double precision , intent(in) :: xj(3)
      double precision  :: smat(9,9)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, norbi, norbj
      double precision :: bi(9), bj(9), xjuc(3), rab
!-----------------------------------------------
!***********************************************************************
!
!  H1ELEC FORMS THE ONE-ELECTRON MATRIX BETWEEN TWO ATOMS.
!
!   ON INPUT    NI   = ATOMIC NO. OF FIRST ATOM.
!               NJ   = ATOMIC NO. OF SECOND ATOM.
!               XI   = COORDINATES OF FIRST ATOM.
!               XJ   = COORDINATES OF SECOND ATOM.
!
!   ON OUTPUT   SMAT = MATRIX OF ONE-ELECTRON INTERACTIONS.
!
!***********************************************************************
      rab = (xi(1)-xj(1))**2+(xi(2)-xj(2))**2+(xi(3)-xj(3))**2
      if (rab > cutofs .or. (rab > 3.24D0 .and. (ni==102 .or. nj==102))) then
        smat = 0.D0
        return
      end if
      xjuc = xj - xi
      call diat (ni, nj, xjuc, smat)
      bi(1) = betas(ni)*0.5D0
      bi(2) = betap(ni)*0.5D0
      bi(3) = bi(2)
      bi(4) = bi(2)
      bi(5) = betad(ni)*0.5D0
      bi(6) = bi(5)
      bi(7) = bi(5)
      bi(8) = bi(5)
      bi(9) = bi(5)
      bj(1) = betas(nj)*0.5D0
      bj(2) = betap(nj)*0.5D0
      bj(3) = bj(2)
      bj(4) = bj(2)
      bj(5) = betad(nj)*0.5D0
      bj(6) = bj(5)
      bj(7) = bj(5)
      bj(8) = bj(5)
      bj(9) = bj(5)
      norbi = natorb(ni)
      norbj = natorb(nj)
      do j = 1, norbj
        smat(:norbi,j) = smat(:norbi,j)*(bi(:norbi)+bj(j))
      end do
      return
      end subroutine h1elec
