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

      subroutine mecih(diag, cimat, nmos, lab, xy)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use meci_C, only : microa, microb, nalmat, ispqr, is, iiloop, jloop
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: nmos
      integer , intent(in) :: lab
      double precision, dimension(*), intent(in)  :: xy
      double precision , intent(in) :: diag(*)
      double precision , intent(out) :: cimat(*)
!-----------------------------------------------
!   L o c a l   V a r iiloop a b l e s
!-----------------------------------------------
      integer :: ik, ix, iy, j
      double precision, external :: aababc, aabacd, aabbcd, babbbc, babbcd
!-----------------------------------------------
!
!     BUILD THE C.I. MATRIX 'CIMAT' IN PACKED CANONICAL FORM.
!
!
      ik = 0
!
!     OUTER LOOP TO FILL C.I. MATRIX.
      do iiloop = 1, lab
        is = 2
!
!     INNER LOOP.
        do jloop = 1, iiloop
          ik = ik + 1
          cimat(ik) = 0.D0
          ix = 0
          iy = 0
          do j = 1, nmos
            ix = ix + abs(microa(j,iiloop)-microa(j,jloop))
            iy = iy + abs(microb(j,iiloop)-microb(j,jloop))
          end do
!
!                              CHECK IF MATRIX ELEMENT HAS TO BE ZERO
!
          if (ix + iy>4 .or. nalmat(iiloop)/=nalmat(jloop)) cycle
          if (ix + iy == 4) then
            if (ix == 0) then
              cimat(ik) = babbcd(microa(1,iiloop),microb(1,iiloop),microa(1,jloop),microb(1,jloop&
                ),nmos,xy)
            else if (ix == 2) then
              cimat(ik) = aabbcd(microa(1,iiloop),microb(1,iiloop),microa(1,jloop),microb(1,jloop&
                ),nmos,xy)
            else
              cimat(ik) = aabacd(microa(1,iiloop),microb(1,iiloop),microa(1,jloop),microb(1,jloop&
                ),nmos,xy)
            end if
          else if (ix == 2) then
            cimat(ik) = aababc(microa(1,iiloop),microb(1,iiloop),microa(1,jloop),nmos,xy)
          else if (iy == 2) then
            cimat(ik) = babbbc(microa(1,iiloop),microb(1,iiloop),microb(1,jloop),nmos,xy)
          else
            cimat(ik) = diag(iiloop)
          end if
        end do
        ispqr(iiloop,1) = is - 1
      end do
      return
      end subroutine mecih
