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

      subroutine sympop(h, i, iskip, deldip)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use symmetry_C, only : nsym, ipo
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i
      integer , intent(out) :: iskip
      double precision  :: h(*)
      double precision  :: deldip(3,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j
!-----------------------------------------------
      do j = 1, nsym
        if (ipo(i,j) >= i) cycle
        call symh (h, deldip, i, j, ipo)
        iskip = 3
!  atom ipo(i,j) is suitable for transition dipole calc'n
!
!            K=I*3-2
!            WRITE(IW,*)' Transition dipole after symmetry operation'
!            WRITE(IW,'(3f12.5)')(deldip(l,k),l=1,3)
!            WRITE(IW,'(3f12.5)')(deldip(l,k+1),l=1,3)
!            WRITE(IW,'(3f12.5)')(deldip(l,k+2),l=1,3)
!
!   INSERT DELDIP ROTATION HERE
!
        go to 20
      end do
      iskip = 0
   20 continue
      return
      end subroutine sympop
