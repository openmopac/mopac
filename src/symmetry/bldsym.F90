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

      subroutine bldsym(ioper, iplace)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE symmetry_C, only : elem, cub
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
      integer , intent(in) :: ioper
      integer  :: iplace
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(3,20) :: j
      integer :: i
      double precision :: twopi, angle

      save j, twopi
!-----------------------------------------------
!***********************************************************************
!
!   BLDSYM  constructs the point-group symmetry operations as
!           3 by 3 matrices, 20 in all.  The operations are, in
!           order
!   1 C2(X)         6 Sigma(YZ)     11 C6      16 S8
!   2 C2(Y)         7 inversion     12 C7      17 S10
!   3 C2(Z)         8 C3            13 C8      18 S12
!   4 Sigma(XY)     9 C4            14 S4      19 S(5) ?
!   5 Sigma(XZ)    10 C5            15 S6      20 C(infinity)
!
!***********************************************************************
      data j(1,1), j(2,1), j(3,1)/ 1, -1, -1/
      data j(1,2), j(2,2), j(3,2)/ -1, 1, -1/
      data j(1,3), j(2,3), j(3,3)/ -1, -1, 1/
      data j(1,4), j(2,4), j(3,4)/ 1, 1, -1/
      data j(1,5), j(2,5), j(3,5)/ 1, -1, 1/
      data j(1,6), j(2,6), j(3,6)/ -1, 1, 1/
      data j(1,7), j(2,7), j(3,7)/ -1, -1, -1/
      data j(1,8), j(2,8), j(3,8)/ 3, 0, 1/
      data j(1,9), j(2,9), j(3,9)/ 4, 0, 1/
      data j(1,10), j(2,10), j(3,10)/ 5, 0, 1/
      data j(1,11), j(2,11), j(3,11)/ 6, 0, 1/
      data j(1,12), j(2,12), j(3,12)/ 7, 0, 1/
      data j(1,13), j(2,13), j(3,13)/ 8, 0, 1/
      data j(1,14), j(2,14), j(3,14)/ 4, 0, -1/
      data j(1,15), j(2,15), j(3,15)/ 6, 0, -1/
      data j(1,16), j(2,16), j(3,16)/ 8, 0, -1/
      data j(1,17), j(2,17), j(3,17)/ 10, 0, -1/
      data j(1,18), j(2,18), j(3,18)/ 12, 0, -1/
      data j(1,19), j(2,19), j(3,19)/ 5, 0, -1/
      data j(1,20), j(2,20), j(3,20)/ 0, 0, -1/
      data twopi/ 6.2831853071796D0/
      do i = 1, 3
        elem(i,:,iplace) = 0.D0
        elem(i,i,iplace) = j(i,ioper)
      end do
      if (ioper /= 20) then
        if (j(1,ioper) >= 2) then
          angle = twopi/dble(j(1,ioper))
          elem(1,1,iplace) = cos(angle)
          elem(2,2,iplace) = elem(1,1,iplace)
          elem(2,1,iplace) = sin(angle)
          elem(1,2,iplace) = -elem(2,1,iplace)
        end if
        if (ioper==8 .or. ioper==15) call mult33 (cub, iplace)
        return
      end if
      elem(1,2,iplace) = 1.D0
      elem(2,1,iplace) = 1.D0
      return
      end subroutine bldsym
