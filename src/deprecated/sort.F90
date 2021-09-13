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

      subroutine sort(val, vec, n) 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real , intent(inout) :: val(*) 
      complex , intent(inout) :: vec(n,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k 
      real :: x 
      complex :: sum 
!-----------------------------------------------
      do i = 1, n 
        k = 0
        x = 1.E9 
        do j = i, n 
          if (val(j) >= x) cycle  
          k = j 
          x = val(j) 
        end do 
        do j = 1, n 
          sum = vec(j,k) 
          vec(j,k) = vec(j,i) 
          vec(j,i) = sum 
        end do 
        val(k) = val(i) 
        val(i) = x 
      end do 
      return  
      end subroutine sort 
