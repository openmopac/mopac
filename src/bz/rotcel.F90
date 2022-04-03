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

subroutine rotcel (t2, isize, jsize, t1)
  implicit none
  integer, intent(in) :: isize, jsize
  double precision, dimension(9,9), intent(inout) :: t2
  double precision, dimension(9,9), intent(in) :: t1
  integer :: i, j, k, k2, l, l2
  double precision :: sum
  double precision, dimension(isize,jsize) :: work 
  do i = 1, jsize
    do j = 1, isize
      sum = 0.d0
      k2 = 0
      do k = 1, isize
        k2 = k2 + 1
        l2 = 0
        do l = 1, jsize
          l2 = l2 + 1
          sum = sum + t1(l2, i)*t2(k, l)*t1(k2, j)
        end do
      end do
      work(j, i) = sum
    end do
  end do
  do i = 1, isize
    do j = 1, jsize
      t2(i, j) = work(i, j)
    end do
  end do
end subroutine rotcel
