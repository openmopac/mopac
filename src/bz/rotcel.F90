! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
