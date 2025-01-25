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

subroutine sort (val, vec, n)
  implicit none
  integer, intent(in) :: n
  real, dimension(*), intent(inout) :: val
  complex, dimension(n, *), intent(inout) :: vec
  integer :: i, j, k
  real :: x
  complex :: sum
  ! 
  ! ... Executable Statements ...
  ! 
  do i = 1, n
    x = 1.e9
    do j = i, n
      if (val(j) < x) then
        k = j
        x = val(j)
      end if
    end do
    do j = 1, n
      sum = vec(j, k)
      vec(j, k) = vec(j, i)
      vec(j, i) = sum
    end do
    val(k) = val(i)
    val(i) = x
  end do
end subroutine sort
