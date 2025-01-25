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

subroutine chrge_for_MOZYME (p, q)
!
! Evaluate the atomic partial charge.
!
    use molkst_C, only: mpack, numat
    use MOZYME_C, only : iorbs
    implicit none
    double precision, dimension (mpack), intent (in) :: p
    double precision, dimension (numat), intent (out) :: q
!
    integer :: i, ii, j
    double precision :: sum
    integer, external :: ijbo
!
    do i = 1, numat
      ii = ijbo (i, i)
      sum = 0.d0
      do j = 1, iorbs(i)
        ii = ii + j
        sum = sum + p(ii)
      end do
      q(i) = sum
    end do
end subroutine chrge_for_MOZYME
