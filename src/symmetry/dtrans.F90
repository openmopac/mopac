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

      subroutine dtrans(d, h, ioper, first, r)
!
!  Perform symmetry operation ioper on the "d" orbital set d
!
!  "r" holds the orientation matrix that rotates the system so that it is orientated
!  along the principal axes
!
      USE symmetry_C, only : nclass, elem
      implicit none
      integer , intent(in) :: ioper
      logical , intent(inout) :: first
      double precision , intent(inout) :: d(5)
      double precision , intent(inout) :: h(5)
      double precision , intent(in) :: r(3,3)
      integer :: i, k
      double precision, dimension(5,5,12), save :: t1
      double precision, dimension(3,3) :: s
!-----------------------------------------------
      if (first) then
        first = .FALSE.
!
!  Copy r into s because dtrans modifies the first argument
!
        s = r
        t1 = 0.d0
        call dtran2 (s, t1, 1)
        do k = 2, nclass
          s = elem(:,:,k)
          call dtran2 (s, t1, k)
        end do
      end if
!
!  Matrix multiply d*t1*d(trans)
!
      do i = 1, 5
        h(i) = 0.D0
        h(i) = h(i) + sum(t1(i,:,1)*d)
      end do
      do i = 1, 5
        d(i) = 0.D0
        d(i) = d(i) + sum(t1(:,i,ioper)*h)
      end do
      return
      end subroutine dtrans
