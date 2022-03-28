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
