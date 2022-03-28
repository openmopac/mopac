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

      subroutine mult(c, s, vecs, n)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
      integer , intent(in) :: n
      double precision , intent(in) :: c(n,n)
      double precision , intent(in) :: s(n,n)
      double precision , intent(out) :: vecs(n,n)
!
!  Local variables
!
      integer :: i, j, k
      double precision :: sum
!-----------------------------------------------
!**********************************************************************
!
!   MULT IS USED IN THE MULLIKEN ANALYSIS ONLY. IT PERFORMS THE
!        OPERATION:-
!                                   VECS=BACK-TRANSFORMED EIGENVECTORS
!        VECS  =  C*S               C   =UN-BACK-TRANSFORMED VECTORS
!                                   S   =1/SQRT(OVERLAP MATRIX)
!
!**********************************************************************
      do i = 1, n
        do j = 1, n
          sum = 0.D0
          do k = 1, n
            sum = sum + c(k,i)*s(j,k)
          end do
          vecs(j,i) = sum
        end do
      end do
      return
      end subroutine mult
