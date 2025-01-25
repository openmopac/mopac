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
