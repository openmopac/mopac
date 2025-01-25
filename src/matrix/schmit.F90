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

      subroutine schmit(u, n, ndim)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: ndim
      double precision , intent(inout) :: u(ndim,ndim)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, k, k1, npass, j
      double precision :: zero, small, one, dot, scale

      save zero, small, one
!-----------------------------------------------
      data zero, small, one/ 0.D0, 0.01D0, 1.D0/
      ii = 0
      do k = 1, n
        k1 = k - 1
!
!     NORMALIZE KTH COLUMN VECTOR
!
        dot = zero
        dot = dot + sum(u(:n,k)*u(:n,k))
        if (Abs(dot) < 1.d-20) go to 100
        scale = one/sqrt(dot)
        u(:n,k) = scale*u(:n,k)
   30   continue
        if (k1 == 0) cycle
        npass = 0
!
!     PROJECT OUT K-1 PREVIOUS ORTHONORMAL VECTORS FROM KTH VECTOR
!
   40   continue
        npass = npass + 1
        do j = 1, k1
          dot = zero
          dot = dot + sum(u(:n,j)*u(:n,k))
          u(:n,k) = u(:n,k) - dot*u(:n,j)
        end do
!
!     SECOND NORMALIZATION (AFTER PROJECTION)
!     IF KTH VECTOR IS SMALL BUT NOT ZERO THEN NORMALIZE
!     AND PROJECT AGAIN TO CONTROL ROUND-OFF ERRORS.
!
        dot = zero
        dot = dot + sum(u(:n,k)*u(:n,k))
        if (Abs(dot) < 1.d-20) go to 100
        if (dot<small .and. npass>2) go to 100
        scale = one/sqrt(dot)
        u(:n,k) = scale*u(:n,k)
        if (dot < small) go to 40
        cycle
!
!     REPLACE LINEARLY DEPENDENT KTH VECTOR BY A UNIT VECTOR.
!
  100   continue
        ii = ii + 1
!     IF(II.GT.N) CALL MOPEND
        u(ii,k) = one
        go to 30
      end do
      return
      end subroutine schmit
