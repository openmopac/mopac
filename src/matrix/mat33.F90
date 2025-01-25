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

      subroutine mat33(a, b, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!  A subroutine that will multiply two 3 by 3 matricies in the following
!     fashion:    C = A(transpose) B A
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(9)
      double precision , intent(in) :: b(9)
      double precision , intent(out) :: c(9)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision, dimension(9) :: t
!-----------------------------------------------
!
!
      t(1) = b(1)*a(1) + b(2)*a(4) + b(3)*a(7)
      t(2) = b(1)*a(2) + b(2)*a(5) + b(3)*a(8)
      t(3) = b(1)*a(3) + b(2)*a(6) + b(3)*a(9)
      t(4) = b(4)*a(1) + b(5)*a(4) + b(6)*a(7)
      t(5) = b(4)*a(2) + b(5)*a(5) + b(6)*a(8)
      t(6) = b(4)*a(3) + b(5)*a(6) + b(6)*a(9)
      t(7) = b(7)*a(1) + b(8)*a(4) + b(9)*a(7)
      t(8) = b(7)*a(2) + b(8)*a(5) + b(9)*a(8)
      t(9) = b(7)*a(3) + b(8)*a(6) + b(9)*a(9)
!
      c(1) = a(1)*t(1) + a(4)*t(4) + a(7)*t(7)
      c(2) = a(1)*t(2) + a(4)*t(5) + a(7)*t(8)
      c(3) = a(1)*t(3) + a(4)*t(6) + a(7)*t(9)
      c(4) = a(2)*t(1) + a(5)*t(4) + a(8)*t(7)
      c(5) = a(2)*t(2) + a(5)*t(5) + a(8)*t(8)
      c(6) = a(2)*t(3) + a(5)*t(6) + a(8)*t(9)
      c(7) = a(3)*t(1) + a(6)*t(4) + a(9)*t(7)
      c(8) = a(3)*t(2) + a(6)*t(5) + a(9)*t(8)
      c(9) = a(3)*t(3) + a(6)*t(6) + a(9)*t(9)
!
      return
      end subroutine mat33
