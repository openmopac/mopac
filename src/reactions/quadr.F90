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

      subroutine quadr(f0, f1, f2, x1, x2, a, b, c)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: f0
      double precision , intent(in) :: f1
      double precision , intent(in) :: f2
      double precision , intent(in) :: x1
      double precision , intent(in) :: x2
      double precision , intent(out) :: a
      double precision , intent(out) :: b
      double precision , intent(out) :: c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!***************************************************
!                                                  *
!    QUADR CALCULATES THE A, B AND C IN THE EQUNS. *
!                                                  *
!     A                   =   F0                   *
!     A + B.X0 + C.X0**2  =   F1                   *
!     A + B.X2 + C.X2**2  =   F2                   *
!                                                  *
!***************************************************
      c = (x2*(f1 - f0) - x1*(f2 - f0))/(x2*x1*(x1 - x2))
      b = (f1 - f0 - c*x1**2)/x1
      a = f0
      return
      end subroutine quadr
