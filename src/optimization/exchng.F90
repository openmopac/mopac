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

      subroutine exchng(a, b, c, d, x, y, t, q, n)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      double precision , intent(in) :: a
      double precision , intent(out) :: b
      double precision , intent(in) :: c
      double precision , intent(out) :: d
      double precision , intent(in) :: t
      double precision , intent(out) :: q
      double precision , intent(in) :: x(*)
      double precision , intent(out) :: y(*)
!-----------------------------------------------
!********************************************************************
!
! THE CONTENTS OF A, C, T, AND X ARE STORED IN B, D, Q, AND Y!
!
!   THIS IS A DEDICATED ROUTINE, IT IS CALLED BY LINMIN AND LOCMIN ONLY.
!
!********************************************************************
      b = a
      d = c
      q = t
      y(:n) = x(:n)
      return
!
      end subroutine exchng
