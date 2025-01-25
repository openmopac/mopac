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

subroutine epseta (eps, eta)
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    double precision, intent (out) :: eps, eta
   !
   !.. Intrinsic Functions ..
    intrinsic Max, Tiny, Epsilon

   ! ... Executable Statements ...
   !
   !     RETURN ETA, THE SMALLEST REPRESENTABLE NUMBER,
   !     AND EPS, THE SMALLEST NUMBER FOR WHICH 1+EPS.NE.1.
   !
    eps = Max (Tiny (1.d0), 1.d-39)
    eta = Max (Epsilon (1.d0), 1.d-17)
end subroutine epseta
