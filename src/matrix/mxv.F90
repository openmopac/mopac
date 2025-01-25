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

      subroutine mxv(a, nar, vecx, nbr, vecy)
      implicit none
      integer, parameter :: incy = 1, incx = 1
      integer  :: nar
      integer  :: nbr
      double precision  :: a(nar,nbr)
      double precision  :: vecx(nbr)
      double precision  :: vecy(nar)
!
!     RECTANGULAR MATRIX-VECTOR PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.

      call dgemv ('N', nar, nbr, 1.0d0, a, nar, vecx, incx, 0.0d0, vecy, incy)
      return
      end subroutine mxv
