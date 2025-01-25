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

      subroutine mxm(a, nar, b, nbr, c, ncc)
      implicit none
      integer  :: nar
      integer  :: nbr
      integer  :: ncc
      double precision  :: a(nar,nbr)
      double precision  :: b(nbr,ncc)
      double precision  :: c(nar,ncc)
!
!     RECTANGULAR MATRIX PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.
!
      call dgemm ('N', 'N', nar, ncc, nbr, 1.0D0, a, nar, b, nbr, 0.0D0, c, nar)
      return
      end subroutine mxm
