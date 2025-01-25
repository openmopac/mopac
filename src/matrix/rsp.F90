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

      subroutine rsp(a, n, root, vect)
      implicit none
      integer  :: n
      double precision  :: a(*)
      double precision, intent (out)  :: root(n)
      double precision  :: vect(n,n)
!
! Trivial case: n = 1
!
      if (n == 1) then
        root(1) = a(1)
        vect(1,1) = 1.d0
        return
      end if
      call eigenvectors_LAPACK(vect, a, root, n)
      return
      end subroutine rsp
