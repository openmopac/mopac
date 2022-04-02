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
