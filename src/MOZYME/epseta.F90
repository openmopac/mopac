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
