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

subroutine symtnn (geo, na)
!**********************************************************************
!
!  SYMTRY COMPUTES THE BOND LENGTHS AND ANGLES THAT ARE FUNCTIONS OF
!         OTHER BOND LENGTHS AND ANGLES.
!
! ON INPUT GEO     = KNOWN INTERNAL COORDINATES
!          NDEP    = NUMBER OF DEPENDENCY FUNCTIONS.
!          IDEPFN  = ARRAY OF DEPENDENCY FUNCTIONS.
!          LOCDEP  = ARRAY OF LABELS OF DEPENDENT ATOMS.
!          LOCPAR  = ARRAY OF LABELS OF REFERENCE ATOMS.
!
!  ON OUTPUT THE ARRAY "GEO" IS FILLED
!***********************************************************************
    use common_systm
    use common_symult
    use common_geosym
    implicit none
    double precision, dimension (3, natoms), intent (inout) :: geo
    integer, dimension (natoms), intent (in) :: na
    integer :: i, j, locn, n
    double precision :: value
!
!     NOW COMPUTE THE DEPENDENT PARAMETERS.
!
    n = 0
    do i = 1, ndep
      if (idepfn(i) == 19 .and. depmul(n+1) > 1.d-3) then
        n = n + 1
        call haddon (value, locn, idepfn(i), locpar(i), geo, na, depmul(n))
      else
        call haddon (value, locn, idepfn(i), locpar(i), geo, na, depmul(1))
      end if
      j = locdep(i)
      if(j <= natoms) geo(locn, j) = value
    end do
end subroutine symtnn
