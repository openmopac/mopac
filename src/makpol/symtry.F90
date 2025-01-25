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
