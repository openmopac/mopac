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

subroutine maksym (loc, xparam)
    use common_systm
    use common_geosym, only : locpar, idepfn, locdep
    implicit none
    integer, dimension (2, nvar), intent (inout) :: loc
    double precision, dimension (nvar), intent (inout) :: xparam
    integer :: i, j, locl, loop
    double precision :: twopi, xref, xstore(maxsym)
    twopi = 4.d0 * Asin (1.d0)
    ndep = 0
    do i = 1, nvar
      if (loc(2, i) == 3) then
!
!  FORCE DIHEDRALS INTO SAME HALF-CIRCLE
!
        j = Int (sign(0.5d0, xparam(i))+xparam(i)/twopi)
        xparam(i) = xparam(i) - j * twopi
      end if
      xstore(i) = xparam(i)
    end do
    do loop = 1, nvar
      if (xstore(loop) >=-1.d4) then
        xref = xstore(loop)
        locl = loc(2, loop)
        do i = loop + 1, nvar
          if (Abs (xref-xstore(i)) < 1.d-3 .and. loc(2, i) == locl) then
            ndep = ndep + 1
            locpar(ndep) = loc(1, loop)
            idepfn(ndep) = locl
            locdep(ndep) = loc(1, i)
            xstore(i) = -1.d5
          end if
        end do
      end if
    end do
    j = 0
    do i = 1, nvar
      if (xstore(i) >-1.d4) then
        j = j + 1
        loc(1, j) = loc(1, i)
        loc(2, j) = loc(2, i)
        xparam(j) = xparam(i)
      end if
      nvar = j
    end do
    return
end subroutine maksym
