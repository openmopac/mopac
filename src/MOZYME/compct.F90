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

subroutine compct (nncnew, ncnew, ncmnew, itop, nc, ic, iws, n01, c, n02, &
& nmos, idone, lb, mb, icref, inref)
    implicit none
    integer, intent (in) :: icref, idone, inref, itop, n01, n02, nmos
    integer, intent (out) :: lb, mb
    integer, dimension (n01), intent (inout) :: ic
    integer, dimension (nmos), intent (inout) :: iws, nc, ncmnew, nncnew
    integer, dimension (nmos), intent (out) :: ncnew
    double precision, dimension (n02), intent (inout) :: c
    integer :: i, ii, inco, inew, inic, ioco, ioic, n, natom, ncoef
   !
   !   Start the compression at the LMO below the current LMO.  Gaps can
   !   only exist in the domain IDONE-1 to 1.
   !
    inic = icref
    inco = inref
   !   The last array element of  C used is INCO
   !   The last array element of IC used is INIC
    inew = itop - 1
    ii = idone
    do i = idone - 1, 1, -1
      if (nc(i) /= 0) then
        ii = ii - 1
        inew = inew + 1
        natom = nc(i)
        ncoef = iws(i)
         !
         !  Decrement counters by number of atoms in LMO and by number of
         !  coefficients in LMO.
         !
        inic = inic - natom
        inco = inco - ncoef
         !
         !  Get starting address of atoms and coefficients.
         !
        ioic = nncnew(inew)
        ioco = ncmnew(inew)
         !
         !  Move atom index addresses.  Count backwards to prevent overwriting.
         !
        do n = natom, 1, -1
          ic(inic+n) = ic(ioic+n)
        end do
         !
         !  Move orbital coefficients. Count backwards to prevent overwriting.
         !
        do n = ncoef, 1, -1
          c(inco+n) = c(ioco+n)
        end do
         !
         !  RESET starting addresses
         !
        ncnew(inew) = natom
        nncnew(inew) = inic
        ncmnew(inew) = inco
        nc(ii) = natom
        iws(ii) = ncoef
        if (inew == nmos) exit
      end if
    end do
   !
   !   ZERO OUT ALL LOWER LMOs
   !
    do i = ii - 1, 1, -1
      nc(i) = 0
    end do
   !
   !   RESET LOWER BOUNDARY OF LMOs. (THIS IS USED BY SELMOS)
   !
    lb = inic
    mb = inco
end subroutine compct
