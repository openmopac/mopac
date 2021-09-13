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

subroutine newflg ()
    use molkst_C, only: natoms, numat, nvar
    use common_arrays_C, only : coord, geo, txtatm, na, nb, nc, loc, xparam
    use chanel_C, only: iw
    implicit none
    integer :: i, j, k
!
    if (numat /= natoms) then
      write (iw,*) " NEWGEO CAN ONLY BE USED IF THERE ARE NO DUMMY ATOMS,", &
     & " OR IF 'RESEQ' IS USED"
      call mopend ("NEWGEO cannot be used here")
    end if
!
!    Force all coordinates to be internal
!
    call xyzint (coord, numat, na, nb, nc, 1.D0, geo)
!
!  Convert coordinates to Cartesian
!
    call gmetry(geo, coord)
    k = 0
    do i = 1, numat
      if (txtatm(i)(:6) /= "ATOM  ") cycle
      if (txtatm(i)(13:15) /= " N " .and. txtatm(i)(13:15) /= " C ") cycle
!
!  CONVERT ATOM COORDINATES FROM INTERNAL TO CARTESIAN
!
      na(i) = 0
      nb(i) = 0
      nc(i) = 0
      do j = 1, 3
        geo(j, i) = coord(j, i)
      end do
      k = 1
    end do
    do i = 1, nvar
      xparam(i) = geo(loc(2, i), loc(1, i))
    end do
    if (k == 0) call mopend("KEYWORD ""NEWGEO"" USED, BUT NO BACKBONE ATOMS IDENTIFIED")
    return
end subroutine newflg
