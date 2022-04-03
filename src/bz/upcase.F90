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

subroutine upcase (keywrd)
  implicit none
  character(len=*), intent(inout) :: keywrd
!
  integer :: i, icapa, iline, ilowa, ilowz
  icapa = Ichar ("A")
  ilowa = Ichar ("a")
  ilowz = Ichar ("z")
  do i = 1, Len (keywrd)
    iline = Ichar (keywrd(i:i))
    if (iline>=ilowa .and. iline<=ilowz) then
      keywrd(i:i) = Char(iline+icapa-ilowa)
    end if
  end do
end subroutine upcase
