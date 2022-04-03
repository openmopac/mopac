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

subroutine lower (a, n)
  implicit none
  character(len=*), intent(inout) :: a
  integer, intent(in) :: n
  integer :: i, j, lowa, lupa, lupz
  intrinsic Char, Ichar
  lowa = Ichar ("a")
  lupa = Ichar ("A")
  lupz = Ichar ("Z")
  do i = 1, n
    if (Ichar (a(i:i))>=lupa .and. Ichar (a(i:i))<=lupz) then
      j = Ichar (a(i:i)) + lowa - lupa
      a(i:i) = Char(j)
    end if
  end do
end subroutine lower
