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

subroutine Details_of_cursor_point(unit, mouseevent, keystate, MouseXpos,MouseYpos)
!
! Generate and print details of a point on a contour-map.  The point is selected by a mouse-click
! on the map
!
  use common_common, only : line, xscale, yscale, xoffset, yoffset, &
  top_l, top_r, bottom_l, grid, isurf, phonon, units
  USE IFQWIN
  INTEGER unit
  INTEGER mouseevent
  INTEGER keystate
  INTEGER MouseXpos
  INTEGER MouseYpos
  real :: xx, yy, xyz(3), sum
  integer :: status, npoint = -1
  character :: letter*1, num*3
  save :: npoint
  npoint = npoint + 1  
  xx = 0.5*(MouseXpos - xoffset)/xscale 
  yy = 0.5*(MouseYpos - yoffset)/yscale 
  if (xx > 1.0 .or. xx < 0.0 .or. yy > 1.0 .or. yy < 0.0 ) then
    line = "EXIT"
    npoint = -1
    return
  end if
  status = SETACTIVEQQ (10) ! unit = 10
  if (npoint > 25) then
    letter = char(ichar("A") + npoint - 26)
  else
    letter = char(ichar("a") + npoint)
  end if
  write(line,'(a)')letter
  call graphics(2.0*xx - 0.01, 2.0*yy - 0.02, 96)
  xyz = top_l - xx*(top_l - top_r) - yy*(top_l - bottom_l) 
  sum = grid(81 - nint(80.0*xx),81 - nint(80.0*yy), isurf)
  if (phonon) then
    num = "7.1"
  else
     num = "9.4"
  end if
  write(line,'(a,f6.2,",",f6.2,",",f6.2,a,f'//num//',a)') &
  "     Point "//letter//" = (", xyz, ")   Value:", sum, units
  status = SETACTIVEQQ (6)
  write(6,*)trim(line)
  return
end subroutine Details_of_cursor_point
