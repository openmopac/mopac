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
    use common_common, only : line, xscale, yscale, xoffset, yoffset, &
      top_l, top_r, bottom_l
   USE IFQWIN
    INTEGER unit
    INTEGER mouseevent
    INTEGER keystate
    INTEGER MouseXpos
    INTEGER MouseYpos
    double precision :: xx, yy
    integer :: status
    xx = (MouseXpos - xoffset)/xscale 
    yy = (MouseYpos - yoffset)/yscale 
    write(line,'(a,f6.2,",",f6.2,",",f6.2,a,f6.2)') &
    "Point (", xx, yy, yy, ") Value:", xx*yy
    status = SETACTIVEQQ (6)
    write(6,*)trim(line)
   return
  end subroutine Details_of_cursor_point
