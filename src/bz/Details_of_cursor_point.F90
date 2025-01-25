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
