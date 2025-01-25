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
