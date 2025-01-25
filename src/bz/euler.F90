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

subroutine euler (px, py, pz, i1, i2, i3, ipen)
  implicit none
  integer, intent(in) :: i1, i2, i3, ipen
  double precision, intent(in) :: px, py, pz
  ! 
  double precision :: ca = 1.d0, cb = 1.d0, sa = 0.d0, sb = 0.d0
  double precision :: r1, x, xy, y, z
  real :: cx, cy
  double precision, dimension(3) :: c
!
  if (i1 == 0) then
    !
    ! SET UP COS AND SIN ANGLES FOR EULARIAN TRANSFORM.
    !
    x = px
    y = py
    z = pz
    xy = x*x + y*y
    if (xy > 1.d-15) then
      r1 = Sqrt (xy + z*z)
      xy = Sqrt (xy)
      ca = x / xy
      cb = z / r1
      sa = y / xy
      sb = xy / r1
    end if
  else
    !
    !  SWAP AROUND THE POINTS TO BE PLOTTED.
    !
    c(i1) = px
    c(i2) = py
    c(i3) = pz
    x = c(1)
    y = c(2)
    z = c(3)
    !
    !  NOW DO THE EULARIAN TRANSFORM ITSELF.
    !
    cx = real(ca*cb*x + sa*cb*y - sb*z)
    cy = real(-sa*x + ca*y)
    call graphics (cx, cy, ipen)
  end if
end subroutine euler
