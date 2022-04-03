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
