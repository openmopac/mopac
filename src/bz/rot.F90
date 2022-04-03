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

subroutine rot (xyz)
  use common_common, only : tvec, id
  implicit none 
  double precision, dimension(3), intent(inout) :: xyz
  integer :: i, j, k
  double precision :: sum, x, y, z
  double precision, dimension(3, 3) :: unit
  do i = 1, 3
    sum = 0.d0
    do j = 1, 3
      sum = sum + tvec(j, i)**2
    end do
    sum = 1.d0 / Sqrt (sum)
    do k = 1, 3
      unit(k, i) = tvec(k, i) * sum
    end do
  end do
  sum = 1.d7
  do i = 1, id
    sum = Min (sum, unit(i, i))
  end do
  sum = 1.d0 / sum
  do i = 1, 3
    do j = 1, 3
      unit(i, j) = sum * unit(i, j)
    end do
  end do
  !
  !  UNIT is now the unitary transform for TVEC
  !
  x = xyz(1)*unit(1, 1) + xyz(2)*unit(2, 1) + xyz(3)*unit(3, 1)
  y = xyz(1)*unit(1, 2) + xyz(2)*unit(2, 2) + xyz(3)*unit(3, 2)
  z = xyz(1)*unit(1, 3) + xyz(2)*unit(2, 3) + xyz(3)*unit(3, 3)
  xyz(1) = x
  xyz(2) = y
  xyz(3) = z
end subroutine rot
