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

subroutine rapid0 (loop)
!
!
    use param_global_C, only : valvar, numvar, xparamp, locvar
!
!
    implicit none
    integer, intent (in) :: loop
!----------------------------------------------------------------
    integer :: i
    double precision :: funct1
!----------------------------------------------------------------
  !
  !   DELTAS will hold the CHANGE IN VALUE of the parameters
  !
    do i = 1, numvar
      xparamp(i) = 0.d0
    end do
  !
  !  Optimize the parameters
  !
    call rapid1 (loop, xparamp, numvar, funct1)
  !
  !  Update the values of the parameters
  !
    do i = 1, numvar
      valvar(i) = valvar(i) - xparamp(i)
      if(locvar(1,i) > 3 .and. locvar(1,i) < 7) valvar(i) = max(0.05d0, valvar(i))
      call update (locvar(1, i), locvar(2, i), valvar(i), 0.d0)
    end do
end subroutine rapid0
