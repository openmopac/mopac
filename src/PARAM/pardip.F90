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

double precision function pardip (coord, nat)
  !**********************************************************************
  !
  !   PARDIP  1:  Computes the dipole moment.
  !
  !**********************************************************************
    use molkst_C, only : numat
    use parameters_c, only : tore
    use common_arrays_C, only : p, q
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    double precision, dimension (3, numat), intent (inout) :: coord
    integer, dimension (numat), intent (in) :: nat
  !
  !.. Local Scalars ..
    integer :: i, l
  !
  !.. Local Arrays ..
    double precision, dimension (3) :: dumy
  !
  !.. External Calls ..
    external chrge
  !
  !.. External Functions ..
    double precision, external :: dipole
  !
  ! ... Executable Statements ...
  !
  !
  !  Calculate the total number of electrons on each atom
  !
    call chrge (p, q)
  !
  !  Calculate the charge on each atom
  !
    do i = 1, numat
      l = nat(i)
      q(i) = tore(l) - q(i)
    end do
  !
  !  Calculate the dipole
  !
    pardip = dipole (p, coord, dumy, 1)
end function pardip
