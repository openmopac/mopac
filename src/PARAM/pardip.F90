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
