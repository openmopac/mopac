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
