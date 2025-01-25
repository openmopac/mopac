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

double precision function pargeo (grad, wtgeo, refgeo, loc, diffs, ndif)
  !***********************************************************************
  !
  !   PARGEO  1:  Computes the gradient error.
  !
  !***********************************************************************
    use molkst_C, only : nvar
    use param_global_C, only : power, wtz
    implicit none
    double precision, dimension (300), intent (in) :: grad
    double precision, intent (in) :: wtgeo
    character (len=12), dimension (300), intent (in) :: refgeo
    integer, dimension (2, 300), intent (in) :: loc
    double precision, dimension (300), intent (inout) :: diffs
    integer, intent (out) :: ndif
    integer :: i, lim
    double precision :: relative_weight, sum
    intrinsic Min
!------------------------------------------------------------------------
    lim = Min (100, nvar)
    ndif = 0
    sum = 0.d0
    do i = 1, lim
      if (refgeo(i) /= " ") then
        ndif = ndif + 1
        if (loc(2, i) == 1) then
                                      !
          relative_weight = 1.d0      !  Relative weight of bond-length
                                      !
        else if (loc(2, i) == 2) then !
                                      !
          relative_weight = 1.d0      !  Relative weight of bond-angle
                                      !
        else if (loc(2,i) == 3) then  !
                                      !
          relative_weight = 1.d0      !  Relative weight of dihedral
                                      !
        end if
        diffs(ndif) = grad(i) * wtgeo * relative_weight * wtz
        sum = sum + Abs(diffs(ndif)) ** power
      end if
    end do
    pargeo = sum
end function pargeo
