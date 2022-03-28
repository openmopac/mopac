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
