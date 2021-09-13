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

subroutine optgeo (xparam, yparam, nvar, refgeo, glim)
  !**********************************************************************
  !
  !   OPTGEO 1:  Optimizes the current geometry
  !              If the geometry is a geometry reference, the reference
  !              data are restored before exit.
  !  On exit, xparam holds the optimized geometry; geometric references are
  !           NOT changed.
  !           yparam holds to optimized geometric references
  !
  !**********************************************************************
    use param_global_C, only : ifiles_8, molnam, contrl
    use molkst_C, only : moperr, errtxt, escf, keywrd
    use common_arrays_C, only : gradients => grad, loc
    implicit none
    integer, intent (inout) :: nvar
    double precision, intent (in) :: glim
    character (len=12), dimension (300), intent (in) :: refgeo
    double precision, dimension (nvar), intent (inout) :: xparam
    double precision, dimension (nvar), intent (out) ::  yparam
    logical :: lsav
    integer :: i, lim, j
    double precision :: sum
    double precision, dimension (nvar) :: x_store
    integer, dimension (2,nvar) :: loc_store
    external compfg, ef
    double precision, external :: dot
    intrinsic Min, Sqrt
!-------------------------------------------------------------------------
  !
  !  Compute the gradient.
  !
       lim = Min (300, nvar)
       j = 0
       do i = 1, lim
        if (refgeo(i) /= " ") j = j + 1
      end do
      if(nvar - j > 0 .or. nvar == 0)then
        escf = 0.d0
        call compfg (xparam, .true., escf, .true., gradients, .true.)
          if (moperr) then
            write (ifiles_8, "(2A)") " Deadly error detected in COMPFG in OPTGEO, in system ", &
           & molnam
            write (ifiles_8, "(2A)") " Error message:", errtxt (1:50)
          end if
  !
  !  Compute the GNORM.  First, set gradients due to geometric references to zero.
  !
          do i = 1, lim
            if (refgeo(i) /= " ") then
            gradients(i) = 0.d0
          end if
        end do
      if (nvar > 0) then
        sum = Sqrt (dot(gradients, gradients, nvar))
      else
        sum = 0.d0
      end if
    else
      sum=0.d0
    end if
  !
  !  Store the geometric variables - they might be used if PARTAB is
  !  called.
  !
    do i = 1, nvar
      x_store(i) = xparam(i)
      loc_store(1,i) = loc(1,i)
      loc_store(2,i) = loc(2,i)
    end do
  !
  !  If the GNORM is large, re-optimize the geometry
  !
    lsav = (sum > glim)
    if (lsav) then
      if( glim > 0.d0) then
    !
    !  Arrange to optimize only the variables that are NOT geometric references
    !  This is for parameter optimization.  If PARTAB is called, then glim < 0
    !  and the whole geometry is optimized.
    !
        j = 0
        do i = 1,nvar
          if (i > lim .or. refgeo(i) == " ") then
            j = j + 1
            loc(1,j) = loc(1,i)
            loc(2,j) = loc(2,i)
            xparam(j) = xparam(i)
          end if
        end do
        i = nvar
        nvar = j
      end if
      if(nvar == 1) then
       call flepo (xparam, nvar, escf)
      else if(nvar > 1) then
        if (index(keywrd,' LBFGS') + index(contrl,' LBFGS')  /= 0) then
          call lbfgs(xparam, escf)
        else
          call ef (xparam, escf)
        end if
      end if
      if (glim < 0) then
!
!  This is used for PARTAB: restore reference data to xparam,
!  put optimized geometry into yparam
!
        do i = 1, lim
          if (refgeo(i) /= " ") then
            yparam(i) = xparam(i)
            xparam(i) = x_store(i)
            if (abs(gradients(i)) < 1.d-4) gradients(i) = 1.d-4
          end if
        end do
      else
!
!  Restore geometry.  x_store contains the reference geometric data plus the
!  un-optimized geometry.  xparam contains the optimized geometry
!
        nvar = i
        j = 0
        do i = 1, nvar
          if (i > lim .or. refgeo(i) == " ") then
            j = j + 1
            x_store(i) = xparam(j)
          end if
        end do
        do i = 1, nvar
          loc(1,i) = loc_store(1,i)
          loc(2,i) = loc_store(2,i)
          xparam(i) = x_store(i)
        end do
        if (moperr) then
          write (ifiles_8, "(2A)") " Deadly error detected in EF, in system ", molnam
          write (ifiles_8, "(2A)") " Error message:", errtxt (1:50)
        end if
      end if
    end if
    return
end subroutine optgeo
