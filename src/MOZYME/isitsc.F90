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

subroutine isitsc (escf, selcon, emin, iemin, iemax, okscf, niter, itrmax)
    use molkst_C, only: iscf
    use MOZYME_C, only : ovmax, energy_diff
    implicit none
    logical, intent (out) :: okscf
    integer, intent (in) :: itrmax, niter
    integer, intent (inout) :: iemax, iemin
    double precision, intent (in) :: emin, escf, selcon
    logical :: scf1 = .false.
    integer :: i, iemax1, iemin1
    double precision :: energy_test, fmo_test
    double precision, dimension (10) :: escf0
    data escf0 / 10 * 0.d0 /
    save
   !
   ! Test the change in energy on successive iterations and the maximum
   ! element of the occ-vir block of the Fock matrix in the LMO basis
   !
    energy_test = selcon
    fmo_test = selcon * 5.0d0

    if (ovmax < fmo_test .and. Abs (energy_diff) < energy_test &
         & .and. scf1 .or. niter > itrmax) then
      okscf = .true.
      iscf = 2
      if (scf1) iscf = 1
    else
      scf1 = (ovmax < fmo_test .and. Abs (energy_diff) < energy_test)
      if (emin /= 0.d0) then
        !*****************************************************************
        !
        !  THE FOLLOWING TESTS ARE INTENDED TO ALLOW A FAST EXIT FROM
        !  ITER IF THE RESULT IS 'GOOD ENOUGH' FOR THE CURRENT STEP IN
        !  THE GEOMETRY OPTIMIZATION
        !
        if (escf < emin) then
          !
          !  THE ENERGY IS LOWER THAN THE PREVIOUS MINIMUM.
          !  NOW CHECK THAT IT IS CONSISTENTLY LOWER.
          !
          iemax = 0
          iemin1 = iemin
          iemin = Min (5, iemin+1)
          if (iemin1 == 5) then
            do i = 2, 5
              escf0(i-1) = escf0(i)
            end do
          end if
          escf0(iemin) = escf
          !
          !  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 10%
          !  OF THE ENERGY GAIN FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
          !  MINIMUM.
          !
          if (iemin > 3) then
            do i = 2, iemin
              if (Abs (escf0(i)-escf0(i-1)) > 0.1d0*(emin-escf)) go to 1000
            end do
               !
               ! IS GOOD ENOUGH -- RAPID EXIT
               !
            okscf = .true.
            iscf = 1
            return
          end if
        else
            !
            !  THE ENERGY HAS RISEN ABOVE THAT OF THE PREVIOUS MINIMUM.
            !  WE NEED TO CHECK WHETHER THIS IS A FLUKE OR IS THIS REALLY
            !  A BAD GEOMETRY.
            !
          iemin = 0
          iemax1 = iemax
          iemax = Min (5, iemax+1)
          if (iemax1 == 5) then
            do i = 2, 5
              escf0(i-1) = escf0(i)
            end do
          end if
          escf0(iemax) = escf
          !
          !  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 10%
          !  OF THE ENERGY LOST FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
          !  MINIMUM.
          !
          if (iemax > 3) then
            do i = 2, iemax
              if (Abs (escf0(i)-escf0(i-1)) > 0.1d0*(escf-emin)) go to 1000
            end do
               !
               ! IS GOOD ENOUGH -- RAPID EXIT
               !
            okscf = .true.
            iscf = 1
            return
          end if
        end if
      end if
1000  okscf = .false.
    end if
  end subroutine isitsc
  logical function PLS_faulty()
!
!  When some systems are run using MOZYME, the DIAGG1 - DIAGG2 combination fails to converge,
!  and the ovmax converges to a non-zero minimum.  If the job is stopped and a <file>.den
!  is generated, then on restarting the same job, the fault is automatically corrected.
!
!  PLS_faulty detects the conditions of the failure, at run time, and silently writes out
!  the <file>.den, then after reading in the same file, it re-runs the SCF calculation.
!  This corrects the fault.
!
    use MOZYME_C, only : ovmax
    use molkst_C, only : escf
    implicit none
    integer :: loop = -1, j = 0
    integer, parameter :: loop_lim = 6
    double precision :: ovmax_old = 0.d0, array_ovmax(loop_lim), escf_old = 0.d0, array_escf(loop_lim)
    save
    if (loop == -1) then
      array_ovmax = 10.d0
      array_escf = 10.d0
      loop = 0
    else
      loop = loop + 1
      if (loop == loop_lim + 1) then
        do loop = 2, loop_lim
          array_ovmax(loop - 1) = array_ovmax(loop)
          array_escf(loop - 1) = array_escf(loop)
        end do
        loop = loop_lim
      end if
      array_ovmax(loop) = abs(ovmax - ovmax_old)
      ovmax_old = ovmax
      array_escf(loop) = abs(escf - escf_old)
      escf_old = escf
      do j = 1, loop
        if (array_ovmax(j) > 0.01d0) exit
      end do
      if (j <= loop) then
        do j = 1, loop
          if (array_escf(j) > 0.1d0) exit
        end do
      end if
    end if
!
!  Check to see if (a) ovmax (PLS) has converged, and (b) that ovmax is not near zero.
!
    PLS_faulty = (j > loop .and. ovmax > 0.1d0)
    if (j > loop  .and. ovmax > 0.1d0) then
!
!  Deliberately overwrite array_ovmax with different numbers.
!  This prevents the appearance of convergance, for loop_lim iterations.
!
      do j = 1, loop_lim
        array_ovmax(j) = j
        array_escf(j) = j
      end do
    end if
    return
  end function PLS_faulty
