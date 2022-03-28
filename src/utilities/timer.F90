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

      subroutine timer(a)
!
!  Write out CPU and wall-clock times for the current step, and total times
!
      Use molkst_C, only : wall_clock_0, wall_clock_1, CPU_1, CPU_0
      USE chanel_C, only : iw
      implicit none
      character , intent(in) :: a*(*)
      character :: b*20
      double precision :: CPU_start, CPU_last, CPU_now
      real :: wall_clock_last, dummy
      logical :: first
      double precision, external :: seconds
      save first, CPU_start, CPU_last, wall_clock_last
      data first/ .TRUE./
      if (first) then
!
!  DEFINE THE ZERO OF TIME
!
        CPU_start = seconds(1)
        CPU_start = CPU_0
        CPU_last = CPU_start
        wall_clock_last = wall_clock_1
        wall_clock_0 = wall_clock_last
        first = .FALSE.
      end if
!
! "dummy" call: wall_clock_1 is set in seconds(1), but call needed to set CPU_1
!
      dummy = real(seconds(1))
      if (dummy < -300.0) return
      CPU_now = CPU_1
      b = a
      write (iw, '(2x, a20, a, f7.2,  a, f11.2, a, f7.2, a, f9.2)') b, 'CPU INTERVAL:', &
      CPU_now - CPU_last, ', INTEGRAL:', CPU_now - CPU_start, &
      '    WALL-CLOCK INTERVAL:', wall_clock_1 - wall_clock_last, ', INTEGRAL:', &
      wall_clock_1 - wall_clock_0
      CPU_last = CPU_now
      wall_clock_last = wall_clock_1
      return
      end subroutine timer
