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

      subroutine timout(nout) 
      use molkst_C, only : wall_clock_0, wall_clock_1, CPU_1, CPU_0
!
!     CONVERT THE TIME FROM SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
!
      implicit none
      integer , intent(in) :: nout  
      integer :: CPU_days, CPU_hours, CPU_mins, &
        wall_clock_days, wall_clock_hours, wall_clock_mins 
      double precision :: mins, days, hours, CPU_secs, wall_clock_secs, tim
      character :: day*5, hour*6, minute*8
      tim = CPU_1 - CPU_0
      days = (wall_clock_1 - wall_clock_0)/86400.0 
      wall_clock_days = idint(days) 
      day = " DAY "
      if (wall_clock_days > 1) day = " DAYS"
      hours = (days - dble(wall_clock_days))*24.0D0 
      wall_clock_hours = idint(hours) 
      hour = " HOUR "
      if (wall_clock_hours > 1) hour = " HOURS"
      mins = (hours - dble(wall_clock_hours))*60.0D0 
      wall_clock_mins = idint(mins)
      minute = " MINUTE "
      if (wall_clock_mins > 1) minute = " MINUTES"
      wall_clock_secs = (mins - dble(wall_clock_mins))*60.0D0 
      
      if (wall_clock_days > 0) then 
        write (nout, "(10x, 'WALL-CLOCK TIME', 9x, '=', i2, a, i3, a, i3, a, ' AND',f7.3, ' SECONDS')") &
        wall_clock_days, trim(day), wall_clock_hours, trim(hour), wall_clock_mins, trim(minute), wall_clock_secs 
      else if (wall_clock_hours > 0) then 
        write (nout, "(10x, 'WALL-CLOCK TIME', 9x, '=',  i3, a, i3, a, ' AND',f7.3, ' SECONDS')") &
        wall_clock_hours, trim(hour), wall_clock_mins, trim(minute), wall_clock_secs 
      else if (wall_clock_mins > 0) then 
         write (nout, "(10x, 'WALL-CLOCK TIME', 9x, '=', i3, a, ' AND',f7.3, ' SECONDS')") &
         wall_clock_mins, trim(minute), wall_clock_secs 
      else 
        write (nout, "(10x, 'WALL-CLOCK TIME', 9x, '=', f15.3, ' SECONDS')") wall_clock_secs 
      end if 
      
      days = tim/86400.0 
      CPU_days = idint(days) 
      day = " DAY "
      if (CPU_days > 1) day = " DAYS"
      hours = (days - dble(CPU_days))*24.0D0 
      CPU_hours = idint(hours) 
      hour = " HOUR "
      if (CPU_hours > 1) hour = " HOURS"
      mins = (hours - dble(CPU_hours))*60.0D0 
      CPU_mins = idint(mins) 
      minute = " MINUTE "
      if (CPU_mins > 1) minute = " MINUTES"
      CPU_secs = (mins - dble(CPU_mins))*60.0D0 
      
      if (CPU_days > 0) then 
        write (nout, "(10x, 'COMPUTATION TIME', 8x, '=', i2, a, i3, a, i3, a, ' AND',f7.3, ' SECONDS')") &
        CPU_days, trim(day), CPU_hours, trim(hour), CPU_mins, trim(minute), CPU_secs 
      else if (CPU_hours > 0) then 
        write (nout, "(10x, 'COMPUTATION TIME', 8x, '=',  i3, a, i3, a, ' AND',f7.3, ' SECONDS')") & 
        CPU_hours, trim(hour), CPU_mins, trim(minute), CPU_secs 
      else if (CPU_mins > 0) then 
         write (nout, "(10x, 'COMPUTATION TIME', 8x, '=', i3, a, ' AND',f7.3, ' SECONDS')") &
         CPU_mins, trim(minute), CPU_secs 
      else 
        write (nout, "(10x, 'COMPUTATION TIME', 8x, '=', f15.3, ' SECONDS')") CPU_secs 
      end if 
      return
      end subroutine timout 
