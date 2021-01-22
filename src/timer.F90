      subroutine timer(a) 
!
!  Write out CPU and wall-clock times for the current step, and total times
!
      USE vast_kind_param, ONLY:  double 
      Use molkst_C, only : wall_clock_0, wall_clock_1, CPU_1, CPU_0
      USE chanel_C, only : iw 
      use second2_I 
      implicit none
      character , intent(in) :: a*(*) 
      character :: b*20
      real(double) :: CPU_start, CPU_last, CPU_now 
      real :: wall_clock_last, dummy
      logical :: first 
      save first, CPU_start, CPU_last, wall_clock_last
      data first/ .TRUE./  
      if (first) then 
!
!  DEFINE THE ZERO OF TIME
!
        CPU_start = second2(1) 
        CPU_start = CPU_0
        CPU_last = CPU_start 
        wall_clock_last = wall_clock_1
        wall_clock_0 = wall_clock_last
        first = .FALSE. 
      endif 
!
! "dummy" call: wall_clock_1 is set in second2(1), but call needed to set CPU_1
!
      dummy = real(second2(1))
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
