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

      double precision function seconds (mode)
!
!  If mode = 1:   Get current CPU time - CPU_0, and return without checking for <file>.end
!     mode = 2:   Get current CPU time - CPU_0, and return after checking for <file>.end,
!                 if <file>.end exists, and contains anything, increase time by 10^8 seconds.
!
      USE chanel_C, only : iw, iend, end_fn, output_fn
      use molkst_C, only : is_PARAM, wall_clock_0, wall_clock_1, CPU_1, CPU_0, prt_coords, use_disk
      implicit none
      integer , intent(in) :: mode
      double precision :: shut
      logical :: setok, exists, first = .true., opend
      character :: x
      integer :: i99, l
      integer :: itim1(8), ihour, imin, isec, imsec, jhour = 100, idays = -1
      save
!******************************************************
!
!   SECONDS, ON EXIT, CONTAINS THE NUMBER OF CPU SECONDS
!   SINCE THE START OF THE CALCULATION.
!
!   NOTE: THIS FUNCTION WAS NAMED 'SECOND' IN OLDER
!         MOPAC VERSIONS, WHICH IS THE NAME OF A GNU
!         FORTRAN INTRINSIC FUNCTION
!
!******************************************************
      data setok/ .TRUE./
      data shut/ 0.D0/
      seconds = 0.D0 ! dummy time upon failure
      if (first) then
        first = .false.
!
!   CHECK TO SEE IF AN OLD '.end' FILE EXISTS.  IF IT DOES, DELETE IT.
!
        if (use_disk) then
          rewind iend
          open(unit=iend, file=end_fn, iostat = i99)
          if (i99 /= 0) then
            write(iw,"(a)")" File '"//trim(end_fn)//"' is unavailable for use"
            write(iw,"(a)")" Correct the error and re-submit"
            call to_screen(" File '"//trim(end_fn)//"' is unavailable for use")
            call to_screen(" Correct the error and re-submit")
            call mopend("File '"//trim(end_fn)//"' is unavailable for use")
            return
          end if
          read (iend, '(A)', end=20, err=20, iostat = i99) x
!
!   FILE EXISTS.  DELETE IT.
!
          if (ichar(x) /= (-1)) close(iend, status='DELETE')
        end if
!
!  Set the zero of time
!
20      call cpu_time(CPU_0)
        call date_and_time(values = itim1)
        ihour = itim1(5)
        imin  = itim1(6)
        isec  = itim1(7)
        imsec = itim1(8)
        wall_clock_0 = float(3600*ihour + 60*imin + isec) + (float(imsec)/1.e3 )
      end if

!**********************************************************************
!
!   NOW TO SEE IF A FILE CALLED <FILENAME>.end EXISTS, IF IT DOES
!   THEN INCREMENT CPU TIME BY 10,000,000 SECONDS.
!
!***********************************************************************
      if (mode == 2 .and. use_disk) then
        shut = 0.D0
        inquire (file=end_fn, exist = exists)
        if (exists) then
          open(unit=iend, file=end_fn, iostat=i99)
          if (i99 /= 0) then
            close(iend, err=99)   ! Do not stop job if shutdown is faulty
    99      return
          end if
          read (iend, '(A)', end=10, err=10) x
  !
  !          FILE EXISTS, THEREFORE INCREMENT TIME
  !
          shut = 1.D7
          if (setok) then
            write (iw, '(3/10X,''****   JOB STOPPED BY OPERATOR   ****'')')
            if (prt_coords) then
              write (iw, '(3/10X,''****   GEOMETRY AT THIS POINT    ****'')')
              call geout (1)
            end if
            setok = .FALSE.
          end if
     10   continue
          close(iend, iostat = i99)
        end if
      end if
!
!  Every call to second executes the following code
!
      call cpu_time(CPU_1)
      call date_and_time(values = itim1)
      ihour = itim1(5)
      if (ihour < jhour) idays = idays + 1
      jhour = ihour
      ihour = ihour + 24*idays
      imin  = itim1(6)
      isec  = itim1(7)
      imsec = itim1(8)
      wall_clock_1 = float(3600*ihour + 60*imin + isec) + (float(imsec)/1.e3 )
      seconds = dble(wall_clock_1 - wall_clock_0) + shut

      if (is_PARAM) return
      if (.not. use_disk) return
 11   endfile (iw, iostat = i99)
      l = 0
      if (i99 /= 0) then
        l = l + 1
        write(0,"(i3,3a)")21 - l," File """,trim(output_fn),""" is unavailable for use"
        write(0,"(a)")"    Correct fault (probably the output file is in use elsewhere)"
        call sleep(5)
        inquire(unit=iw, opened=opend)
        if (opend) close (iw)
        open(unit=iw, file=output_fn)
        if (l < 20) goto 11
        write(0,"(a)") "Job abandoned due to output file being unavailable"
        call mopend("Job abandoned due to output file being unavailable")
        call sleep(5)
        return
      end if
      if (l > 0) then
       write(0,"(a)")" Fault successfully corrected.  Job continuing"
      end if
      backspace (iw)
      return
      end function seconds
