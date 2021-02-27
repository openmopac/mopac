program mopac
  use chanel_C, only : iw0
  use molkst_C, only : program_name, gui, line, verson !, site_no, academic
  implicit none
  program_name = "Standalone MOPAC "
  call getdatestamp(line, verson)
#ifdef MOPAC_OS
  verson(7:7) = MOPAC_OS
#else
  verson(7:7) = "X"
#endif
  gui = .false.
  iw0 = -1
!  iw0 = 0
!
! The call to "password" checks that the password is valid.  If it's not valid, the run will be stopped.
! if it is valid, ijulian will be incremented to, e.g. ijulian = 365.
! To by-pass "password" replace the following line with: "site_no = 999".   
!  call password 
!  academic = .false.
!  site_no = 999 
  call run_mopac 
end program mopac
  
