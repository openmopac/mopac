      SUBROUTINE FINISH
!
!   MOPEND SHUTS ALL FILES WHICH MAY HAVE BEEN OPENED
!        AND THEN STARTS A RAPID RETURN TO THE MAIN SEGMENT
!            
      use chanel_C, only: iend, end_fn
      use param_global_C, only : ifiles_8
      implicit none
      logical :: exists
      integer :: i
      inquire (file = end_fn, exist = exists)
      if (exists) then
        open(unit=iend, file=end_fn, status='UNKNOWN', position='asis', iostat=i)
        if (i == -100) return
        close(iend, status = 'delete', iostat=i)
        if (i == -100) return
      end if
      end_fn = end_fn(:len_trim(end_fn) - 3)//"res"
      inquire (file = end_fn, exist = exists)
      if (exists) then
        open(unit=iend, file=end_fn, status='UNKNOWN', position='asis', iostat=i)
        if (i == -100) return
        close(iend, status = 'delete', iostat=i)
        if (i == -100) return
      end if
      write (ifiles_8, '(/,'' == PARAM DONE =='')') 
      stop
      END
