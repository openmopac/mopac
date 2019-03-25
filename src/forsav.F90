      subroutine forsav(time, deldip, ipt, fmatrx, coord, nvar, refh, evecs, &
        jstart, fconst) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE chanel_C, only : iw, ires, restart_fn
      use molkst_C, only : numat, norbs
!...Translated by Pacific-Sierra Research 77to90  4.4G  07:59:25  03/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mopend_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ipt 
      integer , intent(in) :: nvar 
      integer  :: jstart 
      real(double)  :: time 
      real(double)  :: refh 
      real(double)  :: deldip(3,*) 
      real(double)  :: fmatrx((nvar*(nvar+1))/2) 
      real(double)  :: coord(nvar) 
      real(double)  :: evecs(nvar**2) 
      real(double)  :: fconst(nvar) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, i99, j, linear, n33, old_numat, old_norbs
      logical :: opend, exists
!-----------------------------------------------
!***********************************************************************
!
!  FORSAV SAVES AND RESTORES DATA USED IN THE FORCE CALCULATION.
!
! ON INPUT TIME = TOTAL TIME ELAPSED SINCE THE START OF THE CALCULATION.
!          IPT  = LINE OF FORCE MATRIX REACHED, IF IN WRITE MODE,
!               = 0 IF IN READ MODE.
!        FMATRX = FORCE MATRIX
!***********************************************************************
      i99 = 0 
   10 continue 
      j = i99 
      inquire(unit=ires, opened=opend) 
      if (opend) close(unit=ires, status='KEEP') 
      if (ipt == 0) then
        inquire (file=restart_fn, exist = exists)
        if ( .not. exists) then
          call mopend ('RESTART file does not exist') 
          return
        end if
      end if
      open(unit=ires, file=restart_fn, status='UNKNOWN', &
      form='UNFORMATTED', iostat=i99, position='asis') 
      if (i99 /= 0) then 
        if (j /= 0) then 
          call mopend ('Fatal error in trying to open RESTART file') 
          return  
        endif 
!
!   A restart file exists, but cannot be mounted.  Delete it.
!   If it can't be deleted, carry on regardless.
!
        open(unit=ires, file=restart_fn, status='OLD', iostat=i99, &
          position='asis') 
        close(ires, status='DELETE', iostat=i99) 
        go to 10 
      endif 
      rewind ires 
      if (ipt == 0) then 
!
!   READ IN FORCE DATA
!
        read (ires, iostat=i99) time, ipt, refh, old_numat, old_norbs
        if (norbs /= old_norbs .or. numat /= old_numat) then
              call mopend("Restart file read in does not match current data set")
              return
        end if
        linear = (nvar*(nvar + 1))/2 
        read (ires, iostat=i99) (coord(i),i=1,nvar) 
        read (ires, iostat=i99) (fmatrx(i),i=1,linear) 
        read (ires, iostat=i99) ((deldip(j,i),j=1,3),i=1,ipt) 
        n33 = nvar*nvar 
        read (ires, iostat=i99) (evecs(i),i=1,n33) 
        read (ires, iostat=i99) jstart, (fconst(i),i=1,nvar) 
        if (i99 /= 0) then
          call mopend ('INSUFFICIENT DATA ON DISK FILES FOR A FORCE CALCULATION RESTART.') 
          write (iw, '(/10X,"PERHAPS THIS STARTED OFF AS A FORCE CALCULATION")') 
          write (iw, '(10X,"BUT THE GEOMETRY HAD TO BE OPTIMIZED FIRST, IN WHICH CASE" &
      &,/10X,"REMOVE THE KEY-WORD ""FORCE"".")') 
        end if
        return  
      else 
!
!    WRITE FORCE DATA
!
        rewind ires 
        if (time > 1.D7) time = time - 1.D7 
        write (ires) time, ipt, refh, numat, norbs
        linear = (nvar*(nvar + 1))/2 
        write (ires) (coord(i),i=1,nvar) 
        write (ires) (fmatrx(i),i=1,linear) 
        write (ires) ((deldip(j,i),j=1,3),i=1,ipt) 
        n33 = nvar*nvar 
        write (ires) (evecs(i),i=1,n33) 
        write (ires) jstart, (fconst(i),i=1,nvar)  
        call den_in_out(1)
        close(ires) 
      endif 
      return  
      end subroutine forsav 
