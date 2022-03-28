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

      subroutine forsav(time, deldip, ipt, fmatrx, coord, nvar, refh, evecs, &
        jstart, fconst)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE chanel_C, only : iw, ires, restart_fn
      use molkst_C, only : numat, norbs
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ipt
      integer , intent(in) :: nvar
      integer  :: jstart
      double precision  :: time
      double precision  :: refh
      double precision  :: deldip(3,*)
      double precision  :: fmatrx((nvar*(nvar+1))/2)
      double precision  :: coord(nvar)
      double precision  :: evecs(nvar**2)
      double precision  :: fconst(nvar)
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
        end if
!
!   A restart file exists, but cannot be mounted.  Delete it.
!   If it can't be deleted, carry on regardless.
!
        open(unit=ires, file=restart_fn, status='OLD', iostat=i99, &
          position='asis')
        close(ires, status='DELETE', iostat=i99)
        go to 10
      end if
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
        close (ires)
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
      end if
      return
      end subroutine forsav
