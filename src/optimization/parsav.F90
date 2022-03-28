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

      subroutine parsav(mode, n, m, q, r, efslst, xlast, iiium)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, ONLY:  keywrd, numat, norbs
      USE chanel_C, only : ires, iw, restart_fn
      use common_arrays_C, only : aicorr, errfn
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mode
      integer  :: n, m, iiium(6)
      double precision  :: q(n,n)
      double precision, dimension(n) :: efslst, xlast
      double precision  :: r(n + 2,n + 2)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  i, j, io_stat, old_numat, old_norbs
      logical :: opend
!-----------------------------------------------
!*********************************************************************
!
!   PARSAV SAVES AND RESTORES DATA USED IN NLLSQ GRADIENT MINIMIZATION.
!
!    IF MODE IS 0 DATA ARE RESTORED, IF 1 THEN SAVED.
!
!*********************************************************************
      inquire(unit=ires, opened=opend)
      if (opend) close(unit=ires, status='KEEP')
      open(unit=ires, file=restart_fn, status='UNKNOWN', form=&
        'UNFORMATTED', position='asis', iostat = io_stat)
      if (io_stat /= 0) then
          write(iw,*)" Restart file either does not exist or is not available for reading"
          call mopend ("Restart file either does not exist or is not available for reading")
          return
        end if
      rewind ires
      if (mode == 0) then
!
!  MODE=0: RETRIEVE DATA FROM DISK.
!
        write (iw, '(2/10X,''RESTORING DATA FROM DISK''/)')
        read (ires, iostat = io_stat) old_numat, old_norbs, (xlast(i),i=1,n), m, iiium, efslst, n
        if (norbs /= old_norbs .or. numat /= old_numat) then
              call mopend("Restart file read in does not match current data set")
              return
        end if
        read (ires, iostat = io_stat) ((q(j,i),j=1,m),i=1,m)
        read (ires, iostat = io_stat) ((r(j,i),j=1,n),i=1,n)
        if (index(keywrd,'AIDER') /= 0) then
          read (ires, iostat = io_stat) (aicorr(i),i=1,n)
          read (ires, iostat = io_stat) (errfn(i),i=1,n)
        end if
        if (io_stat /= 0) then
          write(iw,*)" Restart file is currupt"
          call mopend ("Restart file is currupt")
        end if
        close(ires)
        return
      end if
      if (mode == 1) then
        write (iw, &
      '(2/10X,                                              ''- - - - - - - TIM&
      &E UP - - - - - - -'',2/)')
        write (iw, '(10X,A)') ' - THE CALCULATION IS BEING DUMPED TO DISK', &
          '   RESTART IT USING THE KEY-WORD "RESTART"'
        write (iw, '(/10X,''CURRENT VALUE OF GEOMETRY'',/)')
        call geout (iw)
      end if
      write (ires) numat, norbs, (xlast(i),i=1,n), m, iiium, efslst, n
      write (ires) ((q(j,i),j=1,m),i=1,m)
      write (ires) ((r(j,i),j=1,n),i=1,n)
      if (index(keywrd,'AIDER') /= 0) write (ires) (aicorr(i),i=1,n)
      if (index(keywrd,'AIDER') /= 0) write (ires) (errfn(i),i=1,n)
!*****
!     The density matrix is required by ITER upon restart .
!
      call den_in_out(1)
!
!*****
      close(ires)
      return
      end subroutine parsav
