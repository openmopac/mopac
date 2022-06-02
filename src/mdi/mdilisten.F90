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

      subroutine mdi_listen()
      use MDI, only : MDI_Init
      USE chanel_C, only : iw

      character(len=1024) :: mdi_options
      integer :: i, ierr
      logical :: use_mdi
!
!  Check for a -mdi command-line option
!
      use_mdi = .false.
      do i = 1, iargc()
        call getarg (i, mdi_options)
        if (mdi_options == '-mdi' .OR. mdi_options == '--mdi') then
          if ( i .eq. iargc() ) then
            WRITE(iw,*)'No argument was provided for the -mdi command-line option'
            stop
          else
            call getarg (i+1, mdi_options)
            write(iw,*) " FOUND -MDI OPTION: ",TRIM(mdi_options)
            CALL MDI_Init( mdi_options, ierr)
            use_mdi = .true.
          endif
        endif
      end do

      return
      end subroutine mdi_listen
