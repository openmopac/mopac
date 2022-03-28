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
