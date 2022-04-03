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

subroutine gettxt
    use common_titles, only : koment, title
    use common_systm, only : ir, iw
    use common_keywrd, only : keywrd
    implicit none
      keywrd = " "
      koment = "    NULL  "
      title = "    NULL  "
      read (ir, "(A)", end=1700, err=1600) keywrd (:120)
      call upcase (keywrd(1:len_trim(keywrd)), len_trim(keywrd))
      if (Index (keywrd, "SETUP") /= 0) then
      else if (Index (keywrd(1:80), " +") /= 0) then
!
!  READ SECOND KEYWORD LINE
!
        read (ir, "(A)", end=1700, err=1600) keywrd (81:160)
        call upcase (keywrd(81:160), 80)
        if (Index (keywrd(81:160), "SETUP") /= 0) then
        else if (Index (keywrd(81:160), " +") /= 0) then
!
!  READ THIRD KEYWORD LINE
!
          read (ir, "(A)", end=1700, err=1600) keywrd (161:240)
          call upcase (keywrd(161:240), 80)
        end if
!
!  READ TITLE LINE
!
      else if (Index (keywrd(:80), "&") /= 0) then
        read (ir, "(A)", end=1700, err=1600) keywrd (81:160)
        call upcase (keywrd(81:160), 80)
        if (Index (keywrd(81:160), "SETUP") /= 0) then
        else if (Index (keywrd(81:160), "&") /= 0) then
          read (ir, "(A)", end=1700, err=1600) keywrd (161:240)
          call upcase (keywrd(161:240), 80)
        else
          read (ir, "(A)", end=1700, err=1600) title
        end if
      else
        read (ir, "(A)", end=1700, err=1600) koment, title
      end if
      
      call upcase (keywrd, 241)
      go to 1800
1600  write (iw, "(A)") " ERROR IN READ OF FIRST THREE LINES"
1700  continue
1800  continue
      return
end subroutine gettxt
