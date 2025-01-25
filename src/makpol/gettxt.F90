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
