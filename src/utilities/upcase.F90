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

      subroutine upcase(keywrd, n)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      character , intent(inout) :: keywrd*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icapa, ilowa, ilowz, i, iline, j
      character :: keybuf*3000
!-----------------------------------------------
!
!  UPCASE WILL TAKE A CHARACTER STRING, IN KEYWRD, AND PUT IT INTO
!  UPPER CASE.  KEYWRD IS LIMITED TO 80 CHARACTERS
!
      icapa = ichar('A')
      ilowa = ichar('a')
      ilowz = ichar('z')
      keybuf = keywrd
      do i = 1, n
        iline = ichar(keywrd(i:i))
        if (iline>=ilowa .and. iline<=ilowz) &
        keywrd(i:i) = char(iline + icapa - ilowa)
        if (iline /= 9) cycle
!
!  Change tabs to spaces.  A tab is ASCII character 9.
!
        keywrd(i:i) = ' '
      end do
!
!   If the word EXTERNAL is present, do NOT change case of the following
!   character string.
!
      i = index(keywrd,'EXTERNAL=')
      if (i /= 0) then
        j = index(keywrd(i+1:),' ') + i
        keywrd(i+9:j) = keybuf(i+9:j)
      end if
      return
      end subroutine upcase
