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

      logical function myword (keywrd, testwd)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(inout) :: keywrd*(*)
      character , intent(in) :: testwd*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, k, len_key
      logical :: quote
      character :: key*3000
!-----------------------------------------------
      myword = .FALSE.
      len_key = len_trim(keywrd)
      key = keywrd(:len_key)
!
!      Remove all quoted text from "key" so that, when searching for a keyword, index will not
!      look at any text that is inside quotation marks.
!
      k = 0
      do
        k = k + 1
        if (k > len_key) exit
        if (key(k:k) == '"') then
          j = k
          do
            j = j + 1
            if (j > len_key) exit
            if (key(j:j) == '"') exit
          end do
          key(k:j) = " "
          k = j + 1
        end if
      end do
!
!
!  If the keyword containes quoted text, then ignore text between quotation marks.
!  quote starts off .FALSE.
!  quote is turned .TRUE. if a quotation mark, '"', is found
!  quote is turned .FALSE. when another quotation mark is found.
      quote = .false.
   10 continue
      j = index(key,testwd)
      if (j /= 0) then
!
!  Keyword found.  Now delete keyword in its entirety
!
        key(j:j+1) = " "
20      continue
        do while(keywrd(j:j) == ' ')
          j = j + 1
        end do
        myword = .TRUE.
        do k = j, len_key
          if (keywrd(k:k) == '"') quote = (.not. quote)
          if (.not. quote) then
            if (keywrd(k:k)=='=' .or. keywrd(k:k)==' ') then
!
!     CHECK FOR ATTACHED '=' SIGN
!
              j = k
              if (keywrd(j:j) == '=') go to 50
!
!     CHECK FOR SEPARATED '=' SIGN
!
              do j = k + 1, len_key
                if (keywrd(j:j) == '=') go to 50
                if (keywrd(j:j) /= ' ') go to 10
              end do
!
!    THERE IS NO '=' SIGN ASSOCIATED WITH THIS KEYWORD
!
              go to 10
     50       continue
              keywrd(j:j) = ' '
!
!   THERE MUST BE A NUMBER AFTER THE '=' SIGN, SOMEWHERE
!
              go to 20
            end if
          end if
          keywrd(k:k) = ' '
        end do
      end if
      return
      end function myword
