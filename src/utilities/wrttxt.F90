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

      subroutine wrttxt(iprt)
      use molkst_C, only : koment, title, refkey, keywrd, line
      implicit none
      integer , intent(in) :: iprt
      integer :: i, j, k
      logical :: l_chains = .false., l_start = .false.
!-----------------------------------------------
!
! Is CHAINS or START_RES present in the data set?
!
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        line = " "//trim(refkey(i))
        call upcase(line, len_trim(line))
        if (.not. l_chains) l_chains = (index(line, " CHAINS") /= 0)
        if (.not. l_start)  l_start  = (index(line, " START_RES") /= 0)
      end do
!
!  Is CHAINS present in the keyword?
!
      i = index(keywrd," CHAINS")
      if(i /= 0 .and. .not. l_chains) then
!
!  CHAINS is present in keyword, but was not present in the data-set,
!  so add CHAINS keyword to refkey(1)
!
        j = index(keywrd(i + 7:), ")") + i + 7
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if
!
!  Is START_RES present in the keyword?
!
      i = index(keywrd," START_RES")
      if(i /= 0 .and. .not.l_start) then
!
!  START_RES is present in keyword, but was not present in the data-set,
!  so add START_RES keyword to refkey(1)
!
        j = index(keywrd(i + 10:), ")") + i + 10
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if
      if (refkey(2) == " ") then
        i = index(refkey(1), " +")
        if (i /= 0) then
          refkey(1)(i:i + 1) = " "
          refkey(2) = " NULL"
        end if
      end if
      k = 0
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        if (index(refkey(i), " +") == 0) k = k + 1
        write(iprt,'(a)', iostat = j)trim(refkey(i))
        if (j /= 0) then
          call mopend("ERROR DETECTED WHILE TRYING TO WRITE KEYWORDS TO A FILE")
          return
        end if
      end do
      if (index(koment, " NULL") == 0 .and. k < 3) write (iprt, '(A)') trim(koment)
      if (index(koment, " NULL") == 0 .and. k < 4) write (iprt, '(A)') trim(title)
      return
      end subroutine wrttxt
