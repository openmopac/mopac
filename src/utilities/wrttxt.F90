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
