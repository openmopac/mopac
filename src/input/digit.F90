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

!     ******************************************************************
      double precision function digit (string, istart)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: istart
      character , intent(in) :: string*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i9, ineg, ipos, idot, ispc, l, idig, i, n, j
      double precision :: c1, c2, deciml
      logical :: sign
!-----------------------------------------------
!     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
!     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
!     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
!     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
!
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      i0 = ichar('0')
      i9 = ichar('9')
      ineg = ichar('-')
      ipos = ichar('+')
      idot = ichar('.')
      ispc = ichar(' ')
!
      c1 = 0.D0
      c2 = 0.D0
      sign = .TRUE.
      l = len(string)
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
      idig = 0
      do i = istart, l
        n = ichar(string(i:i))
        if (n>=i0 .and. n<=i9) then
          idig = idig + 1
          c1 = c1*1.D1 + n - i0
        else if (n==ineg .or. n==ipos .or. n==ispc) then
          if (n == ineg) sign = .FALSE.
        else if (n == idot) then
          exit
        else
          go to 40
        end if
      end do
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
      deciml = 1.D0
      do j = i + 1, l
        n = ichar(string(j:j))
        if (n>=i0 .and. n<=i9) then
          deciml = deciml/1.D1
          c2 = c2 + (n - i0)*deciml
        else if (n /= ispc) then
          exit
        end if
      end do
!
!     PUT THE PIECES TOGETHER
   40 continue
      digit = c1 + c2
      if (.not.sign) digit = -digit
      return
      end function digit
