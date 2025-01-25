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
