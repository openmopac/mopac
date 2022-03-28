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

      double precision function reada (string, istart)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: istart
      character  :: string*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i9, idot, ineg, ipos, icapd, icape, ismld, ismle, l, i, &
        iadd, n, j
      logical :: expnnt
      double precision, external :: digit
!-----------------------------------------------
!     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
!
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      i0 = ichar('0')
      i9 = ichar('9')
      idot = ichar('.')
      ineg = ichar('-')
      ipos = ichar('+')
      icapd = ichar('D')
      icape = ichar('E')
      ismld = ichar('d')
      ismle = ichar('e')
!
      l = len(string)
!
!     FIND THE START OF THE NUMERIC FIELD
      do i = istart, l
        iadd = 0
        n = ichar(string(i:i))
!
!       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
        if (n>=i0 .and. n<=i9) go to 20
!
!       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
        if (n==ineg .or. n==ipos) then
          iadd = iadd + 1
          if (i + iadd > l) go to 50
          n = ichar(string(i+iadd:i+iadd))
          if (n>=i0 .and. n<=i9) go to 20
        end if
!
!       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
        if (n /= idot) cycle
        iadd = iadd + 1
        if (i + iadd > l) go to 50
        n = ichar(string(i+iadd:i+iadd))
        if (n>=i0 .and. n<=i9) go to 20
      end do
      go to 50
!
!     FIND THE END OF THE NUMERIC FIELD
   20 continue
      expnnt = .FALSE.
      do j = i + 1, l
        iadd = 0
        n = ichar(string(j:j))
!
!       CONTINUE SEARCH FOR END IF DIGIT FOUND
        if (n>=i0 .and. n<=i9) cycle
!
!       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
        if (n==ineg .or. n==ipos) then
          if (.not.expnnt) go to 40
          iadd = iadd + 1
          if (j + iadd > l) go to 40
          n = ichar(string(j+iadd:j+iadd))
          if (n>=i0 .and. n<=i9) cycle
        end if
        if (n == idot) then
          iadd = iadd + 1
          if (j + iadd > l) go to 40
          n = ichar(string(j+iadd:j+iadd))
          if (n>=i0 .and. n<=i9) cycle
          if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) cycle
        end if
        if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) then
          if (expnnt) go to 40
          expnnt = .TRUE.
          cycle
        end if
        go to 40
      end do
      j = l + 1
   40 continue
      n = ichar(string(j-1:j-1))
      if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) j = j - 1
!
!     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
      n = 0
      n = n + index(string(i:j-1),'e')
      n = n + index(string(i:j-1),'E')
      n = n + index(string(i:j-1),'d')
      n = n + index(string(i:j-1),'D')
      if (n == 0) then
        reada = digit(string(i:j-1),1)
      else
        reada = digit(string(:i+n-2),i)*1.D1**digit(string(:j-1),i+n)
      end if
      return
!
!     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
   50 continue
      reada = 0.D0
      return
      end function reada
