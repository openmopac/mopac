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

double precision function digit (string, istart)
  !     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
  !     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
  !     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
  !     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
  !
  ! 
  !.. Implicit Declarations .. 
  implicit none
  ! 
  !.. Formal Arguments .. 
  character(len=*), intent(in) :: string
  integer, intent(in) :: istart
  ! 
  !.. Local Scalars .. 
  logical :: sign
  integer :: i, i0, i9, idig, idot, ineg, ipos, ispc, j, l, n
  double precision :: c1, c2, deciml
  ! 
  !.. Intrinsic Functions .. 
  intrinsic Ichar, Len
  ! 
  ! ... Executable Statements ...
  ! 
  !
  !     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
  i0 = Ichar ("0")
  i9 = Ichar ("9")
  ineg = Ichar ("-")
  ipos = Ichar ("+")
  idot = Ichar (".")
  ispc = Ichar (" ")
  !
  c1 = 0.d0
  c2 = 0.d0
  sign = .true.
  l = Len (string)
  !
  !     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
  idig = 0
  do i = istart, l
    n = Ichar (string(i:i))
    if (n>=i0 .and. n<=i9) then
      idig = idig + 1
      c1 = c1*1.d1 + n - i0
    elseif (n==ineg .or. n==ipos .or. n==ispc) then
      if (n == ineg) then
        sign = .false.
      end if
    else
      goto 1000
    end if
  end do
  !
  !     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
  1100 deciml = 1.d0
  do j = i+1, l
    n = Ichar (string(j:j))
    if (n>=i0 .and. n<=i9) then
      deciml = deciml / 1.d1
      c2 = c2 + (n-i0)*deciml
    elseif (n /= ispc) then
      exit
    end if
  end do
  goto 1200
  1000 if (n == idot) goto 1100
  !
  !     PUT THE PIECES TOGETHER
  1200 digit = c1 + c2
  if (.not. sign) then
    digit = -digit
  end if
end function digit
