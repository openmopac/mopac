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
  !     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
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
  logical :: expnnt
  integer :: i, i0, i9, iadd, icapd, icape, idot, ineg, ipos, ismld, ismle, &
       & j, l, n
  double precision, external :: digit
  ! 
  ! ... Executable Statements ...
  ! 
  !
  !     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
  i0 = Ichar ("0")
  i9 = Ichar ("9")
  idot = Ichar (".")
  ineg = Ichar ("-")
  ipos = Ichar ("+")
  icapd = Ichar ("D")
  icape = Ichar ("E")
  ismld = Ichar ("d")
  ismle = Ichar ("e")
  !
  l = Len (string)
  !
  !     FIND THE START OF THE NUMERIC FIELD
  do i = istart, l
    iadd = 0
    n = Ichar (string(i:i))
    !
    !       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
    if (n>=i0 .and. n<=i9) goto 1000
    !
    !       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
    if (n==ineg .or. n==ipos) then
      iadd = iadd + 1
      if (i+iadd > l) exit
      n = Ichar (string(i+iadd:i+iadd))
      if (n>=i0 .and. n<=i9) goto 1000
    end if
    !
    !       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
    if (n == idot) then
      iadd = iadd + 1
      if (i+iadd > l) exit
      n = Ichar (string(i+iadd:i+iadd))
      if (n>=i0 .and. n<=i9) goto 1000
    end if
  end do
  !
  !     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
  reada = 0.d0
  return
  !
  !     FIND THE END OF THE NUMERIC FIELD
  1000 expnnt = .false.
  do j = i+1, l
    iadd = 0
    n = Ichar (string(j:j))
    !
    !       CONTINUE SEARCH FOR END IF DIGIT FOUND
    if (n<i0 .or. n>i9) then
      !
      !       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
      if (n==ineg .or. n==ipos) then
        if (.not. expnnt) goto 1100
        iadd = iadd + 1
        if (j+iadd > l) goto 1100
        n = Ichar (string(j+iadd:j+iadd))
        if (n>=i0 .and. n<=i9) cycle
      end if
      if (n == idot) then
        iadd = iadd + 1
        if (j+iadd > l) goto 1100
        n = Ichar (string(j+iadd:j+iadd))
        if (n>=i0 .and. n<=i9) then
          cycle
        elseif (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) then
          cycle
        end if
      end if
      if (n/=icape .and. n/=ismle .and. n/=icapd .and. n/=ismld) goto 1100
      if (expnnt) goto 1100
      expnnt = .true.
    end if
  end do
  j = l + 1
  1100 n = Ichar (string(j-1:j-1))
  if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) then
    j = j - 1
  end if
  !
  !     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
  n = 0
  n = n + Index (string(i:j-1), "e")
  n = n + Index (string(i:j-1), "E")
  n = n + Index (string(i:j-1), "d")
  n = n + Index (string(i:j-1), "D")
  if (n == 0) then
    reada = digit (string(i:j-1), 1)
  else
    reada = digit (string(:i+n-2), i) * 1.d1**digit (string(:j-1), i+n)
  end if
end function reada
