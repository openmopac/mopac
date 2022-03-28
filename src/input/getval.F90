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

      subroutine getval(line, x, t)
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
      double precision , intent(out) :: x
      character  :: line*80
      character , intent(out) :: t*12
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      character :: ch1, ch2
      double precision, external :: reada
!-----------------------------------------------
      ch1 = line(1:1)
      ch2 = line(2:2)
      if ((ichar(ch1)<ichar('A') .or. ichar(ch1)>ichar('Z')) .and. (ichar(ch2)<&
        ichar('A') .or. ichar(ch2)>ichar('Z'))) then
!
!   IS A NUMBER
!
        x = reada(line,1)
        t = ' '
      else
        i = index(line,' ')
        t = line(:i)
        x = -999.D0
      end if
      return
      end subroutine getval
