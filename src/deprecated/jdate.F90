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

function jdate() result (line)
    !
    !   Simulate the function jdate
    !   On return, contains the year and Julian date in the character*8 form: yyddd
    !   with the last 3 characters unused
    !
    implicit none

    ! Result value
    character (len=8) :: line

    ! External routines
    intrinsic Date_and_Time

    ! Local variables
    integer :: iday, imonth, iyear, ijulian, i
    integer :: dim(12) = (/ 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30/)
    character(len=8) :: date

    line = ' '

    ! This intrinsic routine returns a date-string containing
    !  CCYYMMDD
    call Date_and_Time(date)

    ! Assign the year (in current century, and with leading zeros)
    read (date(3:4), '(i2)')   iyear     
    write(line(1:2), '(i2.2)') iyear

    ! Setup the absolute day in the year (ddd)
    read(date(5:6),'(i2)') imonth 
    read(date(7:8),'(i2)') iday 
    ijulian = iday 
    do i = 1,imonth
        ijulian = ijulian + dim(i)
    end do
    write(line(3:5),'(i3.3)') ijulian

end function jdate

