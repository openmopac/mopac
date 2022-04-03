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

subroutine getsym (locpar, idepfn, locdep, depmul, na)
    use common_systm, only : line, natoms, ifiles_1,iw, ir, ndep
    implicit none
    integer, dimension (3*natoms), intent (inout) :: idepfn, locpar
    integer, dimension (3*natoms), intent (out) :: locdep
    double precision, dimension (natoms), intent (out) :: depmul
    integer, dimension (natoms + 3), intent (in) :: na
    integer :: i, j, k, l, ll, loop, n, nerror, nvalue
    double precision :: sum
    character (len=60), dimension (19, 2), save :: text
    integer, dimension (40) :: ivalue
    double precision, dimension (40) :: value
    data text (1:19, 1) &
         &/ " BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ", &
         &" BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    " , &
         &" DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE", &
         &" DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS 270 DEGREES + REFERENCE DIHEDRAL  ", &
         &" DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ", &
         &" BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ", &
         &" BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ", &
         &" BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ", &
         &" NOT USED - USE SYMMETRY FUNCTION 19 INSTEAD                ", &
         &" BOND LENGTH IS A MULTIPLE OF THE REFERENCE BOND LENGTH     " /
    data text (1:19, 2) &
         &/ " X COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  " , &
         &" Z COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ", &
         &" X COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ", &
         &" Z COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ", &
         &" X COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ", &
         &" Z COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ", &
         &" X COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ", &
         &" Z COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ", &
         &" X COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ", &
         &" Z COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ", &
         &" X COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ", &
         &" Y COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ", &
         &" Z COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ", &
         &" NOT USED                                                   " /
!
! TITLE OUTPUT
!
    if (ifiles_1 == 1) then
10000 format (///5 x, "PARAMETER DEPENDENCE DATA"//&
           &"        REFERENCE ATOM      FUNCTION NO.    DEPENDENT ATOM(S)")
      write (iw, 10000)
    end if
!
! INPUT SYMMETRY : FUNCTION, REFERANCE PARAMETER, AND DEPENDENT ATOMS
!
    n = 0
    ndep = 0
    nerror = 0
    depmul(1) = 0.d0
    do
      read (ir, "(A)", end=1000) line
      call nuchar (line, value, nvalue)
!   INTEGER VALUES
      do i = 1, nvalue
        ivalue(i) = Nint (value(i))
      end do
!   FILL THE LOCDEP ARRAY
      if (nvalue == 0 .or. ivalue(3) == 0) exit
      if (ivalue(2) == 19) then
        if (na(ivalue(1)) == 0) then
!
!  Not allowed: a Cartesian coordinate cannot use function 19
!
          write (iw,*) "Atom ", ivalue (1), " is Cartesian.  " // &
               &"Function 19 cannot be used here."
          return
        end if
        do i = 4, nvalue
          if (ivalue(i) == 0) exit
          ndep = ndep + 1
          locdep(ndep) = ivalue(i)
          locpar(ndep) = ivalue(1)
          idepfn(ndep) = 19
          n = n + 1
          sum = value(3)**2
          do l = 1, 16
            j = Nint (sum*l)
            if (Abs (j-sum*l) < l*1.d-3) then
              value(3) = Sqrt ((1.d0*j)/l)
            end if
          end do
          depmul(n) = value(3)
        end do
      else
        if (na(ivalue(1)) /= 0 .and. ivalue(2) == 18) then
!
!  Not allowed: an internal coordinate cannot use function 18
!
          write (iw,*) "Atom ", ivalue (1), " is internal.  " // &
               &"Function 18 cannot be used here."
          return
        end if
        do i = 3, nvalue
          if (ivalue(i) == 0) exit
          ndep = ndep + 1
          locdep(ndep) = ivalue(i)
          locpar(ndep) = ivalue(1)
          idepfn(ndep) = ivalue(2)
          if (ivalue(i) > natoms) then
            nerror = 1
          end if
        end do
      end if
      ll = i - 1
      if (ivalue(2) == 19) then
        if (ifiles_1 == 1) then
10010     format (i13, i13, f13.8, i9, 6 i5, 10(/, 43 x, 7 i5))
          write (iw, 10010) ivalue(1), ivalue(2), value(3), &
               &(ivalue(j), j=4, ll)
        end if
      else if (ifiles_1 == 1) then
10020   format (i13, i19, i16, 6 i5, 10(/, 43 x, 7 i5))
        write (iw, 10020) ivalue(1), ivalue(2), (ivalue(j), j=3, ll)
      end if
 !     if (nerror == 1) then
 !       write (iw,*) " A SYMMETRY FUNCTION IS USED TO ", &
 !            &"DEFINE A NON-EXISTENT ATOM"
 !       return
 !     end if
    end do
!
! CLEAN UP
1000 if (ifiles_1 == 1) then
10030 format (/ 10 x, "   DESCRIPTIONS OF THE FUNCTIONS USED", /)
      write (iw, 10030)
    end if
    do loop = 1, 2
      do j = 1, 19
        do i = 1, ndep
          if (idepfn(i) == j) go to 1100
        end do
        cycle
1100    do k = 1, natoms
          if (locpar(i) == k) then
            if (loop == 1 .and. na(k) /= 0) then
              go to 1200
            else if (loop == 2 .and. na(k) == 0) then
              go to 1200
            end if
          end if
        end do
        cycle
1200    if (ifiles_1 == 1) then
10040     format (i4, 5 x, a)
          write (iw, 10040) j, text(j, loop)
        end if
      end do
    end do
end subroutine getsym
