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

subroutine haddon (w, l, m, loc, a, na, fact)
!**********************************************************************
!
!   HADDON CALCULATES THE VALUE OF A SYMMETRY-DEPENDENT VARIABLE
!
!  ON INPUT: M   = NUMBER SPECIFYING THE SYMMETRY OPERATION
!            LOC = ADDRESS OF REFERENCE ATOM
!            A   = ARRAY OF INTERNAL COORDINATES
!  ON OUTPUT W   = VALUE OF DEPENDENT FUNCTION
!            L   = 1 (FOR BOND LENGTH), 2 (ANGLE), OR 3 (DIHEDRAL)
!**********************************************************************
    use common_systm
    implicit none
    integer, intent (in) :: loc, m, na(natoms + 3)
    integer, intent (out) :: l
    double precision, intent (in) :: fact, a(3,natoms + 3)
    double precision, intent (out) :: w
    integer :: i
    double precision :: pi
    pi = 2.d0 * Asin (1.d0)
    if (m > 19 .or. m < 1) then
      write (iw, "(///10X,'UNDEFINED SYMMETRY FUNCTION',I3, ' USED')") m
    end if
    i = loc
    if (na(i) == 0) then
      select case (m)
      case (1)
!
!   Start of Cartesian relationships.
!
!                          X = X
        l = 1
        w = a(1, i)
        return
      case (2)
!                          Y = Y
        l = 2
        w = a(2, i)
        return
      case (3)
!                          Z = Z
        l = 3
        w = a(3, i)
        return
      case (4)
!                          X = -X
        l = 1
        w = -a(1, i)
        return
      case (5)
!                          Y = -Y
        l = 2
        w = -a(2, i)
        return
      case (6)
!                          Z = -Z
        l = 3
        w = -a(3, i)
        return
      case (7)
!                          X = Y
        l = 1
        w = a(2, i)
        return
      case (8)
!                          Y = Z
        l = 2
        w = a(3, i)
        return
      case (9)
!                          Z = X
        l = 3
        w = a(1, i)
        return
      case (10)
!                          X = -Y
        l = 1
        w = -a(2, i)
        return
      case (11)
!                          Y = -Z
        l = 2
        w = -a(3, i)
        return
      case (12)
!                          Z = -X
        l = 3
        w = -a(1, i)
        return
      case (13)
!                          X = Z
        l = 1
        w = a(3, i)
        return
      case (14)
!                          Y = X
        l = 2
        w = a(1, i)
        return
      case (15)
!                          Z = Y
        l = 3
        w = a(2, i)
        return
      case (16)
!                          X = -Z
        l = 1
        w = -a(3, i)
        return
      case (17)
!                          Y = -X
        l = 2
        w = -a(1, i)
        return
      case (18, 19)
!                          Z = -Y
        l = 3
        w = -a(2, i)
        return
      case default
      end select
    else
      select case (m)
      case (2)
!
!  Set bond-angles equal
!
        l = 2
        w = a(2, i)
        return
      case (3)
!
!  Set dihedrals equal
!
        w = a(3, i)
      case (4)
        w = (pi/2.0d00) - a(3, i)
      case (5)
        w = (pi/2.0d00) + a(3, i)
      case (6)
        w = (2.0d00*pi/3.0d00) - a(3, i)
      case (7)
        w = (2.0d00*pi/3.0d00) + a(3, i)
      case (8)
        w = pi - a(3, i)
      case (9)
        w = pi + a(3, i)
      case (10)
        w = (4.0d00*pi/3.0d00) - a(3, i)
      case (11)
        w = (4.0d00*pi/3.0d00) + a(3, i)
      case (12)
        w = (3.0d00*pi/2.0d00) - a(3, i)
      case (13)
        w = (3.0d00*pi/2.0d00) + a(3, i)
      case (14)
        w = -a(3, i)
      case (15)
        l = 1
        w = a(1, i) / 2.0d00
        return
      case (16)
        l = 2
        w = a(2, i) / 2.0d00
        return
      case (17)
        l = 2
        w = pi - a(2, i)
        return
      case (18, 19)
        l = 1
        w = a(1, i) * fact
        return
      case default
        go to 1000
      end select
      l = 3
      return
    end if
!
!  Set bond-lengths equal
!
1000 l = 1
    w = a(1, i)
end subroutine haddon
