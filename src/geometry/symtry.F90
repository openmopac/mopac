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

      subroutine symtry ()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : ndep
      USE symmetry_C, ONLY: locpar, idepfn, locdep, depmul
      use common_arrays_C, only : geo
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, i, locn, j
      double precision :: value
!-----------------------------------------------
!**********************************************************************
!
!  SYMTRY COMPUTES THE BOND LENGTHS AND ANGLES THAT ARE FUNCTIONS OF
!         OTHER BOND LENGTHS AND ANGLES.
!
! ON INPUT GEO     = KNOWN INTERNAL COORDINATES
!          NDEP    = NUMBER OF DEPENDENCY FUNCTIONS.
!          IDEPFN  = ARRAY OF DEPENDENCY FUNCTIONS.
!          LOCDEP  = ARRAY OF LABELS OF DEPENDENT ATOMS.
!          LOCPAR  = ARRAY OF LABELS OF REFERENCE ATOMS.
!
!  ON OUTPUT THE ARRAY "GEO" IS FILLED
!***********************************************************************
!
!     NOW COMPUTE THE DEPENDENT PARAMETERS.
!
      n = 0
      do i = 1, ndep
        if (idepfn(i) == 19 .and. depmul(n + 1) > 1.d-20) then
          n = n + 1
          call haddon (value, locn, idepfn(i), locpar(i), geo, depmul(n))
        else
          call haddon (value, locn, idepfn(i), locpar(i), geo, depmul(i))
        end if
        j = locdep(i)
        geo(locn,j) = value
      end do
      return
      end subroutine symtry
      subroutine haddon(w, l, m, loc, a,  fact)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw
      use molkst_C, only : numcal
      use common_arrays_C, only : na
      use funcon_C, only : pi
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: l
      integer , intent(in) :: m
      integer , intent(in) :: loc
      double precision , intent(out) :: w
      double precision , intent(in) :: fact
      double precision , intent(in) :: a(3,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icalcn, i
      save icalcn
!-----------------------------------------------
      data icalcn/ 0/
      if (numcal /= icalcn) then
        icalcn = numcal
      end if
     !**********************************************************************
    if (m > 19 .or. m < 1) then
      write (iw, "(///10X,'UNDEFINED SYMMETRY FUNCTION',I3, ' USED')") m
      call mopend ("UNDEFINED SYMMETRY FUNCTION USED")
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
