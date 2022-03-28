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

      subroutine xyzint(xyz, numat, na, nb, nc, degree, geo)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use funcon_C, only : pi
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision  :: degree
      integer :: numat
      integer  :: na(numat)
      integer  :: nb(numat)
      integer  :: nc(numat)
      double precision  :: xyz(3,numat)
      double precision  :: geo(3,numat)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, im1, k, l
      double precision :: sum, r
!-----------------------------------------------
!**********************************************************************
!
! XYZINT WORKS OUT THE INTERNAL COORDINATES OF A MOLECULE.
!        THE "RULES" FOR THE CONNECTIVITY ARE AS FOLLOWS:
!        ATOM I IS DEFINED AS BEING AT A DISTANCE FROM THE NEAREST
!        ATOM J, ATOM J ALREADY HAVING BEEN DEFINED.
!        ATOM I MAKES AN ANGLE WITH ATOM J AND THE ATOM K, WHICH HAS
!        ALREADY BEEN DEFINED, AND IS THE NEAREST ATOM TO J
!        ATOM I MAKES A DIHEDRAL ANGLE WITH ATOMS J, K, AND L. L HAVING
!        BEEN DEFINED AND IS THE NEAREST ATOM TO K, AND J, K AND L
!        HAVE A CONTAINED ANGLE IN THE RANGE 15 TO 165 DEGREES,
!        IF POSSIBLE.
!
!        IF(NA(2).EQ.1 THEN THE ORIGINAL CONNECTIVITY IS USED.
!
!        NOTE THAT GEO AND XYZ MUST NOT BE THE SAME IN THE CALL.
!
!   ON INPUT XYZ    = CARTESIAN ARRAY OF NUMAT ATOMS
!            DEGREE = 1 IF ANGLES ARE TO BE IN RADIANS
!            DEGREE = 57.29578 IF ANGLES ARE TO BE IN DEGREES
!
!**********************************************************************
!
!  Use original connectivity, if present
!
        geo = 0.d0
        do i = 1, numat
          im1 = i - 1
          j = na(i)
          if (j == 0) then
            na(i) = 2
            nb(i) = 3
            nc(i) = 4
            if (i == 1) cycle
            sum = 1.d30
            do j = 1, im1
              r = (xyz(1,i) - xyz(1,j))**2 + (xyz(2,i) - xyz(2,j))**2 + (xyz(3,i) - xyz(3,j))**2
              if (.not.(r < sum .and. na(j) /= j .and. nb(j) /=j )) cycle
              sum = r
              k = j
            end do
  !
  !   ATOM I IS NEAREST TO ATOM K
  !
            na(i) = k
            j = k
            if (i > 2) nb(i) = na(k)
            nc(i) = nb(k)
!
!   FIND ANY ATOM TO RELATE TO NA(I)
!
          end if
          if (i > 3) then
!
!  Check that the angle na(i) - nb(i) - nc(i) is meaningful
!
            call bangle (xyz, na(i), nb(i), nc(i), sum)
            if (sum < 1.d-2 .or. Abs(pi - sum) < 1.d-2) then
!
!  Angle is zero or 180 degrees.  Search for an atom nearest to 90 degrees
!
              r = 2.d0
              l = 0
              do k = 1, im1
                if (k == na(i) .or. k == nb(i)) cycle
                call bangle (xyz, na(i), nb(i), k, sum)
                if (Abs(pi*0.5d0 - sum)  < r) then
                  r =  Abs(pi*0.5d0 - sum)
                  l = k
                  end if
                  if (r < 0.5d0) exit ! angle is acceptable
              end do
              nc(i) = l
            end if
            call dihed (xyz, i, j, nb(i), nc(i), geo(3,i))
          end if
          geo(3,i) = geo(3,i)*degree
          if (i > 2) call bangle (xyz, i, j, nb(i), geo(2,i))
          geo(2,i) = geo(2,i)*degree
          geo(1,i) = sqrt((xyz(1,i)-xyz(1,j))**2 + &
                          (xyz(2,i)-xyz(2,j))**2 + &
                          (xyz(3,i)-xyz(3,j))**2)
        end do
        na(1) = 0
        nb(1) = 0
        nc(1) = 0
        if (numat > 1) then
          nb(2) = 0
          nc(2) = 0
          if (numat > 2) nc(3) = 0
          na(2) = 1
        end if
        return
      end subroutine xyzint
