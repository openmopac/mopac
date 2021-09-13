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

      double precision function dist2 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!
!     DETERMINE DISTANCES BETWEEN NEIGHBORING ATOMS
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(3) 
      double precision , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dist2 = (a(1)-b(1))**2 + (a(2)-b(2))**2 + (a(3)-b(3))**2 
      return  
      end function dist2 


      double precision function dot1 (a, b) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: a(3) 
      double precision , intent(in) :: b(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      dot1 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3) 
      return  
      end function dot1 
