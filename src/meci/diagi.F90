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

      double precision function diagi (ialpha, ibeta, eiga, xy, nmos)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmos
      integer , intent(in) :: ialpha(nmos)
      integer , intent(in) :: ibeta(nmos)
      double precision , intent(in) :: eiga(*)
      double precision , intent(in) :: xy(nmos,nmos,nmos,nmos)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      double precision :: x
!-----------------------------------------------
!***********************************************************************
!
!  CALCULATES THE ENERGY OF A MICROSTATE DEFINED BY IALPHA AND IBETA
!
!***********************************************************************
      x = 0.0D0
      do i = 1, nmos
        if (ialpha(i) == 0) cycle
        x = x + eiga(i)
        do j = 1, nmos
          x = x + ((xy(i,i,j,j)-xy(i,j,i,j))*ialpha(j)*0.5D0+xy(i,i,j,j)*ibeta(&
            j))
        end do
      end do
      do i = 1, nmos
        if (ibeta(i) == 0) cycle
        x = x + eiga(i)
        do j = 1, i - 1
          x = x + (xy(i,i,j,j)-xy(i,j,i,j))*ibeta(j)
        end do
      end do
      diagi = x
      return
      end function diagi
