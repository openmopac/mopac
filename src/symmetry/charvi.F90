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

      double precision function charvi (vects, jorb, ioper, r, nvecs)
!-----------------------------------------------
!   M o d u l elem s
!-----------------------------------------------
      use symmetry_C, only : elem, jelem
      use molkst_C, only : numat
!***********************************************************************
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m elem t elem r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m elem n t s
!-----------------------------------------------
      integer , intent(in) :: jorb
      integer , intent(in) :: ioper
      integer , intent(in) :: nvecs
      double precision , intent(in) :: vects(nvecs,nvecs)
      double precision , intent(in) :: r(3,3)
!-----------------------------------------------
!   L o c a l   V a r i a b l elem s
!-----------------------------------------------
      integer :: iatom, jatom
      double precision, dimension(:), allocatable :: vect1, vect2
      double precision, dimension(5) :: h
      double precision, dimension(3) :: p
      double precision :: sum
! For Mopac BLAS
      double precision, external :: ddot
!

!-----------------------------------------------
      charvi = 1.D0
      if (ioper == 1) return
!
!   Non-trivial case
!
      allocate(vect1(nvecs), vect2(nvecs))
      vect1(:nvecs) = 0.D0
      vect2(:nvecs) = 0.D0
!#      WRITE(IW,*)' Vibration, raw'
!#      WRITE(IW,'(3f12.6)')(VECTS(i,jorb),i=1,3*NUMAT)
      do iatom = 1, numat
        jatom = jelem(ioper,iatom)
        p(1) = vects(iatom*3-2,jorb)
        p(2) = vects(iatom*3-1,jorb)
        p(3) = vects(iatom*3-0,jorb)
        h(1) = r(1,1)*p(1) + r(2,1)*p(2) + r(3,1)*p(3)
        h(2) = r(1,2)*p(1) + r(2,2)*p(2) + r(3,2)*p(3)
        h(3) = r(1,3)*p(1) + r(2,3)*p(2) + r(3,3)*p(3)
        p(1) = elem(1,1,ioper)*h(1) + elem(1,2,ioper)*h(2) + elem(1,3,ioper)*h(3)
        p(2) = elem(2,1,ioper)*h(1) + elem(2,2,ioper)*h(2) + elem(2,3,ioper)*h(3)
        p(3) = elem(3,1,ioper)*h(1) + elem(3,2,ioper)*h(2) + elem(3,3,ioper)*h(3)
        vect1(iatom*3-2:iatom*3) = h(:3)
        vect2(jatom*3-2:jatom*3) = p
      end do
!#      WRITE(IW,*)' Vibration, coordinate independent '
!#      WRITE(IW,'(3f12.6)')(VECT2(i),i=1,3*NUMAT)
!#      WRITE(IW,*)' Vibration, transformed            '
!#      WRITE(IW,'(3f12.6)')(VECT1(i),i=1,3*NUMAT)
!
!   Evaluate the term <psi(i)|R(theta)|psi(i)>
!
      sum = ddot(nvecs,vect1(:nvecs),1,vect2(:nvecs),1)
      charvi = sum
      return
      end function charvi
