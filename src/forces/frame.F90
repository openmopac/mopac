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

      subroutine frame(fmat, numat, mode)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : coord, atmass
      use molkst_C, only : id, n_trivial => itemp_2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: numat
      integer  :: mode
      double precision , intent(inout) :: fmat(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, n3, l
      double precision, dimension(6,3*numat) :: vib
      double precision, dimension(3,3) :: rot
      double precision :: shift(6)
      double precision, dimension(3,numat) :: coord1
      double precision :: sum, wtmass, x, y, z, sums(6)
!-----------------------------------------------
!**********************************************************************
!
!   FRAME APPLIES AN RIGID ORIENTATION TO THE MOLECULE IN A FORCE
!         CALCULATION. THE TRANSLATIONS ARE GIVEN A 'FORCE CONSTANT'
!         OF T(X)=500 MILLIDYNES/ANGSTROM
!            T(Y)=600 MILLIDYNES/ANGSTROM
!            T(Z)=700 MILLIDYNES/ANGSTROM
!         AND THE ROTATIONS ARE GIVEN A 'FORCE CONSTANT' OF
!            R(X)=800 MILLIDYNES/ANGSTROM
!            R(Y)=900 MILLIDYNES/ANGSTROM
!            R(Z)=1000 MILLIDYNES/ANGSTROM,
!    THE ROTATIONS ARE MADE ABOUT AXES DETERMINED BY THE MOMENTS
!    OF INERTIA, WHICH IN TURN DEPEND ON THE ISOTOPIC MASSES. FOR
!    THE NORMAL FREQUENCY CALCULATION THESE ARE THE REAL MASSES,
!    FOR THE FORCE CALCULATION THEY ARE ALL UNITY.
!**********************************************************************
      call axis (x, y, z, rot)
      do i = 1, numat
        do j = 1, 3
          sum = 0.D0
          do k = 1, 3
            sum = sum + coord(k,i)*rot(k,j)
          end do
          coord1(j,i) = sum
        end do
      end do
      n3 = numat*3
      j = 0
      wtmass = 1.D0
      if (mode == 1) then
        do i = 1, numat
          wtmass = sqrt(atmass(i))
          j = j + 1
          vib(1,j) = wtmass
          vib(2,j) = 0.D0
          vib(3,j) = 0.D0
          vib(4,j) = 0.D0
          vib(5,j) = coord1(3,i)*wtmass
          vib(6,j) = coord1(2,i)*wtmass
          j = j + 1
          vib(1,j) = 0.D0
          vib(2,j) = wtmass
          vib(3,j) = 0.D0
          vib(4,j) = coord1(3,i)*wtmass
          vib(5,j) = 0.D0
          vib(6,j) = -coord1(1,i)*wtmass
          j = j + 1
          vib(1,j) = 0.D0
          vib(2,j) = 0.D0
          vib(3,j) = wtmass
          vib(4,j) = -coord1(2,i)*wtmass
          vib(5,j) = -coord1(1,i)*wtmass
          vib(6,j) = 0.D0
        end do
      else
! Translation along "x"
        vib(1,j+1:numat*3-2+j:3) = wtmass
        vib(1,j+2:numat*3-1+j:3) = 0.D0
        vib(1,j+3:numat*3+j:3) = 0.D0
! Translation along "y"
        vib(2,j+1:numat*3-2+j:3) = 0.D0
        vib(2,j+2:numat*3-1+j:3) = wtmass
        vib(2,j+3:numat*3+j:3) = 0.D0
! Translation along "z"
        vib(3,j+1:numat*3-2+j:3) = 0.D0
        vib(3,j+2:numat*3-1+j:3) = 0.D0
        vib(3,j+3:numat*3+j:3) = wtmass
! Rotation about "x"
        vib(4,j+1:numat*3-2+j:3) = 0.D0
        vib(4,j+2:numat*3-1+j:3) =  coord1(3,:numat)*wtmass
        vib(4,j+3:numat*3+j:3)   = -coord1(2,:numat)*wtmass
! Rotation about "y"
        vib(5,j+1:numat*3-2+j:3) =  coord1(3,:numat)*wtmass
        vib(5,j+2:numat*3-1+j:3) =  0.D0
        vib(5,j+3:numat*3+j:3)   = -coord1(1,:numat)*wtmass
! Rotation about "z"
        vib(6,j+1:numat*3-2+j:3) =  coord1(2,:numat)*wtmass
        vib(6,j+2:numat*3-1+j:3) = -coord1(1,:numat)*wtmass
        vib(6,j+3:numat*3+j:3)   =  0.D0
      end if
!
! Unitary transform to rotate vibrations into frame of moments of inertia.
!
      j = 1
      do i = 1, numat
        do k = 4, 6
          x = vib(k,j)
          y = vib(k,j+1)
          z = vib(k,j+2)
          vib(k,j)   = x*rot(1,1) + y*rot(1,2) + z*rot(1,3)
          vib(k,j+1) = x*rot(2,1) + y*rot(2,2) + z*rot(2,3)
          vib(k,j+2) = x*rot(3,1) + y*rot(3,2) + z*rot(3,3)
        end do
        j = j + 3
      end do
      sums = 0.D0
      do i = 1,n3
        sums(:) = sums(:) + vib(:,i)**2
      end do
      if (n_trivial == 5) then
        sum = 0.01d0 + Min(sums(1), sums(2), sums(3), sums(4), sums(5), sums(6))  ! System is linear, so don't use one mode
        do i = 1,6
          if (sums(i) > sum) sums(i) = sqrt(1.D0/sums(i))
        end do
      else
        sums = sqrt(1.d0/sums)
      end if

      if (id /= 0) sums(4:6) = 0.D0

      do i = 1, 6
        vib(i,:n3) = vib(i,:n3)*sums(i)
        shift(i) = 40000.D0 + i*100.D0
      end do
      l = 0
      do i = 1, n3
        do j = 1, i
          l = l + 1
          sum = 0.D0
          do k = 1, 6
            sum = sum + vib(k,i)*shift(k)*vib(k,j)
          end do
          fmat(l) = fmat(l) + sum
        end do
      end do
      return
      end subroutine frame
