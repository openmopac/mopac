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

      subroutine xyzcry(tvec, numat, dxyz, iw)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: numat
      integer , intent(in) :: iw
      double precision , intent(inout) :: tvec(3,3)
      double precision , intent(inout) :: dxyz(3,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      double precision :: sum, ca, sa, sum1
!-----------------------------------------------
!
!   Convert Cartesian derivatives into fractional unit cell derivatives.
!
!      WRITE(IW,*)' Atoms 1 and 2'
!      WRITE(IW,'(3F12.6)')((DXYZ(I,J),I=1,3),J=1,2)
      sum = sqrt(tvec(2,1)**2+tvec(3,1)**2)
!      WRITE(IW,*)' TVEC'
!      WRITE(IW,'(3f12.4)')TVEC
      if (sum > 1.D-6) then
!
!    Rotate to eliminate TVEC(3,1)
!
        ca = tvec(2,1)/sum
        sa = tvec(3,1)/sum
        do i = 1, 3
          sum1 = tvec(2,i)*ca + tvec(3,i)*sa
          tvec(3,i) = (-tvec(2,i)*sa) + tvec(3,i)*ca
          tvec(2,i) = sum1
        end do
!      WRITE(IW,*)' TVEC'
!      WRITE(IW,'(3f12.4)')TVEC
        do i = 1, numat
          sum1 = dxyz(2,i)*ca + dxyz(3,i)*sa
          dxyz(3,i) = (-dxyz(2,i)*sa) + dxyz(3,i)*ca
          dxyz(2,i) = sum1
        end do
!
!    Rotate to eliminate TVEC(2,1)
!
        sum = sqrt(tvec(1,1)**2+tvec(2,1)**2)
        ca = tvec(1,1)/sum
        sa = tvec(2,1)/sum
        do i = 1, 3
          sum1 = tvec(1,i)*ca + tvec(2,i)*sa
          tvec(2,i) = (-tvec(1,i)*sa) + tvec(2,i)*ca
          tvec(1,i) = sum1
        end do
        do i = 1, numat
          sum1 = dxyz(1,i)*ca + dxyz(2,i)*sa
          dxyz(2,i) = (-dxyz(1,i)*sa) + dxyz(2,i)*ca
          dxyz(1,i) = sum1
        end do
!      WRITE(IW,*)' TVEC'
!      WRITE(IW,'(3f12.4)')TVEC
      end if
!
!    Rotate to eliminate TVEC(3,2)
!
      sum = sqrt(tvec(2,2)**2+tvec(3,2)**2)
      if (sum > 1.D-6) then
        ca = tvec(2,2)/sum
        sa = tvec(3,2)/sum
        do i = 2, 3
          sum1 = tvec(2,i)*ca + tvec(3,i)*sa
          tvec(3,i) = (-tvec(2,i)*sa) + tvec(3,i)*ca
          tvec(2,i) = sum1
        end do
!      WRITE(IW,*)' TVEC'
!      WRITE(IW,'(3f12.4)')TVEC
        do i = 1, numat
          sum1 = dxyz(2,i)*ca + dxyz(3,i)*sa
          dxyz(3,i) = (-dxyz(2,i)*sa) + dxyz(3,i)*ca
          dxyz(2,i) = sum1
        end do
      end if
!
!  Convert unit cell into it's reciprocal
!
      do i = 1, 3
        sum = 0.D0
        do j = 1, i
          sum = sum + tvec(j,i)**2
        end do
        tvec(:i,i) = tvec(:i,i)/sum
      end do
!      WRITE(IW,*)' TVEC'
!      WRITE(IW,'(3f12.4)')TVEC
!         WRITE(IW,'(A)')' Fractional Unit Cell Derivatives'
!         WRITE(IW,'(I4,3F12.5)')(I,(DXYZ(J,I),J=1,3),I=1,NUMAT)
      dxyz(3,:numat) = dxyz(3,:numat)/tvec(3,3)
      dxyz(2,:numat) = dxyz(2,:numat) - dxyz(3,:numat)*tvec(2,3)
      dxyz(1,:numat) = dxyz(1,:numat) - dxyz(3,:numat)*tvec(1,3)
      dxyz(2,:numat) = dxyz(2,:numat)/tvec(2,2)
      dxyz(1,:numat) = dxyz(1,:numat) - dxyz(2,:numat)*tvec(1,2)
      dxyz(1,:numat) = dxyz(1,:numat)/tvec(1,1)
      write (iw, '(A)') ' Fractional Unit Cell Derivatives'
      write (iw, '(I4,3F12.5)') (i,(dxyz(j,i),j=1,3),i=1,numat)
      return
      end subroutine xyzcry
