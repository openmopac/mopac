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

      subroutine intfc(fmatrx, xparam, georef, nar, nbr, ncr)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use molkst_C, only : natoms, nvar, numat, l123, l1u, l2u, l3u
      use common_arrays_C, only : na, nb, nc, geo, loc, labels, nat
      use chanel_C, only : iw
      use elemts_C, only : elemnt
      use to_screen_C, only : fcint
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nar(natoms)
      integer , intent(in) :: nbr(natoms)
      integer , intent(in) :: ncr(natoms)
      double precision , intent(in) :: fmatrx((3*numat*(3*numat + 1))/2)
      double precision  :: xparam(nvar)
      double precision , intent(in) :: georef(3,natoms)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, ilim, n3, k, ncol
      double precision :: step, stepi, sum, dumy(3*numat*l123)
      logical :: precis, lxyz
!-----------------------------------------------
      nvar = 0
      numat = 0
      l1u = 0
      l2u = 0
      l3u = 0
      j = 0
      do i = 1, natoms
      if (nar(i) /= 0) j = j + 1
      end do
      lxyz = (j == 0) ! All coordinates are Cartesian
      na(:natoms) = nar
      nb(:natoms) = nbr
      nc(:natoms) = ncr
      if (allocated(fcint)) deallocate (fcint)
      allocate(fcint(3,natoms))
      fcint(:,:natoms) = 0.D0
      geo(:,:natoms) = georef
      do i = 1, natoms
        if (labels(i) == 99) cycle
        numat = numat + 1
        if (lxyz) then
          ilim = 3
        else
          ilim = min(3,i - 1)
        end if
        do j = 1, ilim
          nvar = nvar + 1
          loc(1,nvar) = i
          loc(2,nvar) = j
          xparam(nvar) = geo(j,i)
        end do
      end do
      n3 = 3*numat
!
!   Calculate force constants over internal coordinates
!
      step = 1.D-7
      stepi = 0.5D0/step
      precis = .TRUE.
      if (lxyz) then
        write (iw, '(/,10X,A,/)') &
          ' FORCE CONSTANT IN CARTESIAN COORDINATES (Millidynes/A)'
        write (iw, '(A)') &
          '    ATOM   CHEMICAL       X               Y               Z'
        write (iw, '(A)') &
          '   NUMBER   SYMBOL  FORCE CONSTANT  FORCE CONSTANT  FORCE CONSTANT'
      else
        write (iw, '(/,10X,A,/)') &
          ' FORCE CONSTANT IN INTERNAL COORDINATES (Millidynes/A)'
        write (iw, '(A)') &
          '    ATOM   CHEMICAL  BOND LENGTH      BOND ANGLE     TWIST ANGLE'
        write (iw, '(A)') &
          '   NUMBER   SYMBOL  FORCE CONSTANT  FORCE CONSTANT  FORCE CONSTANT'
      end if
      write (iw, '(A)')
      l123 = 1
      do i = 1, nvar
        j = i
        call jcarin (xparam, step, precis, dumy, ncol, i, j)
        dumy(:n3) = dumy(:n3)*stepi
!
!   Calculate internal force constant
!
        sum = 0.D0
        do j = 1, n3
          do k = 1, n3
            if (j >= k) then
              sum = sum + dumy(j)*fmatrx((j*(j-1))/2+k)*dumy(k)
            else
              sum = sum + dumy(j)*fmatrx((k*(k-1))/2+j)*dumy(k)
            end if
          end do
        end do
!
!   1.D-5 is to correct the multiplier in FREQCY of the force constants
!
        fcint(loc(2,i),loc(1,i)) = sum*1.d-5
      end do
      if (ncol < 0) return ! dummy statement, just to use ncol
      write (iw, '(I7,6X,A2,3F16.6)') (i,elemnt(nat(i)),(fcint(j,i),j=1,3),i=1,&
        numat)
      call to_screen("To_file: Internal Force Constants")
      deallocate(fcint)
      return
      end subroutine intfc
