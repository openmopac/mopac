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

      subroutine set(s1, s2, na, nb, rab, ii)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE overlaps_C, only : isp, ips, sa, sb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: na
      integer , intent(in) :: nb
      integer , intent(in) :: ii
      double precision , intent(in) :: s1
      double precision , intent(in) :: s2
      double precision , intent(in) :: rab
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jcall
      double precision :: alpha, beta
!-----------------------------------------------
!***********************************************************************
!
!     SET IS PART OF THE OVERLAP CALCULATION, CALLED BY OVERLP.
!         IT CALLS AINTGS AND BINTGS
!
!***********************************************************************
      if (na <= nb) then
        isp = 1
        ips = 2
        sa = s1
        sb = s2
      else
        isp = 2
        ips = 1
        sa = s2
        sb = s1
      end if
      j = ii + 2
      if (ii > 3) j = j - 1
      alpha = 0.5D00*rab*(sa + sb)
      beta = 0.5D00*rab*(sb - sa)
      jcall = j - 1
      call aintgs (alpha, jcall)
      call bintgs (beta, jcall)
      return
!
      end subroutine set
      subroutine aintgs(x, k)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE overlaps_C, only : a
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: k
      double precision , intent(in) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision :: c
!-----------------------------------------------
!***********************************************************************
!
!    AINTGS FORMS THE "A" INTEGRALS FOR THE OVERLAP CALCULATION.
!
!***********************************************************************
      c = exp((-x))
      a(1) = c/x
      do i = 1, k
        a(i+1) = (a(i)*i+c)/x
      end do
      return
!
      end subroutine aintgs
      subroutine bintgs(x, k)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE overlaps_C, only : b, fact
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: k
      double precision , intent(in) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: io, last, i, m
      double precision :: absx, expx, expmx, y, xf
!-----------------------------------------------
!**********************************************************************
!
!     BINTGS FORMS THE "B" INTEGRALS FOR THE OVERLAP CALCULATION.
!
!**********************************************************************
      io = 0
      absx = abs(x)
      if (absx > 3.D00) go to 40
      if (absx > 2.D00) then
        if (k <= 10) go to 40
        last = 15
        go to 60
      end if
      if (absx > 1.D00) then
        if (k <= 7) go to 40
        last = 12
        go to 60
      end if
      if (absx > 0.5D00) then
        if (k <= 5) go to 40
        last = 7
        go to 60
      end if
      if (absx <= 1.D-6) go to 90
      last = 6
      go to 60
   40 continue
      expx = exp(x)
      expmx = 1.D00/expx
      b(1) = (expx - expmx)/x
      do i = 1, k
        b(i+1) = (i*b(i)+(-1.D00)**i*expx-expmx)/x
      end do
      go to 110
   60 continue
      do i = io, k
        y = 0.0D00
        do m = io, last
          xf = 1.0D00
          if (m /= 0) xf = fact(m)
          y = y + (-x)**m*(2*mod(m + i + 1,2))/(xf*(m + i + 1))
        end do
        b(i+1) = y
      end do
      go to 110
   90 continue
      do i = io, k
        b(i+1) = (2*mod(i + 1,2))/(i + 1.D0)
      end do
  110 continue
      return
!
      end subroutine bintgs
