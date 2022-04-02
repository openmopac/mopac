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

      subroutine sp_two_electron
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameters_C, only : zsn, zpn, gss, gsp, gpp, gp2, hsp, &
      main_group
      use mndod_C, only : iii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ni, ns
      double precision :: es, ep, r033, r233
      double precision, external :: rsc
!
      do ni = 1, 80
        ns = iii(ni)
        es = zsn(ni)
        ep = zpn(ni)
        if (es < 1.d-4 .or. ep < 1.d-4 .or. main_group(ni)) cycle
        gss(ni) = rsc(0,ns,es,ns,es,ns,es,ns,es)
        gsp(ni) = rsc(0,ns,es,ns,es,ns,ep,ns,ep)
        hsp(ni) = rsc(1,ns,es,ns,ep,ns,es,ns,ep)/3.D0
        r033 = rsc(0,ns,ep,ns,ep,ns,ep,ns,ep)
        r233 = rsc(2,ns,ep,ns,ep,ns,ep,ns,ep)
        gpp(ni) = r033 + 0.16D0*r233
        gp2(ni) = r033 - 0.08D0*r233
      end do
      return
      end subroutine sp_two_electron
