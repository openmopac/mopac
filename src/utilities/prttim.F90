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

      subroutine prttim(tleft, tprt, txt)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision , intent(in) :: tleft
      double precision , intent(out) :: tprt
      character , intent(out) :: txt
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      tprt = tleft
      txt = 'S'
      if (tprt >= 604800D0) then
        tprt = tprt/604800.D0
        txt = 'W'
      else if (tprt >= 86400.D0) then
        tprt = tprt/86400.D0
        txt = 'D'
      else if (tprt >= 3600.D0) then
        tprt = tprt/3600.D0
        txt = 'H'
      else if (tprt >= 60.D0) then
        tprt = tprt/60.D0
        txt = 'M'
      end if
      return
      end subroutine prttim
