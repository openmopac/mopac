! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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
