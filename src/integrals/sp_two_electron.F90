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
