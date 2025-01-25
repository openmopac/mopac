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

      subroutine dfield()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : tore
      USE molkst_C, only :numat, efield
      USE funcon_C, only : ev, a0, fpc_9
      use common_arrays_C, only : nat, p, dxyz
!***********************************************************************
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
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      double precision, dimension(numat) :: q2
      double precision :: fldcon
!-----------------------------------------------
      call chrge (p, q2)
      q2(:numat) = tore(nat(:numat)) - q2(:numat)
!
!   FLDCON=(h/Ao)*(eV to Kcal/mol)
!
      fldcon = ev/a0*fpc_9
      do i = 1, numat
        dxyz(3*(i-1) + 1) = dxyz(3*(i-1) + 1) + efield(1)*q2(i)*fldcon
        dxyz(3*(i-1) + 2) = dxyz(3*(i-1) + 2) + efield(2)*q2(i)*fldcon
        dxyz(3*(i-1) + 3) = dxyz(3*(i-1) + 3) + efield(3)*q2(i)*fldcon
      end do
      return
      end subroutine dfield
