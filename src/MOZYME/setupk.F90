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

subroutine setupk (nocc1)
    use MOZYME_C, only : icocc, ncf, nncf, kopt
    use molkst_C, only : numat
    implicit none
    integer, intent (in) :: nocc1
    integer :: i, j, k, l
   !********************************************************************
   !
   !   SETUPK DETERMINES WHICH ATOMS NEED TO BE CONSIDERED IN
   !   RE-CONSTRUCTING THE FOCK MATRIX.
   !
   !   THE SET OF ATOMS TO BE USED IS PUT INTO KOPT.  ALL ATOMS USED IN
   !   THE LMOs IN THE SCF CALCULATION ARE PUT INTO KOPT.
   !
   !   (Note:  Perhaps a better test would be 'All atoms that make a
   !   significant contribution to the LMOs in the SCF calculation?'
   !
   !********************************************************************
    kopt = 0
    do i = 1, nocc1
      j = nncf(i)
      do k = 1, ncf(i)
        kopt(icocc(k+j)) = 1
      end do
    end do
!
!   If atom i is to be used in constructing the Fock matrix, then kopt(i) = 1
!   otherwise kopt(i) = 0
!
!   Now compress kopt, so that all the atoms used for the Fock matrix are in
!   order.
!
    l = 0
    do i = 1, numat
      if (kopt(i) == 1) then
        l = l + 1
        kopt(l) = i
      end if
    end do
!
!  Force a zero in after the last atom.  This can be used in finding the end of the atom list.
!
    if (l /= numat) kopt(l+1) = 0
end subroutine setupk
