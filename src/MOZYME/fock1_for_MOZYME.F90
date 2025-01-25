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

subroutine fock1_for_MOZYME (f, ptot, w, kr, iab, ilim)
    implicit none
    integer, intent (in) :: iab, ilim
    integer, intent (inout) :: kr
    double precision, dimension ((iab*(iab+1))/2), intent (in) :: ptot
    double precision, dimension ((iab*(iab+1))/2), intent (inout) :: f
    double precision, dimension (ilim, ilim), intent (in) :: w
!
    integer :: i, ij, ijp, ijw, ikw, im, ip, j, jlw, jm, jp, k, klw, l
    double precision :: sum
   ! *********************************************************************
   !
   ! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTER ELEMENTS.
   !
   ! *********************************************************************
   !
   !   One-center coulomb and exchange terms for atom II.
   !
   !  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
   !                        -(PA(k,l)        )*<i,k|j,l>), k,l on atom II.
   !
    do i = 1, iab
      do j = 1, i
         !
         !    Address in 'F'
         !
        ij = (i*(i-1)) / 2 + j
         !
         !    'J' Address IJ in W
         !
        ijw = (i*(i-1)) / 2 + j
        sum = 0.d0
        do k = 1, iab
          do l = 1, iab
            ip = Max (k, l)
            jp = Min (k, l)
               !
               !    Address in 'P'
               !
            ijp = (ip*(ip-1)) / 2 + jp
               !
               !    'J' Address KL in W
               !
            im = Max (k, l)
            jm = Min (k, l)
            klw = (im*(im-1)) / 2 + jm
               !
               !    'K' Address IK in W
               !
            im = Max (k, j)
            jm = Min (k, j)
            ikw = (im*(im-1)) / 2 + jm
               !
               !    'K' Address JL in W
               !
            im = Max (l, i)
            jm = Min (l, i)
            jlw = (im*(im-1)) / 2 + jm
               !
               !   The term itself
               !
            sum = sum + ptot(ijp) * w(ijw, klw) - 0.5d0 * ptot(ijp) * w &
           & (ikw, jlw)
          end do
        end do
        f(ij) = f(ij) + sum
      end do
    end do
    kr = kr + ilim ** 2
end subroutine fock1_for_MOZYME
