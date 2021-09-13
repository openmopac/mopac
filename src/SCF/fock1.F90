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

subroutine fock1(f, ptot, pa, mpack, w, kr, ia, ib, ilim)
        implicit none
        integer, intent (in) :: ia, ib, ilim, mpack
        integer, intent (inout) :: kr
        double precision, dimension (mpack), intent (in) :: pa, ptot
        double precision, dimension (mpack), intent (inout) :: f
        double precision, dimension (ilim, ilim), intent (in) :: w
        integer :: i, ij, ijp, ijw, ikw, im, ip, iw, j, jlw, jm, jp, jw, k, klw, &
       & kw, l, lw
        double precision :: sum
        intrinsic Max, Min
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
        do i = ia, ib
          iw = i - ia + 1
          do j = ia, i
            jw = j - ia + 1
             !
             !    Address in 'F'
             !
            ij = (i*(i - 1))/2 + j
             !
             !    'J' Address IJ in W
             !
            ijw = (iw*(iw - 1))/2 + jw
            sum = 0.d0
            do k = ia, ib
              kw = k - ia + 1
              do l = ia, ib
                lw = l - ia + 1
                ip = Max (k, l)
                jp = Min (k, l)
                   !
                   !    Address in 'P'
                   !
                ijp = (ip*(ip - 1))/2 + jp
                   !
                   !    'J' Address KL in W
                   !
                im = Max (kw, lw)
                jm = Min (kw, lw)
                klw = (im*(im - 1))/2 + jm
                   !
                   !    'K' Address IK in W
                   !
                im = Max (kw, jw)
                jm = Min (kw, jw)
                ikw = (im*(im - 1))/2 + jm
                   !
                   !    'K' Address JL in W
                   !
                im = Max (lw, iw)
                jm = Min (lw, iw)
                jlw = (im*(im - 1))/2 + jm
                   !
                   !   The term itself
                   !
                sum = sum + ptot(ijp) * w(ijw, klw) - pa(ijp) * w(ikw, jlw)
              end do
            end do
            f(ij) = f(ij) + sum
          end do
        end do
        kr = kr + ilim**2
end subroutine fock1
