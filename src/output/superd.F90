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

      subroutine superd(c, eigs, norbs, nelecs, numat, nat)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw
      use elemts_C, only : elemnt
      use parameters_C, only: natorb
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbs
      integer , intent(in) :: nelecs
      integer , intent(in) :: numat
      integer , intent(in) :: nat(numat)
      double precision , intent(in) :: c(norbs,norbs)
      double precision , intent(in) :: eigs(norbs)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ihomo, isize, i, ii, k, j, l, ll
      double precision, dimension(4) :: orb
      double precision :: alpha, totloc, totnuc, cdens, deloc, denuc, seltot, &
        selpol, sum
!-----------------------------------------------
      ihomo = nelecs/2
      isize = norbs
      alpha = (eigs(ihomo)+eigs(ihomo+1))/2.D0
      write (iw, '(A,F12.6)') ' Mulliken electronegativity:    ', (-alpha)
      write (iw, '(A,F12.6)') ' Parr & Pople absolute hardness:', (eigs(ihomo+1&
        )-eigs(ihomo))/2.D0
      write (iw, '(A,F12.6)') ' Schuurmann MO shift alpha:     ', alpha
      write (iw, *)
      write (iw, '(A,F12.6)') ' Ehomo:                         ', eigs(ihomo)
      write (iw, '(A,F12.6)') ' Elumo:                         ', eigs(ihomo+1)
      write (iw, '(2/,A,/)') '  a   n        Dn(r)        De(r)   q(r) - Z(r)'
      totloc = 0.D0
      totnuc = 0.D0
      i = 1
      do ii = 1, numat
        if (nat(ii) == 1) then
          i = i + 1
        else
!
!   CALCULATE DELOCALIZABILITY
!
          cdens = 0.D0
          deloc = 0.D0
          denuc = 0.D0
          do k = i, i + natorb(nat(ii)) - 1
            do j = 1, ihomo
              cdens = cdens + c(k,j)*c(k,j)
              deloc = deloc + c(k,j)*c(k,j)/(eigs(j)-alpha)
            end do
            do j = ihomo + 1, isize
              denuc = denuc - c(k,j)*c(k,j)/(eigs(j)-alpha)
            end do
          end do
          write (iw, '(1X,A2,1X,I3,3F13.6)') elemnt(nat(ii)), ii, 2.D0*denuc, &
            2.D0*deloc, (-2.D0*cdens)
          i = i + 4
          totnuc = totnuc + 2.D0*denuc
          totloc = totloc + 2.D0*deloc
        end if
      end do
      write (iw, '(/,A,2F13.6)') ' Total:', totnuc, totloc
      write (iw, '(/,A,/)') '  a   n        piS(r)'
      seltot = 0.D0
      i = 1
      do ii = 1, numat
        if (nat(ii) == 1) then
          i = i + 1
        else
          selpol = 0.D0
          do l = i, i + natorb(nat(ii)) - 1
            do j = 1, ihomo
              do k = ihomo + 1, isize
                selpol = selpol + c(l,j)*c(l,j)*c(l,k)*c(l,k)/(eigs(k)-eigs(j))
              end do
            end do
          end do
          selpol = -4.D0*selpol
          write (iw, '(1X,A2,1X,I3,F13.6)') elemnt(nat(ii)), ii, selpol
          seltot = seltot + selpol
          i = i + natorb(nat(ii))
        end if
      end do
      write (iw, '(/,A,F13.6,/)') ' Total:', seltot
      write (iw, '(A,/)') &
        '  a   n       homo-1         homo         lumo       lumo+1'
      i = 1
      do ii = 1, numat
        if (nat(ii) == 1) then
          i = i + 1
        else
          ll = 0
          do j = ihomo - 1, ihomo + 2
            ll = ll + 1
            sum = 0.D0
            do l = i, i + natorb(nat(ii)) - 1
              sum = sum + c(l,j)*c(l,j)
            end do
            orb(ll) = sum*100.D0
          end do
          write (iw, '(1X,A2,1X,I3,4F13.6)') elemnt(nat(ii)), ii, orb
          i = i + 4
        end if
      end do
      return
      end subroutine superd
