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

      subroutine resolv(c, cold, mdim, eig, nocc, size_of_LMO)
      use common_arrays_C, only : nfirst, nlast
      use molkst_C, only : norbs, numat, keywrd
      implicit none
      integer , intent(in) :: mdim
      integer , intent(in) :: nocc
      double precision , intent(inout) :: c(mdim,mdim)
      double precision , intent(in) :: cold(mdim,mdim)
      double precision , intent(in) :: eig(nocc), size_of_LMO(nocc)
!
      integer , dimension(4) :: idegen
      integer :: loop, j, k, nsec, i, ij, ii, jj, l, mo1, mo2, mo3, atom_a, atom_b
      double precision, dimension(10) :: sec
      double precision, dimension(16) :: vec
      double precision, dimension(4) :: eigs
      double precision :: thresh, sum1, sum, coi, coj, sum2, sum3, &
        sqrt_2, sqrt_3, sqrt_23, sqrt_6
      logical :: l_banana, l_rabbit, lone_pair, multiple
      logical, allocatable :: used_MO(:)
!***********************************************************************
!
!   RESOLVE removes any ill-definition in the LMOs.
!
!   If two or more LMOs have the same atomic contributions, then
!   any linear combination of these LMOs is an acceptable solution.
!   This is undesirable, in that the energies of these LMOs is therefore
!   not defined.  RESOLVE will identify such sets of LMOs, and resolve
!   them so that they have a zero energy interaction, that is, the
!   integral <psi(1)|F|psi(2)> is zero.
!
!***********************************************************************
      atom_a = 0
      atom_b = 0
      thresh = 0.05d0
      lone_pair = .false.
      l_banana = (index(keywrd, "BANANA") /= 0)
      l_rabbit = (index(keywrd, "RABBIT") /= 0)
      allocate(used_MO(norbs))
      used_MO = .false.
      sqrt_2 = sqrt(0.5d0)
      sqrt_3 = sqrt(1/3.d0)
      sqrt_23 = sqrt(2/3.d0)
      sqrt_6 = sqrt(1/6.d0)

      l160: do loop = 1, nocc
!
!  Test the LMO to see if it is a potential candidate for degeneracy.
!
!  If the LMO is entirely on one atom (a lone pair) or
!  if the LMO is over 30% on one atom and entirely on two atoms (a double or triple bond)
!
        if (used_MO(loop)) cycle
        do j = 1, numat
          if (nlast(j) == nfirst(j)) cycle
          sum1 = 0.D0
          do k = nfirst(j), nlast(j)
            sum1 = sum1 + c(k,loop)**2
          end do
!
!  Only LMOs that are near 100% on an atom (a lone pair) or at least 30% on one atom and entirely
!  on two atoms (a bond) are potential candidates.
!
          if (sum1 > 1.d0 - thresh) then
            lone_pair = .true.
            atom_a = j
            exit
          else if (sum1 > 0.3d0) then
            do l = j + 1, numat
              sum2 = 0.D0
              do k = nfirst(l), nlast(l)
                sum2 = sum2 + c(k,loop)**2
              end do
              if (sum1 + sum2 > 1.d0 - thresh) exit
            end do
            if (l <= numat) then
              atom_a = j
              atom_b = l
              exit
            end if
          end if
        end do
        if (j > numat) then
!
! No atoms qualify as belonging to either a lone pair or to an atom in a diatomic bond
! so go on to the next LMO
!
          continue
          cycle  l160
        end if
!
! If a lone pair, then "atom_a" is the atom involved
! If a bond, then "atom_a" is the first atom, and "atom_b" is the second atom in the bond
!
        nsec = 1
        idegen(nsec) = loop
!
!   LMO 'LOOP' is a candidate.  Now identify any related LMOs
!
        do i = loop + 1, nocc
!
! Check that LMO "i" has not already been used.
!
          if (used_MO(i)) cycle
!
! Check that LMO "i" has roughly the same number of centers as LMO "loop".
! If it does not, then it is not a candidate for rabbit ears or banana bonds
!
          if (abs(size_of_LMO(i) - size_of_LMO(loop)) > 0.3d0) cycle
          multiple = .false.
          sum1 = 0.d0
          do k = nfirst(atom_a), nlast(atom_a)
              sum1 = sum1 + c(k,i)**2
          end do
          if (lone_pair) then
            multiple = (sum1 > 1.d0 - thresh)
          else
            sum2 = 0.d0
            do k = nfirst(atom_b), nlast(atom_b)
                sum2 = sum2 + c(k,i)**2
            end do
            multiple = (sum1 + sum2 > 1.d0 - thresh)
          end if
          if (multiple) then
            nsec = nsec + 1
            idegen(nsec) = i
            used_MO(i) = .true.
          end if
        end do
        if (nsec /= 1 .and. nsec < 4) then
!
!   Build small secular determinant.
!
          ij = 0
          do ii = 1, nsec
            i = idegen(ii)
            do jj = 1, ii
              j = idegen(jj)
              sum = 0.D0
              do l = 1, nocc
                coi = 0.D0
                coj = 0.D0
                do k = 1, norbs
                  coi = coi + cold(k,l)*c(k,i)
                  coj = coj + cold(k,l)*c(k,j)
                end do
                sum = sum + coi*eig(l)*coj
              end do
              ij = ij + 1
              sec(ij) = sum
            end do
          end do
!
!   Diagonalize, to identify LCMO
!
          call rsp (sec, nsec, eigs, vec)
          sum = eigs(1) ! dummy use of eigs
!
!    Crude, but fast, way of rotating LMOs
!
          select case (nsec)
          case default
            mo1 = idegen(1)
            mo2 = idegen(2)
            do i = 1, norbs
              sum1 = vec(1)*c(i,mo1) + vec(2)*c(i,mo2)
              sum2 = vec(3)*c(i,mo1) + vec(4)*c(i,mo2)
              c(i,mo1) = sum1
              c(i,mo2) = sum2
            end do
            if (lone_pair .and. l_rabbit .or. .not. lone_pair .and. l_banana) then
              sum3 = sqrt(0.5d0)
              do i = 1, norbs
                sum1 = sum3*c(i,mo1) + sum3*c(i,mo2)
                sum2 = sum3*c(i,mo1) - sum3*c(i,mo2)
                c(i,mo1) = sum1
                c(i,mo2) = sum2
              end do
            end if
          case (3)
            mo1 = idegen(1)
            mo2 = idegen(2)
            mo3 = idegen(3)
            do i = 1, norbs
              sum1 = vec(1)*c(i,mo1) + vec(2)*c(i,mo2) + vec(3)*c(i,mo3)
              sum2 = vec(4)*c(i,mo1) + vec(5)*c(i,mo2) + vec(6)*c(i,mo3)
              sum3 = vec(7)*c(i,mo1) + vec(8)*c(i,mo2) + vec(9)*c(i,mo3)
              c(i,mo1) = sum1
              c(i,mo2) = sum2
              c(i,mo3) = sum3
            end do
            if (lone_pair .and. l_rabbit .or. .not. lone_pair .and. l_banana) then
              do i = 1, norbs
                sum1 = sqrt_3*c(i,mo1) - sqrt_2*c(i,mo2) +  sqrt_6*c(i,mo3)
                sum2 = sqrt_3*c(i,mo1) + sqrt_2*c(i,mo2) +  sqrt_6*c(i,mo3)
                sum3 = sqrt_3*c(i,mo1)                   - sqrt_23*c(i,mo3)
                c(i,mo1) = sum1
                c(i,mo2) = sum2
                c(i,mo3) = sum3
              end do
              continue
            end if
          end select
        end if
      end do l160
      return
      end subroutine resolv
