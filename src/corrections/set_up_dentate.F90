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

  subroutine set_up_dentate()
!
!   Work out which atoms are connected to which other atoms.
!   Atoms are assumed to be connected if they are within a certain distance of each other.
!
!  On exit:  nbonds(i) = number of atoms attached to atom "i"
!            ibonds(j,i) = atom numbers of atoms attached to atom i (there are nbonds(i) of these)
    use molkst_C, only : numat, l11, l21, l31, id, pdb_label
    use mozyme_C, only : nijbo, tyres
    use common_arrays_C, only : nat, coord, tvec, nbonds, ibonds, txtatm
    use atomradii_C, only: atom_radius_covalent, radius
    implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    integer :: i, ik, j, k, jk, kl, il, iu, l
    double precision :: rmin, safety
   !
   !.. Local Arrays ..
    double precision, dimension (3) :: coord1
    nbonds = 0
    ibonds = 0
    if (allocated(radius)) deallocate (radius)
    allocate (radius(numat))
    call extvdw_for_MOZYME (radius, atom_radius_covalent)
   !
   !  ATOMS ARE ASSUMED ATTACHED IF THEY ARE WITHIN
   !  1.1  TIMES THE SUM OF THEIR COVALENT RADII.
   !
    do i = 1, numat
      do j = 1, i - 1
        rmin = 1.d6
        if (id == 0) then
          rmin = (coord(1, i)-coord(1, j)) ** 2 + &
                 (coord(2, i)-coord(2, j)) ** 2 + &
                 (coord(3, i)-coord(3, j)) ** 2
        else
          do ik = -l11, l11
            do jk = -l21, l21
              do kl = -l31, l31
                coord1(:) = coord(:, j) + tvec(:, 1)*ik + tvec(:, 2)*jk + tvec(:, 3)*kl
                rmin = Min (rmin, (coord1(1)-coord(1, i))**2 &
                                + (coord1(2)-coord(2, i))**2 &
                                + (coord1(3)-coord(3, i))**2)
              end do
            end do
          end do
        end if
!
!  Apply safety criteria to specific diatomic pairs
!
        il = min(nat(i), nat(j))
        iu = max(nat(i), nat(j))
        safety = 1.1d0
        select case (il)
          case (1)
            if (iu == 6) safety = 1.25d0  ! C-H
          case (5)
            if (iu ==7) safety = 1.0d0   ! B-N bonds are unusually short
          case (6)
            if (iu < 8) safety = 1.2d0   ! C-C, C-N, C-O
          case (16)
            if (iu == 16) safety = 1.2d0 ! S-S
          end select
        if (pdb_label) then
          if (il /= 1 .and. iu /= 1 .and. txtatm(i)(18:27) == txtatm(j)(18:27)) then
            do k = 1, 20
              if (txtatm(i)(18:20) == tyres(k)) exit
            end do
            if (k < 21) then
!
!  Both atoms are heavy atoms in the same standard amino acid residue, so be more tolerant.
!
              safety = safety*1.05d0
            end if
          end if
        end if
        if (rmin < (safety*(radius(i) + radius(j)))**2) then
          if (nbonds(i) < 15 .and. nbonds(j) < 15) then
            nbonds(i) = nbonds(i) + 1
            nbonds(j) = nbonds(j) + 1
            ibonds(nbonds(i), i) = j
            ibonds(nbonds(j), j) = i
          end if
        end if
      end do
    end do
!
!  Check for H attached to H
!
    do i = 1, numat
      if (nat(i) /= 1) cycle
      if (nbonds(i) < 2) cycle
      k = 0
      do j = 1, nbonds (i)
        if (nat(ibonds(j,i)) /= 1) then
          k = k +1
          ibonds(k,i) = ibonds(j,i)
        end if
      end do
      nbonds(i) = k
    end do
    if (allocated(nijbo)) then
!
!  Sanity check - don't allow a bond if the atoms are too far apart for nijbo to be positive
!
      do i = 1, numat
        k = 0
        do j = 1, nbonds(i)
          l = ibonds(j,i)
          if (nijbo(l,i) > -1) then
            k = k + 1
            ibonds(k,i) = l
          end if
          nbonds(i) = k
        end do
      end do
    end if
    return
  end subroutine set_up_dentate
!
!
!
  double precision function nsp2_correction()
!
!   Add a molecular mechanics correction to all nitrogen atoms that have exactly
!   three ligands.
!
    use common_arrays_C, only : nat, coord, nbonds, ibonds
    use molkst_C, only : numat, method_pm6, method_pm7, method_pm6_org, method_pm8
    implicit none
    integer :: i, j
    double precision :: correction, sum
    double precision, external :: nsp2_atom_correction
    if ( .not. (method_pm6 .or. method_pm7 .or. method_pm6_org .or. method_pm8)) then
      nsp2_correction = 0.d0
      return
    end if
    correction = 0.d0
    do i = 1, numat
      if (nat(i) == 7 .and. nbonds(i) == 3) then
        j = 0
        if (nat(ibonds(1,i)) == 1) j = 1
        if (nat(ibonds(2,i)) == 1) j = j + 1
        if (nat(ibonds(3,i)) == 1) j = j + 1
        if ( j < 2) then
          sum = nsp2_atom_correction(coord, i, ibonds(1,i), ibonds(2,i), ibonds(3,i))
          correction = correction + sum
         end if
      end if
    end do
    nsp2_correction = correction
  end function nsp2_correction
!
!
!
  double precision function nsp2_atom_correction(vectors,n,i,j,k)
  implicit none
  double precision :: vectors(3,*)
  integer :: n,i,j,k

  double precision :: a, b, c, ab, ac, bc, tot, cosa, cosb, cosc
!
!  Evaluate the penalty for non-planarity
!  - done by working out the three angles about the central atom
!  (here "n") subtended by the lines to atoms "i", "j", and "k".
!
  a = sqrt((vectors(1,n) - vectors(1,i))**2 + &
           (vectors(2,n) - vectors(2,i))**2 + &
           (vectors(3,n) - vectors(3,i))**2)
  b = sqrt((vectors(1,n) - vectors(1,j))**2 + &
           (vectors(2,n) - vectors(2,j))**2 + &
           (vectors(3,n) - vectors(3,j))**2)
  c = sqrt((vectors(1,n) - vectors(1,k))**2 + &
           (vectors(2,n) - vectors(2,k))**2 + &
           (vectors(3,n) - vectors(3,k))**2)
  ab = sqrt((vectors(1,j) - vectors(1,i))**2 + &
            (vectors(2,j) - vectors(2,i))**2 + &
            (vectors(3,j) - vectors(3,i))**2)
  ac = sqrt((vectors(1,k) - vectors(1,i))**2 + &
            (vectors(2,k) - vectors(2,i))**2 + &
            (vectors(3,k) - vectors(3,i))**2)
  bc = sqrt((vectors(1,j) - vectors(1,k))**2 + &
            (vectors(2,j) - vectors(2,k))**2 + &
            (vectors(3,j) - vectors(3,k))**2)
!
  cosa = dacos((b**2 +c**2 -bc**2)/(2.d0*b*c))
  cosb = dacos((a**2 +c**2 -ac**2)/(2.d0*a*c))
  cosc = dacos((b**2 +a**2 -ab**2)/(2.d0*b*a))
!
!  tot = difference between the sum of the three angles and
!        360 degrees, expressed as radians.
  tot = 4.d0*asin(1.d0) - (cosa + cosb + cosc)
  nsp2_atom_correction = -0.5d0*exp(-10.d0*tot)
  end function nsp2_atom_correction
!
!
!
  double precision function C_triple_bond_C()
!
! Evaluate energy contribution from acetylenic bonds - this is
! a correction to account for the extra stabilization of yne bonds
!
  use common_arrays_C, only : nat, coord, nbonds, ibonds
  use molkst_C, only : numat, method_pm6, method_pm7, method_pm6_org, method_pm8
  implicit none
  integer :: i, j, k, isum
  double precision :: rab
    if ( .not. (method_pm6 .or. method_pm7 .or. method_pm6_org .or. method_pm8)) then
      C_triple_bond_C = 0.d0
      return
    end if
    isum = 0
    do i = 1, numat
      rab = 10.d0
      if (nat(i) == 6 .and. nbonds(i) == 2) then
!
!  Possible triple bond
!
        do k = 1, nbonds(i)
          j = ibonds(k,i)
          if (j > i) cycle
          if(nat(j) == 6 .and. nbonds(j) == 2) then
            rab = (coord(1,i) - coord(1,j))**2 + (coord(2,i) - coord(2,j))**2 + (coord(3,i) - coord(3,j))**2
!
! C-C double bond length = 1.34 Angstroms
! C-C triple bond length = 1.20 Angstroms
!
            if (rab < 1.27d0**2) goto 99
          end if
        end do
        if (i < j) cycle
!
!   Carbon atom "i" is attached to two other atoms, and is attached to carbon atom "j"
!   and the i-j distance indicates that the bond is acetylenic.
!
      end if
      cycle
99    isum = isum + 1
    end do
    C_triple_bond_C = isum*12.d0  !  (The value "12" was determined empirically
  end function C_triple_bond_C
!
!
!
   double precision function Si_O_H_Correction()
!
!   If an Si-O-H structure, add in a bending perturbation
!
      use common_arrays_C, only : nat, nbonds, ibonds, coord
!
      use molkst_C, only : numat
!
      implicit none
      integer :: i, j, k, Si, O, H
      double precision :: sum
      double precision, external :: Si_O_H_bond_correction
!
      sum = 0.d0
      do i = 1, numat
        if (nat(i) == 8) then
          O = i
          Si = 0
          H = 0
          do j = 1, nbonds(i)
            k = ibonds(j,i)
            if (nat(k) == 14) Si = k
            if (nat(k) == 1) H = k
          end do
          if (Si /= 0 .and. H /= 0) then
          sum = sum + Si_O_H_bond_correction(coord, Si, O, H)
          end if
        end if
      end do
      Si_O_H_Correction = sum
      return
    end function Si_O_H_Correction
!
!
!
    double precision function Si_O_H_bond_correction(coord, Si, O, H)
      use funcon_C, only : pi
      implicit none
      integer :: Si, O, H
      double precision :: coord(3,*)
      double precision :: r_Si_O, r_O_H, ref_angle, angle
!
!  Assuming that a Si - O - H structure has been identified, now apply a correction.
!  The angle should be 115 degrees
!  Apply a Gaussian when Si - O > 1.7 Angstroms and when O - H > 1.0 Angstroms
!  Gaussian drops to 0.05 when atoms are no longer considered connected
!  = 2.02 A for Si -O, and 1.21 A for O - H
!
        r_Si_O = (coord(1,O) - coord(1,Si))**2 + (coord(2,O) - coord(2,Si))**2 + (coord(3,O) - coord(3,Si))**2 - &
              1.7d0**2
        r_O_H  = (coord(1,O) - coord(1,H ))**2 + (coord(2,O) - coord(2,H ))**2 + (coord(3,O) - coord(3,H ))**2 - &
              1.d0
        ref_angle = 125.d0*pi/180.d0
        call bangle (coord, Si, O, H, angle)
        Si_O_H_bond_correction = 15.d0*(angle - ref_angle)**2*exp(-33.d0*max(0.d0, r_Si_O))*exp(-68.d0*max(0.d0, r_O_H ))
        return
      end function Si_O_H_bond_correction
