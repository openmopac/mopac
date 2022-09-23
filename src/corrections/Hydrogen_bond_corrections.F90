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

double precision function Hydrogen_bond_corrections(l_grad, prt)
!
!    Add in hydrogen-bond corrections
!
  use molkst_C, only : numat, E_hb, N_Hbonds, method_pm7, l123, l1u, l2u, l3u, &
    method_pm6_dh_plus, method_pm6_dh2, method_pm6_d3h4, method_pm6_d3h4x, method_pm6_dh2x, &
    line, id, numcal, method_pm6_org
  use parameters_C, only : tore
  use common_arrays_C, only: coord, nat, q, p, hblist, dxyz, cell_ijk
  implicit none
  logical, intent (in) :: l_grad, prt
!
! Local variables
!
  integer :: i, j, k, l, D, H, A, nd_list, d_list(20), d_l(20), i_cell, ii, iii, kkkk, &
    nrpairs, max_h_bonds, icalcn = -1
  logical :: n_h_bonds, first = .true.
  integer, allocatable :: nrbondsa(:), nrbondsb(:)
  double precision :: EC, ER, vector(numat), sum, sum1, delta = 1.d-5, covrad(94)
  double precision, external :: EC_plus_ER, EH_plus, h_bonds4
  logical, external :: connected

  save
   data covrad /&
  &0.32d0,  0.46d0,  1.20d0,  0.94d0,  0.77d0,  0.75d0,  0.71d0,  0.63d0,  0.64d0,  0.67d0, &
  &1.40d0,  1.25d0,  1.13d0,  1.04d0,  1.10d0,  1.02d0,  0.99d0,  0.96d0,  1.76d0,  1.54d0, &
  &1.33d0,  1.22d0,  1.21d0,  1.10d0,  1.07d0,  1.04d0,  1.00d0,  0.99d0,  1.01d0,  1.09d0, &
  &1.12d0,  1.09d0,  1.15d0,  1.10d0,  1.14d0,  1.17d0,  1.89d0,  1.67d0,  1.47d0,  1.39d0, &
  &1.32d0,  1.24d0,  1.15d0,  1.13d0,  1.13d0,  1.08d0,  1.15d0,  1.23d0,  1.28d0,  1.26d0, &
  &1.26d0,  1.23d0,  1.32d0,  1.31d0,  2.09d0,  1.76d0,  1.62d0,  1.47d0,  1.58d0,  1.57d0, &
  &1.56d0,  1.55d0,  1.51d0,  1.52d0,  1.51d0,  1.50d0,  1.49d0,  1.49d0,  1.48d0,  1.53d0, &
  &1.46d0,  1.37d0,  1.31d0,  1.23d0,  1.18d0,  1.16d0,  1.11d0,  1.12d0,  1.13d0,  1.32d0, &
  &1.30d0,  1.30d0,  1.36d0,  1.31d0,  1.38d0,  1.42d0,  2.01d0,  1.81d0,  1.67d0,  1.58d0, &
  &1.52d0,  1.53d0,  1.54d0,  1.55d0 /
!
  if (first) then
    covrad = 4.d0/3.d0*covrad !  This must be done once per entire run
    first = .false.
  end if
  if (icalcn /= numcal) then
    max_h_bonds = 0
    do j = 1, numat
      if (nat(j) == 8 .or. nat(j) == 7) max_h_bonds = max_h_bonds + 1
    end do
    if (id == 0) then
      max_h_bonds = max_h_bonds*125
    else
      max_h_bonds = max_h_bonds*250
    end if
  end if
  if (max_h_bonds == 0) then
    Hydrogen_bond_corrections = 0.d0
    return
  end if
  if (allocated(nrbondsa))   deallocate(nrbondsa)
  if (allocated(nrbondsb))   deallocate(nrbondsb)
  if (allocated(hblist))     deallocate(hblist)
  allocate (hblist(max_h_bonds,10), nrbondsa(max_h_bonds), nrbondsb(max_h_bonds), stat=i)
    if (i /= 0) then
    line = " Cannot allocate arrays for PM6-DH+"
    call to_screen(trim(line))
    call mopend(trim(line))
    Hydrogen_bond_corrections = 0.d0
    return
  end if
  hblist(:,:) = 0
  nrpairs = 0
  if (method_pm6_dh_plus .or. method_PM7) then
    call all_h_bonds(hblist(1,1), hblist(1,9), hblist(1,5), max_h_bonds, nrpairs)
    call setup_DH_Plus(nrpairs, nrbondsa, nrbondsb, n_h_bonds, covrad)
  else
    call all_h_bonds(hblist(1,1), hblist(1,2), hblist(1,3), max_h_bonds, nrpairs)
  end if
  if (method_pm6_dh2 .or. method_pm6_dh2x) then
    call chrge (p, vector)  ! PM6-DH2 needs partial charges
    do i = 1, numat
      j = nat(i)
      q(i) = tore(j) - vector(i)
    end do
  end if
!
!  Calculate Hydrogen bond energy correction
!
!    The hydrogen-bond acceptor types are as follows:
!
!    (1) nitrogen with no hydrogens bonded to it (mostly in aromatic rings) interacting with any hydrogen,
!    (2) nitrogen with one hydrogen (secondary amines) interacting with any hydrogen,
!    (3) nitrogen with two or more hydrogens (primary amines, ammonia) interacting with any hydrogen,
!    (4) oxygen except carbonyl interacting with HN,
!    (5) carbonyl oxygen interacting with HN,
!    (6) oxygen interacting with HO hydrogen different from 7 and 8,
!    (7) oxygen interacting with H in a water molecule, and
!    (8) oxygen interacting with H in a carboxyl group.
!

  N_Hbonds = 0
  E_hb = 0.d0
  do ii = 1, nrpairs
    if (method_pm6_dh_plus .or. method_pm7) then
      H = hblist(ii,9)
      A = hblist(ii,1)
      D = hblist(ii,5)
      sum = EH_plus(ii, hblist, max_h_bonds, nrbondsa, nrbondsb)
      E_hb = E_hb + sum
      nd_list=0
      do i = 1, 9
        if (hblist(ii,i) > 0) then
          j = hblist(ii,i)
          do l = 1, nd_list
            if (d_list(l) == j) exit
          end do
          if (l > nd_list) then
            nd_list = nd_list + 1
            d_list(nd_list) = hblist(ii,i)
          end if
        end if
      end do
    else
      H = hblist(ii,2)
      A = hblist(ii,3)
      D = hblist(ii,1)
      if (method_pm6_d3h4 .or. method_pm6_d3h4x .or. method_pm6_org) then
        sum = H_bonds4(H, A, D, d_list, nd_list)
        E_hb = E_hb + sum
      else
        sum = EC_plus_ER(hblist(ii,1), hblist(ii,2), hblist(ii,3), q(hblist(ii,2)), q(hblist(ii,3)), EC, ER, d_list, nd_list)
        E_hb = E_hb + EC + ER
      end if

    end if
    if (sum < -1.d0) N_Hbonds = N_Hbonds + 1
    if (prt) call prt_hbonds(D, H, A, sum)
!
! Add in gradient contribution, if requested
!
    if (l_grad .and. sum < -0.01d0) then
      do j = 1, nd_list
        k = d_list(j)
        iii = l123*(k - 1)
        if (connected(H, k, 8.d0**2)) then
!
!   kkkk is the cell that atom k is in, relative to atom H
!
          kkkk = (l3u - cell_ijk(3)) + (2*l3u + 1)*(l2u - cell_ijk(2) + (2*l2u + 1)*(l1u - cell_ijk(1)))
          i_cell = iii + kkkk
!
!  Delta is 10^(-5).  During the evaluation of simple functions such as sqrt and acos, half of the precision is lost
!  This means that, instead of 16 digits of precision, there are only 8.  Tests with delta = 10^(-8) showed severe errors
!  in the derivatives.  The "ideal" value of delta should be 10^(-4), but 10^(-5) was selected by JJPS as the best compromise,
!  based on finding the middle of the plateau of constant derivatives.
!
          do i = 1, 3
            coord(i,k) = coord(i,k) + delta
            if (method_pm6_dh_plus .or. method_pm7) then
              sum1 = EH_plus(ii, hblist, max_h_bonds, nrbondsa, nrbondsb)
            else if (method_pm6_d3h4 .or. method_pm6_d3h4x .or. method_pm6_org) then
              sum1 = H_bonds4(H, A, D, d_l, l)
            else
              sum1 = EC_plus_ER(D, H, A, q(H), q(A), EC, ER, d_l, l)
            end if
            if (abs(sum1) > 1.d-5) then
              sum1 = (sum1 - sum)/delta
              dxyz(i_cell*3 + i) = dxyz(i_cell*3 + i) + sum1
            end if
            coord(i,k) = coord(i,k) - delta
          end do
        end if
      end do
    end if
  end do
  Hydrogen_bond_corrections = E_hb
  return
  end function Hydrogen_bond_corrections
  subroutine setup_DH_Plus(nrpairs, nrbondsa, nrbondsb, l_h_bonds, covrad)
  use molkst_C, only : numat
  use common_arrays_C, only: hblist
  implicit none
  integer :: nrpairs, nrbondsa(nrpairs), nrbondsb(nrpairs)
  logical :: l_h_bonds
  double precision covrad(94)
!
!   Local variables
!
  integer :: i, j, i1, k, l, bondlist(10,3), nrbondsc
  double precision :: xa_dist, xb_dist, xh_dist, xc_dist, sum, sum1, old_dist
  double precision, external :: distance, bonding
  logical :: hbs1_ok, hbs2_ok
  hbs1_ok = .false.
  hbs2_ok = .false.
!
!  Identify atoms associated with the hydrogen bonds
!
  do i = 1, nrpairs
    nrbondsa(i) = 0
    nrbondsb(i) = 0
    bondlist(:,  :) = 0
    do j = 1, numat
      xa_dist = distance(j, hblist(i, 1))
      xb_dist = distance(j, hblist(i, 5))
      if (xa_dist < bonding(j, hblist(i, 1), covrad) .and. hblist(i, 1) /= j) then
!
!  Atom j is covalently bonded to atom hblist(i,1) (either the H bond donor or accepter)
!
        nrbondsa(i) = nrbondsa(i) + 1
        bondlist(nrbondsa(i), 1) = j
        if (nrbondsa(i) == 5) then
!
!  Too many bonds - eliminate the longest bond
!
          i1 = hblist(i, 1)
          sum = 0.d0
          l = 0
          do k = 1, 5
            sum1 = distance(i1, bondlist(k,1))
            if (sum < sum1) then
              sum = sum1
              l = k
            end if
          end do
          do k = l, 4
            bondlist(k,1) = bondlist(k + 1,1)
          end do
          nrbondsa(i) = 4
        end if
      end if
      if (xb_dist < bonding(j, hblist(i, 5), covrad) .and. hblist(i, 5) /= j) then
!
!  Atom j is covalently bonded to atom hblist(i,5) (either the H bond donor or accepter)
!
        nrbondsb(i) = nrbondsb(i) + 1
        bondlist(nrbondsb(i), 2) = j
        if (nrbondsb(i) == 5) then
!
!  Too many bonds - eliminate the longest bond
!
          i1 = hblist(i, 5)
          sum = 0.d0
          l = 0
          do k = 1, 5
            sum1 = distance(i1, bondlist(k,2))
            if (sum < sum1) then
              sum = sum1
              l = k
            end if
          end do
          do k = l, 4
            bondlist(k,2) = bondlist(k + 1,2)
          end do
          nrbondsb(i) = 4
        end if
      end if
    end do
! check for 1-3 and 1-4 case - very seldom needed,  but ...
!
!  Is an atom, attached to hblist(i,1) equal to atom hblist(i,5)
!  That is, are the donor and accepter atoms the same.
!
    do j = 1, nrbondsa(i)
      if (bondlist(j, 1) == hblist(i, 5)) hblist(i, 10) = -666     !  1-3
      do k = 1, nrbondsb(i)
!
!  Is an atom, attached to hblist(i,1) equal to atom attached to hblist(i,5)
!  and not the hydrogen atom
!  That is, are both donor and accepter atoms attached to the same atom, other than the hydrogen
!
        if (bondlist(j, 1) == bondlist(k, 2) .and. bondlist(k, 2) /= hblist(i,9)) hblist(i, 10) = -666 !  1-4
      end do
    end do
    if (hblist(i, 10) == -666) cycle
! hbs1 part
    hbs1_ok = .true.
    if (nrbondsa(i) == 3 .or. nrbondsa(i) == 4) then
      old_dist = -1
      do k = 1, nrbondsa(i)
        xh_dist = distance(bondlist(k, 1), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 2) = bondlist(k, 1)
        end if
      end do
      old_dist = -1
      do k = 1, nrbondsa(i)
        xh_dist = distance(bondlist(k, 1), hblist(i, 9))
        if (xh_dist > old_dist .and. bondlist(k, 1) /= hblist(i, 2)) then
          old_dist = xh_dist
          hblist(i, 3) = bondlist(k, 1)
        end if
      end do
      old_dist = -1
      do k = 1, nrbondsa(i)
        xh_dist = distance(bondlist(k, 1), hblist(i, 9))
        if (xh_dist > old_dist .and. bondlist(k, 1) /= hblist(i, 2) .and. bondlist(k, 1) /= hblist(i, 3)) then
          old_dist = xh_dist
          hblist(i, 4) = bondlist(k, 1)
        end if
      end do    ! possible forth atom should be h
    else if (nrbondsa(i) == 2) then
      old_dist = -1
      do k = 1, nrbondsa(i)
        xh_dist = distance(bondlist(k, 1), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 2) = bondlist(k, 1)
        end if
      end do
      do k = 1, nrbondsa(i)
        if (bondlist(k, 1) /= hblist(i, 2)) then
          hblist(i, 3) = bondlist(k, 1)
        end if
      end do
      if (distance(hblist(i, 1), hblist(i, 9)) < bonding(hblist(i, 1), hblist(i, 9), covrad)) then ! not needed
          hblist(i, 4) = hblist(i, 9)
      else
          hblist(i, 4) = hblist(i, 1)
      end if
    else if (nrbondsa(i) == 1) then
      hblist(i, 2) = bondlist(1, 1)
! need bonds on first bonded atom here
      nrbondsc = 0
      do k = 1, numat
        xc_dist = distance(k, hblist(i, 2))
        if (xc_dist < bonding(k, hblist(i, 2), covrad) .and. hblist(i, 2) /= k) then
          nrbondsc = nrbondsc + 1
          bondlist(nrbondsc, 3) = k
        end if
      end do
      old_dist = -1
      do k = 1, nrbondsc
        xh_dist = distance(bondlist(k, 3), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 3) = bondlist(k, 3)
        end if
      end do
      if (distance(hblist(i, 1), hblist(i, 9)) < bonding(hblist(i, 1), hblist(i, 9), covrad)) then ! not needed
          hblist(i, 4) = hblist(i, 9)
      else
          hblist(i, 4) = hblist(i, 1)
      end if
    else if (nrbondsa(i) == 0) then
      if (distance(hblist(i, 1), hblist(i, 9)) < bonding(hblist(i, 1), hblist(i, 9), covrad)) then ! not needed
        hblist(i, 2) = hblist(i, 9)
        hblist(i, 3) = hblist(i, 9)
        hblist(i, 4) = hblist(i, 9)
      else
        hblist(i, 2) = hblist(i, 1)
        hblist(i, 3) = hblist(i, 1)
        hblist(i, 4) = hblist(i, 1)
      end if
    else
      hbs1_ok = .false.
    end if
! hbs2 part
    hbs2_ok = .true.
    if (nrbondsb(i) == 3 .or. nrbondsb(i) == 4) then
      old_dist = -1
      do k = 1, nrbondsb(i)
        xh_dist = distance(bondlist(k, 2), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 6) = bondlist(k, 2)
        end if
      end do
      old_dist = -1
      do k = 1, nrbondsb(i)
        xh_dist = distance(bondlist(k, 2), hblist(i, 9))
        if (xh_dist > old_dist .and. bondlist(k, 2) /= hblist(i, 6)) then
          old_dist = xh_dist
          hblist(i, 7) = bondlist(k, 2)
        end if
      end do
      old_dist = -1
      do k = 1, nrbondsb(i)
        xh_dist = distance(bondlist(k, 2), hblist(i, 9))
        if (xh_dist > old_dist .and. bondlist(k, 2) /= hblist(i, 6) .and. bondlist(k, 2) /= hblist(i, 7)) then
          old_dist = xh_dist
          hblist(i, 8) = bondlist(k, 2)
        end if
      end do      ! possible forth atom should be h
    else if (nrbondsb(i) == 2) then
      old_dist = -1
      do k = 1, nrbondsb(i)
        xh_dist = distance(bondlist(k, 2), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 6) = bondlist(k, 2)
        end if
      end do
      do k = 1, nrbondsb(i)
        if (bondlist(k, 2) /= hblist(i, 6)) then
          hblist(i, 7) = bondlist(k, 2)
        end if
      end do
      if (distance(hblist(i, 5), hblist(i, 9)) < bonding(hblist(i, 5), hblist(i, 9), covrad)) then ! not needed
          hblist(i, 8) = hblist(i, 9)
      else
          hblist(i, 8) = hblist(i, 5)
      end if
    else if (nrbondsb(i) == 1) then
      hblist(i, 6) = bondlist(1, 2)
! need bonds on first bonded atom here
      nrbondsc = 0
      do k = 1, numat
!
! Check only H, N, and DO DO NOT USE!
!
      !   select case (nat(k))
      !   case (1, 7, 8)
          xc_dist = distance(k, hblist(i, 6))
          if (xc_dist < bonding(k, hblist(i, 6), covrad) .and. hblist(i, 6) /= k) then
            nrbondsc = nrbondsc + 1
            bondlist(nrbondsc, 3) = k
          end if
      !  end select
      end do
      old_dist = -1
      do k = 1, nrbondsc
        xh_dist = distance(bondlist(k, 3), hblist(i, 9))
        if (xh_dist > old_dist) then
          old_dist = xh_dist
          hblist(i, 7) = bondlist(k, 3)
        end if
      end do
      if (distance(hblist(i, 5), hblist(i, 9)) < bonding(hblist(i, 5), hblist(i, 9), covrad)) then ! not needed
          hblist(i, 8) = hblist(i, 9)
      else
          hblist(i, 8) = hblist(i, 5)
      end if
    else if (nrbondsb(i) == 0) then
      if (distance(hblist(i, 5), hblist(i, 9)) < bonding(hblist(i, 5), hblist(i, 9), covrad)) then ! not needed
          hblist(i, 6) = hblist(i, 9)
          hblist(i, 7) = hblist(i, 9)
          hblist(i, 8) = hblist(i, 9)
      else
          hblist(i, 6) = hblist(i, 5)
          hblist(i, 7) = hblist(i, 5)
          hblist(i, 8) = hblist(i, 5)
      end if
    else
      hbs2_ok = .false.
    end if
  end do
! check for problems
  if ( .not. hbs1_ok .or. .not. hbs2_ok ) then
!
!  No hydrogen bonds found
!
    l_h_bonds = .false.
  end if
end subroutine setup_DH_Plus
