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

double precision function PM6_DH_H_bond_corrections(l_grad, prt)
!
!    Add a dispersion E_disp, a coulombic, EC, and a repulsive correction, ER,
!    to improve intermolecular interaction energies.
!
  use molkst_C, only : numat, E_hb, N_Hbonds, method_pm7, l123, l1u, l2u, l3u, &
    method_pm6_dh_plus, method_pm6_dh2, method_pm6_dh2x, &
    line, id, numcal
  use parameters_C, only : tore
  use common_arrays_C, only: coord, nat, q, p, hblist, dxyz, cell_ijk
  implicit none
  logical, intent (in) :: l_grad, prt
!
! Local variables
!
  integer :: i, j, k, l, D, H, A, nd_list, d_list(9), d_l(9), i_cell, ii, iii, kkkk, &
    nrpairs, max_h_bonds, icalcn = -1
  logical :: n_h_bonds, first = .true.
  integer, allocatable :: nrbondsa(:), nrbondsb(:)
  double precision, allocatable :: vector(:)
  double precision :: EC, ER, sum, sum1, delta = 1.d-5, covrad(94)
  double precision, external :: EC_plus_ER, EH_plus
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
    PM6_DH_H_bond_corrections = 0.d0
    return
  end if
  if (allocated(nrbondsa))   deallocate(nrbondsa)
  if (allocated(nrbondsb))   deallocate(nrbondsb)
  if (allocated(hblist))     deallocate(hblist)
  if (allocated(vector))     deallocate(vector)
  allocate (hblist(max_h_bonds,10), nrbondsa(max_h_bonds), nrbondsb(max_h_bonds), vector(numat), stat=i)
    if (i /= 0) then
    line = " Cannot allocate arrays for PM6-DH+"
    call to_screen(trim(line))
    call mopend(trim(line))
    PM6_DH_H_bond_corrections = 0.d0
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
      sum = EC_plus_ER(hblist(ii,1), hblist(ii,2), hblist(ii,3), q(hblist(ii,2)), q(hblist(ii,3)), EC, ER, d_list, nd_list)
      E_hb = E_hb + EC + ER
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
            else
              sum1 = EC_plus_ER(D, H, A, q(H), q(A), EC, ER, d_l, l)
            end if
            sum1 = (sum1 - sum)/delta
            if (Abs(sum1) < 50.d0) then
              dxyz(i_cell*3 + i) = dxyz(i_cell*3 + i) + sum1
            end if
            coord(i,k) = coord(i,k) - delta
          end do
        end if
      end do
    end if
  end do
  PM6_DH_H_bond_corrections = E_hb
  return
  end function PM6_DH_H_bond_corrections
