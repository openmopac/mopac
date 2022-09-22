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

   module radii_C
   double precision :: covalent_radii(118)
   data covalent_radii /0.37d0, 0.32d0, 1.34d0, 0.9d0, 0.82d0, 0.77d0, &
     0.75d0, 0.73d0, 0.71d0, 0.69d0, 1.54d0, 1.3d0, 1.18d0, 1.11d0, 1.06d0, 1.02d0, 0.99d0, 0.97d0, &
     1.96d0, 1.74d0, 1.44d0, 1.36d0, 1.25d0, 1.27d0, 1.39d0, 1.25d0, 1.26d0, 1.21d0, 1.38d0, 1.31d0, &
     1.26d0, 1.22d0, 1.19d0, 1.16d0, 1.14d0, 1.1d0, 2.11d0, 1.92d0, 1.62d0, 1.48d0, 1.37d0, 1.45d0, &
     1.56d0, 1.26d0, 1.35d0, 1.31d0, 1.53d0, 1.48d0, 1.44d0, 1.41d0, 1.38d0, 1.35d0, 1.33d0, 1.3d0, &
     2.25d0, 1.98d0, 1.69d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 1.6d0, 1.5d0, 1.38d0, 1.46d0, 1.59d0, 1.28d0, 1.37d0, 1.28d0, 1.44d0, 1.49d0, 0.0d0, &
     0.0d0, 1.46d0, 0.0d0, 0.0d0, 1.45d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0/
  end module radii_C

  double precision function H_bonds4(H_i, i, j, d_list, nd_list)
!
! Rezac J., Hobza P., "Advanced Corrections of Hydrogen Bonding and
! Dispersion for Semiempirical Quantum Mechanical Methods", J. Chem.
! Theory and Comp 8, 141-151 (2012).
!
  use common_arrays_C, only: coord, nat
  use molkst_C, only : numat, N_Hbonds, numcal, method_PM6_D3H4, method_PM7, method_PM8, &
    method_pm6_d3h4x, method_pm6_org
  implicit none
  integer, intent (in) :: H_i, i, j
  integer, intent (out) :: d_list(20), nd_list
!
!  Local variables
!
  integer :: d_i, a_i, k, o1, o2, cc, icalcn = -1
  double precision :: rda, rih, rjh, angle, rdh, rah, e_radial, a, x, e_angular, &
    e_para, e_bond_switch, e_scale_w, e_scale_chd, e_scale_cha, rdhs, ravgs, &
    hydrogens, others, slope, v, fv, fv2, f_o1, f_o2, f_cc, cdist, &
    cv_o1, odist, cv_cc, cv_o2, e_corr_sum, e_corr, sum
  double precision :: M_PI = 3.141592653589793d0, &
    para_oh_o = 2.32d0, para_oh_n = 3.10d0, para_nh_o = 1.07d0, para_nh_n = 2.01d0, &
    multiplier_wh_o = 0.42d0, multiplier_nh4 = 3.61d0, multiplier_coo = 1.41d0
  double precision, external :: distance, cvalence_contribution
  save :: icalcn
    cv_cc = 0.d0
    if (icalcn /= numcal) then
!
!   Variables that change depending on the method
!
      if (method_PM6_D3H4 .or. method_pm6_d3h4x .or. method_pm6_org) then
        icalcn = numcal
        para_oh_o = 2.32d0
        para_oh_n = 3.10d0
        para_nh_o = 1.07d0
        para_nh_n = 2.01d0
        multiplier_wh_o = 0.42d0
        multiplier_nh4 = 3.61d0
        multiplier_coo = 1.41d0
      else if (method_PM7 .or. method_PM8) then
        para_oh_o = 2.32d0
        para_oh_n = 3.10d0
        para_nh_o = 1.07d0
        para_nh_n = 2.01d0
        multiplier_wh_o = 0.42d0
        multiplier_nh4 = 3.61d0
        multiplier_coo = 1.41d0
      else
        para_oh_o = 2.32d0
        para_oh_n = 3.10d0
        para_nh_o = 1.07d0
        para_nh_n = 2.01d0
        multiplier_wh_o = 0.42d0
        multiplier_nh4 = 3.61d0
        multiplier_coo = 1.41d0
      end if
    end if
    d_list(1) = i
    d_list(2) = j
    d_list(3) = H_i
    nd_list = 3
!
! Iterate over donor/acceptor pairs
!
    e_corr_sum = 0.d0
    N_Hbonds = 0
! Calculate donor-acceptor distance
!
    rda = distance(i, j)
!
! Iterate over hydrogens
!

!
! Distances to hydrogen
!
    rih = distance(i, h_i)
    rjh = distance(j, h_i)
    call bangle(coord, i, h_i, j, angle)
    angle = M_PI - angle
    if (angle >= M_PI/2.d0) then
      continue
    end if
!
! Here, we have filtered out everything but correct H-bonds
! Determine donor and acceptor - donor is the closer one
!
    if (rih < rjh) then
      d_i = i
      a_i = j
      rdh = rih
      rah = rjh
    else
      d_i = j
      a_i = i
      rdh = rjh
      rah = rih
    end if
!
! Radial term
!
    e_radial = -0.00303407407407313510d0*rda**7 &
                +0.07357629629627092382d0*rda**6 &
                -0.70087111111082800452d0*rda**5 &
                +3.25309629629461749545d0*rda**4 &
                -7.20687407406838786983d0*rda**3 &
                +5.31754666665572184314d0*rda**2 &
                +3.40736000001102778967d0*rda  &
                -4.68512000000450434811d0
!
!  The e_radial function has a minimum at rda = 5.5 Angstroms:
!  rda  e_radial   rda    e_radial
!  4.50 -0.143     6.75   -1.531
!  4.75 -0.056     7.00   -3.601
!  5.00 -0.013     7.25   -7.527
!  5.25 -0.001     7.50  -14.418
!  5.50  0.000     7.75  -25.819
!  5.75 -0.001     8.00  -43.815
!  6.00 -0.026     8.25  -71.155
!  6.25 -0.152     8.50 -111.391
!  6.50 -0.550     8.75 -169.028
!
! This is set in "cutoff" in find_H__Y_bonds
!
! Angular term
!
    a = angle/(M_PI/2.d0)
    x = -20.d0*a**7 + 70.d0*a**6 - 84.d0*a**5 + 35.d0*a**4
    e_angular = 1.d0 - x*x
!
! Energy coefficient
!
    e_para = 0.d0
    if (nat(d_i) == 8 .and. nat(a_i) == 8) e_para = para_oh_o
    if (nat(d_i) == 8 .and. nat(a_i) == 7) e_para = para_oh_n
    if (nat(d_i) == 7 .and. nat(a_i) == 8) e_para = para_nh_o
    if (nat(d_i) == 7 .and. nat(a_i) == 7) e_para = para_nh_n
!
! Bond switching
!
    if (rdh > 1.15d0) then
      rdhs = rdh - 1.15d0
      ravgs = 0.5d0*rdh + 0.5d0*rah - 1.15d0
      x = rdhs/ravgs
      e_bond_switch = 1.d0 - (-20.d0*x**7 + 70.d0*x**6 -84.d0*x**5 + 35.d0*x**4)
    else
!
! No switching, no gradient
!
      e_bond_switch = 1.d0
    end if
!
! Water scaling
!
    e_scale_w = 1.d0
    if (nat(d_i) == 8 .and. nat(a_i) == 8) then
      hydrogens = 0.d0
      others = 0.d0
      do k = 1, numat
        if (nat(k) == 1) then
          hydrogens = hydrogens + cvalence_contribution(d_i, k, d_list, nd_list)
        else
          others = others + cvalence_contribution(d_i, k, d_list, nd_list)
        end if
      end do
!
! If it is water
!
      if (hydrogens >= 1.d0) then
        slope = multiplier_wh_o -1.d0
        v = hydrogens
        fv = 0.d0
        if (v > 1.d0 .and. v <= 2.d0) then
          fv = v - 1.d0
        end if
        if (v > 2.d0 .and. v < 3.d0) then
          fv = 3.d0 -v
        end if
        fv2 = max(1.d0 - others, 0.d0)
        e_scale_w = 1.d0 + slope*fv*fv2
      end if
    end if
!
! Charged groups
!
    e_scale_chd = 1.d0
    e_scale_cha = 1.d0
!
! Scaled groups: NR4+
!
    if (nat(d_i) == 7) then
      slope = multiplier_nh4 - 1.d0
      v = 0.d0
      do k = 1, numat
        v = v + cvalence_contribution(d_i, k, d_list, nd_list)
      end do
      if (v > 3.d0) then
        v = v - 3.d0
      else
        v = 0.d0
      end if
      e_scale_chd = 1.d0 + slope*v
    end if
!
! Scaled groups: COO-
!
    f_o1 = 0.d0
    f_o2 = 0.d0
    f_cc = 0.d0
    o1 = a_i
    o2 = -1
    cc = -1
    if (nat(a_i) == 8) then
      slope = multiplier_coo - 1.d0
!
! Search for closest C atom
!
      cdist = 9.9d9
      cv_o1 = 0.d0
      do k = 1, numat
        v = cvalence_contribution(o1, k, d_list, nd_list)
        cv_o1 = cv_o1 + v! Sum O1 valence
        sum = distance(o1, k)
        if (v > 0.d0 .and. nat(k) == 6 .and. sum < cdist) then
          cdist =  sum
          cc = k
        end if
      end do
!
! If C found, look for the second O
!
      if (cc /= -1) then
        odist = 9.9d9
        cv_cc = 0.d0
        do k = 1, numat
          v = cvalence_contribution(cc, k, d_list, nd_list)
          cv_cc = cv_cc + v
          sum = distance(cc, k)
          if ((v > 0.d0) .and. (k /= o1) .and. (nat(k) == 8) .and. (sum < odist)) then
            odist =  sum
            o2 = k
          end if
        end do
      end if
!
! O1-C-O2 triad:
!
      if (o2 /= -1) then
!
! Get O2 valence
!
        cv_o2 = 0.d0
        do k = 1, numat
          cv_o2 = cv_o2 +cvalence_contribution(o2, k, d_list, nd_list)
        end do
        f_o1 = max(1.d0 - abs(1.d0 - cv_o1), 0.d0)
        f_o2 = max(1.d0 - abs(1.d0 - cv_o2), 0.d0)
        f_cc = max(1.d0 - abs(3.d0 - cv_cc), 0.d0)
        e_scale_cha = 1.d0 + slope*f_o1*f_o2*f_cc
      end if
    end if
!
! Final energy
!
    e_corr = e_para*e_radial*e_angular*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
    if (e_corr < -1.d0) then
      N_Hbonds = N_Hbonds + 1
    end if
    e_corr_sum = e_corr_sum  + e_corr
    H_bonds4 = e_corr_sum
  end function H_bonds4




  function energy_corr_hh_rep(l_grad, dxyz)
    use common_arrays_C, only: nat, cell_ijk, Vab
    use molkst_C, only : numat, keywrd, e_hh, l123, l1u, l2u, l3u
    use chanel_C, only : iw
    use elemts_C, only: elemnt
    implicit none
    logical, intent (in) :: l_grad
    double precision, intent (inout) :: dxyz(3, numat*l123)
    double precision :: grad_hh(3, numat*l123)
    double precision :: e_corr, r, d_rad, g(3), energy_corr_hh_rep
    integer :: i, j, iii, jjj, i_cell, j_cell, kkkk
    double precision, external :: distance, poly
    e_hh = 0.d0
    grad_hh = 0.d0
!
! Iterate over H atoms twice
!
    do i = 1, numat
      if (nat(i) /= 1) cycle
      do j = 1, i - 1
        if (nat(j) /= 1) cycle
!
! Calculate distance
!
        r = distance(j, i)
        if (cell_ijk(1) /= 0) then
          continue
        end if

        e_corr = poly(r, l_grad, d_rad)
        e_hh = e_hh + e_corr
        if (l_grad) then
!
! Cartesian components of the gradient
!
          g(:) = Vab(:)/r*d_rad
!
! Add pair contribution to the global gradient
!
          iii = l123*(i - 1)
          jjj = l123*(j - 1)
          kkkk = (l3u - cell_ijk(3)) + (2*l3u + 1)*(l2u - cell_ijk(2) + (2*l2u + 1)*(l1u - cell_ijk(1))) + 1
          i_cell = iii + kkkk
          j_cell = jjj + kkkk
          grad_hh(:,i_cell) = grad_hh(:,i_cell) - g(:)
          grad_hh(:,j_cell) = grad_hh(:,j_cell) + g(:)
        end if
      end do
    end do
    energy_corr_hh_rep = e_hh
    if (l_grad) then
      dxyz = dxyz + grad_hh
      if (index(keywrd, " DERIV") > 0) then
        write (iw, '(/25X,a)')"HH REPULSION"
        write (iw, '(" NUMBER  ATOM  ",5X,"X",12X,"Y",12X,"Z",/)')
        write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i, elemnt(nat(i)), grad_hh(:,i), i = 1,numat)
      end if
    end if
    return
  end function energy_corr_hh_rep
  double precision function poly(r, l_grad, dpoly)
    implicit none
    double precision :: r, dpoly
    logical :: l_grad
!
!  poly is the hydrogen-hydrogen correction for two atoms separated by "r" Angstroms
!  dpoly is the gradient of the internal coordinate
!
    if (r <= 1.d0) then
      poly = 25.46293603147693d0
      dpoly = 0.d0
    else if ( r > 1.d0 .and. r < 1.5d0) then
      poly =  -2714.952351603469651d0 * r**5 &
             +17103.650110591705015d0 * r**4 &
             -42511.857982217959943d0 * r**3 &
             +52063.196799138342612d0 * r**2 &
             -31430.658335972289933d0 * r    &
              +7516.084696095140316d0
      if (l_grad) &
        dpoly = -2714.952351603469651d0*5.d0 * r**4 &
               +17103.650110591705015d0*4.d0 * r**3 &
               -42511.857982217959943d0*3.d0 * r**2 &
               +52063.196799138342612d0*2.d0 * r    &
               -31430.658335972289933d0
    else
      poly = 118.7326d0*exp(-1.53965d0*(r**1.72905d0))
      if (l_grad) &
      dpoly = -1.53965d0*1.72905d0*r**0.72905d0*118.7326d0*exp(-1.53965d0*(r**1.72905d0))
    end if
    return
  end function poly


   function cvalence_contribution(atom_a, atom_b, d_list, nd_list)
    use common_arrays_C, only: nat
    use radii_C, only : covalent_radii
    implicit none
    integer, intent (in) :: atom_a, atom_b
    integer, intent (inout) :: d_list(20), nd_list
    integer :: i
    double precision :: r, ri, rj, r0, r1, x, cvalence_contribution
    double precision, external :: distance
    ri = covalent_radii(nat(atom_a))
    rj = covalent_radii(nat(atom_b))
    r0 = ri + rj
    r1 = r0*1.6d0
    r = distance(atom_a, atom_b)
    if (r == 0.d0 .or. r >= r1) then
      cvalence_contribution = 0.d0
    else if (r <= r0) then
      cvalence_contribution = 1.d0
    else
      x = (r - r0)/(r1 - r0)
      cvalence_contribution = 1.d0 - (-20.d0*x**7 + 70.d0*x**6 -84.d0*x**5 + 35.d0*x**4)
      do i = 1, nd_list
        if (d_list(i) == atom_b) exit
      end do
      if (i > nd_list) then
        nd_list = nd_list + 1
        d_list(nd_list) = atom_b
      end if
    end if
    return
  end function cvalence_contribution
