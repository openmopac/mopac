  double precision function H_bonds4(l_grad, dxyz)
    use molkst_C, only : numat, keywrd, e_hb
    use elemts_C, only: elemnt
    use chanel_C, only : iw
    use common_arrays_C, only: nat
    implicit none
    double precision, intent (inout) :: dxyz(3,numat)
    logical, intent (in) :: l_grad
    double precision, external :: energy_corr_h4
    integer :: i
    double precision :: grad_h4(3,numat)
!
!  H4 Correction
!
    e_hb = energy_corr_h4(l_grad, grad_h4)
    H_bonds4 = e_hb
    if (l_grad) then
      if (index(keywrd, " DERIV") > 0) then
        write (iw, '(/25X,a)')"H4 CORRECTIONS" 
        write (iw, '(3X,       ''NUMBER  ATOM '',5X,''X'',12X,''Y'',12X,''Z'',/)') 
        write (iw, '(I6,4x,a2,F13.6,2F13.6)') (i, elemnt(nat(i)), grad_h4(:,i), i = 1,numat) 
      end if      
      dxyz = dxyz + grad_h4
    end if
    return
  end function H_bonds4
  
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
  
  function energy_corr_h4(l_grad, grad_h4)
  use common_arrays_C, only: coord, nat
  use molkst_C, only : numat, N_Hbonds, keywrd, numcal, method_PM6_D3H4, method_PM7, method_PM8
  implicit none
  double precision :: energy_corr_h4
  double precision, intent (out) :: grad_h4(3, numat)
  logical, intent (in) :: l_grad
!
!  Local variables
!
  integer :: i, j, h_i, d_i, a_i, k, o1, o2, cc, icalcn = -1
  double precision :: rda, rih, rjh, angle, rdh, rah, e_radial, d_radial, a, x, d, e_angular, &
    d_radial_d(3), d_radial_a(3), xd, d_angular, d_angular_d(3), d_angular_a(3), d_angular_h(3), &
    e_para, e_bond_switch, e_scale_w, e_scale_chd, e_scale_cha, rdhs, ravgs, d_bs, xd2, d_bs_d(3), d_bs_a(3), &
    d_bs_h(3), hydrogens, others, sign_wat, slope, v, fv, fv2, f_o1, f_o2, f_cc, cdist, &
    cv_o1, odist, cv_cc, cv_o2, e_corr_sum, g(3), e_corr, sum
  double precision :: HB_R_0 = 1.5d0, HR_R_CUTOFF = 5.5d0, M_PI = 3.141592653589793d0, &
    para_oh_o = 2.32d0, para_oh_n = 3.10d0, para_nh_o = 1.07d0, para_nh_n = 2.01d0, &
    multiplier_wh_o = 0.42d0, multiplier_nh4 = 3.61d0, multiplier_coo = 1.41d0
  double precision, external :: distance, cvalence_contribution, cvalence_contribution_d
  logical :: prt
  save :: icalcn, prt
    if (icalcn /= numcal) then  
!
!   Variables that change depending on the method
!
      if (method_PM6_D3H4) then
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
    prt = (index(keywrd," 0SCF ") + index(keywrd," PRT ") /= 0 .and. index(keywrd," DISP(") /= 0)
    grad_h4 = 0.d0
!
! Iterate over donor/acceptor pairs
!
    e_corr_sum = 0.d0
    N_Hbonds = 0
    do i = 1, numat
      if (nat(i) /= 7 .and. nat(i) /= 8) cycle
      do j = 1, i - 1
        if (nat(j) /= 7 .and. nat(j) /= 8) cycle
!
! Calculate donor-acceptor distance
!
        rda = distance(i, j)
        if (rda < HB_R_0 .or. rda > HR_R_CUTOFF) cycle
!
! Iterate over hydrogens
!
        do h_i = 1, numat
          if (nat(h_i) /= 1) cycle
!
! Distances to hydrogen
!
          rih = distance(i, h_i)
          rjh = distance(j, h_i)
          call bangle(coord, i, h_i, j, angle)
          angle = M_PI - angle
          if (angle >= M_PI/2.d0) cycle
!
! Here, we have filterd out everything but correct H-bonds
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
          if (l_grad) then
!
! Radial gradient
!
            d_radial = -0.02123851851851194655d0*rda**6 &
                        +0.44145777777762551519d0*rda**5 &
                        -3.50435555555413991158d0*rda**4 &
                      +13.01238518517846998179d0*rda**3 &
                      -21.62062222220516360949d0*rda**2 &
                      +10.63509333331144368628d0*rda &
                        +3.40736000001102778967d0
!
! Cartesian gradients on D and A atoms
!
            d_radial_d(:) = (coord(:, d_i) - coord(:, a_i))/rda *d_radial
            d_radial_a = -d_radial_d                  
          end if
!
! Angular term
!
          a = angle/(M_PI/2.d0)
          x = -20.d0*a**7 + 70.d0*a**6 - 84.d0*a**5 + 35.d0*a**4  
          e_angular = 1.d0 - x*x
          if (l_grad) then
            xd = (-140.d0*a**6 + 420.d0*a**5 - 420.d0*a**4 + 140.d0*a**3)/(M_PI/2.d0)
            d_angular = -2.d0*xd*x
!
! Dot product of bond vectors
!
            d = (coord(1, d_i) - coord(1, h_i))*(coord(1, a_i) - coord(1, h_i)) + &
                (coord(2, d_i) - coord(2, h_i))*(coord(2, a_i) - coord(2, h_i)) + &
                (coord(3, d_i) - coord(3, h_i))*(coord(3, a_i) - coord(3, h_i)) 
            x = -d_angular/sqrt(1.d0 - (d*d)/(rdh*rdh)/(rah*rah))
!
! Donor atom
!
            d_angular_d(:) = -x*((coord(:, a_i) - coord(:, h_i))/(rdh*rah) - (coord(:, d_i) - coord(:, h_i))*d/(rdh**3*rah))
! 
! Acceptor atom
!
            d_angular_a(:) = -x*((coord(:, d_i) - coord(:, h_i))/(rdh*rah) - (coord(:, a_i) - coord(:, h_i))*d/(rah**3*rdh))
!
! Hydrogen
!
            d_angular_h(:) =  -d_angular_d(:) -  d_angular_a(:)
          end if
!
! Energy coefficient
!
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
!
! Gradient
!
            if (l_grad) then
              d_bs = -(-140.d0*x**6 + 420.d0*x**5 - 420.d0*x**4 + 140.d0*x**3)
              xd = d_bs/ravgs
              xd2 = -0.5d0*d_bs*x/ravgs
              d_bs_d(:) = (coord(:, d_i) - coord(:, h_i))*(xd + xd2)/rdh 
              d_bs_a(:) = (coord(:, a_i) - coord(:, h_i))*xd2/rah
              d_bs_h(:) = -d_bs_d(:) - d_bs_a(:)                       
            end if
          else
!
! No switching, no gradient
!
            e_bond_switch = 1.d0
            if (l_grad) then
              d_bs_d(:) = 0.d0
              d_bs_a(:) = 0.d0
              d_bs_h(:) = 0.d0
            end if
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
                hydrogens = hydrogens + cvalence_contribution(d_i, k)
              else
                others = others + cvalence_contribution(d_i, k)
              end if
            end do
!
! If it is water
!
            if (hydrogens >= 1.d0) then
              sign_wat = 1.d0
              slope = multiplier_wh_o -1.d0
              v = hydrogens
              fv = 0.d0
              if (v > 1.d0 .and. v <= 2.d0) then
                fv = v - 1.d0
                sign_wat = 1.d0
              end if
              if (v > 2.d0 .and. v < 3.d0) then
                fv = 3.d0 -v
                sign_wat = -1.d0
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
              v = v + cvalence_contribution(d_i, k)
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
              v = cvalence_contribution(o1, k)
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
                v = cvalence_contribution(cc, k)
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
                cv_o2 = cv_o2 +cvalence_contribution(o2, k)
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
          if (prt) call prt_hbonds(D_i, H_i, A_i, e_corr)
          e_corr_sum = e_corr_sum + e_corr
!
! Total gradient
! radial
!
          grad_h4(:, d_i) = grad_h4(:, d_i) + d_radial_d(:)*e_para*e_angular*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
          grad_h4(:, a_i) = grad_h4(:, a_i) + d_radial_a(:)*e_para*e_angular*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
!
! angular
!
          grad_h4(:, d_i) = grad_h4(:, d_i) + d_angular_d(:)*e_para*e_radial*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
          grad_h4(:, a_i) = grad_h4(:, a_i) + d_angular_a(:)*e_para*e_radial*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
                grad_h4(:, h_i) = grad_h4(:, h_i) + d_angular_h(:)*e_para*e_radial*e_bond_switch*e_scale_w*e_scale_chd*e_scale_cha
!
! bond_switch
!
          grad_h4(:, d_i) = grad_h4(:, d_i) + d_bs_d(:)*e_para*e_radial*e_angular*e_scale_w*e_scale_chd*e_scale_cha
          grad_h4(:, a_i) = grad_h4(:, a_i) + d_bs_a(:)*e_para*e_radial*e_angular*e_scale_w*e_scale_chd*e_scale_cha
          grad_h4(:, h_i) = grad_h4(:, h_i) + d_bs_h(:)*e_para*e_radial*e_angular*e_scale_w*e_scale_chd*e_scale_cha
!
! water scaling
!
          if (l_grad .and. abs(e_scale_w - 1.d0) > 1.d-14) then
            slope = multiplier_wh_o - 1.d0
            do k = 1, numat
              if (k == d_i) cycle
              x = distance(d_i,k)
              if (nat(k) == 1) then
                xd = cvalence_contribution_d(d_i, k)*sign_wat
                g(:) = -(coord(:, d_i) - coord(:, k))*xd*slope/x
                grad_h4(:, d_i) = grad_h4(:, d_i) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_cha
                grad_h4(:, k  ) = grad_h4(:, k  ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_cha
              else 
                xd = cvalence_contribution_d(d_i, k)
                g(:) = -(coord(:, d_i) - coord(:, k))*xd*slope/x
                grad_h4(:, d_i) = grad_h4(:, d_i) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_cha
                grad_h4(:, k  ) = grad_h4(:, k  ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_cha
              end if
            end do
          end if        
!
! scaled groups: NR4+
!
          if (l_grad .and. abs(e_scale_chd - 1.d0) > 1.d-14) then
            slope = multiplier_nh4 - 1.d0
            do k = 1, numat
              if (k /= d_i) then
                x = distance(d_i, k)
                xd = cvalence_contribution_d(d_i, k)
                g(:) = -(coord(:, d_i) - coord(:, k))*xd*slope/x
                grad_h4(:, d_i) = grad_h4(:, d_i) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_cha*e_scale_w
                grad_h4(:, k  ) = grad_h4(:, k  ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_cha*e_scale_w
              end if
            end do
          end if
!
! scaled groups: COO-
!
          if (.not. l_grad .or. f_o1*f_o2*f_cc == 0.d0) cycle
          slope = multiplier_coo - 1.d0
!
! Atoms around O1
!
          do k = 1, numat
            if (k == o1) cycle
            xd = cvalence_contribution_d(o1, k)
            if (xd == 0.d0) cycle
            x =  distance(o1, k)
            if (cv_o1 > 1.d0) xd = -xd
            xd = xd*(f_o2*f_cc)
            g(:) = -(coord(:, o1) - coord(:, k))*xd*slope/x
            grad_h4(:, o1) = grad_h4(:, o1) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w
            grad_h4(:, k ) = grad_h4(:, k ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w             
          end do
          slope = multiplier_coo - 1.d0
!
! Atoms around O2
!
          do k = 1, numat
            if (k == o2) cycle
            xd = cvalence_contribution_d(o2, k)
            if (xd == 0.d0) cycle
            x =  distance(o2, k)
            if (cv_o2 > 1.d0) xd = -xd
            xd = xd*(f_o1*f_cc)
            g(:) = -(coord(:, o2) - coord(:, k))*xd*slope/x
            grad_h4(:, o2) = grad_h4(:, o2) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w
            grad_h4(:, k ) = grad_h4(:, k ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w
          end do
          slope = multiplier_coo - 1.d0
          do k = 1, numat
            if (k == cc) cycle
            xd = cvalence_contribution_d(cc, k)
            if (xd == 0.d0) cycle
            x =  distance(cc, k)
            if (cv_cc > 3.d0) xd = -xd
            xd = xd*(f_o1*f_o2)
            g(:) = -(coord(:, cc) - coord(:, k))*xd*slope/x
            grad_h4(:, cc) = grad_h4(:, cc) - g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w
            grad_h4(:, k ) = grad_h4(:, k ) + g(:)*e_para*e_radial*e_angular*e_bond_switch*e_scale_chd*e_scale_w
          end do                                    
        end do
      end do
    end do
    energy_corr_h4 = e_corr_sum
  end function energy_corr_h4
  
  function cvalence_contribution(atom_a, atom_b)
    use common_arrays_C, only: nat
    use radii_C, only : covalent_radii
    implicit none
    integer, intent (in) :: atom_a, atom_b
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
    end if
    return     
  end function cvalence_contribution
  
  function cvalence_contribution_d(atom_a, atom_b)
    use common_arrays_C, only: nat
    use radii_C, only : covalent_radii
    implicit none
    integer, intent (in) :: atom_a, atom_b
    double precision :: cvalence_contribution_d
    double precision :: r, ri, rj, r0, r1, x
    double precision, external :: distance
    ri = covalent_radii(nat(atom_a))
    rj = covalent_radii(nat(atom_b))
    r0 = ri + rj
    r1 = r0*1.6d0
    r = distance(atom_a,atom_b)
    if (r == 0.d0 .or. r >= r1 .or. r <= r0) then
      cvalence_contribution_d = 0.d0
    else
      x = (r - r0)/(r1 - r0)
      cvalence_contribution_d = -(-140.d0*x**6 + 420.d0*x**5 - 420.d0*x**4 + 140.d0*x**3)/(r1 - r0)
    end if
  end function cvalence_contribution_d
  
  function energy_corr_hh_rep(l_grad, dxyz) 
    use common_arrays_C, only: nat, coord
    use molkst_C, only : numat, keywrd, e_hh
    use chanel_C, only : iw
    use elemts_C, only: elemnt
    implicit none
    logical, intent (in) :: l_grad
    double precision, intent (inout) :: dxyz(3, numat)
    double precision :: grad_hh(3, numat)
    double precision :: e_corr, r, d_rad, g(3), energy_corr_hh_rep
    integer :: i, j ! iteration counters
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
        r = distance(i,j)
        e_corr = poly(r, l_grad, d_rad)
        e_hh = e_hh + e_corr
        if (l_grad) then
!
! Cartesian components of the gradient
!
          g(:) = -(coord(:,i) - coord(:,j))/r*d_rad
!
! Add pair contribution to the global gradient
!  
          grad_hh(:,i) = grad_hh(:,i) - g(:)
          grad_hh(:,j) = grad_hh(:,j) + g(:)
        end if
      end do
    end do
    energy_corr_hh_rep = e_hh
    if (l_grad) then
      dxyz = dxyz + grad_hh
      if (index(keywrd, " DERIV") > 0) then
        write (iw, '(/25X,a)')"HH REPULSION" 
        write (iw, '(3X,       ''NUMBER  ATOM '',5X,''X'',12X,''Y'',12X,''Z'',/)') 
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

