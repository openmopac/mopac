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

  double precision function EH_plus(i, hblist, max_h_bonds, nrbondsa, nrbondsb)
!
!  Third-Generation Hydrogen-Bonding Corrections for Semiempirical QM Methods and Force Fields
!  Martin Korth, J. Chem. Theory Comput., 2010, 6 (12), pp 3808-3816
!
    use common_arrays_C, only : nat
    use molkst_C, only : method_PM7
    use funcon_C, only : eV, fpc_9, a0, pi
  !  use parameters_C, only : par8, par9
    implicit none
    integer, intent (in) :: max_h_bonds, i, hblist(max_h_bonds,10), nrbondsa(max_h_bonds), nrbondsb(max_h_bonds)
    double precision :: angle_cos, angle, torsion_check, & ! weight_956, &
      angle2_shift, angle2_shift_2, torsion_shift, angle2, torsion_check_bac, &
      angle2_cos, angle2_cos_2, torsion_correct, torsion_value, torsion_cos, torsion_value_2, &
      torsion_cos_2, angle2_cos_new, angle2_cos_2_new, torsion_cos_new, torsion_cos_2_new, &
      scale_a, scale_nsp3, scale_nsp2, scale_osp3, scale_osp2, scale_b, scale_c, &
      hb_dist, xc_dist, ha_dist, damping, hartree2kcal, XY_dist, short
    double precision :: shortcut = 2.4d0    !  2.4 Angstrom
    double precision :: longcut  = 7.0d0    !  7.0 Angstrom
    double precision :: covcut   = 1.2d0    !  1.2 Angstrom (hb mid point)
    logical :: torsion_check_set, torsion_check_set2
    double precision, external :: torsion, distance
    torsion_value = 0.d0
    torsion_shift = 0.d0
    angle2_shift = 0.d0
    angle2_shift_2 = 0.d0
    if (method_PM7) then
       scale_nsp3 = -0.171271D0*a0**2
       scale_osp3 = -0.098822D0*a0**2
       scale_nsp2 = -0.171271D0*a0**2
    else
       scale_nsp3 = -0.16d0*a0**2
       scale_osp3 = -0.12d0*a0**2
       scale_nsp2 = scale_nsp3
    end if
    scale_osp2 = scale_osp3
    hartree2kcal = eV*fpc_9
    EH_plus = 0.d0
    if (hblist(i, 10) /= -666) then
! first angle
      angle_cos = -cos(angle(hblist(i, 1), hblist(i, 9), hblist(i, 5)))   !  cos(pi-angle)
      if (angle_cos <= 0) return
! get target angles for hblist(i, 1)
! - no explicit lone pair representation to keep gradient simpler
! - no coordinate dependent target shifts for the same reason (except NR3)
      torsion_check = 0.d0
      torsion_check_set = .false.
      torsion_check_set2 = .false.
      if (nat(hblist(i, 1)) == 8) then ! sulphur?
        if (nrbondsa(i) == 1) then
!
!  >C=O - - H-X
!
          angle2_shift = pi
          angle2_shift_2 = pi/(180.d0/120.d0)
          torsion_shift = 0.d0
          torsion_check_set2 = .true.! H...O = C -R1
        else
!
! >O - - H-X
!
          angle2_shift = pi/(180.d0/109.48)
          angle2_shift_2 = angle2_shift
          torsion_shift = pi/(180.d0/54.74)
        end if
      else if (nat(hblist(i, 1)) == 7) then
        if (nrbondsa(i) == 2) then
!
! >N - - H-X
!
          angle2_shift = pi/(180.d0/120.d0)
          angle2_shift_2 = angle2_shift
          torsion_shift = 0.d0
        else
!
! >N- - - H-X   WHAT ABOUT Cyanide?
!
          angle2_shift = pi/(180.d0/109.48)
          angle2_shift_2 = angle2_shift
          torsion_shift = pi/(180.d0/54.74)
          torsion_check_set = .true.! NR3 group
        end if
      end if
!
! extrapolation between tetragonal and planar NR3 group
!
      if (torsion_check_set) then
        torsion_check = torsion(hblist(i, 3), hblist(i, 2), hblist(i, 1), hblist(i, 4))   !  torsion2
        if (torsion_check <= -pi) torsion_check = torsion_check  +  2.d0 * pi
        if (torsion_check > pi) torsion_check = torsion_check - 2.d0 * pi
        if (torsion_check < 0) then
          torsion_check = -pi -torsion_check
        else
          torsion_check = pi -torsion_check
        end if
        torsion_check_bac = torsion_check ! save sign
        if (torsion_check < 0) torsion_check = -1.d0 * torsion_check
        torsion_check = torsion_check*180.d0/pi
        torsion_shift = torsion_shift + pi/(180.d0/((54.74d0 - torsion_check)/54.74d0*35.26d0))
        angle2_shift = angle2_shift - pi/(180.d0/((54.74d0 - torsion_check)/54.74d0*19.48d0))
        angle2_shift_2 = angle2_shift
        torsion_check = torsion_check_bac ! restore sign
      end if
! second angle
      angle2 = angle(hblist(i, 2), hblist(i, 1), hblist(i, 9))  !angle2
      angle2_cos = cos(angle2_shift - angle2)
      angle2_cos_2 = cos(angle2_shift_2 - angle2)
      if (angle2_cos_2 > angle2_cos) angle2_cos = angle2_cos_2
      if (angle2_cos <= 0.d0) return
! torsion angle
      torsion_correct = torsion(hblist(i, 3), hblist(i, 2), hblist(i, 1), hblist(i, 9))   !  torsion1
      if (torsion_correct <= -pi) torsion_correct = torsion_correct + 2.d0 * pi
      if (torsion_correct  >  pi) torsion_correct = torsion_correct - 2.d0 * pi
      if ((.not.torsion_check_set2) .or. (abs(torsion_correct*180.d0/pi) > 90.d0)) then
        if (torsion_correct < 0.d0) then
          torsion_correct = -pi - torsion_correct
        else
          torsion_correct =  pi - torsion_correct
        end if
      end if
!
! correction of NR3 torsion angle for through-bond case
!
      if (torsion_check < 0) then ! negative torsion angle occupied by -NR3 r3
        torsion_value = torsion_shift - torsion_correct
        if (torsion_value <= -pi) torsion_value = torsion_value + 2.d0 * pi
        if (torsion_value  >  pi) torsion_value = torsion_value - 2.d0 * pi
        torsion_cos = cos(torsion_value)
      else if (torsion_check > 0) then ! positive torsion angle occupied by -NR3 r3
        torsion_value = -torsion_shift - torsion_correct
        if (torsion_value <= -pi) torsion_value = torsion_value + 2.d0 * pi
        if (torsion_value  >  pi) torsion_value = torsion_value - 2.d0 * pi
        torsion_cos = cos(torsion_value)
      else ! planar -NR3 or general case
        torsion_value = torsion_shift - torsion_correct
        torsion_value_2 = -torsion_shift - torsion_correct
        if (torsion_value   <= -pi) torsion_value   = torsion_value   + 2.d0 * pi
        if (torsion_value   >   pi) torsion_value   = torsion_value   - 2.d0 * pi
        if (torsion_value_2 <= -pi) torsion_value_2 = torsion_value_2 + 2.d0 * pi
        if (torsion_value_2 >   pi) torsion_value_2 = torsion_value_2 - 2.d0 * pi
        torsion_cos = cos(torsion_value)
        torsion_cos_2 = cos(torsion_value_2)
        if (torsion_cos_2 > torsion_cos) torsion_cos = torsion_cos_2
      end if
      if (distance(hblist(i, 9), hblist(i, 1)) > distance(hblist(i, 9), hblist(i, 2)) .and. torsion_check_set2) &
        torsion_cos = 0.d0
      if (hblist(i, 3) == hblist(i, 4) .or. hblist(i, 9) == hblist(i, 4)) torsion_cos = 1.d0
      if (torsion_cos < 0) return
!
! get target angles for hblist(i, 5) - see comment above
!
      torsion_check      = 0.d0
      torsion_check_set  = .false.
      torsion_check_set2 = .false.
      if (nat(hblist(i, 5)) == 8) then ! sulphur?
        if (nrbondsb(i) == 1) then
          angle2_shift = pi
          angle2_shift_2 = pi/(180.d0/120.d0)
          torsion_shift = 0.d0
          torsion_check_set2 = .true.! H...O = C -R1
        else
          angle2_shift = pi/(180.d0/109.48d0)
          angle2_shift_2 = angle2_shift
          torsion_shift = pi/(180.d0/54.74d0)
        end if
      else if (nat(hblist(i, 5)) == 7) then
        if (nrbondsb(i) == 2) then
          angle2_shift = pi/(180.d0/120.d0)
          angle2_shift_2 = angle2_shift
          torsion_shift = 0.d0
        else
          angle2_shift = pi/(180.d0/109.48d0)
          angle2_shift_2 = angle2_shift
          torsion_shift = pi/(180.d0/54.74d0)
          torsion_check_set = .true.! NR3 group
        end if
      end if
! extrapolation between tetragonal and planar NR3 group
      if (torsion_check_set) then
        torsion_check = torsion(hblist(i, 7), hblist(i, 6), hblist(i, 5), hblist(i, 8))  !torsion2_new
        if (torsion_check <= -pi) torsion_check = torsion_check + 2.d0 * pi
        if (torsion_check  >  pi) torsion_check = torsion_check - 2.d0 * pi
        if (torsion_check < 0) then
          torsion_check = -pi -torsion_check
        else
          torsion_check = pi -torsion_check
        end if
        torsion_check_bac = torsion_check ! save sign
        if (torsion_check < 0) torsion_check = -1.d0 * torsion_check
        torsion_check = torsion_check*180.d0/pi
        torsion_shift = torsion_shift + pi/(180.d0/((54.74d0 -torsion_check)/54.74d0*35.26d0))
        angle2_shift = angle2_shift -pi/(180.d0/((54.74d0 -torsion_check)/54.74d0*19.48d0))
        angle2_shift_2 = angle2_shift
        torsion_check = torsion_check_bac ! restore sign
      end if
! second angle
      angle2 = angle(hblist(i, 6), hblist(i, 5), hblist(i, 9))   ! angle2_new
!
!  Patch - as angle approaches 180 degrees, reduce the weight of this term.
!
!     weight_956 = 1.d0 - exp(-10.d0*(angle2 - pi)**2)
      angle2_cos_new = cos(angle2_shift -angle2)
      angle2_cos_2_new = cos(angle2_shift_2 -angle2)
      if (angle2_cos_2_new > angle2_cos_new) angle2_cos_new = angle2_cos_2_new
      if (angle2_cos_new <= 0.d0) return
! torsion angle
      torsion_correct = torsion(hblist(i, 7), hblist(i, 6), hblist(i, 5), hblist(i, 9))   ! torsion1_new
      if (torsion_correct <= -pi) torsion_correct = torsion_correct  +  2.d0 * pi
      if (torsion_correct > pi) torsion_correct = torsion_correct - 2.d0 * pi
      if ((.not.torsion_check_set2) .or. (abs(torsion_correct*180.d0/pi) > 90.d0)) then
        if (torsion_correct < 0.d0) then
          torsion_correct = -pi -torsion_correct
        else
          torsion_correct = pi - torsion_correct
        end if
      end if
! correction of NR3 torsion angle for through-bond case
      if (torsion_check < 0.d0) then ! negative torsion angle occupied by -NR3 r3
        torsion_value = torsion_shift - torsion_correct
        if (torsion_value <= -pi) torsion_value = torsion_value + 2.d0 * pi
        if (torsion_value  >  pi) torsion_value = torsion_value - 2.d0 * pi
        torsion_cos_new = cos(torsion_value)
      else if (torsion_check > 0) then ! positive torsion angle occupied by -NR3 r3
        torsion_value = -torsion_shift - torsion_correct
        if (torsion_value <= -pi) torsion_value = torsion_value + 2.d0 * pi
        if (torsion_value  >  pi) torsion_value = torsion_value - 2.d0 * pi
        torsion_cos_new = cos(torsion_value)
      else ! planar -NR3 or general case
        torsion_value   =  torsion_shift - torsion_correct
        torsion_value_2 = -torsion_shift - torsion_correct
        if (torsion_value   <= -pi) torsion_value   = torsion_value   + 2.d0 * pi
        if (torsion_value    >  pi) torsion_value   = torsion_value   - 2.d0 * pi
        if (torsion_value_2 <= -pi) torsion_value_2 = torsion_value_2 + 2.d0 * pi
        if (torsion_value_2  >  pi) torsion_value_2 = torsion_value_2 - 2.d0 * pi
        torsion_cos_new   = cos(torsion_value)
        torsion_cos_2_new = cos(torsion_value_2)
        if (torsion_cos_2_new > torsion_cos_new)torsion_cos_new = torsion_cos_2_new
      end if
      if (distance(hblist(i, 9), hblist(i, 5)) > distance(hblist(i, 9), hblist(i, 6)) .and. torsion_check_set2) &
        torsion_cos_new = 0.d0
      if (hblist(i, 7) == hblist(i, 8) .or. hblist(i, 9) == hblist(i, 8)) torsion_cos_new = 1.d0
   !   if (torsion_cos_new < 0.d0) return
      torsion_cos_new=abs(torsion_cos_new)
! parameters
      if (nat(hblist(i, 1)) == 7) then ! nitrogen
        if (nrbondsa(i) >= 3) then
          scale_a = scale_nsp3
        else
          scale_a = scale_nsp2
        end if
      else ! oxygen
        if (nrbondsa(i) >= 2) then
          scale_a = scale_osp3
        else
          scale_a = scale_osp2
        end if
      end if ! sulphur?
      if (nat(hblist(i, 5)) == 7) then ! nitrogen
        if (nrbondsb(i) >= 3) then
          scale_b = scale_nsp3
        else
          scale_b = scale_nsp2
        end if
      else ! oxygen
        if (nrbondsb(i) >= 2) then
          scale_b = scale_osp3
        else
          scale_b = scale_osp2
        end if
      end if ! sulphur?
      scale_c = (scale_a + scale_b)/2.d0
!
! damping 1.d0/(1.d0 + exp( -60.d0*(x/2.65d0 - 1.d0)))
!
      ha_dist = distance(hblist(i, 9), hblist(i, 1))
      hb_dist = distance(hblist(i, 9), hblist(i, 5))
      xc_dist = min(ha_dist, hb_dist)

! energies
      if (method_PM7) then
        XY_dist = max(ha_dist, hb_dist) - xc_dist
        if (XY_dist > 0.5d0) then
!
! Make damping a function of the shorter X - H distance.
! Distances greater than covcut = 1.2 Angstroms reduce damping
! distances less than covcut are not affected
!
          damping = 1.d0 - 1.d0/(1.d0 + exp( -60.d0*(xc_dist/covcut - 1.0d0)))
        else
!
! In very strong X - H - Y systems the X - Y distance is small, so an increased
! X - H minimum distance is expected, therefore do not dampen the function.
!
          damping = 1.d0
        end if
! y-x damping
        xc_dist = distance(hblist(i, 1), hblist(i, 5))

!
! Make damping a function of the O - O distance
! Damping is decreased as the O - O distances approaches shortcut Angstroms
!
        damping = damping/(1.d0 + exp( -100.d0*(xc_dist/shortcut - 1.d0)))
!
! At very large O - O distances, make damping go to zero
!
        damping = damping*(1.d0 - 1.d0/(1.d0 + exp( -10.d0*(xc_dist/longcut - 1.d0))))
        XY_dist = distance(hblist(i, 1), hblist(i, 5))
!
!   angle_cos       = cos(angle between donor-hydrogen-acceptor)
!   angle2_cos      = cos(difference between optimum DA-D-H and actual DA-D-H angle)
!   torsion_cos     =
!   angle2_cos_new  = cos(difference between optimum AA-A-H and actual AA-A-H angle)
!   torsion_cos_new = cos(difference between optimum AA'-AA-A-H and actual AA'-AA-A-H torsion angle)
!
        EH_plus = scale_c/XY_dist**2.d0*angle_cos**2* &
          (1.d0 - (1.d0 - angle2_cos*torsion_cos*angle2_cos_new*torsion_cos_new)**2)*hartree2kcal*damping
  !      (1.d0 - (1.d0 - weight_956*angle2_cos*torsion_cos*angle2_cos_new*torsion_cos_new)**2)*hartree2kcal*damping
        if (nat(hblist(i, 1)) == 8 .and. nat(hblist(i, 5)) == 8) then
!
!  Add in extra stabilization as the O - O distance approaches 2.41 Angstroms.
!  As a guide, formic acid dimer has O - O: 2.67 Angstroms, stabilization should total 18.61 kcal/mol
!              the H5O2(+) ion has O - O: 2.38 - 2.43 Angstroms
!
          short = -2.5d0*exp(-80.d0*(max(XY_dist - 2.67d0, 0.d0))**2)*angle_cos**4
          EH_plus = EH_plus + short
        end if
      else
! y-h damping
      damping = 1.d0 - 1.d0/(1.d0 + exp( -60.d0*(xc_dist/covcut - 1.d0)))
! y-x damping
      xc_dist = distance(hblist(i, 1), hblist(i, 5))

! short range
      damping = damping/(1.d0 + exp( -100.d0*(xc_dist/shortcut - 1.d0)))
! long range
      damping = damping*(1.d0 - 1.d0/(1.d0 + exp( -10.d0*(xc_dist/longcut - 1.d0))))
        EH_plus = scale_c/distance(hblist(i, 1), hblist(i, 5))**2.d0* &
        angle_cos**2*angle2_cos**2*torsion_cos**2*angle2_cos_new**2*torsion_cos_new**2*hartree2kcal*damping
      end if
    end if
    return
  end function EH_plus
