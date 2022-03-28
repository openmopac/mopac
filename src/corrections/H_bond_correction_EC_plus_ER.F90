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

  double precision function EC_plus_ER(D, H, A, q1, q2, EC, ER, d_list, nd_list)
!
! Based on "A Transferable H-bonding Correction For Semiempirical Quantum-Chemical Methods",
! by Martin Korth, Michal Pitonak, Jan Rezac and Pavel Hobza, J. Chem. Theory Comput., 2010, 6 (1), pp 344:352
!
  use common_arrays_C, only: coord, nat, nbonds, ibonds
  use molkst_C, only : Rab
  use parameters_C, only : dh2_a_parameters
  use funcon_C, only : pi, a0, eV, fpc_9
  implicit none
  integer :: D, H, A, R1, R2, R3, C_of_CO, i, nH, nC, nO, nN, N_of_HNCO, ii, nd_list, d_list(8)
  double precision :: q1, q2, EC, ER
  double precision :: angle_cos, r, torsion_check, angle, angle2_shift, &
    angle2_shift_2, torsion_shift, angle2_cos, angle2, angle2_cos_2, &
    torsion, torsion_cos, torsion2, torsion_2, torsion_ref, &
    torsion_cos_2, unit_part, attraction, torsion_correct, &
    expo, repulsion, rep_pre, rep_exp, c, sum, multiplier_a, sum_max
  logical :: torsion_check_set, peptide, torsion_check_set2, first = .true.
  double precision, external :: truncation
  logical, external :: connected
  save :: first
  if (first) then
    first = .false.
    dh2_a_parameters(1) = 1.48d0    ! Nitrogen
    dh2_a_parameters(2) = 1.56d0    ! Oxygen, generic
    dh2_a_parameters(3) = 1.55d0    ! Oxygen, acid
    dh2_a_parameters(4) = 0.96d0    ! Oxygen, peptide
    dh2_a_parameters(5) = 0.76d0    ! Oxygen, water
    dh2_a_parameters(6) = 0.85d0    ! Sulfur
  end if
  EC_plus_ER = 0.d0
  EC = 0.d0
  ER = 0.d0
  torsion_shift = 0.d0
  multiplier_a = 0.d0
  angle2_shift = 0.d0
  angle2_shift_2 = 0.d0
!
!  D is the atom that hydrogen atom H is attached to, and A is the distant acceptor
!
!
!  First angle
!
    call bangle (coord, D, H, A, angle)
    d_list(1) = D
    d_list(2) = H
    d_list(3) = A
    nd_list = 3
    angle_cos = -cos(angle)
    if (angle_cos < 0.d0) return
!
! get target angles
!
    torsion_check = 0.0
    torsion_check_set = .false.
    torsion_check_set2 = .false.
    if (nat(A) == 8 .or. nat(A) == 16) then
      if (nbonds(A) == 1) then
        angle2_shift = pi
        angle2_shift_2 = pi/(180.d0/120.d0)
        torsion_shift = 0.0
        torsion_check_set2 = .true.
      else
        angle2_shift = pi/(180.d0/109.48d0)
        angle2_shift_2 = angle2_shift
        torsion_shift = pi/(180.d0/54.74d0)
      end if
      if (nat(A) == 8) then
        multiplier_a = 0.d0
        if (nbonds(A) == 2) then
          if (nat(ibonds(1,A)) == 1 .and. nat(ibonds(2,A)) == 1) multiplier_a = dh2_a_parameters(5)   ! Water
        else if (nat(ibonds(1,A)) == 6) then
          C_of_CO =  ibonds(1,A)
          if (nbonds(C_of_CO) == 3) then
            nH = 0
            nC = 0
            nO = 0
            nN = 0
            do i = 1, 3
              if (nat(ibonds(i,C_of_CO)) == 1) nH = nH +1
              if (nat(ibonds(i,C_of_CO)) == 6) nC = nC +1
              if (nat(ibonds(i,C_of_CO)) == 7) then
                nN = nN +1
                N_of_HNCO = ibonds(i,C_of_CO)
              end if
              if (nat(ibonds(i,C_of_CO)) == 8) nO = nO +1
            end do
            if (nC + nH == 1 .and. nN == 1 .and. nO == 1) then
              peptide = .false.
!
!  Check that H-N-C-O exists and is trans
!
              nH = 0
              do i = 1, nbonds(N_of_HNCO)
                if (nat(ibonds(i,N_of_HNCO)) == 1) then
                  nH = nH + 1
                  call dihed (coord, A, C_of_CO, N_of_HNCO, ibonds(i,N_of_HNCO), sum)
                  sum = min(sum, 2*pi - sum)
                  if (sum > 0.5d0*pi) peptide = .true.
                end if
              end do
              if (peptide .and. nH > 0) multiplier_a = dh2_a_parameters(4)   ! peptide oxygen
            else if (nO == 2) then
              multiplier_a = dh2_a_parameters(3)   ! acid (carboxylate) oxygen
            end if
          end if
        end if
        if (multiplier_a < 1.d-20) multiplier_a  = dh2_a_parameters(2)   ! generic oxygen
      else
        multiplier_a = dh2_a_parameters(6)   ! sulfur
      end if
    else if (nat(A) == 7) then
      multiplier_a = dh2_a_parameters(1)   ! nitrogen
      if (nbonds(A) == 2) then
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
!
!  extrapolation between tetrahedral and planar NR3 group
!
    R1 = 0
    if (nbonds(A) == 1) then
      R1 = ibonds(1,A)
      sum_max = 0.d0
      do ii = 1, nbonds(R1)
        i = ibonds(ii,R1)
!  if (i == A) cycle
       if ( .not. connected(i, H, 1000.d0**2)) return
        if (sum_max < Rab) then
          R2 = i
          sum_max = Rab
        end if
      end do
      R3 = H
    else if (nbonds(A) == 2) then
      sum_max = 0.d0
      do ii = 1, nbonds(A)
        i = ibonds(ii,A)
        if ( .not. connected(i, H, 1000.d0**2)) return
        if (sum_max < Rab) then
          sum_max = Rab
          R1 = i
        end if
      end do
      do ii = 1, nbonds(A)
        i = ibonds(ii,A)
        if (i == R1) cycle
        R2 = i
      end do
      R3 = H
    else if (nbonds(A) == 3) then
      sum_max = 0.d0
      do ii = 1, nbonds(A)
        i = ibonds(ii,A)
        if ( .not. connected(i, H, 1000.d0**2)) return
        if (sum_max < Rab) then
          sum_max = Rab
          R1 = i
        end if
      end do
      sum_max = 0.d0
      do ii = 1, nbonds(A)
        i = ibonds(ii,A)
        if (i == R1) cycle
        if ( .not. connected(i, H, 1000.d0**2)) return
        if (sum_max < Rab) then
          sum_max = Rab
          R2 = i
        end if
      end do
      do ii = 1, nbonds(A)
        i = ibonds(ii,A)
        if (i == R1 .or. i == R2) cycle
        R3 = i
      end do
      if (torsion_check_set) call dihed (coord, R2, R1, A, R3, torsion2)
      torsion_check = torsion2
      if (torsion_check <   -pi) torsion_check = torsion_check + 2.d0*pi
      if (torsion_check >    pi) torsion_check = torsion_check - 2.d0*pi
      if (torsion_check < 0.d0) then
        torsion_check = -pi - torsion_check
      else
        torsion_check =  pi - torsion_check
      end if
      sum = torsion_check
      if (sum <  0.d0) sum = -sum
      sum = 180.d0/pi*sum
      torsion_shift = torsion_shift + pi/(180.d0/((54.74d0 - sum)/54.74d0*35.26d0))
      angle2_shift=angle2_shift - pi/(180.d00/((54.74d0 - sum)/54.74d0*19.48d0))
      angle2_shift_2 = angle2_shift
    end if
    if (R1 == 0) return
!
! second angle
!
    call bangle (coord, R1, A, H, angle2)
    nd_list = nd_list + 1
    d_list(nd_list) = R1
    angle2_cos = cos(angle2_shift - angle2)
    angle2_cos_2 = cos(angle2_shift_2 - angle2)
    if (angle2_cos_2 > angle2_cos) angle2_cos = angle2_cos_2
    if (angle2_cos <= 0.d0) return
!
! torsion angle
! correction of NR3 torsion angle for through-bond case
!
      call dihed (coord, R2, R1, A, H, torsion_ref)
      nd_list = nd_list + 1
      d_list(nd_list) = R2
      torsion_correct = torsion_ref
      if (torsion_correct < -pi) torsion_correct = torsion_correct + 2.d0*pi
      if (torsion_correct >  pi) torsion_correct = torsion_correct - 2.d0*pi
      if (.not. torsion_check_set2 .or. Abs(torsion_correct) > 0.5d0*pi) then
        if (torsion_correct < 0.d0) then
          torsion_correct = -pi - torsion_correct
        else
          torsion_correct =  pi - torsion_correct
        end if
      end if


      if (torsion_check < 0.0) then ! negative torsion angle occupied by -NR3 r3
        torsion = torsion_shift - torsion_correct
        if (torsion <  -pi) torsion = torsion + 2.d0*pi
        if (torsion >   pi) torsion = torsion - 2.d0*pi
        torsion_cos = cos(torsion)
      else if (torsion_check > 0.0) then ! positive torsion angle occupied by -NR3 r3
        torsion = -torsion_shift - torsion_correct
        if (torsion <  -pi) torsion = torsion + 2.d0*pi
        if (torsion >   pi) torsion = torsion - 2.d0*pi
        torsion_cos = cos(torsion)
      else ! planar -NR3 or general case
        torsion = torsion_shift - torsion_correct
        torsion_2 = -torsion_shift - torsion_correct
        if (torsion <  -pi) torsion = torsion + 2.d0*pi
        if (torsion >   pi) torsion = torsion - 2.d0*pi
        if (torsion_2 <  -pi) torsion_2 = torsion_2 + 2.d0*pi
        if (torsion_2 >   pi) torsion_2 = torsion_2 - 2.d0*pi
        torsion_cos = cos(torsion)
        torsion_cos_2 = cos(torsion_2)
        if (torsion_cos_2 > torsion_cos) torsion_cos = torsion_cos_2
      end if
      if ( .not. connected(A, H, 1000.d0**2)) return
      r = Rab
      if ( .not. connected(R1, H, 1000.d0**2)) return
      if (torsion_check_set2 .and. r > Rab) torsion_cos = -1.d0
      if (torsion_cos <= 0.0) return
!
! distance cutoff as ad-hoc solution ...
!
!  To prevent a discontinuity in the gradients, make r a function of distance:
!
!  Above 1.85 Angstroms, r = r
!  Below 1.75 Angstroms r is increased to 1.80 Angstroms
!
    r = truncation(r, 1.80d0, 0.05d0)
    r = r/a0
    expo    = 3.d0    ! "b" in equation 2
    rep_pre = 0.65d0    ! "c" in equation 2
    rep_exp = 5.d0    ! "d" in equation 2
    unit_part = fpc_9*eV
    attraction = multiplier_a*q1*q2/r**expo*unit_part
    repulsion = rep_pre*rep_exp**(-r)*unit_part
    c = (attraction + repulsion)*angle_cos*angle2_cos*torsion_cos
    EC =       attraction*angle_cos*angle2_cos*torsion_cos
    ER =       repulsion*angle_cos*angle2_cos*torsion_cos
    EC_plus_ER = c
  end function EC_plus_ER
