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

double precision function PM6_DH_Dispersion(l_grad)
  use molkst_C, only : numat, E_disp
  use common_arrays_C, only: coord, dxyz
  implicit none
  logical, intent (in) :: l_grad
  double precision, external :: PM6_DH_Disp
!
!  Local
!
  double precision :: sum, sum1, sum2, delta = 1.d-5
  integer :: i, k
  E_disp = PM6_DH_Disp(0, numat)
  PM6_DH_Dispersion = E_disp
  if (l_grad) then
!
!   Dispersion contribution is small, so only calculate the CUC.
!
    do k = 1, numat
      sum2 = PM6_DH_Disp(k, 1)
      do i = 1, 3
        coord(i,k) = coord(i,k) + delta
        sum = PM6_DH_Disp(k, 1)
        sum1 = (sum2 - sum)/delta
        if (Abs(sum1) < 50.d0) then
          dxyz((k - 1)*3 + i) = dxyz((k - 1)*3 + i) - sum1
        end if
        coord(i,k) = coord(i,k) - delta
      end do
    end do
  end if
  return
end function PM6_DH_Dispersion

double precision function PM6_DH_Disp(set_a, nsa)
!
!  Based on materials provided by Jan Rezak, as described in:
!
!  "A Transferable H-bonding Correction For Semiempirical Quantum-Chemical Methods"
!  Martin Korth, Michal Pitonak, Jan Rezac and Pavel Hobza,
!
!  J Chem Theory and Computation 6:344-352 (2010)
!
  use common_arrays_C, only: nat, nbonds, tvec, coord
  use molkst_C, only : Rij => Rab, method_PM7, numcal, line, id, l11, l21, l31, numat
  use chanel_C, only : iw
  use elemts_C, only : atom_names
  implicit none
  integer :: nsa
  integer :: set_a
  logical :: first_1 = .true., l_eles(110), first_2 = .true.
  double precision :: alpha = 20.d0, s = 1.04d0, cscale = 0.89d0, C6i, C6j, V_ab(3)
  integer :: ii, i, j, ni, nj, j_start, icalcn = 0, i1, j1, k1
  double precision :: C(86),R(86),N(86)
  double precision :: C6, R0, E_disp_tmp, damp, disp_limit = 6.5d0, sum, E_disp
  logical, external :: connected
  save :: icalcn, first_1, first_2

! C6 parameters (J.nm^6/mol)
  data C/ &
    0.16D0, 0.084D0,0.00D0, 0.00D0, 5.79D0, 1.65D0, 1.11D0, 0.70D0,   &  !  H He Li Be  B  C  N  O
    0.57D0, 0.45D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 3.25D0, 5.79D0,   &  !  F Ne Na Mg Al Si  P  S
    5.97D0, 3.71D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Cl Ar  K Ca Sc Ti  V Cr
    0.00D0, 0.00D0, 0.04D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Mn Fe Co Ni Cu Zn Ga Ge
    0.00D0, 0.00D0,11.60D0, 4.47D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! As Se Br Kr Rb Sr  Y Zr
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Nb Mo Tc Ru Rh Pd Ag Cd
    0.00D0, 0.00D0, 0.00D0, 0.00D0,25.80D0,16.50D0, 0.00D0, 0.00D0,   &  ! In Sn Sb Te  I Xe Cs Ba
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! La Ce Pr Nd Pm Sm Eu Gd
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Tb Dy Ho Er Tm Yb Lu Hf
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Ta  W Re Os Ir Pt Au Hg
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0/                      ! Tl Pb Bi Po At Rn

! R0 parameters (pm)
  data R/ &
    156.d0,  140.d0,  000.d0,  000.d0,  180.d0,  170.d0,  155.d0,  152.d0,   &  !  H He Li Be  B  C  N  O
    147.d0,  154.d0,  000.d0,  000.d0,  000.d0,  000.d0,  180.d0,  180.d0,   &  !  F Ne Na Mg Al Si  P  S
    175.d0,  188.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! Cl Ar  K Ca Sc Ti  V Cr
    000.d0,  000.d0,  140.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! Mn Fe Co Ni Cu Zn Ga Ge
    000.d0,  000.d0,  185.d0,  202.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! As Se Br Kr Rb Sr  Y Zr
    000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! Nb Mo Tc Ru Rh Pd Ag Cd
    000.d0,  000.d0,  000.d0,  000.d0,  198.d0,  216.d0,  000.d0,  000.d0,   &  ! In Sn Sb Te  I Xe Cs Ba
    000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! La Ce Pr Nd Pm Sm Eu Gd
    000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! Tb Dy Ho Er Tm Yb Lu Hf
    000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0,   &  ! Ta  W Re Os Ir Pt Au Hg
    000.d0,  000.d0,  000.d0,  000.d0,  000.d0,  000.d0/                        ! Tl Pb Bi Po At Rn

! Neff (electrons) - Slater-Kirkwood effective number of electrons
  data N/ &
    0.80D0, 1.42D0, 0.00D0, 0.00D0, 2.16D0, 2.50D0, 2.82D0, 3.15D0,   &  !  H He Li Be  B  C  N  O
    3.48D0, 3.81D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 4.50D0, 4.80D0,   &  !  F Ne Na Mg Al Si  P  S
    5.10D0, 5.40D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Cl Ar  K Ca Sc Ti  V Cr
    0.00D0, 0.00D0, 2.90D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Mn Fe Co Ni Cu Zn Ga Ge
    0.00D0, 0.00D0, 6.00D0, 6.30D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! As Se Br Kr Rb Sr  Y Zr
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Nb Mo Tc Ru Rh Pd Ag Cd
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 6.95D0, 7.25D0, 0.00D0, 0.00D0,   &  ! In Sn Sb Te  I Xe Cs Ba
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! La Ce Pr Nd Pm Sm Eu Gd
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Tb Dy Ho Er Tm Yb Lu Hf
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0,   &  ! Ta  W Re Os Ir Pt Au Hg
    0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0, 0.00D0/                      ! Tl Pb Bi Po At Rn
!
! dispersion energy
!
  E_disp = 0.d0
  if (method_PM7) then
    alpha  = 15.450118D0
    s      = 1.226593D0
    cscale = 2.286419D0
  else
    alpha  = 20.d0
    s      = 1.04d0
    cscale = 0.89d0
  end if
  l_eles = .true.
  do ii = 1, nsa
    if (nsa == numat) then
      i = ii
    else
      i = set_a
    end if
    ni = nat(i)
    if (ni > 86) cycle
    if (R(ni) == 0.d0 .or. C(ni) == 0.d0 .or. N(ni) == 0.d0) then
      if (method_PM7) cycle
      if (icalcn /= numcal .and. l_eles(ni)) then
        l_eles(ni) = .false.
        line = atom_names(ni)
        do
          if (line(1:1) /= " ") exit
          line(1:) = line(2:)
        end do
        if (first_1) write(iw,'(///,6(a,/))') "          *********************************************", &
                                            "          *                                           *", &
                                            "          *                 WARNING                   *", &
                                            "          *                                           *", &
                                            "          *********************************************"
        if (first_2) write(iw,'(10x,a)') "Dispersion parameters missing for "//trim(line)
      end if
      first_1 = .false.
      cycle
    end if
    if (nsa == numat) then
      j_start = ii + 1
    else
      j_start = 1
    end if
    do j  = j_start, numat
      if (i == j) cycle
      nj = nat(j)
      if (nj > 86) cycle
      if (ni == 6) then
        if (nbonds(i) == 4) then
          C6i = 0.95d0
        else
          C6i = 1.65d0
        end if
      else
        C6i = C(ni)
      end if
      if (nj == 6) then
        if (nbonds(j) == 4) then
          C6j = 0.95d0
        else
          C6j = 1.65d0
        end if
      else
        C6j = C(nj)
      end if
      if (R(nj).eq.0.0.or.C(nj).eq.0.0 .or. N(nj).eq.0.0) cycle
      C6 = 2.d0*(C6i**2*C6j**2*N(ni)*N(nj))**(1.d0/3.d0)/((C6i*N(nj)**2)**(1.d0/3.d0) + &
        (C6j*N(ni)**2)**(1.d0/3.d0))
      R0 = (R(ni)**3+R(nj)**3)/(R(ni)**2+R(nj)**2) /1000.d0*2.d0    ! pm to nm
      if (connected(i, j, 100.d0**2)) then
        if (id == 0) then
! atomic distance Rij, in nm
! i-j dispersion energy
          Rij = Rij*0.1d0
          damp = 1.d0 / (1.d0 + exp(-alpha * (Rij/(s*R0) - 1.d0)))
          E_disp_tmp = C6/Rij**6*damp/(1000.d0*4.184d0)
          E_disp = E_disp - E_disp_tmp
        else
          sum = 0.d0
          do i1 = -l11, l11
            do j1 = -l21, l21
              do k1 = -l31, l31
                V_ab = coord(:,i) - coord(:,j) + tvec(:,1)*i1 + tvec(:,2)*j1 + tvec(:,3)*k1
                Rij = sqrt(V_ab(1)**2 + V_ab(2)**2 + V_ab(3)**2)*0.1d0
                if (Rij*10.d0 < disp_limit) then
                  damp = 1.d0 / (1.d0 + exp(-alpha * (Rij/(s*R0) - 1.d0)))
!
!  The damping function *(1.d0 - exp(-(Rij*10.d0 - disp_limit)**2)) is used in solids because the number
!  of atoms increases as Rij**2.  This makes long-range interactions (6 Angstroms to infinity) more important.
!  In molecules, there is no equivalent increase.
!

                  E_disp_tmp = C6/Rij**6*damp/(1000.d0*4.184d0)*(1.d0 - exp(-(Rij*10.d0 - disp_limit)**2))
                    sum = sum + E_disp_tmp
                end if
              end do
            end do
          end do
          E_disp = E_disp - sum
        end if
      end if
    end do
  end do
  if (first_2 .and. .not. first_1) &
  write(iw,'(/,6(a,/))')   "          *********************************************", &
                            "          *                                           *", &
                            "          *              END OF WARNING               *", &
                            "          *                                           *", &
                            "          *********************************************"
  first_2 = .false.
  icalcn = numcal
  PM6_DH_Disp = E_disp*cscale
  return
end function PM6_DH_Disp
