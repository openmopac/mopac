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

double precision function disp_DnX(l_grad)
   use common_arrays_C, only: nat, dxyz, cell_ijk, vab
   use molkst_C, only : numat, l123, l1u, l2u, l3u, method_PM6_DH2X, e_disp
   implicit none
   logical, intent (in) :: l_grad
!
!  Local variables
!
   integer :: i, j, k, l, iii, jjj, kkkk, i_cell, j_cell, iX(100), jOorN(100)
   double precision :: Rab, sum, sum2, sum3,  a_X(3,3), b_X(3,3)
   logical :: first = .true.
   logical, external :: connected
   double precision, external :: distance
   save
      if (first) then
          first = .false.
          if (method_PM6_DH2X) then
!
! Add in Rezac and Hobza's correction: "A halogen-bonding correction for the semiempirical PM6 method"
! Chem. Phys. Lett. 506 286-289 (2011)
!
            a_X(1,1) = 1.0489d12  ! Cl - N
            a_X(2,1) = 1.0226d5   ! Br - N
            a_X(3,1) = 1.2751d12  !  I - N
            a_X(1,2) = 4.6783d8   ! Cl - O
            a_X(2,2) = 9.6021d3   ! Br - O
            a_X(3,2) = 6.0912d5   !  I - O
            b_X(1,1) = -9.946d0   ! Cl - N
            b_X(2,1) = -3.236d0   ! Br - N
            b_X(3,1) = -9.534d0   !  I - N
            b_X(1,2) = -6.867d0   ! Cl - O
            b_X(2,2) = -2.900d0   ! Br - O
            b_X(3,2) = -4.154d0   !  I - O
          else
!
! Parameters for the "X" part in the "D3H4X" method
!
! Use Brahmkshatriya, et al.: "Quantum Mechanical Scoring: Structural and Energetic Insights into
! Cyclin-Dependent Kinase 2 Inhibition by Pyrazolo[1,5-a]pyrimidines" Current Computer-Aided Drug Design, 2013, 9, 118-129
! Table 2
            a_X(1,1) = 1.049d12   ! Cl - N
            a_X(2,1) = 5.560d4    ! Br - N
            a_X(3,1) = 5.237d8    !  I - N
            a_X(1,2) = 1.871d9    ! Cl - O
            a_X(2,2) = 2.160d4    ! Br - O
            a_X(3,2) = 2.436d6    !  I - O
            a_X(3,3) = 1.051d6    !  I - S
            b_X(1,1) = -9.95d0    ! Cl - N
            b_X(2,1) = -3.04d0    ! Br - N
            b_X(3,1) = -6.77d0    !  I - N
            b_X(1,2) = -7.44d0    ! Cl - O
            b_X(2,2) = -3.30d0    ! Br - O
            b_X(3,2) = -4.71d0    !  I - O
            b_X(3,3) = -3.82d0    !  I - S
          end if
          iX(17) = 1
          iX(35) = 2
          iX(53) = 3
          jOorN(7)  = 1
          jOorN(8)  = 2
          jOorN(16) = 3
      end if
      sum = 0.d0
      do i = 1, numat
        if (nat(i) /= 17 .and. nat(i) /= 35 .and. nat(i) /= 53) cycle ! crude, but fast
        k = nat(i)
        do j = 1, numat
          select case (nat(j))
          case (7, 8, 16)
            if (k /= 53 .and. nat(j) == 16) cycle  ! If sulfur, only select iodine
            l = nat(j)
            Rab = distance(i, j)
            sum2 = a_X(iX(k),jOorN(l))
            sum3 = b_X(iX(k),jOorN(l))
            sum = sum + sum2*exp(sum3*Rab)
            if (l_grad) then
              if (connected(i,j, 8.d0**2)) then
  !
  !   kkkk is the cell that atom j is in, relative to atom i
  !
                iii = l123*(i - 1)
                jjj = l123*(j - 1)
                kkkk = (l3u - cell_ijk(3)) + (2*l3u + 1)*(l2u - cell_ijk(2) + (2*l2u + 1)*(l1u - cell_ijk(1)))
                i_cell = iii + kkkk
                j_cell = jjj - kkkk
                do l = 1,3
                  dxyz(i_cell*3 + l) = dxyz(i_cell*3 + l) + Vab(l)*sum2*sum3*exp(sum3*Rab)/Rab
                  dxyz(j_cell*3 + l) = dxyz(j_cell*3 + l) - Vab(l)*sum2*sum3*exp(sum3*Rab)/Rab
                end do
              end if
            end if
          end select
        end do
      end do
      disp_DnX = sum
      e_disp = e_disp + sum
      return
  end function disp_DnX
  subroutine print_post_scf_corrections
  use molkst_C, only : keywrd, E_disp, E_hb, P_Hbonds
  use common_arrays_C, only: H_energy, H_txt
  use chanel_C, only : iw
  implicit none
  double precision :: sum, sum1
  double precision, external :: reada
  integer :: i, j, k
    if (index(keywrd," DISP(") > 0) then
      write(iw,'(/47x,a)')" List of hydrogen bonds found"
      write(iw,'(3x,a,12x,a,16x,a,11x,a,23x,a,17x,a)')"No.", "Donor", &
      "R(D-H)", "Hydrogen",  "Acceptor", "H-bond energy"
      sum1 = -abs(reada(keywrd, index(keywrd," DISP(") + 5))
      k = 0
      do
        sum = 0.d0
        j = 0
        do i = 1, P_hbonds
          if (sum > H_energy(i)) then
            sum = H_energy(i)
            j = i
          end if
        end do
        if (sum > sum1) exit
        k = k + 1
        write(iw,'(i5,3x,a)')k, trim(H_txt(j))
        H_energy(j) = 10.d0
      end do
    end if
    if (index(keywrd, "0SCF") /= 0) then
      write(iw,'(/10x,"DISPERSION ENERGY       =", f17.5, a)') e_disp, " KCAL/MOL"
      write(iw,'(10x,"H-BOND ENERGY           =", f17.5, a,/)') e_hb, " KCAL/MOL"
    end if
  end subroutine print_post_scf_corrections
