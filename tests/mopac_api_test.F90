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

program mopac_api_test
    use mopac_api_f
    implicit none
    integer :: nfail
    nfail = 0

    call test_mopac_scf1(nfail)
    write(*,*) "Passed API test 1"
    call test_mopac_relax1(nfail)
    write(*,*) "Passed API test 2"
    call test_mopac_vibe1(nfail)
    write(*,*) "Passed API test 3"
    call test_mozyme_scf1(nfail)
    write(*,*) "Passed API test 4"
    call test_mozyme_relax1(nfail)
    write(*,*) "Passed API test 5"
    call test_mozyme_vibe1(nfail)
    write(*,*) "Passed API test 6"
    call test_mopac_restart1(nfail)
    write(*,*) "Passed API test 7"
    call test_mozyme_restart1(nfail)
    write(*,*) "Passed API test 8"
    call test_cosmo1(nfail)
    write(*,*) "Passed API test 9"
    call test_crystal1(nfail)
    write(*,*) "Passed API test 10"
    call exit(nfail)
end program mopac_api_test

subroutine test_mopac_scf1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O SCF'

    ! SCF calculation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.76d0
    test_in%coord(2) = 0.59d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.76d0
    test_in%coord(5) = 0.59d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.76975d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 2.307865d0
    test_target%coord_deriv(2) = 2.742432d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -2.307865d0
    test_target%coord_deriv(5) = 2.711610d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.0d0
    test_target%coord_deriv(8) = -5.454042d0
    test_target%coord_deriv(9) = 0.0d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.322260d0
    test_target%charge(2) = 0.322260d0
    test_target%charge(3) = -0.644520d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.147d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.896d0
    test_target%bond_order(2) = 0.895d0
    test_target%bond_order(3) = 0.896d0
    test_target%bond_order(4) = 0.895d0
    test_target%bond_order(5) = 0.895d0
    test_target%bond_order(6) = 0.895d0
    test_target%bond_order(7) = 1.791d0
    test_target%nerror = 0

    call mopac_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_scf1

subroutine test_mopac_relax1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O relax'

    ! 1 - geometry relaxation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.8d0
    test_in%coord(2) = 0.4d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.8d0
    test_in%coord(5) = 0.4d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.79952d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update(1) = 0.758862835d0
    test_target%coord_update(2) = 0.486805290d0
    test_target%coord_update(3) = 0.0d0
    test_target%coord_update(4) = -0.759709263d0
    test_target%coord_update(5) = 0.485142519d0
    test_target%coord_update(6) = 0.0d0
    test_target%coord_update(7) = 0.000846434d0
    test_target%coord_update(8) = -0.093004709d0
    test_target%coord_update(9) = 0.0d0
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = -0.565914d0
    test_target%coord_deriv(2) = -0.294814d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = 0.063505d0
    test_target%coord_deriv(5) = 0.058243d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.502409d0
    test_target%coord_deriv(8) = 0.236571d0
    test_target%coord_deriv(9) = 0.0d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.324544d0
    test_target%charge(2) = 0.324720d0
    test_target%charge(3) = -0.649264d0
    test_target%dipole(1) = -0.005d0
    test_target%dipole(2) = 2.129d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.895d0
    test_target%bond_order(2) = 0.894d0
    test_target%bond_order(3) = 0.895d0
    test_target%bond_order(4) = 0.894d0
    test_target%bond_order(5) = 0.894d0
    test_target%bond_order(6) = 0.894d0
    test_target%bond_order(7) = 1.788d0
    test_target%nerror = 0

    call mopac_relax_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_relax1

subroutine test_mopac_vibe1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    test_name = 'H2O vibe'

    ! 1 - geometry relaxation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.759568602d0
    test_in%coord(2) = 0.486889564d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.760025949d0
    test_in%coord(5) = 0.485109057d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.000457356d0
    test_in%coord(8) = -0.093055521d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.79986d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 0.000042d0
    test_target%coord_deriv(2) = 0.013515d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -0.006056d0
    test_target%coord_deriv(5) = -0.013157d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.006014d0
    test_target%coord_deriv(8) = -0.000359d0
    test_target%coord_deriv(9) = 0.0d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.324672d0
    test_target%charge(2) = 0.324697d0
    test_target%charge(3) = -0.649369d0
    test_target%dipole(1) = -0.003d0
    test_target%dipole(2) = 2.129d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.895d0
    test_target%bond_order(2) = 0.894d0
    test_target%bond_order(3) = 0.895d0
    test_target%bond_order(4) = 0.894d0
    test_target%bond_order(5) = 0.894d0
    test_target%bond_order(6) = 0.894d0
    test_target%bond_order(7) = 1.788d0
    test_target%nerror = 0
    allocate(test_target%freq(9))
    test_target%freq(1) = -18.0d0
    test_target%freq(2) = -14.6d0
    test_target%freq(3) = -4.0d0
    test_target%freq(4) = 0.0d0
    test_target%freq(5) = 2.3d0
    test_target%freq(6) = 4.7d0
    test_target%freq(7) = 1395.0d0
    test_target%freq(8) = 2810.2d0
    test_target%freq(9) = 2856.4d0
    allocate(test_target%disp(9,9))
    test_target%disp(1,1) = 0.1312d0
    test_target%disp(2,1) = -0.1922d0
    test_target%disp(3,1) = 0.8602d0
    test_target%disp(4,1) = 0.1308d0
    test_target%disp(5,1) = 0.1926d0
    test_target%disp(6,1) = 0.3710d0
    test_target%disp(7,1) = -0.0624d0
    test_target%disp(8,1) = 0.0000d0
    test_target%disp(9,1) = -0.1011d0
    test_target%disp(1,2) = -0.2659d0
    test_target%disp(2,2) = 0.4129d0
    test_target%disp(3,2) = -0.0185d0
    test_target%disp(4,2) = -0.2649d0
    test_target%disp(5,2) = -0.4137d0
    test_target%disp(6,2) = 0.6915d0
    test_target%disp(7,2) = 0.1976d0
    test_target%disp(8,2) = -0.0001d0
    test_target%disp(9,2) = -0.0029d0
    test_target%disp(1,3) = 0.3319d0
    test_target%disp(2,3) = -0.3100d0
    test_target%disp(3,3) = -0.4197d0
    test_target%disp(4,3) = 0.3312d0
    test_target%disp(5,3) = 0.3120d0
    test_target%disp(6,3) = 0.5078d0
    test_target%disp(7,3) = 0.3768d0
    test_target%disp(8,3) = 0.0028d0
    test_target%disp(9,3) = 0.1034d0
    test_target%disp(1,4) = -0.0016d0
    test_target%disp(2,4) = 0.2370d0
    test_target%disp(3,4) = 0.0003d0
    test_target%disp(4,4) = -0.0016d0
    test_target%disp(5,4) = 0.2361d0
    test_target%disp(6,4) = 0.0001d0
    test_target%disp(7,4) = -0.0048d0
    test_target%disp(8,4) = 0.9424d0
    test_target%disp(9,4) = 0.0007d0
    test_target%disp(1,5) = 0.0966d0
    test_target%disp(2,5) = 0.1446d0
    test_target%disp(3,5) = 0.1956d0
    test_target%disp(4,5) = 0.0970d0
    test_target%disp(5,5) = -0.1428d0
    test_target%disp(6,5) = -0.3284d0
    test_target%disp(7,5) = 0.8218d0
    test_target%disp(8,5) = 0.0041d0
    test_target%disp(9,5) = -0.3442d0
    test_target%disp(1,6) = 0.0123d0
    test_target%disp(2,6) = 0.0684d0
    test_target%disp(3,6) = 0.2130d0
    test_target%disp(4,6) = 0.0125d0
    test_target%disp(5,6) = -0.0682d0
    test_target%disp(6,6) = -0.1359d0
    test_target%disp(7,6) = 0.2567d0
    test_target%disp(8,6) = 0.0005d0
    test_target%disp(9,6) = 0.9277d0
    test_target%disp(1,7) = -0.3143d0
    test_target%disp(2,7) = 0.5969d0
    test_target%disp(3,7) = -0.0000d0
    test_target%disp(4,7) = 0.3129d0
    test_target%disp(5,7) = 0.5976d0
    test_target%disp(6,7) = -0.0000d0
    test_target%disp(7,7) = 0.0003d0
    test_target%disp(8,7) = -0.2998d0
    test_target%disp(9,7) = -0.0000d0
    test_target%disp(1,8) = 0.5406d0
    test_target%disp(2,8) = 0.4131d0
    test_target%disp(3,8) = -0.0000d0
    test_target%disp(4,8) = 0.5419d0
    test_target%disp(5,8) = -0.4119d0
    test_target%disp(6,8) = -0.0000d0
    test_target%disp(7,8) = -0.2717d0
    test_target%disp(8,8) = -0.0003d0
    test_target%disp(9,8) = 0.0000d0
    test_target%disp(1,9) = -0.6335d0
    test_target%disp(2,9) = -0.2964d0
    test_target%disp(3,9) = -0.0000d0
    test_target%disp(4,9) = 0.6340d0
    test_target%disp(5,9) = -0.2947d0
    test_target%disp(6,9) = -0.0000d0
    test_target%disp(7,9) = -0.0001d0
    test_target%disp(8,9) = 0.1483d0
    test_target%disp(9,9) = 0.0000d0

    call mopac_vibe_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_vibe1

subroutine test_mozyme_scf1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mozyme_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O SCF MOZYME'

    ! SCF calculation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.76d0
    test_in%coord(2) = 0.59d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.76d0
    test_in%coord(5) = 0.59d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.76837d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 2.237d0
    test_target%coord_deriv(2) = 2.318d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -2.222d0
    test_target%coord_deriv(5) = 2.295d0
    test_target%coord_deriv(6) = 0.001d0
    test_target%coord_deriv(7) = -0.014d0
    test_target%coord_deriv(8) = -4.614d0
    test_target%coord_deriv(9) = -0.001d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.320470d0
    test_target%charge(2) = 0.320480d0
    test_target%charge(3) = -0.640950d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.137d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.897d0
    test_target%bond_order(2) = 0.896d0
    test_target%bond_order(3) = 0.897d0
    test_target%bond_order(4) = 0.896d0
    test_target%bond_order(5) = 0.896d0
    test_target%bond_order(6) = 0.896d0
    test_target%bond_order(7) = 1.793d0
    test_target%nerror = 0

    call mozyme_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_scf1

subroutine test_mozyme_relax1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mozyme_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O relax MOZYME'

    ! 1 - geometry relaxation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.8d0
    test_in%coord(2) = 0.4d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.8d0
    test_in%coord(5) = 0.4d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.79952d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update(1) = 0.760052823d0
    test_target%coord_update(2) = 0.459205345d0
    test_target%coord_update(3) = 0.000297581d0
    test_target%coord_update(4) = -0.760101400d0
    test_target%coord_update(5) = 0.459419431d0
    test_target%coord_update(6) = 0.000222660d0
    test_target%coord_update(7) = 0.000048582d0
    test_target%coord_update(8) = -0.118624776d0
    test_target%coord_update(9) = -0.000520243d0
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = -0.103d0
    test_target%coord_deriv(2) = -0.132d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -0.027d0
    test_target%coord_deriv(5) = -0.08d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.13d0
    test_target%coord_deriv(8) = 0.212d0
    test_target%coord_deriv(9) = 0.001d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.325153d0
    test_target%charge(2) = 0.325224d0
    test_target%charge(3) = -0.650376d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.129d0
    test_target%dipole(3) = 0.003d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.895d0
    test_target%bond_order(2) = 0.894d0
    test_target%bond_order(3) = 0.895d0
    test_target%bond_order(4) = 0.894d0
    test_target%bond_order(5) = 0.894d0
    test_target%bond_order(6) = 0.894d0
    test_target%bond_order(7) = 1.788d0
    test_target%nerror = 0

    call mozyme_relax_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_relax1


subroutine test_mozyme_vibe1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mozyme_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    test_name = 'H2O vibe MOZYME'

    ! 1 - geometry relaxation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.759568602d0
    test_in%coord(2) = 0.486889564d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.760025949d0
    test_in%coord(5) = 0.485109057d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.000457356d0
    test_in%coord(8) = -0.093055521d0
    test_in%coord(9) = 0.0d0
    test_in%tolerance = 0.001d0
    test_target%heat = -57.79986d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 0.000042d0
    test_target%coord_deriv(2) = 0.013515d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -0.006056d0
    test_target%coord_deriv(5) = -0.013157d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.006014d0
    test_target%coord_deriv(8) = -0.000359d0
    test_target%coord_deriv(9) = 0.0d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.324672d0
    test_target%charge(2) = 0.324697d0
    test_target%charge(3) = -0.649369d0
    test_target%dipole(1) = -0.003d0
    test_target%dipole(2) = 2.129d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.895d0
    test_target%bond_order(2) = 0.894d0
    test_target%bond_order(3) = 0.895d0
    test_target%bond_order(4) = 0.894d0
    test_target%bond_order(5) = 0.894d0
    test_target%bond_order(6) = 0.894d0
    test_target%bond_order(7) = 1.788d0
    test_target%nerror = 0
    allocate(test_target%freq(9))
    test_target%freq(1) = -18.0d0
    test_target%freq(2) = -14.6d0
    test_target%freq(3) = -4.0d0
    test_target%freq(4) = 0.0d0
    test_target%freq(5) = 2.3d0
    test_target%freq(6) = 4.7d0
    test_target%freq(7) = 1395.0d0
    test_target%freq(8) = 2810.2d0
    test_target%freq(9) = 2856.4d0
    allocate(test_target%disp(9,9))
    test_target%disp(1,1) = 0.1312d0
    test_target%disp(2,1) = -0.1922d0
    test_target%disp(3,1) = 0.8602d0
    test_target%disp(4,1) = 0.1308d0
    test_target%disp(5,1) = 0.1926d0
    test_target%disp(6,1) = 0.3710d0
    test_target%disp(7,1) = -0.0624d0
    test_target%disp(8,1) = 0.0000d0
    test_target%disp(9,1) = -0.1011d0
    test_target%disp(1,2) = -0.2659d0
    test_target%disp(2,2) = 0.4129d0
    test_target%disp(3,2) = -0.0185d0
    test_target%disp(4,2) = -0.2649d0
    test_target%disp(5,2) = -0.4137d0
    test_target%disp(6,2) = 0.6915d0
    test_target%disp(7,2) = 0.1976d0
    test_target%disp(8,2) = -0.0001d0
    test_target%disp(9,2) = -0.0029d0
    test_target%disp(1,3) = 0.3319d0
    test_target%disp(2,3) = -0.3100d0
    test_target%disp(3,3) = -0.4197d0
    test_target%disp(4,3) = 0.3312d0
    test_target%disp(5,3) = 0.3120d0
    test_target%disp(6,3) = 0.5078d0
    test_target%disp(7,3) = 0.3768d0
    test_target%disp(8,3) = 0.0028d0
    test_target%disp(9,3) = 0.1034d0
    test_target%disp(1,4) = -0.0016d0
    test_target%disp(2,4) = 0.2370d0
    test_target%disp(3,4) = 0.0003d0
    test_target%disp(4,4) = -0.0016d0
    test_target%disp(5,4) = 0.2361d0
    test_target%disp(6,4) = 0.0001d0
    test_target%disp(7,4) = -0.0048d0
    test_target%disp(8,4) = 0.9424d0
    test_target%disp(9,4) = 0.0007d0
    test_target%disp(1,5) = 0.0966d0
    test_target%disp(2,5) = 0.1446d0
    test_target%disp(3,5) = 0.1956d0
    test_target%disp(4,5) = 0.0970d0
    test_target%disp(5,5) = -0.1428d0
    test_target%disp(6,5) = -0.3284d0
    test_target%disp(7,5) = 0.8218d0
    test_target%disp(8,5) = 0.0041d0
    test_target%disp(9,5) = -0.3442d0
    test_target%disp(1,6) = 0.0123d0
    test_target%disp(2,6) = 0.0684d0
    test_target%disp(3,6) = 0.2130d0
    test_target%disp(4,6) = 0.0125d0
    test_target%disp(5,6) = -0.0682d0
    test_target%disp(6,6) = -0.1359d0
    test_target%disp(7,6) = 0.2567d0
    test_target%disp(8,6) = 0.0005d0
    test_target%disp(9,6) = 0.9277d0
    test_target%disp(1,7) = -0.3143d0
    test_target%disp(2,7) = 0.5969d0
    test_target%disp(3,7) = -0.0000d0
    test_target%disp(4,7) = 0.3129d0
    test_target%disp(5,7) = 0.5976d0
    test_target%disp(6,7) = -0.0000d0
    test_target%disp(7,7) = 0.0003d0
    test_target%disp(8,7) = -0.2998d0
    test_target%disp(9,7) = -0.0000d0
    test_target%disp(1,8) = 0.5406d0
    test_target%disp(2,8) = 0.4131d0
    test_target%disp(3,8) = -0.0000d0
    test_target%disp(4,8) = 0.5419d0
    test_target%disp(5,8) = -0.4119d0
    test_target%disp(6,8) = -0.0000d0
    test_target%disp(7,8) = -0.2717d0
    test_target%disp(8,8) = -0.0003d0
    test_target%disp(9,8) = 0.0000d0
    test_target%disp(1,9) = -0.6335d0
    test_target%disp(2,9) = -0.2964d0
    test_target%disp(3,9) = -0.0000d0
    test_target%disp(4,9) = 0.6340d0
    test_target%disp(5,9) = -0.2947d0
    test_target%disp(6,9) = -0.0000d0
    test_target%disp(7,9) = -0.0001d0
    test_target%disp(8,9) = 0.1483d0
    test_target%disp(9,9) = 0.0000d0

    call mozyme_vibe_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_vibe1

subroutine test_mopac_restart1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O SCF restart'

    ! SCF calculation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.76d0
    test_in%coord(2) = 0.59d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.76d0
    test_in%coord(5) = 0.59d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.76975d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 2.307865d0
    test_target%coord_deriv(2) = 2.742432d0
    test_target%coord_deriv(3) = 0.0d0
    test_target%coord_deriv(4) = -2.307865d0
    test_target%coord_deriv(5) = 2.711610d0
    test_target%coord_deriv(6) = 0.0d0
    test_target%coord_deriv(7) = 0.0d0
    test_target%coord_deriv(8) = -5.454042d0
    test_target%coord_deriv(9) = 0.0d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.322260d0
    test_target%charge(2) = 0.322260d0
    test_target%charge(3) = -0.644520d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.147d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.896d0
    test_target%bond_order(2) = 0.895d0
    test_target%bond_order(3) = 0.896d0
    test_target%bond_order(4) = 0.895d0
    test_target%bond_order(5) = 0.895d0
    test_target%bond_order(6) = 0.895d0
    test_target%bond_order(7) = 1.791d0
    test_target%nerror = 0

    call mopac_scf_f(test_in, test_restore, test_out)
    call mopac_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_restart1

subroutine test_mozyme_restart1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mozyme_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O SCF MOZYME restart'

    ! SCF calculation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    ! API-based save/load is captured at slightly different point in SCF cycle,
    ! needs tigher tolerance for consistent results
    test_in%tolerance = 0.1d0
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.76d0
    test_in%coord(2) = 0.59d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.76d0
    test_in%coord(5) = 0.59d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.76970d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = 2.296d0
    test_target%coord_deriv(2) = 2.662d0
    test_target%coord_deriv(3) = 0.001d0
    test_target%coord_deriv(4) = -2.295d0
    test_target%coord_deriv(5) = 2.631d0
    test_target%coord_deriv(6) = 0.001d0
    test_target%coord_deriv(7) = -0.001d0
    test_target%coord_deriv(8) = -5.293d0
    test_target%coord_deriv(9) = -0.001d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.321963d0
    test_target%charge(2) = 0.321963d0
    test_target%charge(3) = -0.643926d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.146d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.896d0
    test_target%bond_order(2) = 0.895d0
    test_target%bond_order(3) = 0.896d0
    test_target%bond_order(4) = 0.895d0
    test_target%bond_order(5) = 0.895d0
    test_target%bond_order(6) = 0.895d0
    test_target%bond_order(7) = 1.791d0
    test_target%nerror = 0

    call mozyme_scf_f(test_in, test_restore, test_out)
    call mozyme_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_restart1

subroutine test_cosmo1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'H2O SCF COSMO'

    ! SCF calculation of H2O
    test_in%natom = 3
    test_in%natom_move = 3
    test_in%epsilon = 78.4d0
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.76d0
    test_in%coord(2) = 0.59d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.76d0
    test_in%coord(5) = 0.59d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 0.0d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -65.21458d0
    allocate(test_target%coord_update(3*3))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*3))
    test_target%coord_deriv(1) = -0.773696d0
    test_target%coord_deriv(2) = -1.302274d0
    test_target%coord_deriv(3) = 0.045812d0
    test_target%coord_deriv(4) = 0.734278d0
    test_target%coord_deriv(5) = -1.261587d0
    test_target%coord_deriv(6) = -0.029041d0
    test_target%coord_deriv(7) = 0.039417d0
    test_target%coord_deriv(8) = 2.563860d0
    test_target%coord_deriv(9) = -0.016772d0
    allocate(test_target%charge(3))
    test_target%charge(1) = 0.374234d0
    test_target%charge(2) = 0.374006d0
    test_target%charge(3) = -0.748240d0
    test_target%dipole(1) = 0.0d0
    test_target%dipole(2) = 2.43d0
    test_target%dipole(3) = 0.0d0
    test_target%stress(:) = 0.0d0
    allocate(test_target%bond_index(4))
    test_target%bond_index(1) = 1
    test_target%bond_index(2) = 3
    test_target%bond_index(3) = 5
    test_target%bond_index(4) = 8
    allocate(test_target%bond_atom(7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 3
    test_target%bond_atom(3) = 2
    test_target%bond_atom(4) = 3
    test_target%bond_atom(5) = 1
    test_target%bond_atom(6) = 2
    test_target%bond_atom(7) = 3
    allocate(test_target%bond_order(7))
    test_target%bond_order(1) = 0.860d0
    test_target%bond_order(2) = 0.859d0
    test_target%bond_order(3) = 0.860d0
    test_target%bond_order(4) = 0.859d0
    test_target%bond_order(5) = 0.859d0
    test_target%bond_order(6) = 0.859d0
    test_target%bond_order(7) = 1.718d0
    test_target%nerror = 0

    call mopac_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_cosmo1

subroutine test_crystal1(nfail)
    use mopac_api_f
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system_f) :: test_in
    type(mopac_state_f) :: test_restore
    type(mopac_properties_f) :: test_target
    type(mopac_properties_f) :: test_out
    character(50) :: test_name
    integer :: i
    test_name = 'LiF SCF'

    ! SCF calculation of LiF
    test_in%natom = 54
    test_in%natom_move = 54
    test_in%nlattice = 3
    test_in%nlattice_move = 3
    allocate(test_in%atom(54))
    test_in%atom(1) = 3
    test_in%atom(2) = 9
    test_in%atom(3) = 3
    test_in%atom(4) = 9
    test_in%atom(5) = 3
    test_in%atom(6) = 9
    test_in%atom(7) = 3
    test_in%atom(8) = 9
    test_in%atom(9) = 3
    test_in%atom(10) = 9
    test_in%atom(11) = 3
    test_in%atom(12) = 9
    test_in%atom(13) = 3
    test_in%atom(14) = 9
    test_in%atom(15) = 3
    test_in%atom(16) = 9
    test_in%atom(17) = 3
    test_in%atom(18) = 9
    test_in%atom(19) = 3
    test_in%atom(20) = 9
    test_in%atom(21) = 3
    test_in%atom(22) = 9
    test_in%atom(23) = 3
    test_in%atom(24) = 9
    test_in%atom(25) = 3
    test_in%atom(26) = 9
    test_in%atom(27) = 3
    test_in%atom(28) = 9
    test_in%atom(29) = 3
    test_in%atom(30) = 9
    test_in%atom(31) = 3
    test_in%atom(32) = 9
    test_in%atom(33) = 3
    test_in%atom(34) = 9
    test_in%atom(35) = 3
    test_in%atom(36) = 9
    test_in%atom(37) = 3
    test_in%atom(38) = 9
    test_in%atom(39) = 3
    test_in%atom(40) = 9
    test_in%atom(41) = 3
    test_in%atom(42) = 9
    test_in%atom(43) = 3
    test_in%atom(44) = 9
    test_in%atom(45) = 3
    test_in%atom(46) = 9
    test_in%atom(47) = 3
    test_in%atom(48) = 9
    test_in%atom(49) = 3
    test_in%atom(50) = 9
    test_in%atom(51) = 3
    test_in%atom(52) = 9
    test_in%atom(53) = 3
    test_in%atom(54) = 9
    allocate(test_in%coord(3*54))
    test_in%coord(1) = 0.0d0
    test_in%coord(2) = 0.0d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = 1.795d0
    test_in%coord(5) = 1.795d0
    test_in%coord(6) = 1.795d0
    test_in%coord(7) = 0.0d0
    test_in%coord(8) = 1.795d0
    test_in%coord(9) = 1.795d0
    test_in%coord(10) = 1.795d0
    test_in%coord(11) = 3.59d0
    test_in%coord(12) = 3.59d0
    test_in%coord(13) = 0.0d0
    test_in%coord(14) = 3.59d0
    test_in%coord(15) = 3.59d0
    test_in%coord(16) = 1.795d0
    test_in%coord(17) = 5.385d0
    test_in%coord(18) = 5.385d0
    test_in%coord(19) = 1.795d0
    test_in%coord(20) = 0.0d0
    test_in%coord(21) = 1.795d0
    test_in%coord(22) = 3.59d0
    test_in%coord(23) = 1.795d0
    test_in%coord(24) = 3.59d0
    test_in%coord(25) = 1.795d0
    test_in%coord(26) = 1.795d0
    test_in%coord(27) = 3.59d0
    test_in%coord(28) = 3.59d0
    test_in%coord(29) = 3.59d0
    test_in%coord(30) = 5.385d0
    test_in%coord(31) = 1.795d0
    test_in%coord(32) = 3.59d0
    test_in%coord(33) = 5.385d0
    test_in%coord(34) = 3.59d0
    test_in%coord(35) = 5.385d0
    test_in%coord(36) = 7.18d0
    test_in%coord(37) = 3.59d0
    test_in%coord(38) = 0.0d0
    test_in%coord(39) = 3.59d0
    test_in%coord(40) = 5.385d0
    test_in%coord(41) = 1.795d0
    test_in%coord(42) = 5.385d0
    test_in%coord(43) = 3.59d0
    test_in%coord(44) = 1.795d0
    test_in%coord(45) = 5.385d0
    test_in%coord(46) = 5.385d0
    test_in%coord(47) = 3.59d0
    test_in%coord(48) = 7.18d0
    test_in%coord(49) = 3.59d0
    test_in%coord(50) = 3.59d0
    test_in%coord(51) = 7.18d0
    test_in%coord(52) = 5.385d0
    test_in%coord(53) = 5.385d0
    test_in%coord(54) = 8.975d0
    test_in%coord(55) = 1.795d0
    test_in%coord(56) = 1.795d0
    test_in%coord(57) = 0.0d0
    test_in%coord(58) = 3.59d0
    test_in%coord(59) = 3.59d0
    test_in%coord(60) = 1.795d0
    test_in%coord(61) = 1.795d0
    test_in%coord(62) = 3.59d0
    test_in%coord(63) = 1.795d0
    test_in%coord(64) = 3.59d0
    test_in%coord(65) = 5.385d0
    test_in%coord(66) = 3.59d0
    test_in%coord(67) = 1.795d0
    test_in%coord(68) = 5.385d0
    test_in%coord(69) = 3.59d0
    test_in%coord(70) = 3.59d0
    test_in%coord(71) = 7.18d0
    test_in%coord(72) = 5.385d0
    test_in%coord(73) = 3.59d0
    test_in%coord(74) = 1.795d0
    test_in%coord(75) = 1.795d0
    test_in%coord(76) = 5.385d0
    test_in%coord(77) = 3.59d0
    test_in%coord(78) = 3.59d0
    test_in%coord(79) = 3.59d0
    test_in%coord(80) = 3.59d0
    test_in%coord(81) = 3.59d0
    test_in%coord(82) = 5.385d0
    test_in%coord(83) = 5.385d0
    test_in%coord(84) = 5.385d0
    test_in%coord(85) = 3.59d0
    test_in%coord(86) = 5.385d0
    test_in%coord(87) = 5.385d0
    test_in%coord(88) = 5.385d0
    test_in%coord(89) = 7.18d0
    test_in%coord(90) = 7.18d0
    test_in%coord(91) = 5.385d0
    test_in%coord(92) = 1.795d0
    test_in%coord(93) = 3.59d0
    test_in%coord(94) = 7.18d0
    test_in%coord(95) = 3.59d0
    test_in%coord(96) = 5.385d0
    test_in%coord(97) = 5.385d0
    test_in%coord(98) = 3.59d0
    test_in%coord(99) = 5.385d0
    test_in%coord(100) = 7.18d0
    test_in%coord(101) = 5.385d0
    test_in%coord(102) = 7.18d0
    test_in%coord(103) = 5.385d0
    test_in%coord(104) = 5.385d0
    test_in%coord(105) = 7.18d0
    test_in%coord(106) = 7.18d0
    test_in%coord(107) = 7.18d0
    test_in%coord(108) = 8.975d0
    test_in%coord(109) = 3.59d0
    test_in%coord(110) = 3.59d0
    test_in%coord(111) = 0.0d0
    test_in%coord(112) = 5.385d0
    test_in%coord(113) = 5.385d0
    test_in%coord(114) = 1.795d0
    test_in%coord(115) = 3.59d0
    test_in%coord(116) = 5.385d0
    test_in%coord(117) = 1.795d0
    test_in%coord(118) = 5.385d0
    test_in%coord(119) = 7.18d0
    test_in%coord(120) = 3.59d0
    test_in%coord(121) = 3.59d0
    test_in%coord(122) = 7.18d0
    test_in%coord(123) = 3.59d0
    test_in%coord(124) = 5.385d0
    test_in%coord(125) = 8.975d0
    test_in%coord(126) = 5.385d0
    test_in%coord(127) = 5.385d0
    test_in%coord(128) = 3.59d0
    test_in%coord(129) = 1.795d0
    test_in%coord(130) = 7.18d0
    test_in%coord(131) = 5.385d0
    test_in%coord(132) = 3.59d0
    test_in%coord(133) = 5.385d0
    test_in%coord(134) = 5.385d0
    test_in%coord(135) = 3.59d0
    test_in%coord(136) = 7.18d0
    test_in%coord(137) = 7.18d0
    test_in%coord(138) = 5.385d0
    test_in%coord(139) = 5.385d0
    test_in%coord(140) = 7.18d0
    test_in%coord(141) = 5.385d0
    test_in%coord(142) = 7.18d0
    test_in%coord(143) = 8.975d0
    test_in%coord(144) = 7.18d0
    test_in%coord(145) = 7.18d0
    test_in%coord(146) = 3.59d0
    test_in%coord(147) = 3.59d0
    test_in%coord(148) = 8.975d0
    test_in%coord(149) = 5.385d0
    test_in%coord(150) = 5.385d0
    test_in%coord(151) = 7.18d0
    test_in%coord(152) = 5.385d0
    test_in%coord(153) = 5.385d0
    test_in%coord(154) = 8.975d0
    test_in%coord(155) = 7.18d0
    test_in%coord(156) = 7.18d0
    test_in%coord(157) = 7.18d0
    test_in%coord(158) = 7.18d0
    test_in%coord(159) = 7.18d0
    test_in%coord(160) = 8.975d0
    test_in%coord(161) = 8.975d0
    test_in%coord(162) = 8.975d0
    allocate(test_in%lattice(3*3))
    test_in%lattice(1) = 5.385d0
    test_in%lattice(2) = 5.385d0
    test_in%lattice(3) = 0.0d0
    test_in%lattice(4) = 5.385d0
    test_in%lattice(5) = 0.0d0
    test_in%lattice(6) = 5.385d0
    test_in%lattice(7) = 0.0d0
    test_in%lattice(8) = 5.385d0
    test_in%lattice(9) = 5.385d0
    test_target%heat = -4588.63517d0
    allocate(test_target%coord_update(3*54))
    test_target%coord_update = test_in%coord
    allocate(test_target%coord_deriv(3*54))
    test_target%coord_deriv(1) = -0.003072d0
    test_target%coord_deriv(2) = -0.001544d0
    test_target%coord_deriv(3) = -0.000184d0
    test_target%coord_deriv(4) = -0.002616d0
    test_target%coord_deriv(5) = -0.000781d0
    test_target%coord_deriv(6) = -0.000527d0
    test_target%coord_deriv(7) = -0.003046d0
    test_target%coord_deriv(8) = -0.003246d0
    test_target%coord_deriv(9) = -0.001886d0
    test_target%coord_deriv(10) = -0.003656d0
    test_target%coord_deriv(11) = -0.002632d0
    test_target%coord_deriv(12) = -0.002377d0
    test_target%coord_deriv(13) = -0.004086d0
    test_target%coord_deriv(14) = -0.004968d0
    test_target%coord_deriv(15) = -0.003608d0
    test_target%coord_deriv(16) = -0.003656d0
    test_target%coord_deriv(17) = -0.004478d0
    test_target%coord_deriv(18) = -0.004224d0
    test_target%coord_deriv(19) = -0.003531d0
    test_target%coord_deriv(20) = 0.002452d0
    test_target%coord_deriv(21) = -0.001774d0
    test_target%coord_deriv(22) = -0.000814d0
    test_target%coord_deriv(23) = -0.000424d0
    test_target%coord_deriv(24) = -0.001538d0
    test_target%coord_deriv(25) = -0.003531d0
    test_target%coord_deriv(26) = 0.000725d0
    test_target%coord_deriv(27) = -0.003496d0
    test_target%coord_deriv(28) = -0.001854d0
    test_target%coord_deriv(29) = -0.002273d0
    test_target%coord_deriv(30) = -0.003385d0
    test_target%coord_deriv(31) = -0.004570d0
    test_target%coord_deriv(32) = -0.000998d0
    test_target%coord_deriv(33) = -0.005222d0
    test_target%coord_deriv(34) = -0.001854d0
    test_target%coord_deriv(35) = -0.004120d0
    test_target%coord_deriv(36) = -0.005232d0
    test_target%coord_deriv(37) = -0.000966d0
    test_target%coord_deriv(38) = 0.000375d0
    test_target%coord_deriv(39) = -0.001107d0
    test_target%coord_deriv(40) = -0.000700d0
    test_target%coord_deriv(41) = 0.003466d0
    test_target%coord_deriv(42) = -0.002285d0
    test_target%coord_deriv(43) = -0.000966d0
    test_target%coord_deriv(44) = -0.001352d0
    test_target%coord_deriv(45) = -0.002832d0
    test_target%coord_deriv(46) = -0.001741d0
    test_target%coord_deriv(47) = 0.001616d0
    test_target%coord_deriv(48) = -0.004131d0
    test_target%coord_deriv(49) = -0.002006d0
    test_target%coord_deriv(50) = -0.003074d0
    test_target%coord_deriv(51) = -0.004558d0
    test_target%coord_deriv(52) = -0.001740d0
    test_target%coord_deriv(53) = -0.000230d0
    test_target%coord_deriv(54) = -0.005979d0
    test_target%coord_deriv(55) = -0.003615d0
    test_target%coord_deriv(56) = -0.003741d0
    test_target%coord_deriv(57) = 0.009705d0
    test_target%coord_deriv(58) = 0.000734d0
    test_target%coord_deriv(59) = 0.001680d0
    test_target%coord_deriv(60) = -0.002221d0
    test_target%coord_deriv(61) = -0.003615d0
    test_target%coord_deriv(62) = -0.005463d0
    test_target%coord_deriv(63) = 0.007978d0
    test_target%coord_deriv(64) = -0.000306d0
    test_target%coord_deriv(65) = -0.000167d0
    test_target%coord_deriv(66) = -0.004070d0
    test_target%coord_deriv(67) = -0.004655d0
    test_target%coord_deriv(68) = -0.007189d0
    test_target%coord_deriv(69) = 0.006255d0
    test_target%coord_deriv(70) = -0.000306d0
    test_target%coord_deriv(71) = -0.002013d0
    test_target%coord_deriv(72) = -0.005917d0
    test_target%coord_deriv(73) = -0.004098d0
    test_target%coord_deriv(74) = 0.000230d0
    test_target%coord_deriv(75) = 0.008089d0
    test_target%coord_deriv(76) = 0.002535d0
    test_target%coord_deriv(77) = 0.002039d0
    test_target%coord_deriv(78) = -0.003233d0
    test_target%coord_deriv(79) = -0.004098d0
    test_target%coord_deriv(80) = -0.001493d0
    test_target%coord_deriv(81) = 0.006365d0
    test_target%coord_deriv(82) = 0.001495d0
    test_target%coord_deriv(83) = 0.000193d0
    test_target%coord_deriv(84) = -0.005079d0
    test_target%coord_deriv(85) = -0.005138d0
    test_target%coord_deriv(86) = -0.003218d0
    test_target%coord_deriv(87) = 0.004640d0
    test_target%coord_deriv(88) = 0.001495d0
    test_target%coord_deriv(89) = -0.001654d0
    test_target%coord_deriv(90) = -0.006926d0
    test_target%coord_deriv(91) = -0.001534d0
    test_target%coord_deriv(92) = -0.001845d0
    test_target%coord_deriv(93) = 0.008754d0
    test_target%coord_deriv(94) = 0.002650d0
    test_target%coord_deriv(95) = 0.005927d0
    test_target%coord_deriv(96) = -0.003979d0
    test_target%coord_deriv(97) = -0.001534d0
    test_target%coord_deriv(98) = -0.003569d0
    test_target%coord_deriv(99) = 0.007028d0
    test_target%coord_deriv(100) = 0.001610d0
    test_target%coord_deriv(101) = 0.004081d0
    test_target%coord_deriv(102) = -0.005826d0
    test_target%coord_deriv(103) = -0.002574d0
    test_target%coord_deriv(104) = -0.005294d0
    test_target%coord_deriv(105) = 0.005303d0
    test_target%coord_deriv(106) = 0.001610d0
    test_target%coord_deriv(107) = 0.002233d0
    test_target%coord_deriv(108) = -0.007674d0
    test_target%coord_deriv(109) = 0.003433d0
    test_target%coord_deriv(110) = 0.003270d0
    test_target%coord_deriv(111) = 0.002281d0
    test_target%coord_deriv(112) = 0.002916d0
    test_target%coord_deriv(113) = 0.001338d0
    test_target%coord_deriv(114) = 0.006969d0
    test_target%coord_deriv(115) = 0.003433d0
    test_target%coord_deriv(116) = 0.001544d0
    test_target%coord_deriv(117) = 0.000555d0
    test_target%coord_deriv(118) = 0.001876d0
    test_target%coord_deriv(119) = -0.000509d0
    test_target%coord_deriv(120) = 0.005122d0
    test_target%coord_deriv(121) = 0.002393d0
    test_target%coord_deriv(122) = -0.000181d0
    test_target%coord_deriv(123) = -0.001169d0
    test_target%coord_deriv(124) = 0.001876d0
    test_target%coord_deriv(125) = -0.002357d0
    test_target%coord_deriv(126) = 0.003275d0
    test_target%coord_deriv(127) = 0.002947d0
    test_target%coord_deriv(128) = 0.007238d0
    test_target%coord_deriv(129) = 0.000669d0
    test_target%coord_deriv(130) = 0.004719d0
    test_target%coord_deriv(131) = 0.001697d0
    test_target%coord_deriv(132) = 0.005959d0
    test_target%coord_deriv(133) = 0.002947d0
    test_target%coord_deriv(134) = 0.005512d0
    test_target%coord_deriv(135) = -0.001055d0
    test_target%coord_deriv(136) = 0.003679d0
    test_target%coord_deriv(137) = -0.000149d0
    test_target%coord_deriv(138) = 0.004112d0
    test_target%coord_deriv(139) = 0.001907d0
    test_target%coord_deriv(140) = 0.003787d0
    test_target%coord_deriv(141) = -0.002780d0
    test_target%coord_deriv(142) = 0.003679d0
    test_target%coord_deriv(143) = -0.001997d0
    test_target%coord_deriv(144) = 0.002265d0
    test_target%coord_deriv(145) = 0.005513d0
    test_target%coord_deriv(146) = 0.005162d0
    test_target%coord_deriv(147) = 0.001333d0
    test_target%coord_deriv(148) = 0.004834d0
    test_target%coord_deriv(149) = 0.005587d0
    test_target%coord_deriv(150) = 0.005213d0
    test_target%coord_deriv(151) = 0.005513d0
    test_target%coord_deriv(152) = 0.003437d0
    test_target%coord_deriv(153) = -0.000392d0
    test_target%coord_deriv(154) = 0.003794d0
    test_target%coord_deriv(155) = 0.003740d0
    test_target%coord_deriv(156) = 0.003366d0
    test_target%coord_deriv(157) = 0.004473d0
    test_target%coord_deriv(158) = 0.001712d0
    test_target%coord_deriv(159) = -0.002117d0
    test_target%coord_deriv(160) = 0.003794d0
    test_target%coord_deriv(161) = 0.001892d0
    test_target%coord_deriv(162) = 0.001519d0
    allocate(test_target%lattice_update(3*3))
    test_target%lattice_update = test_in%lattice
    allocate(test_target%lattice_deriv(3*3))
    test_target%lattice_deriv(1) = -2.329702d0
    test_target%lattice_deriv(2) = -2.324697d0
    test_target%lattice_deriv(3) = 2.296600d0
    test_target%lattice_deriv(4) = -2.314736d0
    test_target%lattice_deriv(5) = 2.301304d0
    test_target%lattice_deriv(6) = -2.299118d0
    test_target%lattice_deriv(7) = 2.302052d0
    test_target%lattice_deriv(8) = -2.296322d0
    test_target%lattice_deriv(9) = -2.281867d0
    allocate(test_target%charge(54))
    test_target%charge(1) = 0.855660d0
    test_target%charge(2) = -0.855642d0
    test_target%charge(3) = 0.855660d0
    test_target%charge(4) = -0.855642d0
    test_target%charge(5) = 0.855660d0
    test_target%charge(6) = -0.855642d0
    test_target%charge(7) = 0.855654d0
    test_target%charge(8) = -0.855644d0
    test_target%charge(9) = 0.855654d0
    test_target%charge(10) = -0.855644d0
    test_target%charge(11) = 0.855654d0
    test_target%charge(12) = -0.855644d0
    test_target%charge(13) = 0.855654d0
    test_target%charge(14) = -0.855642d0
    test_target%charge(15) = 0.855654d0
    test_target%charge(16) = -0.855642d0
    test_target%charge(17) = 0.855654d0
    test_target%charge(18) = -0.855642d0
    test_target%charge(19) = 0.855642d0
    test_target%charge(20) = -0.855647d0
    test_target%charge(21) = 0.855642d0
    test_target%charge(22) = -0.855647d0
    test_target%charge(23) = 0.855642d0
    test_target%charge(24) = -0.855647d0
    test_target%charge(25) = 0.855637d0
    test_target%charge(26) = -0.855649d0
    test_target%charge(27) = 0.855637d0
    test_target%charge(28) = -0.855649d0
    test_target%charge(29) = 0.855637d0
    test_target%charge(30) = -0.855649d0
    test_target%charge(31) = 0.855636d0
    test_target%charge(32) = -0.855647d0
    test_target%charge(33) = 0.855636d0
    test_target%charge(34) = -0.855647d0
    test_target%charge(35) = 0.855636d0
    test_target%charge(36) = -0.855647d0
    test_target%charge(37) = 0.855642d0
    test_target%charge(38) = -0.855641d0
    test_target%charge(39) = 0.855642d0
    test_target%charge(40) = -0.855641d0
    test_target%charge(41) = 0.855642d0
    test_target%charge(42) = -0.855641d0
    test_target%charge(43) = 0.855636d0
    test_target%charge(44) = -0.855643d0
    test_target%charge(45) = 0.855636d0
    test_target%charge(46) = -0.855643d0
    test_target%charge(47) = 0.855636d0
    test_target%charge(48) = -0.855643d0
    test_target%charge(49) = 0.855636d0
    test_target%charge(50) = -0.855641d0
    test_target%charge(51) = 0.855636d0
    test_target%charge(52) = -0.855641d0
    test_target%charge(53) = 0.855636d0
    test_target%charge(54) = -0.855641d0
    test_target%dipole(:) = 0.0d0
    test_target%stress(1) = -0.550d0
    test_target%stress(2) = -0.553d0
    test_target%stress(3) = -0.554d0
    test_target%stress(4) = 0.001d0
    test_target%stress(5) = 0.000d0
    test_target%stress(6) = 0.002d0
    allocate(test_target%bond_index(55))
    do i = 1, 55
        test_target%bond_index(i) = 1 + 7*(i-1)
    end do
    allocate(test_target%bond_atom(54*7))
    test_target%bond_atom(1) = 1
    test_target%bond_atom(2) = 6
    test_target%bond_atom(3) = 14
    test_target%bond_atom(4) = 18
    test_target%bond_atom(5) = 38
    test_target%bond_atom(6) = 42
    test_target%bond_atom(7) = 50
    test_target%bond_atom(8) = 2
    test_target%bond_atom(9) = 3
    test_target%bond_atom(10) = 7
    test_target%bond_atom(11) = 9
    test_target%bond_atom(12) = 19
    test_target%bond_atom(13) = 21
    test_target%bond_atom(14) = 25
    test_target%bond_atom(15) = 2
    test_target%bond_atom(16) = 3
    test_target%bond_atom(17) = 14
    test_target%bond_atom(18) = 16
    test_target%bond_atom(19) = 38
    test_target%bond_atom(20) = 40
    test_target%bond_atom(21) = 52
    test_target%bond_atom(22) = 4
    test_target%bond_atom(23) = 5
    test_target%bond_atom(24) = 9
    test_target%bond_atom(25) = 11
    test_target%bond_atom(26) = 21
    test_target%bond_atom(27) = 23
    test_target%bond_atom(28) = 27
    test_target%bond_atom(29) = 4
    test_target%bond_atom(30) = 5
    test_target%bond_atom(31) = 16
    test_target%bond_atom(32) = 18
    test_target%bond_atom(33) = 40
    test_target%bond_atom(34) = 42
    test_target%bond_atom(35) = 54
    test_target%bond_atom(36) = 1
    test_target%bond_atom(37) = 6
    test_target%bond_atom(38) = 7
    test_target%bond_atom(39) = 11
    test_target%bond_atom(40) = 19
    test_target%bond_atom(41) = 23
    test_target%bond_atom(42) = 29
    test_target%bond_atom(43) = 2
    test_target%bond_atom(44) = 6
    test_target%bond_atom(45) = 7
    test_target%bond_atom(46) = 12
    test_target%bond_atom(47) = 38
    test_target%bond_atom(48) = 44
    test_target%bond_atom(49) = 48
    test_target%bond_atom(50) = 8
    test_target%bond_atom(51) = 9
    test_target%bond_atom(52) = 13
    test_target%bond_atom(53) = 15
    test_target%bond_atom(54) = 25
    test_target%bond_atom(55) = 27
    test_target%bond_atom(56) = 31
    test_target%bond_atom(57) = 2
    test_target%bond_atom(58) = 4
    test_target%bond_atom(59) = 8
    test_target%bond_atom(60) = 9
    test_target%bond_atom(61) = 40
    test_target%bond_atom(62) = 44
    test_target%bond_atom(63) = 46
    test_target%bond_atom(64) = 10
    test_target%bond_atom(65) = 11
    test_target%bond_atom(66) = 15
    test_target%bond_atom(67) = 17
    test_target%bond_atom(68) = 27
    test_target%bond_atom(69) = 29
    test_target%bond_atom(70) = 33
    test_target%bond_atom(71) = 4
    test_target%bond_atom(72) = 6
    test_target%bond_atom(73) = 10
    test_target%bond_atom(74) = 11
    test_target%bond_atom(75) = 42
    test_target%bond_atom(76) = 46
    test_target%bond_atom(77) = 48
    test_target%bond_atom(78) = 7
    test_target%bond_atom(79) = 12
    test_target%bond_atom(80) = 13
    test_target%bond_atom(81) = 17
    test_target%bond_atom(82) = 25
    test_target%bond_atom(83) = 29
    test_target%bond_atom(84) = 35
    test_target%bond_atom(85) = 8
    test_target%bond_atom(86) = 12
    test_target%bond_atom(87) = 13
    test_target%bond_atom(88) = 18
    test_target%bond_atom(89) = 44
    test_target%bond_atom(90) = 50
    test_target%bond_atom(91) = 54
    test_target%bond_atom(92) = 1
    test_target%bond_atom(93) = 3
    test_target%bond_atom(94) = 14
    test_target%bond_atom(95) = 15
    test_target%bond_atom(96) = 19
    test_target%bond_atom(97) = 31
    test_target%bond_atom(98) = 33
    test_target%bond_atom(99) = 8
    test_target%bond_atom(100) = 10
    test_target%bond_atom(101) = 14
    test_target%bond_atom(102) = 15
    test_target%bond_atom(103) = 46
    test_target%bond_atom(104) = 50
    test_target%bond_atom(105) = 52
    test_target%bond_atom(106) = 3
    test_target%bond_atom(107) = 5
    test_target%bond_atom(108) = 16
    test_target%bond_atom(109) = 17
    test_target%bond_atom(110) = 21
    test_target%bond_atom(111) = 33
    test_target%bond_atom(112) = 35
    test_target%bond_atom(113) = 10
    test_target%bond_atom(114) = 12
    test_target%bond_atom(115) = 16
    test_target%bond_atom(116) = 17
    test_target%bond_atom(117) = 48
    test_target%bond_atom(118) = 52
    test_target%bond_atom(119) = 54
    test_target%bond_atom(120) = 1
    test_target%bond_atom(121) = 5
    test_target%bond_atom(122) = 13
    test_target%bond_atom(123) = 18
    test_target%bond_atom(124) = 23
    test_target%bond_atom(125) = 31
    test_target%bond_atom(126) = 35
    test_target%bond_atom(127) = 2
    test_target%bond_atom(128) = 6
    test_target%bond_atom(129) = 14
    test_target%bond_atom(130) = 19
    test_target%bond_atom(131) = 24
    test_target%bond_atom(132) = 32
    test_target%bond_atom(133) = 36
    test_target%bond_atom(134) = 20
    test_target%bond_atom(135) = 21
    test_target%bond_atom(136) = 25
    test_target%bond_atom(137) = 27
    test_target%bond_atom(138) = 37
    test_target%bond_atom(139) = 39
    test_target%bond_atom(140) = 43
    test_target%bond_atom(141) = 2
    test_target%bond_atom(142) = 4
    test_target%bond_atom(143) = 16
    test_target%bond_atom(144) = 20
    test_target%bond_atom(145) = 21
    test_target%bond_atom(146) = 32
    test_target%bond_atom(147) = 34
    test_target%bond_atom(148) = 22
    test_target%bond_atom(149) = 23
    test_target%bond_atom(150) = 27
    test_target%bond_atom(151) = 29
    test_target%bond_atom(152) = 39
    test_target%bond_atom(153) = 41
    test_target%bond_atom(154) = 45
    test_target%bond_atom(155) = 4
    test_target%bond_atom(156) = 6
    test_target%bond_atom(157) = 18
    test_target%bond_atom(158) = 22
    test_target%bond_atom(159) = 23
    test_target%bond_atom(160) = 34
    test_target%bond_atom(161) = 36
    test_target%bond_atom(162) = 19
    test_target%bond_atom(163) = 24
    test_target%bond_atom(164) = 25
    test_target%bond_atom(165) = 29
    test_target%bond_atom(166) = 37
    test_target%bond_atom(167) = 41
    test_target%bond_atom(168) = 47
    test_target%bond_atom(169) = 2
    test_target%bond_atom(170) = 8
    test_target%bond_atom(171) = 12
    test_target%bond_atom(172) = 20
    test_target%bond_atom(173) = 24
    test_target%bond_atom(174) = 25
    test_target%bond_atom(175) = 30
    test_target%bond_atom(176) = 26
    test_target%bond_atom(177) = 27
    test_target%bond_atom(178) = 31
    test_target%bond_atom(179) = 33
    test_target%bond_atom(180) = 43
    test_target%bond_atom(181) = 45
    test_target%bond_atom(182) = 49
    test_target%bond_atom(183) = 4
    test_target%bond_atom(184) = 8
    test_target%bond_atom(185) = 10
    test_target%bond_atom(186) = 20
    test_target%bond_atom(187) = 22
    test_target%bond_atom(188) = 26
    test_target%bond_atom(189) = 27
    test_target%bond_atom(190) = 28
    test_target%bond_atom(191) = 29
    test_target%bond_atom(192) = 33
    test_target%bond_atom(193) = 35
    test_target%bond_atom(194) = 45
    test_target%bond_atom(195) = 47
    test_target%bond_atom(196) = 51
    test_target%bond_atom(197) = 6
    test_target%bond_atom(198) = 10
    test_target%bond_atom(199) = 12
    test_target%bond_atom(200) = 22
    test_target%bond_atom(201) = 24
    test_target%bond_atom(202) = 28
    test_target%bond_atom(203) = 29
    test_target%bond_atom(204) = 25
    test_target%bond_atom(205) = 30
    test_target%bond_atom(206) = 31
    test_target%bond_atom(207) = 35
    test_target%bond_atom(208) = 43
    test_target%bond_atom(209) = 47
    test_target%bond_atom(210) = 53
    test_target%bond_atom(211) = 8
    test_target%bond_atom(212) = 14
    test_target%bond_atom(213) = 18
    test_target%bond_atom(214) = 26
    test_target%bond_atom(215) = 30
    test_target%bond_atom(216) = 31
    test_target%bond_atom(217) = 36
    test_target%bond_atom(218) = 19
    test_target%bond_atom(219) = 21
    test_target%bond_atom(220) = 32
    test_target%bond_atom(221) = 33
    test_target%bond_atom(222) = 37
    test_target%bond_atom(223) = 49
    test_target%bond_atom(224) = 51
    test_target%bond_atom(225) = 10
    test_target%bond_atom(226) = 14
    test_target%bond_atom(227) = 16
    test_target%bond_atom(228) = 26
    test_target%bond_atom(229) = 28
    test_target%bond_atom(230) = 32
    test_target%bond_atom(231) = 33
    test_target%bond_atom(232) = 21
    test_target%bond_atom(233) = 23
    test_target%bond_atom(234) = 34
    test_target%bond_atom(235) = 35
    test_target%bond_atom(236) = 39
    test_target%bond_atom(237) = 51
    test_target%bond_atom(238) = 53
    test_target%bond_atom(239) = 12
    test_target%bond_atom(240) = 16
    test_target%bond_atom(241) = 18
    test_target%bond_atom(242) = 28
    test_target%bond_atom(243) = 30
    test_target%bond_atom(244) = 34
    test_target%bond_atom(245) = 35
    test_target%bond_atom(246) = 19
    test_target%bond_atom(247) = 23
    test_target%bond_atom(248) = 31
    test_target%bond_atom(249) = 36
    test_target%bond_atom(250) = 41
    test_target%bond_atom(251) = 49
    test_target%bond_atom(252) = 53
    test_target%bond_atom(253) = 20
    test_target%bond_atom(254) = 24
    test_target%bond_atom(255) = 32
    test_target%bond_atom(256) = 37
    test_target%bond_atom(257) = 42
    test_target%bond_atom(258) = 50
    test_target%bond_atom(259) = 54
    test_target%bond_atom(260) = 1
    test_target%bond_atom(261) = 3
    test_target%bond_atom(262) = 7
    test_target%bond_atom(263) = 38
    test_target%bond_atom(264) = 39
    test_target%bond_atom(265) = 43
    test_target%bond_atom(266) = 45
    test_target%bond_atom(267) = 20
    test_target%bond_atom(268) = 22
    test_target%bond_atom(269) = 34
    test_target%bond_atom(270) = 38
    test_target%bond_atom(271) = 39
    test_target%bond_atom(272) = 50
    test_target%bond_atom(273) = 52
    test_target%bond_atom(274) = 3
    test_target%bond_atom(275) = 5
    test_target%bond_atom(276) = 9
    test_target%bond_atom(277) = 40
    test_target%bond_atom(278) = 41
    test_target%bond_atom(279) = 45
    test_target%bond_atom(280) = 47
    test_target%bond_atom(281) = 22
    test_target%bond_atom(282) = 24
    test_target%bond_atom(283) = 36
    test_target%bond_atom(284) = 40
    test_target%bond_atom(285) = 41
    test_target%bond_atom(286) = 52
    test_target%bond_atom(287) = 54
    test_target%bond_atom(288) = 1
    test_target%bond_atom(289) = 5
    test_target%bond_atom(290) = 11
    test_target%bond_atom(291) = 37
    test_target%bond_atom(292) = 42
    test_target%bond_atom(293) = 43
    test_target%bond_atom(294) = 47
    test_target%bond_atom(295) = 20
    test_target%bond_atom(296) = 26
    test_target%bond_atom(297) = 30
    test_target%bond_atom(298) = 38
    test_target%bond_atom(299) = 42
    test_target%bond_atom(300) = 43
    test_target%bond_atom(301) = 48
    test_target%bond_atom(302) = 7
    test_target%bond_atom(303) = 9
    test_target%bond_atom(304) = 13
    test_target%bond_atom(305) = 44
    test_target%bond_atom(306) = 45
    test_target%bond_atom(307) = 49
    test_target%bond_atom(308) = 51
    test_target%bond_atom(309) = 22
    test_target%bond_atom(310) = 26
    test_target%bond_atom(311) = 28
    test_target%bond_atom(312) = 38
    test_target%bond_atom(313) = 40
    test_target%bond_atom(314) = 44
    test_target%bond_atom(315) = 45
    test_target%bond_atom(316) = 9
    test_target%bond_atom(317) = 11
    test_target%bond_atom(318) = 15
    test_target%bond_atom(319) = 46
    test_target%bond_atom(320) = 47
    test_target%bond_atom(321) = 51
    test_target%bond_atom(322) = 53
    test_target%bond_atom(323) = 24
    test_target%bond_atom(324) = 28
    test_target%bond_atom(325) = 30
    test_target%bond_atom(326) = 40
    test_target%bond_atom(327) = 42
    test_target%bond_atom(328) = 46
    test_target%bond_atom(329) = 47
    test_target%bond_atom(330) = 7
    test_target%bond_atom(331) = 11
    test_target%bond_atom(332) = 17
    test_target%bond_atom(333) = 43
    test_target%bond_atom(334) = 48
    test_target%bond_atom(335) = 49
    test_target%bond_atom(336) = 53
    test_target%bond_atom(337) = 26
    test_target%bond_atom(338) = 32
    test_target%bond_atom(339) = 36
    test_target%bond_atom(340) = 44
    test_target%bond_atom(341) = 48
    test_target%bond_atom(342) = 49
    test_target%bond_atom(343) = 54
    test_target%bond_atom(344) = 1
    test_target%bond_atom(345) = 13
    test_target%bond_atom(346) = 15
    test_target%bond_atom(347) = 37
    test_target%bond_atom(348) = 39
    test_target%bond_atom(349) = 50
    test_target%bond_atom(350) = 51
    test_target%bond_atom(351) = 28
    test_target%bond_atom(352) = 32
    test_target%bond_atom(353) = 34
    test_target%bond_atom(354) = 44
    test_target%bond_atom(355) = 46
    test_target%bond_atom(356) = 50
    test_target%bond_atom(357) = 51
    test_target%bond_atom(358) = 3
    test_target%bond_atom(359) = 15
    test_target%bond_atom(360) = 17
    test_target%bond_atom(361) = 39
    test_target%bond_atom(362) = 41
    test_target%bond_atom(363) = 52
    test_target%bond_atom(364) = 53
    test_target%bond_atom(365) = 30
    test_target%bond_atom(366) = 34
    test_target%bond_atom(367) = 36
    test_target%bond_atom(368) = 46
    test_target%bond_atom(369) = 48
    test_target%bond_atom(370) = 52
    test_target%bond_atom(371) = 53
    test_target%bond_atom(372) = 5
    test_target%bond_atom(373) = 13
    test_target%bond_atom(374) = 17
    test_target%bond_atom(375) = 37
    test_target%bond_atom(376) = 41
    test_target%bond_atom(377) = 49
    test_target%bond_atom(378) = 54
    allocate(test_target%bond_order(54*7))
    test_target%bond_order(1) = 0.268d0
    test_target%bond_order(2) = 0.037d0
    test_target%bond_order(3) = 0.037d0
    test_target%bond_order(4) = 0.037d0
    test_target%bond_order(5) = 0.037d0
    test_target%bond_order(6) = 0.037d0
    test_target%bond_order(7) = 0.037d0
    test_target%bond_order(8) = 0.283d0
    test_target%bond_order(9) = 0.037d0
    test_target%bond_order(10) = 0.037d0
    test_target%bond_order(11) = 0.037d0
    test_target%bond_order(12) = 0.037d0
    test_target%bond_order(13) = 0.037d0
    test_target%bond_order(14) = 0.037d0
    test_target%bond_order(15) = 0.037d0
    test_target%bond_order(16) = 0.268d0
    test_target%bond_order(17) = 0.037d0
    test_target%bond_order(18) = 0.037d0
    test_target%bond_order(19) = 0.037d0
    test_target%bond_order(20) = 0.037d0
    test_target%bond_order(21) = 0.037d0
    test_target%bond_order(22) = 0.283d0
    test_target%bond_order(23) = 0.037d0
    test_target%bond_order(24) = 0.037d0
    test_target%bond_order(25) = 0.037d0
    test_target%bond_order(26) = 0.037d0
    test_target%bond_order(27) = 0.037d0
    test_target%bond_order(28) = 0.037d0
    test_target%bond_order(29) = 0.037d0
    test_target%bond_order(30) = 0.268d0
    test_target%bond_order(31) = 0.037d0
    test_target%bond_order(32) = 0.037d0
    test_target%bond_order(33) = 0.037d0
    test_target%bond_order(34) = 0.037d0
    test_target%bond_order(35) = 0.037d0
    test_target%bond_order(36) = 0.037d0
    test_target%bond_order(37) = 0.283d0
    test_target%bond_order(38) = 0.037d0
    test_target%bond_order(39) = 0.037d0
    test_target%bond_order(40) = 0.037d0
    test_target%bond_order(41) = 0.037d0
    test_target%bond_order(42) = 0.037d0
    test_target%bond_order(43) = 0.037d0
    test_target%bond_order(44) = 0.037d0
    test_target%bond_order(45) = 0.268d0
    test_target%bond_order(46) = 0.037d0
    test_target%bond_order(47) = 0.037d0
    test_target%bond_order(48) = 0.037d0
    test_target%bond_order(49) = 0.037d0
    test_target%bond_order(50) = 0.283d0
    test_target%bond_order(51) = 0.037d0
    test_target%bond_order(52) = 0.037d0
    test_target%bond_order(53) = 0.037d0
    test_target%bond_order(54) = 0.037d0
    test_target%bond_order(55) = 0.037d0
    test_target%bond_order(56) = 0.037d0
    test_target%bond_order(57) = 0.037d0
    test_target%bond_order(58) = 0.037d0
    test_target%bond_order(59) = 0.037d0
    test_target%bond_order(60) = 0.268d0
    test_target%bond_order(61) = 0.037d0
    test_target%bond_order(62) = 0.037d0
    test_target%bond_order(63) = 0.037d0
    test_target%bond_order(64) = 0.283d0
    test_target%bond_order(65) = 0.037d0
    test_target%bond_order(66) = 0.037d0
    test_target%bond_order(67) = 0.037d0
    test_target%bond_order(68) = 0.037d0
    test_target%bond_order(69) = 0.037d0
    test_target%bond_order(70) = 0.037d0
    test_target%bond_order(71) = 0.037d0
    test_target%bond_order(72) = 0.037d0
    test_target%bond_order(73) = 0.037d0
    test_target%bond_order(74) = 0.268d0
    test_target%bond_order(75) = 0.037d0
    test_target%bond_order(76) = 0.037d0
    test_target%bond_order(77) = 0.037d0
    test_target%bond_order(78) = 0.037d0
    test_target%bond_order(79) = 0.283d0
    test_target%bond_order(80) = 0.037d0
    test_target%bond_order(81) = 0.037d0
    test_target%bond_order(82) = 0.037d0
    test_target%bond_order(83) = 0.037d0
    test_target%bond_order(84) = 0.037d0
    test_target%bond_order(85) = 0.037d0
    test_target%bond_order(86) = 0.037d0
    test_target%bond_order(87) = 0.268d0
    test_target%bond_order(88) = 0.037d0
    test_target%bond_order(89) = 0.037d0
    test_target%bond_order(90) = 0.037d0
    test_target%bond_order(91) = 0.037d0
    test_target%bond_order(92) = 0.037d0
    test_target%bond_order(93) = 0.037d0
    test_target%bond_order(94) = 0.283d0
    test_target%bond_order(95) = 0.037d0
    test_target%bond_order(96) = 0.037d0
    test_target%bond_order(97) = 0.037d0
    test_target%bond_order(98) = 0.037d0
    test_target%bond_order(99) = 0.037d0
    test_target%bond_order(100) = 0.037d0
    test_target%bond_order(101) = 0.037d0
    test_target%bond_order(102) = 0.268d0
    test_target%bond_order(103) = 0.037d0
    test_target%bond_order(104) = 0.037d0
    test_target%bond_order(105) = 0.037d0
    test_target%bond_order(106) = 0.037d0
    test_target%bond_order(107) = 0.037d0
    test_target%bond_order(108) = 0.283d0
    test_target%bond_order(109) = 0.037d0
    test_target%bond_order(110) = 0.037d0
    test_target%bond_order(111) = 0.037d0
    test_target%bond_order(112) = 0.037d0
    test_target%bond_order(113) = 0.037d0
    test_target%bond_order(114) = 0.037d0
    test_target%bond_order(115) = 0.037d0
    test_target%bond_order(116) = 0.268d0
    test_target%bond_order(117) = 0.037d0
    test_target%bond_order(118) = 0.037d0
    test_target%bond_order(119) = 0.037d0
    test_target%bond_order(120) = 0.037d0
    test_target%bond_order(121) = 0.037d0
    test_target%bond_order(122) = 0.037d0
    test_target%bond_order(123) = 0.283d0
    test_target%bond_order(124) = 0.037d0
    test_target%bond_order(125) = 0.037d0
    test_target%bond_order(126) = 0.037d0
    test_target%bond_order(127) = 0.037d0
    test_target%bond_order(128) = 0.037d0
    test_target%bond_order(129) = 0.037d0
    test_target%bond_order(130) = 0.268d0
    test_target%bond_order(131) = 0.037d0
    test_target%bond_order(132) = 0.037d0
    test_target%bond_order(133) = 0.037d0
    test_target%bond_order(134) = 0.283d0
    test_target%bond_order(135) = 0.037d0
    test_target%bond_order(136) = 0.037d0
    test_target%bond_order(137) = 0.037d0
    test_target%bond_order(138) = 0.037d0
    test_target%bond_order(139) = 0.037d0
    test_target%bond_order(140) = 0.037d0
    test_target%bond_order(141) = 0.037d0
    test_target%bond_order(142) = 0.037d0
    test_target%bond_order(143) = 0.037d0
    test_target%bond_order(144) = 0.037d0
    test_target%bond_order(145) = 0.268d0
    test_target%bond_order(146) = 0.037d0
    test_target%bond_order(147) = 0.037d0
    test_target%bond_order(148) = 0.283d0
    test_target%bond_order(149) = 0.037d0
    test_target%bond_order(150) = 0.037d0
    test_target%bond_order(151) = 0.037d0
    test_target%bond_order(152) = 0.037d0
    test_target%bond_order(153) = 0.037d0
    test_target%bond_order(154) = 0.037d0
    test_target%bond_order(155) = 0.037d0
    test_target%bond_order(156) = 0.037d0
    test_target%bond_order(157) = 0.037d0
    test_target%bond_order(158) = 0.037d0
    test_target%bond_order(159) = 0.268d0
    test_target%bond_order(160) = 0.037d0
    test_target%bond_order(161) = 0.037d0
    test_target%bond_order(162) = 0.037d0
    test_target%bond_order(163) = 0.283d0
    test_target%bond_order(164) = 0.037d0
    test_target%bond_order(165) = 0.037d0
    test_target%bond_order(166) = 0.037d0
    test_target%bond_order(167) = 0.037d0
    test_target%bond_order(168) = 0.037d0
    test_target%bond_order(169) = 0.037d0
    test_target%bond_order(170) = 0.037d0
    test_target%bond_order(171) = 0.037d0
    test_target%bond_order(172) = 0.037d0
    test_target%bond_order(173) = 0.037d0
    test_target%bond_order(174) = 0.268d0
    test_target%bond_order(175) = 0.037d0
    test_target%bond_order(176) = 0.283d0
    test_target%bond_order(177) = 0.037d0
    test_target%bond_order(178) = 0.037d0
    test_target%bond_order(179) = 0.037d0
    test_target%bond_order(180) = 0.037d0
    test_target%bond_order(181) = 0.037d0
    test_target%bond_order(182) = 0.037d0
    test_target%bond_order(183) = 0.037d0
    test_target%bond_order(184) = 0.037d0
    test_target%bond_order(185) = 0.037d0
    test_target%bond_order(186) = 0.037d0
    test_target%bond_order(187) = 0.037d0
    test_target%bond_order(188) = 0.037d0
    test_target%bond_order(189) = 0.268d0
    test_target%bond_order(190) = 0.283d0
    test_target%bond_order(191) = 0.037d0
    test_target%bond_order(192) = 0.037d0
    test_target%bond_order(193) = 0.037d0
    test_target%bond_order(194) = 0.037d0
    test_target%bond_order(195) = 0.037d0
    test_target%bond_order(196) = 0.037d0
    test_target%bond_order(197) = 0.037d0
    test_target%bond_order(198) = 0.037d0
    test_target%bond_order(199) = 0.037d0
    test_target%bond_order(200) = 0.037d0
    test_target%bond_order(201) = 0.037d0
    test_target%bond_order(202) = 0.037d0
    test_target%bond_order(203) = 0.268d0
    test_target%bond_order(204) = 0.037d0
    test_target%bond_order(205) = 0.283d0
    test_target%bond_order(206) = 0.037d0
    test_target%bond_order(207) = 0.037d0
    test_target%bond_order(208) = 0.037d0
    test_target%bond_order(209) = 0.037d0
    test_target%bond_order(210) = 0.037d0
    test_target%bond_order(211) = 0.037d0
    test_target%bond_order(212) = 0.037d0
    test_target%bond_order(213) = 0.037d0
    test_target%bond_order(214) = 0.037d0
    test_target%bond_order(215) = 0.037d0
    test_target%bond_order(216) = 0.268d0
    test_target%bond_order(217) = 0.037d0
    test_target%bond_order(218) = 0.037d0
    test_target%bond_order(219) = 0.037d0
    test_target%bond_order(220) = 0.283d0
    test_target%bond_order(221) = 0.037d0
    test_target%bond_order(222) = 0.037d0
    test_target%bond_order(223) = 0.037d0
    test_target%bond_order(224) = 0.037d0
    test_target%bond_order(225) = 0.037d0
    test_target%bond_order(226) = 0.037d0
    test_target%bond_order(227) = 0.037d0
    test_target%bond_order(228) = 0.037d0
    test_target%bond_order(229) = 0.037d0
    test_target%bond_order(230) = 0.037d0
    test_target%bond_order(231) = 0.268d0
    test_target%bond_order(232) = 0.037d0
    test_target%bond_order(233) = 0.037d0
    test_target%bond_order(234) = 0.283d0
    test_target%bond_order(235) = 0.037d0
    test_target%bond_order(236) = 0.037d0
    test_target%bond_order(237) = 0.037d0
    test_target%bond_order(238) = 0.037d0
    test_target%bond_order(239) = 0.037d0
    test_target%bond_order(240) = 0.037d0
    test_target%bond_order(241) = 0.037d0
    test_target%bond_order(242) = 0.037d0
    test_target%bond_order(243) = 0.037d0
    test_target%bond_order(244) = 0.037d0
    test_target%bond_order(245) = 0.268d0
    test_target%bond_order(246) = 0.037d0
    test_target%bond_order(247) = 0.037d0
    test_target%bond_order(248) = 0.037d0
    test_target%bond_order(249) = 0.283d0
    test_target%bond_order(250) = 0.037d0
    test_target%bond_order(251) = 0.037d0
    test_target%bond_order(252) = 0.037d0
    test_target%bond_order(253) = 0.037d0
    test_target%bond_order(254) = 0.037d0
    test_target%bond_order(255) = 0.037d0
    test_target%bond_order(256) = 0.268d0
    test_target%bond_order(257) = 0.037d0
    test_target%bond_order(258) = 0.037d0
    test_target%bond_order(259) = 0.037d0
    test_target%bond_order(260) = 0.037d0
    test_target%bond_order(261) = 0.037d0
    test_target%bond_order(262) = 0.037d0
    test_target%bond_order(263) = 0.283d0
    test_target%bond_order(264) = 0.037d0
    test_target%bond_order(265) = 0.037d0
    test_target%bond_order(266) = 0.037d0
    test_target%bond_order(267) = 0.037d0
    test_target%bond_order(268) = 0.037d0
    test_target%bond_order(269) = 0.037d0
    test_target%bond_order(270) = 0.037d0
    test_target%bond_order(271) = 0.268d0
    test_target%bond_order(272) = 0.037d0
    test_target%bond_order(273) = 0.037d0
    test_target%bond_order(274) = 0.037d0
    test_target%bond_order(275) = 0.037d0
    test_target%bond_order(276) = 0.037d0
    test_target%bond_order(277) = 0.283d0
    test_target%bond_order(278) = 0.037d0
    test_target%bond_order(279) = 0.037d0
    test_target%bond_order(280) = 0.037d0
    test_target%bond_order(281) = 0.037d0
    test_target%bond_order(282) = 0.037d0
    test_target%bond_order(283) = 0.037d0
    test_target%bond_order(284) = 0.037d0
    test_target%bond_order(285) = 0.268d0
    test_target%bond_order(286) = 0.037d0
    test_target%bond_order(287) = 0.037d0
    test_target%bond_order(288) = 0.037d0
    test_target%bond_order(289) = 0.037d0
    test_target%bond_order(290) = 0.037d0
    test_target%bond_order(291) = 0.037d0
    test_target%bond_order(292) = 0.283d0
    test_target%bond_order(293) = 0.037d0
    test_target%bond_order(294) = 0.037d0
    test_target%bond_order(295) = 0.037d0
    test_target%bond_order(296) = 0.037d0
    test_target%bond_order(297) = 0.037d0
    test_target%bond_order(298) = 0.037d0
    test_target%bond_order(299) = 0.037d0
    test_target%bond_order(300) = 0.268d0
    test_target%bond_order(301) = 0.037d0
    test_target%bond_order(302) = 0.037d0
    test_target%bond_order(303) = 0.037d0
    test_target%bond_order(304) = 0.037d0
    test_target%bond_order(305) = 0.283d0
    test_target%bond_order(306) = 0.037d0
    test_target%bond_order(307) = 0.037d0
    test_target%bond_order(308) = 0.037d0
    test_target%bond_order(309) = 0.037d0
    test_target%bond_order(310) = 0.037d0
    test_target%bond_order(311) = 0.037d0
    test_target%bond_order(312) = 0.037d0
    test_target%bond_order(313) = 0.037d0
    test_target%bond_order(314) = 0.037d0
    test_target%bond_order(315) = 0.268d0
    test_target%bond_order(316) = 0.037d0
    test_target%bond_order(317) = 0.037d0
    test_target%bond_order(318) = 0.037d0
    test_target%bond_order(319) = 0.283d0
    test_target%bond_order(320) = 0.037d0
    test_target%bond_order(321) = 0.037d0
    test_target%bond_order(322) = 0.037d0
    test_target%bond_order(323) = 0.037d0
    test_target%bond_order(324) = 0.037d0
    test_target%bond_order(325) = 0.037d0
    test_target%bond_order(326) = 0.037d0
    test_target%bond_order(327) = 0.037d0
    test_target%bond_order(328) = 0.037d0
    test_target%bond_order(329) = 0.268d0
    test_target%bond_order(330) = 0.037d0
    test_target%bond_order(331) = 0.037d0
    test_target%bond_order(332) = 0.037d0
    test_target%bond_order(333) = 0.037d0
    test_target%bond_order(334) = 0.283d0
    test_target%bond_order(335) = 0.037d0
    test_target%bond_order(336) = 0.037d0
    test_target%bond_order(337) = 0.037d0
    test_target%bond_order(338) = 0.037d0
    test_target%bond_order(339) = 0.037d0
    test_target%bond_order(340) = 0.037d0
    test_target%bond_order(341) = 0.037d0
    test_target%bond_order(342) = 0.268d0
    test_target%bond_order(343) = 0.037d0
    test_target%bond_order(344) = 0.037d0
    test_target%bond_order(345) = 0.037d0
    test_target%bond_order(346) = 0.037d0
    test_target%bond_order(347) = 0.037d0
    test_target%bond_order(348) = 0.037d0
    test_target%bond_order(349) = 0.283d0
    test_target%bond_order(350) = 0.037d0
    test_target%bond_order(351) = 0.037d0
    test_target%bond_order(352) = 0.037d0
    test_target%bond_order(353) = 0.037d0
    test_target%bond_order(354) = 0.037d0
    test_target%bond_order(355) = 0.037d0
    test_target%bond_order(356) = 0.037d0
    test_target%bond_order(357) = 0.268d0
    test_target%bond_order(358) = 0.037d0
    test_target%bond_order(359) = 0.037d0
    test_target%bond_order(360) = 0.037d0
    test_target%bond_order(361) = 0.037d0
    test_target%bond_order(362) = 0.037d0
    test_target%bond_order(363) = 0.283d0
    test_target%bond_order(364) = 0.037d0
    test_target%bond_order(365) = 0.037d0
    test_target%bond_order(366) = 0.037d0
    test_target%bond_order(367) = 0.037d0
    test_target%bond_order(368) = 0.037d0
    test_target%bond_order(369) = 0.037d0
    test_target%bond_order(370) = 0.037d0
    test_target%bond_order(371) = 0.268d0
    test_target%bond_order(372) = 0.037d0
    test_target%bond_order(373) = 0.037d0
    test_target%bond_order(374) = 0.037d0
    test_target%bond_order(375) = 0.037d0
    test_target%bond_order(376) = 0.037d0
    test_target%bond_order(377) = 0.037d0
    test_target%bond_order(378) = 0.283d0
    test_target%nerror = 0

    call mopac_scf_f(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_crystal1

subroutine test_output(name, input, target, output, nfail)
    use mopac_api_f
    implicit none
    character(50), intent(in) :: name
    type(mopac_system_f), intent(in) :: input
    type(mopac_properties_f), intent(in) :: target, output
    integer, intent(inout) :: nfail
    integer :: i, j
    double precision :: heat_tol, coord_tol, deriv_tol, freq_tol, charge_tol, stress_tol
    double precision, dimension (:,:), allocatable :: fmat1, fmat2, fmat3
    heat_tol = 1d-3
    coord_tol = 1d-3
    deriv_tol = 1d-1
    freq_tol = 3d1
    charge_tol = 1d-3
    stress_tol = 1d-3
    ! compare heat
    if(abs(target%heat - output%heat) > heat_tol) then
        nfail = nfail + 1
        write(*,*) "heat mismatch in test '", trim(name), "':", target%heat, "vs", output%heat
    end if
    ! compare dipole
    do i = 1, 3
        if(abs(target%dipole(i) - output%dipole(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "dipole(", i, ") mismatch in test '", trim(name), "':", target%dipole(i), &
            "vs", output%dipole(i)
        end if
    end do
    ! compare charge
    do i = 1, input%natom
        if(abs(target%charge(i) - output%charge(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "charge(", i, ") mismatch in test '", trim(name), "':", target%charge(i), &
            "vs", output%charge(i)
        end if
    end do
    ! compare updated coordinates
    do i = 1, 3*input%natom_move
        if(abs(target%coord_update(i) - output%coord_update(i)) > coord_tol) then
            nfail = nfail + 1
            write(*,*) "coord_update(", i, ") mismatch in test '", trim(name), "':", &
            target%coord_update(i), "vs", output%coord_update(i)
        end if
    end do
    ! compare coordinate gradients
    do i = 1, 3*input%natom_move
        if(abs(target%coord_deriv(i) - output%coord_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "coord_deriv(", i, ") mismatch in test '", trim(name), "':", &
            target%coord_deriv(i), "vs", output%coord_deriv(i)
        end if
    end do
    ! compare vibrational properties
    if(allocated(target%freq) .neqv. allocated(output%freq)) then
        nfail = nfail + 1
        write(*,*) "freq allocation mistmatch in test '", trim(name), "':", allocated(target%freq), &
        "vs", allocated(output%freq)
    else if(allocated(target%freq) .eqv. .true.) then
        do i = 1, 3*input%natom_move
            if(abs(target%freq(i) - output%freq(i)) > freq_tol) then
                nfail = nfail + 1
                write(*,*) "freq(", i, ") mismatch in test '", trim(name), "':", target%freq(i), &
                "vs", output%freq(i)
            end if
        end do
        allocate(fmat1(3*input%natom_move,3*input%natom_move))
        allocate(fmat2(3*input%natom_move,3*input%natom_move))
        allocate(fmat3(3*input%natom_move,3*input%natom_move))
        do i=1, 3*input%natom_move
            fmat3(:,i) = target%disp(:,i)*target%freq(i)
        end do
        fmat1 = matmul(fmat3,transpose(fmat3))
        do i=1, 3*input%natom_move
            fmat3(:,i) = output%disp(:,i)*output%freq(i)
        end do
        fmat2 = matmul(fmat3,transpose(fmat3))
        do i=1, 3*input%natom_move
            do j=1, 3*input%natom_move
                if(abs(fmat1(i,j) - fmat2(i,j)) > freq_tol*freq_tol) then
                    nfail = nfail + 1
                    write(*,*) "fmat(", i, ",", j, ") mismatch in test '", trim(name), "':", fmat1(i,j), &
                    "vs", fmat2(i,j)
                end if
            end do
        end do
        deallocate(fmat1, fmat2, fmat3)
    end if
    ! compare bond-order matrices
    do i = 1, input%natom+1
        if(target%bond_index(i) /= output%bond_index(i)) then
            nfail = nfail + 1
            write(*,*) "bond_index(", i, ") mismatch in test '", trim(name), "':", target%bond_index(i), "vs", output%bond_index(i)
        end if
    end do
    do i = 1, target%bond_index(input%natom+1)-1
        if(target%bond_atom(i) /= output%bond_atom(i)) then
            nfail = nfail + 1
            write(*,*) "bond_atom(", i, ") mismatch in test '", trim(name), "':", target%bond_atom(i), "vs", output%bond_atom(i)
        end if
        if(abs(target%bond_order(i) - output%bond_order(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "bond_order(", i, ") mismatch in test '", trim(name), "':", target%bond_order(i), "vs", output%bond_order(i)
        end if
    end do
    ! compare updated coordinates
    do i = 1, 3*input%nlattice_move
        if(abs(target%lattice_update(i) - output%lattice_update(i)) > coord_tol) then
            nfail = nfail + 1
            write(*,*) "lattice_update(", i, ") mismatch in test '", trim(name), "':", &
            target%lattice_update(i), "vs", output%lattice_update(i)
        end if
    end do
    ! compare coordinate gradients
    do i = 1, 3*input%nlattice_move
        if(abs(target%lattice_deriv(i) - output%lattice_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "lattice_deriv(", i, ") mismatch in test '", trim(name), "':", &
            target%lattice_deriv(i), "vs", output%lattice_deriv(i)
        end if
    end do
    ! compare stress
    do i = 1, 6
        if(abs(target%stress(i) - output%stress(i)) > stress_tol) then
            nfail = nfail + 1
            write(*,*) "stress(", i, ") mismatch in test '", trim(name), "':", target%stress(i), &
            "vs", output%stress(i)
        end if
    end do
    ! compare error status
    if(target%nerror /= output%nerror) then
        nfail = nfail + 1
        write(*,*) "nerror mismatch in test '", trim(name), "':", target%nerror, "vs", output%nerror
    end if
end subroutine test_output
