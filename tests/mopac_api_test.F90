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
    use mopac_api
    implicit none
    integer :: nfail
    nfail = 0

    call test_mopac_scf1(nfail)
    call test_mopac_relax1(nfail)
    call test_mopac_vibe1(nfail)
    call test_mozyme_scf1(nfail)
    call test_mozyme_relax1(nfail)
    call test_mozyme_vibe1(nfail)
    call exit(nfail)
end program mopac_api_test

subroutine test_mopac_scf1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mopac_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .false.
    test_target%nerror = 0

    call mopac_scf(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_scf1

subroutine test_mopac_relax1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mopac_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .false.
    test_target%nerror = 0

    call mopac_relax(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_relax1

subroutine test_mopac_vibe1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mopac_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .true.
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

    call mopac_vibe(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_vibe1

subroutine test_mozyme_scf1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mozyme_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .false.
    test_target%nerror = 0

    call mozyme_scf(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_scf1

subroutine test_mozyme_relax1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mozyme_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .false.
    test_target%nerror = 0

    call mozyme_relax(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_relax1


subroutine test_mozyme_vibe1(nfail)
    use mopac_api
    implicit none
    integer, intent(inout) :: nfail
    type(mopac_system) :: test_in
    type(mozyme_state) :: test_restore
    type(mopac_properties) :: test_target
    type(mopac_properties) :: test_out
    character(20) :: test_name
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
    test_target%calc_vibe = .true.
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

    call mozyme_vibe(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mozyme_vibe1

subroutine test_output(name, input, target, output, nfail)
    use mopac_api
    implicit none
    character(20), intent(in) :: name
    type(mopac_system), intent(in) :: input
    type(mopac_properties), intent(in) :: target, output
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
        write(*,*) "heat mismatch in test '", name, "':", target%heat, "vs", output%heat
    end if
    ! compare dipole
    do i = 1, 3
        if(abs(target%dipole(i) - output%dipole(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "dipole(", i, ") mismatch in test '", name, "':", target%dipole(i), &
            "vs", output%dipole(i)
        end if
    end do
    ! compare charge
    do i = 1, input%natom
        if(abs(target%charge(i) - output%charge(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "charge(", i, ") mismatch in test '", name, "':", target%charge(i), &
            "vs", output%charge(i)
        end if
    end do
    ! compare updated coordinates
    do i = 1, 3*input%natom_move
        if(abs(target%coord_update(i) - output%coord_update(i)) > coord_tol) then
            nfail = nfail + 1
            write(*,*) "coord_update(", i, ") mismatch in test '", name, "':", &
            target%coord_update(i), "vs", output%coord_update(i)
        end if
    end do
    ! compare coordinate gradients
    do i = 1, 3*input%natom_move
        if(abs(target%coord_deriv(i) - output%coord_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "coord_deriv(", i, ") mismatch in test '", name, "':", &
            target%coord_deriv(i), "vs", output%coord_deriv(i)
        end if
    end do
    ! compare vibrational properties
    if(target%calc_vibe .neqv. output%calc_vibe) then
        nfail = nfail + 1
        write(*,*) "calc_vibe mistmatch in test '", name, "':", target%calc_vibe, &
        "vs", output%calc_vibe
    else if(target%calc_vibe .eqv. .true.) then
        do i = 1, 3*input%natom_move
            if(abs(target%freq(i) - output%freq(i)) > freq_tol) then
                nfail = nfail + 1
                write(*,*) "freq(", i, ") mismatch in test '", name, "':", target%freq(i), &
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
                    write(*,*) "fmat(", i, ",", j, ") mismatch in test '", name, "':", fmat1(i,j), &
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
            write(*,*) "bond_index(", i, ") mismatch in test '", name, "':", target%bond_index(i), "vs", output%bond_index(i)
        end if
    end do
    do i = 1, target%bond_index(input%natom+1)-1
        if(target%bond_atom(i) /= output%bond_atom(i)) then
            nfail = nfail + 1
            write(*,*) "bond_atom(", i, ") mismatch in test '", name, "':", target%bond_atom(i), "vs", output%bond_atom(i)
        end if
        if(abs(target%bond_order(i) - output%bond_order(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "bond_order(", i, ") mismatch in test '", name, "':", target%bond_order(i), "vs", output%bond_order(i)
        end if
    end do
    ! compare updated coordinates
    do i = 1, 3*input%nlattice_move
        if(abs(target%lattice_update(i) - output%lattice_update(i)) > coord_tol) then
            nfail = nfail + 1
            write(*,*) "lattice_update(", i, ") mismatch in test '", name, "':", &
            target%lattice_update(i), "vs", output%lattice_update(i)
        end if
    end do
    ! compare coordinate gradients
    do i = 1, 3*input%nlattice_move
        if(abs(target%lattice_deriv(i) - output%lattice_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "lattice_deriv(", i, ") mismatch in test '", name, "':", &
            target%lattice_deriv(i), "vs", output%lattice_deriv(i)
        end if
    end do
    ! compare stress
    do i = 1, 6
        if(abs(target%stress(i) - output%stress(i)) > stress_tol) then
            nfail = nfail + 1
            write(*,*) "stress(", i, ") mismatch in test '", name, "':", target%stress(i), &
            "vs", output%stress(i)
        end if
    end do
    ! compare error status
    if(target%nerror /= output%nerror) then
        nfail = nfail + 1
        write(*,*) "nerror mismatch in test '", name, "':", target%nerror, "vs", output%nerror
    end if
end subroutine test_output
