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
!    call test_mopac_vibe1(nfail)
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
    allocate(test_in%move_atom(3))
    test_in%move_atom = .true.
    test_target%heat = -57.76975d0
    test_target%natom_move = 3
    allocate(test_target%atom_move(3))
    do i=1, 3
      test_target%atom_move(i) = i
    end do
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
    test_target%nlattice_move = 0
    test_target%status = 0

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
    allocate(test_in%move_atom(3))
    test_in%move_atom = .true.
    test_target%heat = -57.79952d0
    test_target%natom_move = 3
    allocate(test_target%atom_move(3))
    do i=1, 3
      test_target%atom_move(i) = i
    end do
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
    test_target%nlattice_move = 0
    test_target%status = 0

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
    double precision, dimension(3) :: frequency_target, frequency_out
    double precision, dimension(9,3) :: displacement_target, displacement_out
    character(20) :: test_name
    test_name = 'H2O vibe'

    ! 1 - geometry relaxation of H2O
    test_in%natom = 3
    allocate(test_in%atom(3))
    test_in%atom(1) = 1
    test_in%atom(2) = 1
    test_in%atom(3) = 8
    allocate(test_in%coord(3*3))
    test_in%coord(1) = 0.758862835d0
    test_in%coord(2) = 0.486805290d0
    test_in%coord(3) = 0.0d0
    test_in%coord(4) = -0.759709263d0
    test_in%coord(5) = 0.485142519d0
    test_in%coord(6) = 0.0d0
    test_in%coord(7) = 0.000846434d0
    test_in%coord(8) = -0.093004709d0
    test_in%coord(9) = 0.0d0
    test_target%heat = -57.79952d0
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
    test_target%status = 0
    frequency_target(1) = 1394.0d0
    frequency_target(2) = 2811.2d0
    frequency_target(3) = 2858.2d0
    displacement_target(1,1) = -0.3145d0
    displacement_target(2,1) = 0.5972d0
    displacement_target(3,1) = 0.0000d0
    displacement_target(4,1) = 0.3123d0
    displacement_target(5,1) = 0.5976d0
    displacement_target(6,1) = 0.0000d0
    displacement_target(7,1) = 0.0006d0
    displacement_target(8,1) = -0.2999d0
    displacement_target(9,1) = 0.0000d0
    displacement_target(1,2) = 0.5259d0
    displacement_target(2,2) = 0.4066d0
    displacement_target(3,2) = 0.0000d0
    displacement_target(4,2) = 0.5561d0
    displacement_target(5,2) = -0.4184d0
    displacement_target(6,2) = 0.0000d0
    displacement_target(7,2) = -0.2716d0
    displacement_target(8,2) = 0.0030d0
    displacement_target(9,2) = 0.0000d0
    displacement_target(1,3) = 0.6456d0
    displacement_target(2,3) = 0.3055d0
    displacement_target(3,3) = 0.0000d0
    displacement_target(4,3) = -0.6218d0
    displacement_target(5,3) = 0.2850d0
    displacement_target(6,3) = 0.0000d0
    displacement_target(7,3) = -0.0060d0
    displacement_target(8,3) = -0.1482d0
    displacement_target(9,3) = 0.0000d0

    call mopac_vibe(test_in, test_restore, test_out)
    call test_output(test_name, test_in, test_target, test_out, nfail)
end subroutine test_mopac_vibe1

subroutine test_output(name, input, target, output, nfail)
    use mopac_api
    implicit none
    character(20), intent(in) :: name
    type(mopac_system), intent(in) :: input
    type(mopac_properties), intent(in) :: target, output
    integer, intent(inout) :: nfail
    integer :: i
    double precision :: heat_tol, coord_tol, deriv_tol, freq_tol, charge_tol, stress_tol
    heat_tol = 1d-3
    coord_tol = 1d-3
    deriv_tol = 1d-1
    freq_tol = 1d1
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
    ! compare moveable atoms
    if(target%natom_move /= output%natom_move) then
        nfail = nfail + 1
        write(*,*) "natom_move mistmatch in test '", name, "':", target%natom_move, &
        "vs", output%natom_move
    else
        ! compare updated coordinates
        do i = 1, 3*target%natom_move
            if(abs(target%coord_update(i) - output%coord_update(i)) > coord_tol) then
                nfail = nfail + 1
                write(*,*) "coord_update(", i, ") mismatch in test '", name, "':", &
                target%coord_update(i), "vs", output%coord_update(i)
            end if
        end do
        ! compare coordinate gradients
        do i = 1, 3*target%natom_move
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
            do i = 1, 3*target%natom_move
                if(abs(target%freq(i) - output%freq(i)) > freq_tol) then
                    nfail = nfail + 1
                    write(*,*) "freq(", i, ") mismatch in test '", name, "':", target%freq(i), &
                    "vs", output%freq(i)
                end if
            end do
! TODO: add test of vibrational modes (disp) ...
        end if
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
    ! compare moveable lattice vectors
    if(target%nlattice_move /= output%nlattice_move) then
        nfail = nfail + 1
        write(*,*) "nlattice_move mistmatch in test '", name, "':", target%nlattice_move, "vs", output%nlattice_move
    else
        ! compare updated coordinates
        do i = 1, 3*target%nlattice_move
            if(abs(target%lattice_update(i) - output%lattice_update(i)) > coord_tol) then
                nfail = nfail + 1
                write(*,*) "lattice_update(", i, ") mismatch in test '", name, "':", &
                target%lattice_update(i), "vs", output%lattice_update(i)
            end if
        end do
        ! compare coordinate gradients
        do i = 1, 3*target%nlattice_move
            if(abs(target%lattice_deriv(i) - output%lattice_deriv(i)) > deriv_tol) then
                nfail = nfail + 1
                write(*,*) "lattice_deriv(", i, ") mismatch in test '", name, "':", &
                target%lattice_deriv(i), "vs", output%lattice_deriv(i)
            end if
        end do
    end if
    ! compare stress
    do i = 1, 6
        if(abs(target%stress(i) - output%stress(i)) > stress_tol) then
            nfail = nfail + 1
            write(*,*) "stress(", i, ") mismatch in test '", name, "':", target%stress(i), &
            "vs", output%stress(i)
        end if
    end do
    ! compare status
    if(target%status /= output%status) then
        nfail = nfail + 1
        write(*,*) "status mismatch in test '", name, "':", target%status, "vs", output%status
    end if
end subroutine test_output
