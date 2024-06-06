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
    integer :: test_type(2)
    type(mopac_system) :: test_in(2), target_coord(2)
    type(mopac_state) :: test_restore(2)
    type(mopac_properties) :: test_target(2)
    type(mopac_properties) :: test_out
    double precision, dimension(:), allocatable :: coord_out
    integer :: i, ntest, nfail
    ntest = 2
    nfail = 0

    ! assign MOPAC testing information:

    ! 1 - SCF calculation of H2O
    test_type(1) = 1
    test_in(1)%natom = 3
    allocate(test_in(1)%atom(3))
    test_in(1)%atom(1) = 1
    test_in(1)%atom(2) = 1
    test_in(1)%atom(3) = 8
    allocate(test_in(1)%coord(3*3))
    test_in(1)%coord(1) = 0.76d0
    test_in(1)%coord(2) = 0.59d0
    test_in(1)%coord(3) = 0.0d0
    test_in(1)%coord(4) = -0.76d0
    test_in(1)%coord(5) = 0.59d0
    test_in(1)%coord(6) = 0.0d0
    test_in(1)%coord(7) = 0.0d0
    test_in(1)%coord(8) = 0.0d0
    test_in(1)%coord(9) = 0.0d0
    test_target(1)%heat = -57.76975d0
    allocate(test_target(1)%coord_deriv(3*3))
    test_target(1)%coord_deriv(1) = 2.307865d0
    test_target(1)%coord_deriv(2) = 2.742432d0
    test_target(1)%coord_deriv(3) = 0.0d0
    test_target(1)%coord_deriv(4) = -2.307865d0
    test_target(1)%coord_deriv(5) = 2.711610d0
    test_target(1)%coord_deriv(6) = 0.0d0
    test_target(1)%coord_deriv(7) = 0.0d0
    test_target(1)%coord_deriv(8) = -5.454042d0
    test_target(1)%coord_deriv(9) = 0.0d0
    allocate(test_target(1)%charge(3))
    test_target(1)%charge(1) = 0.322260d0
    test_target(1)%charge(2) = 0.322260d0
    test_target(1)%charge(3) = -0.644520d0
    test_target(1)%dipole(1) = 0.0d0
    test_target(1)%dipole(2) = 2.147d0
    test_target(1)%dipole(3) = 0.0d0
    test_target(1)%stress(:) = 0.0d0
    allocate(test_target(1)%bond_index(4))
    test_target(1)%bond_index(1) = 1
    test_target(1)%bond_index(2) = 3
    test_target(1)%bond_index(3) = 5
    test_target(1)%bond_index(4) = 8
    allocate(test_target(1)%bond_atom(7))
    test_target(1)%bond_atom(1) = 1
    test_target(1)%bond_atom(2) = 3
    test_target(1)%bond_atom(3) = 2
    test_target(1)%bond_atom(4) = 3
    test_target(1)%bond_atom(5) = 1
    test_target(1)%bond_atom(6) = 2
    test_target(1)%bond_atom(7) = 3
    allocate(test_target(1)%bond_order(7))
    test_target(1)%bond_order(1) = 0.896
    test_target(1)%bond_order(2) = 0.895
    test_target(1)%bond_order(3) = 0.896
    test_target(1)%bond_order(4) = 0.895
    test_target(1)%bond_order(5) = 0.895
    test_target(1)%bond_order(6) = 0.895
    test_target(1)%bond_order(7) = 1.791
    test_target(1)%status = 0

    ! 2 - geometry relaxation of H2O
    test_type(2) = 2
    test_in(2)%natom = 3
    allocate(test_in(2)%atom(3))
    test_in(2)%atom(1) = 1
    test_in(2)%atom(2) = 1
    test_in(2)%atom(3) = 8
    allocate(test_in(2)%coord(3*3))
    test_in(2)%coord(1) = 0.8d0
    test_in(2)%coord(2) = 0.4d0
    test_in(2)%coord(3) = 0.0d0
    test_in(2)%coord(4) = -0.8d0
    test_in(2)%coord(5) = 0.4d0
    test_in(2)%coord(6) = 0.0d0
    test_in(2)%coord(7) = 0.0d0
    test_in(2)%coord(8) = 0.0d0
    test_in(2)%coord(9) = 0.0d0
    target_coord(2)%natom = 3
    allocate(target_coord(2)%coord(3*3))
    target_coord(2)%coord(1) = 0.758862835d0
    target_coord(2)%coord(2) = 0.486805290d0
    target_coord(2)%coord(3) = 0.0d0
    target_coord(2)%coord(4) = -0.759709263d0
    target_coord(2)%coord(5) = 0.485142519d0
    target_coord(2)%coord(6) = 0.0d0
    target_coord(2)%coord(7) = 0.000846434d0
    target_coord(2)%coord(8) = -0.093004709d0
    target_coord(2)%coord(9) = 0.0d0
    test_target(2)%heat = -57.79952d0
    allocate(test_target(2)%coord_deriv(3*3))
    test_target(2)%coord_deriv(1) = -0.565914d0
    test_target(2)%coord_deriv(2) = -0.294814d0
    test_target(2)%coord_deriv(3) = 0.0d0
    test_target(2)%coord_deriv(4) = 0.063505d0
    test_target(2)%coord_deriv(5) = 0.058243d0
    test_target(2)%coord_deriv(6) = 0.0d0
    test_target(2)%coord_deriv(7) = 0.502409d0
    test_target(2)%coord_deriv(8) = 0.236571d0
    test_target(2)%coord_deriv(9) = 0.0d0
    allocate(test_target(2)%charge(3))
    test_target(2)%charge(1) = 0.324544d0
    test_target(2)%charge(2) = 0.324720d0
    test_target(2)%charge(3) = -0.649264d0
    test_target(2)%dipole(1) = -0.005d0
    test_target(2)%dipole(2) = 2.129d0
    test_target(2)%dipole(3) = 0.0d0
    test_target(2)%stress(:) = 0.0d0
    allocate(test_target(2)%bond_index(4))
    test_target(2)%bond_index(1) = 1
    test_target(2)%bond_index(2) = 3
    test_target(2)%bond_index(3) = 5
    test_target(2)%bond_index(4) = 8
    allocate(test_target(2)%bond_atom(7))
    test_target(2)%bond_atom(1) = 1
    test_target(2)%bond_atom(2) = 3
    test_target(2)%bond_atom(3) = 2
    test_target(2)%bond_atom(4) = 3
    test_target(2)%bond_atom(5) = 1
    test_target(2)%bond_atom(6) = 2
    test_target(2)%bond_atom(7) = 3
    allocate(test_target(2)%bond_order(7))
    test_target(2)%bond_order(1) = 0.895
    test_target(2)%bond_order(2) = 0.894
    test_target(2)%bond_order(3) = 0.895
    test_target(2)%bond_order(4) = 0.894
    test_target(2)%bond_order(5) = 0.894
    test_target(2)%bond_order(6) = 0.894
    test_target(2)%bond_order(7) = 1.788
    test_target(2)%status = 0

    ! run tests
    do i = 1, ntest
        allocate(coord_out(3*test_in(i)%natom))
        select case (test_type(i))
            case (1)
                call mopac_scf(test_in(i), test_restore(i), test_out)
            case (2)
                call mopac_relax(test_in(i), test_restore(i), test_out, coord_out)
            case (3)
            case (4)
            case (5)
        end select
        ! test common output
        call test_output(i, test_in(i), test_target(i), test_out, nfail)
        ! test geometry
        if(test_type(i) == 2 .or. test_type(i) == 5) then
            call test_coord(i, target_coord(i), coord_out, nfail)
        end if
        deallocate(coord_out)
    end do

    call exit(nfail)

end program mopac_api_test

subroutine test_output(num, input, target, output, nfail)
    use mopac_api
    implicit none
    integer, intent(in) :: num
    type(mopac_system), intent(in) :: input
    type(mopac_properties), intent(in) :: target, output
    integer, intent(inout) :: nfail
    integer :: i
    double precision :: heat_tol, deriv_tol, charge_tol, stress_tol
    heat_tol = 1d-3
    deriv_tol = 1d-1
    charge_tol = 1d-3
    stress_tol = 1d-3
    ! compare heats
    if(abs(target%heat - output%heat) > heat_tol) then
        nfail = nfail + 1
        write(*,*) "heat mismatch in test", num, ":", target%heat, "vs", output%heat
    end if
    ! compare coordinate gradients
    do i = 1, 3*input%natom
        if(abs(target%coord_deriv(i) - output%coord_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "coord_deriv(", i, ") mismatch in test", num, ":", target%coord_deriv(i), "vs", output%coord_deriv(i)
        end if
    end do
    ! compare lattice gradients
    do i = 1, 3*input%nlattice
        if(abs(target%lattice_deriv(i) - output%lattice_deriv(i)) > deriv_tol) then
            nfail = nfail + 1
            write(*,*) "lattice_deriv(", i, ") mismatch in test", num, ":", target%lattice_deriv(i), "vs", output%lattice_deriv(i)
        end if
    end do
    ! compare charge
    do i = 1, input%natom
        if(abs(target%charge(i) - output%charge(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "charge(", i, ") mismatch in test", num, ":", target%charge(i), "vs", output%charge(i)
        end if
    end do
    ! compare dipole
    do i = 1, 3
        if(abs(target%dipole(i) - output%dipole(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "dipole(", i, ") mismatch in test", num, ":", target%dipole(i), "vs", output%dipole(i)
        end if
    end do
    ! compare stress
    do i = 1, 6
        if(abs(target%stress(i) - output%stress(i)) > stress_tol) then
            nfail = nfail + 1
            write(*,*) "stress(", i, ") mismatch in test", num, ":", target%stress(i), "vs", output%stress(i)
        end if
    end do
    ! compare bond-order matrices
    do i = 1, input%natom+1
        if(target%bond_index(i) /= output%bond_index(i)) then
            nfail = nfail + 1
            write(*,*) "bond_index(", i, ") mismatch in test", num, ":", target%bond_index(i), "vs", output%bond_index(i)
        end if
    end do
    do i = 1, target%bond_index(input%natom+1)-1
        if(target%bond_atom(i) /= output%bond_atom(i)) then
            nfail = nfail + 1
            write(*,*) "bond_atom(", i, ") mismatch in test", num, ":", target%bond_atom(i), "vs", output%bond_atom(i)
        end if
        if(abs(target%bond_order(i) - output%bond_order(i)) > charge_tol) then
            nfail = nfail + 1
            write(*,*) "bond_order(", i, ") mismatch in test", num, ":", target%bond_order(i), "vs", output%bond_order(i)
        end if
    end do
    ! compare status
    if(target%status /= output%status) then
        nfail = nfail + 1
        write(*,*) "status mismatch in test", num, ":", target%status, "vs", output%status
    end if
end subroutine test_output

subroutine test_coord(num, target, output, nfail)
    use mopac_api
    implicit none
    integer, intent(in) :: num
    type(mopac_system), intent(in) :: target
    double precision, dimension(3*target%natom), intent(in) :: output
    integer, intent(inout) :: nfail
    integer :: i
    double precision :: coord_tol
    coord_tol = 1d-3
    ! compare coordinates
    do i = 1, 3*target%natom
        if(abs(target%coord(i) - output(i)) > coord_tol) then
            nfail = nfail + 1
            write(*,*) "coord(", i, ") mismatch in test", num, ":", target%coord(i), "vs", output(i)
        end if
    end do
end subroutine test_coord