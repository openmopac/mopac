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
    integer :: test_type(1)
    type(mopac_system) :: test_in(1)
    type(mopac_state) :: test_restore(1)
    type(mopac_properties) :: test_target(1)
    type(mopac_properties) :: test_out
    integer :: i, ntest, nfail
    ntest = 1
    nfail = 0

    ! assign MOPAC testing information:

    ! SCF calculation of H2O
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
    test_target(1)%coord_deriv(1) = 2.277043d0
    test_target(1)%coord_deriv(2) = 2.711610d0
    test_target(1)%coord_deriv(3) = 0.0d0
    test_target(1)%coord_deriv(4) = -2.277043d0
    test_target(1)%coord_deriv(5) = 2.742432d0
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
    test_target(1)%status = 0

    ! bond-order matrix in compressed sparse column (CSC) matrix format
    ! with insignificant bond orders (<0.001) truncated
    ! > first index of each atom in CSC bond-order matrix [natom+1]
!    integer, dimension (:), allocatable :: bond_index
    ! > list of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
!    integer, dimension (:), allocatable :: bond_atom
    ! > bond order of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
!    double precision, dimension (:), allocatable :: bond_order

    ! run tests
    do i = 1, ntest
        select case (test_type(i))
            case (1)
                call mopac_scf(test_in(i), test_restore(i), test_out)
            case (2)
            case (3)
            case (4)
            case (5)
        end select
        ! test common output
        call test_output(i, test_target(i), test_out, nfail)
        ! test task-specific output ...
    end do

    call exit(nfail)

end program mopac_api_test

subroutine test_output(num, target, output, nfail)
    use mopac_api
    use molkst_C, only : keywrd
    implicit none
    integer, intent(in) :: num
    type(mopac_properties), intent(in) :: target, output
    integer, intent(inout) :: nfail
    logical :: fail
    fail = .false.
    if(abs(target%heat - output%heat) > 1d-3) then
        fail = .true.
        write(*,*) "heat mismatch in test", num, ":", target%heat, "vs", output%heat
    end if
    if(target%status /= output%status) then
        fail = .true.
        write(*,*) "status mismatch in test", num, ":", target%status, "vs", output%status
    end if
    if (fail) nfail = nfail + 1
end subroutine test_output