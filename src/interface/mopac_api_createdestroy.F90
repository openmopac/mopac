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

! Diskless/stateless Application Programming Interface (API) to core MOPAC operations
submodule (mopac_api) mopac_api_createdestroy
  implicit none

contains

  ! allocate memory & initialize mopac_system
  module subroutine create_mopac_system(natom, natom_move, charge, spin, model, & 
                                        epsilon, atom, coord, nlattice, nlattice_move, &
                                        pressure, lattice, tolerance, max_time, &
                                        system) bind(c)
    integer(c_int), intent(in) :: natom
    integer(c_int), intent(in) :: natom_move
    integer(c_int), intent(in) :: charge
    integer(c_int), intent(in) :: spin
    integer(c_int), intent(in) :: model
    real(c_double), intent(in) :: epsilon
    integer(c_int), dimension(natom), intent(in) :: atom
    real(c_double), dimension(3*natom), intent(in) :: coord
    integer(c_int), intent(in) :: nlattice
    integer(c_int), intent(in) :: nlattice_move
    real(c_double), intent(in) :: pressure
    real(c_double), dimension(3*nlattice), intent(in) :: lattice
    real(c_double), intent(in) :: tolerance
    integer(c_int), intent(in) :: max_time
    type(mopac_system), intent(out) :: system
    integer :: status
    integer(c_int), pointer :: iptr(:)
    real(c_double), pointer :: rptr(:)
    system%natom = natom
    system%natom_move = natom_move
    system%charge = charge
    system%spin = spin
    system%model = model
    system%epsilon = epsilon
    if (natom > 0) then
      allocate(iptr(natom), stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to allocate memory in CREATE_MOPAC_SYSTEM"
        stop 1
      end if  
      iptr(1:natom) = atom(1:natom)
      system%atom = c_loc(iptr)
      allocate(rptr(3*natom), stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to allocate memory in CREATE_MOPAC_SYSTEM"
        stop 1
      end if  
      rptr(1:3*natom) = coord(1:3*natom)
      system%coord = c_loc(rptr)
    end if
    system%nlattice = nlattice
    system%nlattice_move = nlattice_move
    system%pressure = pressure
    if (nlattice > 0) then
      allocate(rptr(3*nlattice), stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to allocate memory in CREATE_MOPAC_SYSTEM"
        stop 1
      end if  
      rptr(1:3*nlattice) = lattice(1:3*nlattice)
      system%lattice = c_loc(rptr)
    end if
    system%tolerance = tolerance
    system%max_time = max_time
  end subroutine create_mopac_system

  ! deallocate memory in mopac_system
  module subroutine destroy_mopac_system(system) bind(c)
    type(mopac_system), intent(in) :: system
    integer(c_int), pointer :: iptr(:)
    real(c_double), pointer :: rptr(:)
    if (system%natom > 0) then
      call c_f_pointer(system%atom, iptr, [system%natom])
      deallocate(iptr)
      call c_f_pointer(system%coord, rptr, [3*system%natom])
      deallocate(rptr)
    end if
    if (system%nlattice > 0) then
      call c_f_pointer(system%lattice, rptr, [3*system%nlattice])
      deallocate(rptr)
    end if
  end subroutine destroy_mopac_system

  ! deallocate memory in mopac_properties
  module subroutine destroy_mopac_properties(system, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mopac_properties), intent(in) :: properties
    integer :: i, j
    integer(c_int), pointer :: iptr(:)
    real(c_double), pointer :: rptr(:)
    character(kind=c_char), pointer :: cptr(:)
    type(c_ptr), pointer :: pptr(:)
    if (system%natom > 0) then
      call c_f_pointer(properties%charge, rptr, [system%natom])
      deallocate(rptr)
      call c_f_pointer(properties%bond_index, iptr, [system%natom+1])
      i = iptr(system%natom+1)
      deallocate(iptr)
      call c_f_pointer(properties%bond_atom, iptr, [i])
      deallocate(iptr)
      call c_f_pointer(properties%bond_order, rptr, [i])
      deallocate(rptr)
    end if
    if (system%natom_move > 0) then
      call c_f_pointer(properties%coord_update, rptr, [3*system%natom_move])
      deallocate(rptr)
      call c_f_pointer(properties%coord_deriv, rptr, [3*system%natom_move])
      deallocate(rptr)
    end if
    if (c_associated(properties%freq)) then
      call c_f_pointer(properties%freq, rptr, [3*system%natom_move])
      deallocate(rptr)
    end if
    if (c_associated(properties%disp)) then
      call c_f_pointer(properties%disp, rptr, [3*system%natom_move*3*system%natom_move])
      deallocate(rptr)
    end if
    if (system%nlattice_move > 0) then
      call c_f_pointer(properties%lattice_update, rptr, [3*system%nlattice_move])
      deallocate(rptr)
      call c_f_pointer(properties%lattice_deriv, rptr, [3*system%nlattice_move])
      deallocate(rptr)
    end if
    if (properties%nerror > 0) then
      call c_f_pointer(properties%error_msg, pptr, [properties%nerror])
      do i=1, properties%nerror
        call c_f_pointer(pptr(i), cptr, [120])
        do j=1, 120 ! using the hard-coded max size of MOPAC's error messages
          if (cptr(j) == c_null_char) exit
        end do
        call c_f_pointer(pptr(i), cptr, [j])
        deallocate(cptr)
      end do
      deallocate(pptr)
    end if
  end subroutine destroy_mopac_properties

  ! deallocate memory in mopac_state
  module subroutine destroy_mopac_state(state) bind(c)
    type(mopac_state), intent(in) :: state
    real(c_double), pointer :: rptr(:)
    if (state%mpack /= 0) then
      call c_f_pointer(state%pa, rptr, [state%mpack])
      deallocate(rptr)
      call c_f_pointer(state%pb, rptr, [state%mpack])
      deallocate(rptr)
    end if
  end subroutine destroy_mopac_state

  ! deallocate memory in mozyme_state
  module subroutine destroy_mozyme_state(state) bind(c)
    type(mozyme_state), intent(in) :: state
    integer(c_int), pointer :: iptr(:)
    real(c_double), pointer :: rptr(:)
    if (state%numat /= 0) then
      call c_f_pointer(state%nbonds, iptr, [state%numat])
      deallocate(iptr)
      call c_f_pointer(state%ibonds, iptr, [9*state%numat])
      deallocate(iptr)
      call c_f_pointer(state%iorbs, iptr, [state%numat])
      deallocate(iptr)
      call c_f_pointer(state%ncf, iptr, [state%noccupied])
      deallocate(iptr)
      call c_f_pointer(state%nce, iptr, [state%nvirtual])
      deallocate(iptr)
      call c_f_pointer(state%icocc, iptr, [state%icocc_dim])
      deallocate(iptr)
      call c_f_pointer(state%icvir, iptr, [state%icvir_dim])
      deallocate(iptr)
      call c_f_pointer(state%cocc, rptr, [state%cocc_dim])
      deallocate(rptr)
      call c_f_pointer(state%cvir, rptr, [state%cvir_dim])
      deallocate(rptr)
    end if
  end subroutine destroy_mozyme_state

end submodule mopac_api_createdestroy
