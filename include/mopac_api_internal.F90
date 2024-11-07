! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University

! C-bound types for diskless/stateless Application Programming Interface (API) to core MOPAC operations
module mopac_api_c
  use, intrinsic :: iso_c_binding
  implicit none

  ! data that defines the atomistic system and MOPAC job options
  type, bind(c) :: mopac_system
    ! number of atoms
    integer(c_int) :: natom
    ! number of atoms that are allowed to move (first natom_move atoms in array)
    integer(c_int) :: natom_move
    ! net charge
    integer(c_int) :: charge
    ! number of spin excitations, floor[(number of alpha electrons)/2 - (number of beta electrons)/2]
    integer(c_int) :: spin
    ! semiempirical model: PM7 = 0, PM6-D3H4 = 1, PM6-ORG = 2, PM6 = 3, AM1 = 4, RM1 = 5
    integer(c_int) :: model
    ! dielectric constant for COSMO implicit solvent, must be 1 (no solvent) for nlattice > 0
    real(c_double) :: epsilon
    ! atomic number of each atom
    type(c_ptr) :: atom ! integer(c_int)[natom]
    ! (x,y,z) coordinates of each atom (Angstroms)
    type(c_ptr) :: coord ! real(c_double)[3*natom]
    ! number of lattice vectors / translation vectors / periodic dimensions
    integer(c_int) :: nlattice
    ! number of lattice vectors that are allowed to move (first nlattice_move vectors in array)
    integer(c_int) :: nlattice_move
    ! external hydrostatic pressure (Gigapascals)
    real(c_double) :: pressure
    ! (x,y,z) coordinates of each lattice vectors (Angstroms)
    type(c_ptr) :: lattice ! real(c_double)[3*nlattice]
    ! relative numerical tolerances (equivalent to GNORM and RELSCF keyword values)
    real(c_double) :: tolerance
    ! time limit for a MOPAC calculation (seconds)
    integer(c_int) :: max_time
  end type

  ! calculated ground-state properties of an atomistic system and MOPAC job info
  type, bind(c) :: mopac_properties
    ! heat of formation (kcal/mol)
    real(c_double) :: heat
    ! dipole moment vector (Debye)
    real(c_double), dimension (3) :: dipole
    ! atomic partial charges
    type(c_ptr) :: charge ! real(c_double)[natom]
    ! (x,y,z) coordinates of each moveable atom (Angstroms)
    type(c_ptr) :: coord_update ! real(c_double)[3*natom_move]
    ! (x,y,z) heat gradients for each moveable atom (kcal/mol/Angstrom)
    type(c_ptr) :: coord_deriv ! real(c_double)[3*natom_move]
    ! vibrational frequencies of normal modes (1/cm)
    type(c_ptr) :: freq ! real(c_double)[3*natom_move], NULL if unavailable
    ! (x,y,z) displacement vectors of normal modes
    type(c_ptr) :: disp ! real(c_double)[3*natom_move,3*natom_move], NULL if unavailable
    ! bond-order matrix in compressed sparse column (CSC) matrix format (0-based indexing)
    ! with insignificant bond orders (<0.01) truncated
    ! diagonal matrix entries are atomic valencies
    ! > first index of each atom in CSC bond-order matrix
    type(c_ptr) :: bond_index ! integer(c_int)[natom+1]
    ! > list of atoms bonded to each atom in CSC format
    type(c_ptr) :: bond_atom ! integer(c_int)[bond_index(natom+1)]
    ! > bond order of atoms bonded to each atom in CSC format
    type(c_ptr) :: bond_order ! real(c_double)[bond_index(natom+1)]
    ! (x,y,z) coordinates of each moveable lattice vectors (Angstroms)
    type(c_ptr) :: lattice_update ! real(c_double)[3*nlattice_move]
    ! (x,y,z) heat gradients for each moveable lattice vector (kcal/mol/Angstrom)
    type(c_ptr) :: lattice_deriv ! real(c_double)[3*nlattice_move]
    ! stress tensor (Gigapascals) in Voigt form (xx, yy, zz, yz, xz, xy)
    real(c_double), dimension (6) :: stress ! 0 if unavailable
    ! number of MOPAC error messages (negative value indicates that allocation of error_msg failed)
    integer(c_int) :: nerror
    ! text of MOPAC error messages
    type(c_ptr) :: error_msg ! type(c_ptr)[nerror] -> character(kind=c_char)[*]
  end type

  ! data that describes the electronic state using standard molecular orbitals
  type, bind(c) :: mopac_state
    ! MOPAC data format is adapted from molkst_C and Common_arrays_C modules
    ! > number of matrix elements in packed lower triangle matrix format
    integer(c_int) :: mpack ! 0 if state is unavailable
    ! > flag for unrestricted Hartree-Fock ground state (0 == restricted, 1 == unrestricted)
    integer(c_int) :: uhf
    ! > alpha density matrix
    type(c_ptr) :: pa ! real(c_double)[mpack]
    ! > beta density matrix
    type(c_ptr) :: pb ! real(c_double)[mpack], NULL if uhf == 0
  end type

  ! data that describes the electronic state using localized molecular orbitals
  type, bind(c) :: mozyme_state
    ! MOZYME data format is adapted from molkst_C, Common_arrays_C, and MOZYME_C modules
    ! > number of real atoms
    integer(c_int) :: numat ! 0 if state is unavailable
    ! > number of Lewis bonds per real atom
    type(c_ptr) :: nbonds ! integer(c_int)[numat]
    ! > list of Lewis-bonded real atoms for each real atom
    type(c_ptr) :: ibonds ! integer(c_int)[9,numat]
    ! > number of orbitals per real atom
    type(c_ptr) :: iorbs ! integer(c_int)[numat]
    ! > number of occupied molecular orbitals
    integer(c_int) :: noccupied
    ! > number of atoms in each occupied LMO
    type(c_ptr) :: ncf ! integer(c_int)[noccupied]
    ! > number of virtual molecular orbitals
    integer(c_int) :: nvirtual
    ! > number of atoms in each virtual LMO
    type(c_ptr) :: nce ! integer(c_int)[nvirtual]
    ! > size of array icocc
    integer(c_int) :: icocc_dim
    ! > index of each real atom in the occupied LMOs
    type(c_ptr) :: icocc ! integer(c_int)[iccoc_dim]
    ! > size of array icvir
    integer(c_int) :: icvir_dim
    ! > index of each real atom in the virtual LMOs
    type(c_ptr) :: icvir ! integer(c_int)[icvir_dim]
    ! > size of array cocc
    integer(c_int) :: cocc_dim
    ! > atomic orbital coefficients of the occupied LMOs
    type(c_ptr) :: cocc ! real(c_double)[cocc_dim]
    ! > size of array cvir
    integer(c_int) :: cvir_dim
    ! > atomic orbital coefficients of the virtual LMOs
    type(c_ptr) :: cvir ! real(c_double)[cvir_dim]
  end type

end module mopac_api_c

! Implementation of Fortran API wrapper
submodule (mopac_api_f) mopac_api_f_internal
  use, intrinsic :: iso_c_binding
  use mopac_api_c
  implicit none

  ! C-bound library subroutines
  interface

    ! MOPAC electronic ground state calculation
    subroutine mopac_scf(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_scf

    ! MOPAC geometry relaxation
    subroutine mopac_relax(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_relax

    ! MOPAC vibrational calculation
    subroutine mopac_vibe(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_vibe

    ! MOZYME electronic ground state calculation
    subroutine mozyme_scf(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_scf

    ! MOZYME geometry relaxation
    subroutine mozyme_relax(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_relax

    ! MOZYME vibrational calculation
    subroutine mozyme_vibe(system, state, properties) bind(c)
      use mopac_api_c
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_vibe

    ! allocate memory for mopac_state
    subroutine create_mopac_state(state) bind(c)
      use mopac_api_c
      type(mopac_state), intent(inout) :: state
    end subroutine create_mopac_state

    ! allocate memory for mozyme_state
    subroutine create_mozyme_state(state) bind(c)
      use mopac_api_c
      type(mozyme_state), intent(inout) :: state
    end subroutine create_mozyme_state

    ! deallocate memory in mopac_properties (associated system is needed for size info)
    subroutine destroy_mopac_properties(properties) bind(c)
      use mopac_api_c
      type(mopac_properties), intent(in) :: properties
    end subroutine destroy_mopac_properties

    ! deallocate memory in mopac_state
    subroutine destroy_mopac_state(state) bind(c)
      use mopac_api_c
      type(mopac_state), intent(in) :: state
    end subroutine destroy_mopac_state

    ! deallocate memory in mozyme_state
    subroutine destroy_mozyme_state(state) bind(c)
      use mopac_api_c
      type(mozyme_state), intent(in) :: state
    end subroutine destroy_mozyme_state

    ! run MOPAC conventionally from an input file
    subroutine run_mopac_from_input(path_to_file) bind(c)
      use iso_c_binding
      character(kind=c_char), dimension(*), intent(in) :: path_to_file
    end subroutine run_mopac_from_input

  end interface

contains

  module subroutine mopac_scf_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mopac_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mopac_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mopac_state_f2c(state, state_c)
    call mopac_scf(system_c, state_c, properties_c)
    call mopac_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mopac_scf_f

  module subroutine mopac_relax_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mopac_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mopac_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mopac_state_f2c(state, state_c)
    call mopac_relax(system_c, state_c, properties_c)
    call mopac_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mopac_relax_f

  module subroutine mopac_vibe_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mopac_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mopac_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mopac_state_f2c(state, state_c)
    call mopac_vibe(system_c, state_c, properties_c)
    call mopac_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mopac_vibe_f

  module subroutine mozyme_scf_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mozyme_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mozyme_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mozyme_state_f2c(state, state_c)
    call mozyme_scf(system_c, state_c, properties_c)
    call mozyme_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mozyme_scf_f

  module subroutine mozyme_relax_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mozyme_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mozyme_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mozyme_state_f2c(state, state_c)
    call mozyme_relax(system_c, state_c, properties_c)
    call mozyme_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mozyme_relax_f

  module subroutine mozyme_vibe_f(system, state, properties)
    type(mopac_system_f), intent(in) :: system
    type(mozyme_state_f), intent(inout) :: state
    type(mopac_properties_f), intent(out) :: properties
    type(mopac_system) :: system_c
    type(mozyme_state) :: state_c
    type(mopac_properties) :: properties_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    call mopac_system_f2c(system, system_c, iwork, rwork)
    call mozyme_state_f2c(state, state_c)
    call mozyme_vibe(system_c, state_c, properties_c)
    call mozyme_state_c2f(state_c, state)
    call mopac_properties_c2f(system_c, properties_c, properties)
    deallocate(iwork)
    deallocate(rwork)
  end subroutine mozyme_vibe_f

  module subroutine run_mopac_from_input_f(path_to_file)
    character(len=240), intent(in) :: path_to_file
    character(kind=c_char), allocatable :: path_to_file_c(:)
    integer :: i, size, status
    size = len_trim(path_to_file)
    allocate(path_to_file_c(size+1), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    do i=1, size
      path_to_file_c(i) = path_to_file(i:i)
    end do
    path_to_file_c(size+1) = c_null_char
    call run_mopac_from_input(path_to_file_c)
  end subroutine run_mopac_from_input_f

  subroutine mopac_system_f2c(system_f, system_c, iwork, rwork)
    type(mopac_system_f), intent(in) :: system_f
    type(mopac_system), intent(out) :: system_c
    integer(c_int), pointer :: iwork(:)
    real(c_double), pointer :: rwork(:)
    integer :: status
    system_c%natom = system_f%natom
    system_c%natom_move = system_f%natom_move
    system_c%charge = system_f%charge
    system_c%spin = system_f%spin
    system_c%model = system_f%model
    system_c%epsilon = system_f%epsilon
    system_c%nlattice = system_f%nlattice
    system_c%nlattice_move = system_f%nlattice_move
    system_c%pressure = system_f%pressure
    system_c%tolerance = system_f%tolerance
    system_c%max_time = system_f%max_time
    allocate(iwork(system_f%natom), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    iwork = system_f%atom
    system_c%atom = c_loc(iwork)
    allocate(rwork(3*system_f%natom+3*system_f%nlattice), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    rwork(:3*system_f%natom) = system_f%coord
    system_c%coord = c_loc(rwork)
    if (system_f%nlattice > 0) then
      rwork(3*system_f%natom+1:) = system_f%lattice
      system_c%lattice = c_loc(rwork(3*system_f%natom+1))
    end if
  end subroutine mopac_system_f2c

  subroutine mopac_properties_c2f(system_c, properties_c, properties_f)
    type(mopac_system), intent(in) :: system_c
    type(mopac_properties), intent(in) :: properties_c
    type(mopac_properties_f), intent(out) :: properties_f
    character(kind=c_char), pointer :: cptr(:)
    type(c_ptr), pointer :: pptr(:)
    integer :: i, j, size, status
    properties_f%nerror = properties_c%nerror
    if (properties_c%nerror == 0) then
      properties_f%heat = properties_c%heat
      properties_f%dipole = properties_c%dipole
      properties_f%stress = properties_c%stress
      call copy_real(properties_f%charge, properties_c%charge, [system_c%natom])
      call copy_real(properties_f%coord_update, properties_c%coord_update, [3*system_c%natom_move])
      call copy_real(properties_f%coord_deriv, properties_c%coord_deriv, [3*system_c%natom_move])
      if (c_associated(properties_c%freq)) then
        call copy_real(properties_f%freq, properties_c%freq, [3*system_c%natom_move])
        call copy_real2(properties_f%disp, properties_c%disp, [3*system_c%natom_move,3*system_c%natom_move])
      end if
      call copy_int(properties_f%bond_index, properties_c%bond_index, [system_c%natom+1])
      do i=1, system_c%natom+1
        properties_f%bond_index(i) = properties_f%bond_index(i) + 1
      end do
      call copy_int(properties_f%bond_atom, properties_c%bond_atom, &
                    [properties_f%bond_index(system_c%natom+1)-1])
      call copy_real(properties_f%bond_order, properties_c%bond_order, &
                    [properties_f%bond_index(system_c%natom+1)-1])
      call copy_real(properties_f%lattice_update, properties_c%lattice_update, [3*system_c%nlattice_move])
      call copy_real(properties_f%lattice_deriv, properties_c%lattice_deriv, [3*system_c%nlattice_move])
    else if (properties_c%nerror > 0) then
      allocate(properties_f%error_msg(properties_c%nerror), stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
        stop 1
      end if
      call c_f_pointer(properties_c%error_msg, pptr, [properties_c%nerror])
      do i=1, properties_c%nerror
        properties_f%error_msg(i) = ' '
        call c_f_pointer(pptr(i), cptr, [120])
        do size=1, 120
          if (cptr(size) == c_null_char) exit
        end do
        size = size - 1
        do j=1, size
          properties_f%error_msg(i)(j:j) = cptr(j)
        end do
      end do
    end if
    call destroy_mopac_properties(properties_c)
  end subroutine mopac_properties_c2f

  subroutine mopac_state_f2c(state_f, state_c)
    type(mopac_state_f), intent(in) :: state_f
    type(mopac_state), intent(out) :: state_c
    real(c_double), pointer :: rptr(:)
    state_c%mpack = state_f%mpack
    if (state_f%uhf) then
      state_c%uhf = 1
    else
      state_c%uhf = 0
    end if
    call create_mopac_state(state_c)
    if (state_f%mpack /= 0) then
      call c_f_pointer(state_c%pa, rptr, [state_c%mpack])
      rptr = state_f%pa
      if (state_f%uhf) then
        call c_f_pointer(state_c%pb, rptr, [state_c%mpack])
        rptr = state_f%pb
      end if
    end if
  end subroutine mopac_state_f2c

  subroutine mopac_state_c2f(state_c, state_f)
    type(mopac_state), intent(in) :: state_c
    type(mopac_state_f), intent(out) :: state_f
    state_f%mpack = state_c%mpack
    state_f%uhf = (state_c%uhf /= 0)
    if (allocated(state_f%pa)) deallocate(state_f%pa)
    if (allocated(state_f%pb)) deallocate(state_f%pb)
    if (state_f%mpack /= 0) then
      call copy_real(state_f%pa, state_c%pa, [state_c%mpack])
      if (state_f%uhf) call copy_real(state_f%pb, state_c%pb, [state_c%mpack])
    end if
    call destroy_mopac_state(state_c)
  end subroutine mopac_state_c2f

  subroutine mozyme_state_f2c(state_f, state_c)
    type(mozyme_state_f), intent(in) :: state_f
    type(mozyme_state), intent(out) :: state_c
    integer(c_int), pointer :: iptr(:), iptr2(:,:)
    real(c_double), pointer :: rptr(:)
    state_c%numat = state_f%numat
    state_c%noccupied = state_f%noccupied
    state_c%nvirtual = state_f%nvirtual
    state_c%icocc_dim = state_f%icocc_dim
    state_c%icvir_dim = state_f%icvir_dim
    state_c%cocc_dim = state_f%cocc_dim
    state_c%cvir_dim = state_f%cvir_dim
    call create_mozyme_state(state_c)
    if (state_f%numat /= 0) then
      call c_f_pointer(state_c%nbonds, iptr, [state_c%numat])
      iptr = state_f%nbonds
      call c_f_pointer(state_c%ibonds, iptr2, [9,state_c%numat])
      iptr2 = state_f%ibonds
      call c_f_pointer(state_c%iorbs, iptr, [state_c%numat])
      iptr = state_f%iorbs
      call c_f_pointer(state_c%ncf, iptr, [state_c%noccupied])
      iptr = state_f%ncf
      call c_f_pointer(state_c%nce, iptr, [state_c%nvirtual])
      iptr = state_f%nce
      call c_f_pointer(state_c%icocc, iptr, [state_c%icocc_dim])
      iptr = state_f%icocc
      call c_f_pointer(state_c%icvir, iptr, [state_c%icvir_dim])
      iptr = state_f%icvir
      call c_f_pointer(state_c%cocc, rptr, [state_c%cocc_dim])
      rptr = state_f%cocc
      call c_f_pointer(state_c%cvir, rptr, [state_c%cvir_dim])
      rptr = state_f%cvir
    end if
  end subroutine mozyme_state_f2c

  subroutine mozyme_state_c2f(state_c, state_f)
    type(mozyme_state), intent(in) :: state_c
    type(mozyme_state_f), intent(out) :: state_f
    state_f%numat = state_c%numat
    state_f%noccupied = state_c%noccupied
    state_f%nvirtual = state_c%nvirtual
    state_f%icocc_dim = state_c%icocc_dim
    state_f%icvir_dim = state_c%icvir_dim
    state_f%cocc_dim = state_c%cocc_dim
    state_f%cvir_dim = state_c%cvir_dim
    if (allocated(state_f%nbonds)) deallocate(state_f%nbonds)
    if (allocated(state_f%ibonds)) deallocate(state_f%ibonds)
    if (allocated(state_f%iorbs)) deallocate(state_f%iorbs)
    if (allocated(state_f%ncf)) deallocate(state_f%ncf)
    if (allocated(state_f%nce)) deallocate(state_f%nce)
    if (allocated(state_f%icocc)) deallocate(state_f%icocc)
    if (allocated(state_f%icvir)) deallocate(state_f%icvir)
    if (allocated(state_f%cocc)) deallocate(state_f%cocc)
    if (allocated(state_f%cvir)) deallocate(state_f%cvir)
    if (state_f%numat /= 0) then
      call copy_int(state_f%nbonds, state_c%nbonds, [state_c%numat])
      call copy_int2(state_f%ibonds, state_c%ibonds, [9,state_c%numat])
      call copy_int(state_f%iorbs, state_c%iorbs, [state_c%numat])
      call copy_int(state_f%ncf, state_c%ncf, [state_c%noccupied])
      call copy_int(state_f%nce, state_c%nce, [state_c%nvirtual])
      call copy_int(state_f%icocc, state_c%icocc, [state_c%icocc_dim])
      call copy_int(state_f%icvir, state_c%icvir, [state_c%icvir_dim])
      call copy_real(state_f%cocc, state_c%cocc, [state_c%cocc_dim])
      call copy_real(state_f%cvir, state_c%cvir, [state_c%cvir_dim])
    end if
    call destroy_mozyme_state(state_c)
  end subroutine mozyme_state_c2f

  subroutine copy_int(data, ptr, size)
    integer, allocatable, intent(out) :: data(:)
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: size(1)
    integer(c_int), pointer :: data_c(:)
    integer :: status
    call c_f_pointer(ptr, data_c, size)
    allocate(data(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    data = data_c
  end subroutine copy_int

  subroutine copy_int2(data, ptr, size)
    integer, allocatable, intent(out) :: data(:,:)
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: size(2)
    integer(c_int), pointer :: data_c(:,:)
    integer :: status
    call c_f_pointer(ptr, data_c, size)
    allocate(data(size(1),size(2)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    data = data_c
  end subroutine copy_int2

  subroutine copy_real(data, ptr, size)
    double precision, allocatable, intent(out) :: data(:)
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: size(1)
    real(c_double), pointer :: data_c(:)
    integer :: status
    call c_f_pointer(ptr, data_c, size)
    allocate(data(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    data = data_c
  end subroutine copy_real

  subroutine copy_real2(data, ptr, size)
    double precision, allocatable, intent(out) :: data(:,:)
    type(c_ptr), intent(in) :: ptr
    integer, intent(in) :: size(2)
    real(c_double), pointer :: data_c(:,:)
    integer :: status
    call c_f_pointer(ptr, data_c, size)
    allocate(data(size(1),size(2)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API wrapper"
      stop 1
    end if
    data = data_c
  end subroutine copy_real2

end submodule mopac_api_f_internal
