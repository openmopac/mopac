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
module mopac_api
  use iso_c_binding
  implicit none

  private
  ! public derived types of the MOPAC API
  public :: mopac_system, mopac_properties, mopac_state, mozyme_state
  ! public subroutines to create or destroy a MOPAC system
  public :: create_mopac_system, destroy_mopac_system
  ! public subroutines to destroy the other derived types of the API
  public :: destroy_mopac_properties, destroy_mopac_state, destroy_mozyme_state
  ! public subroutines of the MOPAC API
  public :: mopac_scf, mopac_relax, mopac_vibe, mozyme_scf, mozyme_relax, mozyme_vibe
  ! public subroutine of the simple, disk-based MOPAC API
  public :: run_mopac_from_input

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
    ! bond-order matrix in compressed sparse column (CSC) matrix format
    ! with insignificant bond orders (<0.01) truncated
    ! diagonal matrix entries are atomic valencies
    ! > first index of each atom in CSC bond-order matrix
    type(c_ptr) :: bond_index ! integer(c_int)[natom+1]
    ! > list of atoms bonded to each atom in CSC format
    type(c_ptr) :: bond_atom ! integer(c_int)[bond_index(natom+1)-1]
    ! > bond order of atoms bonded to each atom in CSC format
    type(c_ptr) :: bond_order ! real(c_double)[bond_index(natom+1)-1]
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
    ! > alpha density matrix
    type(c_ptr) :: pa ! real(c_double)[mpack]
    ! > beta density matrix
    type(c_ptr) :: pb ! real(c_double)[mpack]
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

  interface

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
    end subroutine create_mopac_system

    ! deallocate memory in mopac_system
    module subroutine destroy_mopac_system(system) bind(c)
      type(mopac_system), intent(in) :: system
    end subroutine destroy_mopac_system

    ! deallocate memory in mopac_properties
    module subroutine destroy_mopac_properties(properties) bind(c)
      type(mopac_properties), intent(in) :: properties
    end subroutine destroy_mopac_properties

    ! deallocate memory in mopac_state
    module subroutine destroy_mopac_state(state) bind(c)
      type(mopac_state), intent(in) :: state
    end subroutine destroy_mopac_state

    ! deallocate memory in mozyme_state
    module subroutine destroy_mozyme_state(state) bind(c)
      type(mozyme_state), intent(in) :: state
    end subroutine destroy_mozyme_state

    ! MOPAC electronic ground state calculation
    module subroutine mopac_scf(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_scf

    ! MOPAC geometry relaxation
    module subroutine mopac_relax(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_relax

    ! MOPAC vibrational calculation
    module subroutine mopac_vibe(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_vibe

    ! MOZYME electronic ground state calculation
    module subroutine mozyme_scf(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_scf

    ! MOZYME geometry relaxation
    module subroutine mozyme_relax(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_relax

    ! MOZYME vibrational calculation
    module subroutine mozyme_vibe(system, state, properties) bind(c)
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_vibe

    ! run MOPAC conventionally from an input file
    module subroutine run_mopac_from_input(path_to_file) bind(c)
      character(kind=c_char,len=*), intent(in) :: path_to_file
    end subroutine run_mopac_from_input
    
  end interface

end module mopac_api
