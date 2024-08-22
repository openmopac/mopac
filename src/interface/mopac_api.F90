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
  implicit none

  private
  ! public derived types of the MOPAC API
  public :: mopac_system, mopac_properties, mopac_state, mozyme_state
  ! public subroutines of the MOPAC API
  public :: mopac_scf, mopac_relax, mopac_vibe, mozyme_scf, mozyme_relax, mozyme_vibe

  ! data that defines the atomistic system and MOPAC job options
  type :: mopac_system
    ! number of atoms
    integer :: natom = 0
    ! net charge
    integer :: charge = 0
    ! number of spin excitations, floor[(number of alpha electrons)/2 - (number of beta electrons)/2]
    integer :: spin = 0
    ! dielectric constant for COSMO implicit solvent, must be 1 (no solvent) for nlattice > 0
    double precision :: epsilon = 1.d0
    ! semiempirical model: PM7 = 0, PM6-D3H4 = 1, PM6-ORG = 2, PM6 = 3, AM1 = 4, RM1 = 5
    integer :: model = 0
    ! atomic number of each atom [natom]
    integer, dimension (:), allocatable :: atom
    ! (x,y,z) coordinates of each atom (Angstroms) [3*natom]
    double precision, dimension (:), allocatable :: coord
    ! flag to determine if each atom is allowed to move [natom]
    logical, dimension (:), allocatable :: move_atom
    ! number of lattice vectors / translation vectors / periodic dimensions
    integer :: nlattice = 0
    ! external hydrostatic pressure (Gigapascals)
    double precision :: pressure = 0.d0
    ! (x,y,z) coordinates of each lattice vectors (Angstroms) [3*nlattice]
    double precision, dimension (:), allocatable :: lattice
    ! flag to determine if each lattice vector is allowed to move [nlattice]
    logical, dimension (:), allocatable :: move_lattice
    ! numerical tolerances (relative to their default values)
    double precision :: tolerance = 1.d0
    ! time limit for a MOPAC calculation (seconds)
    integer :: max_time = 3600
  end type

  ! calculated ground-state properties of an atomistic system and MOPAC job info
  type :: mopac_properties
    ! heat of formation (kcal/mol)
    double precision :: heat
    ! dipole moment vector (Debye)
    double precision, dimension (3) :: dipole
    ! atomic partial charges [natom]
    double precision, dimension (:), allocatable :: charge
    ! number of moveable atoms
    integer :: natom_move
    ! index of each moveable atom [natom_move]
    integer, dimension (:), allocatable :: atom_move
    ! (x,y,z) coordinates of each moveable atom (Angstroms) [3*natom_move]
    double precision, dimension (:), allocatable :: coord_update
    ! (x,y,z) heat gradients for each moveable atom (kcal/mol/Angstrom) [3*natom_move]
    double precision, dimension (:), allocatable :: coord_deriv
    ! availability of vibrational properties
    logical :: calc_vibe
    ! vibrational frequencies of normal modes (1/cm) [3*natom_move]
    double precision, dimension (:), allocatable :: freq
    ! displacement vectors of normal modes [3*natom_move,3*natom_move]
    double precision, dimension (:,:), allocatable :: disp
    ! bond-order matrix in compressed sparse column (CSC) matrix format
    ! with insignificant bond orders (<0.001) truncated
    ! diagonal matrix entries are atomic valencies
    ! > first index of each atom in CSC bond-order matrix [natom+1]
    integer, dimension (:), allocatable :: bond_index
    ! > list of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
    integer, dimension (:), allocatable :: bond_atom
    ! > bond order of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
    double precision, dimension (:), allocatable :: bond_order
    ! number of moveable lattice vectors
    integer :: nlattice_move
    ! index of each moveable lattice vector [nlattice_move]
    integer, dimension (:), allocatable :: lattice_move
    ! (x,y,z) coordinates of each moveable lattice vectors (Angstroms) [3*nlattice_move]
    double precision, dimension (:), allocatable :: lattice_update
    ! (x,y,z) heat gradients for each moveable lattice vector (kcal/mol/Angstrom) [3*nlattice_move]
    double precision, dimension (:), allocatable :: lattice_deriv
    ! stress tensor (Gigapascals) in Voigt form (xx, yy, zz, yz, xz, xy) for nlattice_move == 3
    double precision, dimension (6) :: stress
    ! status of MOPAC job
    integer :: status
    ! TO DO: compile list of status values & their meaning
  end type

  ! data that describes the electronic state using standard molecular orbitals
  type :: mopac_state
    ! availability of a saved state
    logical :: save_state = .false.
    ! MOPAC data format is copied from molkst_C and Common_arrays_C modules
    ! > number of matrix elements in packed lower triangle matrix format
    integer :: mpack
    ! > alpha density matrix [mpack]
    double precision, dimension (:), allocatable :: pa
    ! > beta density matrix [mpack]
    double precision, dimension (:), allocatable :: pb
  end type

  ! data that describes the electronic state using localized molecular orbitals
  type :: mozyme_state
    ! availability of a saved state
    logical :: save_state = .false.
    ! MOZYME data format is copied from molkst_C, Common_arrays_C, and MOZYME_C modules
    ! > number of real atoms
    integer :: numat
    ! > number of Lewis bonds per real atom [numat]
    integer, dimension (:), allocatable :: nbonds
    ! > list of Lewis-bonded real atoms for each real atom [9,numat]
    integer, dimension (:,:), allocatable :: ibonds
    ! > number of orbitals per real atom [numat]
    integer, dimension (:), allocatable :: iorbs
    ! > number of occupied molecular orbitals
    integer :: noccupied
    ! > number of atoms in each occupied LMO [noccupied]
    integer, dimension (:), allocatable :: ncf
    ! > number of virtual molecular orbitals
    integer :: nvirtual
    ! > number of atoms in each virtual LMO [nvirtual]
    integer, dimension (:), allocatable :: nce
    ! > size of array icocc
    integer :: icocc_dim
    ! > index of each real atom in the occupied LMOs [iccoc_dim]
    integer, dimension (:), allocatable :: icocc
    ! > size of array icvir
    integer :: icvir_dim
    ! > index of each real atom in the virtual LMOs [icvir_dim]
    integer, dimension (:), allocatable :: icvir
    ! > size of array cocc
    integer :: cocc_dim
    ! > atomic orbital coefficients of the occupied LMOs [cocc_dim]
    integer, dimension (:), allocatable :: cocc
    ! > size of array cvir
    integer :: cvir_dim
    ! > atomic orbital coefficients of the virtual LMOs [cvir_dim]
    integer, dimension (:), allocatable :: cvir
  end type

  interface

    ! MOPAC electronic ground state calculation
    module subroutine mopac_scf(system, state, properties)
    !dec$ attributes dllexport :: mopac_scf
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_scf

    ! MOPAC geometry relaxation
    module subroutine mopac_relax(system, state, properties)
    !dec$ attributes dllexport :: mopac_relax
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_relax

    ! MOPAC vibrational calculation
    module subroutine mopac_vibe(system, state, properties)
    !dec$ attributes dllexport :: mopac_vibe
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_vibe

    ! MOZYME electronic ground state calculation
    module subroutine mozyme_scf(system, state, properties)
    !dec$ attributes dllexport :: mozyme_scf
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_scf

    ! MOZYME geometry relaxation
    module subroutine mozyme_relax(system, state, properties)
    !dec$ attributes dllexport :: mozyme_relax
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
    end subroutine mozyme_relax

    ! MOZYME vibrational calculation
    module subroutine mozyme_vibe(system, state, properties)
      !dec$ attributes dllexport :: mopac_vibe
        type(mopac_system), intent(in) :: system
        type(mozyme_state), intent(inout) :: state
        type(mopac_properties), intent(out) :: properties
      end subroutine mozyme_vibe
  end interface

end module mopac_api
