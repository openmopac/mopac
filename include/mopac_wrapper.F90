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

! Fortran wrapper for diskless/stateless Application Programming Interface (API) to core MOPAC operations
module mopac_api_f
  implicit none

  private
  ! public derived types of the MOPAC API
  public :: mopac_system_f, mopac_properties_f, mopac_state_f, mozyme_state_f
  ! public subroutines of the MOPAC API
  public :: mopac_scf_f, mopac_relax_f, mopac_vibe_f, mozyme_scf_f, mozyme_relax_f, mozyme_vibe_f
  ! public subroutine of the simple, disk-based MOPAC API
  public :: run_mopac_from_input_f

  ! data that defines the atomistic system and MOPAC job options (Fortran wrapper)
  type :: mopac_system_f
    ! number of atoms
    integer :: natom = 0
    ! number of atoms that are allowed to move (first natom_move atoms in array)
    integer :: natom_move = 0
    ! net charge
    integer :: charge = 0
    ! number of spin excitations, floor[(number of alpha electrons)/2 - (number of beta electrons)/2]
    integer :: spin = 0
    ! semiempirical model: PM7 = 0, PM6-D3H4 = 1, PM6-ORG = 2, PM6 = 3, AM1 = 4, RM1 = 5
    integer :: model = 0
    ! dielectric constant for COSMO implicit solvent, must be 1 (no solvent) for nlattice > 0
    double precision :: epsilon = 1.d0
    ! atomic number of each atom [natom]
    integer, dimension (:), allocatable :: atom
    ! (x,y,z) coordinates of each atom (Angstroms) [3*natom]
    double precision, dimension (:), allocatable :: coord
    ! number of lattice vectors / translation vectors / periodic dimensions
    integer :: nlattice = 0
    ! number of lattice vectors that are allowed to move (first nlattice_move vectors in array)
    integer :: nlattice_move = 0
    ! external hydrostatic pressure (Gigapascals)
    double precision :: pressure = 0.d0
    ! (x,y,z) coordinates of each lattice vectors (Angstroms) [3*nlattice]
    double precision, dimension (:), allocatable :: lattice
    ! relative numerical tolerances (equivalent to GNORM and RELSCF keyword values)
    double precision :: tolerance = 1.d0
    ! time limit for a MOPAC calculation (seconds)
    integer :: max_time = 3600
  end type

  ! calculated ground-state properties of an atomistic system and MOPAC job info (Fortran wrapper)
  type :: mopac_properties_f
    ! heat of formation (kcal/mol)
    double precision :: heat
    ! dipole moment vector (Debye)
    double precision, dimension (3) :: dipole
    ! atomic partial charges [natom]
    double precision, dimension (:), allocatable :: charge
    ! (x,y,z) coordinates of each moveable atom (Angstroms) [3*natom_move]
    double precision, dimension (:), allocatable :: coord_update
    ! (x,y,z) heat gradients for each moveable atom (kcal/mol/Angstrom) [3*natom_move]
    double precision, dimension (:), allocatable :: coord_deriv
    ! vibrational frequencies of normal modes (1/cm) [3*natom_move], if available
    double precision, dimension (:), allocatable :: freq
    ! (x,y,z) displacement vectors of normal modes [3*natom_move,3*natom_move], if available
    double precision, dimension (:,:), allocatable :: disp
    ! bond-order matrix in compressed sparse column (CSC) matrix format
    ! with insignificant bond orders (<0.01) truncated
    ! diagonal matrix entries are atomic valencies
    ! > first index of each atom in CSC bond-order matrix [natom+1]
    integer, dimension (:), allocatable :: bond_index
    ! > list of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
    integer, dimension (:), allocatable :: bond_atom
    ! > bond order of atoms bonded to each atom in CSC format [bond_index(natom+1)-1]
    double precision, dimension (:), allocatable :: bond_order
    ! (x,y,z) coordinates of each moveable lattice vectors (Angstroms) [3*nlattice_move]
    double precision, dimension (:), allocatable :: lattice_update
    ! (x,y,z) heat gradients for each moveable lattice vector (kcal/mol/Angstrom) [3*nlattice_move]
    double precision, dimension (:), allocatable :: lattice_deriv
    ! stress tensor (Gigapascals) in Voigt form (xx, yy, zz, yz, xz, xy), 0 if available
    double precision, dimension (6) :: stress
    ! number of MOPAC error messages (negative value indicates that allocation of error_msg failed)
    integer :: nerror
    ! text of MOPAC error messages [nerror,120]
    character(len=120), dimension (:), allocatable :: error_msg
  end type

  ! data that describes the electronic state using standard molecular orbitals (Fortran wrapper)
  type :: mopac_state_f
    ! MOPAC data format is adapted from molkst_C and Common_arrays_C modules
    ! > number of matrix elements in packed lower triangle matrix format
    integer :: mpack = 0 ! 0 if state is unavailable
    ! > flag for unrestricted Hartree-Fock ground state (0 == restricted, 1 == unrestricted)
    logical :: uhf
    ! > alpha density matrix [mpack]
    double precision, dimension (:), allocatable :: pa
    ! > beta density matrix [mpack], empty if uhf == 0
    double precision, dimension (:), allocatable :: pb
  end type

  ! data that describes the electronic state using localized molecular orbitals (Fortran wrapper)
  type :: mozyme_state_f
    ! MOZYME data format is adapted from molkst_C, Common_arrays_C, and MOZYME_C modules
    ! > number of real atoms
    integer :: numat = 0 ! 0 if state is unavailable
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
    double precision, dimension (:), allocatable :: cocc
    ! > size of array cvir
    integer :: cvir_dim
    ! > atomic orbital coefficients of the virtual LMOs [cvir_dim]
    double precision, dimension (:), allocatable :: cvir
  end type

  ! library subroutines
  interface

    ! MOPAC electronic ground state calculation
    module subroutine mopac_scf_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mopac_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mopac_scf_f

    ! MOPAC geometry relaxation
    module subroutine mopac_relax_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mopac_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mopac_relax_f

    ! MOPAC vibrational calculation
    module subroutine mopac_vibe_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mopac_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mopac_vibe_f

    ! MOZYME electronic ground state calculation
    module subroutine mozyme_scf_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mozyme_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mozyme_scf_f

    ! MOZYME geometry relaxation
    module subroutine mozyme_relax_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mozyme_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mozyme_relax_f

    ! MOZYME vibrational calculation
    module subroutine mozyme_vibe_f(system, state, properties)
      type(mopac_system_f), intent(in) :: system
      type(mozyme_state_f), intent(inout) :: state
      type(mopac_properties_f), intent(out) :: properties
    end subroutine mozyme_vibe_f

    ! run MOPAC conventionally from an input file
    module function run_mopac_from_input_f(path_to_file)
      logical :: run_mopac_from_input_f
      character(len=240), intent(in) :: path_to_file
    end function run_mopac_from_input_f

  end interface

end module mopac_api_f
