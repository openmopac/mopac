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

submodule (mopac_api) mopac_api_operations
  implicit none

contains

  ! MOPAC electronic ground state calculation
  subroutine mopac_scf(system, state, properties)
    !dec$ attributes dllexport :: mopac_scf
      use molkst_C, only: keywrd
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
  
      keywrd = "1SCF"
      mopac_initialize(system, status)
      ! TO DO: parse status
      mopac_load(state, status)
      ! TO DO: parse status
      ! TO DO: prepare arguments for compfg
      call compfg (xparam, .true., escf, .true., grad, .true.)
      mopac_save(state, status)
      ! TO DO: parse status
      mopac_finalize(properties)
    end subroutine mopac_scf
  
    ! MOPAC geometry relaxation
    subroutine mopac_relax(system, state, properties, relax_coord, relax_lattice)
    !dec$ attributes dllexport :: mopac_relax
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: relax_coord
      double precision, optional, dimension(:), intent(out) :: relax_lattice
  
      keywrd = "LBFGS"
      mopac_initialize(system, status)
      ! TO DO: parse status
      mopac_load(state, status)
      ! TO DO: parse status
      ! TO DO: prepare arguments for lbfgs
      call lbfgs (xparam, escf)
      ! TO DO: extract relax_coord & relax_lattice
      mopac_save(state, status)
      ! TO DO: parse status
      mopac_finalize(properties)
    end subroutine mopac_relax
  
    ! MOPAC vibrational calculation
    subroutine mopac_vibe(system, state, properties, frequency, displacement)
    !dec$ attributes dllexport :: mopac_vibe
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: frequency
      double precision, optional, dimension(:,:), intent(out) :: displacement
  
      keywrd = "FORCE"
      mopac_initialize(system, status)
      ! TO DO: parse status
      mopac_load(state, status)
      ! TO DO: parse status
      call force ()
      ! TO DO: extract frequency & displacement
      mopac_save(state, status)
      ! TO DO: parse status
      mopac_finalize(properties)
    end subroutine mopac_vibe
  
    ! MOZYME electronic ground state calculation
    subroutine mozyme_scf(system, state, properties)
    !dec$ attributes dllexport :: mozyme_scf
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
  
      keywrd = "MOZYME 1SCF"
      mopac_initialize(system, status)
      ! TO DO: parse status
      mozyme_load(state, status)
      ! TO DO: parse status
      ! TO DO: prepare arguments for compfg
      call compfg (xparam, .true., escf, .true., grad, .true.)
      mozyme_save(state, status)
      ! TO DO: parse status
      mopac_finalize(properties)
    end subroutine mozyme_scf
  
    ! MOZYME geometry relaxation
    subroutine mozyme_relax(system, state, properties, relax_coord, relax_lattice)
    !dec$ attributes dllexport :: mozyme_relax
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: relax_coord
      double precision, optional, dimension(:), intent(out) :: relax_lattice
  
      keywrd = "MOZYME LBFGS"
      mopac_initialize(system, status)
      ! TO DO: parse status
      mozyme_load(state, status)
      ! TO DO: parse status
      ! TO DO: prepare arguments for lbfgs
      call lbfgs (xparam, escf)
      ! TO DO: extract relax_coord & relax_lattice
      mopac_save(state, status)
      ! TO DO: parse status
      mozyme_finalize(properties)
    end subroutine mozyme_relax

end submodule mopac_api_operations