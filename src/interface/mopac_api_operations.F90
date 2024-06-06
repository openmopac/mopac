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
  use Common_arrays_C, only: xparam, grad, lopt
  use molkst_C, only: keywrd, escf, nvar, moperr
  use to_screen_C, only : freq, cnorml
  implicit none

  interface

    ! initial MOPAC modules before calculations
    module subroutine mopac_initialize(system, status)
      type(mopac_system), intent(in) :: system
      integer, intent(out) :: status
    end subroutine mopac_initialize

    ! extract data from MOPAC modules after calculations
    module subroutine mopac_finalize(properties)
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_finalize

    ! save MOPAC density matrices
    module subroutine mopac_save(state, status)
      type(mopac_state), intent(out) :: state
      integer, intent(out) :: status
    end subroutine mopac_save

    ! load MOPAC density matrices, or construct initial guesses
    module subroutine mopac_load(state, status)
      type(mopac_state), intent(in) :: state
      integer, intent(out) :: status
    end subroutine mopac_load

    ! save MOZYME density matrix
    module subroutine mozyme_save(state, status)
      type(mozyme_state), intent(out) :: state
      integer, intent(out) :: status
    end subroutine mozyme_save

    ! load MOZYME density matrix, or construct initial guess
    module subroutine mozyme_load(state, status)
      type(mozyme_state), intent(in) :: state
      integer, intent(out) :: status
    end subroutine mozyme_load

  end interface

contains

  ! MOPAC electronic ground state calculation
  module subroutine mopac_scf(system, state, properties)
    !dec$ attributes dllexport :: mopac_scf
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      integer :: status

      keywrd = " 1SCF PULAY BONDS"
      call mopac_initialize(system, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_load(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      ! call computational routine for SCF calculations
      call compfg (xparam, .true., escf, .true., grad, .true.)
      ! TO DO: compfg error handling
      call mopac_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mopac_scf
  
    ! MOPAC geometry relaxation
    module subroutine mopac_relax(system, state, properties, relax_coord, relax_lattice)
    !dec$ attributes dllexport :: mopac_relax
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: relax_coord
      double precision, optional, dimension(:), intent(out) :: relax_lattice
      integer :: status

      keywrd = " LBFGS PULAY BONDS"
      call mopac_initialize(system, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_load(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      ! turn off lattice relaxation as needed
      if (system%nlattice > 0 .and. .not. present(relax_lattice)) then
        lopt(:,system%natom+1:system%natom+system%nlattice) = 0
        nvar = 3*system%natom
      end if
      ! call computational routine for LBFGS geometry optimization
      call lbfgs (xparam, escf)
      ! TO DO: lbfgs error handling
      ! save relaxed geometry
      if(size(relax_coord) < 3*system%natom) then
        ! TO DO: error handling
      end if
      relax_coord = xparam(:3*system%natom)
      if(system%nlattice > 0 .and. present(relax_lattice)) then
        if(size(relax_lattice) < 3*system%nlattice) then
          ! TO DO: error handling
        end if
        relax_lattice = xparam(3*system%natom+1:3*system%natom+3*system%nlattice)
      end if
      call mopac_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mopac_relax
  
    ! MOPAC vibrational calculation
    module subroutine mopac_vibe(system, state, properties, frequency, displacement)
    !dec$ attributes dllexport :: mopac_vibe
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: frequency
      double precision, optional, dimension(:,:), intent(out) :: displacement
      integer :: status

      keywrd = " FORCE NOREOR PULAY BONDS"
      call mopac_initialize(system, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_load(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      ! turn off lattice relaxation
!      lopt(:,system%natom+1:system%natom+system%nlattice) = 0
!      nvar = 3*system%natom
      ! call computational routine to evaluate vibrational matrix
      call force ()
      ! TO DO: force error handling
      ! save frequencies and displacement vectors
      frequency = freq
      displacement = reshape(cnorml,[nvar, nvar])
      call mopac_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mopac_vibe
  
    ! MOZYME electronic ground state calculation
    module subroutine mozyme_scf(system, state, properties)
    !dec$ attributes dllexport :: mozyme_scf
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      integer :: status

      keywrd = " MOZYME 1SCF ALLBONDS"
      call mopac_initialize(system, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mozyme_load(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      ! call computational routine for SCF calculations
      call compfg (xparam, .true., escf, .true., grad, .true.)
      ! TO DO: compfg error handling
      call mozyme_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mozyme_scf
  
    ! MOZYME geometry relaxation
    module subroutine mozyme_relax(system, state, properties, relax_coord, relax_lattice)
    !dec$ attributes dllexport :: mozyme_relax
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
      double precision, dimension(:), intent(out) :: relax_coord
      double precision, optional, dimension(:), intent(out) :: relax_lattice
      integer :: status

      keywrd = " MOZYME LBFGS ALLBONDS"
      call mopac_initialize(system, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mozyme_load(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      ! turn off lattice relaxation as needed
      if (system%nlattice > 0 .and. .not. present(relax_lattice)) then
        lopt(:,system%natom+1:system%natom+system%nlattice) = 0
        nvar = 3*system%natom
      end if
      ! call computational routine for LBFGS geometry optimization
      call lbfgs (xparam, escf)
      ! TO DO: lbfgs error handling
      ! save relaxed geometry
      if(size(relax_coord) < 3*system%natom) then
        ! TO DO: error handling
      end if
      relax_coord = xparam(:3*system%natom)
      if(system%nlattice > 0 .and. present(relax_lattice)) then
        if(size(relax_lattice) < 3*system%nlattice) then
          ! TO DO: error handling
        end if
        relax_lattice = xparam(3*system%natom+1:3*system%natom+3*system%nlattice)
      end if
      call mozyme_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mozyme_relax

end submodule mopac_api_operations