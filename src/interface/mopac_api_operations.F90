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
  use Common_arrays_C, only: xparam, grad
  use molkst_C, only: keywrd, escf, moperr
  implicit none

  interface

    ! initialize MOPAC modules before calculations
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
    module subroutine mopac_relax(system, state, properties)
    !dec$ attributes dllexport :: mopac_relax
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
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
      ! call computational routine for LBFGS geometry optimization
      call lbfgs (xparam, escf)
      ! TO DO: lbfgs error handling
      call mopac_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mopac_relax
  
    ! MOPAC vibrational calculation
    module subroutine mopac_vibe(system, state, properties)
    !dec$ attributes dllexport :: mopac_vibe
      type(mopac_system), intent(in) :: system
      type(mopac_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
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
      ! call computational routine to evaluate vibrational matrix
      call force ()
      ! TO DO: force error handling
write(*,*) "done"
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

      keywrd = " MOZYME 1SCF GRADIENTS ALLBONDS"
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
    module subroutine mozyme_relax(system, state, properties)
    !dec$ attributes dllexport :: mozyme_relax
      type(mopac_system), intent(in) :: system
      type(mozyme_state), intent(inout) :: state
      type(mopac_properties), intent(out) :: properties
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
      ! call computational routine for LBFGS geometry optimization
      call lbfgs (xparam, escf)
      ! TO DO: lbfgs error handling
      call mozyme_save(state, status)
      if (status /= 0) then
        properties%status = status
        return
      end if
      call mopac_finalize(properties)
    end subroutine mozyme_relax

    ! MOPAC vibrational calculation
    module subroutine mozyme_vibe(system, state, properties)
      !dec$ attributes dllexport :: mopac_vibe
        type(mopac_system), intent(in) :: system
        type(mozyme_state), intent(inout) :: state
        type(mopac_properties), intent(out) :: properties
        integer :: status
  
        keywrd = " FORCE NOREOR PULAY BONDS"
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
        ! call computational routine to evaluate vibrational matrix
        call force ()
        ! TO DO: force error handling
        call mozyme_save(state, status)
        if (status /= 0) then
          properties%status = status
          return
        end if
        call mopac_finalize(properties)
      end subroutine mozyme_vibe

end submodule mopac_api_operations
