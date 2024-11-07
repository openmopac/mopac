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
#ifdef WIN32
!dec$ attributes dllexport :: mopac_scf, mopac_relax, mopac_vibe, mozyme_scf, mozyme_relax, mozyme_vibe
#endif
  use Common_arrays_C, only: xparam, grad, lopt
  use molkst_C, only: keywrd, escf, moperr, nvar, gui, jobnam, run
  implicit none

  interface

    ! initialize MOPAC modules before calculations
    module subroutine mopac_initialize(system)
      type(mopac_system), intent(in) :: system
    end subroutine mopac_initialize

    ! extract data from MOPAC modules after calculations
    module subroutine mopac_finalize(properties)
      type(mopac_properties), intent(out) :: properties
    end subroutine mopac_finalize

    ! save MOPAC density matrices
    module subroutine mopac_save(state)
      type(mopac_state), intent(out) :: state
    end subroutine mopac_save

    ! load MOPAC density matrices, or construct initial guesses
    module subroutine mopac_load(state)
      type(mopac_state), intent(in) :: state
    end subroutine mopac_load

    ! save MOZYME density matrix
    module subroutine mozyme_save(state)
      type(mozyme_state), intent(out) :: state
    end subroutine mozyme_save

    ! load MOZYME density matrix, or construct initial guess
    module subroutine mozyme_load(state)
      type(mozyme_state), intent(in) :: state
    end subroutine mozyme_load

  end interface

contains

  ! MOPAC electronic ground state calculation
  module subroutine mopac_scf(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mopac_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " 1SCF PULAY GRADIENTS BONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mopac_load(state)
    ! call computational routine for SCF calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mopac_save(state)
    call mopac_finalize(properties)
  end subroutine mopac_scf

  ! MOPAC geometry relaxation
  module subroutine mopac_relax(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mopac_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " LBFGS PULAY BONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mopac_load(state)
    ! call computational routine for LBFGS geometry optimization
    if (.not. moperr) call lbfgs (xparam, escf)
    if (.not. moperr) call mopac_save(state)
    call mopac_finalize(properties)
  end subroutine mopac_relax

  ! MOPAC vibrational calculation
  module subroutine mopac_vibe(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mopac_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties
    integer :: i

    keywrd = " FORCETS PULAY BONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mopac_load(state)
    ! call computational routine to evaluate vibrational matrix
    if (.not. moperr) call force()
    ! recompute nvar because it is zero'd after vibrational calculations
    if (.not. moperr) then
      nvar = 0
      do i=1, system%natom+system%nlattice
        if (lopt(1,i) == 1) nvar = nvar + 3
      end do
    end if
    ! call computational routine for standard properties calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mopac_save(state)
    call mopac_finalize(properties)
  end subroutine mopac_vibe

  ! MOZYME electronic ground state calculation
  module subroutine mozyme_scf(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mozyme_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " MOZYME 1SCF GRADIENTS ALLBONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mozyme_load(state)
    ! call computational routine for SCF calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mozyme_save(state)
    call mopac_finalize(properties)
  end subroutine mozyme_scf

  ! MOZYME geometry relaxation
  module subroutine mozyme_relax(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mozyme_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " MOZYME LBFGS ALLBONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mozyme_load(state)
    ! call computational routine for LBFGS geometry optimization
    if (.not. moperr) call lbfgs (xparam, escf)
    if (.not. moperr) call mozyme_save(state)
    call mopac_finalize(properties)
  end subroutine mozyme_relax

  ! MOZYME vibrational calculation
  module subroutine mozyme_vibe(system, state, properties) bind(c)
    type(mopac_system), intent(in) :: system
    type(mozyme_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties
    integer :: i

    keywrd = " MOZYME FORCETS PULAY BONDS"
    call mopac_initialize(system)
    if (.not. moperr) call mozyme_load(state)
    ! call computational routine to evaluate vibrational matrix
    if (.not. moperr) call force()
    ! recompute nvar because it is zero'd after vibrational calculations
    if (.not. moperr) then
      nvar = 0
      do i=1, system%natom+system%nlattice
        if (lopt(1,i) == 1) nvar = nvar + 3
      end do
    end if
    ! call computational routine for standard properties calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mozyme_save(state)
    call mopac_finalize(properties)
  end subroutine mozyme_vibe

  ! Run MOPAC conventionally from an input file
  module subroutine run_mopac_from_input(path_to_file) bind(c)
    character(kind=c_char), dimension(*), intent(in) :: path_to_file
    integer :: i
    i = 1
    do
      if(path_to_file(i) == ' ' .or. path_to_file(i) == c_null_char) exit
      jobnam(i:i) = path_to_file(i)
      i = i + 1
    end do
    gui = .false.
    run = 2
    call run_mopac
    run = 1
    gui = .true.
    jobnam = ' '
  end subroutine run_mopac_from_input
  
end submodule mopac_api_operations
