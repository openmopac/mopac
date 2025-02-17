! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

submodule (mopac_api) mopac_api_operations
  use Common_arrays_C, only: xparam, grad, lopt
  use molkst_C, only: keywrd, escf, moperr, nvar, jobnam, run, verson
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
      type(mopac_state), intent(inout) :: state
    end subroutine mopac_save

    ! load MOPAC density matrices, or construct initial guesses
    module subroutine mopac_load(state)
      type(mopac_state), intent(in) :: state
    end subroutine mopac_load

    ! save MOZYME density matrix
    module subroutine mozyme_save(state)
      type(mozyme_state), intent(inout) :: state
    end subroutine mozyme_save

    ! load MOZYME density matrix, or construct initial guess
    module subroutine mozyme_load(state)
      type(mozyme_state), intent(in) :: state
    end subroutine mozyme_load

  end interface

contains

  ! MOPAC electronic ground state calculation
  module subroutine mopac_scf(system, state, properties) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: mopac_scf
#endif
    type(mopac_system), intent(in) :: system
    type(mopac_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " 1SCF PULAY BONDS"
    if (system%natom_move + system%nlattice_move > 0) trim(keywrd) // " GRADIENTS"
    call mopac_initialize(system)
    if (.not. moperr) call mopac_load(state)
    ! call computational routine for SCF calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mopac_save(state)
    call mopac_finalize(properties)
  end subroutine mopac_scf

  ! MOPAC geometry relaxation
  module subroutine mopac_relax(system, state, properties) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: mopac_relax
#endif
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
#ifdef WIN32
!dec$ attributes dllexport :: mopac_vibe
#endif
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
#ifdef WIN32
!dec$ attributes dllexport :: mozyme_scf
#endif
    type(mopac_system), intent(in) :: system
    type(mozyme_state), intent(inout) :: state
    type(mopac_properties), intent(out) :: properties

    keywrd = " MOZYME 1SCF ALLBONDS"
    if (system%natom_move + system%nlattice_move > 0) trim(keywrd) // " GRADIENTS"
    call mopac_initialize(system)
    if (.not. moperr) call mozyme_load(state)
    ! call computational routine for SCF calculations
    if (.not. moperr) call compfg (xparam, .true., escf, .true., grad, .true.)
    if (.not. moperr) call mozyme_save(state)
    call mopac_finalize(properties)
  end subroutine mozyme_scf

  ! MOZYME geometry relaxation
  module subroutine mozyme_relax(system, state, properties) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: mozyme_relax
#endif
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
#ifdef WIN32
!dec$ attributes dllexport :: mozyme_vibe
#endif
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
  module function run_mopac_from_input(path_to_file) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: run_mopac_from_input
#endif
    interface
      subroutine run_mopac() bind(c)
      end subroutine run_mopac
    end interface
    integer(c_int) :: run_mopac_from_input
    character(kind=c_char), dimension(*), intent(in) :: path_to_file
    integer :: i
    i = 1
    do
      if(path_to_file(i) == ' ' .or. path_to_file(i) == c_null_char) exit
      jobnam(i:i) = path_to_file(i)
      i = i + 1
    end do
    run = 2
    call run_mopac
    if (moperr) then
      run_mopac_from_input = 1
    else
      run_mopac_from_input = 0
    end if
    run = 1
    jobnam = ' '
  end function run_mopac_from_input

  ! Get MOPAC version string
  module subroutine get_mopac_version(version) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: get_mopac_version
#endif
    character(kind=c_char), dimension(*), intent(out) :: version
    integer :: i
    i = 1
    do
      if(verson(i:i) == ' ') exit
      version(i) = verson(i:i)
      i = i + 1
    end do
    version(i) = c_null_char
  end subroutine get_mopac_version

end submodule mopac_api_operations
