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

submodule (mopac_api:mopac_api_operations) mopac_api_saveload
  use Common_arrays_C, only: pa, pb, nbonds, ibonds
  use molkst_C, only: keywrd, uhf, mpack, numat
  use MOZYME_C, only: iorbs, noccupied, ncf, nvirtual, nce, icocc_dim, &
    icocc, icvir_dim, icvir, cocc_dim, cocc, cvir_dim, cvir
  implicit none

contains

  ! save MOPAC density matrices
  module subroutine mopac_save(state, status)
    type(mopac_state), intent(out) :: state
    integer, intent(out) :: status

    state%save_state = .true.
    state%mpack = mpack

    if (allocated(state%pa)) deallocate(state%pa)
    allocate(state%pa(mpack), stat=status)
    if (status /= 0) return
    state%pa = pa
    if (uhf) then
      if (allocated(state%pb)) deallocate(state%pb)
      allocate(state%pb(mpack), stat=status)
      if (status /= 0) return
      state%pb = pb
    end if
  end subroutine mopac_save

  ! load MOPAC density matrices, or construct initial guesses
  module subroutine mopac_load(state, status)
    type(mopac_state), intent(in) :: state
    integer, intent(out) :: status

    if(state%save_state) then
      ! TO DO: compatibility tests
      keywrd = trim(keywrd) // " OLDENS"
      mpack = state%mpack
      if (allocated(pa)) deallocate(pa)
      allocate(pa(mpack), stat=status)
      if (status /= 0) return
      pa = state%pa
      if(uhf) then
        if (allocated(pb)) deallocate(pb)
        allocate(pb(mpack), stat=status)
        if (status /= 0) return
        pb = state%pb
      end if
    end if
  end subroutine mopac_load

  ! save MOZYME density matrix
  module subroutine mozyme_save(state, status)
    type(mozyme_state), intent(out) :: state
    integer, intent(out) :: status

    state%save_state = .true.
    state%numat = numat
    state%noccupied = noccupied
    state%nvirtual = nvirtual
    state%icocc_dim = icocc_dim
    state%icvir_dim = icvir_dim
    state%cocc_dim = cocc_dim
    state%cvir_dim = cvir_dim
    if (allocated(state%nbonds)) deallocate(state%nbonds)
    if (allocated(state%ibonds)) deallocate(state%ibonds)
    if (allocated(state%iorbs)) deallocate(state%iorbs)
    if (allocated(state%ncf)) deallocate(state%ncf)
    if (allocated(state%nce)) deallocate(state%nce)
    if (allocated(state%icocc)) deallocate(state%icocc)
    if (allocated(state%icvir)) deallocate(state%icvir)
    if (allocated(state%cocc)) deallocate(state%cocc)
    if (allocated(state%cvir)) deallocate(state%cvir)
    allocate(state%nbonds(numat), stat=status)
    if (status /= 0) return
    allocate(state%ibonds(9,numat), stat=status)
    if (status /= 0) return
    allocate(state%iorbs(numat), stat=status)
    if (status /= 0) return
    allocate(state%ncf(noccupied), stat=status)
    if (status /= 0) return
    allocate(state%nce(nvirtual), stat=status)
    if (status /= 0) return
    allocate(state%icocc(icocc_dim), stat=status)
    if (status /= 0) return
    allocate(state%icvir(icvir_dim), stat=status)
    if (status /= 0) return
    allocate(state%cocc(cocc_dim), stat=status)
    if (status /= 0) return
    allocate(state%cvir(cvir_dim), stat=status)
    if (status /= 0) return
    state%nbonds = nbonds
    state%ibonds = ibonds
    state%iorbs = iorbs
    state%ncf = ncf
    state%nce = nce
    state%icocc = icocc
    state%icvir = icvir
    state%cocc = cocc
    state%cvir = cvir
  end subroutine mozyme_save

  ! load MOZYME density matrix, or construct initial guess
  module subroutine mozyme_load(state, status)
    type(mozyme_state), intent(in) :: state
    integer, intent(out) :: status

    if(state%save_state) then
      ! TO DO: compatibility tests
      keywrd = trim(keywrd) // " OLDENS"
      numat = state%numat
      noccupied = state%noccupied
      nvirtual = state%nvirtual
      icocc_dim = state%icocc_dim
      icvir_dim = state%icvir_dim
      cocc_dim = state%cocc_dim
      cvir_dim = state%cvir_dim
      if (allocated(nbonds)) deallocate(nbonds)
      if (allocated(ibonds)) deallocate(ibonds)
      if (allocated(iorbs)) deallocate(iorbs)
      if (allocated(ncf)) deallocate(ncf)
      if (allocated(nce)) deallocate(nce)
      if (allocated(icocc)) deallocate(icocc)
      if (allocated(icvir)) deallocate(icvir)
      if (allocated(cocc)) deallocate(cocc)
      if (allocated(cvir)) deallocate(cvir)
      allocate(nbonds(numat), stat=status)
      if (status /= 0) return
      allocate(ibonds(9,numat), stat=status)
      if (status /= 0) return
      allocate(iorbs(numat), stat=status)
      if (status /= 0) return
      allocate(ncf(noccupied), stat=status)
      if (status /= 0) return
      allocate(nce(nvirtual), stat=status)
      if (status /= 0) return
      allocate(icocc(icocc_dim), stat=status)
      if (status /= 0) return
      allocate(icvir(icvir_dim), stat=status)
      if (status /= 0) return
      allocate(cocc(cocc_dim), stat=status)
      if (status /= 0) return
      allocate(cvir(cvir_dim), stat=status)
      if (status /= 0) return
      nbonds = state%nbonds
      ibonds = state%ibonds
      iorbs = state%iorbs
      ncf = state%ncf
      nce = state%nce
      icocc = state%icocc
      icvir = state%icvir
      cocc = state%cocc
      cvir = state%cvir
    else
      call geochk()
    end if
  end subroutine mozyme_load

end submodule mopac_api_saveload