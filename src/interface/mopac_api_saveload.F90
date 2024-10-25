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
    icocc, icvir_dim, icvir, cocc_dim, cocc, cvir_dim, cvir, nnce, nncf, ncocc, ncvir
  implicit none

contains

  ! save MOPAC density matrices
  module subroutine mopac_save(state)
    type(mopac_state), intent(out) :: state
    integer :: status
    real(c_double), pointer :: rptr(:)

    if (state%mpack > 0) then
      call c_f_pointer(state%pa, rptr, [state%mpack])
      deallocate(rptr)
      if (uhf) then
        call c_f_pointer(state%pb, rptr, [state%mpack])
        deallocate(rptr)
      end if
    end if

    state%mpack = mpack
    allocate(rptr(mpack), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_SAVE")
      return
    end if
    rptr = pa
    state%pa = c_loc(rptr)
    if (uhf) then
      allocate(rptr(mpack), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_SAVE")
        return
      end if
      rptr = pb
      state%pb = c_loc(rptr)
    end if
  end subroutine mopac_save

  ! load MOPAC density matrices, or construct initial guesses
  module subroutine mopac_load(state)
    type(mopac_state), intent(in) :: state
    integer :: status
    real(c_double), pointer :: rptr(:)

    if(state%mpack > 0) then
      if (state%mpack /= mpack) then
        call mopend("Attempting to load incompatible MOPAC state")
        return
      end if
      keywrd = trim(keywrd) // " OLDENS"
      mpack = state%mpack
      if (allocated(pa)) deallocate(pa)
      allocate(pa(mpack), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_LOAD")
        return
      end if
      call c_f_pointer(state%pa, rptr, [mpack])
      pa = rptr
      if(uhf) then
        if (allocated(pb)) deallocate(pb)
        allocate(pb(mpack), stat=status)
        if (status /= 0) then
          call mopend("Failed to allocate memory in MOPAC_LOAD")
          return
        end if
        call c_f_pointer(state%pb, rptr, [mpack])
        pb = rptr
        end if
    end if
  end subroutine mopac_load

  ! save MOZYME density matrix
  module subroutine mozyme_save(state)
    type(mozyme_state), intent(out) :: state
    integer :: status
    integer(c_int), pointer :: iptr(:), iptr2(:,:)
    real(c_double), pointer :: rptr(:)

    if (state%numat > 0) then
      call c_f_pointer(state%nbonds, iptr, [state%numat])
      deallocate(iptr)
      call c_f_pointer(state%ibonds, iptr2, [9,state%numat])
      deallocate(iptr2)
      call c_f_pointer(state%iorbs, iptr, [state%numat])
      deallocate(iptr)
      call c_f_pointer(state%ncf, iptr, [state%noccupied])
      deallocate(iptr)
      call c_f_pointer(state%nce, iptr, [state%nvirtual])
      deallocate(iptr)
      call c_f_pointer(state%icocc, iptr, [state%icocc_dim])
      deallocate(iptr)
      call c_f_pointer(state%icvir, iptr, [state%icvir_dim])
      deallocate(iptr)
      call c_f_pointer(state%cocc, rptr, [state%cocc_dim])
      deallocate(rptr)
      call c_f_pointer(state%cvir, rptr, [state%cvir_dim])
      deallocate(rptr)
    end if

    state%numat = numat
    state%noccupied = noccupied
    state%nvirtual = nvirtual
    state%icocc_dim = icocc_dim
    state%icvir_dim = icvir_dim
    state%cocc_dim = cocc_dim
    state%cvir_dim = cvir_dim

    allocate(iptr(numat), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = nbonds
    state%nbonds = c_loc(iptr)
    
    allocate(iptr2(9,numat), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr2 = ibonds
    state%ibonds = c_loc(iptr2)

    allocate(iptr(numat), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = iorbs
    state%iorbs = c_loc(iptr)

    allocate(iptr(noccupied), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = ncf
    state%ncf = c_loc(iptr)

    allocate(iptr(nvirtual), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = nce
    state%nce = c_loc(iptr)

    allocate(iptr(icocc_dim), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = icocc
    state%icocc = c_loc(iptr)

    allocate(iptr(icvir_dim), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    iptr = icvir
    state%icvir = c_loc(iptr)

    allocate(rptr(cocc_dim), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    rptr = cocc
    state%cocc = c_loc(rptr)

    allocate(rptr(cvir_dim), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOZYME_SAVE")
      return
    end if
    rptr = cvir
    state%cvir = c_loc(rptr)
  end subroutine mozyme_save

  ! load MOZYME density matrix, or construct initial guess
  module subroutine mozyme_load(state)
    type(mozyme_state), intent(in) :: state
    integer :: status, i, j, k, l
    integer(c_int), pointer :: iptr(:), iptr2(:,:)
    real(c_double), pointer :: rptr(:)

    if(state%numat > 0) then
      if (state%numat /= numat) then
        call mopend("Attempting to load incompatible MOZYME state")
        return
      end if
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
      if (allocated(nncf)) deallocate(nncf)
      if (allocated(nnce)) deallocate(nnce)
      if (allocated(ncocc)) deallocate(ncocc)
      if (allocated(ncvir)) deallocate(ncvir)
      allocate(nbonds(numat), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(ibonds(9,numat), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(iorbs(numat), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(ncf(noccupied), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(nce(nvirtual), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(icocc(icocc_dim), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(icvir(icvir_dim), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(cocc(cocc_dim), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(cvir(cvir_dim), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(nncf(noccupied), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(nnce(nvirtual), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(ncocc(noccupied), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      allocate(ncvir(nvirtual), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOZYME_LOAD")
        return
      end if
      call c_f_pointer(state%nbonds, iptr, [numat])
      nbonds = iptr
      call c_f_pointer(state%ibonds, iptr2, [9,numat])
      ibonds = iptr2
      call c_f_pointer(state%iorbs, iptr, [numat])
      iorbs = iptr
      call c_f_pointer(state%ncf, iptr, [noccupied])
      ncf = iptr
      call c_f_pointer(state%nce, iptr, [nvirtual])
      nce = iptr
      call c_f_pointer(state%icocc, iptr, [icocc_dim])
      icocc = iptr
      call c_f_pointer(state%icvir, iptr, [icvir_dim])
      icvir = iptr
      call c_f_pointer(state%cocc, rptr, [cocc_dim])
      cocc = rptr
      call c_f_pointer(state%cvir, rptr, [cvir_dim])
      cvir = rptr
      ! reconstruct nncf, nnce, ncocc, & ncvir
      j = 0
      do i = 1, noccupied
        nncf(i) = j
        j = j + ncf(i)
      end do
      j = 0
      do i = 1, nvirtual
        nnce(i) = j
        j = j + nce(i)
      end do
      k = 0
      do i = 1, noccupied
        ncocc(i) = k
        l = 0
        do j = nncf(i) + 1, nncf(i) + ncf(i)
          l = l + iorbs(icocc(j))
        end do
        k = k + l
      end do
      k = 0
      do i = 1, nvirtual
        ncvir(i) = k
        l = 0
        do j = nnce(i) + 1, nnce(i) + nce(i)
          l = l + iorbs(icvir(j))
        end do
        k = k + l
      end do
    else
      call geochk()
    end if
  end subroutine mozyme_load

end submodule mopac_api_saveload
