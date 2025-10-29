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

submodule (mopac_api:mopac_api_operations) mopac_api_saveload
  use Common_arrays_C, only: p, pa, pb, nbonds, ibonds
  use molkst_C, only: keywrd, uhf, mpack, numat
  use MOZYME_C, only: iorbs, noccupied, ncf, nvirtual, nce, icocc_dim, &
    icocc, icvir_dim, icvir, cocc_dim, cocc, cvir_dim, cvir, nnce, nncf, ncocc, ncvir
  implicit none

contains

  ! save MOPAC density matrices
  module subroutine mopac_save(state)
    type(mopac_state), intent(inout) :: state

    if (state%mpack > 0) call destroy_mopac_state(state)

    state%mpack = mpack
    if (mpack > 0) then
      state%pa = create_copy(pa, [mpack])
      if (uhf) then
        state%uhf = 1
        state%pb = create_copy(pb, [mpack])
      else
        state%uhf = 0
        state%pb = c_null_ptr
      end if
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
      call c_f_pointer(state%pa, rptr, [mpack])
      pa = rptr
      if (uhf) then
        if (state%uhf == 1) then
          call c_f_pointer(state%pb, rptr, [mpack])
          pb = rptr
        else
          pb = pa
        end if
        p = pa + pb
      else
        p = pa*2.d0
      end if
    end if
  end subroutine mopac_load

  ! save MOZYME density matrix
  module subroutine mozyme_save(state)
    type(mozyme_state), intent(inout) :: state

    if (state%numat > 0) call destroy_mozyme_state(state)

    state%numat = numat
    if (numat > 0) then
      state%noccupied = noccupied
      state%nvirtual = nvirtual
      state%icocc_dim = icocc_dim
      state%icvir_dim = icvir_dim
      state%cocc_dim = cocc_dim
      state%cvir_dim = cvir_dim

      state%nbonds = create_copy(nbonds, [numat])
      state%ibonds = create_copy(ibonds, [9,numat])
      state%iorbs = create_copy(iorbs, [numat])
      state%ncf = create_copy(ncf, [noccupied])
      state%nce = create_copy(nce, [nvirtual])
      state%icocc = create_copy(icocc, [icocc_dim])
      state%icvir = create_copy(icvir, [icvir_dim])
      state%cocc = create_copy(cocc, [cocc_dim])
      state%cvir = create_copy(cvir, [cvir_dim])
    end if
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
