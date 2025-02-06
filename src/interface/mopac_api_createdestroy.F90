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

! Diskless/stateless Application Programming Interface (API) to core MOPAC operations
submodule (mopac_api) mopac_api_createdestroy
  implicit none

#ifdef MOPAC_API_MALLOC_C
  interface
    function malloc(num) bind(c)
      use iso_c_binding
      integer(c_size_t), value :: num
      integer(c_ptr) :: malloc
    end function malloc
    subroutine free(ptr) bind(c)
      use iso_c_binding
      integer(c_ptr), value :: ptr
    end subroutine free
  end interface
#endif

contains

  ! allocate memory for mopac_state
  module subroutine create_mopac_state(state) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: create_mopac_state
#endif
    type(mopac_state), intent(inout) :: state
    if (state%mpack /= 0) then
      state%pa = create_real(state%mpack)
      if (state%uhf == 1) state%pb = create_real(state%mpack)
    end if
  end subroutine create_mopac_state

  ! allocate memory for mozyme_state
  module subroutine create_mozyme_state(state) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: create_mozyme_state
#endif
    type(mozyme_state), intent(inout) :: state
    if (state%numat /= 0) then
      state%nbonds = create_int(state%numat)
      state%ibonds = create_int2(9,state%numat)
      state%iorbs = create_int(state%numat)
      state%ncf = create_int(state%noccupied)
      state%nce = create_int(state%nvirtual)
      state%icocc = create_int(state%icocc_dim)
      state%icvir = create_int(state%icvir_dim)
      state%cocc = create_real(state%cocc_dim)
      state%cvir = create_real(state%cvir_dim)
    end if
  end subroutine create_mozyme_state

  ! deallocate memory in mopac_properties
  module subroutine destroy_mopac_properties(properties) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: destroy_mopac_properties
#endif
    type(mopac_properties), intent(in) :: properties
    integer :: i
    type(c_ptr), pointer :: pptr(:)
    call destroy_real(properties%charge)
    call destroy_int(properties%bond_index)
    call destroy_int(properties%bond_atom)
    call destroy_real(properties%bond_order)
    call destroy_real(properties%coord_update)
    call destroy_real(properties%coord_deriv)
    call destroy_real(properties%freq)
    call destroy_real(properties%disp)
    call destroy_real(properties%lattice_update)
    call destroy_real(properties%lattice_deriv)
    if (properties%nerror > 0) then
      call c_f_pointer(properties%error_msg, pptr, [properties%nerror])
      do i=1, properties%nerror
        call destroy_char(pptr(i))
      end do
      call destroy_ptr(properties%error_msg)
    end if
  end subroutine destroy_mopac_properties

  ! deallocate memory in mopac_state
  module subroutine destroy_mopac_state(state) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: destroy_mopac_state
#endif
    type(mopac_state), intent(in) :: state
    if (state%mpack /= 0) then
      call destroy_real(state%pa)
      if (state%uhf == 1) call destroy_real(state%pb)
    end if
  end subroutine destroy_mopac_state

  ! deallocate memory in mozyme_state
  module subroutine destroy_mozyme_state(state) bind(c)
#ifdef WIN32
!dec$ attributes dllexport :: destroy_mozyme_state
#endif
    type(mozyme_state), intent(in) :: state
    if (state%numat /= 0) then
      call destroy_int(state%nbonds)
      call destroy_int(state%ibonds)
      call destroy_int(state%iorbs)
      call destroy_int(state%ncf)
      call destroy_int(state%nce)
      call destroy_int(state%icocc)
      call destroy_int(state%icvir)
      call destroy_real(state%cocc)
      call destroy_real(state%cvir)
    end if
  end subroutine destroy_mozyme_state

  ! Unfortunately, the Intel Fortran compiler does not adhere to the Fortran standard for pointers.
  ! Memory allocated to a Fortran pointer erroneously cannot be deallocated if its memory is passed
  ! through a C pointer in a C-bound interface and then reassigned to a Fortran pointer because hidden
  ! information about memory allocation is contained within the original Fortran pointer and is not
  ! retained by the C pointer. To get around this, MOPAC's API will support both C and Fortran memory
  ! managers through the presence/absence of the MOPAC_API_MALLOC preprocessor variable.

  ! allocate memory (C or Fortran memory manager, depending on compiler)
  module function create_int(size)
    integer, intent(in) :: size
    type(c_ptr) :: create_int
#ifndef MOPAC_API_MALLOC
    integer(c_int), pointer :: ptr(:)
    integer :: status
    allocate(ptr(size), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_int = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    integer(c_int) :: mold
    dummy = malloc(c_sizeof(mold)*size)
    create_int = transfer(dummy, create_int)
#endif
  end function create_int
  module function create_int2(size, size2)
    integer, intent(in) :: size
    integer, intent(in) :: size2
    type(c_ptr) :: create_int2
#ifndef MOPAC_API_MALLOC
    integer(c_int), pointer :: ptr(:,:)
    integer :: status
    allocate(ptr(size,size2), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_int2 = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    integer(c_int) :: mold
    dummy = malloc(c_sizeof(mold)*size*size2)
    create_int2 = transfer(dummy, create_int2)
#endif
  end function create_int2
  module function create_real(size)
    integer, intent(in) :: size
    type(c_ptr) :: create_real
#ifndef MOPAC_API_MALLOC
    real(c_double), pointer :: ptr(:)
    integer :: status
    allocate(ptr(size), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_real = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    real(c_double) :: mold
    dummy = malloc(c_sizeof(mold)*size)
    create_real = transfer(dummy, create_real)
#endif
  end function create_real

  ! allocate memory & copy data (C or Fortran memory manager, depending on compiler)
  module function create_copy_int(array, size)
    integer, intent(in) :: array(:)
    integer, intent(in) :: size(1)
    type(c_ptr) :: create_copy_int
    integer(c_int), pointer :: ptr(:)
#ifndef MOPAC_API_MALLOC
    integer :: status
    allocate(ptr(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_copy_int = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    integer(c_int) :: mold
    dummy = malloc(c_sizeof(mold)*size(1))
    create_copy_int = transfer(dummy, create_copy_int)
    call c_f_pointer(create_copy_int, ptr, size)
#endif
    ptr = array(:size(1))
  end function create_copy_int
  module function create_copy_int2(array, size)
    integer, intent(in) :: array(:,:)
    integer, intent(in) :: size(2)
    type(c_ptr) :: create_copy_int2
    integer(c_int), pointer :: ptr(:,:)
#ifndef MOPAC_API_MALLOC
    integer :: status
    allocate(ptr(size(1),size(2)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_copy_int2 = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    integer(c_int) :: mold
    dummy = malloc(c_sizeof(mold)*size(1)*size(2))
    create_copy_int2 = transfer(dummy, create_copy_int2)
    call c_f_pointer(create_copy_int2, ptr, size)
#endif
    ptr = array(:size(1),:size(2))
  end function create_copy_int2
  module function create_copy_real(array, size)
    double precision, intent(in) :: array(:)
    integer, intent(in) :: size(1)
    type(c_ptr) :: create_copy_real
    real(c_double), pointer :: ptr(:)
#ifndef MOPAC_API_MALLOC
    integer :: status
    allocate(ptr(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_copy_real = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    real(c_double) :: mold
    dummy = malloc(c_sizeof(mold)*size(1))
    create_copy_real = transfer(dummy, create_copy_real)
    call c_f_pointer(create_copy_real, ptr, size)
#endif
    ptr = array(:size(1))
  end function create_copy_real
  module function create_copy_char(array, size)
    character(len=*), intent(in) :: array
    integer, intent(in) :: size(1)
    type(c_ptr) :: create_copy_char
    character(kind=c_char), pointer :: ptr(:)
    integer :: i
#ifndef MOPAC_API_MALLOC
    integer :: status
    allocate(ptr(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_copy_char = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    character(kind=c_char) :: mold
    dummy = malloc(c_sizeof(mold)*size(1))
    create_copy_char = transfer(dummy, create_copy_char)
    call c_f_pointer(create_copy_char, ptr, size)
#endif
    do i=1, size(1)-1
      ptr(i) = array(i:i)
    end do
    ptr(size(1)) = c_null_char
  end function create_copy_char
  module function create_copy_ptr(array, size)
    type(c_ptr), intent(in) :: array(:)
    integer, intent(in) :: size(1)
    type(c_ptr) :: create_copy_ptr
    type(c_ptr), pointer :: ptr(:)
#ifndef MOPAC_API_MALLOC
    integer :: status
    allocate(ptr(size(1)), stat=status)
    if (status /= 0) then
      write(*,*) "ERROR: Failed to allocate memory in MOPAC API I/O"
      stop 1
    end if
    create_copy_ptr = c_loc(ptr)
#else
    integer(c_intptr_t) :: dummy
    type(c_ptr) :: mold
    dummy = malloc(c_sizeof(mold)*size(1))
    create_copy_ptr = transfer(dummy, create_copy_ptr)
    call c_f_pointer(create_copy_ptr, ptr, size)
#endif
    ptr = array(:size(1))
  end function create_copy_ptr

  ! deallocate memory (C or Fortran memory manager, depending on compiler)
  module subroutine destroy_int(copy)
    type(c_ptr), intent(in) :: copy
#ifndef MOPAC_API_MALLOC
    integer(c_int), pointer :: ptr
    integer :: status
    if (c_associated(copy)) then
      call c_f_pointer(copy, ptr)
      deallocate(ptr, stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to deallocate memory in MOPAC API I/O"
        stop 1
      end if
    end if
#else
    integer(c_intptr_t) :: copy2
    copy2 = transfer(copy, copy2)
    call free(copy2)
#endif
  end subroutine destroy_int
  module subroutine destroy_real(copy)
    type(c_ptr), intent(in) :: copy
#ifndef MOPAC_API_MALLOC
    real(c_double), pointer :: ptr
    integer :: status
    if (c_associated(copy)) then
      call c_f_pointer(copy, ptr)
      deallocate(ptr, stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to deallocate memory in MOPAC API I/O"
        stop 1
      end if
    end if
#else
    integer(c_intptr_t) :: copy2
    copy2 = transfer(copy, copy2)
    call free(copy2)
#endif
  end subroutine destroy_real
  module subroutine destroy_char(copy)
    type(c_ptr), intent(in) :: copy
#ifndef MOPAC_API_MALLOC
    character(kind=c_char), pointer :: ptr
    integer :: status
    if (c_associated(copy)) then
      call c_f_pointer(copy, ptr)
      deallocate(ptr, stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to deallocate memory in MOPAC API I/O"
        stop 1
      end if
    end if
#else
    integer(c_intptr_t) :: copy2
    copy2 = transfer(copy, copy2)
    call free(copy2)
#endif
  end subroutine destroy_char
  module subroutine destroy_ptr(copy)
    type(c_ptr), intent(in) :: copy
#ifndef MOPAC_API_MALLOC
    type(c_ptr), pointer :: ptr
    integer :: status
    if (c_associated(copy)) then
      call c_f_pointer(copy, ptr)
      deallocate(ptr, stat=status)
      if (status /= 0) then
        write(*,*) "ERROR: Failed to deallocate memory in MOPAC API I/O"
        stop 1
      end if
    end if
#else
    integer(c_intptr_t) :: copy2
    copy2 = transfer(copy, copy2)
    call free(copy2)
#endif
  end subroutine destroy_ptr

end submodule mopac_api_createdestroy
