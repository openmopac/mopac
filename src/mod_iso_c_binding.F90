module iso_c_binding
!
! Do NOT use this module with the compiler - it is intended for use by FORCHECK only!
!
  implicit none
  integer, parameter :: pointer_len           = int_ptr_kind()  
  integer (kind = 4), parameter :: c_intptr_t = pointer_len  
  integer (kind = 4), parameter :: c_char     = 1
  integer (kind = 4), parameter :: c_bool     = 1
  integer (kind = 4), parameter :: c_int      = 4
  integer (kind = 4), parameter :: c_float    = 4
  integer (kind = 4), parameter :: c_long     = 4
  integer (kind = 4), parameter :: c_double   = 8 
  integer (kind = 4), parameter :: c_size_t   = pointer_len  
  type, bind(c) :: c_ptr
    private
    integer(c_intptr_t) :: ptr
  end type c_ptr
end module iso_c_binding
