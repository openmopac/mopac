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
