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

      module vast_kind_param                                        
         integer, parameter :: byte_log = selected_int_kind(2)      
         integer, parameter :: short_log = selected_int_kind(4)     
         integer, parameter :: long_log = selected_int_kind(18)     
         integer, parameter :: byte = selected_int_kind(2)          
         integer, parameter :: short = selected_int_kind(4)         
         integer, parameter :: long = selected_int_kind(18)         
         integer, parameter :: double = selected_real_kind(8)      
         integer, parameter :: extended = selected_real_kind(30)    
         integer, parameter :: double_ext = selected_real_kind(50)  
         integer, parameter :: dble_complex = selected_real_kind(14)
         integer, parameter :: ext_complex = selected_real_kind(30) 
      end module vast_kind_param                                    
