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

        module gpu_info
            interface
                subroutine gpuInfo(hasGpu, hasDouble, nDevices, name,name_size, totalMem, clockRate, major, minor) bind(c, name="getGPUInfo")
                    use iso_c_binding
                    implicit none
                    
                    logical(c_bool)		   :: hasGpu
                    integer(c_int)         :: nDevices
                    logical(c_bool),dimension(6)   :: hasDouble
                    integer(c_int),dimension(6)	   :: clockRate, major, minor, name_size
                    character(c_char),dimension(6) :: name
	                integer(c_size_t),dimension(6) :: totalMem
                    
                end subroutine
            end interface
        end module
        
! ******************************

        module settingGPUcard
            interface setDevice_C
                subroutine setGPU(idevice, stat) bind(c, name='setDevice')
                    use iso_c_binding
                    implicit none
                    
                    logical(c_bool)		   :: stat
                    integer(c_int), value 	:: idevice
                    
                end subroutine
            end interface
        end module
        
