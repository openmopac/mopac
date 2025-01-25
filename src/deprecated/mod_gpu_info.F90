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
        
