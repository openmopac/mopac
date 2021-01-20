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
        
