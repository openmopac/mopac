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


! For GPU MOPAC 

!
!==============================================================================
!   
! Interface for (t)asum CUBLAS
!
  module call_asum_cublas
    interface asum_cublas_gpu
        subroutine asum_cublas(n,vecx,incx,res) bind(c, &
            & name='call_asum_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx            
            real(c_double), dimension(n) :: vecx
            real(c_double):: res           
        end subroutine asum_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)axpy CUBLAS
!
  module call_axpy_cublas
    interface axpy_cublas_gpu
        subroutine axpy_cublas(n,alpha,vecx,incx,vecy,incy) bind(c, &
            & name='call_axpy_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx,incy          
            real(c_double), dimension(n) :: vecx,vecy
            real(c_double), value :: alpha      
        end subroutine axpy_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)copy CUBLAS
!
  module call_copy_cublas
    interface copy_cublas_gpu
        subroutine copy_cublas(n,vecx,incx,vecy,incy) bind(c, &
            & name='call_copy_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx,incy          
            real(c_double), dimension(n) :: vecx,vecy       
        end subroutine copy_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)dot CUBLAS
!
  module call_dot_cublas
    interface dot_cublas_gpu
        subroutine dot_cublas(n,vecx,incx,vecy,incy,res) bind(c, &
            & name='call_dot_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx,incy        
            real(c_double), dimension(n) :: vecx,vecy
            real(c_double) :: res       
        end subroutine dot_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)gemm CUBLAS
!
  module call_gemm_cublas
    interface gemm_cublas_gpu
        subroutine gemm_cublas(tra, trb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(c, &
            & name='call_gemm_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: tra,trb
	        integer(c_int),value     :: m,n,k,lda,ldb,ldc        
            real(c_double), dimension(m,k) :: a
            real(c_double), dimension(k,n) :: b
            real(c_double), dimension(m,n) :: c
            real(c_double), value :: alpha,beta                     
        end subroutine gemm_cublas
    end interface

    interface gemm_cublas_multigpu
        subroutine gemm_cublas_mgpu(tra, trb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(c, &
            & name='call_gemm_cublas_mgpu')
            use iso_c_binding
            implicit none
            character(c_char), value :: tra,trb
            integer(c_int),value     :: m,n,k,lda,ldb,ldc
#if CC12
            real(c_float), dimension(m,k) :: a
            real(c_float), dimension(k,n) :: b
            real(c_float), dimension(m,n) :: c
            real(c_float), value :: alpha,beta
#else
            real(c_double), dimension(m,k) :: a
            real(c_double), dimension(k,n) :: b
            real(c_double), dimension(m,n) :: c
            real(c_double), value :: alpha,beta
#endif
        end subroutine gemm_cublas_mgpu
    end interface

  end module


  module call_gemm_phigemm
    interface gemm_phigemm_gpu
        subroutine gemm_phigemm(tra, trb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(c, &
            & name='phigemm')
            use iso_c_binding
            implicit none
            character(c_char), value :: tra,trb
	          integer(c_int),value     :: m,n,k,lda,ldb,ldc       
            real(c_double), dimension(m,k) :: a
            real(c_double), dimension(k,n) :: b
            real(c_double), dimension(m,n) :: c
            real(c_double), value :: alpha,beta                  
        end subroutine gemm_phigemm
    end interface
  end module

! Interface for (t)gemm CUBLAS THRUST
!
  module call_gemm_cublas_thrust
    interface gemm_cublas_gpu_thrust
        subroutine gemm_cublas_thrust(tra, trb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc) bind(c, &
            & name='call_gemm_cublas_thrust')
            use iso_c_binding
            implicit none
            character(c_char), value :: tra,trb
	          integer(c_int),value     :: m,n,k,lda,ldb,ldc        
            real(c_double), dimension(m,k) :: a
            real(c_double), dimension(k,n) :: b
            real(c_double), dimension(m,n) :: c
            real(c_double), value :: alpha,beta                  
        end subroutine gemm_cublas_thrust
    end interface
  end module

!
!==============================================================================
!   
! Interface for (t)rot CUBLAS
!
  module call_rot_cublas
    interface rot_cublas_gpu
        subroutine rot_cublas(n,vecj,k,veci,l,alpha,beta) bind(c, &
            & name='call_rot_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, k,l           
            real(c_double), dimension(n) :: vecj,veci
            real(c_double), value :: alpha,beta           
        end subroutine
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)gemv CUBLAS
!
  module call_gemv_cublas
    interface gemv_cublas_gpu
        subroutine gemv_cublas(tra,m,n,alpha,a,lda,vecx,incx,beta,vecy,incy) bind(c, &
            & name='call_gemv_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: tra
	        integer(c_int),value     :: m,n,lda,incx,incy	        
            real(c_double), dimension(lda,*) :: a
            real(c_double), dimension(*) :: vecx,vecy
            real(c_double), value :: alpha,beta
        end subroutine gemv_cublas
    end interface
  end module
             
!
!==============================================================================
!   
! Interface for (t)ger CUBLAS
!
  module call_ger_cublas
    interface ger_cublas_gpu
        subroutine ger_cublas(m,n,k,alpha,vecx,incx,vecy,incy,a,lda) bind(c, &
            & name='call_ger_cublas')
            use iso_c_binding
            implicit none
	          integer(c_int),value     :: m,n,k,lda,incx,incy        
            real(c_double), dimension(lda,n) :: a
            real(c_double), dimension(m) :: vecx
            real(c_double), dimension(n) :: vecy
            real(c_double), value :: alpha
        end subroutine ger_cublas
    end interface
  end module

!
!==============================================================================
!   
! Interface for (t)nrm2 CUBLAS
!
  module call_nrm2_cublas
    interface nrm2_cublas_gpu
        subroutine nrm2_cublas(n,vecx,incx,res) bind(c, &
            & name='call_nrm2_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx           
            real(c_double), dimension(n) :: vecx
            real(c_double) :: res          
        end subroutine nrm2_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)scal CUBLAS
!
  module call_scal_cublas
    interface scal_cublas_gpu
        subroutine scal_cublas(n,alpha,vecx,incx) bind(c, &
            & name='call_scal_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx           
            real(c_double), dimension(n) :: vecx
            real(c_double), value :: alpha           
        end subroutine scal_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)swap CUBLAS
!
  module call_swap_cublas
    interface swap_cublas_gpu
        subroutine swap_cublas(n,vecx,incx,vecy,incy) bind(c, &
            & name='call_swap_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx,incy          
            real(c_double), dimension(n) :: vecx,vecy          
        end subroutine swap_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)trmm CUBLAS
!
  module call_trmm_cublas
    interface trmm_cublas_gpu
        subroutine trmm_cublas(side,uplo,tra,diag,m,n,alpha,a,lda,b,ldb,c,ldc) bind(c, &
            & name='call_trmm_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: side,uplo,tra,diag
	          integer(c_int),value     :: m,n,lda,ldb,ldc       
            real(c_double), dimension(lda,*) :: a
            real(c_double), dimension(ldb,n) :: b
            real(c_double), dimension(ldc,n) :: c
            real(c_double), value :: alpha
        end subroutine trmm_cublas
    end interface
  end module

!
!==============================================================================
!   
! Interface for (t)trmv CUBLAS
!
  module call_trmv_cublas
    interface trmv_cublas_gpu
        subroutine trmv_cublas(uplo,tra,diag,n,a,lda,vecx,incx) bind(c, &
            & name='call_trmv_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: uplo,tra,diag
	          integer(c_int),value     :: n,lda,incx      
            real(c_double), dimension(lda,n) :: a
            real(c_double), dimension(n) :: vecx
        end subroutine trmv_cublas
    end interface
  end module

!
!==============================================================================
!   
! Interface for (t)trsm CUBLAS
!
  module call_trsm_cublas
    interface trsm_cublas_gpu
        subroutine trsm_cublas(side,uplo,tra,diag,m,n,alpha,a,lda,b,ldb) bind(c, &
            & name='call_trsm_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: side, uplo,tra,diag
	          integer(c_int),value     :: m,n,lda,ldb        
            real(c_double), dimension(lda,*) :: a
            real(c_double), dimension(ldb,n) :: b
            real(c_double), value :: alpha
        end subroutine trsm_cublas
    end interface
  end module

!
!==============================================================================
!   
! Interface for i(t)amax CUBLAS
!
  module call_iamax_cublas
    interface iamax_cublas_gpu
        subroutine iamax_cublas(n,vecx,incx,res) bind(c, &
            & name='call_iamax_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx            
            real(c_double), dimension(n) :: vecx
            real(c_double):: res           
        end subroutine iamax_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for i(t)amin CUBLAS
!
  module call_iamin_cublas
    interface iamin_cublas_gpu
        subroutine iamin_cublas(n,vecx,incx,res) bind(c, &
            & name='call_iamin_cublas')
            use iso_c_binding
            implicit none
            integer(c_int), value :: n, incx           
            real(c_double), dimension(n) :: vecx
            real(c_double):: res          
        end subroutine iamin_cublas
    end interface
end module

!
!==============================================================================
!   
! Interface for (t)syrk CUBLAS
!
  module call_syrk_cublas
    interface syrk_cublas_gpu
        subroutine syrk_cublas(uplo,tra,n,k,alpha,a,lda,beta,c,ldc) bind(c, &
            & name='call_syrk_cublas')
            use iso_c_binding
            implicit none
            character(c_char), value :: uplo,tra
	          integer(c_int),value     :: n,k,lda,ldc       
            real(c_double), dimension(lda,*) :: a
            real(c_double), dimension(ldc,n) :: c
            real(c_double), value :: alpha,beta
        end subroutine syrk_cublas
    end interface
  end module

module call_syrk_cublas_thrust
    interface syrk_cublas_gpu_thrust
        subroutine syrk_cublas_thrust(uplo,tra,n,k,alpha,a,lda,beta,c,ldc) bind(c, &
            & name='call_syrk_cublas_thrust')
            use iso_c_binding
            implicit none
            character(c_char), value :: uplo,tra
	          integer(c_int),value     :: n,k,lda,ldc        
            real(c_double), dimension(lda,*) :: a
            real(c_double), dimension(ldc,n) :: c
            real(c_double), value :: alpha,beta
        end subroutine syrk_cublas_thrust
    end interface
  end module

!
!==============================================================================
!   

  module handle_cublas
    interface handle_init
        subroutine init_handle() bind(c, &
            & name='create_handle')
            implicit none                  
        end subroutine init_handle
    end interface

    interface handle_destroy
        subroutine destroy_handle() bind(c, name='destroy_handle')
            implicit none
        end subroutine destroy_handle
    end interface
  end module

  module cuda_alloc
    interface cudaMallocHost
      integer (C_INT) function cudaMallocHost(buffer, size)  bind(C,name="cudaMallocHost")
        use iso_c_binding
        implicit none
        type (C_PTR)  :: buffer
        integer (C_LONG), value :: size
      end function cudaMallocHost
    end interface

    interface cudaFreeHost
      integer (C_INT) function cudaFreeHost(buffer)  bind(C,name="cudaFreeHost")
        use iso_c_binding
        implicit none
        type (C_PTR), value :: buffer
      end function cudaFreeHost
    end interface
  end module cuda_alloc

!
!==============================================================================
!   
! Interface for MAGMA
!
  
#if (MAGMA .and. LINUX)

  module initMagma
	interface magmaDsyevd
		
		subroutine magma_dsyevd_Driver1(ngpus,opt1, opt2,n,eigenvecs,m,eigvals, &
			& work_tmp, lwork, iwork_tmp, liwork, info) bind(c,name="MagmaDsyevd_Driver1")
			use iso_c_binding
			implicit none
		
			character(C_CHAR),value :: opt1,opt2
			integer(c_int),value :: n,m,lwork,liwork,ngpus
			integer(c_int) :: iwork_tmp(10),info
			real(c_double) :: eigenvecs(n*n), eigvals(n),work_tmp(10)
		end subroutine
		
		subroutine magma_dsyevd_Driver2(ngpus,opt1, opt2,n,eigenvecs,m,eigvals, &
            & work, lwork, iwork, liwork, info) bind(c,name="MagmaDsyevd_Driver2")
			use iso_c_binding
			implicit none
		
			character(C_CHAR),value :: opt1,opt2
            integer(c_int),value :: n,m,lwork,liwork,ngpus
            integer(c_int) :: iwork(liwork),info
            real(c_double) :: eigenvecs(n*n), eigvals(n),work(lwork)
		end subroutine
		
	end interface
  end module
      
#endif

