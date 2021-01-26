! ---------------------------------------------------------------------------------------
      subroutine mult_symm_AB(a, b, alpha, ndim, mdim, c, beta, iopc)

        Use vast_kind_param, only: double
!        Use mod_vars_cuda, only: ngpus
        Use mamult_I  
        Use iso_c_binding
#if GPU
        Use call_gemm_cublas
        Use mamult_cuda_i 
        use common_arrays_C, only : ifact
#endif
        implicit none
        Integer :: iopc,ndim,mdim           
        Integer :: i
#if GPU
        integer :: igrid, iblock
        real :: tt      
#endif
#if CC12
        real(c_float), dimension(mdim) :: a, b, c          
#else
        real(c_double), dimension(mdim) :: a, b, c  
#endif
        real(double), allocatable, dimension(:,:) :: xa, xb, xc
        
        real(double) :: alpha,beta
        
        continue    
                
! here, performs matrix multiplications

        Select case (iopc) 
        
          case (1) ! mamult
            call mamult (a, b, c, ndim, beta)      
          case (2) ! mamult_gpu
#if GPU
            igrid = 512 ; iblock = 512 
            tt = 0.0
            call mamult_gpu(a, b, c, ndim, mdim, ifact, beta, igrid, iblock, tt, 0)
#endif
          case (3) ! dgemm
            allocate (xa(ndim,ndim), xb(ndim,ndim), xc(ndim,ndim),stat=i)
        
!            forall (i=1:ndim,j=1:ndim) 
!              xa(i,j) = 0.d0
!              xb(i,j) = 0.d0
!            endforall
        
            call dtpttr( 'u', ndim, a, xa, ndim, i )        
            if (i /= 0) stop 'error in dtpttr'          
          
            call dtpttr( 'u', ndim, b, xb, ndim, i )        
            if (i /= 0) stop 'error in dtpttr'          
            
            do i = 1,ndim-1
               call dcopy(ndim-i,xa(i,i+1),ndim,xa(i+1,i),1)
               call dcopy(ndim-i,xb(i,i+1),ndim,xb(i+1,i),1)
            enddo

            if (.not.(beta == 0.d0)) then
!              forall (i=1:ndim,j=1:ndim) 
!                xc(i,j) = 0.d0
!              endforall
             
              call dtpttr( 'u', ndim, c, xc, ndim, i )        
              if (i /= 0) stop 'error in dtpttr'
              do i = 1,ndim-1
                 call dcopy(ndim-i,xc(i,i+1),ndim,xc(i+1,i),1)
              enddo
            endif
                                   
            call dgemm ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
                       & ndim)

            call dtrttp('u', ndim, xc, ndim, c, i )
            
            deallocate (xa,xb,xc,stat=i)         
          case (4) ! dgemm_gpu

            allocate (xa(ndim,ndim), xb(ndim,ndim), xc(ndim,ndim),stat=i)
        
!            forall (i=1:ndim,j=1:ndim) 
!              xa(i,j) = 0.d0
!              xb(i,j) = 0.d0
!            endforall
        
            call dtpttr( 'u', ndim, a, xa, ndim, i )        
            if (i /= 0) stop 'error in dtpttr'          
          
            call dtpttr( 'u', ndim, b, xb, ndim, i )        
            if (i /= 0) stop 'error in dtpttr'          
            
            do i = 1,ndim-1
               call dcopy(ndim-i,xa(i,i+1),ndim,xa(i+1,i),1)
               call dcopy(ndim-i,xb(i,i+1),ndim,xb(i+1,i),1)
            enddo

            if (beta /= 0.d0) then
!              forall (i=1:ndim,j=1:ndim) 
!                xc(i,j) = 0.d0
!              endforall
             
              call dtpttr( 'u', ndim, c, xc, ndim, i )        
              if (i /= 0) stop 'error in dtpttr'
              do i = 1,ndim-1
                 call dcopy(ndim-i,xc(i,i+1),ndim,xc(i+1,i),1)
              enddo
            endif

!            if (ngpus > 1 .and. ndim > 100) then
!
!               call gemm_cublas_mgpu ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
!                           & ndim)
!            
!            else
               call gemm_cublas ("N", "N", ndim, ndim, ndim, alpha, xa, ndim, xb, ndim, beta, xc, &
                           & ndim)
            
!            endif            

            call dtrttp('u', ndim, xc, ndim, c, i )
            
            deallocate (xa,xb,xc,stat=i) 
        end select
                  
        continue
        
      return  
      end subroutine mult_symm_AB

! ---------------------------------------------------------------------------------------      
