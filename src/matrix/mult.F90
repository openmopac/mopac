      subroutine mult(c, s, vecs, n)  
                                      
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
#if GPU
      Use call_gemm_cublas
      Use mod_vars_cuda, only: lgpu, real_cuda, prec
#endif
      implicit none
      integer , intent(in) :: n 
      double precision , intent(in) :: c(n,n) 
      double precision , intent(in) :: s(n,n) 
      double precision , intent(out) :: vecs(n,n)   
!
!  Local variables
!
      integer :: i, j, k
      double precision :: sum      
!-----------------------------------------------
!**********************************************************************
!
!   MULT IS USED IN THE MULLIKEN ANALYSIS ONLY. IT PERFORMS THE
!        OPERATION:-
!                                   VECS=BACK-TRANSFORMED EIGENVECTORS
!        VECS  =  C*S               C   =UN-BACK-TRANSFORMED VECTORS
!                                   S   =1/SQRT(OVERLAP MATRIX)
!
!**********************************************************************

#if GPU
      if (lgpu) then
        call gemm_cublas ('T', 'T', n, n, n, 1.0_prec, c,n, s, n, 0.0_prec,  vecs, n)                                     
        call mkl_dimatcopy('c', 't', n, n, 1.d0, vecs, n, n) ! TODO !?
      else
#endif
        do i = 1, n 
          do j = 1, n 
            sum = 0.D0 
            do k = 1, n 
              sum = sum + c(k,i)*s(j,k) 
            end do 
            vecs(j,i) = sum 
          end do 
        end do
#if GPU
      end if
#endif
      return  
      end subroutine mult  
