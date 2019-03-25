      subroutine mult(c, s, vecs, n)  
                                      
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      Use mod_vars_cuda, only: lgpu
#if GPU
      Use call_gemm_cublas
      Use mod_vars_cuda, only: real_cuda, prec
#endif
      implicit none
      integer , intent(in) :: n 
      real(double) , intent(in) :: c(n,n) 
      real(double) , intent(in) :: s(n,n) 
      real(double) , intent(out) :: vecs(n,n)   
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

      if (lgpu) then
#if GPU
        call gemm_cublas ('T', 'T', n, n, n, 1.0_prec, c,n, s, n, 0.0_prec,  vecs, n)                                     
        call mkl_dimatcopy('c', 't', n, n, 1.d0, vecs, n, n) ! TODO !?
#endif
      else
        do i = 1, n 
          do j = 1, n 
            sum = 0.D0 
            do k = 1, n 
              sum = sum + c(k,i)*s(j,k) 
            end do 
            vecs(j,i) = sum 
          end do 
        end do                  
      endif          
      return  
      end subroutine mult  
