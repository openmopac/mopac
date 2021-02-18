      subroutine mult(c, s, vecs, n)  
                                      
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
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
      do i = 1, n 
        do j = 1, n 
          sum = 0.D0 
          do k = 1, n 
            sum = sum + c(k,i)*s(j,k) 
          end do 
          vecs(j,i) = sum 
        end do 
      end do
      return  
      end subroutine mult  
