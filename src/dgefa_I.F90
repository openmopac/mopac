      MODULE dgefa_I   
      INTERFACE
      SUBROUTINE dgefa (A, LDA, N, IPVT, INFO) 
      integer, INTENT(IN) :: LDA 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, DIMENSION(LDA,N), INTENT(INOUT) :: A      
      integer, DIMENSION(N), INTENT(OUT) :: IPVT 
      integer, INTENT(OUT) :: INFO 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
