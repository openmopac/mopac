      MODULE me08a_I   
      INTERFACE
      SUBROUTINE me08a (A, ALPHA, BETA, N, IA, Q)
      integer, INTENT(IN) :: N 
      integer, INTENT(IN) :: IA  
      complex, DIMENSION(IA,*), INTENT(INOUT) :: A 
      complex, DIMENSION(*), INTENT(INOUT) :: ALPHA 
      complex, DIMENSION(*), INTENT(INOUT) :: BETA       
      complex, DIMENSION(*), INTENT(INOUT) :: Q 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
