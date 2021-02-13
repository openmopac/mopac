      MODULE etred3_I   
      INTERFACE
      SUBROUTINE etred3 (N, A, D, E, E2) 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: A 
      DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: D, e, e2
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
