      MODULE sort_I   
      INTERFACE
      SUBROUTINE sort (VAL, VEC, N) 
      INTEGER, INTENT(IN) :: N 
      REAL, DIMENSION(*), INTENT(INOUT) :: VAL 
      COMPLEX, DIMENSION(N,*), INTENT(INOUT) :: VEC       
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
