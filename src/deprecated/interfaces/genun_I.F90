      MODULE genun_I   
      INTERFACE
      SUBROUTINE genun (U, N) 
      INTEGER, INTENT(INOUT) :: N 
      DOUBLE PRECISION, DIMENSION(3,N), INTENT(OUT) :: U       
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
