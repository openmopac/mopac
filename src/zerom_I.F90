      MODULE zerom_I   
      INTERFACE
      SUBROUTINE zerom (X, M) 
      integer, INTENT(IN) :: M
      DOUBLE PRECISION, DIMENSION(M,M), INTENT(OUT) :: X        
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
