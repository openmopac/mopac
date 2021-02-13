      MODULE fsub_I   
      INTERFACE
      SUBROUTINE fsub (N, X, FVAL) 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, INTENT(IN) :: X 
      DOUBLE PRECISION, INTENT(OUT) :: FVAL 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
