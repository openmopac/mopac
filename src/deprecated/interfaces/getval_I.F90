      MODULE getval_I   
      INTERFACE
      SUBROUTINE getval (LINE, X, T) 
      CHARACTER (LEN = 80), INTENT(IN) :: LINE 
      DOUBLE PRECISION, INTENT(OUT) :: X 
      CHARACTER (LEN = 12), INTENT(OUT) :: T 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
