      MODULE daread_I   
      INTERFACE
      SUBROUTINE daread (V, LEN, NREC) 
      integer, INTENT(IN) :: LEN 
      integer, INTENT(IN) :: NREC 
      DOUBLE PRECISION, DIMENSION(LEN) :: V 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
