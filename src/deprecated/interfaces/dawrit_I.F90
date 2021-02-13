      MODULE dawrit_I   
      INTERFACE
      SUBROUTINE dawrit (V, LEN, NREC) 
      integer, INTENT(IN) :: LEN 
      integer, INTENT(IN) :: NREC 
      DOUBLE PRECISION, DIMENSION(LEN) :: V   
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
