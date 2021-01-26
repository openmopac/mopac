      MODULE worder_I   
      INTERFACE
      SUBROUTINE worder (CC, LM7, IOUT2, W) 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: CC 
      integer, INTENT(IN) :: LM7 
      integer, INTENT(IN) :: IOUT2 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: W 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
