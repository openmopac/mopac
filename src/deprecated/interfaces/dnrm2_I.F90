      MODULE dnrm2_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION dnrm2 (N, DX) 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: DX 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
