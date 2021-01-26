      MODULE ddot_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION ddot (N, DX, INCX, DY, INCY) 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DX 
      integer, INTENT(IN) :: INCX 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DY 
      integer, INTENT(IN) :: INCY 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
