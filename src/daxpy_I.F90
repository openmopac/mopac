      MODULE daxpy_I   
      INTERFACE
      SUBROUTINE daxpy (N, DA, DX, INCX, DY, INCY) 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, INTENT(IN) :: DA 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DX 
      integer, INTENT(IN) :: INCX 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: DY 
      integer, INTENT(IN) :: INCY 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
