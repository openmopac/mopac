      MODULE spcore_I   
      INTERFACE
      SUBROUTINE spcore (NI, NJ, R, CORE) 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      DOUBLE PRECISION, INTENT(IN) :: R 
      DOUBLE PRECISION, DIMENSION(10,2), INTENT(OUT) :: CORE 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
