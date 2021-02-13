      MODULE dvfill_I   
      INTERFACE
      SUBROUTINE dvfill (NPPA, DIRVEC) 
      integer, INTENT(IN) :: NPPA 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(INOUT) :: DIRVEC 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
