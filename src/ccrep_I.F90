      MODULE ccrep_I   
      INTERFACE
      SUBROUTINE ccrep (NI, NJ, R, ENUCLR, gab) 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      DOUBLE PRECISION, INTENT(INOUT) :: R 
      DOUBLE PRECISION, INTENT(IN) :: GAB 
      DOUBLE PRECISION, INTENT(OUT) :: ENUCLR 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
