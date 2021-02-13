      MODULE packp_I   
      INTERFACE
      SUBROUTINE packp (P, PP, MN) 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: P 
      DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: PP 
      integer, INTENT(OUT) :: MN 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
