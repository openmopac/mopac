      MODULE reppd_I   
      INTERFACE
      SUBROUTINE reppd (NI, NJ, RIJ, RI, GAB) 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      DOUBLE PRECISION, INTENT(IN) :: RIJ
      DOUBLE PRECISION, INTENT(OUT) :: gab
      DOUBLE PRECISION, DIMENSION(22), INTENT(OUT) :: RI 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
