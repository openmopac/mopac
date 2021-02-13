      MODULE copym_I   
      INTERFACE
      SUBROUTINE copym (H, F, M) 
      integer, INTENT(IN) :: M  
      DOUBLE PRECISION, DIMENSION(M,M), INTENT(IN) :: H 
      DOUBLE PRECISION, DIMENSION(M,M), INTENT(OUT) :: F 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
