      MODULE tx_I   
      INTERFACE
      SUBROUTINE tx (II, KK, REP, LOGV, V) 
      integer, INTENT(IN) :: II 
      integer, INTENT(IN) :: KK 
      DOUBLE PRECISION, DIMENSION(491), INTENT(IN) :: REP 
      logical, DIMENSION(45,45), INTENT(OUT) :: LOGV 
      DOUBLE PRECISION, DIMENSION(45,45), INTENT(OUT) :: V 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
