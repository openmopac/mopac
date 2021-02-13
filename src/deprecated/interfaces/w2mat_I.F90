      MODULE w2mat_I   
      INTERFACE
      SUBROUTINE w2mat (WW, W, LMW, LIMIJ, LIMKL) 
      integer, INTENT(INOUT) :: LMW 
      integer, INTENT(IN) :: LIMIJ 
      integer, INTENT(IN) :: LIMKL 
      DOUBLE PRECISION, DIMENSION(LIMKL,LIMIJ), INTENT(IN) :: WW 
      DOUBLE PRECISION, DIMENSION(LMW,LMW), INTENT(OUT) :: W       
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
