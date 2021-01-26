      MODULE mxv_I   
      INTERFACE
      SUBROUTINE mxv (a, nar, vecx, nbr, vecy)  
      INTEGER, INTENT(IN) :: NAR, NBR
      DOUBLE PRECISION, DIMENSION(NAR,NBR) :: A       
      DOUBLE PRECISION, DIMENSION(NBR) :: vecx 
      DOUBLE PRECISION, DIMENSION(NAR) :: vecy 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE

