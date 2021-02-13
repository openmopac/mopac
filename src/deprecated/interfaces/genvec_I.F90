      MODULE genvec_I   
      INTERFACE
      SUBROUTINE genvec (U, N) 
      integer, INTENT(INOUT) :: N
      DOUBLE PRECISION, DIMENSION(3,N), INTENT(OUT) :: U        
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
