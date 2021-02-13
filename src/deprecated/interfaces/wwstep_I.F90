      MODULE wwstep_I   
      INTERFACE
      SUBROUTINE wwstep (C12, CC, WW, KMAX, LMAX) 
      DOUBLE PRECISION, DIMENSION(*) :: C12 
      DOUBLE PRECISION, DIMENSION(*) :: CC 
      DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: WW 
      integer, INTENT(IN) :: KMAX 
      integer, INTENT(IN) :: LMAX 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
