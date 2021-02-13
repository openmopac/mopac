      MODULE dihed_I   
      INTERFACE
      SUBROUTINE dihed (XYZ, I, J, K, L, ANGLE) 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(IN) :: XYZ 
      INTEGER, INTENT(IN) :: I, j, k, l 
      DOUBLE PRECISION, INTENT(OUT) :: ANGLE 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
