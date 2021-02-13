      MODULE rotat_I   
      INTERFACE
      SUBROUTINE rotat (COORD, I, J, IX, RIJ, DEL1, IDX) 
      DOUBLE PRECISION, DIMENSION(3,25), INTENT(IN) :: COORD 
      integer, INTENT(IN) :: I 
      integer, INTENT(IN) :: J 
      integer, INTENT(IN) :: IX 
      DOUBLE PRECISION, INTENT(IN) :: RIJ 
      DOUBLE PRECISION, INTENT(IN) :: DEL1 
      integer, INTENT(IN) :: IDX 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
