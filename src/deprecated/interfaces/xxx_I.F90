      MODULE xxx_I   
      INTERFACE
      SUBROUTINE xxx (TYPE, I, J, K, L, R) 
      CHARACTER (LEN = 1), INTENT(IN) :: TYPE 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: K 
      INTEGER, INTENT(IN) :: L 
      CHARACTER (LEN = *), INTENT(OUT) :: R 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
