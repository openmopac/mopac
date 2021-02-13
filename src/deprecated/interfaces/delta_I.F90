      MODULE delta_I   
      INTERFACE
      INTEGER FUNCTION delta (A, B, C, Y) 
      integer, INTENT(IN) :: A 
      integer, INTENT(IN) :: B 
      integer, INTENT(IN) :: C 
      integer, INTENT(IN) :: Y 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
