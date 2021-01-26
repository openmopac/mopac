      MODULE lsame_I   
      INTERFACE
      LOGICAL FUNCTION lsame (CA, CB) 
      character (LEN = 1), INTENT(IN) :: CA 
      character (LEN = 1), INTENT(IN) :: CB 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
