      MODULE wrdkey_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION wrdkey (KEYWRD, KEY, NK, REFKEY, NR, DEF) 
      character (LEN = 241), INTENT(IN) :: KEYWRD 
      character (LEN = *), INTENT(IN) :: KEY 
      integer, INTENT(IN) :: NK 
      character (LEN = *), INTENT(IN) :: REFKEY 
      integer, INTENT(IN) :: NR 
      DOUBLE PRECISION, INTENT(IN) :: DEF 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
