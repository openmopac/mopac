      MODULE upcase_I   
      INTERFACE
      SUBROUTINE upcase (KEYWRD, N) 
      CHARACTER (LEN = *), INTENT(INOUT) :: KEYWRD 
      INTEGER, INTENT(IN) :: N 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
