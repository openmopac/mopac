      MODULE printp_I   
      INTERFACE
      SUBROUTINE printp (I, PARA, VALUE, TXT) 
      integer, INTENT(IN) :: I 
      character (LEN = *), INTENT(IN) :: PARA 
      DOUBLE PRECISION, INTENT(IN) :: VALUE 
      character (LEN = *), INTENT(IN) :: TXT 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
