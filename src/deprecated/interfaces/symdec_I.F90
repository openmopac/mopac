      MODULE symdec_I   
      INTERFACE
      LOGICAL FUNCTION symdec (N1, IELEM) 
      INTEGER, INTENT(IN) :: N1 
      INTEGER, DIMENSION(20), INTENT(IN) :: IELEM 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
