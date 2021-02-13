      MODULE darea1_I   
      INTERFACE
      SUBROUTINE darea1 (V, LEN, IDAF, NS) 
      integer, INTENT(IN) :: LEN 
      integer, INTENT(IN) :: IDAF 
      integer, INTENT(IN) :: NS 
      DOUBLE PRECISION, DIMENSION(LEN) :: V 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
