      MODULE dawrt1_I   
      INTERFACE
      SUBROUTINE dawrt1 (V, LEN, IDAF, NS) 
      integer, INTENT(IN) :: LEN 
      integer, INTENT(IN) :: IDAF 
      integer, INTENT(IN) :: NS 
      DOUBLE PRECISION, DIMENSION(LEN), INTENT(IN) :: V 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
