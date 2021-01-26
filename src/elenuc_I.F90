      MODULE elenuc_I   
      INTERFACE
      SUBROUTINE elenuc (IA, IB, JA, JB, H) 
      integer, INTENT(IN) :: IA 
      integer, INTENT(IN) :: IB 
      integer, INTENT(IN) :: JA 
      integer, INTENT(IN) :: JB 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: H 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
