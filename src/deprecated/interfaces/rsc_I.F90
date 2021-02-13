      MODULE rsc_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION rsc (K, NA, EA, NB, EB, NC, EC, ND, ED) 
      integer, INTENT(IN) :: K 
      integer, INTENT(IN) :: NA 
      DOUBLE PRECISION, INTENT(IN) :: EA 
      integer, INTENT(IN) :: NB 
      DOUBLE PRECISION, INTENT(IN) :: EB 
      integer, INTENT(IN) :: NC 
      DOUBLE PRECISION, INTENT(IN) :: EC 
      integer, INTENT(IN) :: ND 
      DOUBLE PRECISION, INTENT(IN) :: ED 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
