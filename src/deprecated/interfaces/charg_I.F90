      MODULE charg_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION charg (R, L1, L2, M, DA, DB, ADD) 
      DOUBLE PRECISION, INTENT(IN) :: R 
      integer, INTENT(IN) :: L1 
      integer, INTENT(IN) :: L2 
      integer, INTENT(IN) :: M 
      DOUBLE PRECISION, INTENT(IN) :: DA 
      DOUBLE PRECISION, INTENT(IN) :: DB 
      DOUBLE PRECISION, INTENT(IN) :: ADD 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
