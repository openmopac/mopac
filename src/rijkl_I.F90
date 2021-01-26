      MODULE rijkl_I   
      INTERFACE
      DOUBLE PRECISION FUNCTION rijkl (NI, NJ, IJ, KL, LI, LJ, LK, LL, IC, R) 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      integer, INTENT(IN) :: IJ 
      integer, INTENT(IN) :: KL 
      integer, INTENT(IN) :: LI 
      integer, INTENT(IN) :: LJ 
      integer, INTENT(IN) :: LK 
      integer, INTENT(IN) :: LL 
      integer, INTENT(IN) :: IC 
      DOUBLE PRECISION :: R 
      END FUNCTION  
      END INTERFACE 
      END MODULE 
