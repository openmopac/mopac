      MODULE ders_I   
      INTERFACE
      SUBROUTINE ders (M, N, RR, DEL1, DEL2, DEL3, IS, IOL) 
      integer, INTENT(IN) :: M 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, INTENT(IN) :: RR 
      DOUBLE PRECISION, INTENT(IN) :: DEL1 
      DOUBLE PRECISION, INTENT(IN) :: DEL2 
      DOUBLE PRECISION, INTENT(IN) :: DEL3 
      integer, INTENT(IN) :: IS 
      integer, INTENT(IN) :: IOL 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
