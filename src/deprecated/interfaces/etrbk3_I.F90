      MODULE etrbk3_I   
      INTERFACE
      SUBROUTINE etrbk3 (NM, N, NV, A, M, Z) 
      integer, INTENT(IN) :: NM 
      integer, INTENT(IN) :: N 
      integer, INTENT(IN) :: NV 
      DOUBLE PRECISION, DIMENSION(NV), INTENT(IN) :: A 
      integer, INTENT(IN) :: M 
      DOUBLE PRECISION, DIMENSION(NM,M) :: Z 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
