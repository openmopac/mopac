      MODULE einvit_I   
      INTERFACE
      SUBROUTINE einvit (NM, N, D, E, E2, M, W, IND, Z, IERR, RV1, RV2, RV3&
        , RV4, RV6) 
      integer, INTENT(IN) :: NM 
      integer, INTENT(IN) :: N 
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: D 
      DOUBLE PRECISION, DIMENSION(N + 1), INTENT(IN) :: E 
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: E2 
      integer, INTENT(IN) :: M 
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: W 
      integer, DIMENSION(M), INTENT(IN) :: IND 
      DOUBLE PRECISION, DIMENSION(NM,M), INTENT(OUT) :: Z 
      integer, INTENT(OUT) :: IERR 
      DOUBLE PRECISION, DIMENSION(N), INTENT(INOUT) :: RV1, rv2, rv3, rv4, rv6
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
