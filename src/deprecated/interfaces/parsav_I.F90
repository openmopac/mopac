      MODULE parsav_I   
      INTERFACE
      SUBROUTINE parsav (MODE, N, M, Q, R, efslst, xlast, iiium) 
      INTEGER, INTENT(IN) :: MODE 
      INTEGER, INTENT(IN) :: N 
      INTEGER, INTENT(IN) :: M 
      DOUBLE PRECISION, DIMENSION(N,N), INTENT(IN) :: Q 
      DOUBLE PRECISION, DIMENSION(N + 2,N + 2), INTENT(IN) :: R
      integer, dimension(6) :: iiium
      double precision, dimension(n) :: efslst, xlast
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
