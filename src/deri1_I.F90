      MODULE deri1_I   
      INTERFACE
      SUBROUTINE deri1 (NUMBER, GRAD, F, MINEAR, fd, SCALAR, work) 
      use molkst_C, only : mpack, norbs   
      INTEGER, INTENT(IN) :: NUMBER 
      DOUBLE PRECISION, INTENT(OUT) :: GRAD 
      DOUBLE PRECISION, DIMENSION(mpack), INTENT(INOUT) :: F 
      double precision  :: fd(mpack) 
      double precision, intent(out) :: work(norbs, norbs)
      INTEGER, INTENT(IN) :: MINEAR 
      DOUBLE PRECISION, DIMENSION(mpack), INTENT(IN) :: SCALAR 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
