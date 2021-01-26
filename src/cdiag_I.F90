      MODULE cdiag_I   
      INTERFACE
      SUBROUTINE cdiag (A, VALUE, VEC, N) 
      INTEGER, INTENT(IN) :: N 
      COMPLEX, DIMENSION(N,*) :: A 
      REAL, DIMENSION(*) :: VALUE 
      COMPLEX, DIMENSION(N,*) :: VEC 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
