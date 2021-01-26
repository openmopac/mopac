      MODULE dtrmv_I   
      INTERFACE
      SUBROUTINE dtrmv (UPLO, TRANS, DIAG, N, A, LDA, X, INCX) 
      character (LEN = *) :: UPLO 
      character (LEN = *) :: TRANS 
      character (LEN = *) :: DIAG 
      integer, INTENT(IN) :: N, LDA 
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: X 
      integer, INTENT(IN) :: INCX 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
