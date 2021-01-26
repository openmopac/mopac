      MODULE dtrmm_I   
      INTERFACE
      SUBROUTINE dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB) 
      character (LEN = *) :: SIDE 
      character (LEN = *) :: UPLO 
      character (LEN = *) :: TRANSA 
      character (LEN = *) :: DIAG 
      integer, INTENT(IN) :: M 
      integer, INTENT(IN) :: N, LDA, LDB
      DOUBLE PRECISION, INTENT(IN) :: ALPHA 
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(IN) :: A 
      DOUBLE PRECISION, DIMENSION(LDB,*), INTENT(INOUT) :: B 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
