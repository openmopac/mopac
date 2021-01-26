      MODULE dtrti2_I   
      INTERFACE
      SUBROUTINE dtrti2 (UPLO, DIAG, N, A, LDA, INFO) 
      character (LEN = 1) :: UPLO 
      character (LEN = 1) :: DIAG 
      integer, INTENT(IN) :: N, LDA
      DOUBLE PRECISION, DIMENSION(LDA,*), INTENT(INOUT) :: A 
      integer, INTENT(OUT) :: INFO 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
