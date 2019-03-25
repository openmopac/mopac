      MODULE mult_symm_ab_I   
      INTERFACE
      SUBROUTINE mult_symm_ab (a, b, alpha, ndim, mdim, c, beta, iopc)
      USE vast_kind_param,ONLY: DOUBLE 
      INTEGER, INTENT(IN) :: ndim,mdim,iopc

      REAL(DOUBLE), DIMENSION(mdim), INTENT(IN) :: a, b
      REAL(DOUBLE), DIMENSION(mdim) :: c
      REAL(DOUBLE), INTENT(IN) :: alpha,beta
      
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
