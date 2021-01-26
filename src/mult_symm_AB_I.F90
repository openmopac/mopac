      MODULE mult_symm_ab_I   
      INTERFACE
      SUBROUTINE mult_symm_ab (a, b, alpha, ndim, mdim, c, beta, iopc)
      INTEGER, INTENT(IN) :: ndim,mdim,iopc

      DOUBLE PRECISION, DIMENSION(mdim), INTENT(IN) :: a, b
      DOUBLE PRECISION, DIMENSION(mdim) :: c
      DOUBLE PRECISION, INTENT(IN) :: alpha,beta
      
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
