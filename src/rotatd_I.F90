      MODULE rotatd_I   
      INTERFACE
      SUBROUTINE rotatd (NI, NJ, CI, CJ, W, KR, ENUC) 
      integer, intent(in) :: NI, NJ 
!GBR
      integer :: kr      
      DOUBLE PRECISION :: CI(3) 
      DOUBLE PRECISION :: CJ(3)
      DOUBLE PRECISION, dimension(*) :: W 
      integer :: LMW 
      DOUBLE PRECISION :: ENUC 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
