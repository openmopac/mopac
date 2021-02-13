      MODULE symh_I   
      INTERFACE
      SUBROUTINE symh (H, DIP, I, N, IPO) 
      USE molkst_C, only : numat 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: H 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(INOUT) :: DIP 
      integer, dimension (numat, 120), intent (in) :: ipo
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: N 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
