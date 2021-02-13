      MODULE delmol_I   
      INTERFACE
      SUBROUTINE delmol (COORD, I, J, NI, NJ, IA, ID, JA, JD, IX, RIJ, TOMB&
        , ISP) 
      DOUBLE PRECISION, DIMENSION(3,25) :: COORD 
      integer :: I 
      integer :: J 
      integer, INTENT(IN) :: NI 
      integer, INTENT(IN) :: NJ 
      integer, INTENT(IN) :: IA 
      integer, INTENT(IN) :: ID 
      integer, INTENT(IN) :: JA 
      integer, INTENT(IN) :: JD 
      integer :: IX 
      DOUBLE PRECISION :: RIJ 
      DOUBLE PRECISION :: TOMB 
      integer, INTENT(INOUT) :: ISP 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
