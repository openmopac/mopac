      MODULE alphaf_I   
      INTERFACE
      SUBROUTINE alphaf (IWFLA, ATOL, MAXITA, U, F, G, UOLD, H1, D, DA, W2, C, W) 
      use molkst_C, only : norbs
      integer, INTENT(IN) :: IWFLA 
      DOUBLE PRECISION :: ATOL 
      integer, INTENT(IN) :: MAXITA 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: F 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UOLD 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: H1 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: D 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: DA 
      DOUBLE PRECISION, DIMENSION(*), intent(out) :: W2 
      DOUBLE PRECISION, DIMENSION(*), intent(in) :: C 
      DOUBLE PRECISION, DIMENSION(*), intent(in) :: W 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
