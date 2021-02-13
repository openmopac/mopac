      MODULE beopor_I   
      INTERFACE
      SUBROUTINE beopor (IWFLB, MAXITU, BTOL, UA, UB, F, GA, GB, T, H1, D, DA&
        , UAB, UOLD1, G, X, C, W) 
      use molkst_C, only : norbs
      integer, INTENT(IN) :: IWFLB 
      integer, INTENT(IN) :: MAXITU 
      DOUBLE PRECISION :: BTOL 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UA 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UB 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: F 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: GA 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: GB 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: T 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: H1 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: D 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: DA 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UAB 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UOLD1 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: X 
      DOUBLE PRECISION, DIMENSION(*), intent(in) :: C 
      DOUBLE PRECISION, DIMENSION(*) :: W 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
