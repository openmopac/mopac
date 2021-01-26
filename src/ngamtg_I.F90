      MODULE ngamtg_I   
      INTERFACE
      SUBROUTINE ngamtg (X, GD3, UD3, G1, U1, GS, USMD, EPS, US) 
      use molkst_C, only : norbs 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: GD3 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: UD3 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: GS 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: USMD 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: EPS 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: US 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
