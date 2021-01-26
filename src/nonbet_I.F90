      MODULE nonbet_I   
      INTERFACE
      SUBROUTINE nonbet (U1X, U1Y, U1Z, U2X, U2Y, U2Z, G1X, G1Y, G1Z, G2X, G2Y&
        , G2Z) 
      use molkst_C, only : norbs
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U2X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U2Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U2Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G2X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G2Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G2Z 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
