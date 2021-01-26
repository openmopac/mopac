      MODULE nonope_I   
      INTERFACE
      SUBROUTINE nonope (U0X, U1Y, U1Z, U1X, U0Y, U0Z, G0X, G1Y, G1Z, G1X, G0Y&
        , G0Z) 
      use molkst_C, only : norbs
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U0X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U1X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U0Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: U0Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G0X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1Z 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G1X 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G0Y 
      DOUBLE PRECISION, DIMENSION(norbs,norbs) :: G0Z 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
