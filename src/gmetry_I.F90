      MODULE gmetry_I   
      INTERFACE
      SUBROUTINE gmetry (GEO, COORD) 
      use molkst_C, only : natoms      
      DOUBLE PRECISION, DIMENSION(3,NATOMS) :: GEO 
      DOUBLE PRECISION, DIMENSION(3,NATOMS), INTENT(OUT) :: COORD 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
