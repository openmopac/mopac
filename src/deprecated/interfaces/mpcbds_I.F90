      MODULE mpcbds_I   
      INTERFACE
      SUBROUTINE mpcbds (BONDAB) 
      use molkst_C, only : NUMAT 
      DOUBLE PRECISION, DIMENSION((NUMAT*(NUMAT + 1))/2), INTENT(IN) :: BONDAB       
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
