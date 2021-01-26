      MODULE dijkl2_I   
      INTERFACE
      SUBROUTINE dijkl2 (DC) 
      use molkst_C, only : norbs
      use meci_C, only : nmos
      DOUBLE PRECISION, DIMENSION(NORBS,NMOS) :: DC 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
