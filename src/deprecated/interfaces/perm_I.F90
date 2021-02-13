      MODULE perm_I   
      INTERFACE
      SUBROUTINE perm (IPERM, NELS, NMOS, NPERMS, LIMCI)  
      use meci_C, only : maxci     
      INTEGER, INTENT(INOUT) :: NPERMS 
      INTEGER, INTENT(IN) :: LIMCI, nels, nmos
      INTEGER, DIMENSION(NMOS,MAXCI*4), INTENT(INOUT) :: IPERM
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
