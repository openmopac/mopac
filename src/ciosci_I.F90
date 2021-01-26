      MODULE ciosci_I   
      INTERFACE
      SUBROUTINE ciosci (VECTS, lroot, OSCIL, CONF) 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: VECTS 
      DOUBLE PRECISION, DIMENSION(3,*), INTENT(INOUT) :: OSCIL 
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: CONF 
      integer, intent(in) :: lroot
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
