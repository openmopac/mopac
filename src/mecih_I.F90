      MODULE mecih_I   
      INTERFACE
      SUBROUTINE mecih (DIAG, CIMAT, NMOS, LAB, XY) 
!GBR
      integer :: nmos      
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: DIAG 
      DOUBLE PRECISION, DIMENSION(*), INTENT(out) :: cimat
      DOUBLE PRECISION, DIMENSION(*), INTENT(in) :: xy 
      INTEGER, INTENT(IN) :: LAB 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
