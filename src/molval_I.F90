      MODULE molval_I   
      INTERFACE
      SUBROUTINE molval (C, P, RHFUHF) 
      use molkst_C, only : norbs, mpack
      DOUBLE PRECISION, DIMENSION(NORBS,NORBS), INTENT(IN) :: C 
      DOUBLE PRECISION, DIMENSION(mpack), INTENT(IN) :: P 
      DOUBLE PRECISION, INTENT(IN) :: RHFUHF 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
