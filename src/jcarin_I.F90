      MODULE jcarin_I   
      INTERFACE
      SUBROUTINE jcarin (XPARAM, STEP, PRECI, B, NCOL, IL, IU) 
      double precision, intent(in) :: step
      logical, intent(in) :: preci
      integer, intent(in) :: il, iu
      integer, intent(out) :: ncol
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: XPARAM 
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: B 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
