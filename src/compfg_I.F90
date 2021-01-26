      MODULE compfg_I   
      INTERFACE
      SUBROUTINE compfg (XPARAM, INT, ESCF, FULSCF, GRAD, LGRAD) 
      DOUBLE PRECISION, DIMENSION(*), intent (in) :: XPARAM 
      logical, intent(in) :: INT, fulscf 
      DOUBLE PRECISION, intent(out) :: ESCF 
      DOUBLE PRECISION, DIMENSION(*) :: GRAD 
      LOGICAL, intent(in) :: LGRAD 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
