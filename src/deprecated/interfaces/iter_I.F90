      MODULE iter_I   
      INTERFACE
      SUBROUTINE iter (EE, FULSCF, RAND) 
      DOUBLE PRECISION, intent(out) :: EE 
      logical, intent(in) :: FULSCF, RAND 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
