      MODULE chrge_I   
      INTERFACE
      SUBROUTINE chrge (p, Q) 
      USE vast_kind_param,ONLY: DOUBLE
      use molkst_C, only : norbs, natoms
      real (double), dimension (natoms), intent(out) :: q 
      real (double), dimension ((norbs*(norbs+1))/2), intent(in) :: p 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
