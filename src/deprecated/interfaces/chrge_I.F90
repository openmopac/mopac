      MODULE chrge_I   
      INTERFACE
      SUBROUTINE chrge (p, Q) 
      use molkst_C, only : norbs, natoms
      double precision, dimension (natoms), intent(out) :: q 
      double precision, dimension ((norbs*(norbs+1))/2), intent(in) :: p 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
