      MODULE diegrd_I   
      INTERFACE
      subroutine diegrd (dxyz)
      use molkst_C, only: numat
      double precision, dimension (3,numat), intent (inout) :: dxyz 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
