      MODULE search_I   
      INTERFACE
      subroutine search(xparam, alpha, sig, nvar, gmin, funct, amin, anext ) 
      integer  :: nvar 
      double precision , intent(inout) :: alpha 
      double precision , intent(inout) :: gmin 
      double precision  :: funct, amin, anext
      double precision  :: xparam(*) 
      double precision , intent(in) :: sig(*) 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
