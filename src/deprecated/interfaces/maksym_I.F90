      MODULE maksym_I   
      INTERFACE
      subroutine maksym (loc, xparam, xstore)
      use molkst_C, only: nvar, ndep
    
     !.. Implicit Declarations ..
      implicit none
     !
     !.. Formal Arguments ..
      integer, dimension (2, nvar), intent (inout) :: loc
      double precision, dimension (nvar), intent (inout) :: xparam, xstore
     !
      end subroutine maksym
      END INTERFACE 
      END MODULE 
