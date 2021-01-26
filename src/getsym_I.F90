
      MODULE getsym_I   
      INTERFACE
      subroutine getsym (locpar, idepfn, locdep, depmul)
      use molkst_C, only: natoms, ndep
      implicit none
      integer, dimension (3*natoms), intent (inout) :: idepfn, locpar
      integer, dimension (3*natoms), intent (inout) :: locdep
      double precision, dimension (natoms), intent (out) :: depmul
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
