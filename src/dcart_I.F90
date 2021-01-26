      MODULE dcart_I   
      INTERFACE
      SUBROUTINE dcart (coord, dxyz) 
      use molkst_C, only : numat
      double precision, intent(in)  :: coord(3,numat) 
      double precision, intent(out)  :: dxyz(3,numat)   
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
