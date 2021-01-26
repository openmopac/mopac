      MODULE wstore_I   
      INTERFACE
      subroutine wstore (w, kr, ni, ilim)
      integer, intent (in) :: ilim
      double precision, dimension (ilim, ilim), intent (out) :: w
      integer, intent (inout) :: kr
      integer, intent (in) :: ni
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
