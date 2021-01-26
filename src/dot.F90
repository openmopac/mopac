      double precision function dot (x, y, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      double precision , intent(in) :: x(n) 
      double precision , intent(in) :: y(n)
!For MOPAC BLAS        
      double precision, external :: ddot
!
      
!-----------------------------------------------
!***********************************************************************
!
!   DOT FORMS THE SCALAR PRODUCT OF TWO VECTORS.
!
!   ON INPUT     X   =    FIRST VECTOR, OF LENGTH N.
!                Y   =    SECOND VECTOR, OF LENGTH N.
!
!   ON RETURN    DOT =    DOT PRODUCT OF X AND Y.
!
!***********************************************************************
!For MOPAC BLAS              
!     dot = dot_product(x(:n),y(:n)) 
      dot = ddot(n,x(1:n),1,y(1:n),1)
!      
      return  
      end function dot 
