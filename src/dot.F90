      real(kind(0.0d0)) function dot (x, y, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:33:41  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: x(n) 
      real(double) , intent(in) :: y(n)
!For MOPAC BLAS        
      real(double), external :: ddot
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
