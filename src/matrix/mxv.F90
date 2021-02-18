      subroutine mxv(a, nar, vecx, nbr, vecy)     
      implicit none
      integer, parameter :: incy = 1, incx = 1
      integer  :: nar 
      integer  :: nbr 
      double precision  :: a(nar,nbr) 
      double precision  :: vecx(nbr) 
      double precision  :: vecy(nar) 
!
!     RECTANGULAR MATRIX-VECTOR PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.

      call dgemv ('N', nar, nbr, 1.0d0, a, nar, vecx, incx, 0.0d0, vecy, incy)
      return  
      end subroutine mxv 

