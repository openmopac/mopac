      subroutine mxm(a, nar, b, nbr, c, ncc) 
      implicit none
      integer  :: nar 
      integer  :: nbr 
      integer  :: ncc 
      double precision  :: a(nar,nbr) 
      double precision  :: b(nbr,ncc) 
      double precision  :: c(nar,ncc) 
!
!     RECTANGULAR MATRIX PRODUCT C=A*B.
!     EACH MATRIX IS ENTIRELY FULLFILLED AND PACKED.
!
      call dgemm ('N', 'N', nar, ncc, nbr, 1.0D0, a, nar, b, nbr, 0.0D0, c, nar)
      return  
      end subroutine mxm 


