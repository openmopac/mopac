      subroutine mxmt(a, nar, b, nbr, c, ncc) 
      implicit none
      integer  :: nar 
      integer  :: nbr 
      integer  :: ncc 
      double precision  :: a(nar,nbr) 
      double precision  :: b(ncc,nbr) 
      double precision  :: c(nar,ncc)            
!
!     MATRIX PRODUCT C(NAR,NCC) = A(NAR,NBR) * (B(NCC,NBR))'
!     ALL MATRICES RECTANGULAR , PACKED.
!    
      call dgemm ('N', 'T', nar, ncc, nbr, 1.0D0, a, nar, b, ncc, 0.0D0, c, nar)
      return  
      end subroutine mxmt 
 
