      subroutine mtxm(a, nar, b, nbr, c, ncc)  
      implicit none
      integer  :: nar 
      integer  :: nbr 
      integer  :: ncc 
      double precision  :: a(nbr,nar) 
      double precision  :: b(nbr,ncc) 
      double precision  :: c(nar,ncc)  
!-----------------------------------------------
!     MATRIX PRODUCT C(NAR,NCC) = (A(NBR,NAR))' * B(NBR,NCC)
!     ALL MATRICES RECTANGULAR , PACKED.

      call dgemm ('T', 'N', nar, ncc, nbr, 1.0D0, a, nbr, b, nbr, 0.0D0, c, nar)
      return  
      end subroutine mtxm   
