      subroutine mxm(a, nar, b, nbr, c, ncc) 
#if GPU
      Use call_gemm_cublas
      Use mod_vars_cuda, only: lgpu, prec
#endif
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
#if GPU
      if (lgpu .and. (nar >= 100 .or. ncc >= 100 .or. nbr >= 100)) then   
         call gemm_cublas ('N', 'N', nar, ncc, nbr, 1.0_prec, a,nar, b, nbr, 0.0_prec, c, nar)
      else            
#endif
            call dgemm ('N', 'N', nar, ncc, nbr, 1.0D0, a, nar, b, nbr, 0.0D0, c, nar)
#if GPU
      end if 
#endif      
      return  
      end subroutine mxm 


