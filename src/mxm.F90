      subroutine mxm(a, nar, b, nbr, c, ncc) 
      Use mod_vars_cuda, only: lgpu
#if GPU
      Use call_gemm_cublas
      Use mod_vars_cuda, only: prec
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
      if (lgpu .and. (nar >= 100 .or. ncc >= 100 .or. nbr >= 100)) then   
#if GPU
      call gemm_cublas ('N', 'N', nar, ncc, nbr, 1.0_prec, a,nar, b, nbr, 0.0_prec, c, nar)
#endif
      else            
         call dgemm ('N', 'N', nar, ncc, nbr, 1.0D0, a, nar, b, nbr, 0.0D0, c, nar)
      end if 
       
      return  
      end subroutine mxm 


