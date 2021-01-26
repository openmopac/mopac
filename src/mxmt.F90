      subroutine mxmt(a, nar, b, nbr, c, ncc) 
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
      double precision  :: b(ncc,nbr) 
      double precision  :: c(nar,ncc)            
!
!     MATRIX PRODUCT C(NAR,NCC) = A(NAR,NBR) * (B(NCC,NBR))'
!     ALL MATRICES RECTANGULAR , PACKED.
!    
      if (lgpu) then
#if GPU
        call gemm_cublas ('N', 'T', nar, ncc, nbr, 1.0_prec, a,nar, b, ncc, 0.0_prec, c, nar) 
#endif           
      else
         call dgemm ('N', 'T', nar, ncc, nbr, 1.0D0, a, nar, b, ncc, 0.0D0, c, nar)
      endif 

      return  
      end subroutine mxmt 
 
