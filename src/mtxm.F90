      subroutine mtxm(a, nar, b, nbr, c, ncc)  
      USE vast_kind_param, ONLY:  double 
      Use mod_vars_cuda, only: lgpu
#if GPU
      Use call_gemm_cublas
      Use mod_vars_cuda, only: prec
#endif
      implicit none
      integer  :: nar 
      integer  :: nbr 
      integer  :: ncc 
      real(double)  :: a(nbr,nar) 
      real(double)  :: b(nbr,ncc) 
      real(double)  :: c(nar,ncc)  
!-----------------------------------------------
!     MATRIX PRODUCT C(NAR,NCC) = (A(NBR,NAR))' * B(NBR,NCC)
!     ALL MATRICES RECTANGULAR , PACKED.
      if (lgpu) then
#if GPU
      call gemm_cublas ('T', 'N', nar, ncc, nbr, 1.0_prec, a,nbr, b, nbr, 0.0_prec, c, nar)   
#endif        
      else
        call dgemm ('T', 'N', nar, ncc, nbr, 1.0D0, a, nbr, b, nbr, 0.0D0, c, nar)
      endif 
      return  
      end subroutine mtxm   
