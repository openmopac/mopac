      subroutine mxv(a, nar, vecx, nbr, vecy) 
#if GPU
      Use mod_vars_cuda, only: lgpu, prec
      Use call_gemv_cublas
#endif      
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

#if GPU
      if (lgpu .and. (nar >= 100 .or. nbr >= 100)) then      
         call gemv_cublas('N',nar,nbr,1.0_prec,a,nar,vecx,incx,0.0_prec,vecy,incy)
      else                  
#endif
         call dgemv ('N', nar, nbr, 1.0d0, a, nar, vecx, incx, 0.0d0, vecy, incy)
#if GPU
      end if 
#endif
      return  
      end subroutine mxv 

