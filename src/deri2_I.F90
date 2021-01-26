      MODULE deri2_I   
      INTERFACE
      subroutine deri2(minear, f, fd, fci, ninear, nvar_nvo, dxyzr,  throld, &
      diag, scalar, work) 
      use molkst_C, only : mpack
      use meci_C, only : nmeci, mmci
      integer  :: minear 
      integer  :: ninear 
      integer  :: nvar_nvo
      double precision , intent(in) :: throld  
      double precision  :: f(minear,nvar_nvo) 
      double precision  :: fd(ninear,nvar_nvo) 
      double precision  :: fci(ninear,nmeci*(nmeci + 1)*10/ninear) 
      double precision  :: work(*) 
      double precision, dimension(:,:), allocatable  :: b, ab, fb
      double precision , intent(inout) :: dxyzr(nvar_nvo) 
      double precision  :: diag(mpack) 
      double precision  :: scalar(mpack) 
      double precision  :: bab(mmci,mmci) 
      double precision  :: babinv(mmci*mmci) 
      double precision  :: bcoef(400) 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
