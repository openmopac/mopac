      MODULE deri2_I   
      INTERFACE
      subroutine deri2(minear, f, fd, fci, ninear, nvar_nvo, dxyzr,  throld, &
      diag, scalar, work) 
      USE vast_kind_param, ONLY:  double
      use molkst_C, only : mpack
      use meci_C, only : nmeci, mmci
      integer  :: minear 
      integer  :: ninear 
      integer  :: nvar_nvo
      real(double) , intent(in) :: throld  
      real(double)  :: f(minear,nvar_nvo) 
      real(double)  :: fd(ninear,nvar_nvo) 
      real(double)  :: fci(ninear,nmeci*(nmeci + 1)*10/ninear) 
      real(double)  :: work(*) 
      real(double), dimension(:,:), allocatable  :: b, ab, fb
      real(double) , intent(inout) :: dxyzr(nvar_nvo) 
      real(double)  :: diag(mpack) 
      real(double)  :: scalar(mpack) 
      real(double)  :: bab(mmci,mmci) 
      real(double)  :: babinv(mmci*mmci) 
      real(double)  :: bcoef(400) 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
