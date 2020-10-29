      subroutine quadr(f0, f1, f2, x1, x2, a, b, c) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:35:46  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: f0 
      real(double) , intent(in) :: f1 
      real(double) , intent(in) :: f2 
      real(double) , intent(in) :: x1 
      real(double) , intent(in) :: x2 
      real(double) , intent(out) :: a 
      real(double) , intent(out) :: b 
      real(double) , intent(out) :: c 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!***************************************************
!                                                  *
!    QUADR CALCULATES THE A, B AND C IN THE EQUNS. *
!                                                  *
!     A                   =   F0                   *
!     A + B.X0 + C.X0**2  =   F1                   *
!     A + B.X2 + C.X2**2  =   F2                   *
!                                                  *
!***************************************************
      c = (x2*(f1 - f0) - x1*(f2 - f0))/(x2*x1*(x1 - x2)) 
      b = (f1 - f0 - c*x1**2)/x1 
      a = f0 
      return  
      end subroutine quadr 
