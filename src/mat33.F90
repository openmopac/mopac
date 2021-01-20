      subroutine mat33(a, b, c) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!  A subroutine that will multiply two 3 by 3 matricies in the following
!     fashion:    C = A(transpose) B A
!
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:18:42  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: a(9) 
      real(double) , intent(in) :: b(9) 
      real(double) , intent(out) :: c(9) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(9) :: t 
!-----------------------------------------------
!
!
      t(1) = b(1)*a(1) + b(2)*a(4) + b(3)*a(7) 
      t(2) = b(1)*a(2) + b(2)*a(5) + b(3)*a(8) 
      t(3) = b(1)*a(3) + b(2)*a(6) + b(3)*a(9) 
      t(4) = b(4)*a(1) + b(5)*a(4) + b(6)*a(7) 
      t(5) = b(4)*a(2) + b(5)*a(5) + b(6)*a(8) 
      t(6) = b(4)*a(3) + b(5)*a(6) + b(6)*a(9) 
      t(7) = b(7)*a(1) + b(8)*a(4) + b(9)*a(7) 
      t(8) = b(7)*a(2) + b(8)*a(5) + b(9)*a(8) 
      t(9) = b(7)*a(3) + b(8)*a(6) + b(9)*a(9) 
!
      c(1) = a(1)*t(1) + a(4)*t(4) + a(7)*t(7) 
      c(2) = a(1)*t(2) + a(4)*t(5) + a(7)*t(8) 
      c(3) = a(1)*t(3) + a(4)*t(6) + a(7)*t(9) 
      c(4) = a(2)*t(1) + a(5)*t(4) + a(8)*t(7) 
      c(5) = a(2)*t(2) + a(5)*t(5) + a(8)*t(8) 
      c(6) = a(2)*t(3) + a(5)*t(6) + a(8)*t(9) 
      c(7) = a(3)*t(1) + a(6)*t(4) + a(9)*t(7) 
      c(8) = a(3)*t(2) + a(6)*t(5) + a(9)*t(8) 
      c(9) = a(3)*t(3) + a(6)*t(6) + a(9)*t(9) 
!
      return  
      end subroutine mat33 
