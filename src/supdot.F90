      subroutine supdot(s, h, g, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:27:50  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(inout) :: s(*) 
      real(double) , intent(in) :: h(*) 
      real(double) , intent(in) :: g(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j 
      real(double) :: sum, gi 
!-----------------------------------------------
!     (S)=(H)*(G) WITH  H  IN PACKED FORM (CANONICAL ORDER).
      k = 0 
      do i = 1, n 
        sum = 0.D0 
        do j = 1, i 
          sum = sum + g(j)*h(k+j) 
        end do 
        s(i) = sum 
        k = k + i 
      end do 
      if (n == 1) return  
      k = 1 
      do i = 2, n 
        gi = g(i) 
        s(:i-1) = s(:i-1) + h(k+1:i-1+k)*gi 
        k = k + i 
      end do 
      return  
      end subroutine supdot 
