      subroutine supdot(s, h, g, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      double precision , intent(inout) :: s(*) 
      double precision , intent(in) :: h(*) 
      double precision , intent(in) :: g(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, i, j 
      double precision :: sum, gi 
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
