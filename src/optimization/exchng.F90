      subroutine exchng(a, b, c, d, x, y, t, q, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      double precision , intent(in) :: a 
      double precision , intent(out) :: b 
      double precision , intent(in) :: c 
      double precision , intent(out) :: d 
      double precision , intent(in) :: t 
      double precision , intent(out) :: q 
      double precision , intent(in) :: x(*) 
      double precision , intent(out) :: y(*) 
!-----------------------------------------------
!********************************************************************
!
! THE CONTENTS OF A, C, T, AND X ARE STORED IN B, D, Q, AND Y!
!
!   THIS IS A DEDICATED ROUTINE, IT IS CALLED BY LINMIN AND LOCMIN ONLY.
!
!********************************************************************
      b = a 
      d = c 
      q = t 
      y(:n) = x(:n) 
      return  
!
      end subroutine exchng 
