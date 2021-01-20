      subroutine exchng(a, b, c, d, x, y, t, q, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:34:14  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: a 
      real(double) , intent(out) :: b 
      real(double) , intent(in) :: c 
      real(double) , intent(out) :: d 
      real(double) , intent(in) :: t 
      real(double) , intent(out) :: q 
      real(double) , intent(in) :: x(*) 
      real(double) , intent(out) :: y(*) 
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
