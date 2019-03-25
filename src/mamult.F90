      subroutine mamult(a, b, c, n, one) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:34:51  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: one 
      real(double) , intent(in) :: a(*) 
      real(double) , intent(in) :: b(*) 
      real(double) , intent(inout) :: c(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, i, ii, j, jj, k, kk 
      real(double) :: sum 
!-----------------------------------------------
!***********************************************************************
!
!   MAMULT MULTIPLIES A BY B AND PUTS THE RESULT IN C
!
!***********************************************************************
      l = 0 
      do i = 1, n 
        ii = ((i - 1)*i)/2 
        do j = 1, i 
          jj = ((j - 1)*j)/2 
          l = l + 1 
          sum = 0.D0 
          do k = 1, j 
            sum = sum + a(ii+k)*b(jj+k) 
          end do 
          do k = j + 1, i 
            sum = sum + a(ii+k)*b(((k-1)*k)/2+j) 
          end do 
          do k = i + 1, n 
            kk = (k*(k - 1))/2 
            sum = sum + a(kk+i)*b(kk+j) 
          end do 
          c(l) = sum + one*c(l) 
        end do 
      end do 
      return  
      end subroutine mamult 
