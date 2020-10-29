      subroutine sort(val, vec, n) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:27:36  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real , intent(inout) :: val(*) 
      complex , intent(inout) :: vec(n,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k 
      real :: x 
      complex :: sum 
!-----------------------------------------------
      do i = 1, n 
        x = 1.E9 
        do j = i, n 
          if (val(j) >= x) cycle  
          k = j 
          x = val(j) 
        end do 
        do j = 1, n 
          sum = vec(j,k) 
          vec(j,k) = vec(j,i) 
          vec(j,i) = sum 
        end do 
        val(k) = val(i) 
        val(i) = x 
      end do 
      return  
      end subroutine sort 
