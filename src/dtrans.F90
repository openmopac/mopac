      subroutine dtrans(d, ioper, first, r)  
!
!  Perform symmetry operation ioper on the "d" orbital set d 
!  
!  "r" holds the orientation matrix that rotates the system so that it is orientated 
!  along the principal axes
! 
      USE vast_kind_param, ONLY:  double 
      USE symmetry_C, only : nclass, elem 
      use dtran2_I 
      implicit none
      integer , intent(in) :: ioper 
      logical , intent(inout) :: first 
      real(double) , intent(inout) :: d(5) 
      real(double) , intent(in) :: r(3,3) 
      integer :: i, k 
      real(double), dimension(5,5,12) :: t1 
      real(double), dimension(3,3) :: s 
      double precision :: h(5)
!-----------------------------------------------
      if (first) then 
        first = .FALSE. 
!
!  Copy r into s because dtrans modifies the first argument
!
        s = r 
        t1 = 0.d0
        call dtran2 (s, t1, 1) 
        do k = 2, nclass 
          s = elem(:,:,k) 
          call dtran2 (s, t1, k) 
        end do 
      endif 
!
!  Matrix multiply d*t1*d(trans)
!
      do i = 1, 5 
        h(i) = 0.D0 
        h(i) = h(i) + sum(t1(i,:,1)*d) 
      end do 
      do i = 1, 5 
        d(i) = 0.D0 
        d(i) = d(i) + sum(t1(:,i,ioper)*h) 
      end do 
      return  
      end subroutine dtrans 
