      subroutine rsp(a, n, root, vect) 
      USE vast_kind_param, ONLY:  double 
      implicit none
      integer  :: n 
      real(double)  :: a(*) 
      real(double), intent (out)  :: root(n) 
      real(double)  :: vect(n,n) 
!
! Trivial case: n = 1
!
      if (n == 1) then
        root(1) = a(1)
        vect(1,1) = 1.d0
        return
      end if
      call eigenvectors_LAPACK(vect, a, root, n)
      return  
      end subroutine rsp 

