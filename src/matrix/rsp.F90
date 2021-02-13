      subroutine rsp(a, n, root, vect) 
      implicit none
      integer  :: n 
      double precision  :: a(*) 
      double precision, intent (out)  :: root(n) 
      double precision  :: vect(n,n) 
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

