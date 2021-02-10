      subroutine deri0(e, n, scalar, diag, fract, nbo) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      double precision , intent(in) :: fract 
      integer , intent(in) :: nbo(3) 
      double precision , intent(in) :: e(n) 
      double precision , intent(out) :: scalar(*) 
      double precision , intent(inout) :: diag(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nopen, l, j, i 
      double precision :: shift, const 
!-----------------------------------------------
!
!     COMPUTE THE DIAGONAL DOMINANT PART OF THE SUPER-MATRIX AND
!     DEFINE THE SCALAR COEFFICIENTS APPLIED ON EACH ROW OF THE
!     SUPER LINEAR SYSTEM IN ORDER TO REDUCE THE EIGENVALUE SPECTRUM OF
!     THE ELECTRONIC HESSIAN,
!     THUS SPEEDING CONVERGENCE OF RELAXATION PROCESS IN 'DERI2'.
!  INPUT
!     E(N)             : EIGENVALUES OF FOCK MATRIX.
!     N                : NUMBER OF M.O.
!     NBO(3)           : OCCUPANCY BOUNDARIES.
!     FRACT            : PARTIAL OCCUPANCY OF 'OPEN' SHELLS.
!     SCALAR(MINEAR)   : SCALE APPLIED ON EACH COLUMN AND ROW OF THE
!                        SYMMETRIC SUPER SYSTEM.
!
      shift = 2.36D0 
!
!     DOMINANT DIAGONAL PART OF THE SUPER-MATRIX.
!     -------------------------------------------
      nopen = nbo(1) + nbo(2) 
      const = 1.D-3 
      l = 1 
      if (nbo(2)>0 .and. nbo(1)>0) then 
!        OPEN-CLOSED
        do j = 1, nbo(1) 
          if (nopen - nbo(1) > 0) then 
            diag(l:nopen-nbo(1)-1+l) = (e(nbo(1)+1:nopen)-e(j))/(2.D0 - fract&
               + const) 
            l = nopen - nbo(1) + l 
          end if 
        end do 
      end if 
      if (nbo(3)>0 .and. nbo(1)>0) then 
!        VIRTUAL-CLOSED
        do j = 1, nbo(1) 
          if (n - nopen > 0) then 
            diag(l:n-nopen-1+l) = (e(nopen+1:n)-e(j))/2.D0 
            l = n - nopen + l 
          end if 
        end do 
      end if 
      if (nbo(3)/=0 .and. nbo(2)/=0) then 
!        VIRTUAL-OPEN
        do j = nbo(1) + 1, nopen 
          if (n - nopen > 0) then 
            diag(l:n-nopen-1+l) = (e(nopen+1:n)-e(j))/(fract + const) 
            l = n - nopen + l 
          end if 
        end do 
      end if 
!
!     TAKE SCALE FACTORS AS (SHIFT-DIAG)**(-0.5) .
!     ------------------------------------------
      do i = 1, l - 1 
        scalar(i) = sqrt(1.D0/max(0.3D0*diag(i),diag(i)-shift)) 
      end do 
      return  
      end subroutine deri0 
