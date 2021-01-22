
            real(kind(0.0d0)) function ddot (n, dx, incx, dy, incy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(in) :: dx(*) 
      real(double) , intent(in) :: dy(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!           DOT = DX(I) * DY(I)
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      ddot = 0.0D+00 
      dtemp = 0.0D+00 
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1 
        iy = 1 
        if (incx < 0) ix = ((-n) + 1)*incx + 1 
        if (incy < 0) iy = ((-n) + 1)*incy + 1 
        dtemp = dot_product(dx(ix:(n-1)*incx+ix:incx),dy(iy:(n-1)*incy+iy:incy)) 
        ddot = dtemp 
        return  
      endif 
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,5) 
      if (m == 0) go to 40 
      dtemp = dot_product(dx(:m),dy(:m)) 
      if (n < 5) go to 60 
   40 continue 
      mp1 = m + 1 
      dtemp = dtemp + dot_product(dx(mp1:((n-mp1+5)/5)*5-1+mp1),dy(mp1:((n-mp1+5)/5)*5-1+mp1)) 
   60 continue 
      ddot = dtemp 
      return  
      end function ddot 

          subroutine dcopy(n, dx, incx, dy, incy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(in) :: dx(*) 
      real(double) , intent(out) :: dy(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1 
!-----------------------------------------------
!
!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1 
        iy = 1 
        if (incx < 0) ix = ((-n) + 1)*incx + 1 
        if (incy < 0) iy = ((-n) + 1)*incy + 1 
        dy(iy:(n-1)*incy+iy:incy) = dx(ix:(n-1)*incx+ix:incx) 
        return  
      endif 
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,7) 
      if (m /= 0) then 
        dy(:m) = dx(:m) 
        if (n < 7) return  
      endif 
      mp1 = m + 1 
      dy(mp1:((n-mp1+7)/7)*7-1+mp1) = dx(mp1:((n-mp1+7)/7)*7-1+mp1) 
      return  
      end subroutine dcopy 

      subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: k 
      integer , intent(in) :: lda 
      integer , intent(in) :: ldb 
      integer , intent(in) :: ldc 
      real(double) , intent(in) :: alpha 
      real(double) , intent(in) :: beta 
      character  :: transa*(*) 
      character  :: transb*(*) 
      real(double) , intent(in) :: a(lda,k) 
      real(double) , intent(in) :: b(ldb,k + n) ! "k + n" is here so that the debugger can be used
      real(double) , intent(out) :: c(ldc,n) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, info, j, l, nrowa, nrowb 
      real(double) :: temp 
      logical :: nota, notb 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
!
!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
!
!  WHERE  OP( X ) IS ONE OF
!
!     OP( X ) = X   OR   OP( X ) = X',
!
!  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
!  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
!
!  PARAMETERS
!  ==========
!
!  TRANSA - CHARACTER*1.
!           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
!           THE MATRIX MULTIPLICATION AS FOLLOWS:
!
!              TRANSA = 'N' OR 'N',  OP( A ) = A.
!
!              TRANSA = 'T' OR 'T',  OP( A ) = A'.
!
!              TRANSA = 'C' OR 'C',  OP( A ) = A'.
!
!           UNCHANGED ON EXIT.
!
!  TRANSB - CHARACTER*1.
!           ON ENTRY, TRANSB SPECIFIES THE FORM OF OP( B ) TO BE USED IN
!           THE MATRIX MULTIPLICATION AS FOLLOWS:
!
!              TRANSB = 'N' OR 'N',  OP( B ) = B.
!
!              TRANSB = 'T' OR 'T',  OP( B ) = B'.
!
!              TRANSB = 'C' OR 'C',  OP( B ) = B'.
!
!           UNCHANGED ON EXIT.
!
!  M      - INTEGER.
!           ON ENTRY,  M  SPECIFIES  THE NUMBER  OF ROWS  OF THE  MATRIX
!           OP( A )  AND OF THE  MATRIX  C.  M  MUST  BE AT LEAST  ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY,  N  SPECIFIES THE NUMBER  OF COLUMNS OF THE MATRIX
!           OP( B ) AND THE NUMBER OF COLUMNS OF THE MATRIX C. N MUST BE
!           AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  K      - INTEGER.
!           ON ENTRY,  K  SPECIFIES  THE NUMBER OF COLUMNS OF THE MATRIX
!           OP( A ) AND THE NUMBER OF ROWS OF THE MATRIX OP( B ). K MUST
!           BE AT LEAST  ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, KA ), WHERE KA IS
!           K  WHEN  TRANSA = 'N' OR 'N',  AND IS  M  OTHERWISE.
!           BEFORE ENTRY WITH  TRANSA = 'N' OR 'N',  THE LEADING  M BY K
!           PART OF THE ARRAY  A  MUST CONTAIN THE MATRIX  A,  OTHERWISE
!           THE LEADING  K BY M  PART OF THE ARRAY  A  MUST CONTAIN  THE
!           MATRIX A.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSA = 'N' OR 'N' THEN
!           LDA MUST BE AT LEAST  MAX( 1, M ), OTHERWISE  LDA MUST BE AT
!           LEAST  MAX( 1, K ).
!           UNCHANGED ON EXIT.
!
!  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, KB ), WHERE KB IS
!           N  WHEN  TRANSB = 'N' OR 'N',  AND IS  K  OTHERWISE.
!           BEFORE ENTRY WITH  TRANSB = 'N' OR 'N',  THE LEADING  K BY N
!           PART OF THE ARRAY  B  MUST CONTAIN THE MATRIX  B,  OTHERWISE
!           THE LEADING  N BY K  PART OF THE ARRAY  B  MUST CONTAIN  THE
!           MATRIX B.
!           UNCHANGED ON EXIT.
!
!  LDB    - INTEGER.
!           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. WHEN  TRANSB = 'N' OR 'N' THEN
!           LDB MUST BE AT LEAST  MAX( 1, K ), OTHERWISE  LDB MUST BE AT
!           LEAST  MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  BETA   - DOUBLE PRECISION.
!           ON ENTRY,  BETA  SPECIFIES THE SCALAR  BETA.  WHEN  BETA  IS
!           SUPPLIED AS ZERO THEN C NEED NOT BE SET ON INPUT.
!           UNCHANGED ON EXIT.
!
!  C      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDC, N ).
!           BEFORE ENTRY, THE LEADING  M BY N  PART OF THE ARRAY  C MUST
!           CONTAIN THE MATRIX  C,  EXCEPT WHEN  BETA  IS ZERO, IN WHICH
!           CASE C NEED NOT BE SET ON ENTRY.
!           ON EXIT, THE ARRAY  C  IS OVERWRITTEN BY THE  M BY N  MATRIX
!           ( ALPHA*OP( A )*OP( B ) + BETA*C ).
!
!  LDC    - INTEGER.
!           ON ENTRY, LDC SPECIFIES THE FIRST DIMENSION OF C AS DECLARED
!           IN  THE  CALLING  (SUB)  PROGRAM.   LDC  MUST  BE  AT  LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 3 BLAS ROUTINE.
!
!  -- WRITTEN ON 8-FEBRUARY-1989.
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!     IAIN DUFF, AERE HARWELL.
!     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
!     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
!
!
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
!     .. LOCAL SCALARS ..
!     .. PARAMETERS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     SET  NOTA  AND  NOTB  AS  TRUE IF  A  AND  B  RESPECTIVELY ARE NOT
!     TRANSPOSED AND SET  NROWA, NCOLA AND  NROWB  AS THE NUMBER OF ROWS
!     AND  COLUMNS OF  A  AND THE  NUMBER OF  ROWS  OF  B  RESPECTIVELY.
!
      nota = lsame(transa,'N') 
      notb = lsame(transb,'N') 
      if (nota) then 
        nrowa = m 
      else 
        nrowa = k 
      endif 
      if (notb) then 
        nrowb = k 
      else 
        nrowb = n 
      endif 
!
!     TEST THE INPUT PARAMETERS.
!
      info = 0 
      if (.not.nota .and. .not.lsame(transa,'C') .and. .not.lsame(transa,'T')) &
        then 
        info = 1 
      else if (.not.notb .and. .not.lsame(transb,'C') .and. .not.lsame(transb,&
          'T')) then 
        info = 2 
      else if (m < 0) then 
        info = 3 
      else if (n < 0) then 
        info = 4 
      else if (k < 0) then 
        info = 5 
      else if (lda < max(1,nrowa)) then 
        info = 8 
      else if (ldb < max(1,nrowb)) then 
        info = 10 
      else if (ldc < max(1,m)) then 
        info = 13 
      endif 
      if (info /= 0) then 
        call xerbla ('DGEMM ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (m==0 .or. n==0 .or. (Abs(alpha) < 1.d-20 .or. k==0) .and. Abs(beta - 1.d0) < 1.d-20) return  
!
!     AND IF  ALPHA.EQ.ZERO.
!
      if (Abs(alpha) < 1.d-20) then 
        if (Abs(beta - 1.d0) < 1.d-20) then 
          c(:m,:n) = zero 
        else 
          c(:m,:n) = beta*c(:m,:n) 
        endif 
        return  
      endif 
!
!     START THE OPERATIONS.
!
      if (notb) then 
        if (nota) then 
!
!           FORM  C := ALPHA*A*B + BETA*C.
!
          do j = 1, n 
            if (beta == zero) then 
              c(:m,j) = zero 
            else if (Abs(beta - 1.d0) > 1.d-20) then 
              c(:m,j) = beta*c(:m,j) 
            endif 
            do l = 1, k 
              if (b(l,j) == zero) cycle  
              temp = alpha*b(l,j) 
              c(:m,j) = c(:m,j) + temp*a(:m,l) 
            end do 
          end do 
        else 
!
!           FORM  C := ALPHA*A'*B + BETA*C
!
          do j = 1, n 
            do i = 1, m 
              temp = zero 
              temp = temp + sum(a(:k,i)*b(:k,j)) 
              if (beta == zero) then 
                c(i,j) = alpha*temp 
              else 
                c(i,j) = alpha*temp + beta*c(i,j) 
              endif 
            end do 
          end do 
        endif 
      else 
        if (nota) then 
!
!           FORM  C := ALPHA*A*B' + BETA*C
!
          do j = 1, n 
            if (beta == zero) then 
              c(:m,j) = zero 
            else if (Abs(beta-1.d0) > 1.d-20) then 
              c(:m,j) = beta*c(:m,j) 
            endif 
            do l = 1, k 
              if (b(j,l) == zero) cycle  
              temp = alpha*b(j,l) 
              c(:m,j) = c(:m,j) + temp*a(:m,l) 
            end do 
          end do 
        else 
!
!           FORM  C := ALPHA*A'*B' + BETA*C
!
          do j = 1, n 
            do i = 1, m 
              temp = zero 
              temp = temp + sum(a(:k,i)*b(j,:k)) 
              if (beta == zero) then 
                c(i,j) = alpha*temp 
              else 
                c(i,j) = alpha*temp + beta*c(i,j) 
              endif 
            end do 
          end do 
        endif 
      endif 
!
      return  
!
!     END OF DGEMM .
!
      end subroutine dgemm 

   subroutine dscal(n, da, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      real(double) , intent(in) :: da 
      real(double) , intent(inout) :: dx(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m, mp1, nincx 
!-----------------------------------------------
!
!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx /= 1) then 
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
        nincx = n*incx 
        dx(:nincx:incx) = da*dx(:nincx:incx) 
        return  
      endif 
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,5) 
      if (m /= 0) then 
        dx(:m) = da*dx(:m) 
        if (n < 5) return  
      endif 
      mp1 = m + 1 
      dx(mp1:((n-mp1+5)/5)*5-1+mp1) = da*dx(mp1:((n-mp1+5)/5)*5-1+mp1) 
      return  
      end subroutine dscal 

           subroutine daxpy(n, da, dx, incx, dy, incy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(in) :: da 
      real(double) , intent(in) :: dx(*) 
      real(double) , intent(inout) :: dy(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix, iy, m, mp1 
!-----------------------------------------------
!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!           DY(I) = DY(I) + DA * DX(I)
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      if (n <= 0) return  
      if (da == 0.0D+00) return  
      if (incx/=1 .or. incy/=1) then 
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
        ix = 1 
        iy = 1 
        if (incx < 0) ix = ((-n) + 1)*incx + 1 
        if (incy < 0) iy = ((-n) + 1)*incy + 1 
        dy(iy:(n-1)*incy+iy:incy) = dy(iy:(n-1)*incy+iy:incy) + da*dx(ix:(n-1)*incx+ix:incx) 
        return  
      endif 
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,4) 
      if (m /= 0) then 
        dy(:m) = dy(:m) + da*dx(:m) 
        if (n < 4) return  
      endif 
      mp1 = m + 1 
      dy(mp1:((n-mp1+4)/4)*4-1+mp1) = dy(mp1:((n-mp1+4)/4)*4-1+mp1) + da*dx(mp1:((n-mp1+4)/4)*4-1+mp1) 
      return  
      end subroutine daxpy 

                subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   !   use dgetrf_I 
      use dgetrs_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: nrhs 
      integer  :: lda 
      integer  :: ldb 
      integer  :: info 
      integer  :: ipiv(*) 
      real(double)  :: a(lda,*) 
      real(double)  :: b(ldb,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DGESV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!
!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor A as
!     A = P * L * U,
!  where P is a permutation matrix, L is unit lower triangular, and U is
!  upper triangular.  The factored form of A is then used to solve the
!  system of equations A * X = B.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices that define the permutation matrix P;
!          row i of the matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                has been completed, but the factor U is exactly
!                singular, so the solution could not be computed.
!
!  =====================================================================
!
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      if (n < 0) then 
        info = -1 
      else if (nrhs < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      else if (ldb < max(1,n)) then 
        info = -7 
      endif 
      if (info /= 0) then 
        call xerbla ('DGESV ', (-info)) 
        return  
      endif 
!
!     Compute the LU factorization of A.
!
      call dgetrf (n, n, a, lda, ipiv, info) 
!
!        Solve the system A*X = B, overwriting B with X.
!
      if (info == 0) call dgetrs ('No transpose', n, nrhs, a, lda, ipiv, b, ldb&
        , info) 
      return  
!
!     End of DGESV
!
      end subroutine dgesv


      subroutine dswap(n, dx, incx, dy, incy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(inout) :: dx(*) 
      real(double) , intent(inout) :: dy(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix, iy, m, mp1 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      if (n <= 0) return  
      if (incx/=1 .or. incy/=1) then 
!
!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1
!
        ix = 1 
        iy = 1 
        if (incx < 0) ix = ((-n) + 1)*incx + 1 
        if (incy < 0) iy = ((-n) + 1)*incy + 1 
        do i = 1, n 
          dtemp = dx(ix) 
          dx(ix) = dy(iy) 
          dy(iy) = dtemp 
          ix = ix + incx 
          iy = iy + incy 
        end do 
        return  
      endif 
!
!       CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!       CLEAN-UP LOOP
!
      m = mod(n,3) 
      if (m /= 0) then 
        do i = 1, m 
          dtemp = dx(i) 
          dx(i) = dy(i) 
          dy(i) = dtemp 
        end do 
        if (n < 3) return  
      endif 
      mp1 = m + 1 
      do i = mp1, n, 3 
        dtemp = dx(i) 
        dx(i) = dy(i) 
        dy(i) = dtemp 
        dtemp = dx(i+1) 
        dx(i+1) = dy(i+1) 
        dy(i+1) = dtemp 
        dtemp = dx(i+2) 
        dx(i+2) = dy(i+2) 
        dy(i+2) = dtemp 
      end do 
      return  
      end subroutine dswap 


          subroutine dgetrf(m, n, a, lda, ipiv, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dgetf2_I 
      use dlaswp_I 
      use dtrsm_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: m 
      integer  :: n 
      integer  :: lda 
      integer  :: info 
      integer  :: ipiv(*) 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jb, nb 
      integer, external :: ilaenv
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, min 
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      if (m < 0) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,m)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DGETRF', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (m==0 .or. n==0) return  
!
!     Determine the block size for this environment.
!
      nb = ilaenv(1,'DGETRF'," ",m,n,-1, -1) 
      if (nb<=1 .or. nb>=min(m,n)) then 
!
!        Use unblocked code.
!
        call dgetf2 (m, n, a, lda, ipiv) 
      else 
!
!        Use blocked code.
!
        do j = 1, min(m,n), nb 
          jb = min(min(m,n) - j + 1,nb) 
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
          call dgetf2 (m - j + 1, jb, a(j,j), lda, ipiv(j)) 
!
!           Adjust INFO and the pivot indices.
!
      !    if (info==0 .and. iinfo>0) info = iinfo + j - 1 
          ipiv(j:min(m,j+jb-1)) = j - 1 + ipiv(j:min(m,j+jb-1)) 
!
!           Apply interchanges to columns 1:J-1.
!
          call dlaswp (j - 1, a, lda, j, j + jb - 1, ipiv, 1) 
!
          if (j + jb > n) cycle  
!
!              Apply interchanges to columns J+JB:N.
!
          call dlaswp (n - j - jb + 1, a(1,j+jb), lda, j, j + jb - 1, ipiv, 1) 
!
!              Compute block row of U.
!
          call dtrsm ('Left', 'Lower', 'No transpose', 'Unit', jb, n - j - jb&
             + 1, one, a(j,j), lda, a(j,j+jb), lda) 
          if (j + jb > m) cycle  
!
!                 Update trailing submatrix.
!
          call dgemm ('No transpose', 'No transpose', m - j - jb + 1, n - j - &
            jb + 1, jb, (-one), a(j+jb,j), lda, a(j,j+jb), lda, one, a(j+jb,j+&
            jb), lda) 
        end do 
      endif 
      return  
!
!     End of DGETRF
!
      end subroutine dgetrf 


      subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use dlaswp_I 
      use dtrsm_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: nrhs 
      integer  :: lda 
      integer  :: ldb 
      integer , intent(out) :: info 
      character  :: trans 
      integer  :: ipiv(*) 
      real(double)  :: a(lda,*) 
      real(double)  :: b(ldb,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: notran 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      notran = lsame(trans,'N') 
      if (.not.notran .and. .not.lsame(trans,'T') .and. .not.lsame(trans,'C')) &
        then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (nrhs < 0) then 
        info = -3 
      else if (lda < max(1,n)) then 
        info = -5 
      else if (ldb < max(1,n)) then 
        info = -8 
      endif 
      if (info /= 0) then 
        call xerbla ('DGETRS', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n==0 .or. nrhs==0) return  
!
      if (notran) then 
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
        call dlaswp (nrhs, b, ldb, 1, n, ipiv, 1) 
!
!        Solve L*X = B, overwriting B with X.
!
        call dtrsm ('Left', 'Lower', 'No transpose', 'Unit', n, nrhs, one, a, &
          lda, b, ldb) 
!
!        Solve U*X = B, overwriting B with X.
!
        call dtrsm ('Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs, one, &
          a, lda, b, ldb) 
      else 
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
        call dtrsm ('Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs, one, a, &
          lda, b, ldb) 
!
!        Solve L'*X = B, overwriting B with X.
!
        call dtrsm ('Left', 'Lower', 'Transpose', 'Unit', n, nrhs, one, a, lda&
          , b, ldb) 
!
!        Apply row interchanges to the solution vectors.
!
        call dlaswp (nrhs, b, ldb, 1, n, ipiv, -1) 
      endif 
!
      return  
!
!     End of DGETRS
!
      end subroutine dgetrs 

         real(kind(0.0d0)) function dasum (n, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      real(double) , intent(in) :: dx(incx*n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, m, mp1, nincx 
      real(double) :: dtemp 
!-----------------------------------------------
!
!     TAKES THE SUM OF THE ABSOLUTE VALUES.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!
      dasum = 0.0D+00 
      dtemp = 0.0D+00 
      if (n <= 0) return  
      if (incx /= 1) then 
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
        nincx = n*incx 
        do i = 1, nincx, incx 
          dtemp = dtemp + abs(dx(i)) 
        end do 
        dasum = dtemp 
        return  
      endif 
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
      m = mod(n,6) 
      if (m == 0) go to 40 
      do i = 1, m 
        dtemp = dtemp + abs(dx(i)) 
      end do 
      if (n < 6) go to 60 
   40 continue 
      mp1 = m + 1 
      do i = mp1, n, 6 
        dtemp = dtemp + abs(dx(i)) + abs(dx(i+1)) + abs(dx(i+2)) + abs(dx(i+3)) &
           + abs(dx(i+4)) + abs(dx(i+5)) 
      end do 
   60 continue 
      dasum = dtemp 
      return  
      end function dasum 

    real(kind(0.0d0)) function dnrm2 (n, dx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use dot_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: dx(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: norm
      norm = dot(dx,dx,n)
      if (norm > 0.d0) then
        norm = sqrt(norm)
      else
        norm=1.d-30
      end if
      dnrm2 = norm
      return  
      end function dnrm2 

   subroutine dpotri(uplo, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: lda 
      integer  :: info 
      character  :: uplo 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
      logical, external :: lsame
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DPOTRI computes the inverse of a real symmetric positive definite
!  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
!  computed by DPOTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T, as computed by
!          DPOTRF.
!          On exit, the upper or lower triangle of the (symmetric)
!          inverse of A, overwriting the input factor U or L.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the (i,i) element of the factor U or L is
!                zero, and the inverse could not be computed.
!
!  =====================================================================
!
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      if (.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DPOTRI', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
!     Invert the triangular Cholesky factor U or L.
!
      call dtrtri (uplo, 'Non-unit', n, a, lda, info) 
      if (info > 0) return  
!
!     Form inv(U)*inv(U)' or inv(L)'*inv(L).
!
      call dlauum (uplo, n, a, lda, info) 
!
      return  
!
!     End of DPOTRI
!
      end subroutine dpotri 

   integer function ilaenv (ispec, name, opts, n1, n2, n3, n4) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ieeeck_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ispec 
      integer , intent(in) :: n1 
      integer , intent(in) :: n2 
      integer  :: n3 
      integer , intent(in) :: n4 
      character , intent(in) :: name*(*) 
      character  :: opts*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ic, iz, nb, nbmin, nx 
      logical :: cname, sname 
      character :: c1, c2*2, c4*2, c3*3, subnam*6 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC char, ichar, int, min, real 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      if (n3 == -10000 .and. opts == " ") return ! dummy use of n3 and opts
      select case (ispec)  
!
!     Invalid value for ISPEC
!
      case default 
        ilaenv = -1 
        return  
!
      case (1:3)  
        ilaenv = 1 
        subnam = name 
        ic = ichar(subnam(1:1)) 
        iz = ichar('Z') 
        if (iz==90 .or. iz==122) then 
!
!        ASCII character set
!
          if (ic>=97 .and. ic<=122) then 
            subnam(1:1) = char(ic - 32) 
            do i = 2, 6 
              ic = ichar(subnam(i:i)) 
              if (ic<97 .or. ic>122) cycle  
              subnam(i:i) = char(ic - 32) 
            end do 
          endif 
!
        else if (iz==233 .or. iz==169) then 
!
!        EBCDIC character set
!
          if (ic>=129 .and. ic<=137 .or. ic>=145 .and. ic<=153 .or. ic>=162&
             .and. ic<=169) then 
            subnam(1:1) = char(ic + 64) 
            do i = 2, 6 
              ic = ichar(subnam(i:i)) 
              if (.not.(ic>=129 .and. ic<=137 .or. ic>=145 .and. ic<=153 .or. &
                ic>=162 .and. ic<=169)) cycle  
              subnam(i:i) = char(ic + 64) 
            end do 
          endif 
!
        else if (iz==218 .or. iz==250) then 
!
!        Prime machines:  ASCII+128
!
          if (ic>=225 .and. ic<=250) then 
            subnam(1:1) = char(ic - 32) 
            do i = 2, 6 
              ic = ichar(subnam(i:i)) 
              if (ic<225 .or. ic>250) cycle  
              subnam(i:i) = char(ic - 32) 
            end do 
          endif 
        endif 
!
        c1 = subnam(1:1) 
        sname = c1=='S' .or. c1=='D' 
        cname = c1=='C' .or. c1=='Z' 
        if (.not.(cname .or. sname)) return  
        c2 = subnam(2:3) 
        c3 = subnam(4:6) 
        c4 = c3(2:3) 
!
        select case (ispec)  
!
        case default 
          nb = 1 
!
          if (c2 == 'GE') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            else if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') &
                then 
              if (sname) then 
                nb = 32 
              else 
                nb = 32 
              endif 
            else if (c3 == 'HRD') then 
              if (sname) then 
                nb = 32 
              else 
                nb = 32 
              endif 
            else if (c3 == 'BRD') then 
              if (sname) then 
                nb = 32 
              else 
                nb = 32 
              endif 
            else if (c3 == 'TRI') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            endif 
          else if (c2 == 'PO') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            endif 
          else if (c2 == 'SY') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            else if (sname .and. c3=='TRD') then 
              nb = 32 
            else if (sname .and. c3=='GST') then 
              nb = 64 
            endif 
          else if (cname .and. c2=='HE') then 
            select case (c3)  
            case ('TRF')  
              nb = 64 
            case ('TRD')  
              nb = 32 
            case ('GST')  
              nb = 64 
            end select 
          else if (sname .and. c2=='OR') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nb = 32 
            else if (c3(1:1) == 'M') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nb = 32 
            endif 
          else if (cname .and. c2=='UN') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nb = 32 
            else if (c3(1:1) == 'M') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nb = 32 
            endif 
          else if (c2 == 'GB') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                if (n4 <= 64) then 
                  nb = 1 
                else 
                  nb = 32 
                endif 
              else 
                if (n4 <= 64) then 
                  nb = 1 
                else 
                  nb = 32 
                endif 
              endif 
            endif 
          else if (c2 == 'PB') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                if (n2 <= 64) then 
                  nb = 1 
                else 
                  nb = 32 
                endif 
              else 
                if (n2 <= 64) then 
                  nb = 1 
                else 
                  nb = 32 
                endif 
              endif 
            endif 
          else if (c2 == 'TR') then 
            if (c3 == 'TRI') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            endif 
          else if (c2 == 'LA') then 
            if (c3 == 'UUM') then 
              if (sname) then 
                nb = 64 
              else 
                nb = 64 
              endif 
            endif 
          else if (sname .and. c2=='ST') then 
            if (c3 == 'EBZ') nb = 1 
          endif 
          ilaenv = nb 
          return  
!
        case (2)  
          nbmin = 2 
          if (c2 == 'GE') then 
            if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') then 
              if (sname) then 
                nbmin = 2 
              else 
                nbmin = 2 
              endif 
            else if (c3 == 'HRD') then 
              if (sname) then 
                nbmin = 2 
              else 
                nbmin = 2 
              endif 
            else if (c3 == 'BRD') then 
              if (sname) then 
                nbmin = 2 
              else 
                nbmin = 2 
              endif 
            else if (c3 == 'TRI') then 
              if (sname) then 
                nbmin = 2 
              else 
                nbmin = 2 
              endif 
            endif 
          else if (c2 == 'SY') then 
            if (c3 == 'TRF') then 
              if (sname) then 
                nbmin = 8 
              else 
                nbmin = 8 
              endif 
            else if (sname .and. c3=='TRD') then 
              nbmin = 2 
            endif 
          else if (cname .and. c2=='HE') then 
            if (c3 == 'TRD') nbmin = 2 
          else if (sname .and. c2=='OR') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nbmin = 2 
            else if (c3(1:1) == 'M') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nbmin = 2 
            endif 
          else if (cname .and. c2=='UN') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nbmin = 2 
            else if (c3(1:1) == 'M') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nbmin = 2 
            endif 
          endif 
          ilaenv = nbmin 
          return  
!
        case (3)  
          nx = 0 
          if (c2 == 'GE') then 
            if (c3=='QRF' .or. c3=='RQF' .or. c3=='LQF' .or. c3=='QLF') then 
              if (sname) then 
                nx = 128 
              else 
                nx = 128 
              endif 
            else if (c3 == 'HRD') then 
              if (sname) then 
                nx = 128 
              else 
                nx = 128 
              endif 
            else if (c3 == 'BRD') then 
              if (sname) then 
                nx = 128 
              else 
                nx = 128 
              endif 
            endif 
          else if (c2 == 'SY') then 
            if (sname .and. c3=='TRD') nx = 32 
          else if (cname .and. c2=='HE') then 
            if (c3 == 'TRD') nx = 32 
          else if (sname .and. c2=='OR') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nx = 128 
            endif 
          else if (cname .and. c2=='UN') then 
            if (c3(1:1) == 'G') then 
              if (c4=='QR' .or. c4=='RQ' .or. c4=='LQ' .or. c4=='QL' .or. c4==&
                'HR' .or. c4=='TR' .or. c4=='BR') nx = 128 
            endif 
          endif 
          ilaenv = nx 
          return  
!
        case (4)  
          ilaenv = 6 
          return  
!
        case (5)  
          ilaenv = 2 
          return  
!
        case (6)  
          ilaenv = int(real(min(n1,n2))*1.6E0) 
          return  
!
        case (7)  
          ilaenv = 1 
          return  
!
        case (8)  
          ilaenv = 50 
          return  
!
        case (9)  
          ilaenv = 25 
          return  
!
        case (10)  
          ilaenv = 1 
          if (ilaenv == 1) ilaenv = ieeeck(0,0.0,1.0) 
          return  
!
        case (11)  
          ilaenv = 1 
          if (ilaenv == 1) ilaenv = ieeeck(1,0.0,1.0) 
          return  
        end select 
      end select 
!
!     End of ILAENV
!
      end function ilaenv 


 
      subroutine dpotrf(uplo, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use lsame_I 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: lda 
      integer  :: info 
      character  :: uplo 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jb, nb 
      logical :: upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, min 
      integer, external :: ilaenv
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DPOTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the block version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DPOTRF', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
!     Determine the block size for this environment.
!
      nb = ilaenv(1,'DPOTRF',uplo,n,-1,-1,-1) 
      if (nb<=1 .or. nb>=n) then 
!
!        Use unblocked code.
!
        call dpotf2 (uplo, n, a, lda, info) 
      else 
!
!        Use blocked code.
!
        if (upper) then 
!
!           Compute the Cholesky factorization A = U'*U.
!
          do j = 1, n, nb 
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            jb = min(nb,n - j + 1) 
            call dsyrk ('Upper', 'Transpose', jb, j - 1, (-one), a(1,j), lda, &
              one, a(j,j), lda) 
            call dpotf2 ('Upper', jb, a(j,j), lda, info) 
            if (info /= 0) go to 30 
            if (j + jb > n) cycle  
!
!                 Compute the current block row.
!
            call dgemm ('Transpose', 'No transpose', jb, n - j - jb + 1, j - 1&
              , (-one), a(1,j), lda, a(1,j+jb), lda, one, a(j,j+jb), lda) 
            call dtrsm ('Left', 'Upper', 'Transpose', 'Non-unit', jb, n - j - &
              jb + 1, one, a(j,j), lda, a(j,j+jb), lda) 
          end do 
!
        else 
!
!           Compute the Cholesky factorization A = L*L'.
!
          do j = 1, n, nb 
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            jb = min(nb,n - j + 1) 
            call dsyrk ('Lower', 'No transpose', jb, j - 1, (-one), a(j,1), lda&
              , one, a(j,j), lda) 
            call dpotf2 ('Lower', jb, a(j,j), lda, info) 
            if (info /= 0) go to 30 
            if (j + jb > n) cycle  
!
!                 Compute the current block column.
!
            call dgemm ('No transpose', 'Transpose', n - j - jb + 1, jb, j - 1&
              , (-one), a(j+jb,1), lda, a(j,1), lda, one, a(j+jb,j), lda) 
            call dtrsm ('Right', 'Lower', 'Transpose', 'Non-unit', n - j - jb&
               + 1, jb, one, a(j,j), lda, a(j+jb,j), lda) 
          end do 
        endif 
      endif 
      go to 40 
!
   30 continue 
      info = info + j - 1 
!
   40 continue 
      return  
!
!     End of DPOTRF
!
      end subroutine dpotrf 


 
      subroutine xerbla(srname, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
     USE chanel_C, only : iw
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: info 
      character , intent(in) :: srname*6 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!  -- LAPACK AUXILIARY ROUTINE (VERSION 1.0B) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992
!
!     .. SCALAR ARGUMENTS ..
! MOPAC CHANGE
      write (iw, fmt=10) srname, info 
! END OF MOPAC CHANGE
!
      call mopend ('Error in BLAS') 
      return  
!
   10 format(' ** ON ENTRY TO ',a6,' PARAMETER NUMBER ',i2,' HAD ',&
        'AN ILLEGAL VALUE') 
!
!     END OF XERBLA
!
      end subroutine xerbla 

      
      subroutine dgetf2(m, n, a, lda, ipiv) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use idamax_I 
      use dger_I 
      use dscal_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer  :: n 
      integer  :: lda 
      integer  :: info   ! WARNING
      integer , intent(out) :: ipiv(*) 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jp 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, min 
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      if (m < 0) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,m)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DGETF2', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (m==0 .or. n==0) return  
!
      do j = 1, min(m,n) 
!
!        Find pivot and test for singularity.
!
        jp = j - 1 + idamax(m - j + 1,a(j,j),1) 
        ipiv(j) = jp 
        if (a(jp,j) /= zero) then 
!
!           Apply the interchange to columns 1:N.
!
          if (jp /= j) call dswap (n, a(j,1), lda, a(jp,1), lda) 
!
!           Compute elements J+1:M of J-th column.
!
          if (j < m) call dscal (m - j, one/a(j,j), a(j+1,j), 1) 
!
        else if (info == 0) then 
!
          info = j 
        endif 
!
        if (j >= min(m,n)) cycle  
!
!           Update trailing submatrix.
!
        call dger (m - j, n - j, (-one), a(j+1,j), 1, a(j,j+1), lda, a(j+1,j+1)&
          , lda) 
      end do 
      return  
!
!     End of DGETF2
!
      end subroutine dgetf2 


      subroutine dlaswp(n, a, lda, k1, k2, ipiv, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: lda 
      integer , intent(in) :: k1 
      integer , intent(in) :: k2 
      integer , intent(in) :: incx 
      integer , intent(in) :: ipiv(*) 
      real(double) , intent(inout) :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, i1, i2, inc, ip, ix, ix0, j, k, n32 
      real(double) :: temp 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      if (incx > 0) then 
        ix0 = k1 
        i1 = k1 
        i2 = k2 
        inc = 1 
      else if (incx < 0) then 
        ix0 = 1 + (1 - k2)*incx 
        i1 = k2 
        i2 = k1 
        inc = -1 
      else 
        return  
      endif 
!
      n32 = (n/32)*32 
      if (n32 /= 0) then 
        do j = 1, n32, 32 
          ix = ix0 
          do i = i1, i2, inc 
            ip = ipiv(ix) 
            if (ip /= i) then 
              do k = j, j + 31 
                temp = a(i,k) 
                a(i,k) = a(ip,k) 
                a(ip,k) = temp 
              end do 
            endif 
            ix = ix + incx 
          end do 
        end do 
      endif 
      if (n32 /= n) then 
        n32 = n32 + 1 
        ix = ix0 
        do i = i1, i2, inc 
          ip = ipiv(ix) 
          if (ip /= i) then 
            do k = n32, n 
              temp = a(i,k) 
              a(i,k) = a(ip,k) 
              a(ip,k) = temp 
            end do 
          endif 
          ix = ix + incx 
        end do 
      endif 
!
      return  
!
!     End of DLASWP
!
      end subroutine dlaswp 



      integer function ieeeck (ispec, zero, one) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:54:31  01/16/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ispec 
      real , intent(in) :: zero 
      real , intent(in) :: one 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real ::  neginf, negzro, newzro, &
        posinf 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
      ieeeck = 1 
!
      posinf = one/zero 
      if (posinf <= one) then 
        ieeeck = 0 
        return  
      endif 
!
      neginf = -one/zero 
      if (neginf >= zero) then 
        ieeeck = 0 
        return  
      endif 
!
      negzro = one/(neginf + one) 
      if (Abs(negzro) > 1.d-20) then 
        ieeeck = 0 
        return  
      endif 
!
      neginf = one/negzro 
      if (neginf > -1.d-20) then 
        ieeeck = 0 
        return  
      endif 
!
      newzro = negzro + zero 
      if (Abs(newzro) > 1.d-20) then 
        ieeeck = 0 
        return  
      endif 
!
      posinf = one/newzro 
      if (posinf <= one) then 
        ieeeck = 0 
        return  
      endif 
!
      neginf = neginf*posinf 
      if (neginf >= zero) then 
        ieeeck = 0 
        return  
      endif 
!
      posinf = posinf*posinf 
      if (posinf <= one) then 
        ieeeck = 0 
        return  
      endif 
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      if (ispec == 0) return  
!
        ieeeck = 0  ! Bypass all the nan tests
!
      return  
      end function ieeeck 



      subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: lda 
      integer , intent(in) :: ldb 
      real(double) , intent(in) :: alpha 
      character  :: side*(*) 
      character  :: uplo*(*) 
      character  :: transa*(*) 
      character  :: diag*(*) 
      real(double) , intent(in) :: a(lda,*) 
      real(double) , intent(inout) :: b(ldb,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, info, j, k, nrowa 
      real(double) :: temp 
      logical :: lside, nounit, upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DTRSM  SOLVES ONE OF THE MATRIX EQUATIONS
!
!     OP( A )*X = ALPHA*B,   OR   X*OP( A ) = ALPHA*B,
!
!  WHERE ALPHA IS A SCALAR, X AND B ARE M BY N MATRICES, A IS A UNIT, OR
!  NON-UNIT,  UPPER OR LOWER TRIANGULAR MATRIX  AND  OP( A )  IS ONE  OF
!
!     OP( A ) = A   OR   OP( A ) = A'.
!
!  THE MATRIX X IS OVERWRITTEN ON B.
!
!  PARAMETERS
!  ==========
!
!  SIDE   - CHARACTER*1.
!           ON ENTRY, SIDE SPECIFIES WHETHER OP( A ) APPEARS ON THE LEFT
!           OR RIGHT OF X AS FOLLOWS:
!
!              SIDE = 'L' OR 'L'   OP( A )*X = ALPHA*B.
!
!              SIDE = 'R' OR 'R'   X*OP( A ) = ALPHA*B.
!
!           UNCHANGED ON EXIT.
!
!  UPLO   - CHARACTER*1.
!           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX A IS AN UPPER OR
!           LOWER TRIANGULAR MATRIX AS FOLLOWS:
!
!              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
!
!              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
!
!           UNCHANGED ON EXIT.
!
!  TRANSA - CHARACTER*1.
!           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
!           THE MATRIX MULTIPLICATION AS FOLLOWS:
!
!              TRANSA = 'N' OR 'N'   OP( A ) = A.
!
!              TRANSA = 'T' OR 'T'   OP( A ) = A'.
!
!              TRANSA = 'C' OR 'C'   OP( A ) = A'.
!
!           UNCHANGED ON EXIT.
!
!  DIAG   - CHARACTER*1.
!           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT TRIANGULAR
!           AS FOLLOWS:
!
!              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
!
!              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
!                                  TRIANGULAR.
!
!           UNCHANGED ON EXIT.
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF B. M MUST BE AT
!           LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF B.  N MUST BE
!           AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY,  ALPHA SPECIFIES THE SCALAR  ALPHA. WHEN  ALPHA IS
!           ZERO THEN  A IS NOT REFERENCED AND  B NEED NOT BE SET BEFORE
!           ENTRY.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, K ), WHERE K IS M
!           WHEN  SIDE = 'L' OR 'L'  AND IS  N  WHEN  SIDE = 'R' OR 'R'.
!           BEFORE ENTRY  WITH  UPLO = 'U' OR 'U',  THE  LEADING  K BY K
!           UPPER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE UPPER
!           TRIANGULAR MATRIX  AND THE STRICTLY LOWER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           BEFORE ENTRY  WITH  UPLO = 'L' OR 'L',  THE  LEADING  K BY K
!           LOWER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE LOWER
!           TRIANGULAR MATRIX  AND THE STRICTLY UPPER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           NOTE THAT WHEN  DIAG = 'U' OR 'U',  THE DIAGONAL ELEMENTS OF
!           A  ARE NOT REFERENCED EITHER,  BUT ARE ASSUMED TO BE  UNITY.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM.  WHEN  SIDE = 'L' OR 'L'  THEN
!           LDA  MUST BE AT LEAST  MAX( 1, M ),  WHEN  SIDE = 'R' OR 'R'
!           THEN LDA MUST BE AT LEAST MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, N ).
!           BEFORE ENTRY,  THE LEADING  M BY N PART OF THE ARRAY  B MUST
!           CONTAIN  THE  RIGHT-HAND  SIDE  MATRIX  B,  AND  ON EXIT  IS
!           OVERWRITTEN BY THE SOLUTION MATRIX  X.
!
!  LDB    - INTEGER.
!           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
!           IN  THE  CALLING  (SUB)  PROGRAM.   LDB  MUST  BE  AT  LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 3 BLAS ROUTINE.
!
!
!  -- WRITTEN ON 8-FEBRUARY-1989.
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!     IAIN DUFF, AERE HARWELL.
!     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
!     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
!
!
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
!     .. LOCAL SCALARS ..
!     .. PARAMETERS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      lside = lsame(side,'L') 
      if (lside) then 
        nrowa = m 
      else 
        nrowa = n 
      endif 
      nounit = lsame(diag,'N') 
      upper = lsame(uplo,'U') 
!
      info = 0 
      if (.not.lside .and. .not.lsame(side,'R')) then 
        info = 1 
      else if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = 2 
      else if (.not.lsame(transa,'N') .and. .not.lsame(transa,'T') .and. .not.&
          lsame(transa,'C')) then 
        info = 3 
      else if (.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then 
        info = 4 
      else if (m < 0) then 
        info = 5 
      else if (n < 0) then 
        info = 6 
      else if (lda < max(1,nrowa)) then 
        info = 9 
      else if (ldb < max(1,m)) then 
        info = 11 
      endif 
      if (info /= 0) then 
        call xerbla ('DTRSM ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (n == 0) return  
!
!     AND WHEN  ALPHA.EQ.ZERO.
!
      if (alpha == zero) then 
        b(:m,:n) = zero 
        return  
      endif 
!
!     START THE OPERATIONS.
!
      if (lside) then 
        if (lsame(transa,'N')) then 
!
!           FORM  B := ALPHA*INV( A )*B.
!
          if (upper) then 
            do j = 1, n 
              if (Abs(alpha-1.d0) > 1.d-20) then 
                b(:m,j) = alpha*b(:m,j) 
              endif 
              do k = m, 1, -1 
                if (b(k,j) == zero) cycle  
                if (nounit) b(k,j) = b(k,j)/a(k,k) 
                b(:k-1,j) = b(:k-1,j) - b(k,j)*a(:k-1,k) 
              end do 
            end do 
          else 
            do j = 1, n 
              if (Abs(alpha-1.d0) > 1.d-20) then 
                b(:m,j) = alpha*b(:m,j) 
              endif 
              do k = 1, m 
                if (b(k,j) == zero) cycle  
                if (nounit) b(k,j) = b(k,j)/a(k,k) 
                b(k+1:m,j) = b(k+1:m,j) - b(k,j)*a(k+1:m,k) 
              end do 
            end do 
          endif 
        else 
!
!           FORM  B := ALPHA*INV( A' )*B.
!
          if (upper) then 
            do j = 1, n 
              do i = 1, m 
                temp = alpha*b(i,j) 
                temp = temp - sum(a(:i-1,i)*b(:i-1,j)) 
                if (nounit) temp = temp/a(i,i) 
                b(i,j) = temp 
              end do 
            end do 
          else 
            do j = 1, n 
              do i = m, 1, -1 
                temp = alpha*b(i,j) 
                temp = temp - sum(a(i+1:m,i)*b(i+1:m,j)) 
                if (nounit) temp = temp/a(i,i) 
                b(i,j) = temp 
              end do 
            end do 
          endif 
        endif 
      else 
        if (lsame(transa,'N')) then 
!
!           FORM  B := ALPHA*B*INV( A ).
!
          if (upper) then 
            do j = 1, n 
              if (Abs(alpha-1.d0) > 1.d-20) then 
                b(:m,j) = alpha*b(:m,j) 
              endif 
              do k = 1, j - 1 
                if (a(k,j) == zero) cycle  
                b(:m,j) = b(:m,j) - a(k,j)*b(:m,k) 
              end do 
              if (.not.nounit) cycle  
              temp = one/a(j,j) 
              b(:m,j) = temp*b(:m,j) 
            end do 
          else 
            do j = n, 1, -1 
              if (Abs(alpha-1.d0) > 1.d-20) then 
                b(:m,j) = alpha*b(:m,j) 
              endif 
              do k = j + 1, n 
                if (a(k,j) == zero) cycle  
                b(:m,j) = b(:m,j) - a(k,j)*b(:m,k) 
              end do 
              if (.not.nounit) cycle  
              temp = one/a(j,j) 
              b(:m,j) = temp*b(:m,j) 
            end do 
          endif 
        else 
!
!           FORM  B := ALPHA*B*INV( A' ).
!
          if (upper) then 
            do k = n, 1, -1 
              if (nounit) then 
                temp = one/a(k,k) 
                b(:m,k) = temp*b(:m,k) 
              endif 
              do j = 1, k - 1 
                if (a(j,k) == zero) cycle  
                temp = a(j,k) 
                b(:m,j) = b(:m,j) - temp*b(:m,k) 
              end do 
              if (Abs(alpha-1.d0) < 1.d-20) cycle  
              b(:m,k) = alpha*b(:m,k) 
            end do 
          else 
            do k = 1, n 
              if (nounit) then 
                temp = one/a(k,k) 
                b(:m,k) = temp*b(:m,k) 
              endif 
              do j = k + 1, n 
                if (a(j,k) == zero) cycle  
                temp = a(j,k) 
                b(:m,j) = b(:m,j) - temp*b(:m,k) 
              end do 
              if (Abs(alpha-1.d0) < 1.d-20) cycle  
              b(:m,k) = alpha*b(:m,k) 
            end do 
          endif 
        endif 
      endif 
!
      return  
!
!     END OF DTRSM .
!
      end subroutine dtrsm 


      integer function idamax (n, dx, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      real(double) , intent(in) :: dx(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ix 
      real(double) :: dmax 
!-----------------------------------------------
!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!     MODIFIED TO CORRECT PROBLEM WITH NEGATIVE INCREMENT, 8/21/90.
!
!
      idamax = 0 
      if (n < 1) return  
      idamax = 1 
      if (n == 1) return  
      if (incx /= 1) then 
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
        ix = 1 
        if (incx < 0) ix = ((-n) + 1)*incx + 1 
        dmax = dabs(dx(ix)) 
        ix = ix + incx 
        do i = 2, n 
          if (dabs(dx(ix)) > dmax) then 
            idamax = i 
            dmax = dabs(dx(ix)) 
          endif 
          ix = ix + incx 
        end do 
        return  
      endif 
!
!        CODE FOR INCREMENT EQUAL TO 1
!
      dmax = dabs(dx(1)) 
      do i = 2, n 
        if (dabs(dx(i)) <= dmax) cycle  
        idamax = i 
        dmax = dabs(dx(i)) 
      end do 
      return  
      end function idamax 


      logical function lsame (ca, cb) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: ca 
      character , intent(in) :: cb 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: inta, intb, zcode 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC ichar 
!-----------------------------------------------
!
!  -- LAPACK AUXILIARY ROUTINE (VERSION 1.0) --
!     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
!     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
!     FEBRUARY 29, 1992
!
!     .. SCALAR ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  LSAME RETURNS .TRUE. IF CA IS THE SAME LETTER AS CB REGARDLESS OF
!  CASE.
!
!  ARGUMENTS
!  =========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CA AND CB SPECIFY THE SINGLE CHARACTERS TO BE COMPARED.
!
!     .. INTRINSIC FUNCTIONS ..
!     ..
!     .. LOCAL SCALARS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST IF THE CHARACTERS ARE EQUAL
!
      lsame = ca == cb 
      if (lsame) return  
!
!     NOW TEST FOR EQUIVALENCE IF BOTH CHARACTERS ARE ALPHABETIC.
!
      zcode = ichar('Z') 
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
      inta = ichar(ca) 
      intb = ichar(cb) 
!
      if (zcode==90 .or. zcode==122) then 
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
        if (inta>=97 .and. inta<=122) inta = inta - 32 
        if (intb>=97 .and. intb<=122) intb = intb - 32 
!
      else if (zcode==233 .or. zcode==169) then 
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
        if (inta>=129 .and. inta<=137 .or. inta>=145 .and. inta<=153 .or. inta&
          >=162 .and. inta<=169) inta = inta + 64 
        if (intb>=129 .and. intb<=137 .or. intb>=145 .and. intb<=153 .or. intb&
          >=162 .and. intb<=169) intb = intb + 64 
!
      else if (zcode==218 .or. zcode==250) then 
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
        if (inta>=225 .and. inta<=250) inta = inta - 32 
        if (intb>=225 .and. intb<=250) intb = intb - 32 
      endif 
      lsame = inta == intb 
      return  
!
!     RETURN
!
!     END OF LSAME
!
      end function lsame 




      subroutine dtrtri(uplo, diag, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use lsame_I 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: lda 
      integer  :: info 
      character  :: uplo 
      character  :: diag 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, jb, nb, nn 
      logical :: nounit, upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, min 
      integer, external :: ilaenv
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DTRTRI computes the inverse of a real upper or lower triangular
!  matrix A.
!
!  This is the Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are
!          assumed to be 1.
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!               matrix is singular and its inverse can not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      nounit = lsame(diag,'N') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (.not.nounit .and. .not.lsame(diag,'U')) then 
        info = -2 
      else if (n < 0) then 
        info = -3 
      else if (lda < max(1,n)) then 
        info = -5 
      endif 
      if (info /= 0) then 
        call xerbla ('DTRTRI', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
!     Check for singularity if non-unit.
!
      if (nounit) then 
        do info = 1, n 
          if (a(info,info) /= zero) cycle  
          return  
        end do 
        info = 0 
      endif 
!
!     Determine the block size for this environment.
!
      nb = ilaenv(1,'DTRTRI',uplo//diag,n,-1,-1,-1) 
      if (nb<=1 .or. nb>=n) then 
!
!        Use unblocked code
!
        call dtrti2 (uplo, diag, n, a, lda, info) 
      else 
!
!        Use blocked code
!
        if (upper) then 
!
!           Compute inverse of upper triangular matrix
!
          do j = 1, n, nb 
            jb = min(nb,n - j + 1) 
!
!              Compute rows 1:j-1 of current block column
!
            call dtrmm ('Left', 'Upper', 'No transpose', diag, j - 1, jb, one, &
              a, lda, a(1,j), lda) 
            call dtrsm ('Right', 'Upper', 'No transpose', diag, j - 1, jb, (-&
              one), a(j,j), lda, a(1,j), lda) 
!
!              Compute inverse of current diagonal block
!
            call dtrti2 ('Upper', diag, jb, a(j,j), lda, info) 
          end do 
        else 
!
!           Compute inverse of lower triangular matrix
!
          nn = ((n - 1)/nb)*nb + 1 
          do j = nn, 1, -nb 
            jb = min(nb,n - j + 1) 
            if (j + jb <= n) then 
!
!                 Compute rows j+jb:n of current block column
!
              call dtrmm ('Left', 'Lower', 'No transpose', diag, n - j - jb + 1&
                , jb, one, a(j+jb,j+jb), lda, a(j+jb,j), lda) 
              call dtrsm ('Right', 'Lower', 'No transpose', diag, n - j - jb + &
                1, jb, (-one), a(j,j), lda, a(j+jb,j), lda) 
            endif 
!
!              Compute inverse of current diagonal block
!
            call dtrti2 ('Lower', diag, jb, a(j,j), lda, info) 
          end do 
        endif 
      endif 
!
      return  
!
!     End of DTRTRI
!
      end subroutine dtrtri 


 
      subroutine dlauum(uplo, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use lsame_I 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: n 
      integer  :: lda 
      integer  :: info 
      character  :: uplo 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ib, nb 
      logical :: upper 
      integer, external :: ilaenv
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, min 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAUUM computes the product U * U' or L' * L, where the triangular
!  factor U or L is stored in the upper or lower triangular part of
!  the array A.
!
!  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
!  overwriting the factor U in A.
!  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
!  overwriting the factor L in A.
!
!  This is the blocked form of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the triangular factor stored in the array A
!          is upper or lower triangular:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the triangular factor U or L.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular factor U or L.
!          On exit, if UPLO = 'U', the upper triangle of A is
!          overwritten with the upper triangle of the product U * U';
!          if UPLO = 'L', the lower triangle of A is overwritten with
!          the lower triangle of the product L' * L.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DLAUUM', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
!     Determine the block size for this environment.
!
      nb = ilaenv(1,'DLAUUM',uplo,n,-1,-1,-1) 
!
      if (nb<=1 .or. nb>=n) then 
!
!        Use unblocked code
!
        call dlauu2 (uplo, n, a, lda, info) 
      else 
!
!        Use blocked code
!
        if (upper) then 
!
!           Compute the product U * U'.
!
          do i = 1, n, nb 
            ib = min(nb,n - i + 1) 
            call dtrmm ('Right', 'Upper', 'Transpose', 'Non-unit', i - 1, ib, &
              one, a(i,i), lda, a(1,i), lda) 
            call dlauu2 ('Upper', ib, a(i,i), lda, info) 
            if (i + ib > n) cycle  
            call dgemm ('No transpose', 'Transpose', i - 1, ib, n - i - ib + 1&
              , one, a(1,i+ib), lda, a(i,i+ib), lda, one, a(1,i), lda) 
            call dsyrk ('Upper', 'No transpose', ib, n - i - ib + 1, one, a(i,i&
              +ib), lda, one, a(i,i), lda) 
          end do 
        else 
!
!           Compute the product L' * L.
!
          do i = 1, n, nb 
            ib = min(nb,n - i + 1) 
            call dtrmm ('Left', 'Lower', 'Transpose', 'Non-unit', ib, i - 1, &
              one, a(i,i), lda, a(i,1), lda) 
            call dlauu2 ('Lower', ib, a(i,i), lda, info) 
            if (i + ib > n) cycle  
            call dgemm ('Transpose', 'No transpose', ib, i - 1, n - i - ib + 1&
              , one, a(i+ib,i), lda, a(i+ib,1), lda, one, a(i,1), lda) 
            call dsyrk ('Lower', 'Transpose', ib, n - i - ib + 1, one, a(i+ib,i&
              ), lda, one, a(i,i), lda) 
          end do 
        endif 
      endif 
!
      return  
!
!     End of DLAUUM
!
      end subroutine dlauum 


 
 
      subroutine dlauu2(uplo, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use ddot_I 
   !   use dgemv_I 
      use dscal_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer  :: lda 
      integer , intent(out) :: info 
      character  :: uplo 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(double) :: aii 
      logical :: upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DLAUU2 computes the product U * U' or L' * L, where the triangular
!  factor U or L is stored in the upper or lower triangular part of
!  the array A.
!
!  If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
!  overwriting the factor U in A.
!  If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
!  overwriting the factor L in A.
!
!  This is the unblocked form of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the triangular factor stored in the array A
!          is upper or lower triangular:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the triangular factor U or L.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular factor U or L.
!          On exit, if UPLO = 'U', the upper triangle of A is
!          overwritten with the upper triangle of the product U * U';
!          if UPLO = 'L', the lower triangle of A is overwritten with
!          the lower triangle of the product L' * L.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DLAUU2', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
      if (upper) then 
!
!        Compute the product U * U'.
!
        do i = 1, n 
          aii = a(i,i) 
          if (i < n) then 
            a(i,i) = ddot(n - i + 1,a(i,i),lda,a(i,i),lda) 
            call dgemv ('No transpose', i - 1, n - i, one, a(1,i+1), lda, a(i,i&
              +1), lda, aii, a(1,i), 1) 
          else 
            call dscal (i, aii, a(1,i), 1) 
          endif 
        end do 
!
      else 
!
!        Compute the product L' * L.
!
        do i = 1, n 
          aii = a(i,i) 
          if (i < n) then 
            a(i,i) = ddot(n - i + 1,a(i,i),1,a(i,i),1) 
            call dgemv ('Transpose', n - i, i - 1, one, a(i+1,1), lda, a(i+1,i)&
              , 1, aii, a(i,1), lda) 
          else 
            call dscal (i, aii, a(i,1), lda) 
          endif 
        end do 
      endif 
!
      return  
!
!     End of DLAUU2
!
      end subroutine dlauu2 


 
      subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: k 
      integer , intent(in) :: lda 
      integer , intent(in) :: ldc 
      real(double) , intent(in) :: alpha 
      real(double) , intent(in) :: beta 
      character  :: uplo 
      character  :: trans 
      real(double) , intent(in) :: a(lda,*) 
      real(double) , intent(out)  :: c(ldc,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, info, j, l, nrowa 
      real(double) :: temp 
      logical :: upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSYRK  performs one of the symmetric rank k operations
!
!     C := alpha*A*A' + beta*C,
!
!  or
!
!     C := alpha*A'*A + beta*C,
!
!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
!  in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
!
!              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrix  A.  K must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
!     .. External Subroutines ..
!     .. Intrinsic Functions ..
!     .. Local Scalars ..
!     .. Parameters ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      if (lsame(trans,'N')) then 
        nrowa = n 
      else 
        nrowa = k 
      endif 
      upper = lsame(uplo,'U') 
!
      info = 0 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = 1 
      else if (.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. .not.&
          lsame(trans,'C')) then 
        info = 2 
      else if (n < 0) then 
        info = 3 
      else if (k < 0) then 
        info = 4 
      else if (lda < max(1,nrowa)) then 
        info = 7 
      else if (ldc < max(1,n)) then 
        info = 10 
      endif 
      if (info /= 0) then 
        call xerbla ('DSYRK ', info) 
        return  
      endif 
!
!     Quick return if possible.
!
      if (n==0 .or. (Abs(alpha) < 1.d-20 .or. k==0) .and. Abs(beta - one) < 1.d-20) return  
!
!     And when  alpha.eq.zero.
!
      if (alpha == zero) then 
        if (upper) then 
          if (beta == zero) then 
            do j = 1, n 
              c(:j,j) = zero 
            end do 
          else 
            do j = 1, n 
              c(:j,j) = beta*c(:j,j) 
            end do 
          endif 
        else 
          if (beta == zero) then 
            do j = 1, n 
              c(j:n,j) = zero 
            end do 
          else 
            do j = 1, n 
              c(j:n,j) = beta*c(j:n,j) 
            end do 
          endif 
        endif 
        return  
      endif 
!
!     Start the operations.
!
      if (lsame(trans,'N')) then 
!
!        Form  C := alpha*A*A' + beta*C.
!
        if (upper) then 
          do j = 1, n 
            if (beta == zero) then 
              c(:j,j) = zero 
            else if (Abs(beta - one) > 1.d-20) then 
              c(:j,j) = beta*c(:j,j) 
            endif 
            do l = 1, k 
              if (a(j,l) == zero) cycle  
              temp = alpha*a(j,l) 
              c(:j,j) = c(:j,j) + temp*a(:j,l) 
            end do 
          end do 
        else 
          do j = 1, n 
            if (beta == zero) then 
              c(j:n,j) = zero 
            else if (Abs(beta - one) > 1.d-20) then 
              c(j:n,j) = beta*c(j:n,j) 
            endif 
            do l = 1, k 
              if (a(j,l) == zero) cycle  
              temp = alpha*a(j,l) 
              c(j:n,j) = c(j:n,j) + temp*a(j:n,l) 
            end do 
          end do 
        endif 
      else 
!
!        Form  C := alpha*A'*A + beta*C.
!
        if (upper) then 
          do j = 1, n 
            do i = 1, j 
              temp = zero 
              temp = temp + sum(a(:k,i)*a(:k,j)) 
              if (beta == zero) then 
                c(i,j) = alpha*temp 
              else 
                c(i,j) = alpha*temp + beta*c(i,j) 
              endif 
            end do 
          end do 
        else 
          do j = 1, n 
            do i = j, n 
              temp = zero 
              temp = temp + sum(a(:k,i)*a(:k,j)) 
              if (beta == zero) then 
                c(i,j) = alpha*temp 
              else 
                c(i,j) = alpha*temp + beta*c(i,j) 
              endif 
            end do 
          end do 
        endif 
      endif 
!
      return  
!
!     End of DSYRK .
!
      end subroutine dsyrk 

      
 
 
      subroutine dpotf2(uplo, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use ddot_I 
   !   use dgemv_I 
      use dscal_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer  :: lda 
      integer , intent(out) :: info 
      character  :: uplo 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j 
      real(double) :: ajj 
      logical :: upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max, sqrt 
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (n < 0) then 
        info = -2 
      else if (lda < max(1,n)) then 
        info = -4 
      endif 
      if (info /= 0) then 
        call xerbla ('DPOTF2', (-info)) 
        return  
      endif 
!
!     Quick return if possible
!
      if (n == 0) return  
!
      if (upper) then 
!
!        Compute the Cholesky factorization A = U'*U.
!
        do j = 1, n 
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
          ajj = a(j,j) - ddot(j - 1,a(1,j),1,a(1,j),1) 
          if (ajj <= zero) then 
            a(j,j) = ajj 
            go to 30 
          endif 
          ajj = sqrt(ajj) 
          a(j,j) = ajj 
!
!           Compute elements J+1:N of row J.
!
          if (j >= n) cycle  
          call dgemv ('Transpose', j - 1, n - j, (-one), a(1,j+1), lda, a(1,j)&
            , 1, one, a(j,j+1), lda) 
          call dscal (n - j, one/ajj, a(j,j+1), lda) 
        end do 
      else 
!
!        Compute the Cholesky factorization A = L*L'.
!
        do j = 1, n 
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
          ajj = a(j,j) - ddot(j - 1,a(j,1),lda,a(j,1),lda) 
          if (ajj <= zero) then 
            a(j,j) = ajj 
            go to 30 
          endif 
          ajj = sqrt(ajj) 
          a(j,j) = ajj 
!
!           Compute elements J+1:N of column J.
!
          if (j >= n) cycle  
          call dgemv ('No transpose', n - j, j - 1, (-one), a(j+1,1), lda, a(j,&
            1), lda, one, a(j+1,j), 1) 
          call dscal (n - j, one/ajj, a(j+1,j), 1) 
        end do 
      endif 
      go to 40 
!
   30 continue 
      info = j 
!
   40 continue 
      return  
!
!     End of DPOTF2
!
      end subroutine dpotf2 
 
 


      subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: lda 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      real(double) , intent(in) :: alpha 
      real(double) , intent(in) :: beta 
      character  :: trans*(*) 
      real(double) , intent(in) :: a(lda,*) 
      real(double) , intent(in) :: x(*) 
      real(double) , intent(inout) :: y(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: info, ix, iy, j, jx, jy, kx, ky, lenx, leny 
      real(double) :: temp 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DGEMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
!
!     Y := ALPHA*A*X + BETA*Y,   OR   Y := ALPHA*A'*X + BETA*Y,
!
!  WHERE ALPHA AND BETA ARE SCALARS, X AND Y ARE VECTORS AND A IS AN
!  M BY N MATRIX.
!
!  PARAMETERS
!  ==========
!
!  TRANS  - CHARACTER*1.
!           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
!           FOLLOWS:
!
!              TRANS = 'N' OR 'N'   Y := ALPHA*A*X + BETA*Y.
!
!              TRANS = 'T' OR 'T'   Y := ALPHA*A'*X + BETA*Y.
!
!              TRANS = 'C' OR 'C'   Y := ALPHA*A'*X + BETA*Y.
!
!           UNCHANGED ON EXIT.
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
!           M MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
!           N MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
!           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
!           CONTAIN THE MATRIX OF COEFFICIENTS.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCX ) ) WHEN TRANS = 'N' OR 'N'
!           AND AT LEAST
!           ( 1 + ( M - 1 )*ABS( INCX ) ) OTHERWISE.
!           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE
!           VECTOR X.
!           UNCHANGED ON EXIT.
!
!  INCX   - INTEGER.
!           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           X. INCX MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!  BETA   - DOUBLE PRECISION.
!           ON ENTRY, BETA SPECIFIES THE SCALAR BETA. WHEN BETA IS
!           SUPPLIED AS ZERO THEN Y NEED NOT BE SET ON INPUT.
!           UNCHANGED ON EXIT.
!
!  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( M - 1 )*ABS( INCY ) ) WHEN TRANS = 'N' OR 'N'
!           AND AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCY ) ) OTHERWISE.
!           BEFORE ENTRY WITH BETA NON-ZERO, THE INCREMENTED ARRAY Y
!           MUST CONTAIN THE VECTOR Y. ON EXIT, Y IS OVERWRITTEN BY THE
!           UPDATED VECTOR Y.
!
!  INCY   - INTEGER.
!           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           Y. INCY MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 2 BLAS ROUTINE.
!
!  -- WRITTEN ON 22-OCTOBER-1986.
!     JACK DONGARRA, ARGONNE NATIONAL LAB.
!     JEREMY DU CROZ, NAG CENTRAL OFFICE.
!     SVEN HAMMARLING, NAG CENTRAL OFFICE.
!     RICHARD HANSON, SANDIA NATIONAL LABS.
!
!
!     .. PARAMETERS ..
!     .. LOCAL SCALARS ..
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      info = 0 
      if (.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. .not.lsame(&
        trans,'C')) then 
        info = 1 
      else if (m < 0) then 
        info = 2 
      else if (n < 0) then 
        info = 3 
      else if (lda < max(1,m)) then 
        info = 6 
      else if (incx == 0) then 
        info = 8 
      else if (incy == 0) then 
        info = 11 
      endif 
      if (info /= 0) then 
        call xerbla ('DGEMV ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (m==0 .or. n==0 .or. Abs(alpha) < 1.d-20 .and. Abs(beta - 1.d0) < 1.d-20) return  
!
!     SET  LENX  AND  LENY, THE LENGTHS OF THE VECTORS X AND Y, AND SET
!     UP THE START POINTS IN  X  AND  Y.
!
      if (lsame(trans,'N')) then 
        lenx = n 
        leny = m 
      else 
        lenx = m 
        leny = n 
      endif 
      if (incx > 0) then 
        kx = 1 
      else 
        kx = 1 - (lenx - 1)*incx 
      endif 
      if (incy > 0) then 
        ky = 1 
      else 
        ky = 1 - (leny - 1)*incy 
      endif 
!
!     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
!     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
!
!     FIRST FORM  Y := BETA*Y.
!
      if (Abs(beta - 1.d0) > 1.d-20) then 
        if (incy == 1) then 
          if (beta == zero) then 
            y(:leny) = zero 
          else 
            y(:leny) = beta*y(:leny) 
          endif 
        else 
          iy = ky 
          if (Abs(beta) < 1.d-20) then 
            y(iy:(leny-1)*incy+iy:incy) = zero 
          else 
            y(iy:(leny-1)*incy+iy:incy) = beta*y(iy:(leny-1)*incy+iy:incy) 
          endif 
        endif 
      endif 
      if (Abs(alpha) < 1.d-20) return  
      if (lsame(trans,'N')) then 
!
!        FORM  Y := ALPHA*A*X + Y.
!
        jx = kx 
        if (incy == 1) then 
          do j = 1, n 
            if (x(jx) /= zero) then 
              temp = alpha*x(jx) 
              y(:m) = y(:m) + temp*a(:m,j) 
            endif 
            jx = jx + incx 
          end do 
        else 
          do j = 1, n 
            if (x(jx) /= zero) then 
              temp = alpha*x(jx) 
              iy = ky 
              y(iy:(m-1)*incy+iy:incy) = y(iy:(m-1)*incy+iy:incy) + temp*a(:m,j&
                ) 
            endif 
            jx = jx + incx 
          end do 
        endif 
      else 
!
!        FORM  Y := ALPHA*A'*X + Y.
!
        jy = ky 
        if (incx == 1) then 
          do j = 1, n 
            temp = zero 
            temp = temp + sum(a(:m,j)*x(:m)) 
            y(jy) = y(jy) + alpha*temp 
            jy = jy + incy 
          end do 
        else 
          do j = 1, n 
            temp = zero 
            ix = kx 
            temp = temp + sum(a(:m,j)*x(ix:(m-1)*incx+ix:incx)) 
            y(jy) = y(jy) + alpha*temp 
            jy = jy + incy 
          end do 
        endif 
      endif 
!
      return  
!
!     END OF DGEMV .
!
      end subroutine dgemv 


      subroutine dger(m, n, alpha, x, incx, y, incy, a, lda) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: incx 
      integer , intent(in) :: incy 
      integer , intent(in) :: lda 
      real(double) , intent(in) :: alpha 
      real(double) , intent(in) :: x(*) 
      real(double) , intent(in) :: y(*) 
      real(double) , intent(inout) :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: info, ix, j, jy, kx 
      real(double) :: temp 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DGER   PERFORMS THE RANK 1 OPERATION
!
!     A := ALPHA*X*Y' + A,
!
!  WHERE ALPHA IS A SCALAR, X IS AN M ELEMENT VECTOR, Y IS AN N ELEMENT
!  VECTOR AND A IS AN M BY N MATRIX.
!
!  PARAMETERS
!  ==========
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF THE MATRIX A.
!           M MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF THE MATRIX A.
!           N MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY, ALPHA SPECIFIES THE SCALAR ALPHA.
!           UNCHANGED ON EXIT.
!
!  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( M - 1 )*ABS( INCX ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE M
!           ELEMENT VECTOR X.
!           UNCHANGED ON EXIT.
!
!  INCX   - INTEGER.
!           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           X. INCX MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!  Y      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCY ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY Y MUST CONTAIN THE N
!           ELEMENT VECTOR Y.
!           UNCHANGED ON EXIT.
!
!  INCY   - INTEGER.
!           ON ENTRY, INCY SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           Y. INCY MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
!           BEFORE ENTRY, THE LEADING M BY N PART OF THE ARRAY A MUST
!           CONTAIN THE MATRIX OF COEFFICIENTS. ON EXIT, A IS
!           OVERWRITTEN BY THE UPDATED MATRIX.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 2 BLAS ROUTINE.
!
!  -- WRITTEN ON 22-OCTOBER-1986.
!     JACK DONGARRA, ARGONNE NATIONAL LAB.
!     JEREMY DU CROZ, NAG CENTRAL OFFICE.
!     SVEN HAMMARLING, NAG CENTRAL OFFICE.
!     RICHARD HANSON, SANDIA NATIONAL LABS.
!
!
!     .. PARAMETERS ..
!     .. LOCAL SCALARS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      info = 0 
      if (m < 0) then 
        info = 1 
      else if (n < 0) then 
        info = 2 
      else if (incx == 0) then 
        info = 5 
      else if (incy == 0) then 
        info = 7 
      else if (lda < max(1,m)) then 
        info = 9 
      endif 
      if (info /= 0) then 
        call xerbla ('DGER  ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (m==0 .or. n==0 .or. alpha==zero) return  
!
!     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
!     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
!
      if (incy > 0) then 
        jy = 1 
      else 
        jy = 1 - (n - 1)*incy 
      endif 
      if (incx == 1) then 
        do j = 1, n 
          if (y(jy) /= zero) then 
            temp = alpha*y(jy) 
            a(:m,j) = a(:m,j) + x(:m)*temp 
          endif 
          jy = jy + incy 
        end do 
      else 
        if (incx > 0) then 
          kx = 1 
        else 
          kx = 1 - (m - 1)*incx 
        endif 
        do j = 1, n 
          if (y(jy) /= zero) then 
            temp = alpha*y(jy) 
            ix = kx 
            a(:m,j) = a(:m,j) + x(ix:(m-1)*incx+ix:incx)*temp 
          endif 
          jy = jy + incy 
        end do 
      endif 
!
      return  
!
!     END OF DGER  .
!
      end subroutine dger 



      subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: m 
      integer , intent(in) :: n 
      integer , intent(in) :: lda 
      integer , intent(in) :: ldb 
      real(double) , intent(in) :: alpha 
      character  :: side*(*) 
      character  :: uplo*(*) 
      character  :: transa*(*) 
      character  :: diag*(*) 
      real(double) , intent(in) :: a(lda,*) 
      real(double) , intent(inout) :: b(ldb,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, info, j, k, nrowa 
      real(double) :: temp 
      logical :: lside, nounit, upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DTRMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
!
!     B := ALPHA*OP( A )*B,   OR   B := ALPHA*B*OP( A ),
!
!  WHERE  ALPHA  IS A SCALAR,  B  IS AN M BY N MATRIX,  A  IS A UNIT, OR
!  NON-UNIT,  UPPER OR LOWER TRIANGULAR MATRIX  AND  OP( A )  IS ONE  OF
!
!     OP( A ) = A   OR   OP( A ) = A'.
!
!  PARAMETERS
!  ==========
!
!  SIDE   - CHARACTER*1.
!           ON ENTRY,  SIDE SPECIFIES WHETHER  OP( A ) MULTIPLIES B FROM
!           THE LEFT OR RIGHT AS FOLLOWS:
!
!              SIDE = 'L' OR 'L'   B := ALPHA*OP( A )*B.
!
!              SIDE = 'R' OR 'R'   B := ALPHA*B*OP( A ).
!
!           UNCHANGED ON EXIT.
!
!  UPLO   - CHARACTER*1.
!           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX A IS AN UPPER OR
!           LOWER TRIANGULAR MATRIX AS FOLLOWS:
!
!              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
!
!              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
!
!           UNCHANGED ON EXIT.
!
!  TRANSA - CHARACTER*1.
!           ON ENTRY, TRANSA SPECIFIES THE FORM OF OP( A ) TO BE USED IN
!           THE MATRIX MULTIPLICATION AS FOLLOWS:
!
!              TRANSA = 'N' OR 'N'   OP( A ) = A.
!
!              TRANSA = 'T' OR 'T'   OP( A ) = A'.
!
!              TRANSA = 'C' OR 'C'   OP( A ) = A'.
!
!           UNCHANGED ON EXIT.
!
!  DIAG   - CHARACTER*1.
!           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT TRIANGULAR
!           AS FOLLOWS:
!
!              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
!
!              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
!                                  TRIANGULAR.
!
!           UNCHANGED ON EXIT.
!
!  M      - INTEGER.
!           ON ENTRY, M SPECIFIES THE NUMBER OF ROWS OF B. M MUST BE AT
!           LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE NUMBER OF COLUMNS OF B.  N MUST BE
!           AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  ALPHA  - DOUBLE PRECISION.
!           ON ENTRY,  ALPHA SPECIFIES THE SCALAR  ALPHA. WHEN  ALPHA IS
!           ZERO THEN  A IS NOT REFERENCED AND  B NEED NOT BE SET BEFORE
!           ENTRY.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, K ), WHERE K IS M
!           WHEN  SIDE = 'L' OR 'L'  AND IS  N  WHEN  SIDE = 'R' OR 'R'.
!           BEFORE ENTRY  WITH  UPLO = 'U' OR 'U',  THE  LEADING  K BY K
!           UPPER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE UPPER
!           TRIANGULAR MATRIX  AND THE STRICTLY LOWER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           BEFORE ENTRY  WITH  UPLO = 'L' OR 'L',  THE  LEADING  K BY K
!           LOWER TRIANGULAR PART OF THE ARRAY  A MUST CONTAIN THE LOWER
!           TRIANGULAR MATRIX  AND THE STRICTLY UPPER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           NOTE THAT WHEN  DIAG = 'U' OR 'U',  THE DIAGONAL ELEMENTS OF
!           A  ARE NOT REFERENCED EITHER,  BUT ARE ASSUMED TO BE  UNITY.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM.  WHEN  SIDE = 'L' OR 'L'  THEN
!           LDA  MUST BE AT LEAST  MAX( 1, M ),  WHEN  SIDE = 'R' OR 'R'
!           THEN LDA MUST BE AT LEAST MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  B      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDB, N ).
!           BEFORE ENTRY,  THE LEADING  M BY N PART OF THE ARRAY  B MUST
!           CONTAIN THE MATRIX  B,  AND  ON EXIT  IS OVERWRITTEN  BY THE
!           TRANSFORMED MATRIX.
!
!  LDB    - INTEGER.
!           ON ENTRY, LDB SPECIFIES THE FIRST DIMENSION OF B AS DECLARED
!           IN  THE  CALLING  (SUB)  PROGRAM.   LDB  MUST  BE  AT  LEAST
!           MAX( 1, M ).
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 3 BLAS ROUTINE.
!
!  -- WRITTEN ON 8-FEBRUARY-1989.
!     JACK DONGARRA, ARGONNE NATIONAL LABORATORY.
!     IAIN DUFF, AERE HARWELL.
!     JEREMY DU CROZ, NUMERICAL ALGORITHMS GROUP LTD.
!     SVEN HAMMARLING, NUMERICAL ALGORITHMS GROUP LTD.
!
!
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!        EXTERNAL           XRBLA
!     .. INTRINSIC FUNCTIONS ..
!     .. LOCAL SCALARS ..
!     .. PARAMETERS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      lside = lsame(side,'L') 
      if (lside) then 
        nrowa = m 
      else 
        nrowa = n 
      endif 
      nounit = lsame(diag,'N') 
      upper = lsame(uplo,'U') 
!
      info = 0 
      if (.not.lside .and. .not.lsame(side,'R')) then 
        info = 1 
      else if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = 2 
      else if (.not.lsame(transa,'N') .and. .not.lsame(transa,'T') .and. .not.&
          lsame(transa,'C')) then 
        info = 3 
      else if (.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then 
        info = 4 
      else if (m < 0) then 
        info = 5 
      else if (n < 0) then 
        info = 6 
      else if (lda < max(1,nrowa)) then 
        info = 9 
      else if (ldb < max(1,m)) then 
        info = 11 
      endif 
      if (info /= 0) then 
        call xerbla ('DTRMM ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (n == 0) return  
!
!     AND WHEN  ALPHA.EQ.ZERO.
!
      if (alpha == zero) then 
        b(:m,:n) = zero 
        return  
      endif 
!
!     START THE OPERATIONS.
!
      if (lside) then 
        if (lsame(transa,'N')) then 
!
!           FORM  B := ALPHA*A*B.
!
          if (upper) then 
            do j = 1, n 
              do k = 1, m 
                if (b(k,j) == zero) cycle  
                temp = alpha*b(k,j) 
                b(:k-1,j) = b(:k-1,j) + temp*a(:k-1,k) 
                if (nounit) temp = temp*a(k,k) 
                b(k,j) = temp 
              end do 
            end do 
          else 
            do j = 1, n 
              do k = m, 1, -1 
                if (b(k,j) == zero) cycle  
                temp = alpha*b(k,j) 
                b(k,j) = temp 
                if (nounit) b(k,j) = b(k,j)*a(k,k) 
                b(k+1:m,j) = b(k+1:m,j) + temp*a(k+1:m,k) 
              end do 
            end do 
          endif 
        else 
!
!           FORM  B := ALPHA*B*A'.
!
          if (upper) then 
            do j = 1, n 
              do i = m, 1, -1 
                temp = b(i,j) 
                if (nounit) temp = temp*a(i,i) 
                temp = temp + sum(a(:i-1,i)*b(:i-1,j)) 
                b(i,j) = alpha*temp 
              end do 
            end do 
          else 
            do j = 1, n 
              do i = 1, m 
                temp = b(i,j) 
                if (nounit) temp = temp*a(i,i) 
                temp = temp + sum(a(i+1:m,i)*b(i+1:m,j)) 
                b(i,j) = alpha*temp 
              end do 
            end do 
          endif 
        endif 
      else 
        if (lsame(transa,'N')) then 
!
!           FORM  B := ALPHA*B*A.
!
          if (upper) then 
            do j = n, 1, -1 
              temp = alpha 
              if (nounit) temp = temp*a(j,j) 
              b(:m,j) = temp*b(:m,j) 
              do k = 1, j - 1 
                if (a(k,j) == zero) cycle  
                temp = alpha*a(k,j) 
                b(:m,j) = b(:m,j) + temp*b(:m,k) 
              end do 
            end do 
          else 
            do j = 1, n 
              temp = alpha 
              if (nounit) temp = temp*a(j,j) 
              b(:m,j) = temp*b(:m,j) 
              do k = j + 1, n 
                if (a(k,j) == zero) cycle  
                temp = alpha*a(k,j) 
                b(:m,j) = b(:m,j) + temp*b(:m,k) 
              end do 
            end do 
          endif 
        else 
!
!           FORM  B := ALPHA*B*A'.
!
          if (upper) then 
            do k = 1, n 
              do j = 1, k - 1 
                if (a(j,k) == zero) cycle  
                temp = alpha*a(j,k) 
                b(:m,j) = b(:m,j) + temp*b(:m,k) 
              end do 
              temp = alpha 
              if (nounit) temp = temp*a(k,k) 
              if (Abs(temp - 1.d0) < 1.d-20) cycle  
              b(:m,k) = temp*b(:m,k) 
            end do 
          else 
            do k = n, 1, -1 
              do j = k + 1, n 
                if (a(j,k) == zero) cycle  
                temp = alpha*a(j,k) 
                b(:m,j) = b(:m,j) + temp*b(:m,k) 
              end do 
              temp = alpha 
              if (nounit) temp = temp*a(k,k) 
              if (Abs(temp - 1.d0) < 1.d-20) cycle  
              b(:m,k) = temp*b(:m,k) 
            end do 
          endif 
        endif 
      endif 
!
      return  
!
!     END OF DTRMM .
!
      end subroutine dtrmm 

      

      subroutine dtrti2(uplo, diag, n, a, lda, info) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  20:49:17  01/01/08  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use dscal_I 
      use dtrmv_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer  :: lda 
      integer , intent(out) :: info 
      character  :: uplo 
      character  :: diag 
      real(double)  :: a(lda,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j 
      real(double) :: ajj 
      logical :: nounit, upper 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
!     ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DTRTI2 computes the inverse of a real upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading n by n upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are
!          assumed to be 1.
!
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0 
      upper = lsame(uplo,'U') 
      nounit = lsame(diag,'N') 
      if (.not.upper .and. .not.lsame(uplo,'L')) then 
        info = -1 
      else if (.not.nounit .and. .not.lsame(diag,'U')) then 
        info = -2 
      else if (n < 0) then 
        info = -3 
      else if (lda < max(1,n)) then 
        info = -5 
      endif 
      if (info /= 0) then 
        call xerbla ('DTRTI2', (-info)) 
        return  
      endif 
!
      if (upper) then 
!
!        Compute inverse of upper triangular matrix.
!
        do j = 1, n 
          if (nounit) then 
            a(j,j) = one/a(j,j) 
            ajj = -a(j,j) 
          else 
            ajj = -one 
          endif 
!
!           Compute elements 1:j-1 of j-th column.
!
          call dtrmv ('Upper', 'No transpose', diag, j - 1, a, lda, a(1,j), 1) 
          call dscal (j - 1, ajj, a(1,j), 1) 
        end do 
      else 
!
!        Compute inverse of lower triangular matrix.
!
        do j = n, 1, -1 
          if (nounit) then 
            a(j,j) = one/a(j,j) 
            ajj = -a(j,j) 
          else 
            ajj = -one 
          endif 
          if (j >= n) cycle  
!
!              Compute elements j+1:n of j-th column.
!
          call dtrmv ('Lower', 'No transpose', diag, n - j, a(j+1,j+1), lda, a(&
            j+1,j), 1) 
          call dscal (n - j, ajj, a(j+1,j), 1) 
        end do 
      endif 
!
      return  
!
!     End of DTRTI2
!
      end subroutine dtrti2 



      subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  22:48:56  03/08/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lsame_I 
      use xerbla_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: lda 
      integer , intent(in) :: incx 
      character  :: uplo*(*) 
      character  :: trans*(*) 
      character  :: diag*(*) 
      real(double) , intent(in) :: a(lda,*) 
      real(double) , intent(inout) :: x(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: info, ix, j, jx, kx 
      real(double) :: temp 
      logical :: nounit 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC max 
!-----------------------------------------------
!     .. SCALAR ARGUMENTS ..
!     .. ARRAY ARGUMENTS ..
!     ..
!
!  PURPOSE
!  =======
!
!  DTRMV  PERFORMS ONE OF THE MATRIX-VECTOR OPERATIONS
!
!     X := A*X,   OR   X := A'*X,
!
!  WHERE X IS AN N ELEMENT VECTOR AND  A IS AN N BY N UNIT, OR NON-UNIT,
!  UPPER OR LOWER TRIANGULAR MATRIX.
!
!  PARAMETERS
!  ==========
!
!  UPLO   - CHARACTER*1.
!           ON ENTRY, UPLO SPECIFIES WHETHER THE MATRIX IS AN UPPER OR
!           LOWER TRIANGULAR MATRIX AS FOLLOWS:
!
!              UPLO = 'U' OR 'U'   A IS AN UPPER TRIANGULAR MATRIX.
!
!              UPLO = 'L' OR 'L'   A IS A LOWER TRIANGULAR MATRIX.
!
!           UNCHANGED ON EXIT.
!
!  TRANS  - CHARACTER*1.
!           ON ENTRY, TRANS SPECIFIES THE OPERATION TO BE PERFORMED AS
!           FOLLOWS:
!
!              TRANS = 'N' OR 'N'   X := A*X.
!
!              TRANS = 'T' OR 'T'   X := A'*X.
!
!              TRANS = 'C' OR 'C'   X := A'*X.
!
!           UNCHANGED ON EXIT.
!
!  DIAG   - CHARACTER*1.
!           ON ENTRY, DIAG SPECIFIES WHETHER OR NOT A IS UNIT
!           TRIANGULAR AS FOLLOWS:
!
!              DIAG = 'U' OR 'U'   A IS ASSUMED TO BE UNIT TRIANGULAR.
!
!              DIAG = 'N' OR 'N'   A IS NOT ASSUMED TO BE UNIT
!                                  TRIANGULAR.
!
!           UNCHANGED ON EXIT.
!
!  N      - INTEGER.
!           ON ENTRY, N SPECIFIES THE ORDER OF THE MATRIX A.
!           N MUST BE AT LEAST ZERO.
!           UNCHANGED ON EXIT.
!
!  A      - DOUBLE PRECISION ARRAY OF DIMENSION ( LDA, N ).
!           BEFORE ENTRY WITH  UPLO = 'U' OR 'U', THE LEADING N BY N
!           UPPER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE UPPER
!           TRIANGULAR MATRIX AND THE STRICTLY LOWER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           BEFORE ENTRY WITH UPLO = 'L' OR 'L', THE LEADING N BY N
!           LOWER TRIANGULAR PART OF THE ARRAY A MUST CONTAIN THE LOWER
!           TRIANGULAR MATRIX AND THE STRICTLY UPPER TRIANGULAR PART OF
!           A IS NOT REFERENCED.
!           NOTE THAT WHEN  DIAG = 'U' OR 'U', THE DIAGONAL ELEMENTS OF
!           A ARE NOT REFERENCED EITHER, BUT ARE ASSUMED TO BE UNITY.
!           UNCHANGED ON EXIT.
!
!  LDA    - INTEGER.
!           ON ENTRY, LDA SPECIFIES THE FIRST DIMENSION OF A AS DECLARED
!           IN THE CALLING (SUB) PROGRAM. LDA MUST BE AT LEAST
!           MAX( 1, N ).
!           UNCHANGED ON EXIT.
!
!  X      - DOUBLE PRECISION ARRAY OF DIMENSION AT LEAST
!           ( 1 + ( N - 1 )*ABS( INCX ) ).
!           BEFORE ENTRY, THE INCREMENTED ARRAY X MUST CONTAIN THE N
!           ELEMENT VECTOR X. ON EXIT, X IS OVERWRITTEN WITH THE
!           TRANFORMED VECTOR X.
!
!  INCX   - INTEGER.
!           ON ENTRY, INCX SPECIFIES THE INCREMENT FOR THE ELEMENTS OF
!           X. INCX MUST NOT BE ZERO.
!           UNCHANGED ON EXIT.
!
!
!  LEVEL 2 BLAS ROUTINE.
!
!  -- WRITTEN ON 22-OCTOBER-1986.
!     JACK DONGARRA, ARGONNE NATIONAL LAB.
!     JEREMY DU CROZ, NAG CENTRAL OFFICE.
!     SVEN HAMMARLING, NAG CENTRAL OFFICE.
!     RICHARD HANSON, SANDIA NATIONAL LABS.
!
!
!     .. PARAMETERS ..
!     .. LOCAL SCALARS ..
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT PARAMETERS.
!
      info = 0 
      if (.not.lsame(uplo,'U') .and. .not.lsame(uplo,'L')) then 
        info = 1 
      else if (.not.lsame(trans,'N') .and. .not.lsame(trans,'T') .and. .not.&
          lsame(trans,'C')) then 
        info = 2 
      else if (.not.lsame(diag,'U') .and. .not.lsame(diag,'N')) then 
        info = 3 
      else if (n < 0) then 
        info = 4 
      else if (lda < max(1,n)) then 
        info = 6 
      else if (incx == 0) then 
        info = 8 
      endif 
      if (info /= 0) then 
        call xerbla ('DTRMV ', info) 
        return  
      endif 
!
!     QUICK RETURN IF POSSIBLE.
!
      if (n == 0) return  
!
      nounit = lsame(diag,'N') 
!
!     SET UP THE START POINT IN X IF THE INCREMENT IS NOT UNITY. THIS
!     WILL BE  ( N - 1 )*INCX  TOO SMALL FOR DESCENDING LOOPS.
!
      if (incx <= 0) then 
        kx = 1 - (n - 1)*incx 
      else if (incx /= 1) then 
        kx = 1 
      endif 
!
!     START THE OPERATIONS. IN THIS VERSION THE ELEMENTS OF A ARE
!     ACCESSED SEQUENTIALLY WITH ONE PASS THROUGH A.
!
      if (lsame(trans,'N')) then 
!
!        FORM  X := A*X.
!
        if (lsame(uplo,'U')) then 
          if (incx == 1) then 
            do j = 1, n 
              if (x(j) == zero) cycle  
              temp = x(j) 
              x(:j-1) = x(:j-1) + temp*a(:j-1,j) 
              if (.not.nounit) cycle  
              x(j) = x(j)*a(j,j) 
            end do 
          else 
            jx = kx 
            do j = 1, n 
              if (x(jx) /= zero) then 
                temp = x(jx) 
                ix = kx 
                x(ix:(j-2)*incx+ix:incx) = x(ix:(j-2)*incx+ix:incx) + temp*a(:j&
                  -1,j) 
                if (nounit) x(jx) = x(jx)*a(j,j) 
              endif 
              jx = jx + incx 
            end do 
          endif 
        else 
          if (incx == 1) then 
            do j = n, 1, -1 
              if (x(j) == zero) cycle  
              temp = x(j) 
              x(n:1+j:(-1)) = x(n:1+j:(-1)) + temp*a(n:1+j:(-1),j) 
              if (.not.nounit) cycle  
              x(j) = x(j)*a(j,j) 
            end do 
          else 
            kx = kx + (n - 1)*incx 
            jx = kx 
            do j = n, 1, -1 
              if (x(jx) /= zero) then 
                temp = x(jx) 
                ix = kx 
                x(ix:incx*(j-n+1)+ix:(-incx)) = x(ix:incx*(j-n+1)+ix:(-incx))&
                   + temp*a(n:1+j:(-1),j) 
                if (nounit) x(jx) = x(jx)*a(j,j) 
              endif 
              jx = jx - incx 
            end do 
          endif 
        endif 
      else 
!
!        FORM  X := A'*X.
!
        if (lsame(uplo,'U')) then 
          if (incx == 1) then 
            do j = n, 1, -1 
              temp = x(j) 
              if (nounit) temp = temp*a(j,j) 
              temp = temp + sum(a(j-1:1:(-1),j)*x(j-1:1:(-1))) 
              x(j) = temp 
            end do 
          else 
            jx = kx + (n - 1)*incx 
            do j = n, 1, -1 
              temp = x(jx) 
              ix = jx 
              if (nounit) temp = temp*a(j,j) 
              temp = temp + sum(a(j-1:1:(-1),j)*x(ix-incx:incx*(1-j)+ix:(-incx)&
                )) 
              x(jx) = temp 
              jx = jx - incx 
            end do 
          endif 
        else 
          if (incx == 1) then 
            do j = 1, n 
              temp = x(j) 
              if (nounit) temp = temp*a(j,j) 
              temp = temp + sum(a(j+1:n,j)*x(j+1:n)) 
              x(j) = temp 
            end do 
          else 
            jx = kx 
            do j = 1, n 
              temp = x(jx) 
              ix = jx 
              if (nounit) temp = temp*a(j,j) 
              temp = temp + sum(a(j+1:n,j)*x(ix+incx:(n-j)*incx+ix:incx)) 
              x(jx) = temp 
              jx = jx + incx 
            end do 
          endif 
        endif 
      endif 
!
      return  
!
!     END OF DTRMV .
!
      end subroutine dtrmv 
