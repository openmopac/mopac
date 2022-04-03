! Molecular Orbital PACkage (MOPAC)
! Copyright (C) 2021, Virginia Polytechnic Institute and State University
!
! MOPAC is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! MOPAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

      subroutine rsp(a, n, nvect, root, vect) 
      USE vast_kind_param, ONLY:  double 
      implicit none
      integer  :: n 
      integer  :: nvect 
      real(double)  :: a(*) 
      real(double), intent (out)  :: root(n) 
      real(double)  :: vect(n,n) 
!
      integer :: i 
!
! Trivial case: n = 1
!
      if (n == 1) then
        root(1) = a(1)
        vect(1,1) = 1.d0
        return
      end if
!
! Perturb secular determinant matrix to split exact degeneracies
! (This is to get around a bug in the diagonalizer that causes eigenvectors to not be orthonormal)
!
      do i = 1, n                                   ! Do NOT go much higher than 1.d-10, otherwise the geometry 
        a((i*(i+1))/2) = a((i*(i+1))/2) + i*1.d-10  ! optimization might go into an endless loop.
      end do
      call evvrsp ((-1), n, nvect, n*n, n, a, root, vect, 0, i)
      return  
      end subroutine rsp 


      subroutine evvrsp(msgfl, n, nvect, lena, nv, a, root, vect, &
        iorder, ierr) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: msgfl 
      integer  :: n 
      integer  :: nvect 
      integer  :: lena 
      integer  :: nv 
      integer , intent(in) :: iorder 
      integer  :: ierr 
      integer  :: ind(n) 
      real(double)  :: a(lena) 
      real(double)  :: b(n,8) 
      real(double)  :: root(n) 
      real(double)  :: vect(nv,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: lmsgfl, i, j, l, jsv, klim, k 
      real(double) :: t 
!-----------------------------------------------
!***********************************************************************
!*       ROUTINE EVVRSP (MSGFL,N,NVECT,LENA,NV,A,B,IND,ROOT,VECT,IORDER
!*   *                 ,IERR)
!*
!*    AUTHOR:  S. T. ELBERT, AMES LABORATORY-USDOE, JUNE 1985
!*
!*    PURPOSE -
!*       FINDS   (ALL) EIGENVALUES    AND    (SOME OR ALL) EIGENVECTORS
!*                     *    *                                   *
!*       OF A REAL SYMMETRIC PACKED MATRIX.
!*            *    *         *
!*
!*    METHOD -
!*       THE METHOD AS PRESENTED IN THIS ROUTINE CONSISTS OF FOUR STEPS:
!*       FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
!*       HOUSEHOLDER TECHNIQUE (ORTHOGONAL SIMILARITY TRANSFORMATIONS).
!*       SECOND, THE ROOTS ARE LOCATED USING THE RATIONAL QL METHOD.
!*       THIRD, THE VECTORS OF THE TRIDIAGONAL FORM ARE EVALUATED BY THE
!*       INVERSE ITERATION TECHNIQUE.  VECTORS FOR DEGENERATE OR NEAR-
!*       DEGENERATE ROOTS ARE FORCED TO BE ORTHOGONAL.
!*       FOURTH, THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
!*       ORIGINAL ARRAY.
!*
!*       THESE ROUTINES ARE MODIFICATIONS OF THE EISPACK 3
!*       ROUTINES TRED3, TQLRAT, TINVIT AND TRBAK3
!*
!*       FOR FURTHER DETAILS, SEE EISPACK USERS GUIDE, B. T. SMITH
!*       ET AL, SPRINGER-VERLAG, LECTURE NOTES IN COMPUTER SCIENCE,
!*       VOL. 6, 2-ND EDITION, 1976.  ANOTHER GOOD REFERENCE IS
!*       THE SYMMETRIC EIGENVALUE PROBLEM BY B. N. PARLETT
!*       PUBLISHED BY PRENTICE-HALL, INC., ENGLEWOOD CLIFFS, N.J. (1980)
!*
!*    ON ENTRY -
!*       MSGFL  - INTEGER (LOGICAL UNIT NO.)
!*                FILE WHERE ERROR MESSAGES WILL BE PRINTED.
!*                IF MSGFL IS 0, ERROR MESSAGES WILL BE PRINTED ON LU 6.
!*                IF MSGFL IS NEGATIVE, NO ERROR MESSAGES PRINTED.
!*       N      - INTEGER
!*                ORDER OF MATRIX A.
!*       NVECT  - INTEGER
!*                NUMBER OF VECTORS DESIRED.  0 .LE. NVECT .LE. N.
!*       LENA   - INTEGER
!*                DIMENSION OF  A  IN CALLING ROUTINE.  MUST NOT BE LESS
!*                THAN (N*N+N)/2.
!*       NV     - INTEGER
!*                ROW DIMENSION OF VECT IN CALLING ROUTINE.   N .LE. NV.
!*       A      - WORKING PRECISION REAL (LENA)
!*                INPUT MATRIX, ROWS OF THE LOWER TRIANGLE PACKED INTO
!*                LINEAR ARRAY OF DIMENSION N*(N+1)/2.  THE PACKED ORDER
!*                IS A(1,1), A(2,1), A(2,2), A(3,1), A(3,2), ...
!*       B      - WORKING PRECISION REAL (N,8)
!*                SCRATCH ARRAY, 8*N ELEMENTS
!*       IND    - INTEGER (N)
!*                SCRATCH ARRAY OF LENGTH N.
!*       IORDER - INTEGER
!*                ROOT ORDERING FLAG.
!*                = 0, ROOTS WILL BE PUT IN ASCENDING ORDER.
!*                = 2, ROOTS WILL BE PUT IN DESCENDING ORDER.
!*
!*    ON EXIT -
!*       A      - DESTORYED.  NOW HOLDS REFLECTION OPERATORS.
!*       ROOT   - WORKING PRECISION REAL (N)
!*                ALL EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
!*                  IF IORDER = 0, ROOT(1) .LE. ... .LE. ROOT(N)
!*                  IF IORDER = 2, ROOT(1) .GE. ... .GE. ROOT(N)
!*       VECT   - WORKING PRECISION REAL (NV,NVECT)
!*                EIGENVECTORS FOR ROOT(1), ..., ROOT(NVECT).
!*       IERR   - INTEGER
!*                = 0 IF NO ERROR DETECTED,
!*                = K IF ITERATION FOR K-TH EIGENVALUE FAILED,
!*                = -K IF ITERATION FOR K-TH EIGENVECTOR FAILED.
!*                (FAILURES SHOULD BE VERY RARE.  CONTACT C. MOLER.)
!*
!*    EXTERNAL ROUTINES -
!*       ETRED3(ELAU, FREDA, DASUM, DNRM2, DSCAL)
!*       EQLRAT(EPSLON)
!*       EINVIT(EPSLON, ESTPI1, DASUM, DAXPY, DDOT, DNRM2, DSCAL)
!*       ETRBK3(DDOT, DAXPY)
!*
!***********************************************************************
!
!                              DECLARATIONS
!
!
!
!
!
      lmsgfl = msgfl 
      if (msgfl == 0) lmsgfl = 6 
      ierr = n - 1 
      if (n <= 0) go to 150 
      ierr = n + 1 
      if ((n*n + n)/2 > lena) go to 160 
!
!        REDUCE REAL SYMMETRIC MATRIX A TO TRIDIAGONAL FORM
!
      call etred3 (n, a, b(1,1), b(1,2), b(1,3)) 
!
!        FIND ALL EIGENVALUES OF TRIDIAGONAL MATRIX
!
      call eqlrat (n, b(1,1), b(1,2), b(1,3), root, ind, ierr, b(1,4)) 
      if (ierr /= 0) go to 170 
!
!         CHECK THE DESIRED ORDER OF THE EIGENVALUES
!
      b(1,3) = dble(iorder) 
      if (iorder /= 0) then 
        if (iorder /= 2) go to 200 
!
!         ORDER ROOTS IN DESCENDING ORDER (LARGEST FIRST)...
!        TURN ROOT AND IND ARRAYS END FOR END
!
        do i = 1, n/2 
          j = n + 1 - i 
          t = root(i) 
          root(i) = root(j) 
          root(j) = t 
          l = ind(i) 
          ind(i) = ind(j) 
          ind(j) = l 
        end do 
!
!           FIND I AND J MARKING THE START AND END OF A SEQUENCE
!           OF DEGENERATE ROOTS
!
        i = 0 
   90   continue 
        i = i + 1 
        if (i > n) go to 130 
        do j = i, n 
          if (Abs(root(j) - root(i)) > 1.d-20) go to 110 
        end do 
        j = n + 1 
  110   continue 
        j = j - 1 
        if (j == i) go to 90 
!
!                    TURN AROUND IND BETWEEN I AND J
!
        jsv = j 
        klim = (j - i + 1)/2 
        do k = 1, klim 
          l = ind(j) 
          ind(j) = ind(i) 
          ind(i) = l 
          i = i + 1 
          j = j - 1 
        end do 
        i = jsv 
        go to 90 
!
      endif 
  130 continue 
      if (nvect <= 0) return  
      if (nv < n) go to 180 
!
!        FIND EIGENVECTORS OF TRI-DIAGONAL MATRIX VIA INVERSE ITERATION
!
      ierr = lmsgfl 
      call einvit (nv, n, b(1,1), b(1,2), b(1,3), nvect, root, ind, vect, ierr&
        , b(1,4), b(1,5), b(1,6), b(1,7), b(1,8)) 
      if (ierr /= 0) go to 190 
!
!        FIND EIGENVECTORS OF SYMMETRIC MATRIX VIA BACK TRANSFORMATION
!
  140 continue 
      call etrbk3 (nv, n, lena, a, nvect, vect) 
      return  
!
!        ERROR MESSAGE SECTION
!
  150 continue 
      if (lmsgfl < 0) return  
      go to 210 
!
  160 continue 
      if (lmsgfl < 0) return  
      go to 210 
!
  170 continue 
      if (lmsgfl < 0) return  
      go to 210 
!
  180 continue 
      if (lmsgfl < 0) return  
      go to 210 
!
  190 continue 
      go to 140 
!
  200 continue 
      ierr = -1 
      if (lmsgfl < 0) return  
!
  210 continue 
      return  
      end subroutine evvrsp 


      subroutine freda(l, d, a, e) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:01  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l 
      real(double) , intent(in) :: d(l) 
      real(double) , intent(inout) :: a(*) 
      real(double) , intent(in) :: e(l) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jk, j
      real(double) :: f, g 
!-----------------------------------------------
!
!
      jk = 1 
!
!     .......... FORM REDUCED A ..........
!
      do j = 1, l 
        f = d(j) 
        g = e(j) 
!
        if (j > 0) then 
          a(jk:j-1+jk) = a(jk:j-1+jk) - f*e(:j) - g*d(:j) 
          jk = j + jk 
        endif 
!
      end do 
      return  
      end subroutine freda 


!***********************************************************************
!*       ROUTINE EINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
!*
!*    AUTHORS-
!*       THIS IS A MODIFICATION OF ROUTINE TINVIT FROM EISPACK EDITION 3
!*       DATED AUGUST 1983.
!*       TINVIT IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE
!*       IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
!*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)
!*
!*    PURPOSE -
!*       THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
!*       SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES.
!*
!*    METHOD -
!*       INVERSE ITERATION.
!*
!*    ON ENTRY -
!*       NM     - INTEGER
!*                MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!*                ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE
!*                DIMENSION STATEMENT.
!*       N      - INTEGER
!*       D      - W.P. REAL (N)
!*                CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!*       E      - W.P. REAL (N)
!*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!*                IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
!*       E2     - W.P. REAL (N)
!*                CONTAINS THE SQUARES OF CORRESPONDING ELEMENTS OF E,
!*                WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
!*                E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
!*                THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE
!*                SUM OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST
!*                CONTAIN 0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER,
!*                OR 2.0 IF THE EIGENVALUES ARE IN DESCENDING ORDER.
!*                IF TQLRAT, BISECT, TRIDIB, OR IMTQLV
!*                HAS BEEN USED TO FIND THE EIGENVALUES, THEIR
!*                OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
!*       M      - INTEGER
!*                THE NUMBER OF SPECIFIED EIGENVECTORS.
!*       W      - W.P. REAL (M)
!*                CONTAINS THE M EIGENVALUES IN ASCENDING
!*                OR DESCENDING ORDER.
!*       IND    - INTEGER (M)
!*                CONTAINS IN FIRST M POSITIONS THE SUBMATRIX INDICES
!*                ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!*                1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX
!*                FROM THE TOP, 2 FOR THOSE BELONGING TO THE SECOND
!*                SUBMATRIX, ETC.
!*       IERR   - INTEGER (LOGICAL UNIT NUMBER)
!*                LOGICAL UNIT FOR ERROR MESSAGES
!*
!*    ON EXIT -
!*       ALL INPUT ARRAYS ARE UNALTERED.
!*       Z      - W.P. REAL (NM,M)
!*                CONTAINS THE ASSOCIATED SET OF ORTHONORMAL
!*                EIGENVECTORS. ANY VECTOR WHICH WHICH FAILS TO CONVERGE
!*                IS LEFT AS IS (BUT NORMALIZED) WHEN ITERATING STOPPED.
!*       IERR   - INTEGER
!*                SET TO
!*                ZERO    FOR NORMAL RETURN,
!*                -R      IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
!*                        EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
!*                        (ONLY LAST FAILURE TO CONVERGE IS REPORTED)
!*
!*       RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
!*
!*       RV1    - W.P. REAL (N)
!*                DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!*       RV2    - W.P. REAL (N)
!*                SUPER(1)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!*       RV3    - W.P. REAL (N)
!*                SUPER(2)-DIAGONAL ELEMENTS OF U FROM LU DECOMPOSITION
!*       RV4    - W.P. REAL (N)
!*                ELEMENTS DEFINING L IN LU DECOMPOSITION
!*       RV6    - W.P. REAL (N)
!*                APPROXIMATE EIGENVECTOR
!*
!*    DIFFERENCES FROM EISPACK 3 -
!*       EPS3 IS SCALED BY  EPSCAL  (ENHANCES CONVERGENCE, BUT
!*          LOWERS ACCURACY)!
!*       ONE MORE ITERATION (MINIMUM 2) IS PERFORMED AFTER CONVERGENCE
!*          (ENHANCES ACCURACY)!
!*       REPLACE LOOP WITH PYTHAG WITH SINGLE CALL TO DNRM2!
!*       IF NOT CONVERGED, USE PERFORMANCE INDEX TO DECIDE ON ERROR
!*          VALUE SETTING, BUT DO NOT STOP!
!*       L.U. FOR ERROR MESSAGES PASSED THROUGH IERR
!*       USE PARAMETER STATEMENTS AND GENERIC INTRINSIC FUNCTIONS
!*       USE LEVEL 1 BLAS
!*       USE IF-THEN-ELSE TO CLARIFY LOGIC
!*       LOOP OVER SUBSPACES MADE INTO DO LOOP.
!*       LOOP OVER INVERSE ITERATIONS MADE INTO DO LOOP
!*       ZERO ONLY REQUIRED PORTIONS OF OUTPUT VECTOR
!*
!*    EXTERNAL ROUTINES -
!*       EPSLON
!*       BLAS(1)--DASUM, DAXPY, DDOT, DNRM2, DSCAL
!*       INTRINSIC FUNCTIONS - ABS, MAX, SQRT
!*
!*    NOTE -
!*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
!*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
!*
!***********************************************************************
      subroutine einvit(nm, n, d, e, e2, m, w, ind, z, ierr, rv1, rv2, rv3, rv4&
        , rv6) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nm 
      integer , intent(in) :: n 
      integer , intent(in) :: m 
      integer , intent(out) :: ierr 
      integer , intent(in) :: ind(m) 
      real(double)  :: d(n) 
      real(double)  :: e(n + 1) 
      real(double) , intent(in) :: e2(n) 
      real(double) , intent(in) :: w(m) 
      real(double)  :: z(nm,m) 
      real(double) , intent(inout) :: rv1(n) 
      real(double) , intent(inout) :: rv2(n) 
      real(double) , intent(inout) :: rv3(n) 
      real(double) , intent(inout) :: rv4(n) 
      real(double)  :: rv6(n) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: epscal = 0.5D+00 
      real(double), parameter :: grptol = 0.001D+00 
      real(double), parameter :: hundrd = 100.0D+00 
      real(double), parameter :: one = 1.0D+00 
      real(double), parameter :: zero = 0.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: group, i, its, j, jj, p, q, r, s, submat, tag 
      real(double) :: anorm, eps2, eps3, eps4, norm, order, rho, u, uk, v, x0, &
        x1, xu 
      double precision, external :: ddot, epslon, dasum, dnrm2, estpi1
      logical :: convgd 
!-----------------------------------------------
!
!                        DECLARATIONS
!
!-----------------------------------------------------------------------
!
      ierr = 0 
      x0 = zero 
      uk = zero 
      norm = zero 
      eps2 = zero 
      eps3 = zero 
      eps4 = zero 
      group = 0 
      tag = 0 
      order = one - e2(1) 
      q = 0 
      do submat = 1, n 
        p = q + 1 
!
!        .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
!
        do q = p, n - 1 
          if (e2(q+1) == zero) go to 30 
        end do 
        q = n 
!
!        .......... FIND VECTORS BY INVERSE ITERATION ..........
!
   30   continue 
        tag = tag + 1 
        anorm = zero 
        s = 0 
!
        do r = 1, m 
          if (ind(r) /= tag) cycle  
          its = 1 
          x1 = w(r) 
          if (s == 0) then 
!
!        .......... CHECK FOR ISOLATED ROOT ..........
!
            xu = one 
            if (p == q) then 
              rv6(p) = one 
              convgd = .TRUE. 
              go to 160 
!
            endif 
            norm = abs(d(p)) 
            do i = p + 1, q 
              norm = max(norm,abs(d(i))+abs(e(i))) 
            end do 
!
!        .......... EPS2 IS THE CRITERION FOR GROUPING,
!                   EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                   ROOTS ARE MODIFIED BY EPS3,
!                   EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW .........
!
            eps2 = grptol*norm 
            eps3 = epscal*epslon(norm) 
            uk = q - p + 1 
            eps4 = uk*eps3 
            uk = eps4/sqrt(uk) 
            s = p 
            group = 0 
          else 
!
!        .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
!
            if (abs(x1 - x0) >= eps2) then 
!
!                 ROOTS ARE SEPERATE
!
              group = 0 
            else 
!
!                 ROOTS ARE CLOSE
!
              group = group + 1 
              if (order*(x1 - x0) <= eps3) x1 = x0 + order*eps3 
            endif 
!
!        .......... ELIMINATION WITH INTERCHANGES AND
!                   INITIALIZATION OF VECTOR ..........
!
          endif 
!
          u = d(p) - x1 
          v = e(p+1) 
          rv6(p) = uk 
          do i = p + 1, q 
            rv6(i) = uk 
            if (abs(e(i)) > abs(u)) then 
!
!                 EXCHANGE ROWS BEFORE ELIMINATION
!
!                  *** WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
!                      E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY .......
!
              xu = u/e(i) 
              rv4(i) = xu 
              rv1(i-1) = e(i) 
              rv2(i-1) = d(i) - x1 
              rv3(i-1) = e(i+1) 
              u = v - xu*rv2(i-1) 
              v = -xu*rv3(i-1) 
!
            else 
!
!                    STRAIGHT ELIMINATION
!
              xu = e(i)/u 
              rv4(i) = xu 
              rv1(i-1) = u 
              rv2(i-1) = v 
              rv3(i-1) = zero 
              u = d(i) - x1 - xu*v 
              v = e(i+1) 
            endif 
          end do 
!
          if (abs(u) <= eps3) u = eps3 
          rv1(q) = u 
          rv2(q) = zero 
          rv3(q) = zero 
!
!              DO INVERSE ITERATIONS
!
          convgd = .FALSE. 
          do its = 1, 5 
            if (its /= 1) then 
!
!                    .......... FORWARD SUBSTITUTION ..........
!
              if (norm == zero) then 
                rv6(s) = eps4 
                s = s + 1 
                if (s > q) s = p 
              else 
                xu = eps4/norm 
                call dscal (q - p + 1, xu, rv6(p), 1) 
              endif 
!
!                     ... ELIMINATION OPERATIONS ON NEXT VECTOR
!
              do i = p + 1, q 
                u = rv6(i) 
!
!                         IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
!                         WAS PERFORMED EARLIER IN THE
!                         TRIANGULARIZATION PROCESS ..........
!
                if (Abs(rv1(i-1) - e(i)) < 1.d-20) then 
                  u = rv6(i-1) 
                  rv6(i-1) = rv6(i) 
                else 
                  u = rv6(i) 
                endif 
                rv6(i) = u - rv4(i)*rv6(i-1) 
              end do 
            endif 
!
!           .......... BACK SUBSTITUTION
!
            rv6(q) = rv6(q)/rv1(q) 
            v = u 
            u = rv6(q) 
            norm = abs(u) 
            do i = q - 1, p, -1 
              rv6(i) = (rv6(i)-u*rv2(i)-v*rv3(i))/rv1(i) 
              v = u 
              u = rv6(i) 
              norm = norm + abs(u) 
            end do 
            if (group /= 0) then 
!
!                 ....... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                         MEMBERS OF GROUP ..........
!
              j = r 
              do jj = 1, group 
                j = j - 1 
                do while(ind(j) /= tag) 
                  j = j - 1 
                end do 
                call daxpy (q - p + 1, (-ddot(q - p + 1,rv6(p),1,z(p,j),1)),z(p&
                  ,j), 1, rv6(p), 1) 
              end do 
              norm = dasum(q - p + 1,rv6(p),1) 
            endif 
!
            if (convgd) exit  
            if (norm < one) cycle  
            convgd = .TRUE. 
          end do 
!
          xu = one/dnrm2(q - p + 1,rv6(p)) 
!
  160     continue 
          z(:p-1,r) = zero 
          z(p:q,r) = rv6(p:q)*xu 
          z(q+1:n,r) = zero 
!
          if (.not.convgd) then 
            rho = estpi1(q - p + 1,x1,d(p),e(p),z(p,r),anorm) 
!
!               *** SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
!
            if (rho > hundrd) ierr = -r 
          endif 
!
          x0 = x1 
        end do 
!
        if (q /= n) cycle  
        exit  
      end do 
      return  
      end subroutine einvit 


      subroutine elau(hinv, l, d, a, e) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:01  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: l 
      real(double) , intent(in) :: hinv 
      real(double) , intent(in) :: d(l) 
      real(double) , intent(in) :: a(*) 
      real(double) , intent(inout) :: e(l) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: half = 0.5D+00 
      real(double), parameter :: zero = 0.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jl, jk, j, jm1
      real(double) :: f, g, hh 
!-----------------------------------------------
!
!
!
      jl = l 
      e(1) = a(1)*d(1) 
      jk = 2 
      do j = 2, jl 
        f = d(j) 
        g = zero 
        jm1 = j - 1 
!
        if (jm1 > 0) then 
          g = g + sum(a(jk:jm1-1+jk)*d(:jm1)) 
          e(:jm1) = e(:jm1) + a(jk:jm1-1+jk)*f 
          jk = jm1 + jk 
        endif 
!
        e(j) = g + a(jk)*f 
        jk = jk + 1 
      end do 
!
!        .......... FORM P ..........
!
      f = zero 
      e = e*hinv 
      f = f + dot_product(e,d) 
!
!     .......... FORM Q ..........
!
      hh = f*half*hinv 
      e = e - hh*d 
!
      return  
      end subroutine elau 


!*MODULE EIGEN   *DECK EPSLON
!***********************************************************************
!*    FUNCT. ROUTINE EPSLON (X)
!*
!*    AUTHORS -
!*       THIS ROUTINE WAS TAKEN FROM EISPACK EDITION 3 DATED 4/6/83
!*       THIS VERSION IS BY S. T. ELBERT, AMES LABORATORY-USDOE NOV 1986
!*
!*    PURPOSE -
!*       ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
!*
!*    ON ENTRY -
!*       X      - WORKING PRECISION REAL
!*                VALUES TO FIND EPSLON FOR
!*
!*    ON EXIT -
!*       EPSLON - WORKING PRECISION REAL
!*                SMALLEST POSITIVE VALUE SUCH THAT X+EPSLON .NE. ZERO
!*
!*    QUALIFICATIONS -
!*       THIS ROUTINE SHOULD PERFORM PROPERLY ON ALL SYSTEMS
!*       SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
!*          1.  THE BASE USED IN REPRESENTING FLOATING POINT
!*              NUMBERS IS NOT A POWER OF THREE.
!*          2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
!*              THE ACCURACY USED IN FLOATING POINT VARIABLES
!*              THAT ARE STORED IN MEMORY.
!*       THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
!*       FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
!*       ASSUMPTION 2.
!*       UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
!*              A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
!*              B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
!*              C  IS NOT EXACTLY EQUAL TO ONE,
!*              EPS  MEASURES THE SEPARATION OF 1.0 FROM
!*                   THE NEXT LARGER FLOATING POINT NUMBER.
!*       THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
!*       ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
!*
!*    DIFFERENCES FROM EISPACK 3 -
!*       USE IS MADE OF PARAMETER STATEMENTS AND INTRINSIC FUNCTIONS
!*       --NO EXECUTEABLE CODE CHANGES--
!*
!*    EXTERNAL ROUTINES - NONE
!*    INTRINSIC FUNCTIONS - ABS
!*
!*    NOTE -
!*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
!*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
!*
!***********************************************************************
      real(kind(0.0d0)) function epslon (x) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:01  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(in) :: x 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+00 
      real(double), parameter :: one = 1.0D+00 
      real(double), parameter :: three = 3.0D+00 
      real(double), parameter :: four = 4.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double) :: a, b, c, eps 
!-----------------------------------------------
!
!                    DECLARATIONS
!
!
!
!-----------------------------------------------------------------------
!
      a = four/three 
      b = a - one 
      c = b + b + b 
      eps = abs(c - one) 
      do while(eps == zero) 
        b = a - one 
        c = b + b + b 
        eps = abs(c - one) 
      end do 
      epslon = eps*abs(x) 
      return  
      end function epslon 


!***********************************************************************
!*       ROUTINE EQLRAT(N,DIAG,E,E2IN,D,IND,IERR,E2)
!*
!*    AUTHORS -
!*       THIS IS A MODIFICATION OF ROUTINE EQLRAT FROM EISPACK EDITION 3
!*       DATED AUGUST 1983.
!*       TQLRAT IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
!*       ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)
!*
!*    PURPOSE -
!*       THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!*       TRIDIAGONAL MATRIX
!*
!*    METHOD -
!*       RATIONAL QL
!*
!*    ON ENTRY -
!*       N      - INTEGER
!*                THE ORDER OF THE MATRIX.
!*       D      - W.P. REAL (N)
!*                CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
!*       E2     - W.P. REAL (N)
!*                CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF
!*                THE INPUT MATRIX IN ITS LAST N-1 POSITIONS.
!*                E2(1) IS ARBITRARY.
!*
!*     ON EXIT -
!*       D      - W.P. REAL (N)
!*                CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!*                ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!*                ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!*                THE SMALLEST EIGENVALUES.
!*       E2     - W.P. REAL (N)
!*                DESTROYED.
!*       IERR   - INTEGER
!*                SET TO
!*                ZERO       FOR NORMAL RETURN,
!*                J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!*                           DETERMINED AFTER 30 ITERATIONS.
!*
!*    DIFFERENCES FROM EISPACK 3 -
!*       G=G+B INSTEAD OF IF(G.EQ.0) G=B ; B=B/4
!*       F77 BACKWARD LOOPS INSTEAD OF F66 CONSTRUCT
!*       GENERIC INTRINSIC FUNCTIONS
!*       ARRARY  IND  ADDED FOR USE BY EINVIT
!*
!*    EXTERNAL ROUTINES -
!*       EPSLON
!*       INTRINSIC--ABS, SIGN, SQRT
!*
!*    NOTE -
!*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
!*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
!*
!***********************************************************************
      subroutine eqlrat(n, diag, e, e2in, d, ind, ierr, e2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(out) :: ierr 
      integer , intent(out) :: ind(n) 
      real(double) , intent(in) :: diag(n) 
      real(double) , intent(in) :: e(n) 
      real(double) , intent(inout) :: e2in(n) 
      real(double) , intent(inout) :: d(n) 
      real(double) , intent(inout) :: e2(n) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: scale = 1.0D+00/64.0D+00 
      real(double), parameter :: zero = 0.0D+00 
      real(double), parameter :: one = 1.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, l, m, l1, k, itag 
      real(double) :: b, c, f, g, h, p, r, s, t 
      double precision, external :: epslon
!-----------------------------------------------
!
!
!
!
!-----------------------------------------------------------------------
      ierr = 0 
      d(1) = diag(1) 
      ind(1) = 1 
      k = 0 
      itag = 0 
      if (n /= 1) then 
!
        d(2:n) = diag(2:n) 
        e2(:n-1) = e2in(2:n) 
!
        f = zero 
        t = zero 
        b = epslon(one) 
        c = b*b 
        b = b*scale 
        e2(n) = zero 
!
        do l = 1, n 
          h = abs(d(l)) + abs(e(l)) 
          if (t < h) then 
            t = h 
            b = epslon(t) 
            c = b*b 
            b = b*scale 
          endif 
!     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
          m = l - 1 
          m = m + 1 
          do while(e2(m) > c) 
            m = m + 1 
          end do 
!     .......... E2(N) IS ALWAYS ZERO, SO THERE IS AN EXIT
!                FROM THE LOOP ..........
!
          if (m > k) then 
            if (m /= n) e2in(m+1) = zero 
            k = m 
            itag = itag + 1 
          endif 
          if (m == l) go to 80 
!
!           ITERATE
!
          do j = 1, 30 
!              .......... FORM SHIFT ..........
            l1 = l + 1 
            s = sqrt(e2(l)) 
            g = d(l) 
            p = (d(l1)-g)/(2.0D+00*s) 
            r = sqrt(p*p + 1.0D+00) 
            d(l) = s/(p + sign(r,p)) 
            h = g - d(l) 
!
            d(l1:n) = d(l1:n) - h 
!
            f = f + h 
!              .......... RATIONAL QL TRANSFORMATION ..........
            g = d(m) + b 
            h = g 
            s = zero 
            do i = m - 1, l, -1 
              p = g*h 
              r = p + e2(i) 
              e2(i+1) = s*r 
              s = e2(i)/r 
              d(i+1) = h + s*(h + d(i)) 
              g = d(i) - e2(i)/g + b 
              h = g*p/r 
              if (h == zero) then
              e2(l) = s * g
              d(l) = zero
              go to 80
            end if
            end do 
!
            e2(l) = s*g 
            d(l) = h 
!              .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST
            if (h == zero) go to 80 
            if (abs(e2(l)) <= abs(c/h)) go to 80 
            e2(l) = h*e2(l) 
            if (e2(l) == zero) go to 80 
          end do 
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
          ierr = l 
          exit  
!
!           CONVERGED
!
   80     continue 
          p = d(l) + f 
!           .......... ORDER EIGENVALUES ..........
          i = 1 
          if (l /= 1) then 
            if (p < d(1)) go to 100 
            i = l 
!           .......... LOOP TO FIND ORDERED POSITION
            i = i - 1 
            do while(p < d(i)) 
              i = i - 1 
            end do 
!
            i = i + 1 
            if (i == l) go to 120 
  100       continue 
            d(l:1+i:(-1)) = d(l-1:i:(-1)) 
            ind(l:1+i:(-1)) = ind(l-1:i:(-1)) 
!
          endif 
  120     continue 
          d(i) = p 
          ind(i) = itag 
        end do 
!
      endif 
      return  
      end subroutine eqlrat 


!*MODULE EIGEN   *DECK ETRBK3
!***********************************************************************
!*       ROUTINE ETRBK3(NM,N,NV,A,M,Z)
!*
!*    AUTHORS-
!*       THIS IS A MODIFICATION OF ROUTINE TRBAK3 FROM EISPACK EDITION 3
!*       DATED AUGUST 1983.
!*       EISPACK TRBAK3 IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
!*       NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!*       THIS VERSION IS BY S. T. ELBERT (AMES LABORATORY-USDOE)
!*
!*    PURPOSE -
!*       THIS ROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
!*       MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!*       SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  ETRED3.
!*
!*    METHOD -
!*       THE CALCULATION IS CARRIED OUT BY FORMING THE MATRIX PRODUCT
!*          Q*Z
!*       WHERE  Q  IS A PRODUCT OF THE ORTHOGONAL SYMMETRIC MATRICES
!*                Q = PROD(I)[1 - U(I)*.TRANSPOSE.U(I)*H(I)]
!*       U  IS THE AUGMENTED SUB-DIAGONAL ROWS OF  A  AND
!*       Z  IS THE SET OF EIGENVECTORS OF THE TRIDIAGONAL
!*       MATRIX  F  WHICH WAS FORMED FROM THE ORIGINAL SYMMETRIC
!*       MATRIX  C  BY THE SIMILARITY TRANSFORMATION
!*                F = Q(TRANSPOSE) C Q
!*       NOTE THAT ETRBK3 PRESERVES VECTOR EUCLIDEAN NORMS.
!*
!*
!*    COMPLEXITY -
!*       M*N**2
!*
!*    ON ENTRY-
!*       NM     - INTEGER
!*                MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!*                ARRAY PARAMETERS AS DECLARED IN THE CALLING ROUTINE
!*                DIMENSION STATEMENT.
!*       N      - INTEGER
!*                THE ORDER OF THE MATRIX  A.
!*       NV     - INTEGER
!*                MUST BE SET TO THE DIMENSION OF THE ARRAY  A  AS
!*                DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT.
!*       A      - W.P. REAL (NV)
!*                CONTAINS INFORMATION ABOUT THE ORTHOGONAL
!*                TRANSFORMATIONS USED IN THE REDUCTION BY  ETRED3  IN
!*                ITS FIRST  NV = N*(N+1)/2 POSITIONS.
!*       M      - INTEGER
!*                THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
!*       Z      - W.P REAL (NM,M)
!*                CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!*                IN ITS FIRST M COLUMNS.
!*
!*    ON EXIT-
!*       Z      - W.P. REAL (NM,M)
!*                CONTAINS THE TRANSFORMED EIGENVECTORS
!*                IN ITS FIRST M COLUMNS.
!*
!*    DIFFERENCES WITH EISPACK 3 -
!*       THE TWO INNER LOOPS ARE REPLACED BY DDOT AND DAXPY.
!*       MULTIPLICATION USED INSTEAD OF DIVISION TO FIND S.
!*       OUTER LOOP RANGE CHANGED FROM 2,N TO 3,N.
!*       ADDRESS POINTERS FOR  A  SIMPLIFIED.
!*
!*    EXTERNAL ROUTINES -
!*       BLAS(1)--DDOT, DAXPY
!*       INTRINSIC FUNCTIONS - NONE
!*
!*    NOTE -
!*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
!*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
!*
!***********************************************************************
      subroutine etrbk3(nm, n, nv, a, m, z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:01  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nm 
      integer , intent(in) :: n 
      integer , intent(in) :: nv 
      integer , intent(in) :: m 
      real(double)  :: a(nv) 
      real(double)  :: z(nm,m) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: zero = 0.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ii, im1, iz, j 
      real(double) :: h, s 
      double precision, external :: ddot
!-----------------------------------------------
!
!                        DECLARATIONS
!
!
!
!
!-----------------------------------------------------------------------
!
      if (m == 0) return  
      if (n <= 2) return  
!
      ii = 3 
      do i = 3, n 
        iz = ii + 1 
        ii = ii + i 
        h = a(ii) 
        if (h == zero) cycle  
        im1 = i - 1 
        do j = 1, m 
          s = -(ddot(im1,a(iz),1,z(1,j),1)*h)*h 
          call daxpy (im1, s, a(iz), 1, z(1,j), 1) 
        end do 
      end do 
      return  
      end subroutine etrbk3 


!*MODULE EIGEN   *DECK ETRED3
!***********************************************************************
!*       ROUTINE ETRED3(N,A,D,E,E2)
!*
!*    AUTHORS -
!*       THIS IS A MODIFICATION OF ROUTINE TRED3 FROM EISPACK EDITION 3
!*       DATED AUGUST 1983.
!*       EISPACK TRED3 IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
!*       NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!*       THIS VERSION IS BY S. T. ELBERT, AMES LABORATORY-USDOE JUN 1986
!*
!*    PURPOSE -
!*       THIS ROUTINE REDUCES A REAL SYMMETRIC (PACKED) MATRIX, STORED
!*       AS A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
!*       USING ORTHOGONAL SIMILARITY TRANSFORMATIONS, PRESERVING THE
!*       INFORMATION ABOUT THE TRANSFORMATIONS IN  A.
!*
!*    METHOD -
!*       THE TRIDIAGONAL REDUCTION IS PERFORMED IN THE FOLLOWING WAY.
!*       STARTING WITH J=N, THE ELEMENTS IN THE J-TH ROW TO THE
!*       LEFT OF THE DIAGONAL ARE FIRST SCALED, TO AVOID POSSIBLE
!*       UNDERFLOW IN THE TRANSFORMATION THAT MIGHT RESULT IN SEVERE
!*       DEPARTURE FROM ORTHOGONALITY.  THE SUM OF SQUARES  SIGMA  OF
!*       THESE SCALED ELEMENTS IS NEXT FORMED.  THEN, A VECTOR  U  AND
!*       A SCALAR
!*                      H = U(TRANSPOSE) * U / 2
!*       DEFINE A REFLECTION OPERATOR
!*                      P = I - U * U(TRANSPOSE) / H
!*       WHICH IS ORTHOGONAL AND SYMMETRIC AND FOR WHICH THE
!*       SIMILIARITY TRANSFORMATION  PAP  ELIMINATES THE ELEMENTS IN
!*       THE J-TH ROW OF  A  TO THE LEFT OF THE SUBDIAGONAL AND THE
!*       SYMMETRICAL ELEMENTS IN THE J-TH COLUMN.
!*
!*       THE NON-ZERO COMPONENTS OF  U  ARE THE ELEMENTS OF THE J-TH
!*       ROW TO THE LEFT OF THE DIAGONAL WITH THE LAST OF THEM
!*       AUGMENTED BY THE SQUARE ROOT OF  SIGMA  PREFIXED BY THE SIGN
!*       OF THE SUBDIAGONAL ELEMENT.  BY STORING THE TRANSFORMED SUB-
!*       DIAGONAL ELEMENT IN  E(J)  AND NOT OVERWRITING THE ROW
!*       ELEMENTS ELIMINATED IN THE TRANSFORMATION, FULL INFORMATION
!*       ABOUT  P  IS SAVE FOR LATER USE IN  ETRBK3.
!*
!*       THE TRANSFORMATION SETS  E2(J)  EQUAL TO  SIGMA  AND  E(J)
!*       EQUAL TO THE SQUARE ROOT OF  SIGMA  PREFIXED BY THE SIGN
!*       OF THE REPLACED SUBDIAGONAL ELEMENT.
!*
!*       THE ABOVE STEPS ARE REPEATED ON FURTHER ROWS OF THE
!*       TRANSFORMED  A  IN REVERSE ORDER UNTIL  A  IS REDUCED TO TRI-
!*       DIAGONAL FORM, THAT IS, REPEATED FOR  J = N-1,N-2,...,3.
!*
!*    COMPLEXITY -
!*       2/3 N**3
!*
!*    ON ENTRY-
!*       N      - INTEGER
!*                THE ORDER OF THE MATRIX.
!*       A      - W.P. REAL (NV)
!*                CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
!*                INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
!*                ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
!*
!*    ON EXIT-
!*       A      - W.P. REAL (NV)
!*                CONTAINS INFORMATION ABOUT THE ORTHOGONAL
!*                TRANSFORMATIONS USED IN THE REDUCTION.
!*       D      - W.P. REAL (N)
!*                CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL
!*                MATRIX.
!*       E      - W.P. REAL (N)
!*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!*                MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO
!*       E2     - W.P. REAL (N)
!*                CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF
!*                E. MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!*
!*    DIFFERENCES FROM EISPACK 3 -
!*       OUTER LOOP CHANGED FROM II=1,N TO I=N,3,-1
!*       PARAMETER STATEMENT AND GENERIC INTRINSIC FUNCTIONS USED
!*       SCALE.NE.0 TEST NOW SPOTS TRI-DIAGONAL FORM
!*       VALUES LESS THAN EPSLON CLEARED TO ZERO
!*       USE BLAS(1)
!*       U NOT COPIED TO D, LEFT IN A
!*       E2 COMPUTED FROM E
!*       INNER LOOPS SPLIT INTO ROUTINES ELAU AND FREDA
!*       INVERSE OF H STORED INSTEAD OF H
!*
!*    EXTERNAL ROUTINES -
!*       ELAU, FREDA
!*       BLAS(1)--DASUM, DNRM2, DSCAL
!*       INTRINSIC--ABS, MAX, SIGN, SQRT
!*
!*    NOTE -
!*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
!*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
!*
!***********************************************************************
      subroutine etred3(n, a, d, e, e2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double)  :: a(*) 
      real(double) , intent(out) :: d(*) 
      real(double)  :: e(*) 
      real(double) , intent(out) :: e2(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+00 
      real(double), parameter :: zero = 0.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, iia, iz0, l 
      real(double) :: aiimax, f, g, h, hroot, scale, scalei 
      double precision, external :: dasum, dnrm2
!-----------------------------------------------
!
!                        DECLARATIONS
!
!
!
!
!-----------------------------------------------------------------------
!
      if (n > 2) then 
        iz0 = (n*n + n)/2 
        aiimax = abs(a(iz0)) 
        do i = n, 3, -1 
          l = i - 1 
          iia = iz0 
          iz0 = iz0 - i 
          aiimax = max(aiimax,abs(a(iia))) 
          scale = dasum(l,a(iz0+1),1) 
          if (Abs(scale - abs(a(iia-1))) < 1.d-20 .or. Abs(aiimax+scale - aiimax) < 1.d-20) then 
!
!           THIS ROW IS ALREADY IN TRI-DIAGONAL FORM
!
            d(i) = a(iia) 
            if (Abs(aiimax + d(i) - aiimax) < 1.d-20) d(i) = zero 
            e(i) = a(iia-1) 
            if (Abs(aiimax + e(i) - aiimax) < 1.d-20) e(i) = zero 
            e2(i) = e(i)*e(i) 
            a(iia) = zero 
            cycle  
!
          endif 
!
          scalei = one/scale 
          call dscal (l, scalei, a(iz0+1), 1) 
          hroot = dnrm2(l,a(iz0+1)) 
!
          f = a(iz0+l) 
          g = -sign(hroot,f) 
          e(i) = scale*g 
          e2(i) = e(i)*e(i) 
          h = hroot*hroot - f*g 
          if (Abs(h) < 1.d-20) exit
          a(iz0+l) = f - g 
          d(i) = a(iia) 
          a(iia) = one/sqrt(h) 
!           .......... FORM P THEN Q IN E(1:L) ..........
          call elau (one/h, l, a(iz0+1), a, e) 
!           .......... FORM REDUCED A ..........
          call freda (l, a(iz0+1), a, e) 
!
        end do 
      endif 
      e(1) = zero 
      e2(1) = zero 
      d(1) = a(1) 
      e(2) = a(2) 
      e2(2) = a(2)*a(2) 
      d(2) = a(3) 
!
      return  
      end subroutine etred3 


!*MODULE EIGEN   *DECK ESTPI1
!***********************************************************************
!*    FUNCT. ROUTINE ESTPI1 (N,EVAL,D,E,X,ANORM)
!*
!*    AUTHOR -
!*       STEPHEN T. ELBERT (AMES LABORATORY-USDOE) DATE: 5 DEC 1986
!*
!*    PURPOSE -
!*       EVALUATE SYMMETRIC TRIDIAGONAL MATRIX PERFORMANCE INDEX
!*       *        *         *                  *           *
!*       FOR 1 EIGENVECTOR
!*           *
!*
!*    METHOD -
!*       THIS ROUTINE FORMS THE 1-NORM OF THE RESIDUAL MATRIX A*X-X*EVAL
!*       WHERE  A  IS A SYMMETRIC TRIDIAGONAL MATRIX STORED
!*       IN THE DIAGONAL (D) AND SUB-DIAGONAL (E) VECTORS, EVAL IS THE
!*       EIGENVALUE OF AN EIGENVECTOR OF  A,  NAMELY  X.
!*       THIS NORM IS SCALED BY MACHINE ACCURACY FOR THE PROBLEM SIZE.
!*       ALL NORMS APPEARING IN THE COMMENTS BELOW ARE 1-NORMS.
!*
!*    ON ENTRY -
!*       N      - INTEGER
!*                THE ORDER OF THE MATRIX  A.
!*       EVAL   - W.P. REAL
!*                THE EIGENVALUE CORRESPONDING TO VECTOR  X.
!*       D      - W.P. REAL (N)
!*                THE DIAGONAL VECTOR OF  A.
!*       E      - W.P. REAL (N)
!*                THE SUB-DIAGONAL VECTOR OF  A.
!*       X      - W.P. REAL (N)
!*                AN EIGENVECTOR OF  A.
!*       ANORM  - W.P. REAL
!*                THE NORM OF  A  IF IT HAS BEEN PREVIOUSLY COMPUTED.
!*
!*    ON EXIT -
!*       ANORM  - W.P. REAL
!*                THE NORM OF  A, COMPUTED IF INITIALLY ZERO.
!*       ESTPI1 - W.P. REAL
!*          !!A*X-X*EVAL!! / (EPSLON(10*N)*!!A!!*!!X!!);
!*          WHERE EPSLON(X) IS THE SMALLEST NUMBER SUCH THAT
!*             X + EPSLON(X) .NE. X
!*
!*          ESTPI1 .LT. 1 == SATISFACTORY PERFORMANCE
!*                 .GE. 1 AND .LE. 100 == MARGINAL PERFORMANCE
!*                 .GT. 100 == POOR PERFORMANCE
!*          (SEE LECT. NOTES IN COMP. SCI. VOL.6 PP 124-125)
!*
!***********************************************************************
      real(kind(0.0d0)) function estpi1 (n, eval, d, e, x, anorm) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      real(double) , intent(in) :: eval 
      real(double) , intent(inout) :: anorm 
      real(double) , intent(in) :: d(n) 
      real(double) , intent(in) :: e(n) 
      real(double) , intent(in) :: x(n) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(double), parameter :: one = 1.0D+00 
      real(double), parameter :: zero = 0.0D+00 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i 
      real(double) :: rnorm, size, xnorm 
      double precision, external :: epslon
!-----------------------------------------------
!
!                    DECLARATIONS
!
!
!
!-----------------------------------------------------------------------
!
      estpi1 = zero 
      if (n <= 1) return  
      size = 10*n 
      if (anorm == zero) then 
!
!              COMPUTE NORM OF  A
!
        anorm = max(abs(d(1))+abs(e(2)),abs(d(n))+abs(e(n))) 
        do i = 2, n - 1 
          anorm = max(anorm,abs(e(i))+abs(d(i))+abs(e(i+1))) 
        end do 
        if (anorm == zero) anorm = one 
      endif 
!
!           COMPUTE NORMS OF RESIDUAL AND EIGENVECTOR
!
      xnorm = abs(x(1)) + abs(x(n)) 
      rnorm = abs((d(1)-eval)*x(1)+e(2)*x(2)) + abs((d(n)-eval)*x(n)+e(n)*x(n-1&
        )) 
      do i = 2, n - 1 
        xnorm = xnorm + abs(x(i)) 
        rnorm = rnorm + abs(e(i)*x(i-1)+(d(i)-eval)*x(i)+e(i+1)*x(i+1)) 
      end do 
!
      estpi1 = rnorm/(epslon(size)*anorm*xnorm) 
      return  
      end function estpi1 
