! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

! HEigensystem.F
! diagonalization of a Hermitian n-by-n matrix using the Jacobi algorithm
! code adapted from the "Handbook" routines for complex A
! (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
! this file is part of the Diag library
! last modified 9 Aug 11 th
!************************************************************************
! HEigensystem diagonalizes a Hermitian n-by-n matrix.
!* Input: n, A = n-by-n matrix, Hermitian
!* (only the upper triangle of A needs to be filled).
!* Output: d = vector of eigenvalues, U = transformation matrix
!* these fulfill diag(d) = U A U^+ = U A U^-1 with U unitary.

	subroutine HEigensystem(n, A, d, U, sort)
	implicit none
	integer n, sort
	Complex A(n,n), U(n,n)
	double precision d(n)

	integer p, q, j
	double precision sym_eps
	Real red, off, thresh
	Real t, delta, invc, s
	Complex x, y, Apq
	Real ev(2,n)

	integer sweep
	common /nsweeps/ sweep

	Real Sq
	Complex c
	Sq(c) = real(c*conjg(c))
	SYM_EPS = 2D0**(-103)

	do p = 1, n
	  ev(1,p) = 0
	  ev(2,p) = real(A(p,p))
	  d(p) = ev(2,p)
	enddo

	do p = 1, n
	  do q = 1, n
	    U(q,p) = (0.0, 0.0)
	  enddo
	  U(p,p) = (1.0, 0.0)
	enddo

	red = 0.04/n**4

	do sweep = 1, 50
	  off = 0
	  do q = 2, n
	    do p = 1, q - 1
	      off = off + Sq(A(p,q))
	    enddo
	  enddo
	  if( .not. off .gt. SYM_EPS ) goto 1

	  thresh = 0
	  if( sweep .lt. 4 ) thresh = off*red

	  do q = 2, n
	    do p = 1, q - 1
	      Apq = A(p,q)
	      off = Sq(Apq)
        if(sweep > 4 .and. off < SYM_EPS*(ev(2,p)**2+ev(2,q)**2) ) then
	        A(p,q) = (0.0, 0.0)
	      else if( off .gt. thresh ) then
	        t = 0.5*(ev(2,p) - ev(2,q))
	        t = 1/(t + sign(sqrt(t**2 + off), t))

	        delta = t*off
	        ev(1,p) = ev(1,p) + delta
	        ev(2,p) = real(d(p)) + ev(1,p)
	        ev(1,q) = ev(1,q) - delta
	        ev(2,q) = real(d(q)) + ev(1,q)

	        invc = sqrt(delta*t + 1)
	        s = t/invc
	        t = delta/(invc + 1)

	        do j = 1, p - 1
	          x = A(j,p)
	          y = A(j,q)
	          A(j,p) = x + s*(conjg(Apq)*y - t*x)
	          A(j,q) = y - s*(Apq*x + t*y)
	        enddo

	        do j = p + 1, q - 1
	          x = A(p,j)
	          y = A(j,q)
	          A(p,j) = x + s*(Apq*conjg(y) - t*x)
	          A(j,q) = y - s*(Apq*conjg(x) + t*y)
	        enddo

	        do j = q + 1, n
	          x = A(p,j)
	          y = A(q,j)
	          A(p,j) = x + s*(Apq*y - t*x)
	          A(q,j) = y - s*(conjg(Apq)*x + t*y)
	        enddo

	        A(p,q) = (0.0, 0.0)

	        do j = 1, n
	          x = U(p,j)
	          y = U(q,j)
	          U(p,j) = x + s*(Apq*y - t*x)
	          U(q,j) = y - s*(conjg(Apq)*x + t*y)
	        enddo
	      endif
	    enddo
	  enddo

	  do p = 1, n
	    ev(1,p) = 0
	    d(p) = ev(2,p)
	  enddo
	enddo

	print *, "Bad convergence in HEigensystem"

1	if( sort .eq. 0 ) return

! sort the eigenvalues

	do p = 1, n - 1
	  j = p
	  t = real(d(p))
	  do q = p + 1, n
	    if( sort*(t - d(q)) .gt. 0 ) then
	      j = q
	      t = real(d(q))
	    endif
	  enddo

	  if( j .ne. p ) then
	    d(j) = d(p)
	    d(p) = t
	    do q = 1, n
	      x = U(p,q)
	      U(p,q) = U(j,q)
	      U(j,q) = x
	    enddo
	  endif
	enddo
	do p = 1, n
	  do q = 1, p
	    x = U(p,q)
	    U(p,q) = U(q,p)
	    U(q,p) = x
	  end do
	end do
	end

