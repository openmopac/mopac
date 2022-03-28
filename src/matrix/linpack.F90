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

subroutine dgedi (a, lda, n, ipvt, det, work, job)
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: job, lda, n
    integer, dimension (n), intent (in) :: ipvt
    double precision, dimension (2), intent (out) :: det
    double precision, dimension (n), intent (out) :: work
    double precision, dimension (lda, n), intent (inout) :: a
   !
   !.. Local Scalars ..
    integer :: i, j, k, kb, kp1, l, nm1
    double precision :: t, ten
   !
   !.. External Calls ..
   !
   !.. Intrinsic Functions ..
    intrinsic Abs, Mod
   !
   ! ... Executable Statements ...
   !
   !     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
   !     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
   !
   !     ON ENTRY
   !
   !        A       DOUBLE PRECISION(LDA, N)
   !                THE OUTPUT FROM DGECO OR DGEFA.
   !
   !        LDA     INTEGER
   !                THE LEADING DIMENSION OF THE ARRAY  A .
   !
   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .
   !
   !        IPVT    INTEGER(N)
   !                THE PIVOT VECTOR FROM DGECO OR DGEFA.
   !
   !        WORK    DOUBLE PRECISION(N)
   !                WORK VECTOR.  CONTENTS DESTROYED.
   !
   !        JOB     INTEGER
   !                = 11   BOTH DETERMINANT AND INVERSE.
   !                = 01   INVERSE ONLY.
   !                = 10   DETERMINANT ONLY.
   !
   !     ON RETURN
   !
   !        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
   !                OTHERWISE UNCHANGED.
   !
   !        DET     DOUBLE PRECISION(2)
   !                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
   !                OTHERWISE NOT REFERENCED.
   !                DETERMINANT = DET(1) * 10.0**DET(2)
   !                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
   !                OR  DET(1) .EQ. 0.0 .
   !
   !     ERROR CONDITION
   !
   !        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
   !        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
   !        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
   !        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
   !        INFO .EQ. 0 .
   !
   !     LINPACK. THIS VERSION DATED 08/14/78 .
   !     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
   !
   !     SUBROUTINES AND FUNCTIONS
   !
   !     BLAS DAXPY,DSCAL,DSWAP
   !     FORTRAN ABS,MOD
   !
   !     COMPUTE DETERMINANT
   !
    if (job/10 /= 0) then
      det(1) = 1.0d0
      det(2) = 0.0d0
      ten = 10.0d0
      do i = 1, n
        if (ipvt(i) /= i) then
          det(1) = -det(1)
        end if
        det(1) = a(i, i) * det(1)
         !        ...EXIT
        if (det(1) == 0.0d0) exit
        do while (Abs (det( 1)) < 1.0d0)
          det(1) = ten * det(1)
          det(2) = det(2) - 1.0d0
        end do
        do while (Abs (det( 1)) >= ten)
          det(1) = det(1) / ten
          det(2) = det(2) + 1.0d0
        end do
      end do
    end if
   !
   !     COMPUTE INVERSE(U)
   !
    if (Mod(job, 10) == 0) return
    do k = 1, n
      a(k, k) = 1.0d0 / a(k, k)
      t = -a(k, k)
      call dscal (k-1, t, a(1, k), 1)
      kp1 = k + 1
      if (n >= kp1) then
        do j = kp1, n
          t = a(k, j)
          a(k, j) = 0.0d0
          call daxpy (k, t, a(1, k), 1, a(1, j), 1)
        end do
      end if
    end do
   !
   !        FORM INVERSE(U)*INVERSE(L)
   !
    nm1 = n - 1
    if (nm1 < 1) return
    do kb = 1, nm1
      k = n - kb
      kp1 = k + 1
      do i = kp1, n
        work(i) = a(i, k)
        a(i, k) = 0.0d0
      end do
      do j = kp1, n
        t = work(j)
        call daxpy (n, t, a(1, j), 1, a(1, k), 1)
      end do
      l = ipvt(k)
      if (l /= k) then
        call dswap (n, a(1, k), 1, a(1, l), 1)
      end if
    end do
end subroutine dgedi

subroutine dgefa (a, lda, n, ipvt, info)
    implicit none
    integer, intent (in) :: lda, n
    integer, intent (out) :: info
    integer, dimension (n), intent (out) :: ipvt
    double precision, dimension (lda, n), intent (inout) :: a
    integer :: j, k, kp1, l, nm1
    double precision :: t
      integer, external :: idamax
!

   !
   ! ... Executable Statements ...
   !
   !     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
   !
   !     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
   !     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
   !     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
   !
   !     ON ENTRY
   !
   !        A       DOUBLE PRECISION(LDA, N)
   !                THE MATRIX TO BE FACTORED.
   !
   !        LDA     INTEGER
   !                THE LEADING DIMENSION OF THE ARRAY  A .
   !
   !        N       INTEGER
   !                THE ORDER OF THE MATRIX  A .
   !
   !     ON RETURN
   !
   !        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
   !                WHICH WERE USED TO OBTAIN IT.
   !                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
   !                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
   !                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
   !
   !        IPVT    INTEGER(N)
   !                AN INTEGER VECTOR OF PIVOT INDICES.
   !
   !        INFO    INTEGER
   !                = 0  NORMAL VALUE.
   !                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
   !                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
   !                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
   !                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
   !                     INDICATION OF SINGULARITY.
   !
   !     LINPACK. THIS VERSION DATED 08/14/78 .
   !     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
   !
   !     SUBROUTINES AND FUNCTIONS
   !
   !     BLAS DAXPY,DSCAL,IDAMAX
   !
   !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
   !
    info = 0
    nm1 = n - 1
    if (nm1 >= 1) then
      do k = 1, nm1
        kp1 = k + 1
         !
         !        FIND L = PIVOT INDEX
         !
        l = idamax (n-k+1, a(k, k), 1) + k - 1
        ipvt(k) = l
         !
         !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
         !
        if (a(l, k) == 0.0d0) then
          info = k
        else
            !
            !           INTERCHANGE IF NECESSARY
            !
          if (l /= k) then
            t = a(l, k)
            a(l, k) = a(k, k)
            a(k, k) = t
          end if
            !
            !           COMPUTE MULTIPLIERS
            !
          t = -1.0d0 / a(k, k)
          call dscal (n-k, t, a(k+1, k), 1)
            !
            !           ROW ELIMINATION WITH COLUMN INDEXING
            !
          do j = kp1, n
            t = a(l, j)
            if (l /= k) then
              a(l, j) = a(k, j)
              a(k, j) = t
            end if
            call daxpy (n-k, t, a(k+1, k), 1, a(k+1, j), 1)
          end do
        end if
      end do
    end if
    ipvt(n) = n
    if (a(n, n) == 0.0d0) then
      info = n
    end if
end subroutine dgefa

