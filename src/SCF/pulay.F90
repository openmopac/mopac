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

      subroutine pulay(f, p, n, fppf, fock, emat, lfock, nfock, msize, start, pl)
      use chanel_C, only : iw
      use molkst_C, only : numcal, keywrd, mpack, npulay
      use Common_arrays_C, only : workmat1, workmat2, workmat3
      implicit none
      integer  :: n
      integer , intent(inout) :: lfock
      integer , intent(inout) :: nfock
      integer , intent(in) :: msize
      double precision , intent(out) :: pl
      logical , intent(inout) :: start
      double precision  :: f(mpack)
      double precision  :: p(mpack)
      double precision  :: fppf(*)
      double precision , intent(inout) :: fock(*)
      double precision , intent(inout) :: emat(npulay+1,npulay+1)
!
      integer :: icalcn, maxlim, linear, mfock, lbase, i, nfock1, j, l, il, ii
      double precision, dimension((npulay+1)**2) :: evec
      double precision, dimension(npulay) :: coeffs
      double precision :: const, d, sum
      logical :: debug
      double precision, external :: ddot
!
      save icalcn, maxlim, debug, linear, mfock
!-----------------------------------------------
!***********************************************************************
!
!   PULAY USES DR. PETER PULAY'S METHOD FOR CONVERGENCE.
!         A MATHEMATICAL DESCRIPTION CAN BE FOUND IN
!         "P. PULAY, J. COMP. CHEM. 3, 556 (1982).
!
! ARGUMENTS:-
!         ON INPUT F      = FOCK MATRIX, PACKED, LOWER HALF TRIANGLE.
!                  P      = DENSITY MATRIX, PACKED, LOWER HALF TRIANGLE.
!                  N      = NUMBER OF ORBITALS.
!                  FPPF   = WORKSTORE OF SIZE MSIZE, CONTENTS WILL BE
!                           OVERWRITTEN.
!                  FOCK   =      "       "              "         "
!                  EMAT   = WORKSTORE OF AT LEAST NPULAY**2 ELEMENTS.
!                  START  = LOGICAL, = TRUE TO START PULAY.
!                  PL     = UNDEFINED ELEMENT.
!      ON OUTPUT   F      = "BEST" FOCK MATRIX, = LINEAR COMBINATION
!                           OF KNOWN FOCK MATRICES.
!                  START  = FALSE
!                  PL     = MEASURE OF NON-SELF-CONSISTENCY
!                         = [F*P] = F*P - P*F.
!
!***********************************************************************
      data icalcn/ 0/
      if (icalcn /= numcal) then
        icalcn = numcal
        maxlim = npulay
        debug = index(keywrd,'DEBUGPULAY') /= 0
      end if
      if (start) then
        linear = (n*(n + 1))/2
        mfock = msize/linear
        mfock = min0(maxlim,mfock)
        if (debug) write (iw, '('' MAXIMUM SIZE:'',I5)') mfock
        nfock = 1
        lfock = 1
        start = .FALSE.
      else
        if (nfock < mfock) nfock = nfock + 1
        if (lfock /= mfock) then
          lfock = lfock + 1
        else
          lfock = 1
        end if
      end if
      lbase = (lfock - 1)*linear
!
!   FIRST, STORE FOCK MATRIX FOR FUTURE REFERENCE.
!
      fock(lfock:(linear-1)*mfock+lfock:mfock) = f(:linear)
!
!   NOW FORM /FOCK*DENSITY-DENSITY*FOCK/, AND STORE THIS IN FPPF
!
!      call mamult (p, f, fppf(lbase+1), n, 0.D0)
!      call mamult (f, p, fppf(lbase+1), n, -1.D0)
      call unpack_matrix(p, workmat1, n)
      call unpack_matrix(f, workmat2, n)
      call sym_commute(workmat1, workmat2, workmat3, n)
      call pack_matrix(workmat3, fppf(lbase+1), n)
!
!   FPPF NOW CONTAINS THE RESULT OF FP - PF.
!
      nfock1 = nfock + 1
      do i = 1, nfock
        emat(nfock1,i) = -1.D0
        emat(i,nfock1) = -1.D0
        emat(lfock,i) = ddot(linear,fppf((i-1)*linear+1),1,fppf(lbase+1),1)
        emat(i,lfock) = emat(lfock,i)
      end do
      pl = emat(lfock,lfock)/linear
      emat(nfock1,nfock1) = 0.D0
      if (emat(lfock, lfock) < 1.d-20) return
      const = 1.D0/emat(lfock,lfock)
      emat(:nfock,:nfock) = emat(:nfock,:nfock)*const
      if (debug) then
        write (iw, '('' EMAT'')')
        do i = 1, nfock1
          write (iw, '(6E13.6)') (emat(j,i),j=1,nfock1)
        end do
      end if
      l = 0
      do i = 1, nfock1
        evec(l+1:nfock1+l) = emat(i,:nfock1)
        l = nfock1 + l
      end do
      const = 1.D0/const
      emat(:nfock,:nfock) = emat(:nfock,:nfock)*const
!********************************************************************
!   THE MATRIX EMAT SHOULD HAVE FORM
!
!      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
!      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
!      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
!      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
!      |     .            .      ...     . |
!      |   -1.0         -1.0     ...    0. |
!
!   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
!   TIMES [F*P] FOR ITERATION J.
!
!********************************************************************
      call osinv (evec, nfock1, d)
      if (abs(d) < 1.D-6) then
        start = .TRUE.
        return
      end if
      if (nfock < 2) return
      il = nfock*nfock1
      coeffs(:nfock) = -evec(1+il:nfock+il)
      if (debug) then
        write (iw, '('' EVEC'')')
        write (iw, '(6F12.6)') (coeffs(i),i=1,nfock)
        write (iw, '(''    LAGRANGIAN MULTIPLIER (ERROR) =''                          ,F13.6)') evec(nfock1*nfock1)
      end if
      do i = 1, linear
        sum = 0.D0
        l = 0
        ii = (i - 1)*mfock
        do j = 1, nfock
          sum = sum + coeffs(j)*fock(j+ii)
        end do
        f(i) = sum
      end do
      return
      end subroutine pulay

      subroutine pack_matrix(unpacked, packed, size)
        implicit none
        integer :: info
        integer , intent(in) :: size
        double precision , intent(in) :: unpacked(size, size)
        double precision , intent(out) :: packed(*)
  !-----------------------------------------------
  !***********************************************************************
  !
  !   CONVERT UNPACKED SYMMETRIC MATRIX INTO A PACKED UPPER TRIANGLE
  !   (LAPACK DTRTTP CALL)
  !
  ! ARGUMENTS:-
  !         ON INPUT UNPACKED = UNPACKED SYMMETRIC MATRIX
  !                  SIZE     = DIMENSION OF MATRIX
  !      ON OUTPUT   PACKED   = PACKED UPPER TRIANGLE MATRIX
  !
  !***********************************************************************
        call dtrttp( 'U', size, unpacked, size, packed, info )
        if (info /= 0) stop 'error in dtrttp'
        return
        end subroutine pack_matrix

        subroutine unpack_matrix(packed, unpacked, size)
          implicit none
          integer :: info, i, j
          integer , intent(in) :: size
          double precision , intent(in) :: packed(*)
          double precision , intent(out) :: unpacked(size, size)
    !-----------------------------------------------
    !***********************************************************************
    !
    !   CONVERT PACKED UPPER TRIANGLE INTO AN UNPACKED SYMMETRIC MATRIX
    !   (LAPACK DTPTTR CALLS & FILLING IN THE REST BY HAND)
    !
    ! ARGUMENTS:-
    !         ON INPUT PACKED   = PACKED UPPER TRIANGLE MATRIX
    !                  SIZE     = DIMENSION OF MATRIX
    !      ON OUTPUT   UNPACKED = UNPACKED SYMMETRIC MATRIX
    !
    !***********************************************************************
          call dtpttr( 'U', size, packed, unpacked, size, info )
          if (info /= 0) stop 'error in dtpttr'
          do i = 1, size
            do j = i+1, size
              unpacked(j,i) = unpacked(i,j)
            end do
          end do
          return
          end subroutine unpack_matrix

          subroutine sym_commute(mat1, mat2, mat3, size)
            implicit none
            integer :: i, j
            integer , intent(in) :: size
            double precision , intent(in) :: mat1(size, size)
            double precision , intent(in) :: mat2(size, size)
            double precision , intent(out) :: mat3(size, size)
    !-----------------------------------------------
    !***********************************************************************
    !
    !   COMPUTE THE COMMUTATOR BETWEEN TWO SYMMETRIC MATRICES
    !   (DGEMM & AN IN-PLACE EVALUATION OF THE SECOND TERM)
    !
    ! ARGUMENTS:-
    !         ON INPUT MAT1   = FIRST SYMMETRIC MATRIX IN COMMUTATOR
    !                  MAT2   = SECOND SYMMETRIC MATRIX IN COMMUTATOR
    !                  SIZE   = DIMENSION OF MATRIX
    !      ON OUTPUT   MAT3   = MAT1*MAT2 - MAT2*MAT1
    !
    !***********************************************************************
            call dsymm('L', 'U', size, size, 1.D0, mat1, size, mat2, size, 0.D0, mat3, size)
            do i = 1, size
              do j = i, size
                mat3(i,j) = mat3(i,j) - mat3(j,i)
                mat3(j,i) = -mat3(i,j)
              end do
            end do
            return
            end subroutine sym_commute
