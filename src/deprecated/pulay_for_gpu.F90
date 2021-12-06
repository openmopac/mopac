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

      subroutine pulay_for_gpu(f, p, n, fppf, fock, emat, &
                        & lfock, nfock, msize, start, pl) 
      use chanel_C, only : iw
      use molkst_C, only : numcal, keywrd, mpack
#ifdef GPU
      Use mult_symm_ab_I 
      Use mod_vars_cuda, only: real_cuda, prec, ngpus
#endif      
      Use mod_vars_cuda, only: lgpu 
      implicit none
      integer  :: n,iopc 
      integer , intent(inout) :: lfock 
      integer , intent(inout) :: nfock 
      integer , intent(in) :: msize 
      double precision , intent(out) :: pl 
      logical , intent(inout) :: start 
      double precision  :: f(mpack) 
      double precision  :: p(mpack) 
      double precision  :: fppf(6*mpack)
      double precision , intent(inout) :: fock(*) 
      double precision , intent(inout) :: emat(20,20) 
      integer :: icalcn, maxlim, linear, mfock, lbase, i, nfock1, j, l, il, ii, idim, kdim
      double precision, dimension(1000) :: evec 
      double precision, dimension(20) :: coeffs 
      double precision :: const, d, sum 
      logical :: debug  
#ifdef GPU
      integer :: igrid, iblock
#endif        
      double precision, external :: ddot       
      
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
!                  EMAT   = WORKSTORE OF AT LEAST 20**2 ELEMENTS.
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
        maxlim = 6
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
      ! TODO: use dcopy
      fock(lfock:(linear-1)*mfock+lfock:mfock) = f(:linear) 
!
!   NOW FORM /FOCK*DENSITY-DENSITY*FOCK/, AND STORE THIS IN FPPF
!      
      idim = linear + lbase
      
      if (lgpu) then
        iopc = 4
      else
        iopc = 3
      end if

      call mult_symm_AB(p, f, 1.d0, n, linear, fppf(lbase+1:idim), 0.d0, iopc)
      call mult_symm_AB(f, p, 1.d0, n, linear, fppf(lbase+1:idim), -1.d0, iopc)
      
!      call mamult (p, f, fppf(lbase+1), n, 0.D0) 
!      call mamult (f, p, fppf(lbase+1), n, -1.D0) 

!
!   FPPF NOW CONTAINS THE RESULT OF FP - PF.
!
      nfock1 = nfock + 1 
      do i = 1, nfock
        kdim = i*linear        
        emat(nfock1,i) = -1.D0 
        emat(i,nfock1) = -1.D0  
        emat(lfock,i) = ddot(linear,fppf((i-1)*linear+1:kdim),1,fppf(lbase+1:idim),1)
        
        emat(i,lfock) = emat(lfock,i) 
      end do 
      pl = emat(lfock,lfock)/linear 
      emat(nfock1,nfock1) = 0.D0 
      if (emat(lfock, lfock) < 1.d-20) return
      const = 1.D0/emat(lfock,lfock) 
      
      emat(:nfock,:nfock) = emat(:nfock,:nfock)*const 
      ! TODO:check call mkl_dimatcopy('c', 't', 20, 20, const, emat, nfock, nfock)
      
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
      ! TODO:check call mkl_dimatcopy('c', 't', 20, 20, const, emat, nfock, nfock)
      
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
        write (iw, '("    LAGRANGIAN MULTIPLIER (ERROR) =", F13.6)') evec(nfock1*nfock1) 
      end if 
     
     ! TODO: make it parallel
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
      end subroutine pulay_for_gpu
      
