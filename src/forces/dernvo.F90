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

      subroutine dernvo()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numcal, numat, nopen, nclose, norbs, gnorm, &
      & fract, keywrd, mpack, moperr, lm61, n2elec
      use meci_C, only : nbo, nmos, nelec, lab
      use common_arrays_C, only : eigs, dxyz
      use chanel_C, only : iw
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
      double precision, dimension(:), allocatable  :: scalar, diag, fmooff, &
      work
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  minear, ninear, nvax, icalcn, i, ilast, ifirst, j, k, l, ll, &
      n1, n2
      double precision, dimension(nmos*norbs + 1) :: fmoon
      double precision, dimension(3*numat) :: dxyzr
      double precision, dimension(nmos*norbs*20) :: eigbb
      double precision :: sum, throld, sumx, sumy, sumz
      logical :: debug, dcar, large, relaxd
      character :: blank*60

      save  debug, dcar, large,  minear, ninear, nvax, icalcn, &
      & diag, scalar, fmooff
!-----------------------------------------------
!**********************************************************************
!
!    IMPLEMENTATION OF ANALYTICAL FORMULATION FOR OPEN SHELL OR CI,
!                      VARIABLES FINITE DIFFERENCE METHODS,
!                      STATISTICAL ESTIMATE OF THE ERRORS,
!                   BY D. LIOTARD
!                      LABORATOIRE DE CHIMIE STRUCTURALE
!                      UNIVERSITE DE PAU ET DES PAYS DE L'ADOUR
!                      AVENUE DE L'UNIVERSITE, 64000, PAU (FRANCE)
!
!
!   MODIFIED BY JJPS TO CONFORM TO MOPAC CONVENTIONS
!   (NOTE BY JJPS:  PROF. LIOTARD'S TECHNIQUE WORKS.  IF THIS
!   IMPLEMENTATION DOES NOT WORK, THE REASON IS A FAULT INTRODUCED
!   BY JJPS, AND DOES NOT REFLECT ON PROF. LIOTARD'S ABILITY)
!
!
!    AS THE WAVE FUNCTION IS NOT VARIATIONALLY OPTIMIZED, I.E.
!    HALF-ELECTRON OR CI, THE DERIVATIVES OF THE 1 AND 2-ELECTRON
!    INTEGRALS IN A.O. BASIS ARE EVALUATED IN CARTESIAN COORDINATES
!    BY A 1 OR 2 POINTS FINITE DIFFERENCE FORMULA AND STORED.
!    THUS ONE GETS THE NON-RELAXED (I.E. FROZEN ELECTRONIC CLOUD)
!    CONTRIBUTION TO THE FOCK EIGENVALUES AND 2-ELECTRON INTEGRALS IN
!    AN M.O. BASIS.  THE NON-RELAXED GRADIENT COMES FROM THE
!    NON-RELAXED C.I. MATRIX DERIVATIVE (SUBROUTINE DERI1).
!    THE DERIVATIVES OF THE M.O. COEFFICIENTS ARE THEN WORKED OUT
!    ITERATIVELY (OK FOR BOTH CLOSED SHELLS AND HALF-ELECTRON CASES)
!    AND STORED. THUS ONE GETS THE ELECTRONIC RELAXATION CONTRIBUTION TO
!    THE FOCK EIGENVALUES AND 2-ELECTRON INTEGRALS IN M.O. BASIS.
!    FINALLY THE RELAXATION CONTRIBUTION TO THE C.I. MATRIX DERIVATIVE
!    GIVES THE RELAXATION CONTRIBUTION TO THE GRADIENT (ROUTINE DERI2).
!
!
!        COORD  HOLDS THE CARTESIAN COORDINATES.
!    INPUT
!        DXYZ   NOT DEFINED.
!    EXIT
!        DXYZ   DERIVATIVES OF ENERGY W.R.T CARTESIAN COORDINATES,
!               IN KCAL/MOL/ANGSTROM (3 * NUMAT OF THESE)
!
!**********************************************************************
      data icalcn/ 0/
!        ACTUAL SIZES FOR C.I. GRADIENT CALCULATION.
        nbo(1) = nclose
        nbo(2) = nopen - nclose
        nbo(3) = norbs - nopen
        minear = nbo(2)*nbo(1) + nbo(3)*nopen
        ninear = nmos
        k = 0
        l = 0
        n2 = 0
        do i = 1, 3
          l = k + 1
          k = k + nbo(i)
          n1 = Max (0, nelec - l)
          n2 = Min (k, nelec + nmos) - l
          if (n2 > n1) then
            ninear = ninear + (n2*(n2 + 1) - n1*(n1 + 1))/2
          end if
        end do
        j = n2 + 1
        n2 = n2 + l
        if (j > 0 .and. n2 < norbs) then
          ninear = ninear + j * (norbs-n2)
        else if (ninear <= 0) then
          ninear = 1
        end if
        i = max(norbs**2 + 45*lm61, (lab*(lab+1))/2, n2elec + mpack)
        allocate(work(i))
        if(allocated(scalar)) deallocate(scalar)
        allocate(scalar(minear))
        if(allocated(diag)) deallocate(diag)
        allocate(diag(minear))
        if(allocated(fmooff)) deallocate(fmooff)
        allocate(fmooff(2*minear))
!
!     SELECT THE REQUIRED OPTION AND READ KEYWORDS
!     --------------------------------------------
!
      if (icalcn /= numcal) then
        icalcn = numcal
        debug = index(keywrd,'DERNVO') /= 0
        large = index(keywrd,'LARGE') /= 0
        dcar = index(keywrd,'FORC') + index(keywrd,'PREC') /= 0
        nvax = 3*numat
      end if
      dxyzr(:nvax) = 0.D0
!        SCALING ROW FACTORS TO SPEED CV OF RELAXATION PROCEDURE.
      call deri0 (eigs, norbs, scalar, diag, fract, nbo)
!
!   BECAUSE DERI2 IS CPU INTENSIVE, AND THE CONTRIBUTION TO THE
!   DERIVATIVE DUE TO RELAXATION OF THE ELECTRON CLOUD IS RELATIVELY
!   INSENSITIVE TO CHANGES IN GEOMETRY, WHERE POSSIBLE ONLY CALCULATE
!   THE DERIVATIVE EVERY 2 CALLS TO DERNVO
!
      eigbb = 0.D0
      sum = 0.D0
      if (gnorm<1.D0 .or. dcar) then
        dxyzr(:nvax) = 0.D0
        relaxd = .FALSE.
      end if
      do i = 1, nvax
        sum = sum + abs(dxyzr(i))
      end do
      relaxd = sum > 1.D-7
!
!  IF DXYZR CONTAINS DATA, USE IT AND FLUSH AFTER USE.
!
      ilast = 0
   50 continue
      ifirst = ilast + 1
      ilast = min(nvax,ilast + 1)
      j = 1 - minear
      k = 1 - ninear
      do i = ifirst, ilast
        k = k + ninear
        j = j + minear
!
!        NON-RELAXED CONTRIBUTION (FROZEN ELECTRONIC CLOUD) IN DXYZ
!        AND NON-RELAXED FOCK MATRICES IN FMOOFF AND FMOON.
!   CONTENTS OF F-MO-OFF: OPEN-CLOSED, VIRTUAL-CLOSED, AND VIRTUAL-OPEN
!   CONTENTS OF F-MO-ON:  CLOSED-CLOSED, OPEN-OPEN AND VIRTUAL-VIRTUAL
!   OVER M.O. INDICES
!
        call deri1 (i, dxyz(i), fmooff(j), minear, fmoon(k), scalar, work)
      end do
      if (debug) then
        if (ifirst==1 .and. large) then
          write (iw, *) ' CONTENTS OF FMOOFF '
          write (iw, *) ' OPEN-CLOSED'
          write (iw, '(7X,I3,5I12)') (j,j=nclose + 1,nopen)
          do i = 1, nclose
            write (iw, '(I3,6F12.6)') i, (fmooff(j),j=(i - 1)*nbo(2) + 1,i*nbo(2))
          end do
!
!
          write (iw, *) ' VIRTUAL-CLOSED'
          k = nclose*nbo(2)
          write (iw, '(7X,I3,5I12)') (j,j=nopen + 1,min(nopen + 6,norbs))
          do i = 1, nclose
            write (iw, '(I3,6F12.6)') i, &
            (fmooff(j+k),j=(i - 1)*nbo(3) + 1,min(6 + (i - 1)*nbo(3),i*nbo(3)))
          end do
          k = nclose*nbo(2) + nbo(3)*nclose
!
!
          write (iw, *) ' VIRTUAL-OPEN'
          write (iw, '(7X,I3,4I12)') (j,j=nclose + 1,nopen)
          do i = 1, min(6,nbo(3))
            write (iw, '(I3,6F12.6)') i + nopen, &
            (fmooff(j+k),j=(i - 1)*nbo(2) + 1,min((i - 1)*nbo(2)+6,i*nbo(2)))
          end do
          write (iw, *) ' CONTENTS OF FMOON (ACTIVE-SPACE -- ACTIVE SPACE)'
          k = (nmos*(nmos - 1))/2
          ll = 1
          blank = ' '
          do i = 1, nmos
            l = ll + nmos - i - 1
            write (iw, '(A,5F12.6)') blank(:12*i), (fmoon(j),j=ll,l), fmoon(k+i)
            ll = l + 1
          end do
        end if
      end if
!        COMPUTE THE ELECTRONIC RELAXATION CONTRIBUTION.
!
!   DERNVO PROVIDES THE FOLLOWING SCRATCH AREAS TO DERI2: EIGB, WORK2,
!          WORK3, CBETA.  THESE ARE DIMENSIONED ON ENTRY TO DERI2
!          WHICH IS WHY THEY ARE NOT DECLARED THERE.  THEY ARE NOT USED
!          AT ALL IN DERNVO.
!
!
!  The following function was chosen as a guide to THROLD.  It is NOT
!  intended to be hard-and-fast.
!
!     throld = max(0.001D0,min(thref,gnorm**3*0.00002D0))
      throld = 0.00001D0
      if (.not.relaxd) call deri2 (minear, fmooff, fmoon, eigbb &
        , ninear, ilast - ifirst + 1, dxyzr(ifirst) &
        , throld,  diag, scalar, work)
      if (moperr) goto 99
      if (ilast < nvax) go to 50
      if (debug) then
        sumx = 0.D0
        sumy = 0.D0
        sumz = 0.D0
        do i = 1, numat
          sumx = sumx + dxyz(i*3-2)
          sumy = sumy + dxyz(i*3-1)
          sumz = sumz + dxyz(i*3)
        end do
        write (iw, *) ' CARTESIAN DERIVATIVES DUE TO FROZEN CORE'
        write (iw, '('' ATOM    X           Y           Z'')')
        do i = 1, numat
          write (iw, '(I4,3F12.7)') i, dxyz(i*3-2), dxyz(i*3-1), dxyz(i*3)
        end do
        write (iw, '(/10X,''RESIDUAL ERROR'')')
        write (iw, '(4X,3F12.7)') sumx, sumy, sumz
        write (iw, *)
        sumx = 0.D0
        sumy = 0.D0
        sumz = 0.D0
        do i = 1, numat
          sumx = sumx + dxyzr(i*3-2)
          sumy = sumy + dxyzr(i*3-1)
          sumz = sumz + dxyzr(i*3)
        end do
        write (iw, *) ' CARTESIAN DERIVATIVES DUE TO RELAXING CORE'
        write (iw, '('' ATOM    X           Y           Z'')')
        do i = 1, numat
          write (iw, '(I4,3F12.7)') i, dxyzr(i*3-2), dxyzr(i*3-1), dxyzr(i*3)
        end do
        write (iw, '(/10X,''RESIDUAL ERROR'')')
        write (iw, '(4X,3F12.7)') sumx, sumy, sumz
        write (iw, *)
      end if
      dxyz(:nvax) = dxyz(:nvax) + dxyzr(:nvax)
      if (relaxd) then
        dxyzr(:nvax) = 0.D0
      end if
      sumx = 0.D0
      sumy = 0.D0
      sumz = 0.D0
      do i = 1, numat
        sumx = sumx + dxyz(i*3-2)
        sumy = sumy + dxyz(i*3-1)
        sumz = sumz + dxyz(i*3)
      end do
      sum = max(1.D-10,abs(sumx) + abs(sumy) + abs(sumz))
      if (debug) then
        write (iw, *) 'CARTESIAN DERIVATIVES FROM ANALYTICAL C.I. CALCULATION'
        write (iw, '('' ATOM    X           Y           Z'')')
        do i = 1, numat
          write (iw, '(I4,3F12.7)') i, dxyz(i*3-2), dxyz(i*3-1), dxyz(i*3)
        end do
        write (iw, '(/10X,''RESIDUAL ERROR'')')
        write (iw, '(4X,3F12.7)') sumx, sumy, sumz
        write (iw, *)
      end if
  99  continue
      if(allocated(scalar)) deallocate(scalar)
      if(allocated(diag))   deallocate(diag)
      if(allocated(fmooff)) deallocate(fmooff)
      return
      end subroutine dernvo
