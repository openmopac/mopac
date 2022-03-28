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

      subroutine deri1(number, grad, f, minear, fd, scalar, work)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE molkst_C, only : numat, nopen, numcal, keywrd, norbs, mpack, &
      & n2elec
      use chanel_C, only : iw
      use common_arrays_C, only : nfirst, nlast, eigb, p, pa, c, coord
      USE funcon_C, only : fpc_9
      use derivs_C, only : wmat, hmat, fmat
      use meci_C, only : nelec, nmos, lab, nbo, nstate, maxci, xy, vectci
      implicit none
      integer , intent(in) :: number
      integer , intent(in) :: minear
      double precision , intent(out) :: grad
      double precision  :: f(mpack)
      double precision  :: fd(mpack)
      double precision , intent(in) :: scalar(mpack)
      double precision, intent(out) :: work(norbs, norbs)
!
      integer ::  icalcn, nati, natx, i, l, j, nend, loop, ninit, n1&
        , n2, lcut, k, ll, ncol
      double precision :: const, step, enucl2, gse, sum
      logical :: debug
      double precision, external :: ddot, helect
      save debug, icalcn, const
!-----------------------------------------------
!********************************************************************
!
!     DERI1 COMPUTE THE NON-RELAXED DERIVATIVE OF THE NON-VARIATIONALLY
!     OPTIMIZED WAVEFUNCTION ENERGY WITH RESPECT TO ONE CARTESIAN
!     COORDINATE AT A TIME
!                             AND
!     COMPUTE THE NON-RELAXED FOCK MATRIX DERIVATIVE IN M.O BASIS AS
!     REQUIRED IN THE RELAXATION SECTION (ROUTINE 'DERI2').
!
!   INPUT
!     C(NORBS,NORBS) : M.O. COEFFICIENTS.
!     COORD  : CARTESIAN COORDINATES ARRAY.
!     NUMBER : LOCATION OF THE REQUIRED VARIABLE IN COORD.
!     WORK   : WORK ARRAY OF SIZE N*N.
!     WMAT     : WORK ARRAYS FOR d<PQ|RS> (2-CENTERS  A.O)
!   OUTPUT
!     C,COORD,NUMBER : NOT MODIFIED.
!     GRAD   : DERIVATIVE OF THE HEAT OF FORMATION WITH RESPECT TO
!              COORD(NUMBER), WITHOUT RELAXATION CORRECTION.
!     F(MINEAR) : NON-RELAXED FOCK MATRIX DERIVATIVE WITH RESPECT TO
!              COORD(NUMBER), EXPRESSED IN M.O BASIS, SCALED AND
!              PACKED, OFF-DIAGONAL BLOCKS ONLY.
!     FD     : IDEM BUT UNSCALED, DIAGONAL BLOCKS, C.I-ACTIVE ONLY.
!
!***********************************************************************
      data icalcn/ 0/
      n2 = 0
!
      if (icalcn /= numcal) then
        if(allocated(wmat)) deallocate(wmat)
        i = Max (n2elec + mpack, (lab*(lab+1))/2, 20*mpack, maxci*(maxci+1)/2)
        allocate(wmat(i))
        if(allocated(hmat)) deallocate(hmat)
        i = Max ((lab*(lab+1))/2, 20*minear, mpack)
        allocate(hmat(i))
        if(allocated(fmat)) deallocate(fmat)
        i = Max ((lab*(lab+1))/2, 20*minear, mpack)
        allocate(fmat(i))
        const = fpc_9
        debug = index(keywrd,'DERI1') /= 0
        icalcn = numcal
      end if
      if (debug) call timer ('BEFORE DERI1')
      step = 1.D-3
!
!     2 POINTS FINITE DIFFERENCE TO GET THE INTEGRAL DERIVATIVES
!     ----------------------------------------------------------
!     STORED IN HMAT AND WMAT, WITHOUT DIVIDING BY THE STEP SIZE.
!
      nati = (number - 1)/3 + 1
      natx = number - 3*(nati - 1)
      call dhcore (coord, hmat, wmat, enucl2, nati, natx, step)
!
! HMAT HOLDS THE ONE-ELECTRON DERIVATIVES OF ATOM NATI FOR DIRECTION
!      NATX W.R.T. ALL OTHER ATOMS
! WMAT HOLDS THE TWO-ELECTRON DERIVATIVES OF ATOM NATI FOR DIRECTION
!      NATX W.R.T. ALL OTHER ATOMS
      step = 0.5D0/step
!
!     NON-RELAXED FOCK MATRIX DERIVATIVE IN A.O BASIS.
!     ------------------------------------------------
!     STORED IN FMAT, DIVIDED BY STEP.
!
      call dcopy (mpack, hmat, 1, fmat, 1)
      call dfock2 (fmat, p, pa, wmat, numat, nfirst, nlast, nati)
!
!  FMAT HOLDS THE ONE PLUS TWO - ELECTRON DERIVATIVES OF ATOM NATI FOR
!       DIRECTION NATX W.R.T. ALL OTHER ATOMS
!
!       DERIVATIVE OF THE SCF-ONLY ENERGY (I.E BEFORE C.I CORRECTION)
!
      grad = (helect(norbs,p,hmat,fmat) + enucl2)*step
!     TAKE STEP INTO ACCOUNT IN FMAT
      fmat(:mpack) = fmat(:mpack)*step
!
!     RIGHT-HAND SIDE SUPER-VECTOR F = C' FMAT C USED IN RELAXATION
!     -----------------------------------------------------------
!     STORED IN NON-STANDARD PACKED FORM IN F(MINEAR) AND FD.
!     THE SUPERVECTOR IS THE NON-RELAXED FOCK MATRIX DERIVATIVE IN
!     M.O BASIS: F(IJ)= ( (C' * FOCK * C)(I,J) )   WITH I.GT.J .
!     F IS SCALED AND PACKED IN SUPERVECTOR FORM WITH
!                THE CONSECUTIVE FOLLOWING OFF-DIAGONAL BLOCKS:
!             1) OPEN-CLOSED  I.E. F(IJ)=F(I,J) WITH I OPEN & J CLOSED
!                                  AND I RUNNING FASTER THAN J,
!             2) VIRTUAL-CLOSED SAME RULE OF ORDERING,
!             3) VIRTUAL-OPEN   SAME RULE OF ORDERING.
!     FD IS PACKED OVER THE C.I-ACTIVE M.O WITH
!                THE CONSECUTIVE DIAGONAL BLOCKS:
!             1) CLOSED-CLOSED   IN CANONICAL ORDER, WITHOUT THE
!                                DIAGONAL ELEMENTS,
!             2) OPEN-OPEN       SAME RULE OF ORDERING,
!             3) VIRTUAL-VIRTUAL SAME RULE OF ORDERING.
!
!     PART 1 : WORK(N,N) = FMAT(N,N) * C(N,N)
      do i = 1, norbs
        call supdot (work(1,i), fmat, c(1,i), norbs)
      end do
!
!     PART 2 : F(IJ) =  (C' * WORK)(I,J) ... OFF-DIAGONAL BLOCKS.
      l = 1
      if (nbo(2)/=0 .and. nbo(1)/=0) then
!        OPEN-CLOSED
        ! TODO: GBR future modifications
        call mtxm (c(1,nbo(1)+1), nbo(2), work, norbs, f(l), nbo(1))
        l = l + nbo(2)*nbo(1)
      end if
      if (nbo(3)/=0 .and. nbo(1)/=0) then
!        VIRTUAL-CLOSED
        ! TODO: GBR future modifications
        call mtxm (c(1,nopen+1), nbo(3), work, norbs, f(l), nbo(1))
        l = l + nbo(3)*nbo(1)
      end if
!        VIRTUAL-OPEN
      ! TODO: GBR future modifications
      if (nbo(3)/=0 .and. nbo(2)/=0) &
      call mtxm (c(1,nopen+1), nbo(3), work(1,nbo(1)+1), norbs, f(l), nbo(2))
!     SCALE F ACCORDING TO THE DIAGONAL METRIC TENSOR 'SCALAR '.
      f(:minear) = f(:minear)*scalar(:minear)
      if (debug) then
        write (iw, *) ' F IN DERI1'
        j = min(20,minear)
        write (iw, '(5F12.6)') (f(i),i=1,j)
      end if
!
!     PART 3 : SUPER-VECTOR FD, C.I-ACTIVE DIAGONAL BLOCKS, UNSCALED.
      l = 1
      nend = 0
      do loop = 1, 3
        ninit = nend + 1
        nend = nend + nbo(loop)
        n1 = max(ninit,nelec + 1)
        n2 = min(nend,nelec + nmos)
        if (n2 < n1) cycle
        do i = n1, n2
          if (i <= ninit) cycle
          call mxm (c(1,i), 1, work(1,ninit), norbs, fd(l), i - ninit)
          l = l + i - ninit
        end do
      end do
      ncol = n2 - ninit + 1
      if (ncol > 0 .and. n2 < norbs) then
        ! TODO: GBR future modifications
        call mtxm (c(1, n2+1), norbs-n2, work(1,ninit), norbs, fd(l), ncol)
        l = l + ncol * (norbs-n2)
      end if
!
!     NON-RELAXED C.I CORRECTION TO THE ENERGY DERIVATIVE.
!     ----------------------------------------------------
!
!     C.I-ACTIVE FOCK EIGENVALUES DERIVATIVES, STORED IN FD(CONTINUED).
      lcut = l
      do i = nelec + 1, nelec + nmos
        fd(l) = ddot(norbs,c(1,i),1,work(1,i),1)
        l = l + 1
      end do
!
!     C.I-ACTIVE 2-ELECTRONS INTEGRALS DERIVATIVES. STORED IN XY.
!   FMAT IS USED HERE AS SCRATCH SPACE
!
      call dijkl1 (c(1,nelec+1), norbs, nati, wmat, fmat, hmat, fmat, xy)
      xy(:nmos,:nmos,:nmos,:nmos) = xy(:nmos,:nmos,:nmos,:nmos)*step
!
!     BUILD THE C.I MATRIX DERIVATIVE, STORED IN WMAT.
      call mecid (fd(lcut), gse, eigb, work, xy)
      if (debug) then
        write (iw, *) ' GSE:', gse
        write (iw, *) ' EIGB:', (eigb(i),i=1,10)
        write (iw, *) ' WORK:', (work(i,1),i=1,10)
      end if
      call mecih (work, wmat, nmos, lab, xy)
!
!     NON-RELAXED C.I CONTRIBUTION TO THE ENERGY DERIVATIVE.
      sum = 0.D0
      do i = 1, nstate
        j = (i - 1)*lab + 1
        call supdot (work, wmat, vectci(j), lab)
        sum = sum + ddot(lab,vectci(j),1,work,1)
      end do
      grad = (grad + sum/nstate)*const
      if (debug) then
        write (iw, '('' * * * GRADIENT COMPONENT NUMBER'',I4)') number
        write (iw, &
      '('' NON-RELAXED C.I-ACTIVE FOCK EIGENVALUES DERIVATIVES (E.V.)'')')
        write (iw, '(8F10.4)') (fd(lcut-1+i),i=1,nmos)
        write (iw, '('' NON-RELAXED 2-ELECTRONS DERIVATIVES (E.V.)''/  &
        &''    I    J    K    L       d<I(1)J(1)|K(2)L(2)>'')')
        do i = 1, nmos
          do j = 1, i
            do k = 1, i
              ll = k
              if (k == i) ll = j
              do l = 1, ll
                write (iw, '(4I5,F20.10)') nelec + i, nelec + j, nelec + k, nelec + l, xy(i,j,k,l)
              end do
            end do
          end do
        end do
        write (iw, '('' NON-RELAXED GRADIENT COMPONENT'',F10.4,'' KCAL/MOLE'')') grad
        call timer ('AFTER DERI1')
      end if
      return
      end subroutine deri1
