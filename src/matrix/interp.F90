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

     subroutine interp(np, nq, mode, e, fp, cp, theta, vec_interp, fock_interp, &
     p_interp, h_interp, vecl, eold_ref)
      use molkst_C, only : keywrd, numcal, norbs
      use common_arrays_C, only :
      use chanel_C, only : iw
      use funcon_C, only : fpc_9, pi
      implicit none
      integer  :: np
      integer  :: nq
      integer , intent(inout) :: mode
      double precision , intent(in) :: e
      double precision , intent(in) :: fp((norbs*(norbs + 1))/2)
      double precision , intent(inout) :: cp(norbs,norbs)
      double precision , dimension (norbs, norbs), intent (inout) :: vec_interp, &
    & fock_interp, p_interp
      double precision , dimension (norbs*norbs), intent (inout) :: h_interp
      double precision, dimension(norbs), intent (inout) :: theta
      double precision, intent (inout)  :: vecl(norbs*norbs)
      double precision, intent(inout) :: eold_ref
!
      integer ::  icalcn, minpq, np1, np2, ii, i, i1, j, k, ik, ipoint, j1, &
        ij, il, k1, k2
      double precision :: xold, zero, ff, dum, dedx, deold, xnow, denow, enow, ck, sk
      integer :: npnts
      double precision, dimension(12) :: x, f, df
      double precision :: xlow, xhigh, xmin, emin
      logical :: debug
      save
!-----------------------------------------------
!*********************************************************************
!
! INTERP: AN INTERPOLATION PROCEDURE FOR FORCING SCF CONVERGANCE
!         ORIGINAL THEORY AND FORTRAN WRITTEN BY R.N. CAMP AND
!         H.F. KING, J. CHEM. PHYS. 75, 268 (1981)
!*********************************************************************
!
! ON INPUT norbs     = NUMBER OF ORBITALS
!          NP    = NUMBER OF FILLED LEVELS
!          NQ    = NUMBER OF EMPTY LEVELS
!          MODE  = 1, DO NOT RESET.
!          E     = ENERGY
!          FP    = FOCK MATRIX, AS LOWER HALF TRIANGLE, PACKED
!          CP    = EIGENVECTORS OF FOCK MATRIX OF ITERATION -1
!                  AS PACKED ARRAY OF NORBS**2 COEFFICIENTS
!
! ON OUTPUT CP   = BEST GUESSED SET OF EIGENVECTORS
!           MODE = 2 OR 3 - USED BY CALLING PROGRAM
!*********************************************************************
      data icalcn/ 0/
      data zero, ff/ 0.0D0, 0.9D0/
      data xold/ 0.D0/
      if (icalcn /= numcal) then
        debug = index(keywrd,'INTERP') /= 0
        icalcn = numcal
        emin = 0.d0
        npnts = 0
      end if
!
!         FF=FACTOR FOR CONVERGENCE TEST FOR 1D SEARCH.
!
      minpq = min0(np,nq)
      np1 = np + 1
      np2 = max0(1,np/2)
      if (mode /= 2) then
!
!     (MODE=1 OR 3 ENTRY)
!     TRANSFORM FOCK MATRIX TO CURRENT MO BASIS.
!     ONLY THE OFF DIAGONAL OCC-VIRT BLOCK IS COMPUTED.
!     STORE IN FOCK ARRAY
!
        ii = 0
        do i = 1, norbs
          i1 = i + 1
          do j = 1, nq
            dum = zero
            dum = dum + sum(fp(ii+1:i+ii)*cp(:i,j+np))
            if (i /= norbs) then
              ik = ii + i + i
              do k = i1, norbs
                dum = dum + fp(ik)*cp(k,j+np)
                ik = ik + k
              end do
            end if
            p_interp(i,j) = dum
          end do
          ii = ii + i
        end do
        do i = 1, np
          do j = 1, nq
            dum = zero
            dum = dum + sum(cp(:norbs,i)*p_interp(:norbs,j))
            fock_interp(i,j) = dum
          end do
        end do
        if (mode /= 3) then
!
!     CURRENT POINT BECOMES OLD POINT (MODE=1 ENTRY)
!
          vec_interp(:norbs,:norbs) = cp(:norbs,:norbs)
          eold_ref = e
          xold = 1.0D0
          mode = 2
          return
        end if
!
!     (MODE=3 ENTRY)
!     F CORRESPONDS TO CURRENT POINT IN CORRESPONDING REPRESENTATION.
!     VEC DOES NOT HOLD CURRENT VECTORS. VEC SET IN LAST MODE=2 ENTRY.
!
        npnts = npnts + 1
        if (debug) write (iw, '(''   INTERPOLATED ENERGY:'',F13.6)') e*fpc_9
        ipoint = npnts
      else
!
!    (MODE=2 ENTRY) CALCULATE THETA, AND U, V, W MATRICES.
!                   U ROTATES CURRENT INTO OLD MO.
!                   V ROTATES CURRENT INTO CORRESPONDING CURRENT MO.
!                   W ROTATES OLD INTO CORRESPONDING OLD MO.
!
        j1 = 1
        do i = 1, norbs
          if (i == np1) j1 = np1
          do j = j1, norbs
            p_interp(i,j) = zero
            p_interp(i,j) = p_interp(i,j) + sum(cp(:norbs,i)*vec_interp(:norbs,j))
          end do
        end do
!
!     U = CP(DAGGER)*VEC IS NOW IN P ARRAY.
!     VEC IS NOW AVAILABLE FOR TEMPORARY STORAGE.
!
        ij = 0
        do i = 1, np
          do j = 1, i
            ij = ij + 1
            h_interp(ij) = 0.D0
            h_interp(ij) = h_interp(ij) + sum(p_interp(i,np1:norbs)*p_interp(j,np1:norbs))
          end do
        end do
        call rsp (h_interp, np, theta, vecl)
        do i = np, 1, -1
          il = i*np - np
          vec_interp(np:1:(-1),i) = vecl(np+il:il+1:(-1))
        end do
        do i = 1, np2
          dum = theta(np1-i)
          theta(np1-i) = theta(i)
          theta(i) = dum
          do j = 1, np
            dum = vec_interp(j,np1-i)
            vec_interp(j,np1-i) = vec_interp(j,i)
            vec_interp(j,i) = dum
          end do
        end do
        do i = 1, minpq
          theta(i) = max(theta(i),zero)
          theta(i) = min(theta(i),1.D0)
          theta(i) = asin(sqrt(theta(i)))
        end do
!
!     THETA MATRIX HAS NOW BEEN CALCULATED, ALSO UNITARY VP MATRIX
!     HAS BEEN CALCULATED AND STORED IN FIRST NP COLUMNS OF VEC MATRIX.
!     NOW COMPUTE WQ
!
        do i = 1, nq
          do j = 1, minpq
            vec_interp(i,np+j) = zero
            vec_interp(i,np+j) = vec_interp(i,np+j) + sum(p_interp(:np,np+i)*vec_interp(:np,j))
          end do
        end do
        call schmit (vec_interp(1,np1), nq, norbs)
!
!     UNITARY WQ MATRIX NOW IN LAST NQ COLUMNS OF VEC MATRIX.
!     TRANSPOSE NP BY NP BLOCK OF U STORED IN P
!
        do i = 1, np
          do j = 1, i
            dum = p_interp(i,j)
            p_interp(i,j) = p_interp(j,i)
            p_interp(j,i) = dum
          end do
        end do
!
!     CALCULATE WP MATRIX AND HOLD IN FIRST NP COLUMNS OF P
!
        do i = 1, np
          h_interp(:np) = p_interp(i,:np)
          do j = 1, np
            p_interp(i,j) = zero
            p_interp(i,j) = p_interp(i,j) + sum(h_interp(:np)*vec_interp(:np,j))
          end do
        end do
        call schmib (p_interp, np, norbs)
!
!     CALCULATE VQ MATRIX AND HOLD IN LAST NQ COLUMNS OF P MATRIX.
!
        do i = 1, nq
          h_interp(:nq) = p_interp(np+i,np+1:nq+np)
          do j = np1, norbs
            p_interp(i,j) = zero
            p_interp(i,j) = p_interp(i,j) + sum(h_interp(:nq)*vec_interp(:nq,j))
          end do
        end do
        call schmib (p_interp(1,np1), nq, norbs)
!
!     CALCULATE (DE/DX) AT OLD POINT
!
        dedx = zero
        do i = 1, np
          do j = 1, nq
            dum = zero
            dum = dum + sum(theta(:minpq)*p_interp(i,:minpq)*vec_interp(j,np+1:minpq+np))
            dedx = dedx + dum*fock_interp(i,j)
          end do
        end do
!
!     STORE OLD POINT INFORMATION FOR SPLINE FIT
!
        deold = -4.0D0*dedx
        x(2) = xold
        f(2) = eold_ref
        df(2) = deold
!
!     MOVE VP OUT OF vec_interp ARRAY INTO FIRST NP COLUMNS OF p_interp MATRIX.
!
        p_interp(:np,:np) = vec_interp(:np,:np)
        k1 = 0
        k2 = np
        do j = 1, norbs
          if (j == np1) then
            k1 = np
            k2 = nq
          end if
          do i = 1, norbs
            dum = zero
            dum = dum + sum(cp(i,k1+1:k2+k1)*p_interp(:k2,j))
            vec_interp(i,j) = dum
          end do
        end do
!
!     CORRESPONDING CURRENT MO VECTORS NOW HELD IN vec_interp.
!     COMPUTE vec_interp(DAGGER)*FP*vec_interp
!     STORE OFF-DIAGONAL BLOCK IN fock_interp ARRAY.
!
        ii = 0
        do i = 1, norbs
          i1 = i + 1
          do j = 1, nq
            dum = zero
            dum = dum + sum(fp(ii+1:i+ii)*vec_interp(:i,j+np))
            if (i /= norbs) then
              ik = ii + i + i
              do k = i1, norbs
                dum = dum + fp(ik)*vec_interp(k,j+np)
                ik = ik + k
              end do
            end if
            p_interp(i,j) = dum
          end do
          ii = ii + i
        end do
        do i = 1, np
          do j = 1, nq
            dum = zero
            dum = dum + sum(vec_interp(:norbs,i)*p_interp(:norbs,j))
            fock_interp(i,j) = dum
          end do
        end do
!
!     SET LIMITS ON RANGE OF 1-D SEARCH
!
        npnts = 2
        ipoint = 1
        xnow = zero
        xhigh = pi/(2.d0*theta(1))
!
!     1.5708 IS MAXIMUM ROTATION ANGLE (90 DEGREE = 3.14159/2 RADIAN).
!
        xlow = -0.5D0*xhigh
      end if
!
!     CALCULATE (DE/DX) AT CURRENT POINT AND
!     STORE INFORMATION FOR SPLINE FIT
!     ***** JUMP POINT FOR MODE=3 ENTRY *****
!
      dedx = zero
      do k = 1, minpq
        dedx = dedx + theta(k)*fock_interp(k,k)
      end do
      denow = -4.0D0*dedx
      enow = e
!     PERFORM 1-D SEARCH AND DETERMINE EXIT MODE.
!
      x(ipoint) = xnow
      f(ipoint) = enow
      df(ipoint) = denow
      call spline(x, f, df, xhigh, xlow,  xmin, npnts)
      if (eold_ref - enow<=ff*(eold_ref - emin) .and. ipoint<=10) then
!
!     (MODE=3 EXIT) RECOMPUTE CP VECTORS AT PREDICTED MINIMUM.
!
        xnow = xmin
        cp = vec_interp
        do k = 1, minpq
          ck = cos(xnow*theta(k))
          sk = sin(xnow*theta(k))
          if (debug) write (iw, '('' ROTATION ANGLE:'',F12.4)') sk*57.29578D0
          cp(:norbs,k) = ck*vec_interp(:norbs,k) - sk*vec_interp(:norbs,np+k)
          cp(:norbs,np+k) = sk*vec_interp(:norbs,k) + ck*vec_interp(:norbs,np+k)
        end do
        mode = 3
        return
      end if
!
!     (MODE=2 EXIT) CURRENT VECTORS GIVE SATISFACTORY ENERGY IMPROVEMENT
!     CURRENT POINT BECOMES OLD POINT FOR THE NEXT 1-D SEARCH.
!
      if (mode /= 2) then
        vec_interp(:norbs,:norbs) = cp(:norbs,:norbs)
        mode = 2
      end if
      eold_ref = enow
      if (npnts <= 200) return
      write (iw, 600)
      do k = 1, npnts
        write (iw, 610) k, x(k), f(k), df(k)
      end do
      write (iw, 620)
      return
  600 format('  K','     X(K) ','       F(K)    ','     DF(K)')
  610 format(i3,f10.5,2f15.10)
  620 format(10x)
      end subroutine interp
      subroutine spline(x, f, df, xhigh, xlow,  xmin, npnts)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!***********************************************************************
      implicit none
      integer :: npnts
      double precision, dimension(12) :: x, f, df
      double precision :: xlow, xhigh, xmin, fmin, dfmin
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k1, j2, n1, k
      double precision :: close, big, huge, ustep, dstep, step, xstart, xstop, dx, &
        x1, x2, dum, a, b, c, bb, ac3, xk, r, fm
      logical :: skip1, skip2

      save close, big, huge, ustep, dstep
!-----------------------------------------------
!
!     FIT F(X) BY A CUBIC SPLINE GIVEN VALUES OF THE FUNCTION
!     AND ITS FIRST DERIVATIVE AT npnts PNTS.
!     SUBROUTINE RETURNS VALUES OF XMIN,FMIN, AND DFMIN
!     AND MAY REORDER THE DATA.
!     XLOW AND XHIGH SET LIMITS ON THE INTERVAL WITHIN WHICH
!     TO SEARCH.  SUBROUTINE MAY FURTHER REDUCE THIS INTERVAL.
!
      data close, big, huge, ustep, dstep/ 1.0D-8, 500.0D0, 1.0D+10, 1.0D0, &
        2.0D0/
!
!     SUBROUTINE ASSUMES THAT THE FIRST npnts-1 DATA PNTS HAVE BEEN
!     PREVIOUSLY ORDERED,  X(I).LT.X(I+1) FOR I=1,2,...,npnts-2
!     NOW MOVE NTH POINT TO ITS PROPER PLACE.
!
      xmin = x(npnts)
      fmin = f(npnts)
      dfmin = df(npnts)
      n1 = npnts - 1
      k = n1
      k1 = k
      j2 = min0(1,k1)
      do k = k1, j2, -1
        if (x(k) < xmin) go to 20
        x(k+1) = x(k)
        f(k+1) = f(k)
        df(k+1) = df(k)
      end do
   20 continue
      x(k+1) = xmin
      f(k+1) = fmin
      df(k+1) = dfmin
      if (df(1) > 0.0D0) then
!
!     DEFINE THE INTERVAL WITHIN WHICH WE TRUST THE SPLINE FIT.
!     USTEP =  UP HILL STEP SIZE FACTOR
!     DSTEP = DOWN HILL STEP SIZE FACTOR
!
        step = dstep
      else
        step = ustep
      end if
      xstart = x(1) - step*(x(2)-x(1))
      xstart = max(xstart,xlow)
      if (df(npnts) > 0.0D0) then
        step = ustep
      else
        step = dstep
      end if
      xstop = x(npnts) + step*(x(npnts)-x(n1))
      xstop = min(xstop,xhigh)
!
!     SEARCH FOR MINIMUM
!
      do k = 1, n1
        skip1 = k /= 1
        skip2 = k /= n1
        if (f(k) < fmin) then
          xmin = x(k)
          fmin = f(k)
          dfmin = df(k)
        end if
        dx = x(k+1) - x(k)
!
!     SKIP INTERVAL IF PNTS ARE TOO CLOSE TOGETHER
!
        if (dx <= close) cycle
        x1 = 0.0D0
        if (k == 1) x1 = xstart - x(1)
        x2 = dx
        if (k == n1) x2 = xstop - x(n1)
!
!     (A,B,C)=COEF OF (CUBIC,QUADRATIC,LINEAR) TERMS
!
        dum = (f(k+1)-f(k))/dx
        a = (df(k)+df(k+1)-dum-dum)/(dx*dx)
        b = (dum + dum + dum - df(k)-df(k)-df(k+1))/dx
        c = df(k)
!
!     XK = X-X(K) AT THE MINIMUM WITHIN THE KTH SUBINTERVAL
!     TEST FOR PATHOLOGICAL CASES.
!
        bb = b*b
        ac3 = (a + a + a)*c
        if (bb < ac3) go to 90
        if (b <= 0.0D0) then
          if (abs(b) > huge*abs(a)) go to 90
        else
          if (bb > big*abs(ac3)) go to 60
        end if
!
!     WELL BEHAVED CUBIC
!
        xk = ((-b) + sqrt(bb - ac3))/(a + a + a)
        go to 70
!
!     CUBIC IS DOMINATED BY QUADRATIC TERM
!
   60   continue
        r = ac3/bb
        xk = -(((0.039063D0*r + 0.0625D0)*r + 0.125D0)*r + 0.5D0)*c/b
   70   continue
        if (xk<x1 .or. xk>x2) go to 90
   80   continue
        fm = ((a*xk + b)*xk + c)*xk + f(k)
        if (fm <= fmin) then
          xmin = xk + x(k)
          fmin = fm
          dfmin = ((a + a + a)*xk + b + b)*xk + c
!
!     EXTRAPOLATE TO END OF INTERVAL IF K=1 AND/OR K=N1
!
        end if
   90   continue
        if (skip1) go to 100
        skip1 = .TRUE.
        xk = x1
        go to 80
  100   continue
        if (skip2) cycle
        skip2 = .TRUE.
        xk = x2
        go to 80
      end do
      return
      end subroutine spline

