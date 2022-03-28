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

        subroutine nllsq()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use chanel_C, only : iw, ilog, igpt, log
      use molkst_C, only : iflepo, keywrd, last, moperr, escf, nvar, &
      & tleft, nscf, tdump, time0, line
      use common_arrays_C, only : xparam, grad
      implicit none
      double precision, dimension (:,:), allocatable  :: q, r
      integer :: icyc, irst, jrst, maxcyc, m, nrst, ifrtl, nsst, ixso, &
        np1, np2, jcyc, i, j, jj, ii, i1, j1, k, jend, ncount
      integer, dimension(6) :: iiium
      double precision, dimension(nvar) :: y, efs, p, xlast
      double precision :: alf, ssq, pn, tol2, tols1, tols2, tols5, tols6, tlast,  &
        eps, t, time1, temp, prt, efsss, sin, cos, work, pnlast, pn2, ttmp, &
        ssqlst, time2, tcycle, yn, const, ytail, bet, tmp, tolx, tprt
      logical :: middle, resfil
      character :: txt
      double precision, external :: ddot, reada, seconds
!***********************************************************************
!
!  NLLSQ IS A NON-DERIVATIVE, NONLINEAR LEAST-SQUARES MINIMIZER. IT USES
!        BARTEL'S PROCEDURE TO MINIMIZE A FUNCTION WHICH IS A SUM OF
!        SQUARES.
!
!    ON INPUT nvar    = NUMBER OF UNKNOWNS
!             xparam    = PARAMETERS OF FUNCTION TO BE MINIMIZED.
!
!    ON EXIT  xparam    = OPTIMIZED PARAMETERS.
!
!    THE FUNCTION TO BE MINIMIZED IS "COMPFG". COMPFG MUST HAVE THE
!    CALLING SEQUENCE
!                  CALL COMPFG(XPARAM,.TRUE.,ESCF,.TRUE.,EFS,.TRUE.)
!                  SSQ=DOT(EFS,EFS,nvar)
!    WHERE   EFS  IS A VECTOR WHICH  COMPFG  FILLS WITH THE nvar INDIVIDUAL
!                 COMPONENTS OF THE ERROR FUNCTION AT THE POINT xparam
!            SSQ IS THE VALUE OF THE SUM OF THE  EFS  SQUARED.
!    IN THIS FORMULATION OF NLLSQ M AND nvar ARE THE SAME.
!    THE PRECISE DEFINITIONS OF THESE TWO QUANTITIES IS:
!
!     nvar = NUMBER OF PARAMETERS TO BE OPTIMIZED.
!     M = NUMBER OF REFERENCE FUNCTIONS. M MUST BE GREATER THEN, OR
!         EQUAL TO, nvar
!***********************************************************************
!     Q = ORTHOGONAL MATRIX   (M BY M)
!     R = RIGHT-TRIANGULAR MATRIX   (M BY nvar)
!     MXCNT(1) = MAX ALLOW OVERALL FUN EVALS
!     MXCNT(2) = MAX ALLOW NO OF FNC EVALS PER LIN SEARCH
!     TOLS1 = RELATIVE TOLERANCE ON xparam OVERALL
!     TOLS2 = ABSOLUTE TOLERANCE ON xparam OVERALL
!     TOLS5 = RELATIVE TOLERANCE ON xparam FOR LINEAR SEARCHES
!     TOLS6 = ABSOLUTE TOLERANCE ON xparam FOR LINEAR SEARCHES
!     NRST = NUMBER OF CYCLES BETWEEN SIDESTEPS
!     **********
      allocate(q(nvar,nvar), r(nvar + 2, nvar + 2), stat = i)
      if (i /= 0) then
        write(iw,*)" Unable to allocate memory for Hessian matrices in NLLSQ"
        call mopend("Unable to allocate memory for Hessian matrices in NLLSQ")
        goto 99
      end if
      middle = index(keywrd,' RESTART') /= 0
      maxcyc = 100000
      if (index(keywrd,' CYCLES') /= 0) maxcyc = nint(reada(keywrd,index(keywrd,' CYCLES')))
      iflepo = 10
      pn = 0.d0
      ssq = 0.d0
      alf = 0.d0
!*
      m = nvar
!*
      tol2 = 4.D-1
      if (index(keywrd,' GNORM') /= 0) then
        tol2 = reada(keywrd,index(keywrd,' GNORM'))
        if (tol2<0.01D0 .and. index(keywrd,' LET')==0) then
          write (iw, '(/,A)') '  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01'
          tol2 = 0.01D0
        end if
      end if
      last = 0
      tols1 = 1.D-12
      tols2 = 1.D-10
      tols5 = 1.D-6
      tols6 = 1.D-3
      nrst = 4
      tlast = tleft
      resfil = .FALSE.
      tleft = tleft - seconds(2) + time0
!     **********
!     SET UP COUNTERS AND SWITCHES
!     **********
      ifrtl = 0
      nsst = 0
      ixso = nvar
      np1 = nvar + 1
      np2 = nvar + 2
      icyc = 0
      jcyc = icyc
      irst = 0
      jrst = 1
      eps = tols5
      t = tols6
      grad = 0.d0
!     **********
!     GET STARTING-POINT FUNCTION VALUE
!     SET UP ESTIMATE OF INITIAL LINE STEP
!     **********
      if (middle) then
        q = 0.d0
        r = 0.d0
        call parsav (0, nvar, m, q, r, grad, xlast, iiium)
        if (moperr) goto 99
        nscf = iiium(1)
        close(igpt)
        ncount = iiium(5)
        xparam(:nvar) = xlast(:nvar)
        time1 = seconds(2)
        if (index(keywrd,' 1SCF') /= 0) then
          iflepo = 13
          last = 1
          goto 99
        end if
        jcyc = icyc
        go to 60
      end if
      call compfg (xparam, .TRUE., escf, .TRUE., grad, .TRUE.)
      ssq = ddot(nvar,grad,1,grad,1)
      ncount = 1
      do i = 1, m
        do j = 1, nvar
          r(i,j) = 0.0D0
          if (i /= j) cycle
          r(i,j) = 1.0D0
        end do
        do j = i, m
          q(i,j) = 0.0D0
          q(j,i) = 0.0D0
          if (i /= j) cycle
          q(i,i) = 1.0D0
        end do
      end do
      temp = 0.0D0
      temp = sum(xparam(:nvar)**2)
      alf = 100.0D0*(eps*sqrt(temp) + t)
!     **********
!     MAIN LOOP
!     **********
      time1 = seconds(2)
   60 continue
      ifrtl = ifrtl + 1
      icyc = icyc + 1
      irst = irst + 1
!     **********
!     SET  PRT,  THE LEVENBERG-MARQUARDT PARAMETER.
!     **********
      prt = sqrt(ssq)
!     **********
!     IF A SIDESTEP IS TO BE TAKEN, GO TO 31
!     **********
      if (irst < nrst) then
!     **********
!     SOLVE THE SYSTEM    Q*R*P = -grad    IN THE LEAST-SQUARES SENSE
!     **********
        nsst = 0
        do i = 1, m
          temp = 0.0D0
          temp = -sum(q(:m,i)*grad(:m))
          efs(i) = temp
        end do
        do j = 1, nvar
          jj = np1 - j
          do i = 1, j
            ii = np2 - i
            r(ii,jj) = r(i,j)
          end do
        end do
        l160: do i = 1, nvar
          i1 = i + 1
          y(i) = prt
          efsss = 0.0D0
          if (i < nvar) then
            y(i1:nvar) = 0.0D0
          end if
          do j = i, nvar
            ii = np2 - j
            jj = np1 - j
            if (abs(y(j)) >= abs(r(ii,jj))) then
              temp = y(j)*sqrt(1.0D0 + (r(ii,jj)/y(j))**2)
            else
              temp = r(ii,jj)*sqrt(1.0D0 + (y(j)/r(ii,jj))**2)
            end if
            sin = r(ii,jj)/temp
            cos = y(j)/temp
            r(ii,jj) = temp
            temp = efs(j)
            efs(j) = sin*temp + cos*efsss
            efsss = sin*efsss - cos*temp
            if (j >= nvar) cycle  l160
            j1 = j + 1
            do k = j1, nvar
              jj = np1 - k
              temp = r(ii,jj)
              r(ii,jj) = sin*temp + cos*y(k)
              y(k) = sin*y(k) - cos*temp
            end do
          end do
        end do l160
        p(nvar) = efs(nvar)/r(2,1)
        i = nvar
  170   continue
        i = i - 1
        if (i <= 0) go to 200
        temp = efs(i)
        k = i + 1
        ii = np2 - i
        temp = temp - dot_product(r(ii,np1-k:np1-nvar:(-1)),p(k:nvar))
        jj = np1 - i
        p(i) = temp/r(ii,jj)
        go to 170
  200   continue
      else
!     **********
!     SIDESTEP SECTION
!     **********
        jrst = jrst + 1
        nsst = nsst + 1
        if (nsst >= ixso) go to 670
        if (jrst > nvar) jrst = 2
        irst = 0
!     **********
!     PRODUCTION OF A VECTOR ORTHOGONAL TO THE LAST P-VECTOR
!     **********
        work = pn*(abs(p(1))+pn)
        temp = p(jrst)
        p(1) = temp*(p(1)+sign(pn,p(1)))
        p(2:nvar) = temp*p(2:nvar)
        p(jrst) = p(jrst) - work
      end if
!     **********
!     COMPUTE NORM AND NORM-SQUARE OF THE P-VECTOR
!     **********
      pnlast = pn
      pn = 0.D0
      pn2 = 0.0D0
      do i = 1, nvar
        pn = pn + abs(p(i))
        pn2 = pn2 + p(i)**2
      end do
      if (pn < 1.D-20) then
        write (iw, &
      '('' SYSTEM DOES NOT APPEAR TO BE OPTIMIZABLE.'',/    ,'' THIS CAN HAPPEN&
      & IF (A) IT WAS OPTIMIZED TO BEGIN WITH'',/     ,'' OR                 (B)&
      & IT IS NEITHER A GROUND NOR A'',        '' TRANSITION STATE'')')
        call geout (1)
        call mopend ('SYSTEM DOES NOT APPEAR TO BE OPTIMIZABLE.')
        goto 99
      end if
      pn2 = dmax1(1.D-20,pn2)
      pn = sqrt(pn2)
      alf = min(1.D20,alf)
      if (icyc > 1) then
        alf = alf*1.D-20*pnlast/pn
        alf = min(1.D10,alf)
        alf = alf*1.D20
      end if
      ttmp = alf*pn
      if (ttmp < 0.0001D0) alf = 0.001D0/pn
!     **********
!     PRINTING SECTION
!     **********
!#      WRITE(IW,501)TLEFT,ICYC,SSQ
      efs(:nvar) = xparam(:nvar)
!     **********
!     PERFORM LINE-MINIMIZATION FROM POINT xparam IN DIRECTION P OR -P
!     **********
      ssqlst = ssq
      efs(:nvar) = 0.D0
      xlast(:nvar) = xparam(:nvar)
      call locmin (m, xparam, nvar, p, ssq, alf, efs, ncount)
      if (ssqlst < ssq) then
     !   if (ierr == 0) ssq = ssqlst
        xparam(:nvar) = xlast(:nvar)
        irst = nrst
        pn = pnlast
        go to 590
      end if
!     **********
!     PRODUCE THE VECTOR   R*P
!     **********
      do i = 1, nvar
        temp = 0.0D0
        temp = sum(r(i,i:nvar)*p(i:nvar))
        y(i) = temp
      end do
!     **********
!     PRODUCE THE VECTOR ...
!                  Y  =    (EFS-grad-ALF*Q*R*P)/(ALF*(NORMSQUARE(P))
!     COMPUTE NORM OF THIS VECTOR AS WELL
!     **********
      work = alf*pn2
      yn = 0.0D0
      do i = 1, m
        temp = sum(q(i,:nvar)*y(:nvar))
        temp = efs(i) - grad(i) - alf*temp
        grad(i) = efs(i)
        yn = yn + temp**2
        efs(i) = temp/work
      end do
      yn = sqrt(yn)/work
!     **********
!     THE BROYDEN UPDATE   NEW MATRIX = OLD MATRIX + Y*(P-TRANS)
!     HAS BEEN FORMED.  IT IS NOW NECESSARY TO UPDATE THE  QR DECOMP.
!     FIRST LET    Y = (Q-TRANS)*Y.
!     **********
      do i = 1, m
        temp = 0.0D0
        temp = sum(q(:m,i)*efs(:m))
        y(i) = temp
      end do
!     **********
!     REDUCE THE VECTOR Y TO A MULTIPLE OF THE FIRST UNIT VECTOR USING
!     A HOUSEHOLDER TRANSFORMATION FOR COMPONENTS nvar+1 THROUGH M AND
!     ELEMENTARY ROTATIONS FOR THE FIRST nvar+1 COMPONENTS.  APPLY ALL
!     TRANSFORMATIONS TRANSPOSED ON THE RIGHT TO THE MATRIX Q, AND
!     APPLY THE ROTATIONS ON THE LEFT TO THE MATRIX R.
!     THIS GIVES    (Q*(V-TRANS))*((V*R) + (V*Y)*(P-TRANS)),    WHERE
!     V IS THE COMPOSITE OF THE TRANSFORMATIONS.  THE MATRIX
!     ((V*R) + (V*Y)*(P-TRANS))    IS UPPER HESSENBERG.
!     **********
      if (m > np1) then
!
! THE NEXT THREE LINES WERE INSERTED TO TRY TO GET ROUND OVERFLOW BUGS.
!
        const = 1.D-12
        do i = np1, m
          const = max(abs(y(np1)),const)
        end do
        ytail = sum((y(np1:m)/const)**2)
        ytail = sqrt(ytail)*const
        bet = (1.0D25/ytail)/(ytail + abs(y(np1)))
        y(np1) = sign(ytail + abs(y(np1)),y(np1))
        do i = 1, m
          tmp = sum(q(i,np1:m)*y(np1:m)*1.D-25)
          tmp = bet*tmp
          q(i,np1:m) = q(i,np1:m) - tmp*y(np1:m)
        end do
        y(np1) = ytail
        i = np1
      else
        i = m
      end if
  410 continue
      j = i
      i = i - 1
      if (i <= 0) go to 480
      do while(y(j) == 0.D0)
        j = i
        i = i - 1
        if (i <= 0) go to 480
      end do
      if (abs(y(i)) >= abs(y(j))) then
        temp = abs(y(i))*sqrt(1.0D0 + (y(j)/y(i))**2)
      else
        temp = abs(y(j))*sqrt(1.0D0 + (y(i)/y(j))**2)
      end if
      cos = y(i)/temp
      sin = y(j)/temp
      y(i) = temp
      do k = 1, m
        temp = cos*q(k,i) + sin*q(k,j)
        work = (-sin*q(k,i)) + cos*q(k,j)
        q(k,i) = temp
        q(k,j) = work
      end do
      if (i > nvar) go to 410
      r(j,i) = -sin*r(i,i)
      r(i,i) = cos*r(i,i)
      if (j > nvar) go to 410
      do k = j, nvar
        temp = cos*r(i,k) + sin*r(j,k)
        work = (-sin*r(i,k)) + cos*r(j,k)
        r(i,k) = temp
        r(j,k) = work
      end do
      go to 410
  480 continue
      r(1,:nvar) = r(1,:nvar) + yn*p(:nvar)
      jend = np1
      if (m == nvar) jend = nvar
      do j = 2, jend
        i = j - 1
        if (r(j,i) == 0.D0) cycle
        if (abs(r(i,i)) >= abs(r(j,i))) then
          temp = abs(r(i,i))*sqrt(1.0D0 + (r(j,i)/r(i,i))**2)
        else
          temp = abs(r(j,i))*sqrt(1.0D0 + (r(i,i)/r(j,i))**2)
        end if
        cos = r(i,i)/temp
        sin = r(j,i)/temp
        r(i,i) = temp
        if (j <= nvar) then
          do k = j, nvar
            temp = cos*r(i,k) + sin*r(j,k)
            work = (-sin*r(i,k)) + cos*r(j,k)
            r(i,k) = temp
            r(j,k) = work
          end do
        end if
        do k = 1, m
          temp = cos*q(k,i) + sin*q(k,j)
          work = (-sin*q(k,i)) + cos*q(k,j)
          q(k,i) = temp
          q(k,j) = work
        end do
      end do
!     **********
!     CHECK THE STOPPING CRITERIA
!     **********
      temp = 0.0D0
      temp = sum(xparam(:nvar)**2)
      tolx = tols1*sqrt(temp) + tols2
      if (sqrt(alf*pn2) <= tolx) go to 650
      if (ssq >= 2.D0*nvar) go to 590
      do i = 1, nvar
!*****
!     The stopping criterion is that no individual gradient be
!         greater than TOL2
!*****
        if (abs(grad(i)) >= tol2) go to 590
      end do
!#      WRITE(IW,730) SSQ
      go to 660
  590 continue
      time2 = time1
      time1 = seconds(2)
      tcycle = time1 - time2
      tleft = tleft - tcycle
      if (tleft < 0.0D0) tleft = -0.1D0
      call prttim (tleft, tprt, txt)
      if (resfil) then
        write (iw, 600) tprt, txt, sqrt(ssq), escf
  600   format('  RESTART FILE WRITTEN,     TIME LEFT:',f6.2,a1,'  GRAD.:',f10.3,&
          ' HEAT:',g14.7)
        resfil = .FALSE.
      else
        write (line, '(" CYCLE:", i6, " TIME:", f8.3, " TIME LEFT:", &
                 & f6.2, a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
                 icyc, min(tcycle,9999.99D0), tprt, txt, &
                 & min(sqrt(ssq),999999.999D0), escf
        write(iw,"(a)")line(:len_trim(line))
        endfile (iw)
        backspace (iw)
        if (log) write (ilog, "(a)")line(:len_trim(line))
        call to_screen(line)
        endfile (iw)
        backspace (iw)
        if (log) then
          endfile (ilog)
          backspace (ilog)
        end if
      end if
      endfile (iw)
      backspace (iw)
      if (log) then
        endfile (ilog)
        backspace (ilog)
      end if
      if (tlast - tleft > tdump) then
        tlast = tleft
        resfil = .TRUE.
        xlast(:nvar) = xparam(:nvar)
        iiium(1) = nscf
        call parsav (2, nvar, m, q, r, grad, xlast, iiium)
        if (moperr) goto 99
      end if
      if (tleft>tcycle*2 .and. icyc-jcyc<maxcyc) go to 60
      iiium(5) = ncount
      xlast(:nvar) = xparam(:nvar)
      iiium(1) = nscf
      call parsav (1, nvar, m, q, r, grad, xlast, iiium)
      if (moperr) goto 99
      iflepo = -1
      goto 99
  650 continue
      write (iw, 760) ncount
      go to 870
  660 continue
      write (iw, 770) ncount
      go to 870
  670 continue
      write (iw, 680) ixso
  680 format(' ',5x,'ATTEMPT TO GO DOWNHILL IS UNSUCCESSFUL AFTER',i5,5x,&
        'ORTHOGONAL SEARCHES')
      go to 870
  760 format(' TEST ON xparam SATISFIED, NUMBER OF FUNCTION CALLS = ',i5)
  770 format(' TEST ON SSQ SATISFIED, NUMBER OF FUNCTION CALLS = ',i5)
  870 continue
      last = 1
      call compfg (xparam, .TRUE., escf, .TRUE., grad, .TRUE.)
  99  if (allocated(q)) deallocate(q, stat = i)
      if (i /= 0) return
      if (allocated(r)) deallocate(r, stat = i)
      if (i /= 0) return
      return
      end subroutine nllsq
