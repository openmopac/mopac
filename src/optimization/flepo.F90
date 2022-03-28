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

      subroutine flepo(xparam, nvar, funct1)
      use molkst_C, only : numcal, gnorm, iflepo, keywrd, emin, tleft, &
      & time0, moperr, nscf, limscf, tdump, last, cosine, line
      use common_arrays_C, only : hesinv, grad
      USE chanel_C, only : iw, ilog, log, input_fn
      implicit none
      integer  :: nvar
      double precision  :: funct1
      double precision  :: xparam(nvar)
!
      integer , dimension(9) :: mdfp
      integer :: icalcn, igg1, nrst, ihdim, itry1, jcyc, lnstop, irepet, jnrst&
        , ncount, maxcyc, i, ireset, icyc, i80, ii, j, k, nto6, &
        iinc1, iinc2, ic
      double precision, dimension(nvar) :: xvar, gvar, xd, gd, glast, xlast, gg, &
        pvect
      double precision, dimension(9) :: xdfp
      double precision :: rst, sfact, dell, einc, del, const, rootv, delhof, tolerf&
        , tolerg, tolrg, tolerx, drop, frepf, cncadd, absmin, alpha, pnorm, &
        cycmx, cos, tx1, tx2, tlast, totime, therb, gnormr, step, funct2, &
        gdnorm, sum, deltag, deltax, ggd, yhy, sy, xvari, ggi, pnlast, dott, &
        beta, smval, dropn, totim, xn, tx, tf, bsmvf, tcycle, tprt
      logical :: restrt, geook, dfp, saddle, minprt, print, thiel, okf, &
        resfil, lgrad
      character :: txt
      double precision, external :: ddot, reada, seconds
      save mdfp, icalcn, rst, sfact, dell, einc, igg1, del, restrt, geook, dfp&
        , const, saddle, minprt, rootv, print, delhof, tolerf, tolerg, nrst, &
        thiel, tolrg, tolerx, drop, frepf, ihdim, cncadd, absmin, itry1, &
        okf, jcyc, lnstop, irepet, alpha, pnorm, jnrst, cycmx, cos, ncount, &
        resfil, tx1, tx2, tlast, totime, maxcyc
!-----------------------------------------------
!
!     *
!     THIS SUBROUTINE ATTEMPTS TO MINIMIZE A REAL-VALUED FUNCTION OF
!     THE N-COMPONENT REAL VECTOR XPARAM ACCORDING TO THE
!     BFGS FORMULA. RELEVANT REFERENCES ARE
!
!     BROYDEN, C.G., JOURNAL OF THE INSTITUTE FOR MATHEMATICS AND
!                     APPLICATIONS, VOL. 6 PP 222-231, 1970.
!     FLETCHER, R., COMPUTER JOURNAL, VOL. 13, PP 317-322, 1970.
!
!
!     SHANNO, D.F. MATHEMATICS OF COMPUTATION, VOL. 24, PP 647-656
!                    1970.
!
!   SEE ALSO SUMMARY IN
!
!    HEAD, J.D.; AND ZERNER, M.C., CHEMICAL PHYSICS LETTERS, VOL. 122,
!          264 (1985).
!    SHANNO, D.F., J. OF OPTIMIZATION THEORY AND APPLICATIONS
!          VOL.46, NO 1 PP 87-94 1985.
!     *
!     THE FUNCTION CAN ALSO BE MINIMIZED USING THE
!     DAVIDON-FLETCHER-POWELL ALGORITHM (COMPUTER JOURNAL, VOL. 6,
!     P. 163).
!
!     THE USER MUST SUPPLY THE SUBROUTINE
!     COMPFG(XPARAM,.TRUE.,FUNCT,.TRUE.,GRAD,LGRAD)
!     WHICH COMPUTES FUNCTION VALUES  FUNCT AT GIVEN VALUES FOR THE
!     VARIABLES XPARAM, AND THE GRADIENT GRAD IF LGRAD=.TRUE.
!     THE MINIMIZATION PROCEEDS BY A SEQUENCE OF ONE-DIMENSIONAL
!     MINIMIZATIONS.  THESE ARE CARRIED OUT WITHOUT GRADIENT COMPUTATION
!     BY THE SUBROUTINE LINMIN, WHICH SOLVES THE SUBPROBLEM OF
!     MINIMIZING THE FUNCTION FUNCT ALONG THE LINE XPARAM+ALPHA*PVECT,
!     WHERE XPARAM
!     IS THE VECTOR OF CURRENT VARIABLE VALUES,  ALPHA IS A SCALAR
!     VARIABLE, AND  PVECT  IS A SEARCH-DIRECTION VECTOR PROVIDED BY THE
!     BFGS OR DAVIDON-FLETCHER-POWELL ALGORITHM.  EACH ITERATION STEP CA
!     OUT BY FLEPO PROCEEDS BY LETTING LINMIN FIND A VALUE FOR ALPHA
!     WHICH MINIMIZES  FUNCT  ALONG  XPARAM+ALPHA*PVECT, BY
!     UPDATING THE VECTOR  XPARAM  BY THE AMOUNT ALPHA*PVECT, AND
!     FINALLY BY GENERATING A NEW VECTOR  PVECT.  UNDER
!     CERTAIN RESTRICTIONS (POWELL, J.INST.MATHS.APPLICS.(1971),
!     V.7,21-36)  A SEQUENCE OF FUNCT VALUES CONVERGING TO SOME
!     LOCAL MINIMUM VALUE AND A SEQUENCE OF
!     XPARAM VECTORS CONVERGING TO THE CORRESPONDING MINIMUM POINT
!     ARE PRODUCED.
!                          CONVERGENCE TESTS.
!
!     HERBERTS TEST: THE ESTIMATED DISTANCE FROM THE CURRENT POINT
!                    POINT TO THE MINIMUM IS LESS THAN TOLERA.
!
!                    "HERBERTS TEST SATISFIED - GEOMETRY OPTIMIZED"
!
!     GRADIENT TEST: THE GRADIENT NORM HAS BECOME LESS THAN TOLERG
!                    TIMES THE SQUARE ROOT OF THE NUMBER OF VARIABLES.
!
!                    "TEST ON GRADIENT SATISFIED".
!
!     XPARAM TEST:  THE RELATIVE CHANGE IN XPARAM, MEASURED BY ITS NORM,
!                   OVER ANY TWO SUCCESSIVE ITERATION STEPS DROPS BELOW
!                   TOLERX.
!
!                    "TEST ON XPARAM SATISFIED".
!
!     FUNCTION TEST: THE CALCULATED VALUE OF THE HEAT OF FORMATION
!                    BETWEEN ANY TWO CYCLES IS WITHIN TOLERF OF
!                    EACH OTHER.
!
!                    "HEAT OF FORMATION TEST SATISFIED"
!
!     FOR THE GRADIENT, FUNCTION, AND XPARAM TESTS A FURTHER CONDITION,
!     THAT NO INDIVIDUAL COMPONENT OF THE GRADIENT IS GREATER
!     THAN TOLERG, MUST BE SATISFIED, IN WHICH CASE THE
!     CALCULATION EXITS WITH THE MESSAGE
!
!                     "PETERS TEST SATISFIED"
!
!     WILL BE PRINTED, AND FUNCT AND XPARAM WILL CONTAIN THE LAST
!     FUNCTION VALUE CUM VARIABLE VALUES REACHED.
!
!
!     ALGORITHMS CHOOSE SEARCH DIRECTIONS
!     ON THE BASIS OF LOCAL PROPERTIES OF THE FUNCTION.  A MATRIX  H,
!     WHICH IN FLEPO IS PRESET WITH THE IDENTITY, IS MAINTAINED AND
!     UPDATED AT EACH ITERATION STEP.  THE MATRIX DESCRIBES A LOCAL
!     METRIC ON THE SURFACE OF FUNCTION VALUES ABOVE THE POINT XPARAM.
!     THE SEARCH-DIRECTION VECTOR  PVECT  IS SIMPLY A TRANSFORMATION
!     OF THE GRADIENT  GRAD  BY THE MATRIX H.
!
      equivalence (mdfp(1), jcyc), (mdfp(2), jnrst), (mdfp(3), ncount), (mdfp(4&
        ), lnstop), (xdfp(1), alpha), (xdfp(2), cos), (xdfp(3), pnorm), (xdfp(4&
        ), drop), (xdfp(5), del), (xdfp(6), frepf), (xdfp(7), cycmx), (xdfp(8)&
        , totime)
      data icalcn/ 0/
!
!   START OF ONCE-ONLY SECTION
!
      emin = 0.D0
      smval = 0.d0
      gnormr = 0.d0
      beta = 0.d0
      therb = 0.d0
      icyc = 0
      if (icalcn /= numcal) then
        if (allocated(hesinv)) deallocate(hesinv)
        allocate (hesinv((nvar*(nvar+1))/2))
        hesinv = 0.0d0
!
!   THE FOLLOWING CONSTANTS SHOULD BE SET BY THE USER.
!
        rst = -1.D0
        nrst = 30
        sfact = 1.5D0
        dell = 0.01D0
        einc = 0.3D0
        igg1 = 3
        del = dell
!
!    THESE CONSTANTS SHOULD BE SET BY THE PROGRAM.
!
        maxcyc = 100000
        if (index(keywrd,' CYCLES') /= 0) maxcyc = nint(reada(keywrd,index(&
          keywrd,' CYCLES')))
        restrt = index(keywrd,'RESTAR') /= 0
        thiel = index(keywrd,'NOTHIE') == 0
        geook = index(keywrd,'GEO-OK') /= 0
        saddle = index(keywrd,'SADDLE') /= 0
        minprt = .not.saddle
        const = 1.D0
!
!      THE DAVIDON-FLETCHER-POWELL METHOD IS NOT RECOMMENDED
!      BUT CAN BE INVOKED BY USING THE KEYWORD 'DFP'
!
        dfp = index(keywrd,'DFP') /= 0
!
!  ORDER OF PRECISION:   'GNORM' TAKES PRECEDENCE OVER 'FORCE', WHICH
!                        TAKES PRECEDENCE OVER 'PRECISE'.
        tolerg = 1.0D0
        if (index(keywrd,'PREC') /= 0) tolerg = 0.2D0
        if (index(keywrd,'FORCE') /= 0) tolerg = 0.1D0
!
!      READ IN THE GRADIENT-NORM LIMIT, IF SPECIFIED
!
        if (index(keywrd,'GNORM=') /= 0) then
          rootv = 1.D0
          const = 1.D-20
          tolerg = reada(keywrd,index(keywrd,'GNORM='))
          if (index(keywrd,' LET')==0 .and. tolerg<1.D-2) then
            write (iw, '(/,A)') '  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01'
            tolerg = 1.D-2
          end if
        else
          rootv = sqrt(nvar + 1.D-5)
        end if
        tolerx = 0.0001D0*const
        delhof = 0.0010D0*const
        tolerf = 0.002D0*const
        tolrg = tolerg
        if (index(keywrd,'PREC') /= 0) then
          tolerx = tolerx*0.01D0
          delhof = delhof*0.01D0
          tolerf = tolerf*0.01D0
          einc = einc*0.01D0
        end if
        therb = 5.D0*tolerg*rootv
!
!  MINOR BOOK-KEEPING
!
        tlast = tleft
        tx2 = seconds(2)
        tleft = tleft - tx2 + time0
        print = index(keywrd,'FLEPO') /= 0
!
!   THE FOLLOWING CONSTANTS SHOULD BE SET TO SOME ARBITARY LARGE VALUE.
!
        drop = 1.D15
        frepf = 1.D15
!
!     AND FINALLY, THE FOLLOWING CONSTANTS ARE CALCULATED.
!
        ihdim = (nvar*(nvar + 1))/2
        cncadd = 1.0D00/rootv
        cncadd = min(0.15D00,cncadd)
        icalcn = numcal
        if (restrt .and. nvar > 0) then
          jnrst = 1
          mdfp(9) = 0
          gd = 0.d0
          xlast = 0.d0
          call dfpsav (totime, xparam, gd, xlast, funct1, mdfp, xdfp)
          if (moperr) return
          jcyc = jcyc - 1
          i = int(totime/1000000.D0)
          totime = totime - i*1000000.D0
          time0 = time0 - totime
          nscf = mdfp(5)
          write (iw, &
      '(2/10X,''TOTAL TIME USED SO FAR:'',                   F13.2,'' SECONDS''&
      &)') totime
          if (index(keywrd,' 1SCF') /= 0) then
            last = 1
            lgrad = index(keywrd,' GRAD') /= 0
            call compfg (xparam, .TRUE., funct1, .TRUE., grad, lgrad)
            if (moperr) return
            iflepo = 13
            emin = 0.D0
            return
          end if
        end if
!
!   END OF ONCE-ONLY SETUP
!
      end if
!
!     FIRST, WE INITIALIZE THE VARIABLES.
!
      ireset = 0
      absmin = 1.D6
      itry1 = 0
      jcyc = 0
      lnstop = 1
      irepet = 1
      limscf = .TRUE.
      alpha = 1.0D00
      pnorm = 1.0D00
      jnrst = 0
      cycmx = 0.D0
      cos = 0.0D00
      totime = 0.D0
      ncount = 1
      resfil = .FALSE.
      if (const>1.D-5 .and. saddle) then
!
!   WE DON'T NEED HIGH PRECISION DURING A SADDLE-POINT CALCULATION.
!
! For Mopac BLAS
!        if (nvar > 0) gnorm = sqrt(dot(grad,grad,nvar)) - 3.D0
        if (nvar > 0) gnorm = dsqrt(ddot(nvar,grad,1,grad,1)) - 3.D0

!
        gnorm = min(10.D0,gnorm)
        if (gnorm > 1.D0) tolerg = tolrg*gnorm
        write (iw, '('' GRADIENT CRITERION IN FLEPO ='',F10.3)') tolerg
      end if
      if (nvar == 1) then
        pvect(1) = 0.01D0
        alpha = 1.D0
        go to 270
      end if
      totime = 0.D0
!
! CALCULATE THE VALUE OF THE FUNCTION -> FUNCT1, AND GRADIENTS -> GRAD.
! NORMAL SET-UP OF FUNCT1 AND GRAD, DONE ONCE ONLY.
!
      call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
      iflepo = 16
      if (nvar == 0) return
      if (moperr) return
      call dcopy (nvar, grad, 1, gd, 1)
      if (nvar /= 0) then
! For Mopac BLAS
!        gnorm = sqrt(dot(grad,grad,nvar))
         gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
!
        gnormr = gnorm
        if (lnstop/=1 .and. cos>rst .and. (jnrst<nrst .or. .not.dfp) .and. &
          restrt) then
          call dcopy (nvar, gd, 1, glast, 1)
        else
          call dcopy (nvar, grad, 1, glast, 1)
        end if
      end if
      if (gnorm<tolerg .or. nvar==0) then
        iflepo = 2
        if (restrt) then
          call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
          if (moperr) return
        else
          call compfg (xparam, .TRUE., funct1, .TRUE., grad, .FALSE.)
          if (moperr) return
        end if
        tx2 = seconds(1)
        emin = 0.D0
        return
      end if
      tx1 = seconds(2)
      tleft = tleft - tx1 + tx2
!     *
!     START OF EACH ITERATION CYCLE ...
!     *
!
! For Mopac BLAS
!        gnorm = sqrt(dot(grad,grad,nvar))
         gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
!

      if (gnormr < 1.D-10) gnormr = gnorm
      icyc = jcyc
      go to 30
   10 continue
      if (cos < rst) then
        gd(:nvar) = 0.5D0
      end if
   30 continue
      jcyc = jcyc + 1
      jnrst = jnrst + 1
      i80 = 0
      if (i80==1 .or. lnstop==1 .or. cos<=rst .or. jnrst>=nrst .and. dfp) then
!
!     *
!     RESTART SECTION
!     *
!
   50   continue
        do i = 1, nvar
!
!  MAKE THE FIRST STEP A WEAK FUNCTION OF THE GRADIENT
!
          step = abs(grad(i))*0.0002D0
          step = max(0.01D0,min(0.04D0,step))
          xd(i) = xparam(i) - sign(del,grad(i))
        end do
!
! THIS CALL OF COMPFG IS USED TO CALCULATE THE SECOND-ORDER MATRIX IN H
! IF THE NEW POINT HAPPENS TO IMPROVE THE RESULT, THEN IT IS KEPT.
! OTHERWISE IT IS SCRAPPED, BUT STILL THE SECOND-ORDER MATRIX IS O.K.
!
        call compfg (xd, .TRUE., funct2, .TRUE., gd, .TRUE.)
        if (moperr) return
        if (.not.geook .and. dsqrt(ddot(nvar,gd,1,gd,1))/gnorm>10.D0 .and. gnorm>20&
           .and. jcyc>2) then
          del = del/10.0D0
          if (del >= 0.00005D0) then
            go to 50
          else
!
!  THE GEOMETRY IS BADLY SPECIFIED IN THAT MINOR CHANGES IN INTERNAL
!  COORDINATES LEAD TO LARGE CHANGES IN CARTESIAN COORDINATES, AND THESE
!  LARGE CHANGES ARE BETWEEN PAIRS OF ATOMS THAT ARE CHEMICALLY BONDED
!  TOGETHER.
            write (iw, '('' GRADIENTS OF OLD GEOMETRY, GNORM='',       F13.6)')&
               gnorm
            write (iw, '(6F12.6)') (grad(i),i=1,nvar)
            gdnorm = dsqrt(ddot(nvar,gd,1,gd,1))
            write (iw, '('' GRADIENTS OF NEW GEOMETRY, GNORM='',       F13.6)')&
               gdnorm
            write (iw, '(6F12.6)') (gd(i),i=1,nvar)
            write (iw, &
      '(3/20X,                                       ''CALCULATION ABANDONED AT&
      & THIS POINT!'')')
            write (iw, &
      '(2/10X,'' SMALL CHANGES IN INTERNAL '',2/       ''COORDINATES ARE   '',/&
      &10X,                                      '' CAUSING A LARGE CHANGE IN TH&
      &E DISTANCE BETWEEN'',              /   10X,'' CHEMICALLY-BOUND ATOMS. THE&
      & GEOMETRY OPTIMIZATION'',/     10X,'' PROCEDURE WOULD LIKELY PRODUCE INCO&
      &RRECT RESULTS'')')
            call geout (1)
            call mopend ('CALCULATION ABANDONED IN FLEPO')
            return
          end if
        end if
        ncount = ncount + 1
        hesinv(:ihdim) = 0.0D00
        sum = 0.D0
        ii = 0
        j = 0
        if (funct2 < funct1) then
          do i = 1, nvar
            ii = ii + i
            deltag = grad(i) - gd(i)
            deltax = xparam(i) - xd(i)
            if (abs(deltag) < 0.001D0) deltag = 0.001D0
            ggd = max(1.D0,abs(gd(i)))
            hesinv(ii) = min(0.1D0/ggd,deltax/deltag)
            if (hesinv(ii) <= 0.0D00) cycle
            j = j + 1
            sum = sum + hesinv(ii)
          end do
        else
          do i = 1, nvar
            ii = ii + i
            deltag = grad(i) - gd(i)
            deltax = xparam(i) - xd(i)
            if (abs(deltag) < 0.001D0) deltag = 0.001D0
            ggd = max(1.D0,abs(grad(i)))
            hesinv(ii) = min(0.1D0/ggd,deltax/deltag)
            if (hesinv(ii) <= 0.0D00) cycle
            j = j + 1
            sum = sum + hesinv(ii)
          end do
        end if
        if (j /= 0) then
          sum = sum/j
          ii = 0
          do i = 1, nvar
            ii = ii + i
            if (hesinv(ii) >= 0) cycle
            hesinv(ii) = sum
          end do
        end if
        jnrst = 0
        if (jcyc < 2) cosine = 1.D0
        if (funct2 >= funct1) then
          if (print) write (iw, 100) funct1, funct2
  100     format(' FUNCTION VALUE=',f13.7,'  WILL NOT BE REPLACED BY VALUE=',&
            f13.7,/,10x,'CALCULATED BY RESTART PROCEDURE',/)
          cosine = 1.D0
        else
          if (print) write (iw, 110) funct1, funct2
  110     format(' FUNCTION VALUE=',f13.7,' IS BEING REPLACED BY VALUE=',f13.7,&
            /,10x,' FOUND IN RESTART PROCEDURE',/,6x,'THE CORRESPONDING',&
            ' X VALUES AND GRADIENTS ARE ALSO BEING REPLACED',/)
          funct1 = funct2
          call dcopy (nvar, xd, 1, xparam, 1)
          call dcopy (nvar, gd, 1, grad, 1)
! For Mopac BLAS
!        gnorm = sqrt(dot(grad,grad,nvar))
         gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
!
          if (gnormr < 1.D-10) gnormr = gnorm
        end if
      else
!
!     *
!     UPDATE VARIABLE-METRIC MATRIX
!     *
!
        xvar(:nvar) = xparam(:nvar) - xlast(:nvar)
        gvar(:nvar) = grad(:nvar) - glast(:nvar)
        call supdot (gg, hesinv, gvar, nvar)

! For Mopac BLAS
!        yhy = dot(gg,gvar,nvar)
!        sy = dot(xvar,gvar,nvar)
        yhy = ddot(nvar,gg,1,gvar,1)
        sy = ddot(nvar,xvar,1,gvar,1)
!

        k = 0
!
!    UPDATE ACCORDING TO DAVIDON-FLETCHER-POWELL
!
        if (dfp) then
          do i = 1, nvar
            xvari = xvar(i)/sy
            ggi = gg(i)/yhy
            if (i > 0) then
              hesinv(k+1:i+k) = hesinv(k+1:i+k) + xvar(:i)*xvari - gg(:i)*ggi
              k = i + k
            end if
          end do
!
!     UPDATE USING THE BFGS FORMALISM
!
        else
          yhy = 1.0D0 + yhy/sy
          do i = 1, nvar
            xvari = xvar(i)/sy
            ggi = gg(i)/sy
            if (i > 0) then
              hesinv(k+1:i+k) = hesinv(k+1:i+k) - gg(:i)*xvari - xvar(:i)*ggi&
                 + yhy*xvar(:i)*xvari
              k = i + k
            end if
          end do
        end if
      end if
!
!     *
!     ESTABLISH NEW SEARCH DIRECTION
!     *
      pnlast = pnorm
      call supdot (pvect, hesinv, grad, nvar)
      pnorm = dsqrt(ddot(nvar,pvect,1,pvect,1))
      if (pnorm > 1.5D0*pnlast) then
!
!  TRIM PVECT BACK
!
        pvect(:nvar) = pvect(:nvar)*1.5D0*pnlast/pnorm
        pnorm = 1.5D0*pnlast
      end if
! For Mopac BLAS
!     dott = -dot(pvect,grad,nvar)
      dott = -ddot(nvar,pvect,1,grad,1)
!
      pvect(:nvar) = -pvect(:nvar)
      cos = -dott/(pnorm*gnorm)
      if (print) write (iw, 170) jcyc, funct1
  170 format(' ','AT THE BEGINNING OF CYCLE',i5,'  THE FUNCTION',' VALUE IS ',&
        f13.6,/,'  THE CURRENT POINT IS ...')
      if (print) write (iw, 180) gnorm, cos
  180 format('  GRADIENT NORM = ',f10.4,/,'  ANGLE COSINE =',f10.4)
      if (print) then
        write (iw, 190)
  190   format('  THE CURRENT POINT IS ...')
        nto6 = (nvar - 1)/6 + 1
        iinc1 = -5
        do i = 1, nto6
          write (iw, '(/)')
          iinc1 = iinc1 + 6
          iinc2 = min(iinc1 + 5,nvar)
          write (iw, 200) (j,j=iinc1,iinc2)
          write (iw, 210) (xparam(j),j=iinc1,iinc2)
          write (iw, 220) (grad(j),j=iinc1,iinc2)
          write (iw, 230) (pvect(j),j=iinc1,iinc2)
  200     format(' ',3x,'I',9x,i3,9(8x,i3))
  210     format(' ',1x,'XPARAM(I)',1x,f9.4,2x,9(f9.4,2x))
  220     format(' ',1x,'GRAD  (I)',f10.4,1x,9(f10.4,1x))
  230     format(' ',1x,'PVECT (I)',2x,f10.6,1x,9(f10.6,1x))
        end do
      end if
      lnstop = 0
      alpha = alpha*pnlast/pnorm
      call dcopy (nvar, grad, 1, glast, 1)
      call dcopy (nvar, xparam, 1, xlast, 1)
      if (jnrst == 0) alpha = 1.0D00
      drop = abs(alpha*dott)
      if (print) write (iw, 250) drop
  250 format(' ',' -ALPHA.P.G =',f18.6,/)
      if (gnorm<=therb .and. jnrst/=0 .and. drop<delhof) then
!
!   HERBERT'S TEST: THE PREDICTED DROP IN ENERGY IS LESS THAN DELHOF
!   IF PASSED, CALL COMPFG TO GET A GOOD SET OF EIGENVECTORS, THEN EXIT
!
        if (minprt) write (iw, 260)
  260   format(/,/,10x,'HERBERTS TEST SATISFIED - GEOMETRY OPTIMIZED')
!
!   FLEPO IS ENDING PROPERLY. THIS IS IMMEDIATELY BEFORE THE RETURN.
!
        last = 1
        call compfg (xparam, .TRUE., funct1, .TRUE., grad, .FALSE.)
        if (moperr) return
        iflepo = 3
        time0 = time0 - totime
        emin = 0.D0
        tx2 = seconds(1)
        return
      end if
      beta = alpha
      smval = funct1
      dropn = -abs(drop/alpha)
!
!    UPDATE GEOMETRY USING THE G-DIIS PROCEDURE
!
      okf = .FALSE.
      ic = 2
  270 continue
      call linmin (xparam, alpha, pvect, nvar, funct1, okf, ic, dropn)
      if (moperr) return
      if (nvar == 1) then
        write (iw, &
      '('' ONLY ONE VARIABLE, THEREFORE ENERGY A MINIMUM'')')
        last = 1
        lgrad = index(keywrd,'GRAD') /= 0
        grad(1) = 0.D0
!
!  This call is essential.  If this call is not made, then
!  the eigenvalues will be incorrect.
!
        call compfg (xparam, .TRUE., funct1, .TRUE., grad, lgrad)
        if (moperr) return
        iflepo = 14
        emin = 0.D0
        tx2 = seconds(1)
        if (tx2 - time0 > tleft) then
!
!  This step is necessary so that a GRID calculation can
!  be re-started, even if there is only one variable.
!
          iflepo = -1
          totim = totime + seconds(1) - time0
          mdfp(9) = 1
          mdfp(5) = nscf
          call dfpsav (totim, xparam, gd, xlast, funct1, mdfp, xdfp)
          if (moperr) return
        end if
        return
      end if
!   WE WANT ACCURATE DERIVATIVES AT THIS POINT
!
!   LINMIN DOES NOT GENERATE ANY DERIVATIVES, THEREFORE COMPFG MUST BE
!   CALLED TO END THE LINE SEARCH
!
!  IF THE DERIVATIVES ARE TO BE CALCULATED USING FULL SCF'S, THEN CHECK
!  WHETHER TO DO FULL SCF'S (CRITERION FROM FLEPO: GRAD IS NULL).
!
      if (ireset>10 .or. gnorm<40.D0 .and. gnorm/gnormr<0.33D0) then
        ireset = 0
        gnormr = 0.D0
        grad(:nvar) = 0.D0
      end if
      ireset = ireset + 1
!
!
!     RESTORE TO STANDARD VALUE BEFORE COMPUTING THE GRADIENT
      if (thiel) then
        call compfg (xparam, ic/=1, sum, .TRUE., grad, .TRUE.) ! sum is not used
        if (moperr) return
      else
        call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
        if (moperr) return
      end if
! For Mopac BLAS
      !gnorm = sqrt(dot(grad,grad,nvar))
      gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
!
      if (gnormr < 1.D-10) gnormr = gnorm
      ncount = ncount + 1
      if (.not.okf) then
        lnstop = 1
        if (minprt) write (iw, &
      '(/,20X, ''NO POINT LOWER IN ENERGY THAN THE STARTING POINT&
      & '',/,20X,''COULD BE FOUND IN THE LINE MINIMIZATION'')')
        funct1 = smval
        alpha = beta
        call dcopy (nvar, glast, 1, grad, 1)
        call dcopy (nvar, xlast, 1, xparam, 1)
        if (jnrst == 0) then
          write (iw, 290)
  290     format(' ',/,/,20x,'SINCE COS WAS JUST RESET,THE SEARCH',&
            ' IS BEING ENDED')
!
!           FLEPO IS ENDING BADLY. THIS IS IMMEDIATELY BEFORE THE RETURN
!
          last = 1
          call compfg (xparam, .TRUE., funct1, .TRUE., grad, .TRUE.)
          if (moperr) return
          iflepo = 4
          time0 = time0 - totime
          tx2 = seconds(1)
          emin = 0.D0
          return
        end if
        if (print) write (iw, 300)
  300   format(' ',20x,'COS WILL BE RESET AND ANOTHER ','ATTEMPT MADE')
        cos = 0.0D00
        go to 430
      end if
! For Mopac BLAS
!      xn = sqrt(dot(xparam,xparam,nvar))
      xn = dsqrt(ddot(nvar,xparam,1,xparam,1))
!
      tx = dabs(alpha*pnorm)
      if (xn /= 0.0D00) tx = tx/xn
      tf = dabs(smval - funct1)
      if (absmin - smval < 1.D-7) then
        itry1 = itry1 + 1
        if (itry1 > 10) then
          write (iw, &
      '(2/,'' HEAT OF FORMATION IS ESSENTIALLY STATIONARY'')')
          go to 420
        end if
      else
        itry1 = 0
        absmin = smval
      end if
      if (print) write (iw, 310) ncount, cos, tx*xn, alpha, (-drop), (-tf), &
        & gnorm
  310 format(/,'           NUMBER OF COUNTS =',i6,'         COS    =',f11.4,/,&
        '  ABSOLUTE  CHANGE IN X     =',f13.6,'  ALPHA  =',f11.4,/,&
        '  PREDICTED CHANGE IN F     =  ',g11.4,'  ACTUAL =  ',g11.4,/,&
        '  GRADIENT NORM             =  ',g11.4,/,/)
      if (tx <= tolerx) then
        if (minprt) write (iw, 320)
  320   format(' TEST ON X SATISFIED')
        go to 350
      end if
      if (tf <= tolerf) then
        if (minprt) write (iw, 330)
  330   format(' HEAT OF FORMATION TEST SATISFIED')
        go to 350
      end if
      if (gnorm <= tolerg*rootv) then
        if (minprt) write (iw, 340)
  340   format(' TEST ON GRADIENT SATISFIED')
        go to 350
      end if
      go to 430
  350 continue
      do i = 1, nvar
        if (abs(grad(i)) <= tolerg) cycle
        irepet = irepet + 1
        if (irepet <= 1) then
          frepf = funct1
          cos = 0.0D00
        end if
        if (minprt) write (iw, 370) tolerg
  370   format(20x,'HOWEVER, A COMPONENT OF GRADIENT IS ','LARGER THAN',f6.2,/)
        if (abs(funct1 - frepf) > einc) irepet = 0
        if (irepet > igg1) then
          write (iw, 380) igg1, einc
  380     format(10x,' THERE HAVE BEEN',i2,' ATTEMPTS TO ',&
            'REDUCE THE GRADIENT.',/,10x,&
            ' DURING THESE ATTEMPTS THE ENERGY DROPPED',' BY LESS THAN',f4.1,&
            ' KCAL/MOLE',/,10x,&
            ' FURTHER CALCULATION IS NOT JUSTIFIED AT THIS TIME.')
          if (index(keywrd,'PREC') == 0) write (iw, 390)
  390     format(10x,' TO CONTINUE, START AGAIN WITH THE WORD "PRECISE"')
          last = 1
          call compfg (xparam, .TRUE., funct1, .TRUE., grad, .FALSE.)
          if (moperr) return
          iflepo = 8
          time0 = time0 - totime
          tx2 = seconds(1)
          emin = 0.D0
          return
        else
          go to 430
        end if
      end do
      if (minprt) write (iw, 410)
  410 format('PETERS TEST SATISFIED')
  420 continue
      last = 1
      call compfg (xparam, .TRUE., funct1, .TRUE., grad, .FALSE.)
      if (moperr) return
      iflepo = 6
      time0 = time0 - totime
      tx2 = seconds(1)
      emin = 0.D0
      return
!
!   ALL TESTS HAVE FAILED, WE NEED TO DO ANOTHER CYCLE.
!
  430 continue
      bsmvf = abs(smval - funct1)
      if (bsmvf > 10.D00) cos = 0.0D00
      del = 0.002D00
      if (bsmvf > 1.0D00) del = dell/2.0D00
      if (bsmvf > 5.0D00) del = dell
      tx2 = seconds(2)
      tcycle = tx2 - tx1
      tx1 = tx2
!
! END OF ITERATION LOOP, EVERYTHING IS STILL O.K. SO GO TO
! NEXT ITERATION, IF THERE IS ENOUGH TIME LEFT.
!
      if (tcycle < 100000.D0) cycmx = max(cycmx,tcycle)
      tleft = tleft - tcycle
      if (tleft < 0) tleft = -0.1D0
      if (tcycle > 1.D5) tcycle = 0.D0
      if (tlast - tleft > tdump) then
        totim = totime + seconds(1) - time0
        tlast = tleft
        mdfp(9) = 2
        resfil = .TRUE.
        mdfp(5) = nscf
        call dfpsav (totim, xparam, gd, xlast, funct1, mdfp, xdfp)
        if (moperr) return
      end if
      call prttim (tleft, tprt, txt)
      if (resfil) then
        write (line, '(" RESTART FILE WRITTEN,      TIME LEFT:", f6.2, &
           & a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
           tprt, txt, Min (gnorm, 999999.999d0), funct1
        write(iw,"(a)")trim(line)
        call to_screen(trim(line))
        endfile (iw)
        backspace (iw)
        if (log) then
          write (ilog, "(a)")trim(line)
          endfile (ilog)
          backspace (ilog)
        end if
        resfil = .false.
      else
        write (line, '(" CYCLE:", i6, " TIME:", f8.3, " TIME LEFT:", &
          & f6.2, a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
          jcyc, Min (tcycle, 9999.99d0), tprt, txt, &
          & Min (gnorm, 999999.999d0), funct1
          if (minprt) then
            write(iw,"(a)")trim(line)
            endfile (iw)
            backspace (iw)
          end if
          if (log) then
            write (ilog, "(a)")trim(line)
            endfile (ilog)
            backspace (ilog)
          end if
          call to_screen(trim(line))
      end if
      if (mod(jcyc,30) == 0) then
        line = trim(input_fn)
        call add_path(line)
        i = len_trim(line) - 5
        call to_screen(line(:i))
      end if
      call to_screen("To_file: Geometry optimizing")
      if (tleft>sfact*cycmx .and. jcyc-icyc<maxcyc) go to 10
      write (iw, 460)
  460 format(20x,'THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE',/,30x,&
        'NOW GOING TO FINAL')
      totim = totime + seconds(1) - time0
      mdfp(9) = 1
      mdfp(5) = nscf
      call dfpsav (totim, xparam, gd, xlast, funct1, mdfp, xdfp)
      if (moperr) return
      iflepo = -1
      return
      end subroutine flepo
