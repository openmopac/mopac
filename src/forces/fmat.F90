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

      subroutine fmat(fmatrx, nreal, tscf, tder, deldip, heat, evecs, ts)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use common_arrays_C, only : coord, na, labels, atmass, p, q, xparam
      use molkst_C, only : numat, nvar, keywrd, tleft, moperr, natoms, tdump, &
      emin, mozyme
      use parameters_C, only : tore
      use funcon_C, only : fpc_10
      use chanel_C, only : iw, ilog, log
      USE elemts_C, only : elemnt
      implicit none
      integer , intent(out) :: nreal
      logical, intent (in) :: ts
      double precision , intent(in) :: tscf
      double precision , intent(in) :: tder
      double precision  :: heat
      double precision  :: fmatrx(((3*numat + 1)*(3*numat))/2)
      double precision  :: deldip(3,3*numat)
      double precision  :: evecs(9*numat*numat)
!
      integer ::  i, lin, maxcyc, istart, jstart, kountf, lu, iskip, j, &
        ii, ll, l, k, kk, iloop, store_numat
      double precision, dimension(3*natoms) :: grad, grold
      double precision, dimension(3) :: del2
      double precision, dimension(3*natoms) :: g2old, eigs, g2rad, fconst, dumy
      double precision :: fact, tlast, totime, time2, estime, delta, &
        time3, tstep, estim, escf
      logical :: debug, restrt, prnt, resfil, precis, big
      character :: line*120
      double precision, external :: dipole_for_MOZYME
      double precision, external :: ddot, dipole, dot, reada, seconds
      character :: XYZ(3)*1, num
      save fact
!**********************************************************************
!
!  VALUE CALCULATES THE SECOND-ORDER OF THE ENERGY WITH
!        RESPECT TO THE CARTESIAN COORDINATES I AND J AND PLACES IT
!        IN FMATRX
!
!  ON INPUT NATOMS  = NUMBER OF ATOMS IN THE SYSTEM.
!           XPARAM  = INTERNAL COORDINATES OF MOLECULE STORED LINEARLY
!
!
!  ON OUTPUT FMATRX = SECOND DERIVATIVE OF THE ENERGY WITH RESPECT TO
!                    CARTESIAN COORDINATES I AND J.
!**********************************************************************
      data del2/ 3*0.D0/
      data xyz /"X", "Y", "Z"/
      num = char(Int(log10(numat*1.0)) + ichar("1") + 1)
!
!    FACT IS THE CONVERSION FACTOR FROM KCAL/MOLE TO ERGS
!
      fact = 4.184D0/fpc_10*1.D21
!          = 4.184/6.0221367d23*10^21
!
! SET UP CONSTANTS AND FLAGS
      na = 0
!
!  SET UP THE VARIABLES IN XPARAM AND LOC,THESE ARE IN
!  CARTESIAN COORDINATES
!
      numat = 0
      do i = 1, natoms
        if (labels(i)==99 .or. labels(i)==107) cycle
        numat = numat + 1
        labels(numat) = labels(i)
      end do
      natoms = numat
      grad = 0.d0
      lin = (nvar*(nvar + 1))/2
      fmatrx(:lin) = 0.D0
      maxcyc = 100000
      if (index(keywrd,' CYCLES') /= 0 .and. index(keywrd, " IRC") == 0) &
        maxcyc = nint(reada(keywrd,index(keywrd,' CYCLES')))
      prnt = index(keywrd,'IRC=') == 0
      precis = index(keywrd,' PREC') /= 0
      restrt = index(keywrd,'RESTART') /= 0
      if (index(keywrd,'NLLSQ') /= 0) restrt = .FALSE.
      debug = index(keywrd,'FMAT') /= 0
      big = (index(keywrd,'LARGE') /= 0 .and. debug)
      if (prnt) write (iw, &
      '(2/4X,''FIRST DERIVATIVES WILL BE USED IN THE CALCULATION OF SECOND DERIVATIVES'')')
      tlast = tleft
      resfil = .FALSE.
      if (restrt) then
        istart = 0
        i = 0
        totime = 0.d0
        jstart = 0
        fconst = 0.d0
        call forsav (totime, deldip, istart, fmatrx, coord, nvar, heat, evecs, jstart, fconst)
        if (moperr) return
        kountf = (istart*(istart + 1))/2
        istart = istart + 1
        jstart = jstart + 1
        time2 = seconds(1)
      else
        kountf = 0
        totime = 0.D0
        if (tscf > 0.D0) tleft = tleft - tscf - tder
        istart = 1
      end if
! CALCULATE FMATRX
      if (istart > 1) then
        estime = (nvar - istart + 1)*totime/(istart - 1.D0)
      else
        estime = nvar*(tscf + tder)*2.D0
        if (precis) estime = estime*2.D0
      end if
! 20190807 removed conditional output to stabilize output file
!      if (tscf > 0)
      write (iw, &
      '(/10X,''ESTIMATED TIME TO COMPLETE CALCULATION ='',F12.2,'' SECONDS'')') estime
      if (restrt) then
        if (istart <= nvar) write (iw, &
          '(/10X,''STARTING AGAIN AT LINE'',18X,I4)') istart
        write (iw, '(/10X,''TIME USED UP TO RESTART ='',F22.2)') totime
      end if
      lu = kountf
      eigs(:nvar) = 0.D0
      if (.not. ts) call symr
      iskip = 0
      do i = istart, nvar
        if (((i - 1)/3)*3 == i - 1) then
!
!  START OF A NEW ATOM.  DOES A SYMMETRY OPERATION RELATE AN ALREADY
!  CALCULATED ATOM TO THIS ONE
!
          j = (i + 2)/3
          call sympop (fmatrx, j, iskip, deldip)
        end if
        if (iskip > 0) then
          write (line, '(" STEP:",I4,23X,"INTEGRAL =",F10.2," TIME LEFT:",F10.2)') i, totime, tleft
          write(iw,"(a)")line(:len_trim(line))
          call to_screen(line)
          endfile (iw)
          backspace (iw)
          iskip = iskip - 1
          lu = lu + i
          cycle
        end if
        time2 = seconds(1)
        delta = 1.D0/120.D0
        g2old = 0.D0
        if (precis) then
!
!   DETERMINE A GOOD STEP SIZE
!
          g2old = 100.D0
          xparam(i) = xparam(i) + delta
          call compfg (xparam, .TRUE., escf, .TRUE., g2old, .TRUE.)
          if (moperr) return
          xparam(i) = xparam(i) - delta
          delta = delta*10.D0/sqrt(ddot(nvar,g2old,1,g2old,1))
!
!   CONSTRAIN DELTA TO A 'REASONABLE' VALUE
!
          delta = min(0.05D0,max(0.005D0,delta))
          if (debug) write (iw, '(A,I3,A,F15.8)') ' STEP:', i, ' DELTA :', delta
          g2old(1) = 100.D0
          xparam(i) = xparam(i) + delta
          emin = 0.d0
          call compfg (xparam, .TRUE., escf, .TRUE., g2old, .TRUE.)
          if (moperr) return
!
          if (debug) write (iw, '(A,F12.5)') ' GNORM +1.0*DELTA',  &
            & dsqrt(ddot(nvar,g2old,1,g2old,1))
          xparam(i) = xparam(i) - delta*2.D0
          g2rad = 100.D0
          emin = 0.d0
          call compfg (xparam, .TRUE., escf, .TRUE., g2rad, .TRUE.)
          if (moperr) return
          xparam(i) = xparam(i) + delta
          if (debug) write (iw, '(A,F12.5)') ' GNORM -1.0*DELTA', &
            dsqrt(ddot(nvar,g2rad,1,g2rad,1))
        else
          if (debug) write (iw, '(A,I3,A,F15.8)') ' STEP:', i, ' DELTA :', delta
        end if
        xparam(i) = xparam(i) + 0.5D0*delta
        grold = 100.D0
        emin = 0.d0
        call compfg (xparam, .TRUE., escf, .TRUE., grold, .TRUE.)
        if (moperr) return
        if (debug) then
          if (precis) then
            write (iw, '(A,F12.5)') ' GNORM +0.5*DELTA', dsqrt(ddot(nvar,grold,1,grold,1))
          else
            write (iw, '(A,F12.5)') ' GNORM +DELTA:', dsqrt(ddot(nvar,grold,1,grold,1))
          end if
        end if
!
        call chrge (p, q)

        q(:numat) = tore(labels(:numat)) - q(:numat)
!
!   ESTIME IS USED HERE AS DIPOLE IS A FUNCTION
!
        store_numat = numat
        numat = nvar/3
        if (mozyme) then
          estime = dipole_for_MOZYME (deldip(1,i),0)
        else
          estime = dipole(p,xparam,deldip(1,i),0)
        end if
        numat = store_numat
        xparam(i) = xparam(i) - delta
        grad(1) = escf ! dummy use of escf and grad
        emin = 0.d0
        call compfg (xparam, .TRUE., escf, .TRUE., grad, .TRUE.)
        if (moperr) return
        if (debug) then
          if (precis) then
            write (iw, '(A,F12.5)') ' GNORM -0.5*DELTA', sqrt(dot(grad,grad,nvar))
          else
            write (iw, '(A,F12.5)') ' GNORM -DELTA:', sqrt(dot(grad,grad,nvar))
          end if
        end if
        call chrge (p, q)
        q(:numat) = tore(labels(:numat)) - q(:numat)
!
!   ESTIME IS USED HERE AS DIPOLE IS A FUNCTION
!
        store_numat = numat
        numat = nvar/3
        if (mozyme) then
          estime = dipole_for_MOZYME (del2, 0)
        else
          estime = dipole(p,xparam,del2, 0)
        end if
        numat = store_numat
        xparam(i) = xparam(i) + delta*0.5D0
        deldip(1,i) = (deldip(1,i)-del2(1))*0.5D0/delta
        deldip(2,i) = (deldip(2,i)-del2(2))*0.5D0/delta
        deldip(3,i) = (deldip(3,i)-del2(3))*0.5D0/delta
        ll = lu + 1
        lu = ll + i - 1
        l = 0
        if (precis) then
          if (lu - ll + 1 > 0) then
!
!       G2OLD = X + 1.0*DELTA
!       GROLD = X + 0.5*DELTA
!       GRAD  = X - 0.5*DELTA
!       G2RAD = X - 1.0*DELTA
!
            dumy(:lu-ll+1) = (8.D0*(grold(:lu-ll+1)-grad(:lu-ll+1)) - &
              (g2old(:lu-ll+1)-g2rad(:lu-ll+1)))/delta*fact/12.D0
            eigs(:lu-ll+1) = (2.D0*(grold(:lu-ll+1)-grad(:lu-ll+1)) - &
              (g2old(:lu-ll+1)-g2rad(:lu-ll+1)))/delta**3*fact/28.D0
!
!  CORRECT FOR 4'TH ORDER CONTAMINATION
!
            fmatrx(ll:lu) = fmatrx(ll:lu) + dumy(:lu-ll+1)
            l = lu - ll + 1
          end if
          l = l - 1
          do k = i, nvar
            l = l + 1
            kk = (k*(k - 1))/2 + i
            dumy(l) = (8.D0*(grold(l)-grad(l))-(g2old(l)-g2rad(l)))/delta*fact/12.D0
            eigs(l) = (2.D0*(grold(l)-grad(l))-(g2old(l)-g2rad(l)))/delta**3*fact/28.D0
!
!  CORRECT FOR 4'TH ORDER CONTAMINATION
!
            fmatrx(kk) = fmatrx(kk) + dumy(l)
          end do
        else
          if (lu - ll + 1 > 0) then
            dumy(l+1:lu-ll+1+l) = (grold(l+1:lu-ll+1+l)-grad(l+1:lu-ll+1+l))*0.5D0/delta*fact
            fmatrx(ll:lu) = fmatrx(ll:lu) + dumy(l+1:lu-ll+1+l)
            l = lu - ll + 1 + l
          end if
          l = l - 1
          do k = i, nvar
            l = l + 1
            kk = (k*(k - 1))/2 + i
            dumy(l) = (grold(l)-grad(l))*0.5D0/delta*fact
            fmatrx(kk) = fmatrx(kk) + dumy(l)
          end do
        end if
        if (big) then
          j = (i - 1)/3 + 1
          l = i - j*3 + 3
          write (iw, '(a, i'//num//', a)') ' CONTRIBUTIONS TO F-MATRIX FOR '//xyz(l)//" COMPONENT OF ATOM", &
            j, " ("//elemnt(labels(j))//")"
          if (precis) then
            write (iw, '(A)') &
      ' ELEMENT  +1.0*DELTA  +0.5*DELTA  -0.5*DELTA  -1.0*DELTA   2''ND ORDER 4TH ORDER'
            write (iw, '(I7,6F12.6)') (l,g2old(l),grold(l),grad(l),g2rad(l),dumy(l),eigs(l),l=1,nvar)
          else
            write (iw, '(A)') '  MATRIX ELEMENT  +Delta      -Delta      Hessian'
            do l = 1, nvar
              j = (l - 1)/3 + 1
              k = l - j*3 + 3
              write (iw, '(I7, a2, a3, i'//num//', 3F12.6)') l, xyz(k), elemnt(labels(j)), j, grold(l) ,grad(l), dumy(l)
            end do
            write (iw, '(a, 2f12.6)')"   SUM:", sum(grold(1:nvar)), sum(grad(1:nvar))
          end if
        end if
        time3 = seconds(2)
        tstep = time3 - time2
        tleft = max(0.1D0,tleft - tstep)
        if (tstep > 1.D7) tstep = tstep - 1.D7
        totime = totime + tstep
        if (resfil) then
          write (line, &
      '('' STEP:'',I4,'' RESTART FILE WRITTEN, INTEGRAL ='',F10.2,'' TIME LEFT:'',F10.2)') &
      i, min(totime,9999999.9D0), tleft
          write(iw,'(a)')trim(line)
          endfile (iw)
          backspace (iw)
          if (log) write (ilog, '(a)')trim(line)
          resfil = .FALSE.
        else
          write (line, &
      '('' STEP:'',I4,'' TIME ='',F9.2,'' SECS, INTEGRAL ='',F10.2,'' TIME LEFT:'',F10.2)') &
      i, tstep, totime, tleft
          write(iw,'(a)')trim(line)
          endfile (iw)
          backspace (iw)
          if (log) write (ilog, '(a)')trim(line)
        end if
        call to_screen(line)
        estim = totime/i
        if (tlast - tleft > tdump) then
          tlast = tleft
          resfil = .TRUE.
          jstart = 1
          ii = i
          call forsav (totime, deldip, ii, fmatrx, coord, nvar, heat, evecs, &
            jstart, fconst)
          if (moperr) return
        end if
        if (.not.(i /= nvar .and. tleft - 10.D0 < estim .or. i - istart >= maxcyc - 1)) &
          cycle
        write (iw, &
      '(2/10X,''- - - - - - - TIME UP - - - - - - -'',2/)')
        write (iw, '(/10X,'' POINT REACHED ='',I4)') i
        write (iw, '(/10X,'' RESTART USING KEY-WORD "RESTART"'')')
        write (iw, &
          '(10X,''ESTIMATED TIME FOR THE NEXT STEP ='',F9.2, '' SECONDS'')') estim
        jstart = 1
        ii = i
        call forsav (totime, deldip, ii, fmatrx, coord, nvar, heat, evecs, &
          jstart, fconst)
        if (moperr) return
        write (iw, '(2/10X,''FORCE MATRIX WRITTEN TO DISK'')')
        nreal = -1
        return
      end do
      do i = 1, natoms
        if (atmass(i)>=1.D-20 .or. labels(i)>=99) cycle
        call forsav (totime, deldip, nvar, fmatrx, coord, nvar, heat, evecs, &
          iloop, fconst)
        if (moperr) return
        write (iw, '(A)') ' AT LEAST ONE ATOM HAS A ZERO MASS. A RESTART'
        write (iw, '(A)') ' FILE HAS BEEN WRITTEN AND THE JOB STOPPED'
        call mopend (&
       'AT LEAST ONE ATOM HAS A ZERO MASS.  A RESTART FILE HAS BEEN WRITTEN AND &
       &THE JOB STOPPED')
        return
      end do
      if (istart<=nvar .and. index(keywrd,'ISOT')/=0) &
        call forsav (totime, deldip, nvar, fmatrx, coord, nvar, heat, evecs, iloop, fconst)
      if (moperr) return
      return
      end subroutine fmat
