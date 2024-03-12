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

subroutine ef (xparam, funct)
    use Common_arrays_C, only: geo, loc, hesinv, grad, atmass, nc
    use molkst_C, only: nvar, numcal, last, gnorm, iflepo, line, &
       & tleft, numat, ndep, time0, tdump, natoms, id, keywrd, moperr
    use chanel_C, only: iw0, iw, ilog, log, input_fn
    use ef_C, only: nstep, negreq, iprnt, ef_mode, ddx, xlamd, &
       & xlamd0, skal, rmin, rmax
    use maps_C, only : latom
    implicit none
    double precision, dimension (nvar), intent (inout) :: xparam
    double precision, intent (inout) :: funct
   !
   !*********************************************************************
   !
   !  EF IS A QUASI NEWTON RAPHSON OPTIMIZATION ROUTINE BASED ON
   !     JACK SIMONS P-RFO ALGORITHM AS IMPLEMENTED BY JON BAKER
   !     (J.COMP.CHEM. 7, 385). STEP SCALING TO KEEP LENGTH WITHIN
   !     TRUST RADIUS IS TAKEN FROM CULOT ET AL. (THEO. CHIM. ACTA 82, 189)
   !     THE TRUST RADIUS CAN BE UPDATED DYNAMICALLY ACCORDING TO FLETCHER
   !     SAFEGUARDS ON VALID STEP FOR TS SEARCHES BASED ON ACTUAL/PREDICTED
   !     FUNCTION CHANGE AND CHANGE IN TS MODE ARE OWN MODIFICATIONS
   !
   !  ON ENTRY XPARAM = VALUES OF PARAMETERS TO BE OPTIMISED.
   !           NVAR   = NUMBER OF PARAMETERS TO BE OPTIMISED.
   !
   !  ON EXIT  XPARAM = OPTIMISED PARAMETERS.
   !           FUNCT  = HEAT OF FORMATION IN KCAL/MOL.
   !
   !  CURRENT VERSION IMPLEMENTING COMBINED NR, P-RFO AND QA ALGORITHM
   !      TOGETHER WITH THRUST RADIUS UPDATE AND STEP REJECTION WAS
   !      MADE OCTOBER 1992 BY F.JENSEN, ODENSE, DK
   !
   !*********************************************************************
    character :: txt
    logical :: lorjk, lrjk, lts, old, l_geo_ok
    logical, save :: let, lupd, scf1, newhes, saddle, rrscal, first_time, donr
    integer, save :: icalcn = 0, igthes = 0, ihess, iloop, ireclc, &
         & iupd = 0, ldump, mxstep, ntime
!
    integer :: i, ij, imode, instep, itry1, ittest, j, k, l, neg, &
         & nflush, ih
    integer, dimension (9) :: ipow
    integer, dimension (:), allocatable :: best_nc
!
    double precision, parameter :: demin = 1.0d-2, gmin = 5.0d0, &
         & one = 1.d0, two = 2.d0,  three = 3.0d0, four = 4.d0, &
         & zero = 0.d0, pt5 = 0.5d0, pt75 = 0.75d0, tmone = 1.0d-1, &
         & tmtwo = 1.0d-2, tmsix = 1.0d-06
    double precision :: ss, absmin, deact, depre, dtmp, &
         & olde, ratio, tprt, tstep, tt0, xtmp, rmx, best_funct, &
         & best_gnorm, cosine, cos_const
    double precision, save :: time1, time2, tol2, dmax, ddmin, ddmax, &
           t_hess1, t_hess2
    double precision, save :: odmax = 0.d0, oolde = 0.d0,  odd = 0.d0
    double precision, dimension(:), allocatable :: pmat, u, hessc, eigval, &
         & tvec, svec, fx, oldfx, oldeig, ooldf, oldf, d, vmode
    double precision, dimension(:), allocatable :: x, rm, dx, best_xparam, &
           best_grad
    double precision, dimension(:,:), allocatable :: p, coord
    double precision, dimension(:,:), allocatable, save :: bmat
    double precision, external :: ddot, seconds
!

    intrinsic Abs, Index, Int, Max, Min, Mod, Sqrt
    i = 0
   !
   !   Create temporary arrays
   !
    ih = (nvar*(nvar+1)) / 2
    if (nvar < 7996) then
      allocate (pmat(nvar**2), bmat(nvar, nvar), u(nvar**2), hessc(ih), &
         & eigval(nvar), tvec(nvar), svec(nvar), fx(nvar), oldfx(nvar), &
         & oldeig(nvar), ooldf(nvar), oldf(nvar), d(nvar), vmode(nvar), &
         & best_xparam(nvar), best_grad(nvar), best_nc(natoms), stat=i)
    end if
    if (nvar > 7995 .or. i /= 0) then
      if (index(keywrd, " TS") /= 0) then
        call mopend("Insufficient memory to run TS")
      else
        call mopend("Insufficient memory to run EF")
        write(iw,'(6x,a)') " ", &
        "The L-BFGS geometry optimizer uses less memory", &
        "To use the L-BFGS geometry optimizer, add LBFGS to the keyword line", " "
      end if
      return
    end if
    tstep = 0.d0
    pmat(:) = 0.d0
    bmat(:,:) = 0.d0
    u(:) = 0.d0
    hessc(:) = 0.d0
    eigval(:) = 0.d0
    tvec(:) = 0.d0
    svec(:) = 0.d0
    fx(:) = 0.d0
    oldfx(:) = 0.d0
    oldeig(:) = 0.d0
    ooldf(:) = 0.d0
    oldf(:) = 0.d0
    d(:) = 0.d0
    vmode(:) = 0.d0
    best_xparam(:) = 0.d0
    best_grad(:) = 0.d0
    best_nc(:) = 0
    best_funct = 1.d10
    best_gnorm = 1.d10
   !
   !     GET ALL INITIALIZATION DATA
   !
    cos_const = 1.d0
    ihess = 0
    nstep = 0
    ipow(:) = 0
    if (allocated(hesinv)) then
      i = size(hesinv)
    else
      i = 0
    end if
    l_geo_ok = (index(keywrd, " GEO-OK") /= 0)
    if (icalcn /= numcal .or. nvar*nvar /= i) then
      newhes = .true.
      old = (index(keywrd, " OLD_HESS") /= 0)
      if (old) then
        old = allocated(hesinv)
        if ( .not. old) then
          write(line,"(a)") " OLD_HESS requested, but old Hessian is missing"
          write(iw,'(//10x,a,//)')trim(line)
          call mopend(trim(line))
        end if
        if (old) then
          old = (nvar*nvar == size(hesinv))
          if ( .not. old) then
            write(line,"(a)") " OLD_HESS requested, but old Hessian is of a different size"
            write(iw,'(//10x,a,//)')trim(line)
            call mopend(trim(line))
          end if
        end if
        if ( .not. old) then
          if (allocated(bmat)) deallocate(bmat)
          return
        end if
      end if
      if (.not. old) then
        if (allocated(hesinv)) deallocate(hesinv)
        allocate (hesinv(nvar*nvar))
        hesinv = 0.0D00
      end if
      call efstr (xparam, funct, ihess, ntime, iloop, igthes, mxstep, &
           & ireclc, iupd, dmax, ddmax, ddmin, tol2, time1, time2, nvar, &
           & scf1, lupd, ldump, rrscal, donr, hesinv, pmat, bmat, grad, &
           & oldf, d, vmode)
      if (moperr) go to 1100
      if (old) iloop = -1
      let = (Index (keywrd, " LET") /= 0)
      saddle = (Index (keywrd, " SADDLE") /= 0)
      first_time = .true.
    else
      first_time = (latom /= 0)
    end if
   !
   !  Assume that this is the last SCF, so that if geometry is optimized
   !  then the results will be OK.
   !
    last = 1
    do i = 1, nvar
      grad(i) = 0.d0
    end do
    time1 = seconds(1)
    call compfg (xparam, .true., funct, .true., grad, .true.)
    if (moperr) then
      go to 1100
    end if
    olde = funct
    nflush = 1
    absmin = 1.d9
    itry1 = 0
    lts = .false.
    if (negreq == 1) then
      lts = .true.
    end if
    lorjk = .false.
    if (scf1) then
      gnorm = dSqrt (ddot(nvar, grad, 1,grad, 1))
      iflepo = 1
      go to 1100
    end if
!     CHECK THAT GEOMETRY IS NOT ALREADY OPTIMIZED
    rmx = dSqrt (ddot(nvar, grad, 1, grad, 1))
!
    if (rmx < tol2) then
      iflepo = 2
      last = 1
      icalcn = numcal
      write (iw, '(/, 5 x, "GRADIENT =", f9.5, " IS LESS THAN CUTOFF =", f9.5,//)') rmx, tol2
      go to 1030
    end if
   !
   !  geometry is not optimized, so set LAST non-zero
   !
    last = 0
   !     GET INITIAL HESSIAN. IF ILOOP IS .LE.0 THIS IS AN OPTIMIZATION
   !     RESTART AND HESSIAN SHOULD ALREADY BE AVAILABLE
    t_hess1 = seconds(1)
    if (newhes .and. iloop > 0) then
      call gethes (xparam, igthes, iloop, hesinv, pmat, bmat, grad, &
           & geo, loc, oldf, d, vmode, funct)
      if (moperr) then
        go to 1100
      end if
      newhes = .false.
    end if
    t_hess2 = seconds(1)
    icalcn = numcal
   !     START OF MAIN LOOP
   !     WE NOW HAVE THE GRADIENT AND A HESSIAN. IF THIS IS THE FIRST
   !     TIME THROUGH DON'T UPDATE THE HESSIAN. FOR LATER LOOPS ALSO
   !     CHECK IF WE NEED TO RECALCULATE THE HESSIAN
    iflepo = 0
    funct = olde
    do
      !     STORE VARIOUS THINGS FOR POSSIBLY OMIN REJECTION
      do i = 1, nvar
        oldfx(i) = fx(i)
        ooldf(i) = oldf(i)
        oldeig(i) = eigval(i)
        do j = 1, nvar
          bmat(i, j) = hesinv(i + (j-1)*nvar)
          pmat(i + (j-1)*nvar) = u(i + (j-1)*nvar)
        end do
      end do
      if (ihess >= ireclc .and. iflepo /= 15) then
        iloop = 1
        ihess = 0
        if (igthes /= 3) then
          igthes = 1
        end if
        call gethes (xparam, igthes, iloop, hesinv, pmat, bmat, grad, &
             & geo, loc, oldf, d, vmode, funct)
        if (moperr) then
          go to 1100
        end if
      end if
      if (ihess > 0) then
        call updhes (svec, tvec, grad, nvar, iupd, hesinv, oldf, d)
      end if
      if (iprnt >= 2) then
        call geout (iw)
      end if
      if (iprnt >= 2) then
        write (iw, "(' XPARAM ')")
        write (iw, "(5(2I3,F10.4))") (loc(1, i), loc(2, i), xparam(i), i=1, &
       & nvar)
        write (iw, "(' grad')")
        write (iw, "(3X,8F11.5)") (grad(i), i=1, nvar)
      end if
      !
      !        PRINT RESULTS IN CYCLE
      gnorm = dSqrt (ddot(nvar, grad, 1, grad, 1))
!
      time2 = seconds (2)
      tstep = time2 - time1
      if (tstep < zero) then
        tstep = zero
      end if
      tleft = tleft - tstep
      if (tleft < 0.0d0) then
        tleft = -0.1d0
      end if
      time1 = time2
!
!  The following construction is designed to delete the time required for building the Hessian
!  from the time of the cycle.  This prevents the time for building the Hessian from distorting
!  the cycle time.
!
      if (tleft < (tstep - (t_hess2 - t_hess1))*two) go to 1030
      call prttim (tleft, tprt, txt)
      if ( .not. saddle .and. tstep > 0.5d0 .or. first_time) then
        if (ldump == 0) then
          if (id == 3) then
            call write_cell(iw)
            call write_cell(iw0)
          end if
          write (line, '(" CYCLE:", i6, " TIME:", f8.3, " TIME LEFT:", &
                & f6.2, a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
                nstep + 1, Min (tstep, 9999.99d0), tprt, txt, &
                & Min (gnorm, 999999.999d0), funct
          write(iw,"(a)")line(:len_trim(line))
          if (mod(nstep + 1,30) == 0) then
            line = trim(input_fn)
            call add_path(line)
            i = len_trim(line) - 5
          call to_screen(line(:i))
        end if
          endfile (iw)
          backspace (iw)
          if (log) write (ilog, "(a)")line(:len_trim(line))
          call to_screen(line)
        end if
        if (nflush /= 0) then
          if (Mod(nstep+1, nflush) == 0) then
              endfile (iw)
              backspace (iw)
            if (log) then
              endfile (ilog)
              backspace (ilog)
            end if
          end if
        end if
      else
        if (id == 3) then
          call write_cell(iw)
          call write_cell(iw0)
        end if
        write (line, '(" RESTART FILE WRITTEN,      TIME LEFT:", f6.2, &
        & a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
        tprt, txt, Min (gnorm, 999999.999d0), funct
        write(iw,"(a)")line(:len_trim(line))
        call to_screen(line)
        endfile (iw)
        backspace (iw)
        if (log) then
          write (ilog, '(" RESTART FILE WRITTEN,      TIME LEFT:", f6.2, &
        & a1, "  GRAD.:", f10.3, " HEAT:", g14.7)', err=1000) &
        tprt, txt, Min (gnorm,999999.999d0), funct
        end if
        if (nflush /= 0) then
          if (Mod(nstep+1, nflush) == 0) then
              endfile (iw)
              backspace (iw)
            if (log) then
                endfile (ilog)
                backspace (ilog)
            end if
          end if
        end if
        if (mod(nstep + 1,30) == 0) then
          line = trim(input_fn)
          call add_path(line)
          i = len_trim(line) - 5
          call to_screen(line(:i))
        end if
      end if
      call to_screen("To_file: Geometry optimizing")
!
!  Store best result up to the present
!
      if (ef_mode /= 0 .and. best_gnorm > gnorm .or. &
        & ef_mode == 0 .and. best_funct > funct) then
        best_gnorm = gnorm
        best_funct = funct
        best_xparam = xparam
        best_grad = grad
        best_nc(:natoms) = nc(:natoms)
      end if
      !
1000  ihess = ihess + 1
      nstep = nstep + 1
      !
      !        TEST FOR CONVERGENCE
      !
      rmx = dSqrt (ddot(nvar, grad, 1, grad, 1))
!
      if (rmx < tol2) then
         !
         !     ****** OPTIMIZATION TERMINATION ******
         !
        write (iw, '(/, 5 x, "GRADIENT =", f9.5, " IS LESS THAN CUTOFF =", f9.5,//)') rmx, tol2
        go to 1020
      else
        if (ef_mode == 0 .and. nstep > 5 .and. .not. l_geo_ok) then
          if (absmin-funct < 1.d-5) then
            if (itry1 > 15 .or. (rmx < 0.1d0 .and. itry1 > 9)) then
              write (iw, &
                 &"(//,' HEAT OF FORMATION IS ESSENTIALLY STATIONARY')")
              go to 1020
            end if
            itry1 = itry1 + 1
          else
            itry1 = 0
            absmin = funct
          end if
        end if
        olde = funct
        if (nstep > 5) then
          cosine = ddot(nvar, grad, 1, oldf, 1)/ &
                 & dsqrt(ddot(nvar, grad, 1, grad, 1)*ddot(nvar, oldf, 1, oldf, 1))
          if (cosine > 0.9d0) then
            cos_const = cos_const*1.5d0
          else
            cos_const = 1.d0
          end if
  !        write(iw,*)"COSINE:", cosine, cos_const
        end if
        do i = 1, nvar
          oldf(i) = grad(i)
        end do
         !
         ! IF THE OPTIMIZATION IS IN CARTESIAN COORDINATES, WE SHOULD REMOVE
         ! TRANSLATION AND ROTATION MODES. POSSIBLE PROBLEM IF RUN IS IN
         ! INTERNAL BUT WITH EXACTLY 3*NATOMS VARIABLE (I.E. DUMMY ATOMS
         ! ARE ALSO OPTIMIZED).
        if (nvar == 3*numat .and. numat > 2) then ! Sometime, add check that system is not linear.
                                                  ! Here NUMAT = 2 is a primitive check
          allocate (p(nvar,nvar), x(nvar), rm(nvar), dx(nvar), coord(3,nvar))

          p(:,:) = 0.d0
          x(:) = 0.d0
          rm(:) = 0.d0
          dx(:) = 0.d0
          coord(:,:) = 0.d0
          call prjfc (hesinv, xparam, nvar, u, p, atmass, x, rm, dx, coord)
          deallocate (p, x, rm, dx, coord)

          if (moperr) then
            go to 1100
          end if
        end if
        ij = 0
        do i = 1, nvar
          do j = 1, i
            ij = ij + 1
            hessc(ij) = hesinv(j + (i-1)*nvar)
          end do
        end do
        call rsp (hessc, nvar, eigval, u)
!
        do i = 1, nvar
          if (Abs (eigval(i)) < tmsix) then
            eigval(i) = zero
          end if
        end do
        ij = nvar ** 2
        do i = nvar, 1, -1
          do j = nvar, 1, -1
            u(j + (i-1)*nvar) = u(ij)
            ij = ij - 1
          end do
        end do
        if (iprnt >= 3) then
          call prthes (eigval, nvar, hesinv, u)
        end if
        if (mxstep == 0) then
          nstep = 0
          go to 1030
        end if
        neg = 0
        do i = 1, nvar
          if (eigval(i) < zero) then
            neg = neg + 1
          end if
        end do
        if (iprnt >= 1) then
10030     format (/, 10 x, "HESSIAN HAS", i3, " NEGATIVE EIGENVALUE(S)", &
               & 6 f7.1, /)
          write (iw, 10030) neg, (eigval(i), i=1, neg)
        end if
         ! IF AN EIGENVALUE HAS BEEN ZERO OUT IT IS PROBABLY ONE OF THE T,R
         ! MODES IN A CARTESIAN OPTIMIZATION. ZERO CORRESPONDING FX TO
         ! ALLOW FORMATION OF STEP WITHOUT THESE CONTRIBUTIONS. A SAFER
         ! CRITERION FOR DECIDING WHETHER THIS ACTUALLY IS A CARTESIAN
         ! OPTIMIZATION SHOULD BE PUT IN SOME DAY...
        do i = 1, nvar
          fx(i) = ddot (nvar, u(1 + (i-1)*nvar), 1, grad, 1)
          if (Abs (eigval(i)) == zero) then
            fx(i) = zero
          end if
        end do
        do
            !     FORM GEOMETRY STEP D
          call formd (eigval, fx, nvar, dmax, ddmin, lts, lorjk, rrscal, &
               & donr, u, d, vmode)
          if (moperr) then
            go to 1099
          end if
            ! IF LORJK IS TRUE, THEN TS MODE OVERLAP IS LESS THAN OMIN,
            ! REJECT PREVIOUS STEP
          if (lorjk) then
            if (iprnt >= 1) then
              write (iw,*) "      NOW UNDOING PREVIOUS STEP"
            end if
            dmax = odmax
            ddx = odd
            olde = oolde
            do i = 1, nvar
              fx(i) = oldfx(i)
              oldf(i) = ooldf(i)
              eigval(i) = oldeig(i)
              do j = 1, nvar
                hesinv(i + (j-1)*nvar) = bmat(i, j)
                u(i + (j-1)*nvar) = pmat(i + (j-1)*nvar)
              end do
            end do
            do i = 1, nvar
              xparam(i) = xparam(i) - d(i)
              k = loc(1, i)
              l = loc(2, i)
              geo(l, k) = xparam(i)
            end do
            if (ndep /= 0) then
              call symtry ()
            end if
            dmax = Min (dmax, ddx) / two
            odmax = dmax
            odd = ddx
            nstep = nstep - 1
            if (dmax < ddmin) then
              go to 1010
            else if (iprnt >= 1) then
              write (iw,*) "      FINISH UNDOING, NOW GOING FOR NEW STEP"
            end if
          else
               !
               !  FORM NEW TRIAL XPARAM AND GEO
               !
            d(:nvar) = d(:nvar)*cos_const
            do i = 1, nvar
              xparam(i) = xparam(i) + d(i)
              k = loc(1, i)
              l = loc(2, i)
              geo(l, k) = xparam(i)
            end do
            if (ndep /= 0) then
              call symtry ()
            end if
               !
               !     COMPARE PREDICTED E-CHANGE WITH ACTUAL
               !
            depre = zero
            imode = 1
            if (ef_mode /= 0) then
              imode = ef_mode
            end if
            do i = 1, nvar
              if (lts .and. i == imode) then
                xtmp = xlamd0
              else
                xtmp = xlamd
              end if
              if (Abs (xtmp-eigval(i)) >= tmtwo) then
                ss = skal * fx(i) / (xtmp-eigval(i))
                depre = depre + ss * (fx(i) + pt5 * ss * eigval(i))
              end if
            end do
               !
               !     GET GRADIENT FOR NEW GEOMETRY
               !
            grad(:nvar) = 0.d0
            call compfg (xparam, .true., funct, .true., grad, .true.)
            if (moperr) then
              go to 1100
            end if
            deact = funct - olde
            if (depre == zero) then
              write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
              call mopend ("in EF")
              go to 1100
            end if
            ratio = deact / depre
         !   if (nstep > 1 .and. Abs(depre) < 1.d0) ratio = 1.d0 ! ignore ratio when step is small
            if (iprnt >= 1) then
		!                   12345678901234567890123456789012345678901234567890
               write (iw, '("       HoF         ACTUAL,  PREDICTED ENERGY CHANGE, RATIO",/ &
                   & 2 f14.7, f20.7, f13.7)') funct, deact, depre, ratio
            end if
               !
               ! POSSIBLY REJECT THE STEP IF THE RATIO BETWEEN ACTUAL AND
               ! PREDICTED CHANGE IN ENERGY IS OUTSIDE RMIN AND RMAX LIMITS
               ! THE DEFAULT VALUES OF RMIN=0.0 FOR MINIMIZATIONS IS
               ! EQUIVALENT
               ! TO NOT ALLOWING THE ENERGY THE RAISE.
               ! DON'T WORRY IS THE ABSOLUTE CHANGES ARE SMALL ( < DEMIN)
               !
            lrjk = .false.
            if ((ratio < rmin .or. ratio > rmax) .and. &
                 & (Abs (depre) > demin .or. Abs (deact) > demin)) then
              dtmp = Min (dmax, ddx) / two
              if (dtmp <= ddmin) then
                dtmp = ddmin
              end if
              if (iprnt >= 1) then
10050           format (1 x, "UNACCEPTABLE RATIO,", &
                     & " REJECTING STEP, REDUCING DMAX TO", f7.4)
                write (iw, 10050) dtmp
              end if
              lrjk = .true.
            end if
               !     IF THE TRUST RADIUS IS EQUAL TO DDMIN, CONTINUE ANYWAY
            if (lrjk .and. Abs(dmax - ddmin) < 1.d-20) then
              if (iprnt >= 1) then
10060           format (1 x, "NEW TRUST RADIUS WOULD BE BELOW DDMIN", &
                     & " ACCEPTING STEP ANYWAY")
                write (iw, 10060)
              end if
              lrjk = .false.
            end if
               !
               !   BYPASS ALL TESTS
               !
            if (let) then
              lrjk = .false.
            end if
               !
            if ( .not. lrjk) exit
            do i = 1, nvar
              xparam(i) = xparam(i) - d(i)
              k = loc(1, i)
              l = loc(2, i)
              geo(l, k) = xparam(i)
            end do
            if (ndep /= 0) then
              call symtry ()
            end if
            dmax = Min (dmax, ddx) / two
            if (dmax < ddmin) then
              go to 1010
            end if
          end if
        end do
        if (iprnt >= 1) then
10070     format (5 x, "STEPSIZE USED IS", f9.5)
          write (iw, 10070) ddx
        end if
        if (iprnt >= 2) then
          write (iw, "(' CALCULATED STEP')")
          write (iw, "(3X,8F10.6)") (d(i), i=1, nvar)
        end if
         !
         !     POSSIBLE USE DYNAMICAL TRUST RADIUS
        odmax = dmax
        odd = ddx
        oolde = olde
        if (lupd .and. ((rmx > gmin) .or. (Abs (depre) > demin .or. &
             & Abs (deact) > demin))) then
            ! FLETCHER RECOMMEND DMAX=DMAX/4 AND DMAX=DMAX*2
            ! THESE ARE A LITTLE MORE CONSERVATIVE SINCE HESSIAN IS BEING
            ! UPDATED
            ! DON'T REDUCE TRUST RADIUS DUE TO RATIO FOR MIN SEARCHES
          if (lts .and. ratio <= tmone .or. ratio >= three) then
            dmax = Min (dmax, ddx) / two
          end if
          if (lts .and. ratio >= pt75 .and. ratio <= (four/three) .and. ddx > &
         & (dmax-tmsix)) then
            dmax = dmax * Sqrt (two)
          end if
            ! ALLOW WIDER LIMITS FOR INCREASING TRUST RADIUS FOR MIN SEARCHES
          if ( .not. lts .and. ratio >= pt5 .and. ddx >(dmax-tmsix)) then
            dmax = dmax * Sqrt (two)
          end if
            ! BE BRAVE IF  0.90 < RATIO < 1.10 ...
          if (Abs (ratio-one) < tmone) then
            dmax = dmax * Sqrt (two)
          end if
          dmax = Max (dmax, ddmin)
          dmax = Min (dmax, ddmax)
        end if
         ! ALLOW STEPSIZE UP TO 0.1 IN THE END-GAME WHERE CHANGES ARE LESS
         ! THAN DEMIN AND GRADIENT IS LESS THAN GMIN
        if (lupd .and. rmx < gmin .and. (Abs (depre) < demin .and. &
             & Abs (deact) < demin)) then
          dmax = Max (dmax, tmone)
        end if
        if (iprnt >= 1) then
10080     format (5 x, "NEW TRUST RADIUS = ", f8.5)
          write (iw, 10080) dmax
        end if
         !frj  this test should never be encountered with the current update...
1010    if (dmax < ddmin) exit
         !     CHECK STEPS AND ENOUGH TIME FOR ANOTHER PASS
        if (nstep >= mxstep) go to 1030
         !     IN USER UNFRIENDLY ENVIROMENT, SAVE RESULTS EVERY 1 CPU HRS
        ittest = Int ((time2-time0)/tdump)
         !     RETURN FOR ANOTHER CYCLE
        if (ittest > ntime) then
          ldump = 1
          ntime = Max (ittest, (ntime+1))
          ipow(1) = ihess
          ipow(2) = nstep
          ipow(9) = 2
          tt0 = seconds (1) - time0
          instep = -nstep
          call efsav (tt0, hesinv, funct, grad, xparam, pmat, instep, &
               bmat, ipow, oldf, d, vmode)

          if (moperr) then
            go to 1100
          end if
        else
          ldump = 0
        end if
      end if
    end do
10090 format (/, 5 x, "TRUST RADIUS NOW LESS THAN ", f8.5, &
           & " OPTIMIZATION TERMINATING",/, 5 x, &
           & " THE GEOMETRY MAY NOT BE COMPLETELY OPTIMIZED",/, 5 x, &
           & " (TO CONTINUE, ADD 'LET DDMIN=0.0' TO THE KEYWORD LINE)")
    write (iw, 10090) ddmin
1099 continue
!
!  If current point is not the best, then load in the best point
!
    if (ef_mode /= 0 .and. best_gnorm < gnorm .or. ef_mode == 0 .and. best_funct < funct) then
          funct = best_funct
          gnorm = best_gnorm
          xparam(:nvar) = best_xparam(:nvar)
          grad(:nvar) = best_grad(:nvar)
          nc(:natoms) = best_nc(:natoms)
      end if
1020 iflepo = 15
    last = 1
    tt0 = seconds (1) - time0
    if ( .not. moperr) then
   !     CALL COMPFG TO CALCULATE ENERGY FOR FIXING MO-VECTOR BUG
      call compfg (xparam, .true., funct, .true., grad, .false.)
    end if
    go to 1100
   !     WE RAN OUT OF TIME OR TOO MANY ITERATIONS. DUMP RESULTS
1030 if (tleft < tstep*two) call mopend("NOT ENOUGH TIME FOR ANOTHER CYCLE")
    if (nstep >= mxstep) call mopend("EXCESS NUMBER OF OPTIMIZATION CYCLES")
    if (tleft < tstep*two .or. nstep >= mxstep) then
      ipow(1) = ihess
      ipow(9) = 1
      ipow(2) = nstep
      tt0 = seconds (1) - time0
      instep = -nstep
      call efsav (tt0, hesinv, funct, grad, xparam, pmat, instep, &
                & bmat, ipow, oldf, d, vmode)
      iflepo = -1
    end if
   !
1100 continue

   !   Delete temporary arrays
   !
    deallocate (pmat, bmat, u, hessc, eigval, tvec, svec, fx, oldfx, &
         & oldeig, ooldf, oldf, d, vmode, best_xparam, best_grad, best_nc)
    return
end subroutine ef

subroutine efsav (tt0, hess, funct, grad, xparam, pmat, il, bmat, ipow, &
     & oldf, d, vmode)
  USE maps_C, only : rxn_coord, latom, kloop
  use molkst_C, only : nvar, keywrd, nscf, is_PARAM, norbs, numat, mozyme, prt_gradients
  use common_arrays_C, only : profil
  USE ef_C, ONLY: ddx, ef_mode, &
   nstep, negreq, alparm, x0, x1, x2, iloop
  USE chanel_C, only : iw, ires, restart_fn
  implicit none
  integer :: il
  double precision  :: tt0
  double precision  :: funct
  integer  :: ipow(9)
  double precision  :: hess(nvar,nvar)
  double precision  :: grad(nvar), oldf(*), d(nvar), vmode(nvar)
  double precision  :: xparam(nvar)
  double precision  :: pmat(*)
  double precision  :: bmat(nvar,nvar)
!
  integer ::  i, j, linear, io_stat, old_numat, old_norbs
  double precision :: funct1
  logical :: opend
  double precision, external :: ddot
!*********************************************************************
!
! EFSAV STORES AND RETRIEVE DATA USED IN THE EF GEOMETRY
!        OPTIMISATION. VERY SIMILAR TO POWSAV.
!
!  ON INPUT HESS   = HESSIAN MATRIX, PARTIAL OR WHOLE.
!           GRAD   = GRADIENTS.
!           XPARAM = CURRENT STATE OF PARAMETERS.
!           IL     = INDEX OF HESSIAN,
!           JL     = CYCLE NUMBER REACHED SO-FAR.
!           BMAT   = "B" MATRIX!
!           IPOW   = INDICES AND FLAGS.
!           IPOW(9)= 0 FOR RESTORE, 1 FOR DUMP, 2 FOR SILENT DUMP
!
!*********************************************************************
  if (is_PARAM) return
  inquire(unit=ires, opened=opend)
  if (opend) close(unit=ires, status='KEEP')
  open(unit=ires, file=restart_fn, status='UNKNOWN', form=&
  'UNFORMATTED', position='asis', iostat = io_stat)
  if (io_stat /= 0) then
    write(iw,*)" Restart file either does not exist or is not available for reading"
    call mopend ("Restart file either does not exist or is not available for reading")
    return
  end if
  rewind ires
  if (ipow(9)==1 .or. ipow(9)==2) then
    funct1 = dsqrt(ddot(nvar, grad, 1, grad, 1))
    if (ipow(9) == 1 .and. index(keywrd,'STEP1') == 0) then
      write (iw, '(2/10X,''CURRENT VALUE OF GRADIENT NORM ='',F12.6)') funct1
      if (prt_gradients .and. index(keywrd," GRADI") /= 0 .and. mozyme) then
        write (iw, '(3/7X,''CURRENT  POINT  AND  DERIVATIVES'',/)')
        call prtgra ()
      end if
      write (iw, '(/10X,''CURRENT VALUE OF GEOMETRY'',/)')
      call geout (iw)
    end if
!
!  IPOW(1) AND IPOW(9) ARE USED ALREADY, THE REST ARE FREE FOR USE
!
    ipow(8) = nscf
    write (ires) numat, norbs, (xparam(i),i=1,nvar)
    if (latom /= 0) then
      if (index(keywrd,' STEP=') /= 0) then
        write (ires) kloop
        write (ires) rxn_coord
        write (ires) (profil(i),i=1,kloop)
      else
        write (ires) ((alparm(j,i),j=1,3),i=1,nvar)
       write (ires) iloop, x0, x1, x2
      end if
    end if
    write (ires) ipow, il, nstep, funct, tt0
    write (ires) (grad(i),i=1,nvar)
    write (ires) ((hess(j,i),j=1,nvar),i=1,nvar)
    write (ires) ((bmat(j,i),j=1,nvar),i=1,nvar)
    write (ires) (oldf(i),i=1,nvar), (d(i),i=1,nvar), (vmode(i),i=1,nvar)
    write (ires) ddx, ef_mode, nstep, negreq
    linear = (nvar*(nvar + 1))/2
    write (ires) (pmat(i),i=1,linear)
    call den_in_out(1)
    if (index(keywrd,'STEP1') /= 0) return
    close(ires)
    return
  else
    read (ires, iostat = io_stat)old_numat, old_norbs
    if (norbs /= old_norbs .or. numat /= old_numat) then
        call mopend("Restart file read in does not match current data set")
        return
    end if
    if (latom /= 0) then
      if (index(keywrd,' STEP=') /= 0) then
        read (ires, iostat = io_stat) kloop
        read (ires, iostat = io_stat) rxn_coord
        read (ires, iostat = io_stat) (profil(i),i=1,kloop)
      else
        read (ires, iostat = io_stat) ((alparm(j,i),j=1,3),i=1,nvar)
        read (ires, iostat = io_stat) iloop, x0, x1, x2
      end if
    end if
    read (ires, iostat = io_stat) ipow, il, nstep, funct, tt0
    nscf = ipow(8)
    i = int(tt0/1000000)
    tt0 = tt0 - i*1000000
    write (iw, '(2/10X,''TOTAL TIME USED SO FAR:'',F13.2,'' SECONDS'')') tt0
    if(abs(funct) > 1.d-20) write (iw, '(  10X,''              FUNCTION:'',F17.6)') funct
    read (ires, iostat = io_stat) (grad(i),i=1,nvar)
    read (ires, iostat = io_stat) ((hess(j,i),j=1,nvar),i=1,nvar)
    read (ires, iostat = io_stat) ((bmat(j,i),j=1,nvar),i=1,nvar)
    read (ires, iostat = io_stat) (oldf(i),i=1,nvar), (d(i),i=1,nvar), (vmode(i),i=1,nvar)
    read (ires, iostat = io_stat) ddx, ef_mode, nstep, negreq
    linear = (nvar*(nvar + 1))/2
    read (ires, iostat = io_stat) (pmat(i),i=1,linear)
    if (index(keywrd,' STEP1') /= 0) return
    close(ires)
    if (io_stat /= 0) then
      call mopend ("Restart file is currupt")
    end if
    return
  end if
  return
  end subroutine efsav
subroutine efstr (xparam, funct, ihess, ntime, iloop, igthes, mxstep, ireclc, &
     & iupd, dmax, ddmax, ddmin, tol2, time1, time2, nvar, scf1, lupd, &
     & ldump, rrscal, donr, hess, bmat, pmat, grad, oldf, d, vmode)
    use molkst_C, only: limscf, time0, numcal, keywrd, moperr, id
    use chanel_C, only: iw
    use ef_C, only: nstep, ef_mode, negreq, rmin, rmax, omin, iprnt
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    logical, intent (out) :: donr, lupd, rrscal, scf1
    integer, intent (inout) :: nvar
    integer, intent (out) :: igthes, ihess, iloop, ireclc, iupd, ldump, &
         & mxstep, ntime
    double precision, intent (inout) :: funct
    double precision, intent (out) :: ddmax, ddmin, dmax, time1, time2, tol2
    double precision, dimension (nvar), intent (inout) :: d, grad, &
         & oldf, vmode, xparam
    double precision, dimension (nvar, nvar), intent (inout) :: &
         & bmat, hess, pmat
   !
   !.. Local Scalars ..
    logical :: restrt
    integer :: i, ip, its, j, k, mtmp, icalcn=0
    double precision :: tt0
   !
   !.. Local Arrays ..
    integer, dimension (9) :: ipow
   !
   !.. External Calls ..
    double precision, external :: reada
    intrinsic Abs, Index, Int, Nint
   !
   ! ... Executable Statements ...
   !
   !     GET ALL INITIALIZATION DATA
    nvar = Abs (nvar)
    ldump = 0
    lupd = (Index (keywrd, " NOUPD") == 0)
    restrt = (Index (keywrd, " RESTART") /= 0)
    scf1 = (Index (keywrd, " 1SCF") /= 0)
    nstep = 0
    ihess = 0
    ntime = 0
    iloop = 1
    ef_mode = 0
    igthes = 0
    iupd = 2
    negreq = 0
    rmin = 0.0d0
    rmax = 1.d3
    dmax = 0.2d0
    ddmax = 0.5d0
    if (icalcn == numcal) then
!
!  This is a reaction path or a grid, etc. therefore set maximum step small,
!  and re-set the Hessian every 10 cycles of geo optimization.
!
      ddmax = 0.1d0
      ireclc = 10
    else
      icalcn = numcal
      ireclc = 999999
    end if
   !
   !   In geometry optimization, let SCF exit if conditions are right,
   !   that is, if the heat has increased or decreased a lot.
   !
    limscf = .true.
    its = Index (keywrd, " TS ")
    if (its /= 0) then
      limscf = .false.
      ef_mode = 1
      igthes = 1
      iupd = 1
      negreq = 1
      rmin = 0.0d0
      rmax = 4.0d0
      omin = 0.8d0
      dmax = 0.1d0
      ddmax = 0.3d0
    end if
    rrscal = .false.
   !
    i = Index (keywrd, " RSCAL")
    if (i /= 0) then
      rrscal = .true.
    end if
    donr = .true.
    i = Index (keywrd, " NONR")
    if (i /= 0) then
      donr = .false.
    end if
    iprnt = 0
    ip = Index (keywrd, " PRNT=")
    if (ip /= 0) then
      iprnt = Nint (reada (keywrd, ip))
    end if
    if (iprnt > 5) then
      iprnt = 5
    end if
    if (iprnt < 0) then
      iprnt = 0
    end if
    mxstep = 2000
    i = Index (keywrd, " CYCLES=")
    if (i /= 0) then
      mxstep = Nint (reada (keywrd, i))
    end if
    if (i /= 0 .and. mxstep == 0 .and. ip == 0) then
      iprnt = 3
    end if
    i = Index (keywrd, " RECALC=")
    if (i /= 0) then
      ireclc = Nint (reada (keywrd, i))
    end if
    i = Index (keywrd, " IUPD=")
    if (i /= 0) then
      iupd = Nint (reada (keywrd, i))
    end if
    i = Index (keywrd, " MODE=")
    if (i /= 0) then
      ef_mode = Nint (reada (keywrd, i))
    end if
    ddmin = 1.d-4
    i = Index (keywrd, " DDMIN=")
    if (i /= 0) then
      ddmin = reada (keywrd, i)
    end if
    i = Index (keywrd, " DMAX=")
    if (i /= 0) then
      dmax = reada (keywrd, i)
    end if
    i = Index (keywrd, " DDMAX=")
    if (i /= 0) then
      ddmax = reada (keywrd, i)
    end if
    tol2 = 1.d+0
    if (id /= 0) tol2 = id*2.d0 - 1.d0
    if (Index (keywrd, " PREC") /= 0) tol2 = tol2*5.d-2
    i = Index (keywrd, " GNORM=")
    if (i /= 0) then
      tol2 = reada (keywrd, i)
    end if
    if (Index (keywrd, " LET") == 0 .and. tol2 < 0.01d0) then
      write (iw, "(/,A)") "  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01", &
     & "  SPECIFY LET AS KEYWORD TO ALLOW GNORM LESS THAN 0.01"
      tol2 = 0.01d0
    end if
    i = Index (keywrd, " HESS=")
    if (i /= 0) then
      igthes = Nint (reada (keywrd, i))
    end if
    i = Index (keywrd, " RMIN=")
    if (i /= 0) then
      rmin = reada (keywrd, i)
    end if
    i = Index (keywrd, " RMAX=")
    if (i /= 0) then
      rmax = reada (keywrd, i)
    end if
    i = Index (keywrd, " OMIN=")
    if (i /= 0) then
      omin = reada (keywrd, i)
    end if
    time1 = time0
    time2 = time1
   !   DONE WITH ALL INITIALIZING STUFF.
   !   CHECK THAT OPTIONS REQUESTED ARE RESONABLE
    if ((its /= 0) .and. (iupd == 2)) then
      write (iw,*) " TS SEARCH AND BFGS UPDATE WILL NOT WORK"
      call mopend ("TS SEARCH AND BFGS UPDATE WILL NOT WORK")
      return
    end if
    if ((its /= 0) .and. (igthes == 0)) then
      write (iw,*) " TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN"
      call mopend ("TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN")
      return
    end if
    if ((igthes < 0) .or. (igthes > 3)) then
      write (iw,*) " UNRECOGNIZED HESS OPTION", igthes
      call mopend ("UNRECOGNIZED HESS OPTION")
      return
    end if
    if ((omin < 0.d0) .or. (omin > 1.d0)) then
      write (iw,*) " OMIN MUST BE BETWEEN 0 AND 1", omin
      call mopend ("OMIN MUST BE BETWEEN 0 AND 1")
      return
    end if
    nstep = 0
    if (restrt) then
      !
      !   RESTORE DATA. I INDICATES (HESSIAN RESTART OR OPTIMIZATION
      !   RESTART). IF I .GT. 0 THEN HESSIAN RESTART AND I IS LAST
      !   STEP CALCULATED IN THE HESSIAN. IF I .LE. 0 THEN J (NSTEP)
      !   IN AN OPTIMIZATION HAS BEEN DONE.
      !
      ipow(9) = 0
      mtmp = ef_mode
      tt0 = 0.d0
      j = 0
      call efsav (tt0, hess, funct, grad, xparam, pmat, i, bmat, ipow, &
           & oldf, d, vmode)
      if (moperr) return
      ef_mode = mtmp
      j = ipow(2)
      k = Int (tt0/1000000.d0)
      time0 = time0 - tt0 + k * 1000000.d0
      iloop = i
      if (i > 0) then
        igthes = 4
        nstep = j
        write (iw, "(10X,'RESTARTING HESSIAN AT POINT',I4)") iloop
        if (nstep /= 0) then
          write (iw, "(10X,'IN OPTIMIZATION STEP',I4)") nstep
        end if
      else
        nstep = j
        write (iw, "(//10X,'RESTARTING OPTIMIZATION AT STEP',I6)") nstep
      end if
    end if
    mxstep = mxstep + nstep
end subroutine efstr
subroutine formd (eigval, fx, nvar, dmax, ddmin, ts, lorjk, rrscal, &
     & donr, u, d, vmode)
   !     THIS VERSION FORMS GEOMETRY STEP BY EITHER PURE NR, P-RFO OR QA
   !     ALGORITHM, UNDER THE CONDITION THAT THE STEPLENGTH IS LESS THAN
   !     DMAX
    use chanel_C, only: iw
    use ef_C, only: skal, ef_mode, iprnt, ddx, xlamd, xlamd0
    use molkst_C, only: numcal, numat
    implicit none
    logical, intent (in) :: donr, rrscal, ts
    logical, intent (inout) :: lorjk
    integer, intent (in) :: nvar
    double precision, intent (in) :: ddmin, dmax
    double precision, dimension (nvar), intent (in) :: eigval, fx
    double precision, dimension (nvar), intent (inout) :: d, vmode
    double precision, dimension (nvar, nvar), intent (in) :: u
    logical :: frodo1, frodo2, rscal
    integer :: i, it, j, jt, ncnt, newmod
    integer, save :: icalcn = 0
    double precision, parameter :: big = 1.0d+3, half = 0.5d0
    double precision, save :: eone, eigit
    double precision, parameter :: tmtwo = 1.0d-2
    double precision :: bl, bu, d2max, fl, fm, fu, ssmax, ssmin, &
         & sstoll, temp, xlambda, lambda, lambda0, sstep, store_ddx
    double precision, parameter :: eps = 1.0d-12, four = 4.0d+00, &
         & one = 1.0d+0, sfix = 1.0d+01, step_loc = 5.0d-02, ten = 1.0d+1, &
         & tmsix = 1.0d-06, toll = 1.0d-8, zero = 0.0d0
    save :: store_ddx

!For MOPAC BLAS
    double precision, external :: ddot
!

   !
   !
   !.. Intrinsic Functions ..
    intrinsic Abs, Max, Sqrt
   !
   ! ... Executable Statements ...
   !
    if (icalcn /= numcal) then
      icalcn = numcal
      store_ddx = 0.d0
      d = 0.d0
    end if
    skal = one
    rscal = rrscal
    it = 0
    jt = 1
    if (ts) then
      if (ef_mode /= 0) then
        call overlp (dmax, ddmin, newmod, nvar, lorjk, u, vmode)
        if (lorjk) return
         !
         !  ON RETURN FROM OVERLP, NEWMOD IS THE TS MODE
         !
        if (newmod /= ef_mode .and. iprnt >= 1) then
10000     format (5 x, "WARNING! MODE SWITCHING. WAS FOLLOWING MODE ", &
         & i3, " NOW FOLLOWING MODE ", i3)
          write (iw, 10000) ef_mode, newmod
        end if
        ef_mode = newmod
        it = ef_mode
      else
        it = 1
      end if
      eigit = eigval(it)
      if (eigit == zero) then
        if (nvar > 3*numat - 5) then
          call mopend ("Too many variables. By definition, at least one force constant is exactly zero" )
          write(iw,"(10x,a)")"and the lowest force constant is not negative."
          write(iw,"(10x,a,i5)")"Number of variables = ", nvar
          write(iw,"(10x,a,i5)")"Number of atoms     = ", numat
          if (nvar == 3*numat) &
            write(iw,"(10x,a,i5)")"(If coordinates are Cartesian, convert to internal coordinates and re-run.)"
        else
          write(iw,"(a)")" At least one force constant is exactly zero"
          call mopend("At least one force constant is exactly zero")
        end if
        return
      else if (iprnt >= 1) then
10010   format (/, 5 x, "TS MODE IS NUMBER", i3, " WITH EIGENVALUE", &
       & f9.1,/, 5 x, "AND COMPONENTS", /)
        write (iw, 10010) it, eigit
10020   format (5 x, 8 f9.4)
        write (iw, 10020) (u(i, it), i=1, nvar)
      end if
    end if
   !     JT SHOULD BE LOWEST MODE WHICH IS NOT THE TS-MODE AND NOT T, R MOD
   !     IN CARTESIAN COORDINATES. JT SHOULD BE ONE OF THE FIRST 8 MODES,
   !     SEARCH ALL TO ALLOW FOR ACCIDENTAL EIGENVALUES = 0.
    do i = 1, nvar
      if (i /= it .and. eigval(i) /= zero) then
        jt = i
        exit
      end if
    end do
    eone = eigval(jt)
    ssmin = Max (Abs (eone)*eps, (ten*eps))
    ssmax = Max (big, Abs (eone))
    ssmax = ssmax * big
    sstoll = toll
    d2max = dmax * dmax
   !  SOLVE ITERATIVELY FOR LAMBDA
   !  INITIAL GUESS FOR LAMBDA IS ZERO EXCEPT NOTE THAT
   !  LAMBDA SHOULD BE LESS THAN EIGVAL(1)
   !  START BY BRACKETING ROOT, THEN HUNT IT DOWN WITH BRUTE FORCE BISECT.
   !
    frodo1 = .false.
    frodo2 = .false.
    lambda = zero
    lambda0 = zero
    if (ts) then
      if (eigit < zero .and. eone >= zero .and. donr) then
        if (iprnt >= 1) then
          write (iw,*) " TS SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP"
        end if
        go to 1010
      end if
    else if (donr) then
      if (eone >= zero) then
        if (iprnt >= 1) then
          write (iw,*) " MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP"
        end if
        go to 1010
      end if
    end if
   !     .. Head of OUTER_LOOP ..
1000 continue
    if (ts) then
      lambda0 = eigval(it) + Sqrt (eigval(it)**2+four*fx(it)**2)
      lambda0 = lambda0 * half
      if (iprnt >= 1) then
         !
10030   format (1 x, "LAMBDA THAT MAXIMIZES ALONG TS MODES =   ", f15.5)
        write (iw, 10030) lambda0
      end if
    end if
    sstep = step_loc
    if (eone <= zero) then
      lambda = eone - sstep
    end if
    if (eone > zero) then
      sstep = eone
    end if
    bl = lambda - sstep
    bu = lambda + sstep * half
    ncnt = 0
    do
      fl = zero
      fu = zero
      do i = 1, nvar
        if (i /= it) then
          if ((bl-eigval(i) == zero) .or. &
            & (bu-eigval(i) == zero)) then
            write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
            call mopend ("in FORMD")
            return
          end if
          fl = fl + (fx(i)*fx(i)) / (bl-eigval(i))
          fu = fu + (fx(i)*fx(i)) / (bu-eigval(i))
        end if
      end do
      fl = fl - bl
      fu = fu - bu
      if (fl*fu < zero) exit
      bl = bl - (eone-bl)
      bu = bu + half * (eone-bu)
      if (bl <=-ssmax) then
        bl = -ssmax
        frodo1 = .true.
      end if
      if (Abs (eone-bu) <= ssmin) then
        bu = eone - ssmin
        frodo2 = .true.
      end if
      endfile (iw)
      backspace (iw)
      if (frodo1 .and. frodo2) then
        write (iw,*) "NUMERICAL PROBLEMS IN BRACKETING LAMBDA", eone, bl, &
             & bu, fl, fu
        write (iw,*) " GOING FOR FIXED STEP SIZE...."
        go to 1020
      end if
      ncnt = ncnt + 1
      if (ncnt > 1000) then
        write (iw,*) "TOO MANY ITERATIONS IN LAMBDA BISECT", bl, bu, &
             & fl, fu
        call mopend ("TOO MANY ITERATIONS IN LAMBDA BISECT IN EF")
        return
      end if
    end do
    ncnt = 0
    xlambda = zero
    do
      fl = zero
      fu = zero
      fm = zero
      lambda = half * (bl+bu)
      do i = 1, nvar
        if (i /= it) then
          if ((bl-eigval(i) == zero) .or. &
            & (bu-eigval(i) == zero) .or. &
            & (lambda-eigval(i) == zero)) then
            write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
            call mopend ("in FORMD")
            return
          end if
          fl = fl + (fx(i)*fx(i)) / (bl-eigval(i))
          fu = fu + (fx(i)*fx(i)) / (bu-eigval(i))
          fm = fm + (fx(i)*fx(i)) / (lambda-eigval(i))
        end if
      end do
      fl = fl - bl
      fu = fu - bu
      fm = fm - lambda
      if (Abs (xlambda-lambda) < sstoll) exit
      ncnt = ncnt + 1
      if (ncnt > 1000) then
        write (iw,*) "TOO MANY ITERATIONS IN LAMBDA BISECT", bl, bu, &
             & lambda, fl, fu
        call mopend ("TOO MANY ITERATIONS IN LAMBDA BISECT IN EF")
        return
      else
        xlambda = lambda
        if (fm*fu < zero) then
          bl = lambda
        end if
        if (fm*fl < zero) then
          bu = lambda
        end if
      end if
    end do
   !     .. Head of LOOP ..
1010 continue
   !
    if (iprnt >= 1) then
10040 format (1 x, "LAMBDA THAT MINIMIZES ALONG ALL MODES =  ", f17.9)
      write (iw, 10040) lambda
    end if
   !
   !  CALCULATE THE STEP
   !
    do i = 1, nvar
      d(i) = zero
    end do
    do i = 1, nvar
      if (lambda == zero .and. Abs (eigval(i)) < tmtwo) then
        temp = zero
      else
        temp = fx(i) / (lambda-eigval(i))
      end if
      if (i == it) then
        if (Abs (lambda0-eigval(it)) < 1.d-9) then
          write (iw,*) " TS FAILED TO LOCATE TRANSITION STATE"
          write (iw, "(/10X,'CURRENT VALUE OF GEOMETRY',/)")
          call geout (iw)
          call mopend ("TS FAILED TO LOCATE TRANSITION STATE")
          return
        else
          temp = fx(it) / (lambda0-eigval(it))
        end if
      end if
      if (iprnt >= 5) then
        write (iw, "(A,I4,F16.11)") "FORMD, DELTA STEP", i, temp
      end if
      do j = 1, nvar
        d(j) = d(j) + temp * u(j, i)
      end do
    end do
    ddx = dSqrt (ddot(nvar, d, 1, d, 1))
    if ( ((lambda /= zero .and. abs(lambda0 + lambda) <= 1.d-20)) .or. &
        & (lambda /= zero .and. abs(lambda0 + lambda) > 1.d-20) .and. store_ddx > 1.d-9) then
      if (ddx > 2*store_ddx) then
        d = d*store_ddx/ddx
        ddx = 2*store_ddx
      end if
    end if
    if (lambda == zero .and. lambda0 == zero .and. iprnt >= 1) then
10050 format (1 x, "PURE NR-STEP HAS LENGTH", f10.5)
      write (iw, 10050) ddx
    end if
    if (lambda /= zero .and. abs(lambda0 + lambda) > 1.d-20 .and. iprnt >= 1) then
10060 format (1 x, "P-RFO-STEP   HAS LENGTH", f10.5)
      write (iw, 10060) ddx
    end if
    if (lambda /= zero .and. abs(lambda0 + lambda) <= 1.d-20 .and. iprnt >= 1) then
10070 format (1 x, "QA/TRIM-STEP HAS LENGTH", f10.5)
      write (iw, 10070) ddx
    end if
    store_ddx = ddx
    if (ddx < (dmax+tmsix)) then
      xlamd = lambda
      xlamd0 = lambda0
      return
    end if
    if (lambda == zero .and. lambda0 == zero) go to 1000
    if (rscal) go to 1040
1020 lambda = zero
    frodo1 = .false.
    frodo2 = .false.
    sstep = step_loc
    if (eone <= zero) then
      lambda = eone - sstep
    end if
    if (ts .and.-eigit < eone) then
      lambda = -eigit - sstep
    end if
    if (eone > zero) then
      sstep = eone
    end if
    bl = lambda - sstep
    bu = lambda + sstep * half
    do
      fl = zero
      fu = zero
      do i = 1, nvar
        if (i /= it) then
          if ((bl-eigval(i) == zero) .or. &
            & (bu-eigval(i) == zero)) then
            write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
            call mopend ("in FORMD")
            return
          end if
          fl = fl + (fx(i)/(bl-eigval(i))) ** 2
          fu = fu + (fx(i)/(bu-eigval(i))) ** 2
        end if
      end do
      if (ts) then
        if ((bl+eigval(it) == zero) .or. &
          & (bu+eigval(it) == zero)) then
          write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
          call mopend ("in FORMD")
          return
        end if
        fl = fl + (fx(it)/(bl+eigval(it))) ** 2
        fu = fu + (fx(it)/(bu+eigval(it))) ** 2
      end if
      fl = fl - d2max
      fu = fu - d2max
      if (fl*fu < zero) go to 1030
      bl = bl - (eone-bl)
      bu = bu + half * (eone-bu)
      if (bl <=-ssmax) then
        bl = -ssmax
        frodo1 = .true.
      end if
      if (Abs (eone-bu) <= ssmin) then
        bu = eone - ssmin
        frodo2 = .true.
      end if
      if (frodo1 .and. frodo2) exit
    end do
    write (iw,*) "NUMERICAL PROBLEMS IN BRACKETING LAMBDA", eone, bl, bu, &
         & fl, fu
    write (iw,*) " GOING FOR FIXED LEVEL SHIFTED NR STEP..."
   !           BOTH LAMBDA SEARCHES FAILED, GO FOR FIXED LEVEL SHIFTED NR
   !           THIS IS UNLIKELY TO PRODUCE ANYTHING USEFUL, BUT MAYBE WE'RE
    lambda = eone - sfix
    lambda0 = eigit + sfix
    rscal = .true.
    go to 1010
1030 continue
    ncnt = 0
    xlambda = zero
    do
      fl = zero
      fu = zero
      fm = zero
      lambda = half * (bl+bu)
      do i = 1, nvar
        if (i /= it) then
          if ((bl-eigval(i) == zero) .or. &
            & (bu-eigval(i) == zero) .or. &
            & (lambda-eigval(i) == zero)) then
            write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
            call mopend ("in FORMD")
            return
          end if
          fl = fl + (fx(i)/(bl-eigval(i))) ** 2
          fu = fu + (fx(i)/(bu-eigval(i))) ** 2
          fm = fm + (fx(i)/(lambda-eigval(i))) ** 2
        end if
      end do
      if (ts) then
        if ((bl+eigval(it) == zero) .or. &
          & (bu+eigval(it) == zero) .or. &
          & (lambda+eigval(it) == zero)) then
          write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
          call mopend ("in FORMD")
          return
        end if
        fl = fl + (fx(it)/(bl+eigval(it))) ** 2
        fu = fu + (fx(it)/(bu+eigval(it))) ** 2
        fm = fm + (fx(it)/(lambda+eigval(it))) ** 2
      end if
      fl = fl - d2max
      fu = fu - d2max
      fm = fm - d2max
      if (Abs (xlambda-lambda) < sstoll) exit
      ncnt = ncnt + 1
      if (ncnt > 1000) then
        write (iw,*) "TOO MANY ITERATIONS IN LAMBDA BISECT", bl, bu, lambda, &
             & fl, fu
        call mopend ("TOO MANY ITERATIONS IN LAMBDA BISECT IN EF")
        return
      else
        xlambda = lambda
        if (fm*fu < zero) then
          bl = lambda
        end if
        if (fm*fl < zero) then
          bu = lambda
        end if
      end if
    end do
   !
    lambda0 = -lambda
    rscal = .true.
   !     .. End of LOOP ..
    go to 1010
1040 if (ddx == zero) then
      write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
      call mopend ("in FORMD")
      return
    end if
    skal = dmax / ddx
    do i = 1, nvar
      d(i) = d(i) * skal
    end do
    ddx = dSqrt (ddot(nvar, d, 1, d, 1))
    if (iprnt >= 1) then
10080 format (5 x, "CALCULATED STEP SIZE TOO LARGE, SCALED WITH", f9.5)
      write (iw, 10080) skal
    end if
    xlamd = lambda
    xlamd0 = lambda0
end subroutine formd
subroutine gethes (xparam, igthes, iloop, hess, pmat, bmat, grad, geo, loc, &
     & oldf, d, vmode, funct0)
    use molkst_C, only: nvar, natoms, id, nalpha, tleft, time0, ndep, gnorm, &
    keywrd, moperr, line
    use chanel_C, only: iw0, iw
    use ef_C, only: ef_mode, nstep, iprnt
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: igthes, iloop
    double precision, intent (inout) :: funct0
    integer, dimension (2, nvar), intent (in) :: loc
    double precision, dimension (nvar*nvar), intent (inout) :: pmat
    double precision, dimension (nvar), intent (inout) :: d, grad, oldf, &
         & vmode, xparam
    double precision, dimension (3, natoms), intent (inout) :: geo
    double precision, dimension (nvar, nvar), intent (inout) :: bmat, hess
   !
   !.. Local Scalars ..
    integer :: i, if, iidum, ij, j, k, l, mtmp, nxxx, percent
    double precision :: dg1, dg2, dg3, dummy, fdmy, funct, funct1, sum, &
         & tdm, time1, time2, tstep, tstore, tt0, test = 0.d0, xinc
    double precision, parameter :: dghsa = 500.d0, dghsd = 200.d0, &
         & dghss = 1000.d0, two = 2.d0,  zzero = 0.d0
    logical :: lpacifier
   !
   !.. Local Arrays ..
    integer, dimension (9) :: ipow = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
    double precision, dimension(:), allocatable :: gnext1, gmin1
    double precision, external :: seconds

    intrinsic Index, Max, Sign
   !
   ! ... Executable Statements ...
   !
    allocate (gnext1(nvar), gmin1(nvar))
    if (index(keywrd, " MOZ") + index(keywrd, " LOCATE-TS")/= 0) then
!
!  MOZYME is not as precise as conventional, so use a bigger step
!
      xinc = 1.d-2
    else
      xinc = 1.d-3
    end if
    do i = 1, nvar
      gnext1(i) = 0.0d0
      gmin1(i) = 0.0d0
    end do
   !
   !     DGHSX IS HESSIAN DIAGONAL FOR IGTHES=0 (STRETCHING, ANGLE,
   !     DIHEDRAL).  THE VALUES SHOULD BE 'OPTIMUM' FOR CYCLOHEXANONE
   !     XINC IS STEPSIZE FOR HESSIAN CALCULATION. TESTS SHOWS THAT IT SHOU
   !     BE IN THE RANGE 10(-2) TO 10(-4). 10(-3) APPEARS TO BE
   !     A REASONABLE COMPROMISE BETWEEN ACCURACY AND NUMERICAL PROBLEMS
    if (igthes == 0) then
10000 format (/, 10 x, "DIAGONAL MATRIX USED AS START HESSIAN", /)
      write (iw, 10000)
      do i = 1, nvar
        do j = 1, nvar
          hess(i, j) = zzero
        end do
      end do
      !
      !   If a polymer, then calculate the diagonal of the Hessian explicitly
      !
      if (id /= 0) then
        do i = 1, nvar
            !
            !  MAKE THE FIRST STEP A WEAK FUNCTION OF THE GRADIENT
            !
          oldf(i) = xparam(i) - Sign (0.001d0, grad(i))
         !
         !  MAKE THE FIRST STEP A WEAK FUNCTION OF THE GRADIENT
         !
        end do
        call compfg (oldf, .true., funct1, .true., gnext1, .true.)
        do i = 1, nvar
          hess(i, i) = Max (100.d0, (grad(i)-gnext1(i))/(xparam(i)-oldf(i)))
        end do
        if (funct1 < funct0) then
            !
            !   The new point is better than the old point
            !
          funct0 = funct1
          do i = 1, nvar
            xparam(i) = oldf(i)
            grad(i) = gnext1(i)
          end do
        end if
      else
        sum = 200.d0
        if (Index (keywrd, " XYZ") /= 0) then
          dg1 = sum
          dg2 = sum
          dg3 = sum
        else
          dg1 = dghss
          dg2 = dghsa
          dg3 = dghsd
        end if
        ij = 1
        outer_loop: do j = 1, natoms
          do i = 1, 3
            if (loc(2, ij) == i .and. loc(1, ij) == j) then
              if (i == 1) then
                hess(ij, ij) = dg1
              end if
              if (i == 2) then
                hess(ij, ij) = dg2
              end if
              if (i == 3) then
                hess(ij, ij) = dg3
              end if
              ij = ij + 1
              if (ij > nvar) exit outer_loop
            end if
          end do
        end do outer_loop
        ij = ij - 1
        if (ij /= nvar) then
          write (iw,*) "ERROR IN IGTHES=0,IJ,NVAR", ij, nvar
        end if
      end if
    end if
   !
    if (igthes == 2) then
10010 format (/, 10 x, "HESSIAN READ FROM DISK", /)
      write (iw, 10010)
      ipow(9) = 0
      ! USE DUMMY ARRAY FOR CALL EXCEPT FOR HESSIAN
      ! TEMPORARY SET NALPHA = 0, THEN WE CAN READ HESSIAN FROM RHF
      ! RUN FOR USE IN SAY UHF RUNS
      ! ALSO SAVE MODE, TO ALLOW FOLLOWING A DIFFERENT MODE THAN THE ON
      ! CURRENTLY ON RESTART FILE
      nxxx = nalpha
      nalpha = 0
      mtmp = ef_mode
      tdm = 0.d0
      fdmy = 0.d0
      iidum = 0
      call efsav (tdm, hess, fdmy, gnext1, gmin1, pmat, iidum, bmat, &
           & ipow, oldf, d, vmode)
      if (moperr) return
      nalpha = nxxx
      ef_mode = mtmp
      nstep = 0
    end if
    if ((igthes == 1) .or. (igthes == 3) .or. (igthes == 4)) then
      ! IF IGTHES IS .EQ. 4, THEN THIS IS A HESSIAN RESTART.
      ! USE GNEXT1 AND DUMMY FOR CALLS TO COMPFG DURING HESSIAN
      ! CALCULATION
      if (igthes == 1) then
10020   format (/, 10 x, "HESSIAN CALCULATED NUMERICALLY", /)
        write (iw, 10020)
      end if
      if (igthes == 3) then
10030   format (/, 10 x, "HESSIAN CALCULATED DOUBLE NUMERICALLY", /)
        write (iw, 10030)
      end if
      if (iprnt >= 5) then
        write (iw, "(I3,12(8F9.4,/3X))") 0, (grad(if), if=1, nvar)
      end if
      time1 = seconds (1)
      tstore = time1
      percent = (100*(iloop-1))/nvar
      lpacifier = .false.
      do i = iloop, nvar
        xparam(i) = xparam(i) + xinc
        call compfg (xparam, .true., dummy, .true., gnext1, .true.)
        if (iprnt >= 5) then
          write (iw, "(I3,12(8F9.4,/3X))") i, (gnext1(if), if=1, nvar)
        end if
        xparam(i) = xparam(i) - xinc
        if (igthes == 3) then
          xparam(i) = xparam(i) - xinc
          gnorm = 0.d0
          call compfg (xparam, .true., dummy, .true., gmin1, .true.)
          if (iprnt >= 5) then
            write (iw, "(I3,12(8F9.4,/3X))") - i, (gmin1(if), if=1, nvar)
          end if
          xparam(i) = xparam(i) + xinc
          do j = 1, nvar
            hess(i, j) = (gnext1(j)-gmin1(j)) / (xinc+xinc)
          end do
        else
          do j = 1, nvar
            hess(i, j) = (gnext1(j)-grad(j)) / xinc
          end do
        end if
        time2 = seconds (1)
        tstep = time2 - time1
        tleft = tleft - tstep
        time1 = time2
        j = (100*i)/nvar
        test = test + tstep
        if (j/5 > percent/5 .or. (j == 1 .and. percent /= 1) .or. i == 1) then
!
!   User pacifier
!
          percent = j
          j = i - iloop + 1
          if (nvar*test/(j*3600) > 0.1d0 .or. lpacifier) then
            lpacifier = .true.
            sum = (nvar - i)*(test/(j*3600))
            write(line,"(a,i3,a,f7.2,a)") &
              "    Hessian",percent,"% complete.  Estimated remaining time required:", sum, " hours"
            write(iw,"(a)")trim(line)
            call to_screen(line)
          end if
        end if
        if (nvar*test/(j*3600) > 0.1d0 .or. lpacifier) then
          write(line,"(i5,a,i4,a)")i," of",nvar," steps completed"
          write(iw,"(a)")trim(line)
          if (iw0 >= 0) call to_screen(line)
        end if

        if (tleft < tstep*two) then
            !
            !  STORE PARTIAL HESSIAN PATRIX
            !  STORE GRADIENT FOR GEOMETRY AND ILOOP AS POSITIVE
          write (iw, "(A)") " NOT ENOUGH TIME TO COMPLETE HESSIAN"
          write (iw, "(A,I4)") " STOPPING IN HESSIAN AT COORDINATE:", i
          ipow(9) = 1
          ipow(2) = nstep
          tt0 = seconds (1) - time0
          j = i
          call efsav (tt0, hess, funct, grad, xparam, pmat, j, bmat, &
               & ipow, oldf, d, vmode)
          tleft = -0.1d0
          call mopend("NOT ENOUGH TIME TO COMPLETE HESSIAN")
          return
        end if
      end do
      !     FIX LAST ENTRY IN GEO ARRAY, THIS IS CURRENTLY AT VALUE-XINC
      k = loc(1, nvar)
      l = loc(2, nvar)
      geo(l, k) = xparam(nvar)
      if (ndep /= 0) then
        call symtry ()
      end if
      !        ADD ALL TIME USED BACK TO TLEFT, THIS WILL THEN BE SUBTRACTED
      !        AGAIN IN MAIN EF ROUTINE
      time2 = seconds (1)
      tstep = time2 - tstore
      tleft = tleft + tstep
    end if
   !
   !     SYMMETRIZE HESSIAN
   !
    do i = 1, nvar
      do j = 1, i - 1
        hess(i, j) = (hess(i, j)+hess(j, i)) / two
        hess(j, i) = hess(i, j)
      end do
    end do
    deallocate (gnext1, gmin1)

end subroutine gethes
subroutine overlp (dmax, ddmin, newmod, nvar, lorjk, u, vmode)
    use molkst_C, only: numcal
    use chanel_C, only: iw
    use ef_C, only: ef_mode, iprnt, omin
    implicit none
    logical, intent (inout) :: lorjk
    integer, intent (in) :: nvar
    integer, intent (out) :: newmod
    double precision, intent (in) :: ddmin, dmax
    double precision, dimension (nvar), intent (inout) :: vmode
    double precision, dimension (nvar, nvar), intent (in) :: u
    integer :: i
    integer, save :: icalcn = 0, it
    double precision :: ovlp, tovlp
    intrinsic Abs
    double precision, external :: dot
   !
   ! ... Executable Statements ...
   !
   !  ON THE FIRST STEP SIMPLY DETERMINE WHICH MODE TO FOLLOW
   !
   !     IF(NSTEP.EQ.1) THEN
    if (icalcn /= numcal) then
      icalcn = numcal !
      if (ef_mode > nvar) then
        write (iw,*) "ERROR!! MODE IS LARGER THAN NVAR", ef_mode
        call mopend ("ERROR!! MODE IS LARGER THAN NVAR")
        return
      else
        it = ef_mode
        if (iprnt >= 1) then
10000     format (5 x, "HESSIAN MODE FOLLOWING SWITCHED ON"/ &
               & "     FOLLOWING MODE ", i3)
          write (iw, 10000) ef_mode
        end if
      end if
    else
      !
      !  ON SUBSEQUENT STEPS DETERMINE WHICH HESSIAN EIGENVECTOR HAS
      !  THE GREATEST OVERLAP WITH THE MODE WE ARE FOLLOWING
      !
      it = 1
      lorjk = .false.
      tovlp = dot (u(1, 1), vmode, nvar)
      tovlp = Abs (tovlp)
      do i = 2, nvar
        ovlp = dot (u(1, i), vmode, nvar)
        ovlp = Abs (ovlp)
        if (ovlp > tovlp) then
          tovlp = ovlp
          it = i
        end if
      end do
      !
      if (iprnt >= 1) then
10010   format (5 x, "OVERLAP OF CURRENT MODE", i3, &
             & " WITH PREVIOUS MODE IS ", f6.3)
        write (iw, 10010) it, tovlp
      end if
      if (tovlp < omin) then
        if (dmax > ddmin) then
          lorjk = .true.
          if (iprnt >= 1) then
10020       format (5 x, "OVERLAP LESS THAN OMIN", f6.3, &
                 & " REJECTING PREVIOUS STEP")
            write (iw, 10020) omin
          end if
          return
        else if (iprnt >= 1) then
10030     format (5 x, "OVERLAP LESS THAN OMIN", f6.3, " BUT TRUST RADIUS", &
               & f6.3, " IS LESS THAN DDMIN", f6.3,/, 5 x, " ACCEPTING STEP")
          write (iw, 10030) omin, dmax, ddmin
        end if
      end if
    end if
   !
   !  SAVE THE EIGENVECTOR IN VMODE
   !
    do i = 1, nvar
      vmode(i) = u(i, it)
    end do
   !
    newmod = it
!
end subroutine overlp
subroutine prjfc (f, xparam, nvar, cof, p, atmass, x, rm, dx, coord)
    use molkst_C, only: numat
    use chanel_C, only: iw
    implicit none
    integer, intent (in) :: nvar
    double precision, dimension (nvar), intent (in) :: xparam
    double precision, dimension (nvar), intent (inout) :: dx, rm, x
    double precision, dimension (numat), intent (in) :: atmass
    double precision, dimension (3, numat), intent (inout) :: coord
    double precision, dimension (nvar, nvar), intent (inout) :: cof, f, p
    double precision, parameter :: zero = 0.0d+00
    double precision, parameter :: one = 1.0d+00
    double precision, parameter :: cut8 = 1.0d-08
    integer :: i, ia, ib, ic, ii, ij, indx, info, ip, j, ja, jb, jc, &
         & jend, jj, jndx, jp, l, natm_loc, nc1
    double precision :: chk, det, sum, tmp, totm, trp
    integer, dimension (6) :: iscr_loc
    double precision, dimension (2) :: det2
    double precision, dimension (3) :: cmass
    double precision, dimension (3, 3) :: rot, scr
    double precision, dimension (3, 3, 3), save :: tens
    intrinsic Abs, Sqrt
    data tens / 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, -1.0d0, &
   & 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d+0, &
   & 0.0d0, -1.0d0, 0.0d0, 0.0d0, 0.0d0, -1.0d0, 0.0d0, &
   & 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /
    natm_loc = nvar / 3
    nc1 = nvar
    ij = 1
    do i = 1, natm_loc
      coord(1, i) = xparam(ij)
      coord(2, i) = xparam(ij+1)
      coord(3, i) = xparam(ij+2)
      ij = ij + 3
    end do
   !     CALCULATE 1/SQRT(MASS)
    l = 0
    do i = 1, natm_loc
      tmp = one / Sqrt (atmass(i))
      do j = 1, 3
        l = l + 1
        rm(l) = tmp
      end do
    end do
   !     PREPARE GRADIENT
    do i = 1, nc1
      dx(i) = zero
    end do
   !     FIND CMS AND CALCULATED MASS WEIGHTED COORDINATES
    totm = zero
    cmass(1) = zero
    cmass(2) = zero
    cmass(3) = zero
    do i = 1, natm_loc
      totm = totm + atmass(i)
      do j = 1, 3
        cmass(j) = cmass(j) + atmass(i) * coord(j, i)
      end do
    end do
    do j = 1, 3
      cmass(j) = cmass(j) / totm
    end do
    l = 0
    do i = 1, natm_loc
      do j = 1, 3
        tmp = Sqrt (atmass(i))
        l = l + 1
        x(l) = tmp * (coord(j, i)-cmass(j))
      end do
    end do
   ! 2. COMPUTE INERTIA TENSOR.
    do i = 1, 3
      do j = 1, 3
        rot(i, j) = zero
      end do
    end do
    do i = 1, natm_loc
      l = 3 * (i-1) + 1
      rot(1, 1) = rot(1, 1) + x(l+1) ** 2 + x(l+2) ** 2
      rot(1, 2) = rot(1, 2) - x(l) * x(l+1)
      rot(1, 3) = rot(1, 3) - x(l) * x(l+2)
      rot(2, 2) = rot(2, 2) + x(l) ** 2 + x(l+2) ** 2
      rot(2, 3) = rot(2, 3) - x(l+1) * x(l+2)
      rot(3, 3) = rot(3, 3) + x(l) ** 2 + x(l+1) ** 2
    end do
    rot(2, 1) = rot(1, 2)
    rot(3, 1) = rot(1, 3)
    rot(3, 2) = rot(2, 3)
   !
   !CHECK THE INERTIA TENSOR.
    chk = rot(1, 1) * rot(2, 2) * rot(3, 3)
    if (Abs (chk) > cut8) then
      !
      ! 4. COMPUTE INVERSION MATRIX OF ROT.
      !     CALL MXLNEQ(ROT,3,3,DET,JRNK,EPS,SCR,+0)
      !     IF(JRNK.LT.3) CALL MOPEND
      !      IF(JRNK.LT.3) RETURN
      info = 0
      call dgefa (rot, 3, 3, iscr_loc, info)
      if (info /= 0) then
        if (numat == 3 .and. natm_loc == 3) then
          call mopend ("ERROR DETECTED IN EF: SYSTEM IS TRIATOMIC AND ALMOST LINEAR")
          write(iw,'(/10x,a)')"Nine geometric parameters marked for optimization."
          write(iw,'(10x,a)')"Reduce this to three or four."
        else
          call mopend ("ERROR DETECTED IN EF: SYSTEM IS ALMOST BUT NOT EXACTLY LINEAR")
          write(iw,'(/10x,a)')"Either define it as being linear or define it as being non-linear."
        end if
      end if
      if (info /= 0) return
      det = zero
      call dgedi (rot, 3, 3, iscr_loc, det2, scr, 1)
    else if (Abs (rot(1, 1)) > cut8) then
      ! X.NE.0
      if (Abs (rot(2, 2)) > cut8) then
         !* 4. X,Y.NE.0 BUT Z=0
        det = rot(1, 1) * rot(2, 2) - rot(1, 2) * rot(2, 1)
        trp = rot(1, 1)
        rot(1, 1) = rot(2, 2) / det
        rot(2, 2) = trp / det
        rot(1, 2) = -rot(1, 2) / det
        rot(2, 1) = -rot(2, 1) / det
      else if (Abs (rot(3, 3)) > cut8) then
         !* 5. X,Z.NE.0 BUT Y=0
        det = rot(1, 1) * rot(3, 3) - rot(1, 3) * rot(3, 1)
        trp = rot(1, 1)
        rot(1, 1) = rot(3, 3) / det
        rot(3, 3) = trp / det
        rot(1, 3) = -rot(1, 3) / det
        rot(3, 1) = -rot(3, 1) / det
      else
         !* 3. Y,Z=0 BUT X.NE.0
        rot(1, 1) = one / rot(1, 1)
      end if
   ! X=0
    else if (Abs (rot(2, 2)) > cut8) then
      ! Y.NE.0
      if (Abs (rot(3, 3)) > cut8) then
         !* 6. Y,Z.NE.0 BUT X=0
        det = rot(3, 3) * rot(2, 2) - rot(3, 2) * rot(2, 3)
        trp = rot(3, 3)
        rot(3, 3) = rot(2, 2) / det
        rot(2, 2) = trp / det
        rot(3, 2) = -rot(3, 2) / det
        rot(2, 3) = -rot(2, 3) / det
      else
         !* 2. X,Z=0 BUT Y.NE.0
        rot(2, 2) = one / rot(2, 2)
      end if
   ! X,Y=0
    else if (Abs (rot(3, 3)) > cut8) then
      !
      !* 1. X,Y=0 BUT Z.NE.0
      rot(3, 3) = one / rot(3, 3)
    else
10000 format (1 x, "EVERY DIAGONAL ELEMENTS ARE ZERO ?", 3 f20.10)
      write (iw, 10000) rot(1, 1), rot(2, 2), rot(3, 3)
      return
    end if
   !
   ! 5. TOTAL MASS ---> TOTM.
   !
   ! 6. COMPUTE P MATRIX
   !    ----------------
    do ip = 1, natm_loc
      indx = 3 * (ip-1)
      do jp = 1, ip
        jndx = 3 * (jp-1)
        do ic = 1, 3
          jend = 3
          if (jp == ip) then
            jend = ic
          end if
          do jc = 1, jend
            sum = zero
            do ia = 1, 3
              do ib = 1, 3
                if (tens(ia, ib, ic) /= 0.d0) then
                  do ja = 1, 3
                    do jb = 1, 3
                      if (tens(ja, jb, jc) /= 0.d0) then
                        sum = sum + tens(ia, ib, ic) * tens(ja, jb, jc) * &
                       & rot(ia, ja) * x(indx+ib) * x(jndx+jb)
                      end if
                    end do
                  end do
                end if
              end do
            end do
            ii = indx + ic
            jj = jndx + jc
            p(ii, jj) = sum + dx(ii) * dx(jj)
            if (ic == jc) then
              p(ii, jj) = p(ii, jj) + one / (rm(ii)*rm(jj)*totm)
            end if
          end do
        end do
      end do
    end do
   !
   ! 7. COMPUTE DELTA(I,J)-P(I,J)
    do i = 1, nc1
      do j = 1, i
        p(i, j) = -p(i, j)
        if (i == j) then
          p(i, j) = one + p(i, j)
        end if
      end do
    end do
   !
   ! 8. NEGLECT SMALLER VALUES THAN 10**-8.
    do i = 1, nc1
      do j = 1, i
        if (Abs (p(i, j)) < cut8) then
          p(i, j) = zero
        end if
        p(j, i) = p(i, j)
      end do
    end do

! For GPU MOPAC
! GBR_new_addition
   !     USE COF FOR SCRATCH.
    call dgemm ("N", "N", nc1, nc1, nc1, 1.0d0, f, nvar, p, nvar, &
         & 0.0d0, cof, nvar)
   !
   ! 11. COMPUTE P*F*P.
   !
    call dgemm ("N", "N", nc1, nc1, nc1, 1.0d0, p, nvar, cof, nvar, &
         & 0.0d0, f, nvar)
    continue
    return
!
end subroutine prjfc
subroutine prthes (eigval, nvar, hess, u)
    use chanel_C, only: iw
    use ef_C, only: iprnt
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, intent (in) :: nvar
    double precision, dimension (nvar), intent (in) :: eigval
    double precision, dimension (nvar, nvar), intent (in) :: hess, u
   !
   !.. Local Scalars ..
    integer :: i, j, low, nup
   !
   !.. Intrinsic Functions ..
    intrinsic Min
   !
   ! ... Executable Statements ...
   !
    if (iprnt >= 4) then
      write (iw,*) " "
      write (iw,*) "              HESSIAN MATRIX"
      low = 1
      nup = 8
      do
        nup = Min (nup, nvar)
10000   format (/, 3 x, 8 i9)
        write (iw, 10000) (i, i=low, nup)
        do i = 1, nvar
10010     format (1 x, i3, 8 f9.1)
          write (iw, 10010) i, (hess(i, j), j=low, nup)
        end do
        nup = nup + 8
        low = low + 8
        if (low > nvar) exit
      end do
    end if
    write (iw,*) " "
    write (iw,*) "              HESSIAN EIGENVALUES AND -VECTORS"
    low = 1
    nup = 8
    do
      nup = Min (nup, nvar)
      write (iw, 10000) (i, i=low, nup)
10020 format (/, 4 x, 8 f9.1, /)
      write (iw, 10020) (eigval(i), i=low, nup)
      do i = 1, nvar
10030   format (1 x, i3, 8 f9.4)
        write (iw, 10030) i, (u(i, j), j=low, nup)
      end do
      nup = nup + 8
      low = low + 8
      if (low > nvar) exit
    end do
end subroutine prthes
subroutine updhes (svec, tvec, grad, nvar, iupd, hess, oldf, d)
    use molkst_C, only: numcal
    use chanel_C, only: iw
    use ef_C, only: iprnt, ddx
    implicit none
    integer, intent (in) :: iupd, nvar
    double precision, dimension (nvar), intent (in) :: d, grad, oldf
    double precision, dimension (nvar), intent (inout) :: svec, tvec
    double precision, dimension (nvar, nvar), intent (inout) :: hess
    integer :: i
    integer, save :: icalcn = 0
    integer :: j
    double precision :: dds, ddtd, temp
    double precision, save :: zero = 0.0d0
    double precision, external :: dot, ddot
!


   !
   ! ... Executable Statements ...
   !
   !  UPDATING OF THE HESSIAN
   !  DEPENDS ON CURRENT grad, OLD grad AND THE
   !  CORRECTION VECTOR USED ON THE LAST CYCLE
   !  SVEC & TVEC ARE FOR TEMPORARY STORAGE
   !
   !  2 UPDATING PROCEDURES ARE POSSIBLE
   !  (I)   THE POWELL UPDATE
   !        THIS PRESERVES THE SYMMETRIC CHARACTER OF THE HESSIAN
   !        WHILST ALLOWING ITS EIGENVALUE STRUCTURE TO CHANGE.
   !        IT IS THE DEFAULT UPDATE FOR A TRANSITION STATE SEARCH
   !  (II)  THE BFGS UPDATE
   !        THIS UPDATE HAS THE IMPORTANT CHARACTERISTIC OF RETAINING
   !        POSITIVE DEFINITENESS (NOTE: THIS IS NOT RIGOROUSLY
   !        GUARANTEED, BUT CAN BE CHECKED FOR BY THE PROGRAM).
   !        IT IS THE DEFAULT UPDATE FOR A MINIMUM SEARCH
   !
   !     SWITCH : IUPD
   !       IUPD = 0  :  SKIP UPDATE
   !       IUPD = 1  :  POWELL
   !       IUPD = 2  :  BFGS
   !
    if (icalcn /= numcal) then
      icalcn = numcal
      if (iprnt >= 2) then
        if (iupd == 0) then
10000     format (/, 5x, "HESSIAN IS NOT BEING UPDATED", /)
          write (iw, 10000)
        end if
        if (iupd == 1) then
            !
10010     format (/, 5x, &
               & "HESSIAN IS BEING UPDATED USING THE POWELL UPDATE", /)
          write (iw, 10010)
        end if
        if (iupd == 2) then
10020     format (/, 5x, "HESSIAN IS BEING UPDATED USING THE BFGS UPDATE", /)
          write (iw, 10020)
        end if
      end if
    end if
    if (iupd == 0) return
    do i = 1, nvar
      tvec(i) = zero
    end do
    do j = 1, nvar
      do i = 1, nvar
        tvec(i) = tvec(i) + hess(i, j) * d(j)
      end do
    end do
   !
    if (iupd == 1) then
      !
      !   (I) POWELL UPDATE
      !
      do i = 1, nvar
        tvec(i) = grad(i) - oldf(i) - tvec(i)
        svec(i) = grad(i) - oldf(i)
      end do
      dds = ddx * ddx

      ddtd = dot (tvec, d, nvar)
      if (Abs(dds - zero) < 1.d-20) then
        write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
        call mopend ("in UPDHES")
        return
      end if
      ddtd = ddtd / dds
      !
      do i = 2, nvar
        do j = 1, i - 1
          temp = tvec(i) * d(j) + d(i) * tvec(j) - d(i) * ddtd * d(j)
          hess(i, j) = hess(i, j) + temp / dds
          hess(j, i) = hess(i, j)
        end do
      end do
      do i = 1, nvar
        temp = d(i) * (2.0d0*tvec(i)-d(i)*ddtd)
        hess(i, i) = hess(i, i) + temp / dds
      end do
    end if
    if (iupd /= 2) return
   !
   !  (II) BFGS UPDATE
   !
    do i = 1, nvar
      svec(i) = grad(i) - oldf(i)
    end do
    dds = ddot (nvar, svec, 1, d, 1)
    ddtd = ddot (nvar, d, 1, tvec, 1)
!
   !
    if (Abs(dds - zero) < 1.d-20 .or. Abs(ddtd) < 1.d-20) then
      write (iw, "(' CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE')")
      call mopend ("in UPDHES")
      return
    end if
    do i = 2, nvar
      do j = 1, i - 1
        temp = (svec(i)*svec(j)) / dds - (tvec(i)*tvec(j)) / ddtd
        hess(i, j) = hess(i, j) + temp
        hess(j, i) = hess(i, j)
      end do
    end do
    do i = 1, nvar
      temp = (svec(i)*svec(i)) / dds - (tvec(i)*tvec(i)) / ddtd
      hess(i, i) = hess(i, i) + temp
    end do
end subroutine updhes
