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

    subroutine lbfgs (xparam, escf)
!
!  Use the limited-memory quasi-Newton Broyden-Fletcher-Goldfarb-Shanno method for unconstrained optimization
!
!    J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980),
!    Mathematics of Computation 35, pp. 773-782.
!    D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale Optimization (1989),
!    Mathematical Programming B, 45, 3, pp. 503-528.
!
      use molkst_C, only: tleft, time0, iflepo, tdump, gnorm, natoms, keywrd, stress, &
      moperr, nvar, id, line, last, mozyme, numat, prt_gradients, keywrd_txt,  prt_coords
!
      use chanel_C, only: iw0, iw, log, ilog, input_fn
!
      use common_arrays_C, only : nc, grad, loc, geo
!
      use ef_C, only : nstep
!
      implicit none
      double precision, intent (inout) :: escf
      double precision, dimension (nvar), intent (inout) ::  xparam
!
      character :: txt, csave*60, task*60, line1*30
      logical :: resfil, times, lsave(4), geo_ref, restrt
      integer :: i, icyc, niwa, nwa, maxcyc, alloc_stat, m, isave(44), nflush = 1, &
        n_bad = 0, max_bad
!
      double precision :: const, cycmx, stepmx, sum, tstep, tt0, rms, slog, &
     & tlast, tolerg, tprt, tx1, tx2, best_funct, best_gnorm, oldstp(12), dsave(29)
!
      double precision, dimension(:), allocatable :: best_xparam, best_gradients, &
      bot, gold, top, xold, wa
      integer, dimension (:), allocatable :: best_nc, iwa, nbd
! For Mopac BLAS
      double precision, external :: ddot, reada, seconds
!

      save :: resfil, tlast
!
      allocate (best_xparam(nvar), best_gradients(nvar), best_nc(natoms), &
         & nbd(nvar), stat=alloc_stat)
      if (alloc_stat /= 0) then
        call memory_error ("mod_lbfgs")
        return
      end if
      m = 12
      grad = 0.d0
      best_xparam(:) = 0.d0
      best_gradients(:) = 0.d0
      best_nc(:) = 0
      best_funct = 1.d10
      best_gnorm = 1.d10
      niwa = 3 * nvar
      nwa = 2 * nvar * m + 4 * nvar + 11 * m ** 2 + 8 * m
      allocate (bot(nvar), gold(nvar), top(nvar), xold(nvar), wa(nwa), &
        iwa(niwa), stat=alloc_stat)
      wa = 0.d0
      iwa = 0
      task = " Unused"
      csave = " Unused"
      lsave = .false.
      isave = 0
      geo_ref = ((Index (keywrd_txt, " GEO_REF") /= 0) .and. (Index(keywrd, "LOCATE-TS") == 0))
      times = (Index (keywrd, " TIMES") /= 0)
      maxcyc = 100000
      slog = 0.2d0
      if (index(keywrd," SLOG=") /= 0) then
        slog = reada(keywrd, index(keywrd," SLOG="))
      else if (index(keywrd," SLOG") /= 0) then
        slog = 0.05d0
      end if
      if (Index (keywrd, " CYCLES") /= 0) then
        maxcyc = Nint (reada (keywrd, Index (keywrd, " CYCLES")))
      end if
      if (index(keywrd, " LET") + index(keywrd, " PREC") /= 0) then
        if (index(keywrd, " LET(") /= 0) then
          max_bad = Nint (reada (keywrd, Index (keywrd, " LET(")))
          write(iw,'(/10x,a,i4)') &
            "Keyword ""LET"" re-set the number of cycles used in identifying the lowest-energy geometry to",max_bad
        else
          max_bad = 60
        end if
      else
        max_bad = 30
      end if
      if (Index (keywrd, "GNORM=") /= 0) then
        tolerg = reada (keywrd, Index (keywrd, "GNORM="))
        if (Index (keywrd, " LET") == 0 .and. tolerg < 1.d-2) then
          write (iw, "(/,A)") "  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01"
          tolerg = 1.d-2
        end if
      else
        tolerg = 1.d0
        if (id /= 0) tolerg = id*2.d0 - 1.d0
        if (Index (keywrd, " PREC") /= 0) then
          tolerg = tolerg*0.2d0
        end if
      end if
      do i = 1, 10
        oldstp(i) = 1.d0
      end do
      tlast = tleft
      tx2 = seconds (2)
      tx1 = tx2
      tleft = tleft - tx2 + time0
    !
    !  Turn OFF all bounds checking.  This saves memory and speeds things up
    !
      do i = 1, nvar
        nbd(i) = 0
      end do
      restrt = (Index (keywrd, " RESTART") /= 0)
      if (restrt) then
        isave = 0
        dsave = 0.d0
        nstep = 0
        tt0 = 0.d0
        call lbfsav (tt0, 0, wa, nwa, iwa, niwa, task, csave, lsave, isave, dsave, &
             & nstep, escf)
        time0 = time0 - tt0
        if (Index (keywrd, " 1SCF") /= 0) then
          i = Index (keywrd, " GRAD") + Index (keywrd,"DERIV")
          call compfg (xparam, .true., escf, .true., grad, i/=0)
          iflepo = 1
          return
        end if
        if (moperr) return
       !
       !  Don't save GOLD - it's only used for advice to the user.
       !
        call dcopy (nvar, grad, 1, gold, 1)
        icyc = nstep
      else
        task = "START"
        nstep = 0
        icyc = 0
      end if
      cycmx = 0.d0
      tlast = tleft
      resfil = .false.
      do
       !
       !  Check: Is there enough time for another cycle?
       !
       !     ------- THE BEGINNING OF THE LOOP ----------
        if (tleft < 1.5d0*cycmx .or. nstep-icyc+1 > maxcyc) then
          if (nstep-icyc+1 > maxcyc) then
            write (iw, "(/20x,a,/)") &
                 & "NUMBER OF CYCLES EXCEEDED.  NOW GOING TO FINAL"
          else
            write (iw, "(20x,a,/30x,a)") &
                 & "THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE", &
                 & "NOW GOING TO FINAL"
          end if
          write (iw, "(//10x,' - THE CALCULATION IS BEING DUMPED TO DISK')")
          write (iw, "(10x,'   RESTART IT USING THE KEYWORD ""RESTART""')")
          tt0 = seconds (1) - time0
          call lbfsav (tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, &
               & dsave, nstep, escf)
          iflepo = -1
          goto 99
        end if
        if (times) then
          call timer (" Before SETULB")
        end if
       !
       !  Yes, there is time for another cycle.
       !
        call dcopy (nvar, xparam, 1, xold, 1)
       !     THIS IS THE CALL TO THE L-BFGS-B CODE.
       !
        call setulb (nvar, m, xparam, bot, top, nbd, escf, grad, 0.d0, &
       & 0.d0, wa, iwa, task,-1, csave, lsave, isave, dsave)
        if (moperr) goto 99
        if (times) then
          call timer (" AFTER SETULB")
        end if
        dsave(2) = dsave(2) + 1.d4
        if (task(1:2) == "FG") then
          if (nstep > 1) then
            !
            !  How big was the step in parameter space?
            !
            sum = 0.d0
            do i = 1, nvar
              sum = sum + (xparam(i)-xold(i)) ** 2
            end do
            sum = Sqrt (sum)
            i = min(nvar,10)
            stepmx = Min (1.d0, dSqrt (ddot(i,oldstp, 1,oldstp, 1)/i))
            if (sum > stepmx*2.d0) then
              !
              ! Step was too big - this was probably due to an error in SETULB.
              ! Reduce the step to two times STEPMX
              !
              const = 2.d0 * stepmx / sum
              do i = 1, nvar
                xparam(i) = const * xparam(i) + (1.d0-const) * xold(i)
              end do
              sum = 2.d0 * stepmx
            end if
            oldstp(Mod(nstep, 10)+1) = sum
          end if
!
!  Limit step to "slog" Angstroms, default: 0.2A.  This should help damp wild swings
!
          do i = 1, nvar
            if (abs(xparam(i) - xold(i)) > slog) then
              xparam(i) = xold(i) + max(-slog, min(slog, xparam(i) - xold(i)))
            end if
          end do
          call compfg (xparam, .true., escf, .true., grad, .true.)
          if (moperr) goto 99
          if (best_funct < escf) then
            n_bad = n_bad + 1
          else
            n_bad = 0
          end if
!
! Allow for quite large regions of instability in the optimization.
!
          if (n_bad > max_bad .or. (gnorm < 1.d0 .and. n_bad > (max_bad*2)/3)) then
            write (iw, "(//10x,' HEAT OF FORMATION IS ESSENTIALLY STATIONARY')")
            iflepo = 3
            exit
          end if
          if (times) then
            call timer ("AFTER COMPFG")
          end if
            !
            !  Write out this cycle
            !
          nstep = nstep + 1
          tx2 = seconds (2)
          tstep = tx2 - tx1
          cycmx = Max (tstep, cycmx)
          tx1 = tx2
          tleft = tleft - tstep
          if (tlast-tleft > tdump) then
            tlast = tleft
            resfil = .true.
            tt0 = seconds (1) - time0
            call lbfsav (tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, &
                 & dsave, nstep, escf)
            if (moperr) goto 99
          end if
          tleft = Max (0.d0, tleft)
          call prttim (tleft, tprt, txt)
          gnorm = dSqrt (ddot(nvar,grad, 1, grad, 1))
          if (best_funct > escf) then
  !
  !  Store best result up to the present
  !
            best_gnorm = gnorm
            best_funct = escf
            best_xparam = xparam
            best_gradients = grad
            best_nc(:natoms) = nc(:natoms)
          end if
  !
  !   Write out current status
  !
          if (id == 3 .and. nstep > 0) then
            nstep = nstep - 1
            call write_cell(iw)
            call write_cell(iw0)
            nstep = nstep + 1
          end if
          if (geo_ref) then
            call  geo_diff(sum, rms, .false.)
            write(iw,'(/1x, a, f8.2, a, f8.4, a, f8.4, a, f11.2, a)') "Difference to Geo-Ref:", sum, &
            " = total,", sum/numat, " = Average,", sqrt(rms/numat)," = RMS movement.    STRESS:", &
            stress
          end if
          if (geo_ref) then
            write(line1,'(a,g15.7)')"  - STRESS:", escf - stress
          else
            line1 = " "
          end if
          if (resfil) then
            write (line, '(" RESTART FILE WRITTEN,      TIME LEFT:", f6.2, &
             & a1, "  GRAD.:", f10.3, " HEAT:", g14.7, a)') &
             tprt, txt, Min (gnorm, 999999.999d0), escf, trim(line1)
            write(iw,"(a)")trim(line)
            call to_screen(trim(line))
            endfile (iw)
            backspace (iw)
            if (log) write (ilog, '(a)', err = 1000)trim(line)
            resfil = .false.
          else
            write (line, '(" CYCLE:", i6, " TIME:", f8.3, " TIME LEFT:", &
                   & f6.2, a1, "  GRAD.:", f10.3, " HEAT:", g14.7, a)') &
                   nstep, Min (tstep, 9999.99d0), tprt, txt, &
                   & Min (gnorm, 999999.999d0), escf, trim(line1)
            write(iw,"(a)")trim(line)
            endfile (iw)
            backspace (iw)
            if (log) write (ilog, "(a)")trim(line)
            call to_screen(trim(line))
          end if
          if (mod(nstep,30) == 0) then
            line = trim(input_fn)
            call add_path(line)
            i = len_trim(line) - 5
            call to_screen(line(:i))
          end if
          if (nflush /= 0) then
            if (Mod(nstep, nflush) == 0) then
              endfile (iw)
              backspace (iw)
              if (log) then
                endfile (ilog)
                backspace (ilog)
              end if
            end if
          end if
          call to_screen("To_file: Geometry optimizing")
          !
          !  Write out the cosine of the angle that the new gradient makes
          !  with the old gradient.  Ideally, this should be small.
          !
  1000    call dcopy (nvar, grad, 1, gold, 1)
          endfile (iw)
          backspace (iw)
          !
          !  EXIT CRITERIA.  (The criteria in SETULB are ignored.)
          if (gnorm < tolerg) then
            iflepo = 3
            exit
          end if
        else if (task(1:5) /= "NEW_X") then
          write (iw, "(2A)") " L-BFGS Message:", task
          iflepo = 9
          exit
        end if
      end do
!
!  If current point is not the best, then load in the best point
!
99    if (best_funct + 5.d-4 < escf) then
        m = 1
        escf   = best_funct
        gnorm  = best_gnorm
        xparam = best_xparam
        grad   = best_gradients
        nc(:natoms)     = best_nc(:natoms)
        if (iflepo == -1) then
!
! Job ran out of time or out of cycles - most likely shut down using the SHUT command.
! so update the geometry before printing it.
!       Update array geo
          do i = 1, nvar
            geo(loc(2,i) ,loc(1,i)) = xparam(i)
          end do
          call symtry
        else
!
!  Call compfg to re-set all calculated quantities
!
          last = 1
          call compfg (xparam, .true., escf, .true., grad, .false.)
          write (iw, "(//10X,'CURRENT BEST VALUE OF HEAT OF FORMATION ='   ,F14.6)") escf - stress
          if ( prt_coords) then
            write(iw,*)
            call geout (iw)
            write(iw,*)
          end if
        end if
      else
        m = 0
        if (iflepo /= -1 .and. .not. mozyme ) then
!
! Job did not run out of time.  The energy levels are incorrect - the SHIFT is still working
! so do a single SCF calculation to re-set the energy levels.
!
          last = 1
          call compfg (xparam, .true., escf, .true., grad, .false.)
        end if
      end if
      if (iflepo == -1) then
        if (m == 1) then
          write (iw, "(//10X,'CURRENT BEST VALUE OF HEAT OF FORMATION ='   ,F14.6)") escf - stress
        else
          write (iw, "(//10X,'CURRENT VALUE OF HEAT OF FORMATION ='   ,F14.6)") escf - stress
        end if
        if (prt_gradients .and. index(keywrd," GRADI") /= 0 .and. mozyme) then
          write (iw, '(3/7X,''CURRENT  POINT  AND  DERIVATIVES'',/)')
          call prtgra ()
        end if
        call geout (iw)
        write(iw,*)
      end if
      deallocate (best_xparam, best_gradients, best_nc, stat=alloc_stat)
      return
    end subroutine lbfgs
    subroutine lbfsav (tt0, mode, wa, nwa, iwa, niwa, task, csave, lsave, isave, &
   & dsave, nstep, escf)
      use molkst_C, only: nscf, numat, norbs, nvar
      use chanel_C, only: ires, iw, restart_fn
      use common_arrays_C, only: xparam, grad
      implicit none
      character (len=60), intent (inout) :: csave, task
      integer, intent (in) :: mode, niwa, nwa
      integer, intent (inout) :: nstep
      double precision, intent (inout) :: escf, tt0
      logical, dimension (4), intent (inout) :: lsave
      integer, dimension (44), intent (inout) :: isave
      integer, dimension (niwa), intent (inout) :: iwa
      double precision, dimension (29), intent (inout) :: dsave
      double precision, dimension (nwa), intent (inout) :: wa
      logical :: opend
      integer :: old_numat, old_norbs, i, j
      inquire (unit=ires, opened=opend)
      if (opend) then
        close (unit=ires, status="KEEP")
      end if
      open (unit=ires, file=restart_fn, status="UNKNOWN", form="UNFORMATTED")
      rewind (ires)
      if (mode == 1) then
        call den_in_out (1)
       !
       !  Write out XPARAM and GRAD
       !
        write (ires) numat, norbs, (xparam(i),i=1,nvar)
        write (ires) grad
       !
       !  Write out restart file
       !
        write (ires) wa, iwa, task, csave, lsave, isave, dsave, nstep, &
       & escf, nscf, tt0
        close (ires)
      else
        write (iw, "(//10X,'RESTORING DATA FROM DISK'/)")
       !
       !  Read in XPARAM and GRAD
       !
        read (ires, iostat = j)old_numat, old_norbs
        if (norbs /= old_norbs .or. numat /= old_numat .or. j /= 0) then
        call mopend("Restart file read in does not match current data set")
        return
    end if
        read (ires, err=1100, END=1100) grad
       !
       !  Read in restart file
       !
        read (ires, err=1000, END=1000) wa, iwa, task, csave, lsave, isave, &
       & dsave, nstep, escf, nscf, tt0
        i = int(tt0/10000000)
    tt0 = tt0 - i*10000000
    write (iw, '(10X,''TOTAL TIME USED SO FAR:'',F13.2,'' SECONDS'',/)') tt0
        return
1000    call mopend ("RESTART FILE EXISTS, BUT IS CORRUPT")
        return
1100    call mopend ("NO RESTART FILE EXISTS!")
      end if
    end subroutine lbfsav
    subroutine setulb (n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, &
   & iprint, csave, lsave, isave, dsave)
      implicit none
      character (len=60), intent (inout) :: csave, task
      integer, intent (in) :: iprint, m, n
      double precision, intent (in) :: factr, pgtol
      double precision, intent (inout) :: f
      logical, dimension (4), intent (inout) :: lsave
      integer, dimension (44), intent (inout) :: isave
      integer, dimension (3*n), intent (inout) :: iwa
      integer, dimension (n), intent (in) :: nbd
      double precision, dimension (2*m*n+4*n+11*m*m+8*m), intent (inout) :: wa
      double precision, dimension (29), intent (inout) :: dsave
      double precision, dimension (n), intent (in) :: l, u
      double precision, dimension (n), intent (inout) :: g, x
      integer :: ld, lr, lsnd, lss, lsy, lt, lwa, lwn, lws, lwt, &
     & lwy, lz
      if (task == "START") then
        isave(1) = m * n
        isave(2) = m ** 2
        isave(3) = 4 * m ** 2
        isave(4) = 1
        isave(5) = isave(4) + isave(1)
        isave(6) = isave(5) + isave(1)
        isave(7) = isave(6) + isave(2)
        isave(8) = isave(7) + isave(2)
        isave(9) = isave(8)
        isave(10) = isave(9) + isave(2)
        isave(11) = isave(10) + isave(3)
        isave(12) = isave(11) + isave(3)
        isave(13) = isave(12) + n
        isave(14) = isave(13) + n
        isave(15) = isave(14) + n
        isave(16) = isave(15) + n
      end if
      lws = isave(4)
      lwy = isave(5)
      lsy = isave(6)
      lss = isave(7)
      lwt = isave(9)
      lwn = isave(10)
      lsnd = isave(11)
      lz = isave(12)
      lr = isave(13)
      ld = isave(14)
      lt = isave(15)
      lwa = isave(16)
    !
    !
      call mainlb (n, m, x, l, u, nbd, f, g, factr, pgtol, wa(lws), wa(lwy), &
     & wa(lsy), wa(lss), wa(lwt), wa(lwn), wa(lsnd), wa(lz), wa(lr), wa(ld), &
     & wa(lt), wa(lwa), iwa(1), iwa(n+1), iwa(2*n+1), task, iprint, csave, &
     & lsave, isave(22), dsave)
    !
    end subroutine setulb
  !
  !======================= The end of setulb =============================
    subroutine mainlb (n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, sy, ss, &
   & wt, wn, snd, z, r, d, t, wa, Index, iwhere, indx2, task, iprint, csave, &
   & lsave, isave, dsave)
    !
    use chanel_C, only: lbfgs_it, iw
      implicit none
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      character (len=60), intent (inout) :: csave, task
      integer, intent (in) :: iprint, m, n
      double precision, intent (in) :: factr, pgtol
      double precision, intent (inout) :: f
      logical, dimension (4), intent (inout) :: lsave
      integer, dimension (23), intent (inout) :: isave
      integer, dimension (n), intent (in) :: nbd
      integer, dimension (n), intent (inout) :: index, indx2, iwhere
      double precision, dimension (29), intent (inout) :: dsave
      double precision, dimension (8*m), intent (inout) :: wa
      double precision, dimension (n), intent (in) :: l, u
      double precision, dimension (n), intent (inout) :: d, g, r, t, x, z
      double precision, dimension (2*m, 2*m), intent (inout) :: snd, wn
      double precision, dimension (m, m), intent (inout) :: ss, sy, wt
      double precision, dimension (n, m), intent (inout) :: ws, wy
      character (len=3) :: word
      logical :: boxed, cnstnd, prjctd, updatd, wrk
      integer :: col, head, i, iback, ifun, ileave, info, itail, iter, itfile, &
     & iupdat, iword, k, nact, nenter, nfgv, nfree, nint, nintol, nskip
      double precision :: cachyt, cpu1, cpu2, ddum, dnorm, dr, dtd, epsmch, &
     & fold, gd, gdold, lnscht, rr, sbgnrm, sbtime, stp, stpmx, theta, time, &
     & time1, time2, tol, xstep
      double precision, external :: dpmeps
! For Parallel MOPAC
      double precision, external :: ddot

      if (task == "START") then
       !
       !         call timer(time1)
       !
       !        Generate the current machine precision.
       !
        epsmch = dpmeps ()!
       !        The end of the initialization.
       !
        fold = 0.0d0
        dnorm = 0.0d0
        time1 = 0.0d0
        time2 = 0.0d0
        cpu1 = 0.0d0
        cpu2 = 0.0d0
        gd = 0.0d0
        sbgnrm = 0.0d0
        stp = 0.0d0
        stpmx = 0.0d0
        gdold = 0.0d0
        dtd = 0.0d0
       !
       !
       !        Initialize counters and scalars when task='START'.
       !
       !           for the limited memory BFGS matrices:
        col = 0
        head = 1
        theta = one
        iupdat = 0
        updatd = .false.
        iback = 0
        itail = 0
        ifun = 0
        iword = 0
        nact = 0
        ileave = 0
        nenter = 0
       !           for operation counts:
        iter = 0
        nfgv = 0
        nint = 0
        nintol = 0
        nskip = 0
        nfree = n
       !
       !           for stopping tolerance:
        tol = factr * epsmch
       !
       !           for measuring running time:
        cachyt = 0
        sbtime = 0
        lnscht = 0
       !           'word' records the status of subspace solutions.
        word = "---"
       !
       !           'info' records the termination information.
        info = 0
        itfile = 0
       !
        if (iprint >= 1) then
          !                                open a summary file
     !     open (lbfgs_it, file=lbfgs_it_fn, status="unknown")
          itfile = lbfgs_it
        end if
       !
       !        Check the input arguments for errors.
       !
        call errclb (n, m, factr, l, u, nbd, task, info, k)
        if (task(1:5) == "ERROR") then
          xstep = 0.0d0
          k = 0
          call prn3lb (n, x, f, task, iprint, info, itfile, iter, nfgv, &
         & nintol, nskip, nact, sbgnrm, zero, nint, word, iback, stp, xstep, &
         & k, cachyt, sbtime, lnscht)
          return
        else
          !
          call prn1lb (n, m, l, u, x, iprint, itfile, epsmch)
          !        Initialize iwhere & project x onto the feasible set.
          call active (n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)
        end if
      else
       !          restore local variables.
       !
        prjctd = lsave(1)
        cnstnd = lsave(2)
        boxed = lsave(3)
        updatd = lsave(4)
       !
        nintol = isave(1)
        itfile = isave(3)
        iback = isave(4)
        nskip = isave(5)
        head = isave(6)
        col = isave(7)
        itail = isave(8)
        iter = isave(9)
        iupdat = isave(10)
        nint = isave(12)
        nfgv = isave(13)
        info = isave(14)
        ifun = isave(15)
        iword = isave(16)
        nfree = isave(17)
        nact = isave(18)
        ileave = isave(19)
        nenter = isave(20)
       !
        theta = dsave(1)
        fold = dsave(2)
        tol = dsave(3)
        dnorm = dsave(4)
        epsmch = dsave(5)
        cpu1 = dsave(6)
        cachyt = dsave(7)
        sbtime = dsave(8)
        lnscht = dsave(9)
        time1 = dsave(10)
        gd = dsave(11)
        stpmx = dsave(12)
        sbgnrm = dsave(13)
        stp = dsave(14)
        gdold = dsave(15)
        dtd = dsave(16)
        cpu2 = 0.d0
       !        After returning from the driver go to the point where execution
       !        is to resume.
       !
        if (task(1:5) /= "FG_LN") then
          if (task(1:5) == "NEW_X") then
            !
            !     Test for termination.
            !
            if (sbgnrm <= pgtol) then
              !             terminate the algorithm.
              task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"
              go to 1100
            else

              ddum = Max (Abs (fold), Abs (f), one)
              if ((fold-f) <= tol*ddum) then
                !           terminate the algorithm.
                task = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
                ! i.e., to issue a warning if iback>10 in the line search.
                if (iback >= 10) then
                  info = -5
                end if
                go to 1100
              else
                !
                !     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
                do i = 1, n
                  r(i) = g(i) - r(i)
                end do
                rr = ddot (n, r, 1, r, 1)
                if (Abs(stp - one) < 1.d-20) then
                  dr = gd - gdold
                  ddum = -gdold
                else
                  dr = (gd-gdold) * stp
                  call dscal (n, stp, d, 1)
                  ddum = -gdold * stp
                end if
                if (dr <= epsmch*ddum) then
                      !                            skip the L-BFGS update.
                  nskip = nskip + 1
                  updatd = .false.
                  if (iprint >= 1) then
10000               format ("  ys=", 1 p, e10.3, "  -gs=", 1 p, e10.3, &
                   & " BFGS update SKIPPED")
                    write (iw, 10000) dr, ddum
                  end if
                else
                  !cccccccccccccccccccccccccccccccccccccccccccccccc
                  !
                  !     Update the L-BFGS matrix.
                  !
                  !cccccccccccccccccccccccccccccccccccccccccccccccc
                  updatd = .true.
                  iupdat = iupdat + 1
                  !
                  !     Update matrices WS and WY and form the middle
                  !     matrix in B.
                  !
                  call matupd (n, m, ws, wy, sy, ss, d, r, itail, iupdat, &
                       & col, head, theta, rr, dr, stp, dtd)
                  !
                  ! Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
                  !    Store T in the upper triangular of the array wt;
                  !    Cholesky factorize T to J*J' with
                  !    J' stored in the upper triangular of wt.
                  !
                  call formt (m, wt, sy, ss, col, theta, info)
                  if (info /= 0) then
                    ! nonpositive definiteness in Cholesky factorization;
                    ! refresh the lbfgs memory and restart the iteration.
                    if (iprint >= 1) then
10010                 format (/, &
                     & " Nonpositive definiteness in Cholesky factorization", &
                     & " in formt;"/ &
                     & "   refresh the lbfgs memory and restart the iteration.")
                      write (iw, 10010)
                    end if
                    info = 0
                    col = 0
                    head = 1
                    theta = one
                    iupdat = 0
                    updatd = .false.
                  end if
                end if
              end if
            end if
          else if (task(1:5) == "FG_ST") then
            nfgv = 1
             !     Compute the infinity norm of the (-) projected gradient.
            call projgr (n, l, u, nbd, x, g, sbgnrm)
            if (iprint >= 1) then
10020         format (/, "At iterate", i5, 4 x, "f= ", 1 p, d12.5, 4 x, &
             & "|proj g|= ", 1 p, d12.5)
              write (iw, 10020) iter, f, sbgnrm
10030         format (2(1 x, i4), 5 x, "-", 5 x, "-", 3 x, "-", 5 x, "-", 5 x, &
             & "-", 8 x, "-", 3 x, 1 p, 2(1 x, d10.3))
              write (itfile, 10030) iter, nfgv, sbgnrm, f
            end if
            if (sbgnrm <= pgtol) then
              !            terminate the algorithm.
              task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"
              go to 1100
            end if
          else if (task(1:4) == "STOP") then
            if (task(7:9) == "CPU") then
                !           restore the previous iterate.
              call dcopy (n, t, 1, x, 1)
              call dcopy (n, r, 1, g, 1)
              f = fold
            end if
            go to 1100
          else
            go to 1000
          end if
          do
            if (iprint >= 99) then
                !
10040         format (/ /, "ITERATION ", i5)
              write (iw, 10040) iter + 1
            end if
            iword = -1
             !
            if ( .not. cnstnd .and. col > 0) then
                !              skip the search for GCP.
              call dcopy (n, x, 1, z, 1)
              wrk = updatd
              nint = 0
            else
              !
              !cccccccccccccccccccccccccccccccccccccccccccccccccccc
              !
              !     Compute the Generalized Cauchy Point (GCP).
              !
              !cccccccccccccccccccccccccccccccccccccccccccccccccccc
              !
              !      call timer(cpu1)
              call cauchy (n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, &
             & ws, sy, wt, theta, col, head, wa(1), wa(2*m+1), wa(4*m+1), &
             & wa(6*m+1), nint, iprint, sbgnrm, info, epsmch)
              if (info /= 0) then
                !singular triangular system detected; refresh the lbfgs memory.
                if (iprint >= 1) then
10050             format (/, " Singular triangular system detected;",/, &
                 & "   refresh the lbfgs memory and restart the iteration.")
                  write (iw, 10050)
                end if
                info = 0
                col = 0
                head = 1
                theta = one
                iupdat = 0
                updatd = .false.
                   !         call timer(cpu2)
                cachyt = cachyt + cpu2 - cpu1
                cycle
              else
                   !      call timer(cpu2)
                cachyt = cachyt + cpu2 - cpu1
                nintol = nintol + nint
                !
                !     Count the entering and leaving variables for iter > 0;
                !     find the index set of free and active variables at
                !     the GCP.
                !
                call freev (n, nfree, index, nenter, ileave, indx2, iwhere, &
               & wrk, updatd, cnstnd, iprint, iter)
                   !
                nact = n - nfree
              end if
            end if
            !     If there are no free variables or B=theta*I, then
            !          skip the subspace minimization.
            if (nfree == 0 .or. col == 0) exit
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Subspace minimization.
            !
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !
             !      call timer(cpu1)
             !
             !     Form  the LEL^T factorization of the indefinite
             !       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
             !                     [L_a -R_z           theta*S'AA'S ]
             !       where     E = [-I  0]
             !                     [ 0  I]
             !
            if (wrk) then
              call formk (n, nfree, index, nenter, ileave, indx2, iupdat, &
             & updatd, wn, snd, m, ws, wy, sy, theta, col, head, info)
            end if
            if (info /= 0) then
                !          nonpositive definiteness in Cholesky factorization;
                !          refresh the lbfgs memory and restart the iteration.
              if (iprint >= 1) then
10060           format (/, &
               & " Nonpositive definiteness in Cholesky factorization in formk;" / &
               & "   refresh the lbfgs memory and restart the iteration.")
                write (iw, 10060)
              end if
              info = 0
              col = 0
              head = 1
              theta = one
              iupdat = 0
              updatd = .false.
                !         call timer(cpu2)
              sbtime = sbtime + cpu2 - cpu1
            else
              !
              !        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
              !                                from 'cauchy').
              call cmprlb (n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, &
                   & col, head, nfree, cnstnd, info)
              if (info == 0) then
                   !       call the direct method.
                call subsm (n, m, nfree, index, l, u, nbd, z, r, ws, wy, &
               & theta, col, head, iword, wa, wn, iprint, info)
              end if
              if (info /= 0) then
                   !    singular triangular system detected;
                   !    refresh the lbfgs memory and restart the iteration.
                if (iprint >= 1) then
                  write (iw, 10050)
                end if
                info = 0
                col = 0
                head = 1
                theta = one
                iupdat = 0
                updatd = .false.
                   !         call timer(cpu2)
                sbtime = sbtime + cpu2 - cpu1
              else
                   !      call timer(cpu2)
                sbtime = sbtime + cpu2 - cpu1
                exit
              end if
            end if
          end do
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !
          !     Line search and optimality tests.
          !
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !     Generate the search direction d:=z-x.
          !
          do i = 1, n
            d(i) = z(i) - x(i)
          end do
        end if
        do
          !      call timer(cpu1)
          call lnsrlb (n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
         & stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, task, &
         & boxed, cnstnd, csave, isave(22), dsave(17))
          if (info == 0 .and. iback < 20) go to 1200
          !          restore the previous iterate.
          call dcopy (n, t, 1, x, 1)
          call dcopy (n, r, 1, g, 1)
          f = fold
          if (col == 0) exit
          !             refresh the lbfgs memory and restart the iteration.
          if (iprint >= 1) then
10070       format (/, " Bad direction in the line search;",/, &
           & "   refresh the lbfgs memory and restart the iteration.")
            write (iw, 10070)
          end if
          if (info == 0) then
            nfgv = nfgv - 1
          end if
          info = 0
          col = 0
          head = 1
          theta = one
          iupdat = 0
          updatd = .false.
          task = "RESTART_FROM_LNSRCH"
          !            call timer(cpu2)
          lnscht = lnscht + cpu2 - cpu1
          do
            ! ----------------- the beginning of the loop ---------------
            if (iprint >= 99) then
              write (iw, 10040) iter + 1
            end if
            iword = -1
             !
            if ( .not. cnstnd .and. col > 0) then
                !                             skip the search for GCP.
              call dcopy (n, x, 1, z, 1)
              wrk = updatd
              nint = 0
            else
              !
              !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
              !
              !     Compute the Generalized Cauchy Point (GCP).
              !
              !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
              !
              !      call timer(cpu1)
              call cauchy (n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, &
             & ws, sy, wt, theta, col, head, wa(1), wa(2*m+1), wa(4*m+1), &
             & wa(6*m+1), nint, iprint, sbgnrm, info, epsmch)
              if (info /= 0) then
                !singular triangular system detected; refresh the lbfgs memory.
                if (iprint >= 1) then
                  write (iw, 10050)
                end if
                info = 0
                col = 0
                head = 1
                theta = one
                iupdat = 0
                updatd = .false.
                   !         call timer(cpu2)
                cachyt = cachyt + cpu2 - cpu1
                cycle
              else
                   !      call timer(cpu2)
                cachyt = cachyt + cpu2 - cpu1
                nintol = nintol + nint
                !
                !     Count the entering and leaving variables for iter > 0;
                !     find the index set of free and active variables at the GCP.
                !
                call freev (n, nfree, index, nenter, ileave, indx2, iwhere, &
               & wrk, updatd, cnstnd, iprint, iter)
                   !
                nact = n - nfree
              end if
            end if
             !     If there are no free variables or B=theta*I, then
             !                skip the subspace minimization.
            if (nfree == 0 .or. col == 0) exit
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !
            !     Subspace minimization.
            !
            !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             !
             !      call timer(cpu1)
             !
             !     Form  the LEL^T factorization of the indefinite
             !       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
             !                     [L_a -R_z           theta*S'AA'S ]
             !       where     E = [-I  0]
             !                     [ 0  I]
             !
            if (wrk) then
              call formk (n, nfree, index, nenter, ileave, indx2, iupdat, &
             & updatd, wn, snd, m, ws, wy, sy, theta, col, head, info)
            end if
            if (info /= 0) then
                !          nonpositive definiteness in Cholesky factorization;
                !          refresh the lbfgs memory and restart the iteration.
              if (iprint >= 1) then
                write (iw, 10060)
              end if
              info = 0
              col = 0
              head = 1
              theta = one
              iupdat = 0
              updatd = .false.
                !         call timer(cpu2)
              sbtime = sbtime + cpu2 - cpu1
            else
              !
              !     compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
              !                       from 'cauchy').
              call cmprlb (n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, &
             & col, head, nfree, cnstnd, info)
              if (info == 0) then
                   !       call the direct method.
                call subsm (n, m, nfree, index, l, u, nbd, z, r, ws, wy, &
               & theta, col, head, iword, wa, wn, iprint, info)
              end if
              if (info /= 0) then
                !          singular triangular system detected;
                !          refresh the lbfgs memory and restart the iteration.
                if (iprint >= 1) then
                  write (iw, 10050)
                end if
                info = 0
                col = 0
                head = 1
                theta = one
                iupdat = 0
                updatd = .false.
                   !         call timer(cpu2)
                sbtime = sbtime + cpu2 - cpu1
              else
                   !      call timer(cpu2)
                sbtime = sbtime + cpu2 - cpu1
                exit
              end if
            end if
          end do
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !
          !     Line search and optimality tests.
          !
          !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !     Generate the search direction d:=z-x.
          !
          do i = 1, n
            d(i) = z(i) - x(i)
          end do
        end do
       !             abnormal termination.
        if (info == 0) then
          info = -9
          !        restore the actual number of f and g evaluations etc.
          nfgv = nfgv - 1
          ifun = ifun - 1
          iback = iback - 1
        end if
        task = "ABNORMAL_TERMINATION_IN_LNSRCH"
        iter = iter + 1
        go to 1100
1200    if (task(1:5) /= "FG_LN") then
          !    calculate and print out the quantities related to the new X.
          !         call timer(cpu2)
          lnscht = lnscht + cpu2 - cpu1
          iter = iter + 1
          !        Compute the infinity norm of the projected (-)gradient.
          call projgr (n, l, u, nbd, x, g, sbgnrm)
          !        Print iteration information.
          !
          call prn2lb (n, x, f, g, iprint, itfile, iter, nfgv, nact, sbgnrm, &
         & nint, word, iword, iback, stp, xstep)
        end if
        go to 1300
1100    time2 = 0
       !      call timer(time2)
        time = time2 - time1
        k = 0
        call prn3lb (n, x, f, task, iprint, info, itfile, iter, nfgv, nintol, &
       & nskip, nact, sbgnrm, time, nint, word, iback, stp, xstep, k, cachyt, &
       & sbtime, lnscht)
        go to 1300
      end if
    !
    !     Compute f0 and g0.
    !
    !          return to the driver to calculate f and g; reenter at 111.
1000  task = "FG_START"
    !
    !     Save local variables.
    !
1300  lsave(1) = prjctd
      lsave(2) = cnstnd
      lsave(3) = boxed
      lsave(4) = updatd
    !
      isave(1) = nintol
      isave(3) = itfile
      isave(4) = iback
      isave(5) = nskip
      isave(6) = head
      isave(7) = col
      isave(8) = itail
      isave(9) = iter
      isave(10) = iupdat
      isave(12) = nint
      isave(13) = nfgv
      isave(14) = info
      isave(15) = ifun
      isave(16) = iword
      isave(17) = nfree
      isave(18) = nact
      isave(19) = ileave
      isave(20) = nenter
    !
      dsave(1) = theta
      dsave(2) = fold
      dsave(3) = tol
      dsave(4) = dnorm
      dsave(5) = epsmch
      dsave(6) = cpu1
      dsave(7) = cachyt
      dsave(8) = sbtime
      dsave(9) = lnscht
      dsave(10) = time1
      dsave(11) = gd
      dsave(12) = stpmx
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd
    !
    end subroutine mainlb
  !======================= The end of mainlb =============================
  !
    subroutine active (n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)
      use chanel_C, only: iw
      implicit none
      double precision, parameter :: zero = 0.0d0
      logical, intent (out) :: boxed, cnstnd, prjctd
      integer, intent (in) :: iprint, n
      integer, dimension (n), intent (in) :: nbd
      integer, dimension (n), intent (out) :: iwhere
      double precision, dimension (n), intent (in) :: l, u
      double precision, dimension (n), intent (inout) :: x
      integer :: i, nbdd
    !
    !     Initialize nbdd, prjctd, cnstnd and boxed.
    !
      nbdd = 0
      prjctd = .false.
      cnstnd = .false.
      boxed = .true.
    !
    !     Project the initial x to the easible set if necessary.
    !
      do i = 1, n
        if (nbd(i) > 0) then
          if (nbd(i) <= 2 .and. x(i) <= l(i)) then
            if (x(i) < l(i)) then
              prjctd = .true.
              x(i) = l(i)
            end if
            nbdd = nbdd + 1
          else if (nbd(i) >= 2 .and. x(i) >= u(i)) then
            if (x(i) > u(i)) then
              prjctd = .true.
              x(i) = u(i)
            end if
            nbdd = nbdd + 1
          end if
        end if
      end do
    !
    !     Initialize iwhere and assign values to cnstnd and boxed.
    !
      do i = 1, n
        if (nbd(i) /= 2) then
          boxed = .false.
        end if
        if (nbd(i) == 0) then
          !                                this variable is always free
          iwhere (i) = -1 !
          !           otherwise set x(i)=mid(x(i), u(i), l(i)).
        else
          cnstnd = .true.
          if (nbd(i) == 2 .and. u(i)-l(i) <= zero) then
             !                   this variable is always fixed
            iwhere(i) = 3
          else
            iwhere(i) = 0
          end if
        end if
      end do
    !
      if (iprint >= 0) then
        if (prjctd) then
          write (iw,*) "The initial X is infeasible.  Restart with its projection."
        end if
        if ( .not. cnstnd) then
          write (iw,*) "This problem is unconstrained."
        end if
      end if
    !
      if (iprint > 0) then
       !
10000   format (/, "At X0 ", i9, " variables are exactly at the bounds")
        write (iw, 10000) nbdd
      end if
    !
    end subroutine active
  !
  !======================= The end of active =============================
    subroutine bmv (m, sy, wt, col, v, p, info)
      implicit none
      integer, intent (in) :: col, m
      integer, intent (inout) :: info
      double precision, dimension (2*col), intent (in) :: v
      double precision, dimension (2*col), intent (inout) :: p
      double precision, dimension (m, m), intent (in) :: sy
      double precision, dimension (m, m), intent (inout) :: wt
      integer :: i, i2, k
      double precision :: sum
      if (col == 0) return
    !     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
    !                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
    !
    !       solve Jp2=v2+LD^(-1)v1.
      p(col+1) = v(col+1)
      do i = 2, col
        i2 = col + i
        sum = 0.0d0
        do k = 1, i - 1
          sum = sum + sy(i, k) * v(k) / sy(k, k)
        end do
        p(i2) = v(i2) + sum
      end do
    !     Solve the triangular system
      call dtrsl (wt, m, col, p(col+1), 11, info)
      if (info /= 0) return
    !       solve D^(1/2)p1=v1.
      do i = 1, col
        p(i) = v(i) / Sqrt (sy(i, i))
      end do
    !     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
    !                    [  0         J'           ] [ p2 ]   [ p2 ].
    !       solve J^Tp2=p2.
      call dtrsl (wt, m, col, p(col+1), 1, info)
      if (info /= 0) return
    !       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
    !                 =-D^(-1/2)p1+D^(-1)L'p2.
      do i = 1, col
        p(i) = -p(i) / Sqrt (sy(i, i))
      end do
      do i = 1, col
        sum = 0.d0
        do k = i + 1, col
          sum = sum + sy(k, i) * p(col+k) / sy(i, i)
        end do
        p(i) = p(i) + sum
      end do
    !
    end subroutine bmv
  !
  !======================== The end of bmv ===============================
  !
    subroutine cauchy (n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, m, wy, &
   & ws, sy, wt, theta, col, head, p, c, wbp, v, nint, iprint, sbgnrm, info, &
   & epsmch)
      use chanel_C, only: iw
      implicit none
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, intent (in) :: head, iprint, m, n
      integer, intent (inout) :: col, info
      integer, intent (out) :: nint
      double precision, intent (in) :: epsmch, sbgnrm
      double precision, intent (inout) :: theta
      integer, dimension (n), intent (in) :: nbd
      integer, dimension (n), intent (inout) :: iorder, iwhere
      double precision, dimension (2*m), intent (inout) :: c, p, v, wbp
      double precision, dimension (n), intent (in) :: g, l, u
      double precision, dimension (n), intent (inout) :: d, t, x, xcp
      double precision, dimension (m, m), intent (in) :: sy
      double precision, dimension (m, m), intent (inout) :: wt
      double precision, dimension (n, col), intent (in) :: ws, wy
      logical :: bnded, xlower, xupper
      integer :: col2, i, ibkmin, ibp, iter, j, nbreak, nfree, nleft, pointr
      double precision :: bkmin, dibp, dibp2, dt, dtm, f1, f2, f2_org, neggi, &
     & tj, tj0, tl, tsum, tu, wmc, wmp, wmw, zibp
      external daxpy, dcopy, dscal
      double precision, external :: ddot
      tl = 0.d0
      tu = 0.d0
    !
    !     Check the status of the variables, reset iwhere(i) if necessary;
    !       compute the Cauchy direction d and the breakpoints t; initialize
    !       the derivative f1 and the vector p = W'd (for theta = 1).
      if (sbgnrm <= zero) then
        if (iprint >= 0) then
          write (iw,*) "Subgnorm = 0.  GCP = X."
        end if
        call dcopy (n, x, 1, xcp, 1)
      else
        bnded = .true.
        nfree = n + 1
        nbreak = 0
        ibkmin = 0
        bkmin = zero
        col2 = 2 * col
        f1 = zero
        if (iprint >= 99) then
10000     format (/, "---------------- CAUCHY entered-------------------")
          write (iw, 10000)
        end if
       !
       !     We set p to zero and build it up as we determine d.
       !
        do i = 1, col2
          p(i) = zero
        end do
       !
       !     In the following loop we determine for each variable its bound
       !        status and its breakpoint, and update p accordingly.
       !        Smallest breakpoint is identified.
       !
        do i = 1, n
          neggi = -g(i)
          if (iwhere(i) /= 3 .and. iwhere(i) /=-1) then
             !             if x(i) is not a constant and has bounds,
             !             compute the difference between x(i) and its bounds.
            if (nbd(i) <= 2) then
              tl = x(i) - l(i)
            end if
            if (nbd(i) >= 2) then
              tu = u(i) - x(i)
            end if
             !
             !           If a variable is close enough to a bound
             !             we treat it as at bound.
            xlower = nbd(i) <= 2 .and. tl <= zero
            xupper = nbd(i) >= 2 .and. tu <= zero
             !
             !              reset iwhere(i).
            iwhere(i) = 0
            if (xlower) then
              if (neggi <= zero) then
                iwhere(i) = 1
              end if
            else if (xupper) then
              if (neggi >= zero) then
                iwhere(i) = 2
              end if
            else if (Abs (neggi) <= zero) then
              iwhere(i) = -3
            end if
          end if
          pointr = head
          if (iwhere(i) /= 0 .and. iwhere(i) /=-1) then
            d(i) = zero
          else
            d(i) = neggi
            f1 = f1 - neggi * neggi
             !             calculate p := p - W'e_i* (g_i).
            do j = 1, col
              p(j) = p(j) + wy(i, pointr) * neggi
              p(col+j) = p(col+j) + ws(i, pointr) * neggi
              pointr = Mod (pointr, m) + 1
            end do
            if (nbd(i) <= 2 .and. nbd(i) /= 0 .and. neggi < zero) then
                !            x(i) + d(i) is bounded; compute t(i).
              nbreak = nbreak + 1
              iorder(nbreak) = i
              t(nbreak) = tl / (-neggi)
              if (nbreak == 1 .or. t(nbreak) < bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
              end if
            else if (nbd(i) >= 2 .and. neggi > zero) then
                !             x(i) + d(i) is bounded; compute t(i).
              nbreak = nbreak + 1
              iorder(nbreak) = i
              t(nbreak) = tu / neggi
              if (nbreak == 1 .or. t(nbreak) < bkmin) then
                bkmin = t(nbreak)
                ibkmin = nbreak
              end if
            else
                !                x(i) + d(i) is not bounded.
              nfree = nfree - 1
              iorder(nfree) = i
              if (Abs (neggi) > zero) then
                bnded = .false.
              end if
            end if
          end if
        end do
       !     The indices of the nonzero components of d are now stored
       !       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
       !       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
        if (abs(theta - one) > 1.d-20) then
          !    complete the initialization of p for theta not= one.
          call dscal (col, theta, p(col+1), 1)
        end if
       !     Initialize GCP xcp = x.
       !
        call dcopy (n, x, 1, xcp, 1)
       !
        if (nbreak == 0 .and. nfree == n+1) then
          !      is a zero vector, return with the initial xcp as GCP.
          if (iprint > 100) then
             !
10010       format ("Cauchy X =  ",/, (4 x, 1 p, 6(1 x, d11.4)))
            write (iw, 10010) (xcp(i), i=1, n)
          end if
        else
          !     Initialize c = W'(xcp - x) = 0.
          do j = 1, col2
            c(j) = zero
          end do
          !     Initialize derivative f2.
          f2 = -theta * f1
          f2_org = f2
          if (col > 0) then
            call bmv (m, sy, wt, col, p, v, info)
            if (info /= 0) return
            f2 = f2 - ddot (col2, v, 1, p, 1)
          end if
          dtm = -f1 / f2
          tsum = zero
          nint = 1
          if (iprint >= 99) then
            write (iw,*) "There are ", nbreak, "  breakpoints "
          end if
          !     If there are no breakpoints, locate the GCP and return.
          if (nbreak /= 0) then
            nleft = nbreak
            iter = 1
            tj = zero
            do
              !     Find the next smallest breakpoint;
              !       compute dt = t(nleft) - t(nleft + 1).
              !------------------- the beginning of the loop ------------
              tj0 = tj
              if (iter == 1) then
                ! Since we already have the smallest breakpoint we need not do
                ! heapsort yet. Often only one breakpoint is used and the
                ! cost of heapsort is avoided.
                tj = bkmin
                ibp = iorder(ibkmin)
              else
                if (iter == 2) then
                  ! Replace the already used smallest breakpoint with the
                  ! breakpoint numbered nbreak > nlast, before heapsort call.
                  if (ibkmin /= nbreak) then
                    t(ibkmin) = t(nbreak)
                    iorder(ibkmin) = iorder(nbreak)
                  end if
                end if
                call hpsolb (nleft, t, iorder, iter-2)
                tj = t(nleft)
                ibp = iorder(nleft)
              end if
              dt = tj - tj0
              if (dt /= zero .and. iprint >= 100) then
10020           format (/, "Piece    ", i3, " --f1, f2 at start point ", 1 p, &
               & 2(1 x, d11.4))
                write (iw, 10020) nint, f1, f2
10030           format ("Distance to the next break point =  ", 1 p, d11.4)
                write (iw, 10030) dt
10040           format ("Distance to the stationary point =  ", 1 p, d11.4)
                write (iw, 10040) dtm
              end if
                !     If a minimizer is within this interval,
                !       locate the GCP and return.
              if (dtm < dt) go to 1100
                !     Otherwise fix one variable and
                !       reset the corresponding component of d to zero.
              tsum = tsum + dt
              nleft = nleft - 1
              iter = iter + 1
              dibp = d(ibp)
              d(ibp) = zero
              if (dibp > zero) then
                zibp = u(ibp) - x(ibp)
                xcp(ibp) = u(ibp)
                iwhere(ibp) = 2
              else
                zibp = l(ibp) - x(ibp)
                xcp(ibp) = l(ibp)
                iwhere(ibp) = 1
              end if
              if (iprint >= 100) then
                write (iw,*) "Variable  ", ibp, "  is fixed."
              end if
              if (nleft == 0 .and. nbreak == n) exit
                !     Update the derivative information.
              nint = nint + 1
              dibp2 = dibp ** 2
                !     Update f1 and f2.
                !        temporarily set f1 and f2 for col=0.
              f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp
              f2 = f2 - theta * dibp2
                !
              if (col > 0) then
                   !                          update c = c + dt*p.
                call daxpy (col2, dt, p, 1, c, 1)
                ! choose wbp,
                ! the row of W corresponding to the breakpoint encountered.
                pointr = head
                do j = 1, col
                  wbp(j) = wy(ibp, pointr)
                  wbp(col+j) = theta * ws(ibp, pointr)
                  pointr = Mod (pointr, m) + 1
                end do
                   !           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
                call bmv (m, sy, wt, col, wbp, v, info)
                if (info /= 0) return
                wmc = ddot (col2, c, 1, v, 1)
                wmp = ddot (col2, p, 1, v, 1)
                wmw = ddot (col2, wbp, 1, v, 1)
                   !           update p = p - dibp*wbp.
                call daxpy (col2,-dibp, wbp, 1, p, 1)
                   !           complete updating f1 and f2 while col > 0.
                f1 = f1 + dibp * wmc
                f2 = f2 + 2.0d0 * dibp * wmp - dibp2 * wmw
              end if
                !
              f2 = Max (epsmch*f2_org, f2)
              if (nleft > 0) then
                dtm = -f1 / f2
              else
                go to 1000
              end if
            end do
            !        all n variables are fixed,
            !        return with xcp as GCP.
            dtm = dt
            go to 1200
             !                 to repeat the loop for unsearched intervals.
1000        if (bnded) then
              f1 = zero
              f2 = zero
              dtm = zero
            else
              dtm = -f1 / f2
            end if
          end if
          !
          !------------------- the end of the loop -------------------
1100      if (iprint >= 99) then
            write (iw,*)
            write (iw,*) "GCP found in this segment"
10050       format ("Piece    ", i3, " --f1, f2 at start point ", 1 p, 2(1 x, &
           & d11.4))
            write (iw, 10050) nint, f1, f2
            write (iw, 10040) dtm
          end if
          if (dtm <= zero) then
            dtm = zero
          end if
          tsum = tsum + dtm
          !     Move free variables (i.e., the ones w/o breakpoints) and
          !       the variables whose breakpoints haven't been reached.
          call daxpy (n, tsum, d, 1, xcp, 1)
          !     Update c = c + dtm*p = W'(x^c - x)
          !       which will be used in computing r = Z'(B(x^c - x) + g).
1200      if (col > 0) then
            call daxpy (col2, dtm, p, 1, c, 1)
          end if
          if (iprint > 100) then
            write (iw, 10010) (xcp(i), i=1, n)
          end if
          if (iprint >= 99) then
10060       format (/, "---------------- exit CAUCHY----------------------", &
           & /)
            write (iw, 10060)
          end if
        end if
      end if
    end subroutine cauchy
  !
  !====================== The end of cauchy ==============================
  !
    subroutine cmprlb (n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, &
   & col, head, nfree, cnstnd, info)
      implicit none
      logical, intent (in) :: cnstnd
      integer, intent (in) :: col, head, m, n, nfree
      integer, intent (inout) :: info
      double precision, intent (in) :: theta
      integer, dimension (n), intent (in) :: index
      double precision, dimension (4*m), intent (inout) :: wa
      double precision, dimension (n), intent (in) :: g, x, z
      double precision, dimension (n), intent (inout) :: r
      double precision, dimension (m, m), intent (in) :: sy
      double precision, dimension (m, m), intent (inout) :: wt
      double precision, dimension (n, m), intent (in) :: ws, wy
      integer :: i, j, k, pointr
      double precision :: a1, a2
      intrinsic Mod
      if ( .not. cnstnd .and. col > 0) then
        do i = 1, n
          r(i) = -g(i)
        end do
      else
        do i = 1, nfree
          k = Index(i)
          r(i) = -theta * (z(k)-x(k)) - g(k)
        end do
        call bmv (m, sy, wt, col, wa(2*m+1), wa(1), info)
        if (info /= 0) then
          info = -8
          return
        end if
        pointr = head
        do j = 1, col
          a1 = wa(j)
          a2 = theta * wa(col+j)
          do i = 1, nfree
            k = Index(i)
            r(i) = r(i) + wy(k, pointr) * a1 + ws(k, pointr) * a2
          end do
          pointr = Mod (pointr, m) + 1
        end do
      end if
    !
    end subroutine cmprlb
  !
  !======================= The end of cmprlb =============================
  !
    subroutine errclb (n, m, factr, l, u, nbd, task, info, k)
      implicit none
      double precision, parameter :: zero = 0.0d0
      character (len=60), intent (inout) :: task
      integer, intent (in) :: m, n
      integer, intent (inout) :: info
      integer, intent (out) :: k
      double precision, intent (in) :: factr
      integer, dimension (n), intent (in) :: nbd
      double precision, dimension (n), intent (in) :: l, u
      integer :: i
    !
    !     Check the input arguments for errors.
    !
      if (n <= 0) then
        task = "ERROR: N .LE. 0"
      end if
      if (m <= 0) then
        task = "ERROR: M .LE. 0"
      end if
      if (factr < zero) then
        task = "ERROR: FACTR .LT. 0"
      end if
    !
    !     Check the validity of the arrays nbd(i), u(i), and l(i).
    !
      do i = 1, n
        if (nbd(i) < 0 .or. nbd(i) > 3) then
          !                                                   return
          task = "ERROR: INVALID NBD"
          info = -6
          k = i
        end if
        if (nbd(i) == 2) then
          if (l(i) > u(i)) then
             !                                    return
            task = "ERROR: NO FEASIBLE SOLUTION"
            info = -7
            k = i
          end if
        end if
      end do
    !
    end subroutine errclb
  !
  !======================= The end of errclb =============================
    subroutine formk (n, nsub, ind, nenter, ileave, indx2, iupdat, updatd, wn, &
   & wn1, m, ws, wy, sy, theta, col, head, info)
    !
      implicit none
      double precision, parameter :: zero = 0.0d0
      logical, intent (in) :: updatd
      integer, intent (in) :: head, ileave, iupdat, m, n, nenter, nsub
      integer, intent (inout) :: col
      integer, intent (inout) :: info
      double precision, intent (in) :: theta
      integer, dimension (n), intent (in) :: ind, indx2
      double precision, dimension (2*m, 2*m), intent (inout) :: wn, wn1
      double precision, dimension (m, m), intent (in) :: sy
      double precision, dimension (n, m), intent (in) :: ws, wy
      integer :: col2, dbegin, dend, i, ipntr, is, is1, iy, jpntr, js, js1, &
     & jy, k, k1, m2, pbegin, pend, upcl
      double precision :: temp1, temp2, temp3, temp4
      external dcopy
      double precision, external :: ddot
    !
    !     Form the lower triangular part of
    !               WN1 = [Y' ZZ'Y   L_a'+R_z']
    !                     [L_a+R_z   S'AA'S   ]
    !        where L_a is the strictly lower triangular part of S'AA'Y
    !              R_z is the upper triangular part of S'ZZ'Y.
      if (updatd) then
        if (iupdat > m) then
          !                                 shift old part of WN1.
          do jy = 1, m - 1
            js = m + jy
            call dcopy (m-jy, wn1(jy+1, jy+1), 1, wn1(jy, jy), 1)
            call dcopy (m-jy, wn1(js+1, js+1), 1, wn1(js, js), 1)
            call dcopy (m-1, wn1(m+2, jy+1), 1, wn1(m+1, jy), 1)
          end do
        end if
       !          put new rows in blocks (1,1), (2,1) and (2,2).
        pbegin = 1
        pend = nsub
        dbegin = nsub + 1
        dend = n
        iy = col
        is = m + col
        ipntr = head + col - 1
        if (ipntr > m) then
          ipntr = ipntr - m
        end if
        jpntr = head
        do jy = 1, col
          js = m + jy
          temp1 = zero
          temp2 = zero
          temp3 = zero
          !             compute element jy of row 'col' of Y'ZZ'Y
          do k = pbegin, pend
            k1 = ind(k)
            temp1 = temp1 + wy(k1, ipntr) * wy(k1, jpntr)
          end do
          !             compute elements jy of row 'col' of L_a and S'AA'S
          do k = dbegin, dend
            k1 = ind(k)
            temp2 = temp2 + ws(k1, ipntr) * ws(k1, jpntr)
            temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr)
          end do
          wn1(iy, jy) = temp1
          wn1(is, js) = temp2
          wn1(is, jy) = temp3
          jpntr = Mod (jpntr, m) + 1
        end do
       !          put new column in block (2,1).
        jy = col
        jpntr = head + col - 1
        if (jpntr > m) then
          jpntr = jpntr - m
        end if
        ipntr = head
        do i = 1, col
          is = m + i
          temp3 = zero
          !             compute element i of column 'col' of R_z
          do k = pbegin, pend
            k1 = ind(k)
            temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr)
          end do
          ipntr = Mod (ipntr, m) + 1
          wn1(is, jy) = temp3
        end do
        upcl = col - 1
      else
        upcl = col
      end if
    !       modify the old parts in blocks (1,1) and (2,2) due to changes
    !       in the set of free variables.
      ipntr = head
      do iy = 1, upcl
        is = m + iy
        jpntr = head
        do jy = 1, iy
          js = m + jy
          temp1 = zero
          temp2 = zero
          temp3 = zero
          temp4 = zero
          do k = 1, nenter
            k1 = indx2(k)
            temp1 = temp1 + wy(k1, ipntr) * wy(k1, jpntr)
            temp2 = temp2 + ws(k1, ipntr) * ws(k1, jpntr)
          end do
          do k = ileave, n
            k1 = indx2(k)
            temp3 = temp3 + wy(k1, ipntr) * wy(k1, jpntr)
            temp4 = temp4 + ws(k1, ipntr) * ws(k1, jpntr)
          end do
          wn1(iy, jy) = wn1(iy, jy) + temp1 - temp3
          wn1(is, js) = wn1(is, js) - temp2 + temp4
          jpntr = Mod (jpntr, m) + 1
        end do
        ipntr = Mod (ipntr, m) + 1
      end do
    !       modify the old parts in block (2,1).
      ipntr = head
      do is = m + 1, m + upcl
        jpntr = head
        do jy = 1, upcl
          temp1 = zero
          temp3 = zero
          do k = 1, nenter
            k1 = indx2(k)
            temp1 = temp1 + ws(k1, ipntr) * wy(k1, jpntr)
          end do
          do k = ileave, n
            k1 = indx2(k)
            temp3 = temp3 + ws(k1, ipntr) * wy(k1, jpntr)
          end do
          if (is <= jy+m) then
            wn1(is, jy) = wn1(is, jy) + temp1 - temp3
          else
            wn1(is, jy) = wn1(is, jy) - temp1 + temp3
          end if
          jpntr = Mod (jpntr, m) + 1
        end do
        ipntr = Mod (ipntr, m) + 1
      end do
    !     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
    !                                     [-L_a +R_z        S'AA'S*theta]
    !
      m2 = 2 * m
      do iy = 1, col
        is = col + iy
        is1 = m + iy
        do jy = 1, iy
          js = col + jy
          js1 = m + jy
          wn(jy, iy) = wn1(iy, jy) / theta
          wn(js, is) = wn1(is1, js1) * theta
        end do
        do jy = 1, iy - 1
          wn(jy, is) = -wn1(is1, jy)
        end do
        do jy = iy, col
          wn(jy, is) = wn1(is1, jy)
        end do
        wn(iy, iy) = wn(iy, iy) + sy(iy, iy)
      end do
    !
    !     Form the upper triangle of
    !          WN= [  LL'            L^-1(-L_a'+R_z')]
    !              [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
    !
    !        first Cholesky factor (1,1) block of wn to get LL'
    !                          with L' stored in the upper triangle of wn.
      call dpofa (wn, m2, col, info)
      if (info /= 0) then
        info = -1
        return
      end if
    !        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2 * col
      do js = col + 1, col2
        call dtrsl (wn, m2, col, wn(1, js), 11, info)
      end do
    !
    !     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
    !        upper triangle of (2,2) block of wn.
    !
      do is = col + 1, col2
        do js = is, col2
          wn(is, js) = wn(is, js) + ddot (col, wn(1, is), 1, wn(1, js), 1)
        end do
      end do
    !
    !     Cholesky factorization of (2,2) block of wn.
    !
      call dpofa (wn(col+1, col+1), m2, col, info)
    !
      if (info /= 0) then
        info = -2
      end if
    !
    end subroutine formk
  !
  !======================= The end of formk ==============================
  !
    subroutine formt (m, wt, sy, ss, col, theta, info)
      implicit none
      double precision, parameter :: zero = 0.0d0
      integer, intent (in) :: col, m
      integer, intent (inout) :: info
      double precision, intent (in) :: theta
      double precision, dimension (m, m), intent (in) :: ss, sy
      double precision, dimension (m, m), intent (inout) :: wt
      integer :: i, j, k, k1
      double precision :: ddum
      intrinsic Min
    !
    !
    !     Form the upper half of  T = theta*SS + L*D^(-1)*L',
    !        store T in the upper triangle of the array wt.
      do j = 1, col
        wt(1, j) = theta * ss (1, j)
      end do
      do i = 2, col
        do j = i, col
          k1 = Min (i, j) - 1
          ddum = zero
          do k = 1, k1
            ddum = ddum + sy(i, k) * sy(j, k) / sy(k, k)
          end do
          wt(i, j) = ddum + theta * ss (i, j)
        end do
      end do
    !     Cholesky factorize T to J*J' with
    !        J' stored in the upper triangle of wt.
      call dpofa (wt, m, col, info)
    !
      if (info /= 0) then
        info = -3
      end if
    !
    end subroutine formt
  !
  !======================= The end of formt ==============================
    subroutine freev (n, nfree, index, nenter, ileave, indx2, iwhere, wrk, &
   & updatd, cnstnd, iprint, iter)
      use chanel_C, only: iw
    !
      implicit none
      logical, intent (in) :: cnstnd, updatd
      logical, intent (out) :: wrk
      integer, intent (in) :: iprint, iter, n
      integer, intent (inout) :: nfree
      integer, intent (out) :: ileave, nenter
      integer, dimension (n), intent (in) :: iwhere
      integer, dimension (n), intent (inout) :: index
      integer, dimension (n), intent (out) :: indx2
      integer :: i, iact, k
    !
    !
      nenter = 0
      ileave = n + 1
      if (iter > 0 .and. cnstnd) then
       !                           count the entering and leaving variables.
        do i = 1, nfree
          k = Index(i)
          if (iwhere(k) > 0) then
            ileave = ileave - 1
            indx2(ileave) = k
            if (iprint >= 100) then
              write (iw,*) "Variable ", k, " leaves the set of free variables"
            end if
          end if
        end do
        do i = 1 + nfree, n
          k = Index(i)
          if (iwhere(k) <= 0) then
            nenter = nenter + 1
            indx2(nenter) = k
            if (iprint >= 100) then
              write (iw,*) "Variable ", k, " enters the set of free variables"
            end if
          end if
        end do
        if (iprint >= 99) then
          write (iw,*) n + 1 - ileave, " variables leave; ", nenter, &
         & " variables enter"
        end if
      end if
      wrk = (ileave < n+1) .or. (nenter > 0) .or. updatd
    !     Find the index set of free and active variables at the GCP.
      nfree = 0
      iact = n + 1
      do i = 1, n
        if (iwhere(i) <= 0) then
          nfree = nfree + 1
          Index(nfree) = i
        else
          iact = iact - 1
          Index(iact) = i
        end if
      end do
    !
      if (iprint >= 99) then
        write (iw,*) nfree, " variables are free at GCP ", iter + 1
      end if
    !
    end subroutine freev
  !
  !======================= The end of freev ==============================
  !
    subroutine hpsolb (n, t, iorder, iheap)
      implicit none
      integer, intent (in) :: n
      double precision, dimension (n), intent (inout) :: t
      integer, dimension (n), intent (inout) :: iorder
      integer, intent (in) :: iheap
      integer :: i, indxin, indxou, j, k
      double precision :: ddum, out
    !
      if (iheap == 0) then
       !
       !        Rearrange the elements t(1) to t(n) to form a heap.
       !
        do k = 2, n
          ddum = t(k)
          indxin = iorder(k)
          !
          !           Add ddum to the heap.
          i = k
          do while (i >  1)
            j = i / 2
            if (ddum < t(j)) then
              t(i) = t(j)
              iorder(i) = iorder(j)
              i = j
            else
              exit
            end if
          end do
          t(i) = ddum
          iorder(i) = indxin
        end do
      end if
    !     Assign to 'out' the value of t(1), the least member of the heap,
    !        and rearrange the remaining members to form a heap as
    !        elements 1 to n-1 of t.
      if (n <= 1) return
    !
      i = 1
      out = t(1)
      indxou = iorder(1)
      ddum = t(n)
      indxin = iorder(n)
      do
       !
       !        Restore the heap
        j = i + i
        if (j > n-1) exit
        if (t(j+1) < t(j)) then
          j = j + 1
        end if
        if (t(j) < ddum) then
          t(i) = t(j)
          iorder(i) = iorder(j)
          i = j
        else
          exit
        end if
      end do
      t(i) = ddum
      iorder(i) = indxin
    !     Put the least member in t(n).
    !
      t(n) = out
      iorder(n) = indxou
    !
    end subroutine hpsolb
  !
  !====================== The end of hpsolb ==============================
  !
    subroutine lnsrlb (n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
   & stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, task, &
   & boxed, cnstnd, csave, isave, dsave)
    !
      implicit none
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      double precision, parameter :: big = 1.d-5
      double precision, parameter :: ftol = 1.0d-3
      double precision, parameter :: gtol = 0.9d0
      double precision, parameter :: xtol = 0.1d0
      character (len=60), intent (inout) :: csave, task
      logical, intent (in) :: boxed, cnstnd
      integer, intent (in) :: iter, n
      integer, intent (inout) :: ifun, nfgv
      integer, intent (inout) :: info, iback
      double precision, intent (in) :: f
      double precision, intent (inout) :: dnorm, stp, stpmx
      double precision, intent (inout) :: dtd, fold, gd, gdold, xstep
      integer, dimension (2), intent (inout) :: isave
      integer, dimension (n), intent (in) :: nbd
      double precision, dimension (13), intent (inout) :: dsave
      double precision, dimension (n), intent (in) :: l, u
      double precision, dimension (n), intent (inout) :: d, g, r, t, x, z
      integer :: i
      double precision :: a1, a2
      external dcopy
      double precision, external :: ddot
    !
    !
      if (task(1:5) /= "FG_LN") then
       !
        dtd = ddot (n, d, 1, d, 1)
        dnorm = Sqrt (dtd)
       !
       !     Determine the maximum step length.
       !
        stpmx = big
        if (cnstnd) then
          if (iter == 0) then
            stpmx = one
          else
            do i = 1, n
              a1 = d(i)
              if (nbd(i) /= 0) then
                if (a1 < zero .and. nbd(i) <= 2) then
                  a2 = l(i) - x(i)
                  if (a2 >= zero) then
                    stpmx = zero
                  else if (a1*stpmx < a2) then
                    stpmx = a2 / a1
                  end if
                else if (a1 > zero .and. nbd(i) >= 2) then
                  a2 = u(i) - x(i)
                  if (a2 <= zero) then
                    stpmx = zero
                  else if (a1*stpmx > a2) then
                    stpmx = a2 / a1
                  end if
                end if
              end if
            end do
          end if
        end if
        if (iter == 0 .and. .not. boxed) then
          stp = Min (one/dnorm, stpmx)
        else
          stp = one
        end if
       !
        call dcopy (n, x, 1, t, 1)
        call dcopy (n, g, 1, r, 1)
        fold = f
        ifun = 0
        iback = 0
        csave = "START"
      end if
      gd = ddot (n, g, 1, d, 1)
      if (ifun == 0) then
        gdold = gd
        if (gd >= zero) then
          !                               the directional derivative >=0.
          !                               Line search is impossible.
          info = -4
          return
        end if
      end if
    !
      call dcsrch (f, gd, stp, ftol, gtol, xtol, zero, stpmx, csave, isave, &
     & dsave)
    !
      xstep = stp * dnorm
    !
      if (csave(1:4) == "CONV" .or. csave(1:4) == "WARN") then
        task = "NEW_X"
        return
      end if
      task = "FG_LNSRCH"
      ifun = ifun + 1
      nfgv = nfgv + 1
      iback = ifun - 1
      if (Abs(stp - one) < 1.d-20) then
        call dcopy (n, z, 1, x, 1)
      else
        do i = 1, n
          x(i) = stp * d(i) + t(i)
        end do
      end if
    !
    end subroutine lnsrlb
  !
  !======================= The end of lnsrlb =============================
  !
    subroutine matupd (n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head, &
   & theta, rr, dr, stp, dtd)
      implicit none
      double precision, parameter :: one = 1.0d0
      integer, intent (in) :: iupdat, m, n
      integer, intent (inout) :: col, head, itail
      double precision, intent (in) :: dr, dtd, rr, stp
      double precision, intent (out) :: theta
      double precision, dimension (n), intent (inout) :: d, r
      double precision, dimension (m, m), intent (inout) :: ss, sy
      double precision, dimension (n, m), intent (inout) :: ws, wy
      integer :: j, pointr
      external dcopy
      double precision, external :: ddot
      intrinsic Mod
    !
    !
    !     Set pointers for matrices WS and WY.
      if (iupdat <= m) then
        col = iupdat
        itail = Mod (head+iupdat-2, m) + 1
      else
        itail = Mod (itail, m) + 1
        head = Mod (head, m) + 1
      end if
    !     Update matrices WS and WY.
    !
      call dcopy (n, d, 1, ws(1, itail), 1)
      call dcopy (n, r, 1, wy(1, itail), 1)
    !     Set theta=yy/ys.
      theta = rr / dr
    !     Form the middle matrix in B.
    !        update the upper triangle of SS,
    !                                         and the lower triangle of SY:
      if (iupdat > m) then
       !                              move old information
        do j = 1, col - 1
          call dcopy (j, ss(2, j+1), 1, ss(1, j), 1)
          call dcopy (col-j, sy(j+1, j+1), 1, sy(j, j), 1)
        end do
      end if
    !        add new information: the last row of SY
    !                                             and the last column of SS:
      pointr = head
      do j = 1, col - 1
        sy(col, j) = ddot (n, d, 1, wy(1, pointr), 1)
        ss (j, col) = ddot (n, ws(1, pointr), 1, d, 1)
        pointr = Mod (pointr, m) + 1
      end do
      if (Abs(stp - one) < 1.d-20) then
        ss (col, col) = dtd
      else
        ss (col, col) = stp * stp * dtd
      end if
      sy(col, col) = dr
    !
    end subroutine matupd
  !
  !======================= The end of matupd =============================
  !
    subroutine prn1lb (n, m, l, u, x, iprint, itfile, epsmch)

      use chanel_C, only: iw
      implicit none
      integer, intent (in) :: iprint, itfile, m, n
      double precision, intent (in) :: epsmch
      double precision, dimension (n), intent (in) :: l, u, x
      integer :: i
    !
    !
      if (iprint < 0) return
10000 format ("RUNNING THE L-BFGS-B CODE",/,/, "           * * *",/,/, &
     & "Machine precision =", 1 p, d10.3)
      write (iw, 10000) epsmch
      write (iw,*) "N = ", n, "    M = ", m
      if (iprint >= 1) then
10010   format ("RUNNING THE L-BFGS-B CODE",/,/, &
       & "it    = iteration number"/ &
       & "nf    = number of function evaluations",/, &
       & "nint  = number of segments explored during the Cauchy search",/, &
       & "nact  = number of active bounds at the generalized Cauchy point",/, &
       & "sub   = manner in which the subspace minimization terminated:",/, &
       & "        con = converged, bnd = a bound was reached",/, &
       & "itls  = number of iterations performed in the line search",/, &
       & "stepl = step length used",/, &
       & "tstep = norm of the displacement (total step)",/, &
       & "projg = norm of the projected gradient",/, &
       & "f     = function value",/,/, "           * * *",/,/, &
       & "Machine precision =", 1 p, d10.3)
        write (itfile, 10010) epsmch
        write (itfile,*) "N = ", n, "    M = ", m
10020   format (/, 3 x, "it", 3 x, "nf", 2 x, "nint", 2 x, "nact", 2 x, "sub", &
       & 2 x, "itls", 2 x, "stepl", 4 x, "tstep", 5 x, "projg", 8 x, "f")
        write (itfile, 10020)
        if (iprint > 100) then
          !
10030     format (/, a4, 1 p, 6(1 x, d11.4), /, (4 x, 1 p, 6(1 x, d11.4)))
          write (iw, 10030) "L =", (l(i), i=1, n)
          write (iw, 10030) "X0 =", (x(i), i=1, n)
          write (iw, 10030) "U =", (u(i), i=1, n)
        end if
      end if
    !
    end subroutine prn1lb
  !
  !======================= The end of prn1lb =============================
  !
    subroutine prn2lb (n, x, f, g, iprint, itfile, iter, nfgv, nact, sbgnrm, &
   & nint, word, iword, iback, stp, xstep)

      use chanel_C, only: iw
      implicit none

      character (len=3), intent (out) :: word
      integer, intent (in) :: iback, iprint, iter, itfile, iword, n, nact, &
     & nfgv, nint
      double precision, intent (in) :: f, sbgnrm, stp, xstep
      double precision, dimension (n), intent (in) :: g, x
      integer :: i, imod
    !
    !           'word' records the status of subspace solutions.
      if (iword == 0) then
       !                            the subspace minimization converged.
        word = "con"
      else if (iword == 1) then
       !                          the subspace minimization stopped at a bound.
        word = "bnd"
      else if (iword == 5) then
       !                             the truncated Newton step has been used.
        word = "TNT"
      else
        word = "---"
      end if
      if (iprint >= 99) then
        write (iw,*) "LINE SEARCH", iback, " times; norm of step = ", xstep
10000   format (/, "At iterate", i5, 4 x, "f= ", 1 p, d12.5, 4 x, "|proj g|= ",&
       &  1 p, d12.5)
        write (iw, 10000) iter, f, sbgnrm
        if (iprint > 100) then
          !
10010     format (/, a4, 1 p, 6(1 x, d11.4), /, (4 x, 1 p, 6(1 x, d11.4)))
          write (iw, 10010) "X =", (x(i), i=1, n)
          write (iw, 10010) "G =", (g(i), i=1, n)
        end if
      else if (iprint > 0) then
        imod = Mod (iter, iprint)
        if (imod == 0) then
          write (iw, 10000) iter, f, sbgnrm
        end if
      end if
      if (iprint >= 1) then
10020   format (2(1 x, i4), 2(1 x, i5), 2 x, a3, 1 x, i4, 1 p, 2(2 x, d8.1), &
       & 1 p, 2(1 x, d10.3))
        write (itfile, 10020) iter, nfgv, nint, nact, word, iback, stp, xstep, &
       & sbgnrm, f
      end if
    !
    end subroutine prn2lb
  !
  !======================= The end of prn2lb =============================
  !
    subroutine prn3lb (n, x, f, task, iprint, info, itfile, iter, nfgv, &
   & nintol, nskip, nact, sbgnrm, time, nint, word, iback, stp, xstep, k, &
   & cachyt, sbtime, lnscht)

      use chanel_C, only: iw
      implicit none

      character (len=3), intent (in) :: word
      character (len=60), intent (in) :: task
      integer, intent (in) :: iback, info, iprint, iter, itfile, k, n, nact, &
     & nfgv, nint, nintol, nskip
      double precision, intent (in) :: cachyt, f, lnscht, sbgnrm, sbtime, stp, &
     & time, xstep
      double precision, dimension (n), intent (in) :: x
      integer :: i
    !
    !
      if (task(1:5) /= "ERROR") then
       !
        if (iprint >= 0) then
10000     format (/, "           * * *",/,/, &
         & "Tit   = total number of iterations",/, &
         & "Tnf   = total number of function evaluations",/, &
         & "Tnint = total number of segments explored during Cauchy searches",/, &
         & "Skip  = number of BFGS updates skipped",/, &
         & "Nact  = number of active bounds at final generalized Cauchy point",/, &
         & "Projg = norm of the final projected gradient",/, &
         & "F     = final function value",/,/, &
         & "           * * *")
          write (iw, 10000)
10010     format (/, 3 x, "N", 3 x, "Tit", 2 x, "Tnf", 2 x, "Tnint", 2 x, "Skip", &
         & 2 x, "Nact", 5 x, "Projg", 8 x, "F")
          write (iw, 10010)
10020     format (i5, 2(1 x, i4), (1 x, i6), (2 x, i4), (1 x, i5), 1 p, 2(2 x, &
         & d10.3))
          write (iw, 10020) n, iter, nfgv, nintol, nskip, nact, sbgnrm, f
          if (iprint >= 100) then
             !
10030       format (/, a4, 1 p, 6(1 x, d11.4), /, (4 x, 1 p, 6(1 x, d11.4)))
            write (iw, 10030) "X =", (x(i), i=1, n)
          end if
          if (iprint >= 1) then
            write (iw,*) " F =", f
          end if
        end if
      end if
      if (iprint < 0) return
10040 format (/, a60)
      write (iw, 10040) task
      if (info /= 0) then
        if (info ==-1) then
10050     format (/, " Matrix in 1st Cholesky factorization in formk is not Pos. Def.")
          write (iw, 10050)
        end if
        if (info ==-2) then
10060     format (/, " Matrix in 2st Cholesky factorization in formk is not Pos. Def.")
          write (iw, 10060)
        end if
        if (info ==-3) then
10070     format (/, " Matrix in the Cholesky factorization in formt is not Pos. Def.")
          write (iw, 10070)
        end if
        if (info ==-4) then
10080     format (/, " Derivative >= 0, backtracking line search impossible.",/, &
         & "   Previous x, f and g restored.",/, &
         & " Possible causes: 1 error in function or gradient evaluation;",/, &
         & "                  2 rounding errors dominate computation.")
          write (iw, 10080)
        end if
        if (info ==-5) then
10090     format (/, " Warning:  more than 10 function and gradient",/, &
         & "   evaluations in the last line search.  Termination",/, &
         & "   may possibly be caused by a bad search direction.")
          write (iw, 10090)
        end if
        if (info ==-6) then
          write (iw,*) " Input nbd(", k, ") is invalid."
        end if
        if (info ==-7) then
          write (iw,*) " l(", k, ") > u(", k, ").  No feasible solution."
        end if
        if (info ==-8) then
10100     format (/, " The triangular system is singular.")
          write (iw, 10100)
        end if
        if (info ==-9) then
10110     format (/, " Line search cannot locate an adequate point after 20 function",/, &
         & "  and gradient evaluations.  Previous x, f and g restored.",/, &
         & " Possible causes: 1 error in function or gradient evaluation;",/, &
         & "                  2 rounding error dominate computation.")
          write (iw, 10110)
        end if
      end if
      if (iprint >= 1) then
10120   format (/, " Cauchy                time", 1 p, e10.3, " seconds.",/ &
       & " Subspace minimization time", 1 p, e10.3, " seconds.",/ &
       & " Line search           time", 1 p, e10.3, " seconds.")
        write (iw, 10120) cachyt, sbtime, lnscht
      end if
10130 format (/, " Total User time", 1 p, e10.3, " seconds.", /)
      write (iw, 10130) time
      if (iprint < 1) return
      if (info ==-4 .or. info ==-9) then
10140   format (2(1 x, i4), 2(1 x, i5), 2 x, a3, 1 x, i4, 1 p, 2(2 x, d8.1), &
       & 6 x, "-", 10 x, "-")
        write (itfile, 10140) iter, nfgv, nint, nact, word, iback, stp, xstep
      end if
      write (itfile, 10040) task
      if (info /= 0) then
        if (info ==-1) then
          write (itfile, 10050)
        end if
        if (info ==-2) then
          write (itfile, 10060)
        end if
        if (info ==-3) then
          write (itfile, 10070)
        end if
        if (info ==-4) then
          write (itfile, 10080)
        end if
        if (info ==-5) then
          write (itfile, 10090)
        end if
        if (info ==-8) then
          write (itfile, 10100)
        end if
        if (info ==-9) then
          write (itfile, 10110)
        end if
      end if
      write (itfile, 10130) time
    !
    end subroutine prn3lb
  !
  !======================= The end of prn3lb =============================
  !
    subroutine projgr (n, l, u, nbd, x, g, sbgnrm)
    !
      implicit none
      double precision, parameter :: zero = 0.0d0
      integer, intent (in) :: n
      double precision, intent (out) :: sbgnrm
      integer, dimension (n), intent (in) :: nbd
      double precision, dimension (n), intent (in) :: g, l, u, x
      integer :: i
      double precision :: gi
      intrinsic Abs, Max, Min
    !
    !
      sbgnrm = zero
      do i = 1, n
        gi = g(i)
        if (nbd(i) /= 0) then
          if (gi < zero) then
            if (nbd(i) >= 2) then
              gi = Max ((x(i)-u(i)), gi)
            end if
          else if (nbd(i) <= 2) then
            gi = Min ((x(i)-l(i)), gi)
          end if
        end if
        sbgnrm = Max (sbgnrm, Abs (gi))
      end do
    !
    end subroutine projgr
  !
  !======================= The end of projgr =============================
  !
    subroutine subsm (n, m, nsub, ind, l, u, nbd, x, d, ws, wy, theta, col, &
   & head, iword, wv, wn, iprint, info)

      use chanel_C, only: iw
      implicit none
      double precision, parameter :: one = 1.0d0
      double precision, parameter :: zero = 0.0d0
      integer, intent (in) :: col, head, iprint, m, n, nsub
      integer, intent (inout) :: info
      integer, intent (out) :: iword
      double precision, intent (in) :: theta
      integer, dimension (n), intent (in) :: nbd
      integer, dimension (nsub), intent (in) :: ind
      double precision, dimension (2*m), intent (inout) :: wv
      double precision, dimension (n), intent (in) :: l, u
      double precision, dimension (n), intent (inout) :: d, x
      double precision, dimension (2*m, 2*m), intent (inout) :: wn
      double precision, dimension (n, m), intent (in) :: ws, wy
      integer :: col2, i, ibd, j, js, jy, k, m2, pointr
      double precision :: alpha, dk, temp1, temp2
    !
      ibd = 0
      if (nsub <= 0) return
      if (iprint >= 99) then
       !
10000   format (/, "----------------SUBSM entered-----------------", /)
        write (iw, 10000)
      end if
    !
    !     Compute wv = W'Zd.
    !
      pointr = head
      do i = 1, col
        temp1 = zero
        temp2 = zero
        do j = 1, nsub
          k = ind(j)
          temp1 = temp1 + wy(k, pointr) * d(j)
          temp2 = temp2 + ws(k, pointr) * d(j)
        end do
        wv(i) = temp1
        wv(col+i) = theta * temp2
        pointr = Mod (pointr, m) + 1
      end do
    !     Compute wv:=K^(-1)wv.
    !
      m2 = 2 * m
      col2 = 2 * col
      call dtrsl (wn, m2, col2, wv, 11, info)
      if (info /= 0) return
      do i = 1, col
        wv(i) = -wv(i)
      end do
      call dtrsl (wn, m2, col2, wv, 1, info)
      if (info /= 0) return
    !     Compute d = (1/theta)d + (1/theta**2)Z'W wv.
      pointr = head
      do jy = 1, col
        js = col + jy
        do i = 1, nsub
          k = ind(i)
          d(i) = d(i) + wy(k, pointr) * wv(jy) / theta + ws(k, pointr) * &
         & wv(js)
        end do
        pointr = Mod (pointr, m) + 1
      end do
      do i = 1, nsub
        d(i) = d(i) / theta
      end do
    !     Backtrack to the feasible region.
      alpha = one
      temp1 = alpha
      do i = 1, nsub
        k = ind(i)
        dk = d(i)
        if (nbd(k) /= 0) then
          if (dk < zero .and. nbd(k) <= 2) then
            temp2 = l(k) - x(k)
            if (temp2 >= zero) then
              temp1 = zero
            else if (dk*alpha < temp2) then
              temp1 = temp2 / dk
            end if
          else if (dk > zero .and. nbd(k) >= 2) then
            temp2 = u(k) - x(k)
            if (temp2 <= zero) then
              temp1 = zero
            else if (dk*alpha > temp2) then
              temp1 = temp2 / dk
            end if
          end if
          if (temp1 < alpha) then
            alpha = temp1
            ibd = i
          end if
        end if
      end do
      if (alpha < one) then
        dk = d(ibd)
        k = ind(ibd)
        if (dk > zero) then
          x(k) = u(k)
          d(ibd) = zero
        else if (dk < zero) then
          x(k) = l(k)
          d(ibd) = zero
        end if
      end if
      do i = 1, nsub
        k = ind(i)
        x(k) = x(k) + alpha * d(i)
      end do
      if (iprint >= 99) then
        if (alpha < one) then
10010     format ("ALPHA = ", f8.5, " backtrack to the BOX")
          write (iw, 10010) alpha
        else
          write (iw,*) "SM solution inside the box"
        end if
        if (iprint > 100) then
10020     format ("Subspace solution X =  ",/, (4 x, 1 p, 6(1 x, d11.4)))
          write (iw, 10020) (x(i), i=1, n)
        end if
      end if
      if (alpha < one) then
        iword = 1
      else
        iword = 0
      end if
      if (iprint >= 99) then
10030   format (/, "----------------exit SUBSM --------------------", /)
        write (iw, 10030)
      end if
    !
    end subroutine subsm
  !====================== The end of subsm ===============================
  !
    subroutine dcsrch (f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, &
   & isave, dsave)
      implicit none
      double precision, parameter :: zero = 0.0d0
      double precision, parameter :: p5 = 0.5d0
      double precision, parameter :: p66 = 0.66d0
      double precision, parameter :: xtrapl = 1.1d0
      double precision, parameter :: xtrapu = 4.0d0
      character (len=*), intent (inout) :: task
      double precision, intent (in) :: f, ftol, g, gtol, stpmax, stpmin, xtol
      double precision, intent (inout) :: stp
      integer, dimension (2), intent (inout) :: isave
      double precision, dimension (13), intent (inout) :: dsave
      logical :: brackt
      integer :: stage
      double precision :: finit, fm, ftest, fx, fxm, fy, fym, ginit, gm, &
     & gtest, gx, gxm, gy, gym, stmax, stmin, stx, sty, width, width1
      intrinsic Abs, Max, Min
    !
    !     Initialization block.
    !
      if (task(1:5) == "START") then
       !
       !        Check the input arguments for errors.
       !
        if (stp < stpmin) then
          task = "ERROR: STP .LT. STPMIN"
        end if !
        if (stp > stpmax) then
          task = "ERROR: STP .GT. STPMAX"
        end if
        if (g >= zero) then
          task = "ERROR: INITIAL G .GE. ZERO"
        end if
        if (ftol < zero) then
          task = "ERROR: FTOL .LT. ZERO"
        end if
        if (gtol < zero) then
          task = "ERROR: GTOL .LT. ZERO"
        end if
        if (xtol < zero) then
          task = "ERROR: XTOL .LT. ZERO"
        end if
        if (stpmin < zero) then
          task = "ERROR: STPMIN .LT. ZERO"
        end if
        if (stpmax < stpmin) then
          task = "ERROR: STPMAX .LT. STPMIN"
        end if
       !
       !        Exit if there are errors on input.
       !
        if (task(1:5) == "ERROR") return
       !
       !        Initialize local variables.
       !
        brackt = .false.
        stage = 1
        finit = f
        ginit = g
        gtest = ftol * ginit
        width = stpmax - stpmin
        width1 = width / p5
       !
       !        The variables stx, fx, gx contain the values of the step,
       !        function, and derivative at the best step.
       !        The variables sty, fy, gy contain the value of the step,
       !        function, and derivative at sty.
       !        The variables stp, f, g contain the values of the step,
       !        function, and derivative at stp.
       !
        stx = 1.d-5
        fx = finit
        gx = ginit
        sty = 1.d-4
        fy = finit
        gy = ginit
        stmin = zero
        stmax = stp + xtrapu * stp
       !
        task = "FG"
      else
       !
       !        Restore local variables.
       !
        if (isave(1) == 1) then
          brackt = .true.
        else
          brackt = .false.
        end if
        stage = isave(2)
        ginit = dsave(1)
        gtest = dsave(2)
        gx = dsave(3)
        gy = dsave(4)
        finit = dsave(5)
        fx = dsave(6)
        fy = dsave(7)
        stx = dsave(8)
        sty = dsave(9)
        stmin = dsave(10)
        stmax = dsave(11)
        width = dsave(12)
        width1 = dsave(13)
       !
       !     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
       !     algorithm enters the second stage.
       !
        ftest = finit + stp * gtest
        if (stage == 1 .and. f <= ftest .and. g >= zero) then
          stage = 2
        end if
       !
       !     Test for warnings.
       !
        if (brackt .and. (stp <= stmin .or. stp >= stmax)) then
          task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS"
        end if
        if (brackt .and. stmax-stmin <= xtol*stmax) then
          task = "WARNING: XTOL TEST SATISFIED"
        end if
        if (abs(stp - stpmax) < 1.d-20 .and. f <= ftest .and. g <= gtest) then
          task = "WARNING: STP = STPMAX"
        end if
        if (Abs(stp - stpmin) < 1.d-20 .and. (f > ftest .or. g >= gtest)) then
          task = "WARNING: STP = STPMIN"
        end if
       !
       !     Test for convergence.
       !
        if (f <= ftest .and. Abs (g) <= gtol*(-ginit)) then
          task = "CONVERGENCE"
        end if
       !
       !     Test for termination.
       !
        if (task(1:4) /= "WARN" .and. task(1:4) /= "CONV") then
          !
          !     A modified function is used to predict the step during the
          !     first stage if a lower function value has been obtained but
          !     the decrease is not sufficient.
          !
          if (stage == 1 .and. f <= fx .and. f > ftest) then
             !
             !        Define the modified function and derivative values.
             !
            fm = f - stp * gtest !
            fxm = fx - stx * gtest
            fym = fy - sty * gtest
            gm = g - gtest
            gxm = gx - gtest
            gym = gy - gtest
             !
             !        Call dcstep to update stx, sty, and to compute the new step.
             !
            call dcstep (stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, &
           & stmin, stmax)
             !
             !        Reset the function and derivative values for f.
             !
            fx = fxm + stx * gtest
            fy = fym + sty * gtest
            gx = gxm + gtest
            gy = gym + gtest
          else
            !
            !  Call dcstep to update stx, sty, and to compute the new step.
            !
            call dcstep (stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, &
           & stmax)
          end if
          !
          !     Decide if a bisection step is needed.
          !
          if (brackt) then
            if (Abs (sty-stx) >= p66*width1) then
              stp = stx + p5 * (sty-stx)
            end if
            width1 = width
            width = Abs (sty-stx)
          end if
          !
          !     Set the minimum and maximum steps allowed for stp.
          !
          if (brackt) then
            stmin = Min (stx, sty)
            stmax = Max (stx, sty)
          else
            stmin = stp + xtrapl * (stp-stx)
            stmax = stp + xtrapu * (stp-stx)
          end if
          !     Force the step to be within the bounds stpmax and stpmin.
          stp = Max (stp, stpmin)
          stp = Min (stp, stpmax)
          !
          !     If further progress is not possible, let stp be the best
          !     point obtained during the search.
          !
          if (brackt .and. (stp <= stmin .or. stp >= stmax) .or. (brackt .and. &
         & stmax-stmin <= xtol*stmax)) then
            stp = stx
          end if
          !
          !     Obtain another function and derivative.
          !
          task = "FG"
        end if
      end if
    !
    !     Save local variables.
    !
    !
      if (brackt) then
        isave(1) = 1
      else
        isave(1) = 0
      end if
      isave(2) = stage
      dsave(1) = ginit
      dsave(2) = gtest
      dsave(3) = gx
      dsave(4) = gy
      dsave(5) = finit
      dsave(6) = fx
      dsave(7) = fy
      dsave(8) = stx
      dsave(9) = sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1
    !
    end subroutine dcsrch
  !====================== The end of dcsrch ==============================
  !
    subroutine dcstep (stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, &
   & stpmax)
      implicit none
      double precision, parameter :: zero = 0.0d0
      double precision, parameter :: p66 = 0.66d0
      double precision, parameter :: two = 2.0d0
      double precision, parameter :: three = 3.0d0
      logical, intent (inout) :: brackt
      double precision, intent (in) :: dp, fp, stpmax, stpmin
      double precision, intent (inout) :: dx, dy, fx, fy, stp, stx, sty
      double precision :: gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta
    !
      sgnd = dp * (dx/Abs(dx))
!
!  Consider pathologic condition: stp = stx
!
      if (Abs (stp - stx) < 1.d-5) stp = stp + 1.d-5
    !
    !     First case: A higher function value. The minimum is bracketed.
    !     If the cubic step is closer to stx than the quadratic step, the
    !     cubic step is taken, otherwise the average of the cubic and
    !     quadratic steps is taken.
    !
      if (fp > fx) then
        theta = three * (fx-fp) / (stp-stx) + dx + dp
        s = Max (Abs (theta), Abs (dx), Abs (dp))
        gamma = s * Sqrt ((theta/s)**2-(dx/s)*(dp/s))
        if (stp < stx) then
          gamma = -gamma
        end if
        p = (gamma-dx) + theta
        q = ((gamma-dx)+gamma) + dp
        r = p / q
        stpc = stx + r * (stp-stx)
        stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/two) * (stp-stx)
        if (Abs (stpc-stx) < Abs (stpq-stx)) then
          stpf = stpc
        else
          stpf = stpc + (stpq-stpc) / two
        end if
        brackt = .true.
       !
       !     Second case: A lower function value and derivatives of opposite
       !     sign. The minimum is bracketed. If the cubic step is farther from
       !     stp than the secant step, the cubic step is taken, otherwise the
       !     secant step is taken.
       !
      else if (sgnd < zero) then
        theta = three * (fx-fp) / (stp-stx) + dx + dp
        s = Max (Abs (theta), Abs (dx), Abs (dp))
        gamma = s * Sqrt ((theta/s)**2-(dx/s)*(dp/s))
        if (stp > stx) then
          gamma = -gamma
        end if
        p = (gamma-dp) + theta
        q = ((gamma-dp)+gamma) + dx
        r = p / q
        stpc = stp + r * (stx-stp)
        stpq = stp + (dp/(dp-dx)) * (stx-stp)
        if (Abs (stpc-stp) > Abs (stpq-stp)) then
          stpf = stpc
        else
          stpf = stpq
        end if
        brackt = .true.
       !
       !     Third case: A lower function value, derivatives of the same sign,
       !     and the magnitude of the derivative decreases.
       !
      else if (Abs (dp) < Abs (dx)) then
       !
       !        The cubic step is computed only if the cubic tends to infinity
       !        in the direction of the step or if the minimum of the cubic
       !        is beyond stp. Otherwise the cubic step is defined to be the
       !        secant step.
       !
        theta = three * (fx-fp) / (stp-stx) + dx + dp !
       !     Fourth case: A lower function value, derivatives of the
       !     same sign, and the magnitude of the derivative does not
       !     decrease. If the minimum is not bracketed, the step is either
       !     stpmin or stpmax, otherwise the cubic step is taken.
       !
        s = Max (Abs (theta), Abs (dx), Abs (dp))
       !
       !        The case gamma = 0 only arises if the cubic does not tend
       !        to infinity in the direction of the step.
       !
        gamma = s * Sqrt (Max(zero, (theta/s)**2-(dx/s)*(dp/s)))
        if (stp > stx) then
          gamma = -gamma
        end if
        p = (gamma-dp) + theta
        q = (gamma+ (dx-dp)) + gamma
        r = p / q
        if (r < zero .and. gamma /= zero) then
          stpc = stp + r * (stx-stp)
        else if (stp > stx) then
          stpc = stpmax
        else
          stpc = stpmin
        end if
        stpq = stp + (dp/(dp-dx)) * (stx-stp)
       !
        if (brackt) then
          !
          !           A minimizer has been bracketed. If the cubic step is
          !           closer to stp than the secant step, the cubic step is
          !           taken, otherwise the secant step is taken.
          !
          if (Abs (stpc-stp) < Abs (stpq-stp)) then
            stpf = stpc
          else
            stpf = stpq
          end if
          if (stp > stx) then
            stpf = Min (stp+p66*(sty-stp), stpf)
          else
            stpf = Max (stp+p66*(sty-stp), stpf)
          end if
        else
          !
          !           A minimizer has not been bracketed. If the cubic step is
          !           farther from stp than the secant step, the cubic step is
          !           taken, otherwise the secant step is taken.
          !
          if (Abs (stpc-stp) > Abs (stpq-stp)) then
            stpf = stpc
          else
            stpf = stpq
          end if
          stpf = Min (stpmax, stpf)
          stpf = Max (stpmin, stpf)
        end if
      else if (brackt) then
        theta = three * (fp-fy) / (sty-stp) + dy + dp
        s = Max (Abs (theta), Abs (dy), Abs (dp))
        gamma = s * Sqrt ((theta/s)**2-(dy/s)*(dp/s))
        if (stp > sty) then
          gamma = -gamma
        end if
        p = (gamma-dp) + theta
        q = ((gamma-dp)+gamma) + dy
        r = p / q
        stpc = stp + r * (sty-stp)
        stpf = stpc
      else if (stp > stx) then
        stpf = stpmax
      else
        stpf = stpmin
      end if
    !
    !     Update the interval which contains a minimizer.
    !
      if (fp > fx) then
        sty = stp
        fy = fp
        dy = dp
      else
        if (sgnd < zero) then
          sty = stx
          fy = fx
          dy = dx
        end if
        stx = stp
        fx = fp
        dx = dp
      end if
    !
    !     Compute the new step.
    !
      stp = stpf
    !
    end subroutine dcstep
  !====================== The end of dcstep ==============================
    double precision function dpmeps ()
    !     **********
    !
    !     Function dpeps
    !
    !     This function computes the machine precision parameter
    !     dpmeps as the smallest floating point number such that
    !     1 + dpmeps differs from 1.
    !
    !     Replaced with F90 function Epsilon
    !     Original fails on COMPAQ Alpha
    !
    !.. Implicit Declarations ..
      implicit none
    !
      dpmeps = Epsilon (1.d0)
    !
    end function dpmeps
  !====================== The end of dpmeps ==============================
    subroutine dpofa (a, lda, n, info)
      implicit none
      integer, intent (in) :: lda
      double precision, dimension (lda,*), intent (inout) :: a
      integer, intent (in) :: n
      integer, intent (inout) :: info
       integer :: j, jm1, k
      double precision :: s, t
      double precision, external :: ddot
    !     begin block with ...exits to 40
    !
    !
      do j = 1, n
        info = j
        s = 0.0d0
        jm1 = j - 1
        if (jm1 >= 1) then
          do k = 1, jm1
            t = a(k, j) - ddot (k-1, a(1, k), 1, a(1, j), 1)
            t = t / a(k, k)
            a(k, j) = t
            s = s + t * t
          end do
        end if
        s = a(j, j) - s
       !     ......exit
        if (s <= 0.0d0) return
        a(j, j) = Sqrt (s)
      end do
      info = 0
    end subroutine dpofa
  !====================== The end of dpofa ===============================
    subroutine dtrsl (t, ldt, n, b, job, info)
      implicit none
      integer, intent (in) :: ldt
      double precision, dimension (ldt,*), intent (inout) :: t
      integer, intent (in) :: n
      double precision, dimension (*), intent (inout) :: b
      integer, intent (in) :: job
      integer, intent (inout) :: info
      integer :: case, j, jj
      double precision :: temp
      external daxpy
      double precision, external :: ddot
    !
    !
    !     begin block permitting ...exits to 150
    !
    !        check for zero diagonal elements.
    !
      do info = 1, n
       !     ......exit
        if (t(info, info) == 0.0d0) return
      end do
      info = 0
    !
    !        determine the task and go to it.
    !
    case = 1
      if (Mod(job, 10) /= 0) then
      case = 2
      end if
      if (Mod(job, 100)/10 /= 0) then
      case = case + 2
      end if
      select case (case)
      case (2)
       !
       !        solve t*x=b for t upper triangular.
       !
        b(n) = b(n) / t(n, n)
        if (n >= 2) then
          do jj = 2, n
            j = n - jj + 1
            temp = -b(j+1)
            call daxpy (j, temp, t(1, j+1), 1, b(1), 1)
            b(j) = b(j) / t(j, j)
          end do
        end if
      case (3)
       !
       !        solve trans(t)*x=b for t lower triangular.
       !
        b(n) = b(n) / t(n, n)
        if (n >= 2) then
          do jj = 2, n
            j = n - jj + 1
            b(j) = b(j) - ddot (jj-1, t(j+1, j), 1, b(j+1), 1)
            b(j) = b(j) / t(j, j)
          end do
        end if
      case (4)
       !
       !        solve trans(t)*x=b for t upper triangular.
       !
        b(1) = b(1) / t(1, 1)
        if (n >= 2) then
          do j = 2, n
            b(j) = b(j) - ddot (j-1, t(1, j), 1, b(1), 1)
            b(j) = b(j) / t(j, j)
          end do
        end if
      case default
       !
       !        solve t*x=b for t lower triangular
       !
        b(1) = b(1) / t(1, 1)
        if (n >= 2) then
          do j = 2, n
            temp = -b(j-1)
            call daxpy (n-j+1, temp, t(j, j-1), 1, b(j), 1)
            b(j) = b(j) / t(j, j)
          end do
        end if
      end select
    end subroutine dtrsl
  !====================== The end of dtrsl ===============================
