subroutine direct(cycle_no)
! 
    use param_global_C, only : valvar, numvar, power, ifiles_8, fnsnew, &
    error, factor, tfns, fns, diffns, nfns, ihrefs, refgeo, is_a_ref, &
    refcer, reftot, refger, refher, refder, wthof, wtdip, wtips, wtgeo, wtpKa, &
    contrl, nmols, locvar, molele, maxmol, maxfns, molnum, molnam, nref, refdir, &
    heats, refhof
!
    use molkst_C, only : koment, title,   tleft, gnorm, last, nvar,  &
    moperr,  atheat, norbs, is_PARAM, id, rjkab1, keywrd, line, time0
    use cosmo_C, only : iseps, useps
!
    use chanel_C, only : iw, ir
    use common_arrays_C, only : nat, nfirst, nlast, uspd, xparam, c, cb, loc, &
    geometry => geo, na, nb, nc, coord
!-----------------------------------------------------------------------
    implicit none
    integer :: cycle_no
    logical :: geopt, full, lderivs, lfact, lfast, used, exists
    integer :: i, itime, j, loop, loopar, ndmols, nel, nerr, nhmols, nimols, &
    & ngmols, big_loop, k
    double precision :: aveder, aveher, aveier, delta, erdmod, erhmod, erimod, &
   & rmsder, rmsher, rmsier, sumder, sumerr, sumger, sumher, sumier, lim_gnorm, &
   & ergmod, rmsger, sum, sumpKa, phase = 1.d0, scale_lim, store
    character (len=20), dimension (300) :: ttype
    integer, dimension (maxmol) :: loch
    double precision, dimension (numvar) :: dstep
    double precision, dimension (300) :: differ, errors, yparam
    logical, dimension (maxfns) :: lfact_here
    double precision, external :: seconds, reada
  !
  !.. Intrinsic Functions ..
    intrinsic Abs, Index, Int
  !
  ! ... Executable Statements ...
  !
  !***********************************************************************
  !
  !   PARAMETRIZATION OPTIMISATION PROCEDURE.
  !
  !   DIRECT IS AN ITERATIVE METHOD.  ON EACH CYCLE THE FUNCTION
  !   SSQ (SUM OF SQUARES OF ERRORS) IS CALCULATED, ALONG WITH THE
  !   PRECISE PARTIAL DIFFERENTIALS OF THE ERRORS WITH RESPECT TO ALL
  !   PARAMETERS.
  !
  !   IN A SECOND STEP THE PARAMETERS ARE CHANGED IN "RAPID" SO AS TO
  !   REDUCE THE SSQ.  THE ASSUMPTION IS MADE THAT ONLY ONE MINIMUM
  !   EXISTS, AND THAT THE LOCAL SPACE IS STRICTLY PARABOLIC.
  !
  !***********************************************************************
    geopt = (Index (contrl, "1SCF") == 0)
    lfact = (Index (contrl, " NOSCALE") == 0)  
    full = (Index (contrl, "FULL") /= 0)
    lderivs = (Index (contrl, "DERI") /= 0)
    lfast = (Index (contrl, " NOFINE") /= 0)
    useps = (Index (contrl, " EPS") + Index (contrl, " PKA")/= 0)
    iseps = useps
    lim_gnorm = 3.d0
    i = Index (contrl, "GNORM") 
    if( i /= 0) lim_gnorm = Reada(contrl, i)
    used = .false.
    j = 0
    k = 0
    sumher = 0.d0
    sumder = 0.d0
    sumier = 0.d0
    sumger = 0.d0
    sumpKa = 0.d0
    nhmols = 0
    ndmols = 0
    nimols = 0
    ngmols = 0
    aveher = 0.d0
    aveder = 0.d0
    aveier = 0.d0
    rmsher = 0.d0
    rmsder = 0.d0
    rmsier = 0.d0
    rmsger = 0.d0
    erhmod = 0.d0
    erdmod = 0.d0
    erimod = 0.d0
    ergmod = 0.d0
    nfns = 0
    phase = - phase
    itime = Int (seconds (1))
    if (cycle_no == 14) then
     i = 1
    end if
    write (ifiles_8, "(2A)") &
    &"   Molecule                                       Time  Heat     "&
    &, "Dipole     I.P.     Geometry     Total"
    if(numvar == -10) open (unit=55, file="CALCD.TXT", status="UNKNOWN",iostat = loop)
    do loop = 1, nmols
      molnum = molnum + 1
      tleft = 10000.d0
      time0 = seconds(1)
    !
    !  Restore all information for molecule number LOOP
    !
      call calpar
      call getmol (loop)
      write(iw,*)trim(koment)
      endfile (iw) 
      backspace (iw) 
      deallocate(c)
      if (allocated(cb)) deallocate(cb)
      allocate(c(norbs, norbs), cb(norbs, norbs))
      line = molnam
      call upcase (line, len_trim(line))
      if (line(:4) == "CORE") then
!
!  Update core-core heat
!
        do i = 1, nref
          line = trim(refdir(i))//trim(molnam)//".mop"
          inquire (file=trim(line), exist= exists)
          if (exists) goto 99
        end do
        goto 98
    99  open(ir,file=trim(line), iostat=i)
        if (i /= 0) goto 98 
        j = 0
        do 
          read(ir,'(a)', iostat=i)line
          if (i /= 0) goto  98
          if (line(1:1) /= "*") j = j + 1
          if (j == 3) goto 97
        end do
    97  call upcase (line, len_trim(line))
        i = index(line," H=")
        store = refhof
        refhof = reada(line,i+3)
        close (ir)
        if (Abs(store - refhof) > 0.001d0) then
           write (ifiles_8, "(a,f9.3,a,f9.3,a,a)") " HoF changed from", &
             store, " to", refhof, " kcal/mol for system '"// trim(molnam)//"'"
           heats(loop) = refhof 
        end if
      end if
   98 moperr = .false.
      gnorm = 0
      last = 0
    !
    !  If needed, optimize the geometry
    !
  !  Write(ifiles_8,"(A,F12.4)")" Tleft:",tleft
 !     if (cycle_no == 1) is_PARAM = .false.
      if (iseps) then
         call gmetry (geometry, coord) 
      ! The following routine constructs the dielectric screening surface 
        call cosini(.false.)
      !  if (moperr) goto 100 
        call coscav 
        call mkbmat
      !  if (moperr) go to 100  
      end if     
      if (geopt) then
        call optgeo (xparam, yparam, nvar, refgeo(1), lim_gnorm)
      else
         call savgeo (loop, geometry, na, nb, nc, xparam, loc)
      end if
      moperr = .false.
    !
    !  Compute the value of all errors
    !
      rjkab1 = 0.d0
      call parfg (errors, ttype, nerr, loop, .true.)
      if (geopt) then
!
!  Write out the geometry so it can be used by a later run
!
        call savgeo (loop, geometry, na, nb, nc, xparam, loc)
!
!  If unable to write geometry, reset MOPERR and continue
!
        moperr = .false.
      end if
      is_PARAM = .true.
      endfile (ifiles_8)
      backspace (ifiles_8)
      if( .not. is_a_ref(loop)) then
        if (index(keywrd," XFAC") == 0) then
          if (refher /= 0.d0) then
            aveher = aveher + refher
            rmsher = rmsher + Abs(refher) ** power
            erhmod = erhmod + Abs (refher)
            nhmols = nhmols + 1
          end if
        end if
        if (refcer /= 0.d0) then
          aveier = aveier + refcer
          rmsier = rmsier + Abs(refcer) ** power
          erimod = erimod + Abs (refcer)
          nimols = nimols + 1
        end if
        if (refder /= 0.d0) then
          aveder = aveder + refder
          rmsder = rmsder + Abs(refder) ** power
          erdmod = erdmod + Abs (refder)
          ndmols = ndmols + 1
        end if
        if (refger > 1.d-5) then
          if (wtgeo > 1.d-4) then
            do i = 1, nerr
              ergmod = ergmod + Abs(errors(i)/wtgeo)
              rmsger = rmsger + Abs(errors(i)/wtgeo) ** power
            end do
          else
            do i = 1, nerr
              ergmod = ergmod + Abs(errors(i))
              rmsger = rmsger + Abs(errors(i)) ** power
            end do
          end if
          ngmols = ngmols + nerr
          end if
        sumher = sumher + Abs(refher*wthof) ** power
        sumpKa = sumpKa + Abs(refher*wtpKa) ** power
        sumder = sumder + Abs(refder*wtdip) ** power
        sumier = sumier + Abs(refcer*wtips) ** power
        sumger = sumger + refger
        reftot = Abs(refher*(wthof + wtpKa)) ** power + Abs(refder*wtdip) ** power + &
        & Abs(refcer*wtips) ** power + refger
      else
!
!  Set errors for some functions to zero. (these functions are data
!  that are not used by the optimizer, but whose derivatives are used by other
!
        do i = 1, nerr
          errors(i) = 0.d0
        end do
      end if
    !
    !  Store all the errors for use by the parameter optimizer
    !
      do i = 1, nerr
        fns(i+nfns) = errors(i)
        tfns (i+nfns) = koment (1:41) // " " // ttype (i)
      end do
!
!  Set lfact_here: .true. if scaling is allowed for this system
!
      lfact_here(nfns+1) = ((Index( title, " NOSCALE") == 0) .and. lfact)
      do i = 2,nerr
        lfact_here(nfns+i) = lfact_here(nfns+1)
      end do
      if(numvar /= 0) then
        do big_loop = 1,5
        !
        !   Move all parameters delta in the negative direction.
        !
          do loopar = 1, numvar
            if(lfast) then
    !
    !  Trying different step sizes to find the best for the derivatives.
    !
    !   2.d-3 results in a lot of numerical instability
    !   4.d-3 previous value - used from 1999 - 2005
    !   I've not tried higher values than 8.d-3
    !
              dstep(loopar) = 2.d-3 
            else
    !
    !   1.d-2 looks good
    !
              dstep(loopar) = 10.d-3*phase !min(max(2.d-2, abs(valvar(loopar)*0.02d0)), 3.d-2 ) * phase
            end if
              
            delta = -dstep(loopar)
    !
    ! Single-sided derivatives are only good for quick rough work.  For accurate
    ! work, double-sided derivatives are needed.
    
            if (lfast .and. (id /= 3 .or. locvar(2,loopar) < 200)) &
               call update (locvar(1, loopar), locvar(2, loopar), delta, 1.d0)
          end do
        !
        !  Update derived parameters, USPD, and ATHEAT
        !
          call calpar
          call getusp (nat, nfirst, nlast, uspd, atheat)
        !
        !  Compute the value of all errors at the point -0.5*delta
        !
          call parfg (errors, ttype, nerr, 0, .true.)
        !
        !  For each parameter in turn, step two delta in the positive direction
        !
    
          do loopar = 1, numvar
            nel = molele(1, loop)
            used = .true.
            do i = 2, nel
              if (mod(locvar(2, loopar),200) == molele(i, loop))  then
                k = 0
                if (locvar(2,loopar) > 200) then
                  k = 1
                  if (molele(i, loop) /= 99) then
                    do j = 2, nel
                      if (locvar(2, loopar)/200 == molele(j, loop) .and. molele(j, loop) /= 99)  k = 0
                    end do
                  end if
                end if
              end if           
            end do
            if (id == 3 .and. locvar(2,loopar) < 200) k = 1
            if (k /= 0) then
          !
          !  Set all derivatives to zero IF:
          !  (A) The element whose parameter is being optimized is not present, OR
          !  (B) The system is a solid, and the parameter is not a diatomic.
          !             
              do i = 1, nerr
                differ(i) = 0.d0
              end do
              used = .false.
              go to 1100
            end if
           ! if (id == 3) write(ifiles_8,"(2a,2i4,a)")" Derivative for Parameter: ", &
          !    partyp(locvar(1, loopar)), &
            !  mod(locvar(2, loopar),200), locvar(2, loopar)/200, " calculated"
            endfile (ifiles_8) 
            backspace (ifiles_8) 
            if (.not. lfast) then
        !
        !  Compute the value of all errors at the point -0.5*delta
        !
              delta = -1.d0*dstep(loopar) 
              call update (locvar(1, loopar), locvar(2, loopar), delta, 1.d0)
              call calpar
              call getusp (nat, nfirst, nlast, uspd, atheat)
              if (loopar == 8)then
                i = 0
              end if
              call parfg (errors, ttype, nerr, 0, .true.)
            end if
           !
           !  Compute the value of all errors at the point +0.5*delta
           !
            delta = 2.d0*dstep(loopar)
            call update (locvar(1, loopar), locvar(2, loopar), delta, 1.d0)
            call calpar
            call getusp (nat, nfirst, nlast, uspd, atheat) ! Update derived parameters, USPD, and ATHEAT
            call parfg (differ, ttype, nerr, 0, .true.)
           !
           !  Evaluate all derivatives
           !
            do i = 1, nerr
              differ(i) = (differ(i)-errors(i)) / delta
            end do
           !
           !  Store all the derivatives for use by the parameter optimizer
           !
     1100   do i = 1, nerr
              diffns(loopar, i+nfns) = differ(i)
            end do
           !
           ! Step -delta.  This resets the value of the parameter.
           !
            delta=-delta
              call update (locvar(1, loopar), locvar(2, loopar),delta, 1.d0)
            end do
            if (.not. lfast) then
     !
     !  Check that derivatives are being correctly evaluated.
     !
            if (used) then
              call calpar
              call getusp (nat, nfirst, nlast, uspd, atheat)
              call parfg (differ, ttype, nerr, 0, .true.)
              do i = 1, nerr
                 differ(i) = (differ(i)-errors(i)) / delta
              end do
            end if
            j = 0
            do i = 1, nerr  !  do i
              if (abs(differ(i)) > 20 .and. used) then
                j = 1
                if (big_loop == 5) then
       !
       !  Zero out the derivatives - they are not accurate, and can't be made accurate.
       !
                  do loopar = 1, numvar
                    diffns(loopar, i+nfns) = 0.d0
                  end do
                  write(ifiles_8,*)" Derivatives inaccurate - not used"
                end if
              end if
            end do ! do i
            if (j == 0) exit
          else
            exit
          end if      
        end do ! big_loop
      end if  ! numvar 
    !
    !  Reset the value of all parameters.  These were stored in VALVAR.
    !  VALVAR is updated in rapid0
    !
      do loopar = 1, numvar
        call update (locvar(1, loopar), locvar(2, loopar), valvar(loopar), 0.d0)
      end do
      
    !
    !  Store location of heat of formation error and derivatives.
    !
      loch(loop) = nfns + 1
      nfns = nfns + nerr
    !
    !  Save the density matrix for use in the next optimization cycle.
    !
      call resetp (1, loop)
    !
    !  Write out information on this system
    !
     
      if (Mod(loop,5) == 0) then
        do i = 80,2,-2
          if(koment(i:i+1) /= "  ") exit
          koment(i:i+1) = ". "
        end do
      end if
       i = Int (seconds (1))
      if(is_a_ref(loop)) then
      write (ifiles_8, "(A49,I4)") koment, i - itime
      else
      write (ifiles_8, "(A49,I4,F9.3,F9.3,F9.2,F12.3,F12.3)") koment, i - &
     & itime, refher, refder, refcer, refger, reftot
!
!   Print the contribution from pKa and the scalar of the derivative vector 
!   of the first reference function.
!
  !    if (wtpKa > 1.d-5) &
  !    & write(ifiles_8,"(a,f12.3,a,f12.3)")"Total pKa:", sumpKa," Gnorm:", &
  !    & sqrt(dot(diffns(1,nfns - nerr + 1), diffns(1,nfns - nerr + 1), numvar - 2))
        if(numvar == -10) then
!
! Write out heats of formation, only, without spaces
!
            if (refher < -999.9950d0) then
             write(55,"(f9.2)",err=95) refher
            else if (refher < -99.9950d0) then
             write(55,"(f8.2)",err=95) refher
             else if (refher < -9.9950d0) then
             write(55,"(f7.2)",err=95) refher
             else if (refher < -0.9950d0) then
             write(55,"(f6.2)",err=95) refher
             else if (refher < 0.0001d0) then
             write(55,"(f5.2)",err=95) refher
             else if (refher < 9.9950d0) then
             write(55,"(f5.2)",err=95) refher
             else if (refher < 99.9950d0) then
             write(55,"(f6.2)",err=95) refher
             else if (refher < 999.9950d0) then
             write(55,"(f7.2)",err=95) refher
            end if
  95  continue
        end if
      end if
      itime = i
      write(iw,*)trim(koment)//" end"
      endfile (iw) 
      backspace (iw) 
    end do
    write(iw,*)" End of cycle"
    do loop = 1, nmols
      if (ihrefs(1, loop) > 0) then
        call depfn (fns, diffns, loop, loch, ihrefs(1, loop), ihrefs(2, loop), &
       & nfns, numvar)
      end if
    end do
    if (numvar == -10) close (55)
  !
  !  Summarize the errors
  !
    sumerr = sumher + sumder + sumier + sumger + sumpKa
    write (ifiles_8, "(25x,A23,F14.3,F9.3,F9.2,F12.3,F12.3)") " TOTALS", sumher, &
   & sumder, sumier, sumger, sumerr
    if (sumpKa > 1.d0) &
    & write (ifiles_8, "(88x,a,f10.3,a)") " (pKa:",sumpKa,")"
    aveher = aveher / (nhmols+0.000001d0)
    erhmod = erhmod / (nhmols+0.000001d0)
    rmsher = (rmsher/(nhmols+0.000001d0))**(1.d0/power)
    aveder = aveder / (ndmols+0.000001d0)
    erdmod = erdmod / (ndmols+0.000001d0)
    rmsder = (rmsder/(ndmols+0.000001d0))**(1.d0/power)
    aveier = aveier / (nimols+1.d-4)
    rmsier = (rmsier/(nimols+1.d-4))**(1.d0/power)
    erimod = erimod / (nimols+1.d-4)
    rmsger = (rmsger/(ngmols+1.d-5))**(1.d0/power)
    ergmod = ergmod / (ngmols + 1.d-4)
    write (ifiles_8, "(21x,i4,'    AVERAGE ERROR       ',i4,f8.2,i3,f6.2,i5,f5.2,i8,/21x,i4,'&
   &    ROOT MEAN SQUARE ERROR =',f8.2,f9.2,f10.2,f11.2,/,21x,i4,'  UNSIGNED AVE. ERROR     &
   & =',f9.3,f9.3,f10.3,f11.3)") cycle_no, nhmols, aveher, ndmols, aveder, nimols, aveier, &
   & ngmols, cycle_no, rmsher, rmsder, rmsier, rmsger, cycle_no, erhmod, erdmod, erimod, ergmod
    if (fnsnew(1) > -1.d6 .and. numvar > 0) then
      j = 0
      do i = 1, nfns
        error(i) = fns(i) - fnsnew(i)
!
!   If the predicted and calculated quantities are very different, then alter 
!   the weighting for that reference datum.
!
!   The weighting factor is just a "best guess".  Currently:
!
        scale_lim = 10.d0
!
!   If the error is less than scale_lim then
!   if the weighting factor is 1.0 ignore it, otherwise increase it
!   If the error is greater than scale_lim, decrease the weighting factor
!
!
        sum = Abs(error(i))
        if (lfact_here(i)) then  
      !  if (lfact) then          
          if (sum < scale_lim) then
            factor(i) = 1.d0 ! Min(1.d0, (factor(i) + 0.01d0)*2.d0) !  Double scaling factor
          else 
            factor(i) =  1.d0/sum**2
          end if
        end if     
        if (Abs (fns(i)-fnsnew(i)) > 5.d0) then 
          if (j == 0) then
            write (ifiles_8, "(A,43X,2A)") "    Function ", "   Calculated ", " &
           & Predicted     Error      Factor"
            j = 1
          end if
          write (ifiles_8, "(1X,A52,3F12.3,f16.7)") tfns (i), fns (i), fnsnew (i), &
          & error(i), factor(i)
        end if
      end do
    else
      factor = 1.d0
    end if

!    if (full .or. nr) then
!
!  Write out the errors and their weighted derivatives
!
    fnsnew(1:numvar) = 0.d0
    sumerr = 0.d0
    if (full) write(ifiles_8,"(/,'    System                Datum      Error',&
    &'     Derivatives, weighted for errors')")
    do i = 1,nfns
    if (full) write(ifiles_8,"(/2x,A52,F12.4,10(/7x,7f12.4))") &
   & tfns(i),fns(i), &
   & (Sign(1.d0,fns(i))*Abs(fns(i))**(power - 1.d0) * diffns(j, i) * power,j = 1,numvar)
   sumerr = sumerr + Abs(fns(i)) ** power
   do j = 1,numvar
   fnsnew(j) = fnsnew(j)+Sign(1.d0,fns(i))*Abs(fns(i))**(power - 1.d0) * diffns(j, i) * power * factor(i)
   end do
    end do
  ! end if
   if (full) write(ifiles_8,"(/2x,A33,F14.4,3F14.2,10(/7x,7f14.2))") &
   & "Total Error:",sumerr,fnsnew(1:numvar) 
   write(ifiles_8,*)
   if (lderivs) then
!
!  Write out the errors and their derivatives
!
    fnsnew(1:numvar) = 0.d0
    sumerr = 0.d0
    write(ifiles_8,"(/,'    System                Datum      Error',&
    &'     Derivatives, NOT weighted for errors')")
    do i = 1,nfns
    write(ifiles_8,"(/2x,A52,F12.4,10(/7x,7f12.4))") &
   & tfns(i),fns(i), (diffns(j,i),j = 1,numvar)
   sumerr = sumerr + Abs(fns(i)) ** power
   do j = 1,numvar
   fnsnew(j) = fnsnew(j)+Sign(1.d0,fns(i))*Abs(fns(i))**(power-1.d0)*diffns(j,i) &
   & * power
   end do
    end do
       write(ifiles_8,"(/2x,A33,F14.4,3F14.4,10(/7x,7f14.4))") &
   & "Total Error:",sumerr,fnsnew(1:numvar) 
   end if
   return
   end subroutine direct
