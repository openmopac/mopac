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

subroutine Locate_TS
!
! Multi-step procedure for locating a transition state involving two stationary points in a
! protein mechanism
!
! Steps involved:
! A: Using lbfgs_TS, and two stationary points, locate an approximation to the transition state
! B: Using the TS option in ef and only those atoms involved in the active site, locate the transition state.
!      (This uses gradient minimization)
! C: Using normal geometry optimization optimize the positions of all other atoms
! D: Repeat B and C until the geometry is stable
!
    use common_arrays_C, only : xparam, loc, geoa, lopt, geo, grad, lopt_store
    USE molkst_C, only : nvar, keywrd, numcal, numat, density, use_ref_geo, mpack, n2elec, &
     line, step_num, escf, moperr, prt_coords
    use MOZYME_C, only : geo_1, geo_2
    use chanel_C, only : iw, iarc, input_fn, iden
    implicit none
    integer :: i, j, k, l, big_nvar, shell, loop, nloop, nset, store_mpack, store_n2elec, ipdb = 14
    double precision :: stresses(100), gradients(100), gnorm_store
    double precision, allocatable :: big_xparam(:)
    integer :: ninsite(3)
    integer, allocatable :: active_site(:)
    logical :: extra_print, exists, l_refine = .true., opend
    character :: line1*200
    double precision, external :: reada
    character :: num
    allocate(active_site(200))
    if (.not. allocated(geoa)) allocate(geoa(3,numat))
    if (allocated(geo_1)) then
      geoa(:,:numat) = geo_1(:,:numat)
      call build_active_site(active_site, ninsite)
      if (moperr) return
    end if
    gnorm_store = 0.d0
!
!   Fill array "lopt"  When two systems are used, "lopt" is the difference in connectivity (topology) of them.
!
    if (.not. allocated(geo_1)) then
!
! If geo_1 does not exist, only one system was input, i.e. GEO_REF was not used.
! so set up conditions to go straight to "Refine_TS"
!
!  Set "lopt" to
!
      shell = 0
      if (.not. allocated(lopt_store)) then
        allocate(lopt_store(3,numat))
        lopt(:3,:numat) = 1 - lopt(:3,:numat)
      end if
      lopt_store = 1
      deallocate(grad)
      allocate(grad(3*numat))
      goto 99
    end if
!
! If geo_1 does exist, so two systems were input, i.e. GEO_REF was used.
! so run LOCATE-TS in the normal way
!
    i = index(keywrd, " GNORM=")
    if (i /= 0) gnorm_store = reada(keywrd, i)
    gradients = max(4.d0,min(20.d0,sqrt(numat*0.5d0)))
    i = index(keywrd,  " LOCATE-TS") + 11
    do j = i + 1, len_trim(keywrd)
      if (keywrd(j:j) == " ") exit
    end do
    k = index(keywrd(i:j), "C:")
    if (k > 0) then
      extra_print = .true.
      i = i + k + 1
      nloop = 0
      do
        k = ichar(keywrd(i:i)) - ichar("0")
        if (k >= 0 .and. k < 10) then
          nloop = nloop + 1
          stresses(nloop) = reada(keywrd,i)
          if (stresses(nloop) > 999.d0) then
            write(line,'(a, i3, a, f8.1, a)')"Stress for point", nloop, " is", stresses(nloop)," This is too large, max. = 999.0."
            call mopend(trim(line))
            return
          end if
          do
            i = i + 1
            if (ichar(keywrd(i:i)) /= ichar(".") .and. ichar(keywrd(i:i)) < ichar("0") .or. ichar(keywrd(i:i)) > ichar("9")) exit
          end do
          i = i + 1
        else
          exit
        end if
        if (i > j) exit
      end do
    else
      if (keywrd(i:i) > "0" .and. keywrd(i:i) <= "9") then
        call mopend("The first term in LOCATE-TS must be either ""C:"" or ""SET:""")
        return
      end if
      extra_print = .false.
      stresses(1) = 3.d0
      stresses(2) = 30.d0
      stresses(3) = 30.d0
      stresses(4) = 30.d0
      stresses(5) = 30.d0
      nloop = 5
    end if
    k = 0
    if (index(keywrd, " LOCATE-TS ") == 0) then
      i = index(keywrd,  " LOCATE-TS") + 11
      do j = i + 1, len_trim(keywrd)
        if (keywrd(j:j) == " ") exit
      end do
      k = index(keywrd(i:j), "SET:")
    end if
    if (k /= 0) then
      nset = nint(reada(keywrd, i + k + 3))
      write(iw,'(/10x, a, i2, a)')"Set", nset, " selected by keyword 'SET' within keyword 'LOCATE-TS'"
      if (nset > 2) then
        write(iw,'(/10x, a)')"(This is greater than 2, so re-set to 2)"
        nset = 2
      end if
      if (active_site(1) == 0) then
        if (index(keywrd, " LET") == 0) then
          write(line,'(a,i1,a)') """SET=""",nset,""" used in LOCATE-TS, but no atoms involved in bond-making or bond-breaking!"
          call mopend( trim(line))
          write(iw,'(/10x,a,i1,a)') "(Either remove "";SET=""",nset,""" from keyword ""LOCATE-TS"""// &
            "or add keyword ""LET"" to allow job to continue)"
          return
        end if
      end if
    else
      nset = 1
      if (index(keywrd(i:j), "C:") == 0) then
        write(iw,'(/10x, a)')"By default, set 1 will be used.  To change default, see keyword 'LOCATE-TS'"
      else
        l_refine = .false.
      end if
    end if
!
!   Use optimization flags from GEO_REF
!
    nvar = 0
    do i = 1, numat
      do j = 1,3
        if (lopt_store(j,i) == 1) then
          nvar = nvar + 1
          loc(1, nvar) = i
          loc(2, nvar) = j
        end if
      end do
    end do
    loc(:,nvar + 1:) = 0
    big_nvar = 2*nvar
    allocate(big_xparam(big_nvar))
    if (nset > 0) then
      shell = ninsite(nset)
      call select_opt(shell, active_site)
      big_nvar = 2*nvar
      do i = 1, nvar
        k = loc(1,i)
        l = loc(2,i)
        big_xparam(i) = geo_1(l,k)
        big_xparam(i + nvar) = geo_2(l,k)
        xparam(i) = geo(l,k)
      end do
    else
      shell = 0
      big_xparam = 0.d0
    end if
!
! Delete GEO_REF and LOCATE-TS from the new .mop, .arc, and other files
!
    call delete_ref_key("GEO_REF=", len_trim("GEO_REF="), '" ', 2)
    call delete_ref_key("LOCATE-TS", len_trim("LOCATE-TS"), ') ', 2)
    call delete_ref_key("XYZ", len_trim("XYZ"), " ", 1)
    call delete_ref_key("OPT", len_trim("OPT"), " ", 1)
    if (nloop == 0) then
!
!  Use average of the two geometries, i.e., the average of the input data-set and reference geometries.
!
      do i = 1, nvar
        xparam(i)=0.5d0*(big_xparam(i) + big_xparam(i + nvar))
        k = loc(1,i)
        l = loc(2,i)
        geo(l,k) = xparam(i)
      end do
      if (extra_print .and. prt_coords) then
        write(iw,'(/10x,a,/)')"Average of the input and reference geometries"
        call geout (iw)
      end if
    else
!
!            Locate the transition state by climbing the barrier. The strength of the pull from the
!            other geometry is in "density"  The pull is increased in steps, with the steps in "stresses"
!
!
! Run reference geometry
!
      density = 0.d0
      call big_swap(1,2)
      numcal = numcal + 1
      step_num = step_num + 1
      call set_up_rapid("ON")
      call set_up_rapid("OFF")
      call big_swap(0,2)
      geoa(:,:numat) = geo(:,:numat)
!
! Run data-set geometry
!
      call big_swap(1,1)
      nvar = 0
      do i = 1, numat
        do j = 1,3
          if (lopt_store(j,i) == 1) then
            nvar = nvar + 1
            xparam(nvar) = geo(j,i)
          end if
        end do
      end do
      numcal = numcal + 1
      step_num = step_num + 1
      call set_up_rapid("ON")
      call set_up_rapid("OFF")
      if (gnorm_store > 0.1d0) gradients(:nloop) = gnorm_store
      do loop = 1, nloop
        density = stresses(loop)
        if(density > 99.9499d0) then
          num = "6"
        else if(density > 9.9499d0) then   !  DELETE
          num = "5"
        else
          num = "4"
        end if
  !      write(iw,'(/10x,a,f0.2,a,/)')"Constraining constant: ",density, " Kcal/mol/Angstrom^2"
        write(line1,'(a,f'//num//'.2)')"C: ",density
        write(line,'(a,S,f0.5)')"GNORM=",gradients(loop)
        write(line1(len_trim(line1) + 4:),'(a)')trim(line)
        call l_control(trim(line), len_trim(line), 1)
        extra_print = (extra_print .or. loop == nloop)
        call lbfgs_TS(big_xparam, big_nvar, escf, extra_print)
      end do
    end if
!
!           The approximate location of the transition state has now been found.
!           Now switch off the two geometry option, and prepare to refine the transition state.
!
    if (shell == 0 .and. ninsite(1) == 0) return
99  continue
    if (.not. l_refine) goto 98
    use_ref_geo = .false.
!
! Switch off the use of a GEO_REF, from here on the geometry is self-referential
!
    density = 0.d0
!
!   Select set of atoms to be used in the active site.
!
    geoa(:,:numat) = geo(:,:numat)
    line = "LET DDMIN=0.D0 GEO-OK"
    if (gnorm_store > 0.1d0) write(line(len_trim(line) + 1:),'(a,f4.1)')" GNORM=",gnorm_store
    call l_control(line, len_trim(line), 1)
    numcal = numcal + 1
    store_mpack = mpack
    store_n2elec = n2elec
    call moldat(1)
    mpack = store_mpack
    n2elec = store_n2elec
    if (shell > 0) then
      lopt(:,:numat) = 1
      do i = 1, shell
        lopt(:,active_site(i)) = 0
      end do
    end if
    call Refine_TS
    close(iarc)
    if (index(keywrd, " PDBOUT") /= 0) then
      line = input_fn(:len_trim(input_fn) - 5)//".pdb"
      call add_path(line)
      inquire(unit=ipdb, opened=opend)
      if (opend) close(ipdb)
      open(unit=ipdb, file=trim(line), status='UNKNOWN', position='asis')
      call pdbout(ipdb)
      close (ipdb)
    end if
98  line = input_fn(:len_trim(input_fn) - 4)//"den"
    inquire (file=trim(line), exist = exists)
    if (exists) then
      open(unit=iden, file=trim(line))
      close(unit=iden, status='DELETE')
    end if
    return
  end subroutine Locate_TS

  subroutine lbfgs_TS (big_xparam, big_nvar, escf_tot, extra_print)
!
!  Use the limited-memory quasi-Newton Broyden-Fletcher-Goldfarb-Shanno method for unconstrained optimization
!
!    J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980),
!    Mathematics of Computation 35, pp. 773-782.
!    D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale Optimization (1989),
!    Mathematical Programming B, 45, 3, pp. 503-528.
!
  use molkst_C, only: tleft, time0, iflepo, tdump, gnorm, keywrd, density, &
  moperr, nvar, id, line, numat, refkey, title
!
  use chanel_C, only: iw0, iw, log, ilog, input_fn, iarc
!
  use common_arrays_C, only : loc, geo, xparam, geoa
!
  use ef_C, only : nstep
!
  implicit none
  integer, intent (in) :: big_nvar
  logical :: extra_print
  double precision, intent (out) :: escf_tot
  double precision, dimension (big_nvar), intent (inout) ::  big_xparam
!
  character :: txt, csave*60, task*60, bias*7, fmt*3
  logical :: times, lsave(4), first = .true., opend
  integer :: i, k, l, itry1, jcyc, niwa, nwa, alloc_stat, m, isave(44), &
    nflush = 1
!
  double precision :: absmin, const, cycmx, stepmx, sum, tstep, tt0, rms, sum1, sum2, &
  & tlast, tolerg, tprt, tx1, tx2, oldstp(12), dsave(29), escf1, escf2, e_stress
!
  double precision, dimension(:), allocatable :: bot, gold, top, xold, wa, big_grad, store_big_grad
  integer, dimension (:), allocatable :: iwa, nbd
! For Mopac BLAS
  double precision, external :: ddot, dot, reada, seconds
!
  save :: tlast, first
!
  allocate (nbd(big_nvar), big_grad(big_nvar), store_big_grad(big_nvar), stat=alloc_stat)
  if (alloc_stat /= 0) then
    call memory_error ("mod_lbfgs")
    return
  end if
  m = 12
  nstep = 0
  niwa = 3 * big_nvar
  nwa = 2 * big_nvar * m + 4 * big_nvar + 11 * m ** 2 + 8 * m
  allocate (bot(big_nvar), gold(big_nvar), top(big_nvar), xold(big_nvar), wa(nwa), &
    iwa(niwa), stat=alloc_stat)
  wa = 0.d0
  iwa = 0
  task = " Unused"
  csave = " Unused"
  lsave = .false.
  isave = 0
  times = (Index (keywrd, " TIMES") /= 0)
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
!
!  Turn OFF all bounds checking.  This saves memory and speeds things up
!
  do i = 1, big_nvar
    nbd(i) = 0
  end do
  task = "START"
  jcyc = 0
  cycmx = 0.d0
  tlast = tleft
  itry1 = 0
  absmin = 1.d6
  do
    if (times) then
      call timer ("Before SETULB")
    end if
    call dcopy (big_nvar, big_xparam, 1, xold, 1)
!
!     THIS IS THE CALL TO THE L-BFGS-B CODE.
!
    call setulb (big_nvar, m, big_xparam, bot, top, nbd, escf_tot, big_grad, 0.d0, &
    & 0.d0, wa, iwa, task,-1, csave, lsave, isave, dsave)
    if (moperr) goto 99
    if (times) then
      call timer ("AFTER SETULB")
    end if
    dsave(2) = dsave(2) + 1.d4
    if (task(1:2) == "FG") then
      if (jcyc > 1) then
!
!  How big was the step in parameter space?
!
        sum = 0.d0
        do i = 1, big_nvar
          sum = sum + (big_xparam(i)-xold(i)) ** 2
        end do
        sum = Sqrt (sum)
        i = min(big_nvar,10)
        stepmx = Min (1.d0, dSqrt (ddot(i,oldstp, 1,oldstp, 1)/i))
        if (sum > stepmx*2.d0) then
!
! Step was too big - this was probably due to an error in SETULB.
! Reduce the step to two times STEPMX
!
          const = 2.d0 * stepmx / sum
          do i = 1, big_nvar
            big_xparam(i) = const * big_xparam(i) + (1.d0-const) * xold(i)
          end do
          sum = 2.d0 * stepmx
        end if
        oldstp(Mod(jcyc, 10)+1) = sum
      end if
!
!  Limit step to 0.2 Angstroms.  This should help damp wild swings
!
      sum = 0.2d0
      do i = 1, big_nvar
        if (abs(big_xparam(i) - xold(i)) > sum) then
          big_xparam(i) = xold(i) + max(-sum, min(sum, big_xparam(i) - xold(i)))
        end if
      end do
      if (first) then
        call  geo_diff(sum, rms, .false.)
        if (density < 0.01d0) then
          write(iw,'(/10x,a,f23.3,a)')"Current value of GEO_REF constraint:", density, " Kcal/mol/Angstrom^2"
        else
          write(iw,'(/10x,a,f22.2,a)')"Current value of GEO_REF constraint:", density, "  Kcal/mol/Angstrom^2"
        end if
        write(iw,'(10x,a,f26.2,a)') "Distance between the geometries:", sum,"  Angstroms"
      end if
      call compfg_TS (big_xparam, (mod(jcyc,222) == 0), escf1, escf2, .true., big_grad, .true.)
      if (first) then
        first = .false.
        store_big_grad(:big_nvar) = big_grad(:big_nvar)
        e_stress = 0.d0
!
!  Evaluate the cosine of the angle between the gradient-vectors for the two geometries
!  This should be need to -1.0, i.e., the angle should be about 180 degrees.
!
        do k = 1, nvar
          i = loc(1,k)
          l = loc(2,k)
          big_grad(k) = big_grad(k) + (geo(l,i) - geoa(l,i))*density*2.d0
          big_grad(k + nvar) = big_grad(k + nvar) + (geoa(l,i) - geo(l,i))*density*2.d0
          e_stress = e_stress + (geo(l,i) - geoa(l,i))**2*density
        end do
        write(iw,'(/10x,a,f19.3,a)') "Heat of formation of the first geometry:", escf1 - e_stress, " Kcal/mol"
        write(iw,'(10x,a,f18.3,a)') "Heat of formation of the second geometry:", escf2 - e_stress, " Kcal/mol"
        write(iw,'(10x,a,f11.3,a)') "Contribution to heat of formation due to stress:", 2.d0*e_stress, "  Kcal/mol"
        sum1 = sqrt(dot(big_grad, big_grad, nvar))
        sum2 = sqrt(dot(big_grad(nvar + 1), big_grad(nvar + 1), nvar))
        write(iw,'(10x,a,f22.3,a)') "Gradient arising from first geometry:", sum1, "  Kcal/mol/Angstrom"
        write(iw,'(10x,a,f21.3,a)') "Gradient arising from second geometry:", sum2,"  Kcal/mol/Angstrom"
        sum = dot(big_grad, big_grad(nvar + 1), nvar)/(sum1*sum2)
        if (sum < 0.d0) then
          write(iw,'(10x,a, f27.2, a, /)')"Angle between gradient vectors:", acos(sum)*57.2957795d0, "  degrees"
        else
          write(iw,'(10x,a, f15.2, a, /)')"WARNING! - Angle between gradient vectors: ", acos(sum)*57.2957795d0, "  degrees"
        end if
        big_grad(:big_nvar) = store_big_grad(:big_nvar)
      end if
      escf_tot = escf1 + escf2
      if (moperr) goto 99
      if (absmin-escf_tot < 1.d-7) then
        itry1 = itry1 + 1
        if (itry1 > 900 .or. (gnorm < 1.d0 .and. itry1 > 9)) then
          write (iw, "(//,' HEAT OF FORMATION IS ESSENTIALLY STATIONARY')")
          iflepo = 3
          exit
        end if
      else
        itry1 = 0
        absmin = escf_tot
      end if
      if (times) then
        call timer ("AFTER COMPFG_TS")
      end if
!
!  Write out this cycle
!
    jcyc = jcyc + 1
    tx2 = seconds (2)
    tstep = tx2 - tx1
    cycmx = Max (tstep, cycmx)
    tx1 = tx2
    tleft = tleft - tstep
    if (tlast-tleft > tdump) then
      tlast = tleft
      tt0 = seconds (1) - time0
      call lbfsav (tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, &
            & dsave, jcyc, escf_tot)
      if (moperr) goto 99
    end if
    tleft = Max (0.d0, tleft)
    call prttim (tleft, tprt, txt)
    gnorm = dSqrt (ddot(big_nvar, big_grad, 1, big_grad, 1))
!
!   Write out current status
!
    if (id == 3) then
      call write_cell(iw)
      call write_cell(iw0)
    end if
    nstep = nstep + 1
    write (line, '(" CYCLE:", i6, " TIME:", f8.3, " TIME LEFT:", &
      & f6.2, a1, "  GRAD.:", f10.3, " HEAT:", g14.7)') &
      jcyc, Min (tstep, 9999.99d0), tprt, txt, &
      & Min (gnorm, 999999.999d0), escf_tot
    write(iw,"(a)")trim(line)
    endfile (iw)
    backspace (iw)
    if (log) write (ilog, "(a)")trim(line)
    call to_screen(line)
    if (mod(jcyc,30) == 0) then
      line = trim(input_fn)
      call add_path(line)
      i = len_trim(line) - 5
      call to_screen(line(:i))
    end if
    if (nflush /= 0) then
      if (Mod(jcyc, nflush) == 0) then
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
      call dcopy (big_nvar, big_grad, 1, gold, 1)
      endfile (iw)
      backspace (iw)
!
!  EXIT CRITERIA.  (The criteria in SETULB are ignored.)
!
      if (gnorm < tolerg) then
        iflepo = 19
        exit
      end if
    else if (task(1:5) /= "NEW_X") then
      write (iw, "(2A)") " L-BFGS Message:", task
      iflepo = 9
      exit
    end if
  end do
   if (gnorm < tolerg) &
     write (iw, '(/, 5 x, "GRADIENT =", f9.5, " IS LESS THAN CUTOFF =", f9.5,//)') gnorm, tolerg
!
!  Prepare to exit: advise user of current status, in case anything goes wrong.
!
99 call  geo_diff(sum, rms, .false.)
  e_stress = 0.d0
!
!  Evaluate the cosine of the angle between the gradient-vectors for the two geometries
!  This should be need to -1.0, i.e., the angle should be about 180 degrees.
!
  do k = 1, nvar
    i = loc(1,k)
    l = loc(2,k)
    big_grad(k) = big_grad(k) + (geo(l,i) - geoa(l,i))*density*2.d0
    big_grad(k + nvar) = big_grad(k + nvar) + (geoa(l,i) - geo(l,i))*density*2.d0
    e_stress = e_stress + (geo(l,i) - geoa(l,i))**2*density
  end do
  if(density > 99.9499d0) then
    fmt = "5.1"
  else if(density > 9.9499d0) then
    fmt = "4.1"
  else if (density > 0.9499d0) then
    fmt = "3.1"
  else if (density > 0.09499d0) then
    fmt = "4.2"
  else
    fmt = "5.3"
  end if
  line = refkey(1)
  call upcase(line, len_trim(line))
!
! Remove SETUP
!
  k = index(line, " SETUP")
  if (k > 0) then
    l = k + 6
    do
      if (line(l:l) == " ") exit
      l = l + 1
    end do
    refkey(1) = refkey(1)(:k)//refkey(1)(l:)
  end if
  line = refkey(1)
  call upcase(line, len_trim(line))
  k = index(line, "GEO_DAT=") + 9
  if (k /= 9) then
    l = index(line(k:),'" ') + k
    refkey(1) = refkey(1)(:k - 10)//refkey(1)(l:)
    line = line(:k - 10)//line(l:)
  end if
  if (extra_print) then
    k = index(line, "GEO_REF=") + 9
    l = index(line(k:),'"') + k
    line = input_fn(:len_trim(input_fn) - 5)

    write(bias,'(f'//fmt//')') density
    title = "(Bias="//trim(bias)//" first) "//trim(title)
    write(line(len_trim(line) + 1:),'(a,f'//fmt//',a)')" bias=",density, " first.mop"
    call add_path(line)
    inquire(unit=iarc, opened=opend)
    if (opend) close(iarc)
    open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
    write(iw,'(/10x,a,/10x,a,f'//fmt//',a,/10x,a,/)')"First geometry (derived from data-set) after optimization subject to ", &
    &"GEO_REF constraint of ", density, " Kcal/mol/Angstrom^2 towards the reference geometry written to file:", &
    &"'"//trim(line)//"'"
    call geout (iarc)
    close(iarc)
    k = index(title, ")") + 2
    title = trim(title(k:))
    geo(:,:numat) = geoa(:,:numat)
    line = refkey(1)
    call upcase(line, len_trim(line))
    k = index(line, "GEO_REF=") + 9
    l = index(line(k:),'"') + k
    line = input_fn(:len_trim(input_fn) - 5)
    title = "(Bias="//trim(bias)//" second) "//trim(title)
    write(line(len_trim(line) + 1:),'(a,f'//fmt//',a)')" bias=",density, " second.mop"
    call add_path(line)
    inquire(unit=iarc, opened=opend)
    if (opend) close(iarc)
    open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
    write(iw,'(/10x,a,/10x,a,f'//fmt//',a,/10x,a,/)')"Second geometry (derived from reference geometry)"// &
    & "after optimization subject to ", "GEO_REF constraint of ", density, &
    " Kcal/mol/Angstrom^2 towards the data-set geometry written to file:", "'"//trim(line)//"'"
    call geout (iarc)
    close(iarc)
    k = index(title, ")") + 2
    title = trim(title(k:))
  end if
  write(iw,'(/10x,a,a)') "Job name: ","'"//input_fn(:len_trim(input_fn) - 5)//"'"
  if (density < 0.01d0) then
    write(iw,'(10x,a,f23.3,a)')"Current value of GEO_REF constraint:", density, " Kcal/mol/Angstrom^2"
  else
    write(iw,'(/10x,a,f22.2,a)')"Current value of GEO_REF constraint:", density, "  Kcal/mol/Angstrom^2"
  end if
  write(iw,'(10x,a,f24.2,a)') "Distance between the geometries:", sum,"  Angstroms"
  do i = 1, nvar
    k = loc(1,i)
    l = loc(2,i)
    geo(l,k) = big_xparam(i)
  end do
  call geo_diff(sum, rms, .true.)
  write(iw,'(/10x,a,f19.3,a)') "Heat of formation of the first geometry:", escf1 - e_stress, " Kcal/mol"
  write(iw,'(10x,a,f18.3,a)') "Heat of formation of the second geometry:", escf2 - e_stress, " Kcal/mol"
  write(iw,'(10x,a,f11.3,a)') "Contribution to heat of formation due to stress:", 2.d0*e_stress, "  Kcal/mol"
  sum1 = sqrt(dot(big_grad, big_grad, nvar))
  sum2 = sqrt(dot(big_grad(nvar + 1), big_grad(nvar + 1), nvar))
  write(iw,'(10x,a,f22.3,a)') "Gradient arising from first geometry:", sum1, "  Kcal/mol/Angstrom"
  write(iw,'(10x,a,f21.3,a)') "Gradient arising from second geometry:", sum2,"  Kcal/mol/Angstrom"
  sum = dot(big_grad, big_grad(nvar + 1), nvar)/(sum1*sum2)
  if (sum < 0.d0) then
    write(iw,'(10x,a, f27.2, a, /)')"Angle between gradient vectors:", acos(sum)*57.2957795d0, "  degrees"
  else
    write(iw,'(10x,a, f15.2, a, /)')"WARNING! - Angle between gradient vectors: ", acos(sum)*57.2957795d0, "  degrees"
  end if
  do i = 1, nvar
    xparam(i)=0.5d0*(big_xparam(i) + big_xparam(i + nvar))
    k = loc(1,i)
    l = loc(2,i)
    geo(l,k) = xparam(i)
  end do
  if (extra_print) then
    line = input_fn(:len_trim(input_fn) - 5)
    title = "(Bias="//trim(bias)//" average) "//trim(title)
    write(line(len_trim(line) + 1:),'(a,f'//fmt//',a)')" bias=",density, " average.mop"
    call add_path(line)
    inquire(unit=iarc, opened=opend)
    if (opend) close(iarc)
    open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
     write(iw,'(/10x,a,/10x,a,f'//fmt//',a,/10x,a,/)')"Average of first and second geometries after optimization subject to ", &
    &"GEO_REF constraint of ", density, " Kcal/mol/Angstrom^2 towards the data-set geometry written to file:", &
    &"'"//trim(line)//"'"
    call geout (iarc)
    close(iarc)
    k = index(title, ")") + 2
    title = trim(title(k:))
  end if
  return
  end subroutine lbfgs_TS

  subroutine compfg_TS(big_xparam, int, escf1, escf2, fulscf, big_grad, lgrad)
!
! compfg_TS evaluates the heat of formation and, if lgrad, the gradients of the two systems,
! one in the data-set and one in the reference data-set.
!
    USE molkst_C, ONLY: nvar, numat
    use common_arrays_C, only : geo, geoa, grad
    use MOZYME_C, only : geo_1, geo_2
    implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
    double precision, intent(out) :: escf1, escf2
    logical , intent(in) :: int
    logical, intent(in)  :: fulscf
    logical , intent(in) :: lgrad
    double precision, intent(in) :: big_xparam(2*nvar)
    double precision, intent (out) :: big_grad(2*nvar)
    double precision, allocatable :: xparam(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
    integer :: i, nparam = -1
    save xparam, nparam
!  Reallocate xparam as needed
   if (nparam < nvar) then
     if(nparam >= 0) deallocate(xparam)
     nparam = nvar
     allocate(xparam(nparam))
   end if
!
!  Save the current point
!
    do i = 1, nvar
      xparam(i) = big_xparam(i)
    end do
!
!  Run data-set geometry, using reference geometry as geo_ref - this is in geoa
!
    call big_swap(1,1)
    call compfg(xparam, int, escf1, fulscf, grad, lgrad)
    call big_swap(0,1)
    big_grad(:nvar) = grad(:nvar)
    geoa(:,:numat) = geo(:,:numat)
!
! Run reference geometry using data-set as reference - this is now in geoa
!
    call big_swap(1,2)
    xparam = big_xparam(nvar + 1: 2*nvar)
    call compfg(xparam, int, escf2, fulscf, grad, lgrad)
    big_grad(nvar + 1: 2*nvar) = grad(:nvar)
    call big_swap(0,2)
!
!  Re-set geo and geoa to data-set and reference geometries
!
    geo(:,:numat) = geo_1(:,:numat)
    geoa(:,:numat) = geo_2(:,:numat)
    return
  end subroutine compfg_TS

  subroutine l_control(txt, nt, mode)
!
!   l_control has two modes:
!
! (A)  If mode =  1 the word(s) in txt is(are) added to keywrd.  If the word in txt already exists
!                   in keywrd, the word in txt would first be deleted from keywrd.
! (B)  If mode = -1 the word(s) in txt is(are) removed from keywrd
!
!  txt can consists of upper and lower case letters, "-" and "_"; the last two allow hyphenated words.
!
  use molkst_C, only : keywrd
  implicit none
  integer :: nt, mode
  character :: txt*(nt), store*2000, line*2000
  integer :: i, j, trim_len, mt
  character :: local_txt*(nt), ch*1
  local_txt = txt
  do
!
!  Parse each word, one at a time.
!
    do
      if (local_txt(1:1) /= " ") exit
      local_txt = local_txt(2:)
    end do
    ch = " "
    mt = 0
    do i = 1, nt
      mt = mt + 1
      if (local_txt(mt:mt) == '"') then ! Quickly run to the closing '"'
        do mt = mt + 1, nt
          if (local_txt(mt:mt) == '"') exit
        end do
      end if
      if (local_txt(mt:mt) == ' ') exit
      if (mt == nt) exit
    end do
    line = local_txt(:mt)
!
!  "line" holds a single keyword
!
    local_txt = local_txt(mt + 1:)
    if (line(1:1) >= "0" .and. line(1:1) <= "9") then
      i = mt + 1
      else
  !
  ! find the end of the defined keyword, e.g., in pi=3.1415, that would be "i"
  !
      do i = 1, mt
        if ((line(i:i) < "A" .or. line(i:i) > "Z") .and. (line(i:i) < "a" .or. line(i:i) > "z") .and. &
        line(i:i) /= "_" .and. line(i:i) /= "-") exit
      end do
    end if
    trim_len = i - 1
!
!  If keyword already exists, then delete the keyword
!  Keyword can have an "=" sign, or be followed by a "(" or a '"'
!  in which case, find the end of the keyword.
!
    do
      if (keywrd == " ") exit
      i = index(keywrd, " "//line(:trim_len))
      if (i > 0) then
        j = i + trim_len + 1
        if (keywrd(j:j) /= " ") then
          if (keywrd(j:j) == "=") j = j + 1
          ch = " "
          if (keywrd(j:j) == "(") then
            ch = ")"
          else if (keywrd(j:j) == '"') then
            ch = '"'
          end if
          do
            j = j + 1
            if (keywrd(j:j) == ch) exit
          end do
        end if
!
!  Clumsy fix to delete existing keywords that contain quotation marks.
!
        if (j > 1) then
          if (keywrd(j - 1:j - 1) == '"') then
            store = keywrd(:i)//'"'//trim(keywrd(j:))
          else
            do
              if (keywrd(j:j) == " ") exit
              j = j + 1
            end do
            store = keywrd(:i)//trim(keywrd(j + 1:))
          end if
        end if
        keywrd = trim(store)
      else
        exit
      end if
    end do
    store = " "
    if (mode == 1) then
      i = index(keywrd, store(:mt + 50))
      keywrd = keywrd(:i)//line(:mt)//trim(keywrd(i + mt + 1:))
    end if
    if (local_txt == " ") exit
  end do
  return
  end subroutine l_control

  subroutine get_pars(stresses, gradients, relscf, cutoff, nloop)
!
!  Used for calibrating quantities
!
  use molkst_C, only : line
  implicit none
  double precision :: stresses(20), gradients(20), relscf(20), cutoff(20)
  integer :: nloop, io_stat
!
!  Local
!
    read(33,'(a)')line
    nloop = 0
    do
      read(33,'(a)', iostat = io_stat) line
      nloop = nloop + 1
      if (io_stat /= 0) exit
      read(line,*, iostat = io_stat)stresses(nloop), gradients(nloop), relscf(nloop), cutoff(nloop)
      if (io_stat /= 0) exit
    end do
    nloop = nloop - 1
    close (33)
    return
  end subroutine get_pars

  Subroutine big_swap(mode, system)
!
!  Given two entire geometries, via 'GEO_REF="data-set.moip"', store or read
!  a geometry and all arrays specific to that geometry, e.g., occupied and virtual LMO's
!
! On input:
!           mode == 0:  Store data on the current system
!           mode == 1:  Extract data on the current system
!         system == 1:  System is based on the first data-set
!         system == 2:  System is based on the second data-set

    use molkst_C, only : numat, mpack
    use cosmo_C, only: solv_energy
!
    use common_arrays_C, only : coord, geo, nbonds, ibonds, xparam, f, p, dxyz, pa, pb
    use MOZYME_C, only : icocc, icvir, ncocc, ncvir, nncf, nnce, cocc, cvir, &
        cocc_dim, icocc_dim, icvir_dim, cvir_dim, ncf, nce,  iijj, &
        partf, partp, nijbo, iij, ijall, numij, iorbs, &
        nbonds_1, ibonds_1, icocc_1, icvir_1, ncocc_1, ncvir_1, &
        nncf_1, nnce_1, nce_1, ncf_1,  norred, nelred, &
        nijbo_1, iijj_1, iij_1, ijall_1, numij_1, iorbs_1, &
        nbonds_2, ibonds_2, icocc_2, icvir_2, ncocc_2, ncvir_2, &
        nncf_2, nnce_2, nce_2, ncf_2,  parth, parth_1, parth_2, &
        nijbo_2, iijj_2, iij_2, ijall_2, numij_2, iorbs_2, &
        geo_1, geo_2, dxyz_1, dxyz_2, &
        cocc_1, cvir_1, cocc_2, cvir_2, xparam_1, xparam_2, &
        partf_1, partp_1, f_1, p_1, refnuc, &
        partf_2, partp_2, f_2, p_2
    implicit none
    integer, intent (in) :: mode, system
    integer :: icocc_dim_1, icvir_dim_1, cvir_dim_1, cocc_dim_1, mpack_1, &
    icocc_dim_2, icvir_dim_2, cvir_dim_2, cocc_dim_2, mpack_2, norred_1, norred_2, &
      nelred_1, nelred_2
    double precision :: refnuc_1, refnuc_2, solv_energy_1 = 0.d0, solv_energy_2 = 0.d0
    save
    if (mode == 0) then
!
!  Store system
!
      if (system == 1) then
        if (allocated(nbonds)) then
          if (.not. allocated(nbonds_1)) allocate (nbonds_1(numat), ibonds_1(15,numat), &
          geo_1(3,numat), dxyz_1(3*numat))
          nbonds_1 = nbonds
          ibonds_1 = ibonds
        end if
        icocc_dim_1 = icocc_dim
        icvir_dim_1 = icvir_dim
        cvir_dim_1  = cvir_dim
        cocc_dim_1  = cocc_dim
        norred_1    = norred
        nelred_1    = nelred
        refnuc_1    = refnuc
        solv_energy_1 = solv_energy
        mpack_1 = mpack
        call copy_r_2(geo,             geo_1, 3)
        call copy_r_1(dxyz,           dxyz_1)
        call copy_i_1(icocc,         icocc_1)
        call copy_i_1(icvir,         icvir_1)
        call copy_i_1(ncocc,         ncocc_1)
        call copy_i_1(ncvir,         ncvir_1)
        call copy_i_1(nncf,           nncf_1)
        call copy_i_1(nnce,           nnce_1)
        call copy_i_1(ncf,             ncf_1)
        call copy_i_1(nce,             nce_1)
        call copy_i_1(iij,             iij_1)
        call copy_i_1(iijj,           iijj_1)
        call copy_i_1(ijall,         ijall_1)
        call copy_i_1(numij,         numij_1)
        call copy_i_1(iorbs,         iorbs_1)
        call copy_i_2(nijbo,         nijbo_1, numat)
        call copy_r_1(cocc,           cocc_1)
        call copy_r_1(cvir,           cvir_1)
        call copy_r_1(xparam,       xparam_1)
        call copy_r_1(partf,         partf_1)
        call copy_r_1(partp,         partp_1)
        call copy_r_1(parth,         parth_1)
        call copy_r_1(f,                 f_1)
        call copy_r_1(p,                 p_1)
      else
        if (allocated(nbonds)) then
          if (.not. allocated(nbonds_2)) allocate (nbonds_2(numat), ibonds_2(15,numat), &
          geo_2(3,numat), dxyz_2(3*numat))
          nbonds_2 = nbonds
          ibonds_2 = ibonds
        end if
        icocc_dim_2 = icocc_dim
        icvir_dim_2 = icvir_dim
        cvir_dim_2  = cvir_dim
        cocc_dim_2  = cocc_dim
        norred_2    = norred
        nelred_2    = nelred
        refnuc_2    = refnuc
        solv_energy_2 = solv_energy
        mpack_2 = mpack
        call copy_r_2(geo,             geo_2, 3)
        call copy_r_1(dxyz,           dxyz_2)
        call copy_i_1(icocc,         icocc_2)
        call copy_i_1(icvir,         icvir_2)
        call copy_i_1(ncocc,         ncocc_2)
        call copy_i_1(ncvir,         ncvir_2)
        call copy_i_1(nncf,           nncf_2)
        call copy_i_1(nnce,           nnce_2)
        call copy_i_1(ncf,             ncf_2)
        call copy_i_1(nce,             nce_2)
        call copy_i_1(iij,             iij_2)
        call copy_i_1(iijj,           iijj_2)
        call copy_i_1(ijall,         ijall_2)
        call copy_i_1(numij,         numij_2)
        call copy_i_1(iorbs,         iorbs_2)
        call copy_i_2(nijbo,         nijbo_2, numat)
        call copy_r_1(cocc,           cocc_2)
        call copy_r_1(cvir,           cvir_2)
        call copy_r_1(xparam,       xparam_2)
        call copy_r_1(partf,         partf_2)
        call copy_r_1(partp,         partp_2)
        call copy_r_1(parth,         parth_2)
        call copy_r_1(f,                 f_2)
        call copy_r_1(p,                 p_2)
      end if
    else
!
! Extract system
!
      if (system == 1) then
        nbonds(:numat) = nbonds_1(:numat)
        ibonds(:,:numat) = ibonds_1(:,:numat)
        if (icocc_dim_1 > 0) icocc_dim = icocc_dim_1
        if (icvir_dim_1 > 0) icvir_dim = icvir_dim_1
        if (cvir_dim_1 > 0)  cvir_dim  = cvir_dim_1
        if (cocc_dim_1 > 0)  cocc_dim  = cocc_dim_1
        if (norred_1 > 0)    norred    = norred_1
        if (nelred_1 > 0)    nelred    = nelred_1
        if (refnuc_1 > 0)    refnuc    = refnuc_1
        if (abs(solv_energy_1) > 0.1d0) solv_energy = solv_energy_1
        if (mpack_1 > 0) mpack = mpack_1
        call copy_r_2(geo_1,             geo, 3)
        call copy_r_1(dxyz_1,           dxyz)
        call copy_i_1(icocc_1,         icocc)
        call copy_i_1(icvir_1,         icvir)
        call copy_i_1(ncocc_1,         ncocc)
        call copy_i_1(ncvir_1,         ncvir)
        call copy_i_1(nncf_1,           nncf)
        call copy_i_1(nnce_1,           nnce)
        call copy_i_1(ncf_1,             ncf)
        call copy_i_1(nce_1,             nce)
        call copy_i_1(iij_1,             iij)
        call copy_i_1(iijj_1,           iijj)
        call copy_i_1(ijall_1,         ijall)
        call copy_i_1(numij_1,         numij)
        call copy_i_1(iorbs_1,         iorbs)
        call copy_i_2(nijbo_1,         nijbo, numat)
        call copy_r_1(cocc_1,           cocc)
        call copy_r_1(cvir_1,           cvir)
        call copy_r_1(xparam_1,       xparam)
        call copy_r_1(partf_1,         partf)
        call copy_r_1(partp_1,         partp)
        call copy_r_1(parth_1,         parth)
        call copy_r_1(f_1,                 f)
        call copy_r_1(p_1,                 p)
        if (allocated(pa)) pa = p*0.5d0
        if (allocated(pb)) pb = pa
      else
        nbonds(:numat) = nbonds_2(:numat)
        ibonds(:,:numat) = ibonds_2(:,:numat)
        if (icocc_dim_2 > 0) icocc_dim = icocc_dim_2
        if (icvir_dim_2 > 0) icvir_dim = icvir_dim_2
        if (cvir_dim_2 > 0)  cvir_dim  = cvir_dim_2
        if (cocc_dim_2 > 0)  cocc_dim  = cocc_dim_2
        if (norred_2 > 0)    norred    = norred_2
        if (nelred_2 > 0)    nelred    = nelred_2
        if (refnuc_2 > 0)    refnuc    = refnuc_2
        if (abs(solv_energy_2) > 0.1d0) solv_energy = solv_energy_2
        if (mpack_2 > 0) mpack = mpack_2
        call copy_r_2(geo_2,             geo, 3)
        call copy_r_1(dxyz_2,           dxyz)
        call copy_i_1(icocc_2,         icocc)
        call copy_i_1(icvir_2,         icvir)
        call copy_i_1(ncocc_2,         ncocc)
        call copy_i_1(ncvir_2,         ncvir)
        call copy_i_1(nncf_2,           nncf)
        call copy_i_1(nnce_2,           nnce)
        call copy_i_1(ncf_2,             ncf)
        call copy_i_1(nce_2,             nce)
        call copy_i_1(iij_2,             iij)
        call copy_i_1(iijj_2,           iijj)
        call copy_i_1(ijall_2,         ijall)
        call copy_i_1(numij_2,         numij)
        call copy_i_1(iorbs_2,         iorbs)
        call copy_i_2(nijbo_2,         nijbo, numat)
        call copy_r_1(cocc_2,           cocc)
        call copy_r_1(cvir_2,           cvir)
        call copy_r_1(xparam_2,       xparam)
        call copy_r_1(partf_2,         partf)
        call copy_r_1(partp_2,         partp)
        call copy_r_1(parth_2,         parth)
        call copy_r_1(f_2,                 f)
        call copy_r_1(p_2,                 p)
        if (allocated(pa)) pa = p*0.5d0
        if (allocated(pb)) pb = pa
      end if
      coord(:,:numat) = geo(:,:numat)
    end if
  contains
!
    subroutine copy_i_1(from, to)
      integer, allocatable :: from(:), to(:)
      if (allocated(from)) then
        if (allocated(to)) deallocate(to)
        allocate(to(size(from)))
        to = from
      end if
    end subroutine copy_i_1

    subroutine copy_i_2(from, to, dim_1)
      integer, allocatable :: from(:,:), to(:,:)
      integer :: dim_1
      if (allocated(from)) then
        if (allocated(to)) deallocate(to)
        allocate(to(dim_1,size(from)/dim_1))
        to = from
      end if
    end subroutine copy_i_2
!
    subroutine copy_r_1(from, to)
      double precision, allocatable :: from(:), to(:)
      if (allocated(from)) then
        if (allocated(to)) deallocate(to)
        allocate(to(size(from)))
        to = from
      end if
    end subroutine copy_r_1

    subroutine copy_r_2(from, to, dim_1)
      double precision, allocatable :: from(:,:), to(:,:)
      integer :: dim_1
      if (allocated(from)) then
        if (allocated(to)) deallocate(to)
        allocate(to(dim_1,size(from)/dim_1))
        to = from
      end if
    end subroutine copy_r_2
!
  end subroutine big_swap
  subroutine build_active_site(active_site, ninsite)
!
!  The set of atoms in the active site is constructed.
!  Given the topology of the reactants and products (in nbonds and ibonds),
!  the active site atoms are identified by the differences in topology.
!  This is typically 2 - 6 atoms.  By using the nearest neighbors and second nearest neighbors,
!  a good approximation to the active site can be constructed.
!
    use molkst_C, only : numat, line, pdb_label, maxtxt
!
    use common_arrays_C, only : nbonds, ibonds, txtatm, nat
!
    USE chanel_C, only : iw
!
    use elemts_C, only: elemnt
    implicit none
    integer, intent (out) :: active_site(200), ninsite(3)
    integer i, ii, j, jj, k, l, nsite
    integer, allocatable :: nbonds_b(:), ibonds_b(:,:)
    character :: ch*2
    active_site = 0
    ninsite = 0
    allocate (nbonds_b(numat), ibonds_b(15,numat))
    call big_swap(1,1)!  Extract system 1 - the input data set
    nbonds_b = nbonds
    ibonds_b = ibonds
    call big_swap(1,2)!  Extract system 2 - the reference geometry
    call set_up_dentate
    call big_swap(0,2)!  Store system 2 - this adds nbonds and ibonds to the store
!
!  Compare topologies
!
    nsite = 0
    do i = 1, numat
      if (nbonds(i) /= nbonds_b(i)) then
        do k = 1, nsite
          if (active_site(k) == i) exit
        end do
        if (k <= nsite) cycle
        nsite = nsite + 1
        active_site(nsite) = i
      else
        do j = 1, nbonds(i)
!
!  At this point, atoms are in exactly the same sequence, so no need to check every permutation
!
          if (ibonds(j,i) /= ibonds_b(j,i)) then
            do k = 1, nsite
              if (active_site(k) == i) exit
            end do
            if (k <= nsite) cycle
            nsite = nsite + 1
            active_site(nsite) = i
          end if
        end do
      end if
    end do
    if (nsite == 0) return
    k = 0
    do i = 1, nsite
      k = max(k,nbonds_b(active_site(i)))
    end do
    write(ch,'(i2)')max(1, k*6 - 6)

    if (pdb_label) then
      write(iw,'(/,a,'//ch//'x,a,/16x,a)') &
      " Set 1: Atoms involved in covalent bond-breaking and bond-making", "Connectivity of atoms"
      write(iw,'(77x,a,'//ch//'x,a)')"Data-set", "Reference"
      write(iw,'(37x,a)')"(PDB labels from GEO-REF)"
    else
      write(iw,'(/6x,a,'//ch//'x,a)')" Set 1: Atoms involved in covalent", "Connectivity of atoms"
      write(iw,'(7x, a,21x,a,/)')"bond-breaking and bond-making            Data-set", "Reference"
    end if
    write(ch,'(i2)')k*6
    do i = 1, nsite
      write(line,'(i12,10i6)')(ibonds_b(j, active_site(i)), j = 1, nbonds_b(active_site(i)))
      write(line(k*6 + 10:),'(i6,10i6)')(ibonds(j, active_site(i)), j = 1, nbonds(active_site(i)))
      if (pdb_label) then
        write(iw,'(i4, 3x, a, i6, a, a)')i, " Atom:  "//elemnt(nat(active_site(i))), active_site(i), &
    "  PDB Label: ("//txtatm(active_site(i))(:maxtxt)//")", trim(line)
      else
         write(iw,'(i12, 3x, a, i6, 6x, a)')i, " Atom:  "//elemnt(nat(active_site(i))), active_site(i), trim(line)
      end if
    end do
    ninsite(1) = nsite
!
! Select neighbors
!
    ii = nsite
    do jj = 1, ii
      i = active_site(jj)
      do j = 1, nbonds(i)
        k = ibonds(j,i)
        do l = 1, nsite
          if (active_site(l) == k) exit
        end do
        if (l <= nsite) cycle
        nsite = nsite + 1
        active_site(nsite) = k
      end do
    end do
    do jj = 1, ii
      i = active_site(jj)
      do j = 1, nbonds_b(i)
        k = ibonds_b(j,i)
        do l = 1, nsite
          if (active_site(l) == k) exit
        end do
        if (l <= nsite) cycle
        nsite = nsite + 1
        active_site(nsite) = k
      end do
    end do
    write(iw,'(/,a,/)')" Set 2 : Atoms involved in covalent bond-breaking "// &
      "and bond-making, plus nearest neighbors"
    if (pdb_label) then
      write(iw,'(i4, 3x, a, i6, a)')(i, " Atom:  "//elemnt(nat(active_site(i))), active_site(i), &
  "  PDB Label: ("//txtatm(active_site(i))(:maxtxt)//")", i = 1, nsite)
    else
      write(iw,'(i12, 3x, a, i6)')(i, " Atom:  "//elemnt(nat(active_site(i))), active_site(i), i = 1, nsite)
    end if
    ninsite(2) = nsite
    return
  end subroutine build_active_site
  subroutine select_opt(shell, active_site)
  use molkst_C, only: nvar
  use common_arrays_C, only: loc, xparam, geo, geoa
  implicit none
  integer, intent (in) :: shell, active_site(200)
  integer :: i, j, k, l, use(shell), tmp(shell)
  double precision :: x_1(3*shell), x_2(3*shell)
  tmp = active_site(:shell)
  do i = 1, shell
    k = 100000
    l = 0
    do j = 1, shell
      if (tmp(j) < k) then
        k = tmp(j)
        l = j
      end if

    end do
    tmp(l) = 200000
    use(i) = k
  end do
  if (.false.) then
  nvar = 0
  do i = 1, shell
    do j = 1, 3
      nvar = nvar + 1
      loc(1,nvar) = use(i)
      loc(2,nvar) = j
      x_1(nvar) = geoa(j,use(i))
      x_2(nvar) = geo(j,use(i))
    end do
  end do
  call big_swap(1,1)
  xparam(:nvar) = x_1(:nvar)
  call big_swap(0,1)
  call big_swap(1,2)
  xparam(:nvar) = x_2(:nvar)
  call big_swap(0,2)
  end if
  return
  end subroutine select_opt

  subroutine Refine_TS
!
!       Refine the transition state by minimizing the heat of formation of everything
!       except the reaction side, followed by gradient minimization of the active site.
!       This is an iterative process, so up to five cycles of
!       (HoF minimization followed by gradient minimization) are used
!
  use common_arrays_C, only : xparam, loc, lopt, geo, f, h, lopt_store
  USE molkst_C, only : nvar, numat, moperr, keywrd
  use ef_C, only: nstep
  use MOZYME_C, only : partf, parth, mode
  implicit none
  integer :: loop, i, j
  logical :: converged, extra_print, l_ts, l_nllsq, l_sigma
    i = index(keywrd,  " LOCATE-TS")
    do j = i + 12, min(240, i + 100)
      if (keywrd(j:j) == " ") exit
    end do
    extra_print = (index(keywrd(i:j), "C:") > 0)
    l_ts = (index(keywrd, " TS") /= 0)
    l_nllsq = (index(keywrd, " NLLSQ") /= 0)
    l_sigma = (index(keywrd, " SIGMA") /= 0)
!
!  Set default gradient minimization method.
!
    if (.not. (l_ts .or. l_nllsq .or. l_sigma)) l_ts = .true.
    if (.not. l_sigma) then
      call l_control("NLLSQ", len("NLLSQ"), -1)
      call l_control("TS", len("TS"), -1)
    end if
    do loop = 1,5
!
! Set up geometry optimization calculation for everything except the active site.
!
      if (loop > 1) then
        do i = 1, numat
          do j = 1, 3
            lopt(j,i) = 1 - lopt(j,i)
          end do
        end do
      end if
!
!  lopt = geometric variables to be optimized.  This consists of all atoms not involved in bond-making or bond-breaking.
!  lopt_store = all geometric variables to be optimized.
!
! At this point, the set to be optimized is (lopt .and. lopt_store)
!
      loc = 0
      nvar = 0
      do i = 1, numat
        do j = 1, 3
          if (lopt(j,i) == 1 .and. lopt_store(j,i) == 1) then
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = j
            xparam(nvar) = geo(j,i)
          end if
        end do
      end do
      h = 0.d0
      f = 0.d0
      parth = 0.d0
      partf = 0.d0
      mode = 0
      call minimize_energy(loop)
      converged = (nstep < 3)
!
! Set up gradient-minimum calculation for the active site, i.e., where the imaginary mode is.
!
      if (moperr) return
      do i = 1, numat
        do j = 1, 3
          lopt(j,i) = 1 - lopt(j,i)
        end do
      end do
!
!  lopt = geometric variables to be optimized.  This consists of all atoms involved in bond-making or bond-breaking.
!  lopt_store = all geometric variables to be optimized.
!
! At this point, the set to be optimized is (lopt .and. lopt_store)
!
      loc = 0
      nvar = 0
      do i = 1, numat
        do j = 1, 3
          if (lopt(j,i) == 1 .and. lopt_store(j,i) == 1) then
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = j
            xparam(nvar) = geo(j,i)
          end if
        end do
      end do
      call minimize_gradient(loop, converged, extra_print, l_ts, l_nllsq, l_sigma)
      converged = (converged .and. nstep < 3)
      if (converged) exit
    end do
  end subroutine Refine_TS
  subroutine minimize_energy(loop)
    USE molkst_C, only : numat, refkey, line, numcal, nvar, time0, escf, gnorm, keywrd
    use common_arrays_C, only : xparam, grad, geo, loc
    use chanel_C, only : iw
    implicit none
    integer, intent (in) :: loop
    double precision :: gnorm_lim
!
!
!
    integer :: i, k, l
    double precision, external :: seconds, reada


      call l_control("TS", len("TS"), -1)
      write(iw,'(a, i4 ,a,/)')"  Loop:", loop, &
      "  Energy minimization, excluding active site, using L-BFGS"
      gnorm_lim = nint(numat**0.25d0*2.d0 + 1.d0)
      line = refkey(1)
      call upcase(line, len_trim(line))
      i = index(line, " GNORM=")
      if (i /= 0) gnorm_lim = reada(line, i)
      write(line,'(a,f0.1)')"DDMIN=0 GNORM=",gnorm_lim
      call l_control(trim(line), len_trim(line), 1)
      numcal = numcal + 1
      time0 = seconds(1)
      if (nvar > 0) then
        call lbfgs(xparam,escf)
        if (gnorm < gnorm_lim) then
          i = index(keywrd, " GNORM")
          write (iw, '(/, 5 x, "GRADIENT =", f9.5, " IS LESS THAN CUTOFF =", f9.5,//)') gnorm, &
          gnorm_lim
        end if
        do i = 1, nvar
          k = loc(1,i)
          l = loc(2,i)
          geo(l,k) = xparam(i)
        end do
      else
        call compfg (xparam, .TRUE., escf, .TRUE., grad, .FALSE.)
        gnorm = 0.d0
      end if
      return
  end subroutine minimize_energy
  subroutine minimize_gradient(loop, converged, extra_print, l_ts, l_nllsq, l_sigma)
    USE molkst_C, only : numat, line, nvar, time0, escf, title, moperr, iflepo
    use common_arrays_C, only : xparam, geo, loc, geoa
    use chanel_C, only : iw, input_fn, iarc
    implicit none
    integer, intent (in) :: loop
    logical, intent (in) :: extra_print, l_ts, l_nllsq, l_sigma
    logical, intent (out) :: converged
!
!  Local
!
    integer :: i
    double precision, external :: seconds
    logical :: opend
    character :: num*1
    if (l_ts) then
      line = "TS"
    else if (l_nllsq) then
      line = "NLLSQ"
    else
      line = "SIGMA"
    end if
    write(iw,'(a, i4 ,a)')"  Loop:", loop, &
    "  Gradient minimization of atoms in the active site using "//trim(line)
    i = int(log10(loop + 0.05))
    num = char(ichar("1") + i)
    do
      if (title(1:5) /= "(Loop") exit
      line = trim(title(8:))
      title = trim(line)
    end do
    write(line,'("(Loop", i'//num//',")", a)')loop, trim(title)
    title = trim(line)
    line = input_fn(:len_trim(input_fn) - 5)
    write(line(len_trim(line) + 1:),'(a, i'//num//', a)')" Loop",loop, ".mop"
    if (extra_print) then
      call add_path(line)
      inquire(unit=iarc, opened=opend)
      if (opend) close(iarc)
      open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
      write(iw,'(/10x,a,i2,a,/10x,a,/)')"Transition state on cycle",loop, " written to file:", &
      &"'"//trim(line)//"'"
      call geout (iarc)
      close(iarc)
      title = trim(title(8 + i:))
    end if
    call l_control("GNORM=3", len("GNORM=3"), 1)
    geoa(:,:numat) = geo(:,:numat)
    time0 = seconds(1)
    if (loop == 2) call l_control("OLD_HESS", len("OLD_HESS"), 1)
    if (l_ts) then
      call l_control("TS", len("TS"), 1)
      call l_control("OLD_SCF", len("OLD_SCF"), 1)
      call l_control("tighten", len_trim("tighten"), 1)
      call ef(xparam,escf)
    else if (l_sigma) then
      call powsq()
    else if (l_nllsq) then
      call nllsq()
    end if
    iflepo = 19
    if (moperr) then
      write(iw,'(//2x,a,//)')" Gradient minimization failed.  The best geometry at this point will be output to ARC file"
      converged = .true.
      iflepo = 8
      line = input_fn(:len_trim(input_fn) - 5)
      write(line(len_trim(line) + 1:),'(a)')".arc"
      call add_path(line)
      inquire(unit=iarc, opened=opend)
      if (opend) close(iarc)
      open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
      do i = 1, nvar
        geo(loc(2,i),loc(1,i)) = xparam(i)
      end do
      call geout (iarc)
      close(iarc)
      return
    end if
    return
  end subroutine minimize_gradient
