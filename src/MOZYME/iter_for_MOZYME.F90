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

subroutine iter_for_MOZYME (ee)
    use molkst_C, only: norbs, step_num, numcal, nscf, escf, &
       & numat,  enuclr, atheat, emin, keywrd, moperr, line
!
    use chanel_C, only: iw, iend, end_fn
!
    use MOZYME_C, only : icocc, icvir, iorbs, ncocc, &
         ncvir, nncf, nnce, cocc, cvir, nvirtual, noccupied, &
         cocc_dim, icocc_dim, icvir_dim, cvir_dim, nelred, &
         norred, numred, ncf, nce, scfref, thresh, ovmax, tiny, shift, &
         energy_diff, sumt, sumb, ijc, use_three_point_extrap, &
         ws, p1, p2, p3, partf, partp, idiag, mode
!
    use funcon_C, only: fpc_9
    use common_arrays_C, only : f, p
    use iter_C, only : pold
    use cosmo_C, only: useps, lpka, solv_energy
    use linear_cosmo, only : c_proc
    implicit none
!
    double precision, intent (out) :: ee
   !
   !***********************************************************************
   !
   !                   The Occupied Set
   !
   !  COCC:     Starting localized filled M.O.s        Atomic Orbitals
   !  NCOCC:    Starting address of the filled LMOs    Atomic Orbitals
   !  ICOCC:    Atom indices for the filled LMOs       Atoms
   !  NCF:      Number of atoms in each filled LMO     Atoms
   !  NNCF:     Starting address of each M.O. in NCF   Atoms
   !
   !                   The Virtual Set
   !
   !  CVIR:     Starting localized empty M.O.s         Atomic Orbitals
   !  NCVIR:    Starting address of the empty LMOs     Atomic Orbitals
   !  ICVIR:    Atom indices for the empty LMOs        Atoms
   !  NCE:      Number of atoms in each empty LMO      Atoms
   !  NNCE:     Starting address of each M.O. in NCE   Atoms
   !
   !***********************************************************************
    character :: xchar = " "
    logical, save :: prtpls, orthog, times, scf1, bigscf, debug, prtden, &
         & prtfok, okscf, opend, panic = .true., store_useps
    integer, save :: i, idnout, nocc1, nvir1, istabl, idiagg, nhb, itrmax, &
      niter, lno, iemax, iemin, lnv, mn, indi, j, k, l, nij = 0, re_local, &
      icalcn = 0, imol = 0, nmol = 0, lstart = 0, add_niter
    double precision, save :: selcon, eold,  sum
    integer, external :: ijbo
    logical, external :: PLS_faulty
    double precision, external :: helecz, reada
    add_niter = 0
!
        80  continue
    if (nmol /= numcal) then
      !
      !  INITIALIZE
      !
      !     THRESH          ERROR IN H.O.F., RELATIVE TO THRESH=1.D-13
      !
      !     1.D-8             0.21??     SCF unstable
      !     1.D-9             0.0294     SCF stable
      !     1.D-10            0.0112     SCF stable
      !     1.D-11            0.0040     SCF stable
      !     1.D-12            0.0011     SCF stable
      !
      !   IF THRESH IS CHANGED, THEN ALSO CHANGE DEFAULT IN WRTKEY
      !
      thresh = 1.d-13
      scfref = 0.d0
      i = Index (keywrd, " RELTHR")
      if (i /= 0) then
        thresh = reada (keywrd, i) * thresh
      end if
      i = Index (keywrd, " THRESH")
      if (i /= 0) then
        thresh = reada (keywrd, i)
      end if!

      thresh = Max (thresh, 1.d-25)
      prtpls = (Index (keywrd, " PL ") /= 0)
      orthog = (Index (keywrd, " REORTH") /= 0)
      times = (Index (keywrd, " TIMES") /= 0)
      i = Index (keywrd, " DENOUT=")
      if (i /= 0) then
        idnout = Nint (reada (keywrd, i+8))
      else
        idnout = 10000000
      end if
      itrmax = 2000
      if (Index (keywrd, " ITRY") /= 0) then
        itrmax = Nint (reada (keywrd, Index (keywrd, " ITRY")))
      end if
      re_local = 100000000
      i = index(keywrd," RE-LOC")
      if (i /= 0) then
        j = index(keywrd(i + 7:), " ") + i + 7
        if (index(keywrd(i:j), "=") /= 0) then ! Allow for RE-LOC=, RE-LOCAL=, etc.
          i = index(keywrd(i:j), "=") + i
          re_local = nint(reada(keywrd,i))
        end if
      end if
      nocc1 = 0
      nvir1 = 0
      scf1 = (Index (keywrd, " 1SCF") /= 0)
      bigscf = (Index (keywrd, " BIGSCF") /= 0)
      if (Index (keywrd, " OLDEN") == 0) then
        bigscf = .true.
      end if
      debug = (Index (keywrd, " DEBUG") /= 0)
      prtden = (Index (keywrd, " DENS") /= 0 .and. debug)
      prtfok = (Index (keywrd, " FOCK") /= 0 .and. debug)
      debug = (Index (keywrd, " ITER") /= 0)
      idiagg = 0
      nocc1 = noccupied
      store_useps = useps
      if (Index (keywrd, " LEWIS") /= 0) then
        call makvec()
        return
      end if
      if (Index (keywrd, " OLDEN") /= 0) then
          call pinout(0, (index(keywrd, "SILENT") == 0))
          if (add_niter /= 0)  call l_control("OLDEN", len_trim("OLDEN"), -1)
          if (add_niter /= 0)  call l_control("SILENT", len_trim("SILENT"), -1)
          if (moperr) return
          call density_for_MOZYME (p, 0, noccupied, partp)
          partp = p
      else
        useps = .false.
        if (index(keywrd, "OLD_SCF") == 0) call makvec()
      end if
      if (moperr) return
      !
      !  A NEW MOLECULE, THEREFORE SEARCH FOR ALL WEAK INTERACTIONS.
      !
      nhb = 0
    else
      !
      !  A MODIFIED GEOMETRY, THEREFORE SEARCH ONLY FOR NEW
      !  OR VERY WEAK INTERACTIONS.
      !
      nhb = 3
    end if
    if (mod(nscf + 1, re_local) == 0) then
      write(iw,"(/10x,a,/)")"  LMOs being Re-Localized"
      call local_for_MOZYME("OCCUPIED")
      call local_for_MOZYME("VIRTUAL")
    end if
    useps = store_useps
    if (lpka) useps = .true.
    if (icalcn /= step_num) then
      istabl = 0
      eold = 0.d0
      ovmax = 0.d0
      call scfcri (selcon)
      idiagg = 0
      tiny = 0.d0
    end if
    if (index(keywrd," tighten")  /= 0) then
      call l_control("tighten", len_trim("tighten"), -1)
      selcon = max( selcon*0.1d0, 1.d-3)
    end if
   !
   ! Fill the array IDIAG pointing to the diagonal elements of P
   !
    l = 0
    do i = 1, numat
      j = ijbo (i, i)
      do k = 1, iorbs(i)
        l = l + 1
        j = j + k
        idiag(l) = j
      end do
    end do
   !
   !  Zero out POLD and P1 so that the new calculation is not affected
   !  by a previous one.
   !
    pold = 0.d0
    p1 = 0.d0
    nscf = nscf + 1
   !
   !  Force IDIAGG to be even - this is to ensure that DIAGG1
   !  re-builds the interaction list, in case the temporary space
   !  has been used between the calls to ITER.
   !
    if (Mod(idiagg, 2) == 1) idiagg = idiagg + 1
    iemin = 0
    iemax = 0
    niter = 0
    shift = 0.0d0
    use_three_point_extrap = .true.
    if (times) then
      call timer (" At start of ITER")
    end if
!***********************************************************
!
!   Everything is now set up to allow the SCF to be run
!
!***********************************************************
    do  !  Big loop to run the SCF
!----------------
      do
        call tidy (noccupied, ncf, icocc, icocc_dim, cocc, cocc_dim, nncf, ncocc, lno, mn, 1)
        if (moperr) then
!
!  During a run of "tidy", the amount of expansion space for the LMO's to use had
!  become small.  The LMOs were stored to disc.  The old arrays will now be deleted
!  and re-created 60% larger than before.  Then the LMOs are read off disc
!
          deallocate (icocc, cocc, icvir, cvir)

!
!  Re-allocate more memory
!
          icocc_dim = Nint(icocc_dim*1.6)
          cocc_dim = Nint(cocc_dim*1.6)
          icvir_dim = Nint(icvir_dim*1.6)
          cvir_dim = Nint(cvir_dim*1.6)
          allocate (icocc(icocc_dim), cocc(cocc_dim), &
                & icvir(icvir_dim), cvir(cvir_dim), stat = i)
          if (i /= 0) then
            call memory_error(" iter_for MOZYME")
            return
          end if
!
!  Read in old density
!
          call pinout (0, .false.)
          moperr = .false.
        else
          exit
        end if
      end do
      do
        call tidy (nvirtual, nce, icvir, icvir_dim, cvir, cvir_dim, nnce, ncvir, lnv, mn, 2)
        if (moperr) then
!  Delete old memory
!
          deallocate (icocc, cocc, icvir, cvir)

!
!  Re-allocate more memory
!
          icocc_dim = Nint(icocc_dim*1.6)
          cocc_dim = Nint(cocc_dim*1.6)
          icvir_dim = Nint(icvir_dim*1.6)
          cvir_dim = Nint(cvir_dim*1.6)
          allocate (icocc(icocc_dim), cocc(cocc_dim), &
                & icvir(icvir_dim), cvir(cvir_dim), stat = i)
          if (i /= 0) then
            call memory_error(" iter_for MOZYME")
            return
          end if
!
!  Read in old density
!
          call pinout (0, .false.)
          moperr = .false.
        else
          exit
        end if
      end do
!----------------
      nocc1 = nelred / 2
      nvir1 = norred - nocc1
      !
      !   REMOVE ELECTRON DENSITY DUE TO LMO'S INVOLVED IN SCF FROM
      !   THE DENSITY MATRIX
      !
      if (imol /= numcal .or. icalcn /= step_num .and. numat > numred+1) then
        !---------------------------------------------------------
        !
        !   THIS PART IS ONLY RUN WHEN ICALCN IS INCREMENTED
        !
        call density_for_MOZYME (p, 0, noccupied, partp) ! Build the whole density matrix
!
        if (times) call timer (" After DENSIT")
        if (prtden) then
          write (iw, "(' DENSITY MATRIX TO GO INTO PARTP')")
          call vecprt_for_MOZYME (p, norbs)
        end if
        call setupk (nocc1) ! Work out the atom list to be used in the SCF
        if (times) call timer (" After SETUPK")
        if (imol == numcal .and. numat > numred+1) then
          call density_for_MOZYME (partp, -1, nocc1, p) ! Remove density due to atoms to be
                                                        ! used in the SCF
          if (prtden) then
            write (iw, "(' DENSITY MATRIX IN PARTP')")
            call vecprt_for_MOZYME (partp, norbs)
          end if
        end if
        call buildf (f, partf, 0)
        if (prtfok) then
          write (iw, "(' FOCK MATRIX AT START OF ITER')")
          call vecprt_for_MOZYME (f, norbs)
        end if
        if (icalcn /= step_num) then
          if (times) call timer (" After BUILDF")
          ee = helecz()
          if (times) call timer (" After HELEC")
          escf = (ee+enuclr) * fpc_9 + atheat
          if (useps) then
                escf = escf + solv_energy * fpc_9
          end if
          if (prtpls)  write (iw, "(/,A,F16.6,A,/)") " PLS ESCF USING THE OLD LMOs:", escf, " KCAL/MOL"
          endfile (iw)
          backspace (iw)
        end if
        if (imol == numcal .and. numat > numred+1) call buildf (partf, f, -1)
        icalcn = step_num
        imol = numcal
      end if
!
!  Correct any small errors in normalization
!
      call check (nocc1, nncf, ncf, icocc, icocc_dim, iorbs, ncocc, cocc, cocc_dim)
      if (moperr) return
      call check (nvir1, nnce, nce, icvir, icvir_dim, iorbs, ncvir, cvir, cvir_dim)
      if (moperr) return
      if (Mod(niter+1, idnout) == 0) then
        write (iw, "(A)") " .den FILE TO BE WRITTEN OUT"
        endfile (iw)
        backspace (iw)
        call pinout (1, .true.)
        write (iw, "(A)") " .den FILE WRITTEN OUT"
        endfile (iw)
        backspace (iw)
      end if
      call eimp ()
      if (nmol == numcal .and. numat > numred+1) then
        indi = 1
      else
        indi = 0
      end if
      if (niter > 10 .and. add_niter == 0) then
        if (PLS_faulty()) then
!
!  When some systems are run using MOZYME, the DIAGG1 - DIAGG2 combination fails to converge,
!  and the ovmax converges to a non-zero minimum.  If the job is stopped and a <file>.den
!  is generated, then on restarting the same job, the fault is automatically corrected.
!
!  PLS_faulty detects the conditions of the failure, at run time, and silently writes out
!  the <file>.den, then after reading in the same file, it re-runs the SCF calculation.
!  This corrects the fault.
!
          numcal = numcal + 1
          add_niter = niter
          call pinout(1, .false.)
          call l_control("OLDEN", len_trim("OLDEN"), 1)
          call l_control("SILENT", len_trim("SILENT"), 1)
          nscf = nscf - 1
          goto 80
        end if
      end if
      if (bigscf .or. numcal /= 1) then
          call diagg (f, nocc1, nvir1,  idiagg,  partp, indi)
        idiagg = idiagg + 1
      else
        call density_for_MOZYME (p, mode, nocc1,  partp)
        bigscf = .true.
      end if
      niter = niter + 1
      if (Mod(niter, 3) == 0 .and. nhb < 4) then
        nhb = nhb + 1
         !
         !   Check for missed hydrogen bonds and other unusual bonds
         !
          call addhb (nocc1, nvir1, idiagg, nij, nhb)
         !
         !   If hydrogen bonds have been made, set IDIAGG even for DIAGG
         !   to make new interactions.
        if (nij /= 0 .and. Mod (idiagg, 2) == 1) then
          idiagg = idiagg + 1
        end if
      end if
      if (times) then
        if (i == 1) then
          call timer (" After DENS+1")
        else
          call timer (" After DENSIT")
        end if
      end if
      if (use_three_point_extrap) call cnvgz (p, pold, p1, p2, p3, niter, idiag)
      if (times) call timer (" After CNVG")
      if (prtden) then
        write (iw, "(' DENSITY MATRIX ON ITERATION',I4)") niter
        call vecprt_for_MOZYME (p, norbs)
      end if
      if (nmol == numcal .and. numat > numred+1) then
        call buildf (f, partf, 1)  !
      else
        call buildf (f, partf, 0)
      end if
      if (itrmax < 3) return
      if (times) call timer (" After BUILDF")
      if (prtfok) then
        write (iw, "(' FOCK    MATRIX ON ITERATION',I4)") niter
        call vecprt_for_MOZYME (f, norbs)
      end if
      ee = helecz()
      escf = (ee+enuclr) * fpc_9 + atheat
      if (times) call timer (" After HELEC")
      if (useps) then
            escf = escf + solv_energy * fpc_9
      end if
      energy_diff = escf - eold
      eold = escf
      if (Abs(ovmax) < 5.d0*selcon) then
        c_proc = 1.d0
      else
        c_proc = 5.d0*selcon/abs(ovmax)
      end if
      escf = max(-999999.d0, min(999999.d0, escf))
      if (abs(energy_diff) > 9999.D0) energy_diff = 0.D0
      if (prtpls .or. debug .and. niter > itrmax - 20) then
        write (line, "(' ITER.',i7,' PLS=', e10.3,10x,' ENERGY ',f13.5,' DELTAE',f13.7)")  &
        niter + add_niter, ovmax,   escf, energy_diff
        write(iw,"(a)")trim(line)
        call to_screen(line)
        endfile (iw)
        backspace (iw)
        if (debug) then
          write (iw, "(A,F9.6,A,F7.1,A,F9.6,A,F8.2,A,F11.3,A,I7)") "TINY:", &
               & tiny, " SUMT:", sumt, " OVMAX:", ovmax, " SUMB:", sumb, &
               & " DIFF:", energy_diff, " IJ:", ijc
        end if
      end if
      if (debug) then
        call chrge_for_MOZYME (p, ws)
        sum = 0.d0
        do i = 1, numat
          sum = sum + Abs (ws(i))
        end do
        write (iw, "(A,I4)") " Atomic Electron Population on Iteration:", &
               & niter
        write (iw, "(10F8.4)") (ws(i), i=1, numat)
        write (iw, "(A,F12.6)") " Variance:", sum
      end if
      endfile (iw)
      backspace (iw)
      call isitsc (escf, selcon, emin, iemin, iemax, okscf, niter, itrmax)
      if ( .not. bigscf .and. numcal == 1) then
        exit
      else if (okscf .and. niter > 1 .and. (emin /= 0.d0 .or. niter > 3)) then
        exit
      end if
      if (use_three_point_extrap) then
        if (mod(niter,3) == 2 .and. Abs (energy_diff) < 0.1d0) then
          use_three_point_extrap = .false.
          lstart = niter
          shift = 0.0d0
        end if
      else
        if (energy_diff > 0.0d0 .and. shift < 11.0d0 .and. &
             & niter > lstart+2) then
          shift = shift + 2.0d0
          lstart = niter
        end if
      end if
    end do
!************************************************************
!
!   The SCF equations are now solved
!
!************************************************************
    if (istabl > 100000 .and. escf-emin > 200.d0 .and. panic) then
      write (iw, "(A)") " Something disastrous has happened.  " // &
                        & "The LMOs are probably corrupt."
      write (iw, "(A)") " The job should be restarted using RESTART but" // &
             & " OLDENS should NOT be used."
      write (iw, "(/10X,A,F16.6,/10X,A,F16.6)") &
             & " Value of previous Heat of Formation:", emin, &
             & " Value of current Heat of Formation: ", escf
        !
        !          Start the SHUT command.
        !
      inquire (file=end_fn, opened=opend)
      if (opend) then
        rewind (iend)
      else
        open (unit=iend, file=end_fn, status="UNKNOWN")
      end if
      write (iend, "(A)", err=1000) xchar
      go to 1010
        !
        !  The SHUT command is faulty.  Stop everything.
        !
1000  call geout (iw)
      call mopend ("Severe fault in RESTART")
1010  panic = .false.
    end if
    if (Abs (escf-emin) < 1.d0) then
      istabl = istabl + 1
    else
      istabl = 0
    end if
    if (escf < emin .or. emin == 0.d0) then
      emin = escf
    end if
    if (.not. scf1) then
      if (numat > numred+1) then
        call density_for_MOZYME (p, 1, nocc1, partp)
      else
        call density_for_MOZYME (p, 0, nocc1, partp)
      end if
    end if
    icalcn = step_num
    imol = numcal
    if (orthog .and. Mod (nscf + 1, 10) == 1) then
      call reorth (ws)               !   Re-orthogonalize the LMO's
      call density_for_MOZYME (p, 0, noccupied, partp)
      call buildf (f, partf, 0)
      ee = helecz ()
      escf = (ee+enuclr) * fpc_9 + atheat
      if (useps) then
            escf = escf + solv_energy * fpc_9
      end if
    end if
    nmol = numcal
    return
    end subroutine iter_for_MOZYME
