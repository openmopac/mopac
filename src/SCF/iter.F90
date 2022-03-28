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

      subroutine iter(ee, fulscf, rand)
      use common_arrays_C, only : eigs, p, pa, pb, cb, h, &
       &  c, nat, nfirst, nlast, eigb, pdiag, f, w, wk, fb
      use iter_C, only : pold, pold2, pbold, pbold2, &
      & pold3, pbold3, vec_ai, vec_bi, fock_ai, fock_bi, p_ai, &
      p_bi, h_ai, h_bi, vecl_ai, vecl_bi
      USE parameters_C, only :
      USE funcon_C, only : fpc_9
      USE maps_C, ONLY: latom
      USE chanel_C, only : iw, ifiles_1
      USE molkst_C, ONLY: numat, norbs, nalpha, nbeta, uhf, &
    &    nclose, nopen, fract, numcal, mpack, iflepo, iscf, &
    &    enuclr, keywrd, gnorm, moperr, last, nscf, emin, &
         limscf, atheat, is_PARAM, id, line, lxfac, nalpha_open, &
         nbeta_open, npulay, method_indo
      USE reimers_C, only: dd, ff, tot, cc0, aa, dtmp, nb2
      use cosmo_C, only : useps
#ifdef GPU
      Use mod_vars_cuda, only: lgpu, real_cuda, prec
      use density_cuda_i
#endif
      implicit none
      double precision , intent(out) :: ee
      logical , intent(in) :: fulscf
      logical , intent(in) :: rand
      double precision :: selcon
      integer :: l, icalcn, itrmax, na2el, na1el, nb1el, ifill, &
        irrr, jalp, ialp, jbet, ibet, ihomo, ihomob, i, j, iemin, &
        iemax, iredy, niter, modea, modeb, nl1, nl2, nu1, nu2
      double precision, dimension(numat) :: q
      double precision, dimension(10) :: escf0
      double precision :: plb, scfcrt, pl, bshift, pltest, trans, w1, w2, random, &
        shift, shiftb = 0.d0, shfmax, ten, tenold, plchek, scorr, shfto, &
        shftbo, titer0, eold, diff, enrgy, titer, escf, &
        sellim, sum, summ, eold_alpha, eold_beta, theta(norbs), ofract, sum1, sum2
      logical :: debug, prtfok, prteig, prtden, prt1el, minprt, newdg, prtpl, &
        prtvec, camkin, ci, okpuly, oknewd, times, force, allcon, &
        halfe, gs, capps, incitr, timitr, frst, bfrst, ready, glow,  &
        makea, makeb, getout, l_param, opendd
      integer :: iopc_calcp
      character, dimension(3) :: abprt*5
      double precision, external :: capcor, helect, meci, reada, seconds

      save icalcn, debug, prtfok, prteig, prtden, prt1el, abprt, plb, &
        minprt, newdg, scfcrt, prtpl, prtvec, pl, bshift, pltest, itrmax, na2el, &
        na1el, nb1el, ifill, camkin, ci, okpuly, oknewd, &
        times, force, allcon, trans, halfe, w1, w2, random, gs, shift, shiftb, &
        shfmax, capps, ten, tenold, incitr, irrr, plchek, timitr, &
        scorr, shfto, shftbo, titer0, jalp, ialp, jbet, ibet, enrgy, &
        iopc_calcp
      data icalcn/ 0/
      data debug/ .FALSE./
      data prtfok/ .FALSE./
      data titer0/ 0.D0/
      data prteig/ .FALSE./
      data prtden/ .FALSE./
      data prt1el/ .FALSE./
      data ten/ 10.D0/
      data tenold/ 10.D0/
      data plb/ 0.D0/
      data scorr/ 0.D0/
      data abprt/ '     ', 'ALPHA', ' BETA'/
!
!  INITIALIZE
!
      ifill = 0
      ihomo = max(1,nclose + nalpha)
      ihomob = max(1,nclose + nbeta)
      eold = 1.D2
      ready = .FALSE.
      diff = 0.D0
      escf0 = 0.D0
      sellim = 0.d0
      opendd = .false.
      glow = .FALSE.
      if (icalcn /= numcal) then
        call delete_iter_arrays
        l_param = .true.
        enrgy = fpc_9
        irrr = 5
        shift = 0.D0
        icalcn = numcal
        shfmax = 20.D0
!
!    DEBUG KEY-WORDS WORKED OUT
!
        debug = index(keywrd,' DEBUG') /= 0
        minprt = index(keywrd,' SADDLE') + latom==0 .or. debug
        prteig = index(keywrd,' EIGS') /= 0
        prtpl = index(keywrd,' PL ') + index(keywrd,' PLS') /= 0
        prt1el = index(keywrd,' 1ELE') /=0 .and. debug
        prtden = index(keywrd,' DENS') /=0 .and. debug
        prtfok = index(keywrd,' FOCK') /=0 .and. debug
        prtvec = index(keywrd,' VEC') + index(keywrd,' ALLVEC') /=0 .and. debug
        debug = index(keywrd,' ITER') /= 0
!
! INITIALIZE SOME LOGICALS AND CONSTANTS
!
        newdg = .FALSE.
        camkin = .false.
        plchek = 0.005D0
        pl = 1.D0
        plb = 0.d0
        bshift = -80.D0
        shift = 1.D0
        shfto = 0.D0
        shftbo = 0.D0
        itrmax = 2000
        na2el = nclose
        na1el = nalpha + nopen
        nb1el = nbeta + nopen
!
!  USE KEY-WORDS TO ASSIGN VARIOUS CONSTANTS
!
        if (index(keywrd,' FILL') /= 0) ifill = -nint(reada(keywrd,index(keywrd,' FILL')))
        if (index(keywrd,' SHIFT') /= 0) bshift = -reada(keywrd,index(keywrd,' SHIFT'))
        if (Abs(bshift) > 1.d-20) ten = bshift
        if (index(keywrd,' ITRY') /= 0) itrmax = nint(reada(keywrd,index(keywrd,' ITRY')))
        ci = index(keywrd,' MICROS') + index(keywrd,' C.I.') /= 0 .and. .not. method_indo
        okpuly = index(keywrd,' PULAY') /= 0
        oknewd = abs(bshift) < 0.001D0
        if (camkin .and. abs(bshift)>1.D-5) bshift = 4.44D0
        times = index(keywrd,' TIMES') /= 0
        timitr = times
        force = index(keywrd,' FORCE') /= 0
        gs = index(keywrd,' TS') + index(keywrd,' NLLSQ') + index(keywrd,' SIGMA') == 0 .and. .not.force
        allcon = okpuly .or. camkin
!
!   DO WE NEED A CAPPED ATOM CORRECTION?
!
        j = 0
        j = j + count(nat(:numat) == 102)
        capps = j > 0
        iscf = 1
        trans = 0.200D0
        if (index(keywrd,' OLDENS') /= 0) then
           call den_in_out(0)
           if (moperr) return
           if (uhf) then
            pold(1:mpack) = pa(1:mpack)
            pbold(1:mpack) = pb(1:mpack)
          else
            pold(1:mpack) = pa(1:mpack)*2.d0
          end if
        else
!
! If a reaction path of some kind and it's a UHF, don't reset the density matrix after the first SCF
!
          if (.not. is_PARAM .and. (latom == 0 .or. nscf == 0) .or. .not. UHF) then
            p(:mpack) = 0.D0
            pa(:mpack) = 0.D0
            pb(:mpack) = 0.D0
            w1 = na1el/(na1el + 1.D-6 + nb1el)
            w2 = 1.D0 - w1
            if (w1 < 1.D-6) w1 = 0.5D0
            if (w2 < 1.D-6) w2 = 0.5D0
!
!  SLIGHTLY PERTURB THE DENSITY MATRIX IN CASE THE SYSTEM IS
!  TRAPPED IN A S**2 = 0 STATE.
!
            random = 1.0D0
            glow = method_indo .or. glow .or. (gnorm<2.D0 .and. gnorm > 1.D-9)
            if (.not. glow .and. uhf .and. na1el == nb1el) random = 1.1D0
            do i = 1, norbs
              j = (i*(i + 1))/2
              p(j) = pdiag(i)
              pa(j) = p(j)*w1*random
              random = 1.D0/random
              pb(j) = p(j)*w2*random
            end do
            if (uhf) then
              do i = 1, norbs
                random = 1.D0/random
                pb((i*(i+1))/2) = p((i*(i+1))/2)*w2*random
              end do
            end if
          end if
          pold(1:mpack) = pa(1:mpack)
          if (uhf) then
            pbold(1:mpack) = pb(1:mpack)
          end if
          do i = 1, norbs
            pold2(i) = pold((i*(i+1))/2)
          end do
        end if
        halfe = (nopen /= nclose .and. Abs(fract - 2.D0) > 1.d-20 .and. Abs(fract) > 1.d-20)
        if (halfe) then
          iopc_calcp = 3            ! DGEMM on CPU
#ifdef GPU
          if (lgpu) iopc_calcp = 2  ! DGEMM on GPU
#endif
        else
          iopc_calcp = 5            ! DSYRK on CPU
#ifdef GPU
          if (lgpu) iopc_calcp = 4  ! DSYRK on GPU
#endif
        end if
!
        if (gs) gs = .not. halfe .and. .not.ci
!
!   DETERMINE THE SELF-CONSISTENCY CRITERION
!
!
! SCFCRT IS MACHINE-PRECISION DEPENDENT
!
        scfcrt = 1.D-4
!
!  INCREASE PRECISION FOR EVERYTHING EXCEPT NORMAL GROUND-STATE
!  CALCULATIONS
!
        if (index(keywrd,' NLLSQ') + index(keywrd,' SIGMA') + index(keywrd,&
          ' TS')/=0 .or. force) then
          scfcrt = scfcrt*0.001D0
        else if (index(keywrd,' PRECISE')/=0 .or. nopen/=nclose) then
          scfcrt = scfcrt*0.01D0
        end if
        if (index(keywrd,' POLAR') /= 0) scfcrt = min(1.d-6, scfcrt)
        scfcrt = max(scfcrt,1.D-12)
!
!  THE USER CAN STATE THE SCF CRITERION, IF DESIRED.
!
        i = index(keywrd,' SCFCRT')
        j = index(keywrd,' RELSCF')
        if (i /= 0) then
          scfcrt = reada(keywrd,i)
        else if (j /= 0) then
          scfcrt = reada(keywrd,j)*scfcrt
        end if
!
!  For solids, reduce the SCF criterion to match a system with ~20 atoms
!
        if (id == 3) scfcrt = scfcrt*numat/20
        if (debug .or. i + j /= 0) write (iw, '(''  SCF CRITERION ='',G14.4)') scfcrt
        if (scfcrt < 1.D-12) write (iw, &
      '(2/2X,'' THERE IS A RISK OF INFINITE LOOPING WITH THE SCFCRT LESS THAN 1.D-12'')')
!
!   END OF INITIALIZATION SECTION.
!
      else if (nscf>0 .and. .not.uhf) then
!
!   RESET THE DENSITY MATRIX IF MECI HAS FORMED AN EXCITED STATE.  THIS
!   PREVENTS THE SCF GETTING TRAPPED ON AN EXCITED STATE, PARTICULARLY
!   IF THE PULAY CONVERGER IS USED.
!
! GBR_new_addition

        call dcopy(mpack, pa, 1, pb, 1)
        forall (i=1:mpack) p(i) = 2.d0*pa(i)

!          pb(:mpack) = pa(:mpack)
!          p(:mpack) = 2.D0*pa(:mpack)
!**
      end if
!
!   INITIALIZATION OPERATIONS DONE EVERY TIME ITER IS CALLED
!
      makea = .TRUE.
      makeb = .TRUE.
      iemin = 0
      iemax = 0
      if (irrr /= 5) then
        if (uhf) then
          call dcopy(mpack, pa, 1, pold, 1)
          call dcopy(mpack, pb, 1, pbold, 1)
          forall (i=1:norbs)
            pold2(i) = pa((i*(i+1))/2)
            pbold2(i) = pb((i*(i+1))/2)
          endforall
        else
          call dcopy(mpack, p, 1, pold, 1)
          forall (i=1:norbs)
            pold2(i) = p((i*(i+1))/2)
          endforall
        end if
      end if
      camkin = index(keywrd,' KING') + index(keywrd,' CAMP') /= 0
!
!  TURN OFF SHIFT IF NOT A FULL SCF.
!
      if (.not.fulscf) shift = 0.D0
      if (newdg) newdg = abs(bshift) < 0.001D0
      if (last == 1) newdg = .FALSE.
!
!   SELF-CONSISTENCY CRITERIA: SELCON IS IN KCAL/MOL, PLTEST IS
!   A LESS IMPORTANT TEST TO MAKE SURE THAT THE SELCON TEST IS NOT
!   PASSED 'BY ACCIDENT'
!                              IF GNORM IS LARGE, MAKE SELCON BIGGER
!
      selcon = scfcrt
!
!  LET SELCON BE DETERMINED BY SCFCRT AND GNORM, BUT IN NO CASE
!  CAN IT BE MORE THAN 100*SELCON OR 0.1
!
      if (gs) selcon = min(min(scfcrt*100.D0,0.1D0),max(scfcrt*gnorm**3*10**(-id*3),scfcrt))
      pltest = 0.05D0*sqrt(abs(selcon))
!
!  SOMETIMES HEAT GOES SCF BUT DENSITY IS STILL FLUCTUATING IN UHF
!  IN WHICH CASE PAY LESS ATTENTION TO DENSITY MATRIX
!
      if (nalpha/=nbeta .and. uhf) pltest = 0.001D0
      if (debug) write (iw, '(''  SELCON, PLTEST'',3G16.7)') selcon, pltest
      if (prt1el) then
        write (iw, '(2/10X,''ONE-ELECTRON MATRIX AT ENTRANCE TO ITER'')')
        call vecprt (h, norbs)
      end if
      iredy = 1
  180 continue
      niter = 0
      frst = .TRUE.
      if (camkin) then
        modea = 1
        modeb = 1
      else
        modea = 0
        modeb = 0
      end if
      bfrst = .TRUE.
!*********************************************************************
!                                                                    *
!                                                                    *
!                START THE SCF LOOP HERE                             *
!                                                                    *
!                                                                    *
!*********************************************************************
  250 continue

      incitr = modea/=3 .and. modeb/=3
      if (incitr) niter = niter + 1
      if (timitr) then
        titer = seconds(1)
        write (iw, *)
        if (niter > 1) write (iw, '(a,f9.2,a,/)') &
          '     TIME FOR ITERATION:', titer - titer0, ' WALL CLOCK SECONDS'
        titer0 = titer
      end if
      if (niter > itrmax - 10 .and. .not.allcon) then
!***********************************************************************
!                                                                      *
!                   SWITCH ON ALL CONVERGERS                           *
!                                                                      *
!***********************************************************************
        okpuly = .true.
        camkin = .not. halfe
          if (itrmax > 2) write (iw, &
      '(2/,'' ALL CONVERGERS ARE NOW FORCED ON'',/,'' SHIFT=10,&
      & PULAY ON, CAMP-KING ON'',/,'' AND ITERATION COUNTER RESET'',2/)')
        allcon = .TRUE.
        bshift = 4.44D0
        iredy = -4
        eold = 100.D0
        newdg = .FALSE.
        if (is_PARAM .and. l_param) then
          write (ifiles_1,'(a)') "ALL CONVERGERS ARE NOW FORCED ON"
          l_param = .false.
        end if
        go to 180
      end if
!***********************************************************************
!                                                                      *
!                        MAKE THE ALPHA FOCK MATRIX                    *
!                                                                      *
!***********************************************************************
      if (abs(shift)>1.D-10 .and. bshift/=0.D0) then
        l = 0
        if (niter > 1) then
          if (newdg .and. .not.(halfe .or. camkin)) then
!
!  SHIFT WILL APPLY TO THE VIRTUAL ENERGY LEVELS USED IN THE
!  PSEUDODIAGONALIIZATION. IF DIFF IS -VE, GOOD, THEN LOWER THE
!  HOMO-LUMO GAP BY 0.1EV, OTHERWISE INCREASE IT.
            if (diff > 0.D0) then
              shift = 1.D0
!
! IF THE PSEUDODIAGONALIZATION APPROXIMATION -- THAT THE WAVEFUNCTION
! IS ALMOST STABLE -- IS INVALID, TURN OFF NEWDG
              if (diff > 1) newdg = .FALSE.
            else
              shift = -0.1D0
            end if
          else
            if (ihomo < norbs) then
              shift = ten + eigs(ihomo+1) - eigs(ihomo) + shift
            else
              shift = 0.D0
            end if
          end if
          if (diff > 0.D0) then
            if (shift > 4.D0) shfmax = 4.5D0
            if (shift > shfmax) shfmax = max(shfmax - 0.5D0,0.D0)
          end if
!
!   IF SYSTEM GOES UNSTABLE, LIMIT SHIFT TO THE RANGE -INFINITY - SHFMAX
!   BUT IF SYSTEM IS STABLE, LIMIT SHIFT TO THE RANGE -INFINITY - +20
!
          shift = max(-20.D0,min(shfmax,shift))
          if (abs(shift - shfmax) < 1.D-5) shfmax = shfmax + 0.01D0
!
!  THE CAMP-KING AND PULAY CONVERGES NEED A CONSTANT SHIFT.
!  IF THE SHIFT IS ALLOWED TO VARY, THESE CONVERGERS WILL NOT
!  WORK PROPERLY.
!
          if (okpuly .or. abs(bshift-4.44D0)<1.D-5) then
            shift = -8.D0
            if (newdg) shift = 0.D0
          end if
          if (uhf) then
            if (newdg .and. .not.(halfe .or. camkin)) then
              shiftb = ten - tenold
            else
              shiftb = ten + eigb(ihomob+1) - eigb(ihomob) + shiftb
            end if
            if (diff > 0.D0) shiftb = min(4.D0,shiftb)
            shiftb = max(-20.D0,min(shfmax,shiftb))
            if (okpuly .or. abs(bshift-4.44D0)<1.D-5) then
              shiftb = -8.D0
              if (newdg) shiftb = 0.D0
            end if
            eigb(ihomob+1:norbs) = eigb(ihomob+1:norbs) + shiftb
          end if
        end if
        tenold = ten
        if (pl > plchek) then
          shftbo = shiftb
          shfto = shift
        else
          shiftb = shftbo
          shift = shfto
        end if
        if (id == 0) eigs(ihomo+1:norbs) = eigs(ihomo+1:norbs) + shift
        if (id /= 0) shift = -80.D0
        if (lxfac) shift=0.d0
        forall (i=1:mpack)
          f(i) = h(i) + shift*pa(i) + 1.D-16*i
        endforall

        do i=1,norbs
           f(i*(i+1)/2) = f(i*(i+1)/2) - shift
        end do
      else if (last==0 .and. niter<2 .and. fulscf) then
!
!  SLIGHTLY PERTURB THE FOCK MATRIX IN CASE THE SYSTEM IS
!  TRAPPED IN A METASTABLE EXCITED ELECTRONIC STATE
!
        random = 0.001D0
        glow = method_indo .or. glow .or. (gnorm < 2.D0 .and. gnorm > 1.D-9)
        if (glow) random = 0.D0
        do i = 1, mpack
          random = -random   ! GBR: This sounds strange. Could Random variable be placed out of the loop?
          f(i) = h(i) + random
        end do
      else
        call dcopy(mpack,h,1,f,1)
      end if
  320 continue
      if (timitr) call timer ('BEFORE FOCKS')
      if (id /= 0) then
        call fock2 (f, p, pa, w, w, wk, numat, nfirst, nlast, 2)
      else if (method_indo) then
        if (.not. allocated(dd)) then
          allocate(dd(mpack,2))
          allocate(cc0(norbs,norbs))
          allocate(aa(mpack))
          allocate(ff(mpack,2))
          allocate(dtmp(mpack,2))
          nb2 = mpack
! First round - fill dd based on occupied fractions
          ofract = dble(nclose)*2/(dble(nclose) + dble(nopen))
          do i=1,mpack
            dd(i,1) = p(i) * ofract
            dd(i,2) = p(i) * (1.D0 - ofract)
          end do
!          write (6,*) 'ofract',ofract,nclose,nopen
          opendd = .False.
        else
! Following rounds - fill dd based on orbitals if open shell
!          write (0,*) 'dd sum', c
          if (nopen > nclose .and. pl < pltest*10 .and. abs(diff) < sellim*10) then
            opendd = .True.
          end if
          if (opendd) then
            l = 0
            nl2 = 1
            nu2 = nclose
            nl1 = nclose + 1
            nu1 = nopen
!            write (0,*) 'nl',nl2,nu2,nl1,nu1,nclose,nopen
            do i = 1, norbs
              do j = 1, i
                l = l + 1
                sum2 = sum(c(i,nl2:nu2)*c(j,nl2:nu2))
                sum1 = sum(c(i,nl1:nu1)*c(j,nl1:nu1))
                dd(l,1) = sum2*2.D0
                dd(l,2) = sum1
              end do
            end do
          else
            do i=1,mpack
              dd(i,1) = p(i) - dd(i,2)
            end do
            do i=1,mpack
              dd(i,1) = p(i)
              dd(i,2) = 0.D0
            end do
          end if
        end if
!        do i = 1,6
!          write (0,*) dd(i,1) + dd(i,2), p(i), dd(i,1) + dd(i,2) - p(i)
!        end do
!        write (6,*) 'p',p, pa, pb
        do i=1,norbs
          do j = 1,norbs
            cc0(j,i) = c(i,j)
          end do
        end do

        call scfmat(tot,shiftb)

        do i=1,mpack
! Use full Fock matrix (same for closed-shell, diff for open-shell)
          f(i)  = aa(i)
          dtmp(i,1) = dd(i,1)
          dtmp(i,2) = dd(i,2)
!          p(i) = dd(i,1) + dd(i,2)
!          pa(i) = p(i)/2.D0
!          pb(i) = p(i)/2.D0
        end do

        do i=1,norbs
          do j = 1,norbs
            c(j,i) = cc0(i,j)
          end do
        end do
! dd and p contain the full density matrix (sum of alpha and beta electrons)
        if (useps) then
          call addfck (f,p)
        end if
      else
        call fock2 (f, p, pa, w, w, w, numat, nfirst, nlast, 2)
      end if
      if (lxfac) then
        ee = helect(norbs,pa,h,f)*2.d0
        if (ci .or. halfe) then
          call eigenvectors_LAPACK(c, f, eigs, norbs)
          summ = meci()
          ee = ee + summ
        end if
        return
      end if
      if (timitr) call timer ('AFTER  FOCKS ')
      if (prtfok) then
        if (uhf) write (iw, "('   ALPHA FOCK MATRIX ON ITERATION',i3)") niter
        if ( .not. uhf) write (iw, "('   FOCK MATRIX ON ITERATION',i3)") niter
        call vecprt (f, norbs)
      end if
!***********************************************************************
!                                                                      *
!                        MAKE THE BETA FOCK MATRIX                     *
!                                                                      *
!***********************************************************************
      if (uhf) then
        if (shiftb /= 0.D0) then
          l = 0
          do i = 1, norbs
            if (i > 0) then
              fb(l+1:i+l) = h(l+1:i+l) + shiftb*pb(l+1:i+l)
              l = i + l
            end if
            fb(l) = fb(l) - shiftb
          end do
        else if (rand .and. last==0 .and. niter<2 .and. fulscf) then
          random = 0.001D0
          if (glow) random = 0.D0
          do i = 1, mpack
            random = -random
            fb(i) = h(i) + random
          end do
        else
           call dcopy(mpack,h,1,fb,1)
        end if
        if (id /= 0) then
          call fock2 (fb, p, pb, w, w, wk, numat, nfirst, nlast, 2)
        else
          call fock2 (fb, p, pb, w, w, w, numat, nfirst, nlast, 2)
        end if
        if (prtfok) then
          write (iw, "('   BETA FOCK MATRIX ON ITERATION',i3)") niter
          call vecprt (fb, norbs)
        end if
      end if
      if (.not.fulscf) go to 600
!
!   CODE THE FOLLOWING LINE IN PROPERLY SOMETIME
!   THIS OPERATION IS BELIEVED TO GIVE RISE TO A BETTER FOCK MATRIX
!   THAN THE CONVENTIONAL GUESS.
!
      if (irrr == 0) then
        do i = 1, norbs
          f((i*(i+1))/2) = f((i*(i+1))/2)*0.5D0
        end do
      end if
      irrr = 2
!***********************************************************************
!                                                                      *
!                        CALCULATE THE ENERGY IN KCAL/MOLE             *
!                                                                      *
!***********************************************************************
      if (niter >= itrmax) then
        if (diff < 1.D-3 .and. pl < 1.D-4 .and. .not. force) then
          if (abs(shift) < 1.D-10) write (iw, &
      '('' """""""""""""""UNABLE TO ACHIEVE'',           '' SELF-CONSISTENCE, J&
      &OB CONTINUING'')')
          incitr = .TRUE.
          getout = .TRUE.
          go to 410
        end if
        if (minprt) write (iw, 390)
  390   format(/,/,10x,'"""""""""""""UNABLE TO ','ACHIEVE SELF-CONSISTENCE',/)
        write (iw, 400) diff, pl
  400   format(/,/,10x,'DELTAE= ',e12.4,5x,'DELTAP= ',e12.4,/,/,/)
        iflepo = 9
        iscf = 2
        call writmo
        call mopend ('UNABLE TO ACHIEVE SELF-CONSISTENCE')
        return
      end if
      ee = helect(norbs,pa,h,f)
      if (uhf) then
        ee = ee + helect(norbs,pb,h,fb)
      else
        ee = ee*2.D0
      end if
      if (capps) ee = ee + capcor(nat,nfirst,nlast,p,h)
      if (uhf) then
        if (bshift /= 0.D0) then
          if (nalpha_open > nalpha) then
            scorr = shift*(nalpha_open - nalpha)*enrgy*0.5d0*(fract*(1.D0 - fract))
          else
            scorr = shift*(nbeta_open - nbeta)*enrgy*0.5d0*(fract*(1.D0 - fract))
          end if
        end if
      else
        if (bshift /= 0.D0) scorr = shift*(nopen - nclose)*enrgy*0.25D0*(fract*(2.D0 - fract))
      end if
      escf = (ee + enuclr)*enrgy + atheat + scorr
      getout = .FALSE.
  410 continue
      if (incitr) then
        if (getout) go to 470
        diff = escf - eold
        if (diff > 0) then
          ten = ten - 1.D0
        else
          ten = ten*0.975D0 + 0.05D0
        end if
        sellim = max(selcon,1.d-15*max(abs(ee),1.D0))
!
! SCF TEST:  CHANGE IN HEAT OF FORMATION IN KCAL/MOL SHOULD BE
!            LESS THAN SELLIM.  THE OTHER TESTS ARE SAFETY MEASURES
!
        if (.not.(niter > 4 .and. (pl == 0.D0 .or. &
        pl < pltest .and. abs(diff) < sellim) .and. ready)) go to 490
!***********************************************************************
!                                                                      *
!          SELF-CONSISTENCY TEST, EXIT MODE FROM ITERATIONS            *
!                                                                      *
!***********************************************************************
  470   continue
        if (abs(shift) < 1.D-10) go to 600
        shift = 0.D0
        shiftb = 0.D0
        f(:mpack) = h(:mpack)
        makea = .TRUE.
        makeb = .TRUE.
        go to 320
  490   continue
        if (limscf .and. emin/=0.D0 .and. .not.(ci .or. halfe)) then
!
!  THE FOLLOWING TESTS ARE INTENDED TO ALLOW A FAST EXIT FROM ITER
!  IF THE RESULT IS 'GOOD ENOUGH' FOR THE CURRENT STEP IN THE GEOMETRY
!  OPTIMIZATION
!
          if (escf < emin) then
!
!  THE ENERGY IS LOWER THAN THE PREVIOUS MINIMUM.  NOW CHECK THAT
!  IT IS CONSISTENTLY LOWER.
!
            iemax = 0
            iemin = min(5,iemin + 1)
            escf0(:iemin-1) = escf0(2:iemin)
            escf0(iemin) = escf
!
!  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 5%
!  OF THE ENERGY GAIN FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
!  MINIMUM.
!
            if (iemin > 3) then
              do i = 2, iemin
                if (abs(escf0(i)-escf0(i-1)) > 0.05D0*(emin - escf)) go to 540
              end do
!
! IS GOOD ENOUGH -- RAPID EXIT
!
              if (debug) write (iw, *) &
                ' RAPID EXIT BECAUSE ENERGY IS CONSISTENTLY LOWER'
              incitr = .TRUE.
              getout = .TRUE.
              go to 410
            end if
          else
!
!  THE ENERGY HAS RISEN ABOVE THAT OF THE PREVIOUS MINIMUM.
!  WE NEED TO CHECK WHETHER THIS IS A FLUKE OR IS THIS REALLY
!  A BAD GEOMETRY.
!
            iemin = 0
            iemax = min(5,iemax + 1)
            escf0(:iemax-1) = escf0(2:iemax)
            escf0(iemax) = escf
!
!  IS THE DIFFERENCE IN ENERGY BETWEEN TWO ITERATIONS LESS THAN 5%
!  OF THE ENERGY LOST FOR THIS GEOMETRY RELATIVE TO THE PREVIOUS
!  MINIMUM.
!
            if (iemax > 3) then
              do i = 2, iemax
                if (abs(escf0(i)-escf0(i-1)) > 0.05D0*(escf - emin)) go to 540
              end do
!
! IS GOOD ENOUGH -- RAPID EXIT
!
              if (debug) write (iw, *) &
                ' RAPID EXIT BECAUSE ENERGY IS CONSISTENTLY HIGHER'
              incitr = .TRUE.
              getout = .TRUE.
              go to 410
            end if
          end if
        end if
  540   continue
        ready = iredy>0 .and. (abs(diff)<sellim*10.D0 .or. pl==0.D0)
        iredy = iredy + 1
      end if
      if (prtpl .or. debug .and. niter > itrmax - 20) then
        if (escf > 999999.D0) then
          escf = 999999.D0
        end if
        if (escf < -999999.D0) escf = -999999.D0
        if (abs(diff) > 9999.D0) diff = 0.D0
        if (incitr) then
          write (line,'('' ITERATION'',I4,'' PLS='',2E10.3,'' ENERGY  '',F14.6,'' DELTAE'',F13.7,f14.3)') &
            niter, pl, plb, escf, diff
          write(iw,'(a)')trim(line)
          call to_screen(line)
          endfile (iw)
          backspace (iw)
        end if
      end if
      if (incitr) eold = escf
!***********************************************************************
!                                                                      *
!                        INVOKE THE CAMP-KING CONVERGER                *
!                                                                      *
!***********************************************************************
      if (niter>2 .and. camkin .and. makea) then
        if (.not. Allocated (vec_ai)) then
        allocate (vec_ai(norbs, norbs), fock_ai(norbs, norbs), &
        & p_ai(norbs, norbs), h_ai(norbs**2), vecl_ai(norbs**2), stat=i)
        if (i /= 0) then
          call memory_error("Camp-King converger in Iter")
          return
        end if
! TODo: make it parallel
        vec_ai(:,:)= 0.0d0
        fock_ai(:,:) = 0.0d0
        p_ai(:,:) = 0.0d0
        h_ai(:) = 0.0d0
        vecl_ai(:) = 0.0d0
      end if
        call interp (na1el, norbs - na1el, modea, escf/enrgy, f, c, &
        theta, vec_ai, fock_ai, p_ai, h_ai, vecl_ai, eold_alpha)
      end if
      makeb = .FALSE.
      if (modea /= 3) then
        makeb = .TRUE.
        if (debug) then
          write (iw, *) ' Diagonal of FOCK Matrix'
          write (iw, '(8F10.6)') (f((i*(i + 1))/2),i = 1,norbs)
        end if
        if (nscf == 2) then
          continue
          end if
        if (newdg) then
!***********************************************************************
!                                                                      *
!                        INVOKE PULAY'S CONVERGER                      *
!                                                                      *
!***********************************************************************
           if (okpuly .and. makea .and. iredy>1) then
#ifdef GPU
              if (lgpu) then
                 call pulay_for_gpu (f, pa, norbs, pold, pold2, pold3, &
                 & jalp, ialp, npulay*mpack, frst, pl)
              else
#endif
                 call pulay (f, pa, norbs, pold, pold2, pold3, &
                 & jalp, ialp, npulay*mpack, frst, pl)
#ifdef GPU
              end if
#endif
          end if

!***********************************************************************
!                                                                      *
!           DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT           *
! WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S  *
!                                                                      *
!***********************************************************************
! JEM NOTE: pseudo-diagonalization does less work, but it is less stable
!           and slower because a large fraction of its work is stuck
!           using BLAS level-1 operations right now. This could be
!           accelerated using a QR diagonalization analog at some point ...
!          if (halfe .or. camkin) then
            if (timitr) call timer ('BEFORE FULL DIAG')
            call eigenvectors_LAPACK(c, f, eigs, norbs)
            if (timitr) call timer ('AFTER  FULL DIAG')
!          else
!            if (lgpu) then
!               if (timitr) call timer ('BEFORE GPU DIAG')
!               call diag_for_GPU (f, c, na1el, eigs, norbs, mpack)
!               if (timitr) call timer ('AFTER  GPU DIAG')
!            else
!              if (timitr) call timer ('BEFORE CPU DIAG')
!               call diag_for_GPU (f, c, na1el, eigs, norbs, mpack)
!               if (timitr) call timer ('AFTER  CPU DIAG')
!            end if
!          end if
        else
          if (timitr) call timer ('BEFORE FULL DIAG')
          call eigenvectors_LAPACK(c, f, eigs, norbs)
          if (timitr) call timer ('AFTER  FULL DIAG')
        end if
        j = 1
        if (prtvec) then
            j = 1
            if (uhf) j = 2
            write (iw, &
         & '(2/10X,A,'' EIGENVECTORS AND EIGENVALUES ON ITERATION'',I3)') abprt(j), niter
            call matout (c, eigs, norbs, norbs, norbs)
        else
          if (prteig) write (iw, 550) abprt(j), niter, (eigs(i),i=1,norbs)
        end if
  550   format(10x,a,'  EIGENVALUES ON ITERATION',i3,/,10(6g13.6,/))
      end if
      if (ifill /= 0) call swap (c, norbs, norbs, na2el, ifill)
!***********************************************************************
!                                                                      *
!            CALCULATE THE ALPHA OR RHF DENSITY MATRIX                 *
!                                                                      *
!***********************************************************************
        if (timitr) call timer ('BEFORE DENSIT')
        if (uhf) then
          call density_for_GPU (c, fract, nalpha, nalpha_open, 1.d0, mpack,norbs, 1, pa, iopc_calcp)
          if (modea /= 3 .and. .not. (newdg .and. okpuly)) then
            i = niter
            if (camkin) i = 7
            call cnvg (pa, pold, pold2,  i, pl)
          end if
        else
          if (halfe) then
            call densit (c, norbs, norbs, na2el, 2.d0, na1el, fract, p, 1)
          else
            call density_for_GPU (c, fract, na2el, na1el, 2.d0, mpack, norbs, 1, p, iopc_calcp)
          end if
          if (modea/=3 .and. .not.(newdg .and. okpuly)) then
            call cnvg (p, pold, pold2,  niter, pl)
        end if
      end if
      if (timitr) call timer ('AFTER  DENSIT')
!***********************************************************************
!                                                                      *
!                       UHF-SPECIFIC CODE                              *
!                                                                      *
!***********************************************************************
      if (uhf) then
!***********************************************************************
!                                                                      *
!                        INVOKE THE CAMP-KING CONVERGER                *
!                                                                      *
!***********************************************************************
        if (niter > 2 .and. camkin .and. makeb) then
          if (.not. Allocated (vec_bi)) then
          allocate (vec_bi(norbs, norbs), fock_bi(norbs, norbs), &
               & p_bi(norbs, norbs), h_bi(norbs**2), vecl_bi(norbs**2), stat=i)
          if (i /= 0) then
            call memory_error("Camp-King converger in Iter")
            return
          end if
        end if
          call interp (nb1el, norbs - nb1el, modeb, escf/enrgy, fb, cb, &
          theta, vec_bi, fock_bi, p_bi, h_bi, vecl_bi, eold_beta)
        end if
        makea = .FALSE.
        if (modeb /= 3) then
          makea = .TRUE.
          if (newdg) then
!***********************************************************************
!                                                                      *
!                        INVOKE PULAY'S CONVERGER                      *
!                                                                      *
!***********************************************************************
            if (okpuly .and. makeb .and. iredy>1) then
#ifdef GPU
              if (lgpu) then
                 call pulay_for_gpu (fb, pb, norbs, pbold, pbold2, pbold3, &
                 & jbet, ibet, npulay*mpack, bfrst, plb)
              else
#endif
                call pulay (fb, pb, norbs, pbold, pbold2, pbold3, &
                 & jbet, ibet, npulay*mpack, bfrst, plb)
#ifdef GPU
              end if
#endif
          end if

!***********************************************************************
!                                                                      *
!           DIAGONALIZE THE ALPHA OR RHF SECULAR DETERMINANT           *
! WHERE POSSIBLE, USE THE PULAY-STEWART METHOD, OTHERWISE USE BEPPU'S  *
!                                                                      *
!***********************************************************************
! JEM NOTE: pseudo-diagonalization does less work, but it is less stable
!           and slower because a large fraction of its work is stuck
!           using BLAS level-1 operations right now. This could be
!           accelerated using a QR diagonalization analog at some point ...
!            if (halfe .or. camkin) then
              if (timitr) call timer ('BEFORE FULL DIAG')
              call eigenvectors_LAPACK(cb, fb, eigb, norbs)
              if (timitr) call timer ('AFTER  FULL DIAG')
!            else
!              if (lgpu) then
!                 if (timitr) call timer ('BEFORE GPU DIAG')
!                 call diag_for_GPU (fb, cb, nb1el, eigb, norbs, mpack)
!                 if (timitr) call timer ('AFTER  GPU DIAG')
!              else
!                if (timitr) call timer ('BEFORE CPU DIAG')
!                call diag_for_GPU (fb, cb, nb1el, eigb, norbs, mpack)
!                if (timitr) call timer ('AFTER  CPU DIAG')
!              end if
!            end if
          else
            if (timitr) call timer ('BEFORE FULL DIAG')
            call eigenvectors_LAPACK(cb, fb, eigb, norbs)
            if (timitr) call timer ('AFTER  FULL DIAG')
          end if

          if (prtvec) then
            write (iw, &
      '(2/10X,A,'' EIGENVECTORS AND EIGENVALUES ON '',   ''ITERATION'',I3)') &
            &  abprt(3), niter
            call matout (cb, eigb, norbs, norbs, norbs)
          else
            if (prteig) write (iw, 550) abprt(3), niter, (eigb(i),i=1,norbs)
          end if

        end if
!***********************************************************************
!                                                                      *
!                CALCULATE THE BETA DENSITY MATRIX                     *
!                                                                      *
!***********************************************************************
        if (timitr) call timer ('BEFORE B-DENS')
        call density_for_GPU (cb, fract, nbeta, nbeta_open, 1.d0, mpack, norbs, 1, pb, iopc_calcp)
        if (.not.(newdg .and. okpuly)) then
          i = niter
          if (camkin) i = 7
          call cnvg (pb, pbold, pbold2, i, plb)
        end if
        if (timitr) call timer ('AFTER  B-DENS')
      end if
!***********************************************************************
!                                                                      *
!                   CALCULATE THE TOTAL DENSITY MATRIX                 *
!                                                                      *
!***********************************************************************
      if (uhf) then
        forall (i=1:mpack)
           p(i) = pa(i) + pb(i)
        endforall
      else
        forall (i=1:mpack)
           pa(i) = p(i)*0.5d0
           pb(i) = pa(i)
        endforall
      end if
      if (debug) then
        call chrge (p, q)
        write (iw, *) ' CHARGES'
        write (iw, '(8F10.7)') (q(i),i=1,numat)
      end if
      if (prtden) then
        write (iw, '('' DENSITY MATRIX ON ITERATION'',I4)') niter
        call vecprt (p, norbs)
      end if
      if (itrmax < 3) return
      oknewd = pl<sellim .or. oknewd
      newdg = pl<trans .and. oknewd .or. newdg
      if (pl < trans*0.3333D0) oknewd = .TRUE.
      go to 250
!*********************************************************************
!                                                                    *
!                                                                    *
!                      END THE SCF LOOP HERE                         *
!                NOW CALCULATE THE ELECTRONIC ENERGY                 *
!                                                                    *
!                                                                    *
!*********************************************************************
!          SELF-CONSISTENCE ACHIEVED.
!
  600 continue
      ee = helect(norbs,pa,h,f)
      if (uhf) then
        ee = ee + helect(norbs,pb,h,fb)
      else
        ee = ee*2.D0
      end if
      if (capps) ee = ee + capcor(nat,nfirst,nlast,p,h)
      if (timitr) call timer ('BEFORE FINAL DIAG')
      if (nscf==0 .or. last==1 .or. ci .or. halfe) then
!
!  PUT F AND FB INTO POLD IN ORDER TO NOT DESTROY F AND FB
!  AND DO EXACT DIAGONALISATIONS
!
        call dcopy(mpack,f,1,pold,1)
        call eigenvectors_LAPACK(c, pold, eigs, norbs)
        if (last == 1) call phase_lock(c, norbs)
        if (uhf) then
           call dcopy(mpack,fb,1,pold,1)
           call eigenvectors_LAPACK(cb, pold, eigb, norbs)
           if (last == 1) call phase_lock(cb, norbs)
           call dcopy(mpack,pa,1,pold,1)
        else
           call dcopy(mpack,p,1,pold,1)
        end if
        if (ci .or. halfe) then
          if (timitr) call timer ('BEFORE MECI')
          summ = meci()
          if (timitr) call timer ('AFTER  MECI')
          if (moperr) return
          ee = ee + summ
          if (prtpl) then
            escf = (ee + enuclr)*enrgy + atheat
            write (iw, '(27X,''AFTER MECI, ENERGY  '',F14.7)') escf
          end if
        end if
      end if
      if (timitr) call timer ('AFTER  FINAL DIAG')
      nscf = nscf + 1
      if (debug) write (iw, '('' NO. OF ITERATIONS ='',I6)') niter
      if (allcon .and. abs(bshift - 4.44D0) < 1.D-7) then
        camkin = .FALSE.
        allcon = .FALSE.
        newdg = .FALSE.
        bshift = -10.D0
        okpuly = .FALSE.
      end if
      shift = 1.D0
      escf = (ee + enuclr)*enrgy + atheat
      if (emin == 0.D0) then
        emin = escf
      else
        emin = min(emin,escf)
      end if
      return
      end subroutine iter
      subroutine delete_iter_arrays
!
!  If these arrays exist, they must be deleted before starting ITER.
!  If they're not deleted, and they're too small, the stack will
!  become corrupt if INTERP is called.
!
      use iter_C, only : vec_ai, vec_bi, fock_ai, fock_bi, p_ai, &
      p_bi, h_ai, h_bi, vecl_ai, vecl_bi
      implicit none
        if (allocated(vec_ai)) deallocate(vec_ai)
        if (allocated(vec_bi)) deallocate(vec_bi)
        if (allocated(fock_ai)) deallocate(fock_ai)
        if (allocated(fock_bi)) deallocate(fock_bi)
        if (allocated(p_ai)) deallocate(p_ai)
        if (allocated(p_bi)) deallocate(p_bi)
        if (allocated(h_ai)) deallocate(h_ai)
        if (allocated(h_bi)) deallocate(h_bi)
        if (allocated(vecl_ai)) deallocate(vecl_ai)
        if (allocated(vecl_bi)) deallocate(vecl_bi)
      end subroutine delete_iter_arrays


      subroutine den_in_out(mode)
      use chanel_C, only : iw, iden, density_fn
      use common_arrays_C, only : p, pa, pb
      use molkst_C, only: uhf, keywrd, norbs, numat, mozyme
!
! Read and write the density matrix
!
!   If mode == 0 read the density matrix in (fill p, pa, and, if uhf, pb)
!   if mode == 1 write out the density matrix
!
      implicit none
      integer, intent (in) :: mode
      integer :: io_stat, icount, old_norbs, old_numat
      logical :: formatted, opend
      if (Index (keywrd, " DENOUT") == 0 .and. mode == 1) return
      if (Index (keywrd, " OLDENS") == 0 .and. mode == 0) return
      call l_control("NEWDEN", len("NEWDEN"), 1)
      if (mozyme) then
        call pinout(mode, .true.)
        return
      end if
      formatted = (Index(keywrd," DENOUTF") /= 0)
      inquire(unit=iden, opened=opend)
      if (opend) close(unit=iden, status='KEEP')
      do icount = 1,2
        if (formatted)then
          open(unit=iden, file=density_fn)
        else
          open(unit=iden, file=density_fn, status='UNKNOWN', &
          form='UNFORMATTED', position='asis')
        end if
        rewind iden
        if (mode == 0) then
!
!  Try reading the file using formatted input if DENOUTF exists.
!  If the read does not work, try the other way
!
          if (formatted)then
             read (iden, *, iostat=io_stat)old_norbs, old_numat, pa
          else
             read (iden, iostat=io_stat)old_norbs, old_numat, pa
          end if
          if (old_norbs > 0 .and. old_norbs < 100000) then   ! Delete this conditional in 2008
          if (norbs /= old_norbs .or. numat /= old_numat) then
            call mopend("Density file read in does not match current data set")
            return
          end if
          end if
          if (icount < 2.and. io_stat /= 0) then
            formatted = (.not. formatted)
            close(iden)
            cycle
          end if
          if (io_stat /= 0) then
            call to_screen(" Density Restart File missing or corrupt")
            call mopend ('Density Restart File missing or corrupt')
            return
          end if
            if (uhf) then
              if (formatted)then
                read (iden, *, iostat=io_stat) pb
              else
                read (iden, iostat=io_stat) pb
              end if
            if (io_stat /= 0) then
              if (index(keywrd, " GEO-OK") == 0) then
                call mopend("Beta Density Restart File missing or corrupt")
                write(iw,'(10x,a)')'(Most likely the previous job did not use UHF)'
                write(iw,'(10x,a)')"To continue, using the RHF density matrix as the starting point, add keyword ""GEO-OK"""
                return
              end if
              pb = pa
            end if
            p = pa + pb
          else
            p = pa*2.d0
          end if
        else
          if (formatted)then
            write (iden, *, iostat=io_stat)norbs, numat, pa
          else
            write (iden, iostat=io_stat)norbs, numat, pa
          end if
          if (uhf) then
            if (formatted)then
              write (iden, *, iostat=io_stat) pb
            else
              write (iden, iostat=io_stat) pb
            end if
          end if
          close (iden)
          exit
        end if
        close (iden)
      end do
      end subroutine den_in_out
