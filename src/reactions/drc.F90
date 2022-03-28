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

     subroutine drc(startv, startk)
      USE molkst_C, only : nvar, tleft, escf, gnorm, ndep, numat, keywrd, &
      line, nopen, nclose, moperr, jloop => itemp_1, prt_velocity
      USE chanel_C, only : iw0, iw, ires, iscr, restart_fn
      use common_arrays_C, only : geo, loc, grad, xparam, atmass, nat, &
      & errfn, coord, na, p, q, labels, geoa, txtatm
      USE parameters_C, only : tore
      USE elemts_C, only : elemnt
      use drc_C, only : georef
      implicit none
      double precision :: startv(9*numat*numat)
      double precision , intent(in) :: startk(3*numat)
!
      integer :: i, j, iskin,   iloop, k, kl, ncoprt, bigcycles, jloop_lim, i_constant = 1, &
        ii, i1, maxcyc, iupper, ilp, io_stat, ilim, iw00, percent, n_escf, n_min, iwd = 9
      integer, dimension(2,3*numat) :: mcoprt
      double precision, dimension(3*numat) :: velo0, velo1, velo2, velo3, gerror, grold2
      double precision, dimension(10) :: past10
      double precision, dimension(3*numat) :: grold
      double precision :: ekin, elost1, etold, dlold2, tnow, oldtim, delold, gtot, &
        accu, gnlim, half, addonk, deltat, quadr, etot, const, one, summ, &
        summas, ams, error, velvec, delta1, elost, sum, dummy, tcycle, &
        old_hof, new_hof, start_hof=1.d20, minstep, escf_old, &
        escf_diff, stepx, escf_min, stepxx, dipvec(3), dip = 0.d0, cnvg
      logical :: addk, letot, let, velred, opend, parmax, debug, l_debug = .false., l_irc, &
        l_dipole, l_pdb
      double precision, external :: ddot, dipole, reada, seconds
!***********************************************************************
!                                                                      *
!    DRC IS DESIGNED TO FOLLOW A REACTION PATH FROM THE TRANSITION     *
!    STATE.  TWO MODES ARE SUPPORTED, FIRST: GAS PHASE:- AS THE SYSTEM *
!    MOVES FROM THE T/S THE MOMENTUM OF THE ATOMS IS STORED AND THE    *
!    POSITION OF THE ATOMS IS RELATED TO THE OLD POSITION BY (A) THE   *
!    CURRENT VELOCITY OF THE ATOM, AND (B) THE FORCES ACTING ON THAT   *
!    ATOM.  THE SECOND MODE IS CONDENSED PHASE, IN WHICH THE ATOMS MOVE*
!    IN RESPONSE TO THE FORCES ACTING ON THEM. I.E. INFINITELY DAMPED  *
!                                                                      *
!***********************************************************************
      addk = .TRUE.
      ekin = 0.d0
      escf_old = 10.d0
      escf_diff = 20.d0
      elost1 = 0.d0
      etold = 0.d0
      dlold2 = 0.d0
      past10 = 0.d0
      old_hof = 0.d0
      new_hof = 0.d0
      tnow = seconds(1)
      oldtim = seconds(1)
      delold = 10.d0
      percent = 0
      gtot = 0.d0
      gnorm = 0.d0
      iw00 = iw0
      iw0 = -1
      iloop = 1
      n_escf = 0
      escf_min = 1.d20
      parmax = .FALSE.
      l_pdb = (txtatm(1)(:7) /= "")
      call l_control("LDRC_FIRST", len("LDRC_FIRST"), 1)
      if (nopen /= nclose .and. Index (keywrd, " IRC") /= 0) then
        minstep = 5.d-16 ! Gradients are not very accurate, so use larger minimum step.
      else
        minstep = 1.d-16
      end if
      if (allocated(grad))  deallocate(grad)
      if (allocated(errfn)) deallocate(errfn)
      if (allocated(georef)) deallocate(georef)
      allocate(grad(3*numat), errfn(3*numat), georef(3,numat))
      errfn = 0.d0
      inquire(unit=iscr, opened=opend)
      if (opend) close(unit=iscr)
      open(unit=iscr, status='SCRATCH', position='asis')
      if (index(keywrd,' PREC') /= 0) then
        accu = 0.25d0
        cnvg = 0.000001d0
      else
        accu = 1.d0
        cnvg = 0.00001d0
      end if
      l_dipole = (index(keywrd,' DIPOLE') /= 0)
      l_irc = (index(keywrd,' DRC') == 0)
      if (l_irc) then
        if (index(keywrd,' X-PRIORITY=') /= 0) then
          stepx = reada(keywrd,index(keywrd,'X-PRIO') + 5)
        else
          stepx = 0.05D0
        end if
        stepx = stepx*0.2d0  !  Use five steps per point
        stepxx = 0.01d0
      else
        stepx = 0.d0
        stepxx = 0.d0
      end if
!
!   COSMO surfaces are relatively rough, so increase number of steps to
!   ensure on the descent side
!
      if (Index (keywrd, " EPS") /= 0) then
        ilim = 300
        accu = accu*10
      else
        ilim = 30
      end if
      gnlim = 1.d0
      past10(5) = 100.d0
      debug = (index(keywrd,' DEBUG') /= 0)
      i = index(keywrd,'GNORM')
      if (i /= 0) gnlim = reada(keywrd,i)
      if (gnlim < 0.9d0) then
        n_min = 80
      else
        n_min = 40
      end if
!
! DETERMINE EXCESS KINETIC ENERGY
!
      iskin = 0
      if (index(keywrd,'KINE') /= 0) then
        iskin = 1
        addonk = reada(keywrd,index(keywrd,'KINE'))
        if (addonk < 0.d0) startv = -startv
        write (iw,'(2/10X,'' EXCESS KINETIC ENERGY ENTERED INTO SYSTEM ='',F12.6)') Abs(addonk)
        if (addonk < 0.d0) then
          addonk = -addonk
          write (iw,'(10X,'' KINETIC ENERGY SUPPLIED WAS NEGATIVE, SO INITIAL VELOCITY IS REVERSED'')')
        end if
      else
        addonk = 0.d0
        if (l_irc) addonk=2.d0  ! Just a guess to get the IRC started
      end if
      velred = index(keywrd,'VELO') /= 0
      if (ddot(3*numat,startv,1,startv,1) > 0.001d0) then
!
!     PRINT OUT INITIAL VELOCITIES
!
        if (prt_velocity) then
          write (iw, '(10x, A)') ' INITIAL VELOCITY IN DRC (Angstroms/Femtosecond)'
          write (iw, '(3F16.5)') (startv(i), i = 1, numat*3)
        end if
        startv(:numat*3) = -startv(:numat*3)
      end if
      let = (velred .and. Index (keywrd, " IRC") == 0)
      if (index(keywrd,' SYMM') /= 0) ndep = 0
   !
   !      CONVERT TO CARTESIAN COORDINATES, IF NOT ALREADY DONE.
   !
      call gmetry (geo, coord)
      geo(:,:numat) = coord(:,:numat)
      na = 0
      do i = 1, numat
        if (atmass(i) >= 1.d-1) cycle
        write (iw, '(A,I3,A)') ' ATOMIC MASS OF ATOM', i, ' TOO SMALL'
        return
      end do
!
!  TRANSFER COORDINATES TO XPARAM AND LOC
!
      if (.not. l_irc) then
        parmax = (loc(1,1) /= 0)
        ncoprt = 0
        if (parmax) then
          mcoprt(:,:nvar) = loc(:,:nvar)
          ncoprt = nvar
        end if
      else
        ncoprt = 0
      end if
      ncoprt = 0  !  Do NOT print turning points. These just mess up the output.
      if (allocated(geoa)) deallocate(geoa)
      allocate( geoa(3,numat))
      geoa(:,:numat) = geo(:,:numat)
      georef = geo
      if (index(keywrd, " NOOPT") == 0) then
        k = 0
        do i = 1, numat
          do j = 1, 3
            k = k + 1
            loc(1,k) = i
            loc(2,k) = j
          end do
        end do
        nvar = k
      end if


      do i = 1, nvar
        xparam(i) = geo(loc(2,i), loc(1,i))
      end do
!
! DETERMINE DAMPING FACTOR
!
      if (index(keywrd,'DRC=') /= 0) then
        half = reada(keywrd,index(keywrd,'DRC='))
        write (iw, &
      '(2/10X,'' DAMPING FACTOR FOR KINETIC ENERGY ='',F12.6)') half
      else if (index(keywrd,' DRC') == 0) then
!
!  IRC
!
        half = 0.d0
!   startv = startv*0.02d0
        if (addonk > 1.d-4) startv(:3*numat) = startv(:3*numat)*addonk
      else
!
! undamped DRC
!
        half = 1.d6
      end if
!
!  LETOT IS TRUE IF CORRECTIONS ARE NOT TO BE MADE PART WAY INTO
!        THE CALCULATION
!
!  USAGE OF LETOT:
! (1) WHILE LETOT IS FALSE, NO DAMPING WILL BE DONE
! (2) WHEN LETOT IS TURNED TRUE,
!     IF AN IRC, THEN ETOT IS RESET SO THE ERROR IS ZERO.
!     IF A  DRC, EXCESS KINETIC ENERGY USED TO START THE RUN IS REMOVED.
!
      letot = index(keywrd,'IRC=') == 0 .and. .not.let
      half = sign(max(0.000001d0,abs(half)),half)

!
!   LOOP OVER TIME-INTERVALS OF DELTAT SECOND
!
      deltat = 1.d-16
      if (index(keywrd,'IRC=') /= 0) deltat = min(1.d-13, deltat*(numat*0.25d0)**4)  !  Big step if starting from a transition state
      quadr = 1.d0
      etot = 0.d0
      escf = 0.d0
      const = 1.d0
      if (index(keywrd,'RESTART') /= 0 .and. index(keywrd,'IRC=') == 0) then
!
!  RESTART FROM A PREVIOUS RUN
!
        inquire(unit=ires, opened=opend)
        if (opend) close(unit=ires, status='KEEP')
        open(unit=ires, file=restart_fn, status='UNKNOWN', &
          form='UNFORMATTED', position='asis', iostat = io_stat)
        if (io_stat /= 0) then
          write(iw,*)" Restart file either does not exist or is not available for reading"
          call mopend ("Restart file either does not exist or is not available for reading")
          return
        end if
        rewind ires
        read (ires, iostat = io_stat) (xparam(i),i=1,nvar)
        read (ires, iostat = io_stat) (velo0(i),i=1,nvar)
        read (ires, iostat = io_stat) (grad(i),i=1,nvar)
        read (ires, iostat = io_stat) (grold(i),i=1,nvar)
        read (ires, iostat = io_stat) (grold2(i),i=1,nvar)
        read (ires, iostat = io_stat) etot, escf, ekin, delold, deltat, dlold2, iloop, gnorm, &
            letot, elost1, gtot
        close (ires)
        if (io_stat /= 0) then
          call mopend ("Restart file is corrupt")
          return
        end if
        if ((abs(etot) > 1.d5 .or. abs(etot) < 1.d-5) .and. (abs(escf) > 1.d5 .or. abs(escf) < 1.d-5)) then
          call mopend ("RESTART FILE DOES NOT APPEARS TO BE AN IRC OR DRC RESTART. ")
          write(iw,'(10x,a)')"COULD IT BE A FORCE RESTART, IF SO USE IRC=n"
          return
        end if
        write (iw, &
      '(2/10X,''CALCULATION RESTARTED, CURRENT KINETIC ENERGY='',F10.5,2/)') ekin
        go to 100
      else
!                         NOT A RESTART

        velo0(:nvar) = 0.d0
        grold2(:nvar) = 0.d0
        grold(:nvar) = 0.d0
        grad(:nvar) = 0.d0
        if (index(keywrd,'IRC=') /= 0 .or. velred) then
!
!  GET HOLD OF VELOCITY VECTOR
!
          if (index(keywrd,'IRC=') /= 0) then
            k = nint(reada(keywrd,index(keywrd,'IRC=')))
          else
            k = 1
          end if
          if (k < 0) then
            k = -k
            one = -1.d0
          else
            one = 1.d0
          end if
!
!   Modify all velocities to set net momentum to zero
!
          kl = (k - 1)*nvar
          summ = 0.d0
          velo1(1) = 0.d0
          velo1(2) = 0.d0
          velo1(3) = 0.d0
          summas = 0.d0
          i = 0
          do ii = 1, numat
            ams = atmass(ii)
            summas = summas + ams
            velo0(i+1) = startv(kl+i+1)*one
            velo1(1) = velo1(1) + velo0(i+1)*ams
            velo0(i+2) = startv(kl+i+2)*one
            velo1(2) = velo1(2) + velo0(i+2)*ams
            velo0(i+3) = startv(kl+i+3)*one
            velo1(3) = velo1(3) + velo0(i+3)*ams
            i = i + 3
          end do
          velo1(:3) = -velo1(:3)/summas
          i = 0
          if (addonk > 1.d-5 .or. .not.velred) then
            do ii = 1, numat
              ams = atmass(ii)
              do i1 = 1, 3
                i = i + 1
                velo0(i) = velo0(i) + velo1(i1)
                summ = summ + velo0(i)**2*ams
              end do
            end do
          else
            do ii = 1, numat
              ams = atmass(ii)
              do i1 = 1, 3
                i = i + 1
                summ = summ + velo0(i)**2*ams
              end do
            end do
          end if

          if (addonk < 1.d-5 .and. velred) addonk = 0.5d0*summ/4.184D10
          if (addonk < 1.d-5 .and. .not. velred) then
            if (abs(half) > 1.d-3 .and. startk(k) > 105.d0) then
              write (iw, '(A,F10.3,A,/,A)') &
                ' BY DEFAULT, ONE QUANTUM OF ENERGY, EQUIVALENT TO', &
                startk(k), ' CM(-1)', ' WILL BE USED TO START THE DRC'
!
!    2.8585086D-3 CONVERTS CM(-1) INTO KCAL/MOLE
!
              addonk = startk(k)*2.8585086D-3
              write (iw, '(A,F7.2,A)') ' THIS REPRESENTS AN ENERGY OF', &
              addonk, ' KCALS/MOLE'
            else if (abs(half) > 1.d-3) then
              write (iw, '(A,F9.2,A)') ' THE VIBRATIONAL FREQUENCY (', &
              startk(k), 'CM(-1)) IS TOO SMALL', ' FOR ONE QUANTUM TO BE USED'
              write (iw, '(A)') &
                ' INSTEAD 0.3KCAL/MOLE WILL BE USED TO START THE IRC'
              addonk = 0.3d0
            else
              addonk = 0.3d0
            end if
          end if
!
!   AT THIS POINT ADDONK IS IN KCAL/MOLE
!   NORMALIZE SO THAT TOTAL K.E. = ONE QUANTUM (DEFAULT) (DRC ONLY)
!                              OR 0.3KCAL/MOLE (IRC ONLY)
!                              OR ADDONK IF KINETIC=NN SUPPLIED
!
          if (summ < 1.d-4) then
            write (iw, '(A)') ' SYSTEM IS APPARENTLY NOT MOVING!'
            return
          end if
!
!  ADDONK IS EXCESS KINETIC ENERGY.  IF THE CALCULATION IS AN IRC,
!  THIS ENERGY MUST BE REMOVED AFTER A SHORT 'TIME'.
!
!  MAKE AN AD-HOC CORRECTION: IF ADDONK IS NON-ZERO AND HALF IS LARGER
!  THAN 0.1, MODIFY ADDONK TO REFLECT ERRORS DUE TO START-UP.
!
          if (half > 0.1d0 .and. half < 10000.d0) &
            addonk = addonk*(1.d0 + 0.06972d0/half)
!
!  MAKE AN AD-HOC CORRECTION: IF ADDONK IS NON-ZERO AND HALF IS LESS
!  THAN -0.1, MODIFY ADDONK TO REFLECT ERRORS DUE TO START-UP.
!
          if (half<(-0.1d0) .and. half>(-10000.d0)) addonk = addonk*(1.d0 + &
            0.06886d0/half)
          summ = sqrt(addonk/(0.5d0*summ/4.184D10))
          addk = .FALSE.
          if (summ > 1.d-10) then
            velo0(:nvar) = velo0(:nvar)*summ
!
!  IF IT IS A DRC, DESTROY ADDONK.  THE KINETIC ENERGY USED WILL COME
!  FROM THE VELOCITY ONLY.
!
            if (half > 1.d-3) addonk = 0.d0
          end if
        end if
      end if
  100 continue
      if (Index (keywrd, " BIGCYCLES") /= 0) then
        bigcycles = Nint (reada (keywrd, Index (keywrd, " BIGCYCLES")))*4
      else
        bigcycles = -1
      end if


      if (index(keywrd,' CYCLES') /= 0) then
        maxcyc = 1000000
        jloop_lim = nint(reada(keywrd,index(keywrd,' CYCLES')))
      else
        maxcyc = 4999
        jloop_lim = 100000
      end if
      iupper = iloop + maxcyc
      ilp = iloop
      one = 0.d0
      if (index(keywrd,'RESTART') /= 0 .and. index(keywrd,'IRC=') == 0) one = 1.d0
      gerror(:nvar) = 0.d0
      if (l_debug) then
!
! Use this block to print diagnostic details of the reaction coordinate
!
        line = "debug.txt"
        call add_path(line)
        open(iwd, file = trim(line))
        write(iwd,'(a)')" Point  Velocity      Coordinate      Delta-x    Gradient     ESCF   Deltat*(10^16)"
      end if
      do iloop = ilp, iupper
!
!  MOVEMENT OF ATOMS WILL BE PROPORTIONAL TO THE AVERAGE VELOCITIES
!  OF THE ATOMS BEFORE AND AFTER TIME INTERVAL
!
!
!  RAPID CHANGE IN GRADIENT IMPLIES SMALL STEP SIZE FOR DELTAT
!
!   KINETIC ENERGY = 1/2 * M * V * V
!                  = 0.5 / (4.184D10) * M * V * V
!   NEW VELOCITY = OLD VELOCITY + GRADIENT * TIME / MASS
!                = KCAL/ANGSTROM*SECOND/(ATOMIC WEIGHT)
!                =4.184*10**10(ERGS)*10**8(PER CM)*DELTAT(SECONDS)
!   NEW POSITION = OLD POSITION - AVERAGE VELOCITY * TIME INTERVAL
!
!
!   ESTABLISH REFERENCE TOTAL ENERGY
!
        error = etot - (ekin + escf)
        if (iloop > 2) then
          quadr = 1.d0 + error/(ekin*const + 0.001d0)*0.5d0
          quadr = min(1.3d0,max(0.8d0,quadr))
        else
          quadr = 1.d0
        end if
        if ((let .or. ekin > 0.2d0) .and. addk) then
!
!   DUMP IN EXCESS KINETIC ENERGY
!
          etot = etot + addonk
          addk = .FALSE.
          addonk = 0.d0
        end if
!
!  CALCULATE THE DURATION OF THE NEXT STEP.
!  STEP SIZE IS THAT REQUIRED TO PRODUCE A CONSTANT CHANGE IN GEOMETRY
!
!
!  IF DAMPING IS USED, CALCULATE THE NEW TOTAL ENERGY AND
!  THE RATIO FOR REDUCING THE KINETIC ENERGY
!
        const = max(1.d-36,0.5d0**(deltat*1.d15/half))
        const = sqrt(const)
        ekin = 0.d0
        delta1 = delold + dlold2
        elost = 0.d0
        startv(:nvar) = xparam(:nvar)
        if (iloop > 3) then
          do i = 1, nvar
!
!   CALCULATE COMPONENTS OF VELOCITY AS
!   V = V(0) + V'*T + V"*T*T
!   WE NEED ALL THREE TERMS, V(0), V' AND V"
!
            velo1(i) = 1.d0/atmass(loc(1,i))*grad(i)   ! Velocity
            velo3(i) = 2.d0/atmass(loc(1,i))* &
                (delta1*(grold(i)-grad(i)) - &
                delold*(grold2(i)-grad(i)))/(delta1*(delold**2*1.d30) - &
                delold*(delta1**2*1.d30))
            velo2(i) = 1.d0/atmass(loc(1,i))* &
            (grad(i)-grold(i)-0.5d0*velo3(i)*(1.d30*delold**2))/(delold*1.d15)
          end do
          if (l_irc) then
            velo2 = 0.d0
            velo3 = 0.d0
          end if
          do
            do i = 1, nvar
!
!  MOVE ATOMS THROUGH DISTANCE EQUAL TO VELOCITY * DELTA-TIME, NOTE
!  VELOCITY CHANGES FROM START TO FINISH, THEREFORE AVERAGE.
!
              startv(i) = xparam(i) - 1.d8*( &
              deltat*velo0(i)*one + &
              0.5d0*deltat**2*velo1(i) + &
              0.16666d0*(deltat**2*1.d15)*deltat*velo2(i) + &
              0.0416666d0*deltat**2*(1.d30*deltat**2)*velo3(i))
            end do
            if (stepxx < 1.d-5 .or. (l_irc .and. iloop > 4)) exit
            sum = 1.d-5
            do i = 1, nvar
              sum = sum + (xparam(i) - startv(i))**2
            end do
            sum = stepxx/sqrt(sum)
            if (sum > 0.8d0 .and. sum < 1.25d0) exit
            deltat = deltat*max(0.99d0, min(1.01d0, sum))
          end do
          sum = 0.d0
          do i = 1, numat
            do j = 1,3
              sum = sum + (geo(j,i) - georef(j,i))**2
            end do
          end do
          if (sum > 0.01d0) gnorm = 4.d0
          if (l_irc) velo0 = 0.d0
          do i = 1, nvar
!
!   CORRECT ERRORS DUE TO CUBIC COMPONENTS IN ENERGY GRADIENT,
!   ALSO TO ADD ON EXCESS ENERGY, IF NECESSARY.
!
!
!   MODIFY VELOCITY IN LIGHT OF CURRENT ENERGY GRADIENTS.
!
!   VELOCITY = OLD VELOCITY + (DELTA-T / ATOMIC MASS) * CURRENT GRADIENT
!                           + 1/2 *(DELTA-T * DELTA-T /ATOMIC MASS) *
!                             (SLOPE OF GRADIENT)
!              SLOPE OF GRADIENT = (GRAD(I)-GROLD(I))/DELOLD
!
!
!   THIS EXPRESSION IS ACCURATE TO SECOND ORDER IN TIME.
!
            velo0(i) =                              velo0(i) + &          ! Velocity
                                             deltat*velo1(i) + &          ! Acceleration
                                    0.5d0*deltat**2*velo2(i)*1.d15 + &    ! Rate of acceleration
                0.166666d0*deltat*(1.d30*deltat**2)*velo3(i)              ! Rate of rate of acceleration
            if (let .or. gnorm > 3.d0) then
              let = .TRUE.
              elost = elost + velo0(i)**2*atmass(loc(1,i))*(1 - const**2)
              velo0(i) = velo0(i)*const*quadr
            end if
!
!  CALCULATE KINETIC ENERGY (IN 2*ERGS AT THIS POINT)
!
            ekin = ekin + velo0(i)**2*atmass(loc(1,i))
          end do
          if (l_debug) then
!
! Use this block to print diagnostic details of the reaction coordinate
!
            if (iloop > 25) then
              write(iwd,'(i6, f10.2, f15.6, f15.6, f10.4, 2f12.6)')iloop, velo0(42)/(const*quadr), &
              startv(42), startv(42) - xparam(42), grad(42)/4.184D18, escf, deltat*1.d16
            end if
          end if
          xparam(:nvar) = startv(:nvar)
        else
!
!   First few steps do not have velo3
!
!
          if (l_irc) then
            velo1 = 0.d0
            velo2 = 0.d0
            velo3 = 0.d0
          end if
          do i = 1, nvar
!
!   CALCULATE COMPONENTS OF VELOCITY AS
!   V = V(0) + V'*T + V"*T*T
!   WE NEED ALL THREE TERMS, V(0), V' AND V"
!
            velo1(i) = 1.d0/atmass(loc(1,i))*grad(i)
            velo2(i) = 1.d0/atmass(loc(1,i))*(grad(i)-grold(i))/(1.d15*delold)
            velo3(i) = 0.d0
!
!  MOVE ATOMS THROUGH DISTANCE EQUAL TO VELOCITY * DELTA-TIME, NOTE
!  VELOCITY CHANGES FROM START TO FINISH, THEREFORE AVERAGE.
!
            xparam(i) = xparam(i) - 1.d8*(deltat*velo0(i)*one + &
            0.5d0*deltat**2*velo1(i) + &
            0.16666d0*(deltat**2*1.d15)*deltat*velo2(i) + &
            0.0416666d0*deltat**2*(1.d30*deltat**2)*velo3(i))
!
!   CORRECT ERRORS DUE TO CUBIC COMPONENTS IN ENERGY GRADIENT,
!   ALSO TO ADD ON EXCESS ENERGY, IF NECESSARY.
!
!
!   MODIFY VELOCITY IN LIGHT OF CURRENT ENERGY GRADIENTS.
!
!   VELOCITY = OLD VELOCITY + (DELTA-T / ATOMIC MASS) * CURRENT GRADIENT
!                           + 1/2 *(DELTA-T * DELTA-T /ATOMIC MASS) *
!                             (SLOPE OF GRADIENT)
!              SLOPE OF GRADIENT = (GRAD(I)-GROLD(I))/DELOLD
!
!
!   THIS EXPRESSION IS ACCURATE TO SECOND ORDER IN TIME.
!
            velo0(i) = velo0(i) + deltat*velo1(i) + 0.5d0*deltat**2*velo2(i)*&
              1.d15 + 0.166666d0*deltat*(1.d30*deltat**2)*velo3(i)
            if (let .or. gnorm>3.d0) then
              let = .TRUE.
              elost = elost + velo0(i)**2*atmass(loc(1,i))*(1 - const**2)
              velo0(i) = velo0(i)*const*quadr
            end if
!
!  CALCULATE KINETIC ENERGY (IN 2*ERGS AT THIS POINT)
!
            ekin = ekin + velo0(i)**2*atmass(loc(1,i))
          end do
        end if
        one = 1.d0
        if (let .or. gnorm > 3.d0) then
          if (.not.letot) then
            if (abs(half) < 1.d-3) then
!
!  IT IS AN IRC, SO RESET THE TOTAL ENERGY
!
              etot = escf + elost1
              addonk = 0.d0
              elost1 = 0.d0
              elost = 0.d0
              deltat = 5.d-16
            else if (iskin == 0) then
!
!  IT IS A DRC AND KINETIC NOT USED, SO REMOVE EXTRA KINETIC ENERGY
!
              etot = etot - addonk
            end if
          end if
          letot = .TRUE.
        end if
!
!  CONVERT ENERGY INTO KCAL/MOLE
!
        ekin = 0.5d0*ekin/4.184D10
!
!  IF IT IS A DAMPED DRC, MODIFY ETOT TO REFLECT LOSS OF KINETIC ENERGY
!
        if (letot .and. abs(half) > 0.00001d0) etot = etot - ekin/const**2 + ekin
        elost1 = elost1 + 0.5d0*elost/4.184D10
!
! STORE OLD GRADIENTS FOR DELTA - VELOCITY CALCULATION
!
        grold2(:nvar) = grold(:nvar)
        grold(:nvar) = grad(:nvar)
        grad(:nvar) = 0.d0
!
!   CALCULATE ENERGY AND GRADIENTS
!
        if (Abs(new_hof) > 1.d-20) old_hof = escf
        call compfg (xparam, .TRUE., escf, .TRUE., grad, .TRUE.)
        if (Abs(new_hof) < 1.d-20) then
        old_hof = escf + 5.d0
        new_hof = escf
        end if
        if (moperr) return
        if (debug) then
          write(iw,"(//,a, i8)")" Point calculated:    ",iloop
          write(iw,"(a, f13.5)")" Heat of formation:    ", escf
          if (start_hof > 1.d19) start_hof = escf
          write(iw,"(a, f10.5)")" Calculated energy change:", elost1
          write(iw,"(a, f11.5)")" Predicted energy change:", etot - escf
          if (Abs(elost1) > 1.d-4) write(iw,"(a, f11.5)")" Cumulative error:       ", escf + elost1 - etot

          write(iw,"(a,  f13.4)")" 'Time' interval (fs):",deltat*1.d15

          write(iw,"(a)")" Geometry supplied to COMPFG in DRC"
          do i = 1, numat
            if (l_pdb) then
              write(iw,'(i6, 2x, a,3F15.6)') i, elemnt(nat(i))//"("//txtatm(i)(:26)//")", (xparam((i - 1)*3 + j), j = 1,3)
            else
              write(iw,"(1x,a2,f15.4,2f17.4)")elemnt(nat(i)), (xparam((i - 1)*3 + j), j = 1,3)
            end if
          end do
          write(iw,"(/,a)")" Forces (gradients) acting on atoms"
          do i = 1, numat
            if (l_pdb) then
              write(iw,'(i6, 2x, a,3F15.6)') i, elemnt(nat(i))//"("//txtatm(i)(:26)//")", (grad((i - 1)*3 + j), j = 1,3)
            else
              write(iw,"(i6, 1x,a2,f15.4,2f17.4)")i, elemnt(nat(i)), (grad((i - 1)*3 + j), j = 1,3)
            end if
          end do
          velvec = 0.d0
          do i = 1, nvar
            velvec = velvec + velo0(i)**2
          end do
          if (velvec > 1.d0) then
            write(iw,"(/,a)")" Velocity of atoms, in cm/sec"
            do i = 1, numat
              if (l_pdb) then
                write(iw,'(i6, 2x, a,3F15.6)') i, elemnt(nat(i))//"("//txtatm(i)(:26)//")", (-velo0((i - 1)*3 + j), j = 1,3)
              else
                write(iw,"(i6, 1x,a2,f15.4,2f17.4)")i, elemnt(nat(i)), (-velo0((i - 1)*3 + j), j = 1,3)
              end if
            end do
          else
            write(iw,"(/,a)")" Acceleration of atoms, in 10^(-14)cm/(sec*sec) = 10^10 Angstroms/(fs**2)"
            do i = 1, numat
              write(iw,"(i6, 1x,a2,f15.4,2f17.4)")i, elemnt(nat(i)),&
              (-1.d-14*velo1((i - 1)*3 + j), j = 1,3)
            end do
          end if
        end if
        new_hof = escf
        if (iloop > 2) then
          gnorm = 0.d0
          do i = 1, nvar, 3
            sum = dsqrt(ddot(3,grad(i),1,grad(i),1)/(ddot(3,velo0(i),1,velo0(i),1)+1.d-20))
            gerror(i:2+i) = gerror(i:2+i) + grad(i:2+i) + velo0(i:2+i)*sum
          end do
          gnorm = dsqrt(ddot(nvar,gerror,1,gerror,1))
          gtot = gnorm
        end if
        gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
!
!   CONVERT GRADIENTS INTO ERGS/CM
!
        grad(:nvar) = grad(:nvar)*4.184D18
!
!   SPECIAL TREATMENT FOR FIRST POINT - SET "OLD" GRADIENTS EQUAL TO
!   CURRENT GRADIENTS.
!
        if (iloop == 1) then
          grold(:nvar) = grad(:nvar)
        end if
        dlold2 = delold
        delold = deltat
        sum = 0.d0
        do i = 1, nvar
          sum = sum + ((grad(i)-grold(i))/4.184D18)**2
        end do
        if (abs(half) < 0.001d0) then
          deltat = deltat*min(2.d0,2.d-4*accu/(abs(escf + elost1 - etold) + 1.d-20))**0.25d0
          etold = escf + elost1
          if (iloop > ilim .and. addonk < 1.d-5 .and. &
          & abs(old_hof - new_hof) < cnvg .or. iloop > 10000) then
!
!  Deciding when to terminate the IRC is tricky.
            i_constant = i_constant + 1
            if (i_constant > 3 .or. iloop > 10000) then
              iw0 = iw00
              write (iw, '(2/,'' IRC CALCULATION COMPLETE '')')
              return
            end if
          else
            i_constant = 1
          end if
          if (escf < escf_min .or. iloop < 10) then
            escf_min = escf
            n_escf = 0
          else
            n_escf = n_escf + 1
            if (n_escf > n_min) then
              write (iw, '(2/10x,A,/)') &
            ' POTENTIAL ENERGY HAS STOPPED DROPPING'
              iw0 = iw00
              return
            end if
          end if
        else
          deltat = deltat*min(1.05d0,10.d0*accu/(sum + 1.d-4))
          deltat = min(deltat,3.d-15*accu)
          past10(10) = gnorm
          sum = 0.d0
          do i = 1, 9
            sum = sum + abs(past10(i)-past10(i+1))
            past10(i) = past10(i+1)
          end do
          if (sum < gnlim) then
            write (iw, '(2/,A)') &
              ' GRADIENT CONSTANT AND SMALL -- ASSUME ALL MOTION STOPPED'
            write (iw, '(A)') ' TO CONTINUE, USE KEYWORD ''GNORM=0'''
            iw0 = iw00
            return
          end if
          deltat = min(deltat,2.d-15)
!***********************************************************************
!
!         TESTING CODE - REMOVE BEFORE FINAL VERSION ASSEMBLED
!#          (ILOOP/400)*400.EQ.ILOOP)DELTAT=-DELTAT
!
!***********************************************************************
        end if
        deltat = max(minstep,deltat)
        if (l_dipole) then
          call chrge (p, q)
          q(:numat) = tore(labels(:numat)) - q(:numat)
          dip = dipole(p, xparam, dipvec,0)
        end if
        iw0 = iw00
        if (abs(half) < 0.00001d0) then
!
!   FOR THE IRC:
!
! ESCF   = POTENTIAL ENERGY
! ELOST1 = ENERGY LOST (IN DRC, THIS WOULD HAVE BEEN THE KINETIC ENERGY)
! ETOT   = COMPUTED TOTAL ENERGY = STARTING POTENTIAL ENERGY
!
!   IN DRCOUT  'TOTAL' = ESCF + ELOST1
!              'ERROR' = ESCF + ELOST1 - ETOT
!
      call prtdrc (deltat, xparam, georef, escf, elost1, gtot, etot, velo0, mcoprt, ncoprt, parmax, &
          l_dipole, dip)
        else
!
!   FOR THE DRC:
!
! ESCF   = POTENTIAL ENERGY
! EKIN   = CURRENT KINETIC ENERGY
! ETOT   = COMPUTED TOTAL ENERGY = STARTING POTENTIAL ENERGY -
!          KINETIC ENERGY LOST THROUGH DAMPING, IF PRESENT.
!
!   IN DRCOUT  'TOTAL' = ESCF + EKIN
!              'ERROR' = ESCF + EKIN - ETOT
!
          call prtdrc (deltat, xparam, georef, escf, ekin, dummy, etot, velo0, mcoprt, ncoprt, parmax, &
          l_dipole, dip)
          if (iloop > 10 .and. (escf - escf_old)/ escf_diff < 0.d0) then
            bigcycles = bigcycles - 1
            if (bigcycles == 0) iupper = iloop
          end if
          escf_diff = escf - escf_old
          escf_old = escf
        end if
        iw0 = -1
        tnow = seconds(2)
        tcycle = tnow - oldtim
        oldtim = tnow
        tleft = tleft - tcycle
        if (iw00 > -1) then
          i = nint((100.0*iloop)/(iupper))
          if (i /= percent) then
            percent = i
            write(line,"(i4,a)")percent, "% of DRC/IRC done"
            call to_screen(line)
          end if
        end if
        if (jloop < jloop_lim .and. iloop /= iupper .and. tleft >= 3*tcycle) cycle
        if (tleft < 3*tcycle) then
          inquire(unit=ires, opened=opend)
          if (opend) close(unit=ires, status='DELETE')
          open(unit=ires, file=restart_fn, status='UNKNOWN', form='UNFORMATTED', position='asis')
          rewind ires
          write (ires) 1,1, (xparam(i),i=1,nvar)
          write (ires) (velo0(i),i=1,nvar)
          write (ires) (grad(i),i=1,nvar)
          write (ires) (grold(i),i=1,nvar)
          write (ires) (grold2(i),i=1,nvar)
          i = iloop + 1
          write (ires) etot, escf, ekin, delold, deltat, dlold2, i, gnorm, letot, elost1, gtot
          close (ires)
          escf = -1.d9
          call prtdrc (deltat, xparam, georef, escf, ekin, elost, etot, velo0, mcoprt, ncoprt, parmax, &
          l_dipole, dip)
          call den_in_out(1)
        end if
        do j = 1, 3
          do i = 1, numat
            geo(j, i) = coord(j, i)
            coord(j, i) = 0.0d0
          end do
        end do
        na = 0
        if ((jloop < jloop_lim .or. iloop == iupper) .and. bigcycles < 0) then
          write (iw, '(/10X,'' NUMBER OF CYCLES EXCEEDED, RESTART FILE WRITTEN'')')
          write (iw, '(10x," Number of cycles allowed in this run =",i7)') maxcyc + 1
          write (iw, '(10x," To increase the number of cycles, use keyword ""CYCLES=n"",")')
          write (iw, '(10x," where ""n"" is the number of cycles you want to run")')
        else
          write (iw, '(/10a,'' RUNNING OUT OF TIME, RESTART FILE WRITTEN'')')
        end if
        iw0 = iw00
        return
      end do
      iw0 = iw00
      return
      end subroutine drc
