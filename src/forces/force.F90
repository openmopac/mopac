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

      subroutine force()
!
      use molkst_C, only : natoms, ndep,  nvar, gnorm, iflepo, keywrd, &
      & last, numat, escf, id, jloop => itemp_1, numcal, n_trivial => itemp_2, &
      moperr, this_point => itemp_2, zpe, mozyme, uhf, prt_force, prt_normal_coords, &
      prt_orientation, maxtxt, l_normal_html
!
      use common_arrays_C, only : xparam, na, nb, nc, geo, geoa, ca => c, cb, eigs, eigb, &
      & labels, coord, loc, grad, errfn, na_store, nat, lopt, fmatrx, p, q, txtatm
!
      use symmetry_C, only : name, igroup
!
      use parameters_C, only : tore
      USE elemts_C, only : elemnt
      USE funcon_C, only : fpc_10, fpc_6, fpc_8, a0, ev, fpc_9
      use to_screen_C, only : dipt, travel, freq, redmas, cnorml, force_const
      USE chanel_C, only : iw, ixyz
      implicit none
!
      integer , dimension(60) :: irot
      integer , dimension(2,3*natoms) :: locold
      integer , dimension(natoms) :: nar, nbr, ncr
      integer :: j, i, l, nvaold, ndeold, iu, il, nvib, ij, &
        im1, ju, jl, ii, jj, ni, k, nto6, nrem6, iinc1, iinc2, store_natoms, store_numat
      double precision, dimension(3,3*numat) :: deldip, trdip
      double precision, dimension(3,3) :: rot
      double precision :: time2, tscf, tder, time1, time3, a, b, c, &
        sum, const, summ, sum1, sym
      double precision, dimension(:), allocatable :: store, ff, oldf, &
        velocity
      double precision, dimension (:,:), allocatable :: store_coord
      logical :: restrt, linear, debug, prnt, large, ts
      double precision, external :: ddot, dipole, reada, seconds
!**********************************************************************
!
!   FORCE CALCULATES THE FORCE CONSTANTS FOR THE MOLECULE, AND THE
!         VIBRATIONAL FREQUENCIES.  ISOTOPIC SUBSTITUTION IS ALLOWED.
!
!**********************************************************************
!  Point-Groups are stored in order of increasing symmetry.  Groups
!   supported are, in order
!
! C1, Cs, Ci, C2, D2, C2v, C2h, D2h, C3, C4, S4, D3, C3v, C3h, C5, C6,
! S6, C7, C8, S8, C4v, D4, D2d, C5v, C6v, D6, C4h, D3h, D3d, C7v, C7h,
! C8v, D8, D4d, D4h, C5h, D5, D5h, D5d, C6h, D6h, D6d, D7, D7h, D7d,
! C8h, D8h, T, Td, O, Th, Oh, I, Ih, Cv, Dh, R3
!
      data irot/ 1, 1, 1, 2, 4, 2, 2, 4, 3, 4, 2, 6, 3, 3, 5, 6, 3, 7, 8, 4, 4&
        , 8, 4, 5, 6, 12, 4, 6, 6, 7, 7, 8, 16, 8, 8, 5, 10, 10, 10, 6, 12, 12&
        , 14, 14, 14, 8, 16, 12, 12, 24, 12, 24, 24, 24, 1, 2, 4*1/
!
      k = 0
      nvaold = 0
      rot = 0.d0
      a = 0.d0
      b = 0.d0
      c = 0.d0
!
!
! TEST GEOMETRY TO SEE IF IT IS OPTIMIZED
!
      time2 = -1.D9
! Save the connectivity and geometry for restarts
      nar(:natoms) = na(:natoms)
      nbr(:natoms) = nb(:natoms)
      ncr(:natoms) = nc(:natoms)
      if (allocated(geoa)) deallocate(geoa)
      allocate(geoa(3,natoms))
      geoa = geo
      ts = (index(keywrd, " FORCETS") + index(keywrd, " MINI") /= 0)
      store_natoms = natoms
      if (id == 0 .and. index(keywrd,' NOREOR') == 0 .and. .not. ts) then
!
!   NEED TO ENSURE THAT XYZINT WILL WORK CORRECTLY BEFORE CALL
!   TO DRC.
!
        l = 0
        do i = 1, natoms
          if (labels(i) == 99) cycle
          l = l + 1
          labels(l) = labels(i)
          txtatm(l) = txtatm(i)
        end do
        numat = l
        natoms = numat
        nvib = 0
        call xyzint (coord, numat, na, nb, nc, 1.D0, geo)
        na_store = na
        call gmetry (geo, coord)
        nvaold = nvar
        locold(1,:nvar) = loc(1,:nvar)
        locold(2,:nvar) = loc(2,:nvar)
      end if
!
!   Unconditionally convert structure into Cartesian coordinates
!   because FORCE works in Cartesian coordinates.
!   Pay careful attention to the translation vectors - they must
!   be at the end of the coordinates, and should not be marked for
!   optimization.
!
      if ( .not. ts) nvar = 0
      ndeold = ndep
      ndep = 0
      do i = 1, numat + id
        do j = 1, 3
          geo(j,i) = coord(j,i)
          if ( .not. ts) then
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = j
            xparam(nvar) = geo(j,i)
          end if
          labels(i) = nat(i)
        end do
      end do
      if (ts) then
        do i = 1, nvar
          xparam(i) = geo(loc(2,i), loc(1,i))
        end do
      end if
      if (id > 0) then
        natoms = numat + id
        nvar = nvar - 3*id ! No need to check for TS - solids can't have transition states.
        labels(numat + 1:numat + id) = 107
      end if
      na = 0
      if (ts) then
        i = nvar
      else
        i = 3*(numat + id)
      end if

!
! Deallocate grad and errfn in case the number of variables set in the input file was smaller
! than the number needed in FORCE
!
      if (allocated(grad))            deallocate (grad)
      if (allocated(errfn))           deallocate (errfn)
      if (allocated(dipt))            deallocate (dipt)
      if (allocated(travel))          deallocate (travel)
      if (allocated(force_const))     deallocate (force_const)
      if (allocated(freq))            deallocate (freq)
      if (allocated(redmas))          deallocate (redmas)
      if (allocated(cnorml))          deallocate (cnorml)
      if (allocated(fmatrx))          deallocate (fmatrx)
      allocate(cnorml(i**2), fmatrx((i*(i+1))/2),grad(i), ff(i**2), errfn(i), &
       & store((i*(i+1))/2), oldf((i*(i + 1))/2), dipt(3*numat), &
       travel(3*numat), force_const(3*numat), freq(3*numat), redmas(3*numat,2), stat = j)
       if (j /= 0) then
         write(iw,*)" Failed to allocate memory in FORCE"
         call mopend("Failed to allocate memory in FORCE")
         return
       end if
      errfn = 0.d0
      grad = 0.d0
!
!   IF A RESTART, THEN TSCF AND TDER WILL BE FAULTY, THEREFORE SET TO -1
!
      tscf = -1.D0
      tder = -1.D0
      prnt = index(keywrd,'RC=') == 0
      debug = index(keywrd,'DFORCE') /= 0
      large = index(keywrd,'LARGE') /= 0
      restrt = index(keywrd,'RESTART') /= 0
      time1 = seconds(1)
      if (.not. restrt) then
        call compfg (xparam, .TRUE., escf, .TRUE., grad, .FALSE.)
        if (numat == 1) goto 98
        if (moperr) goto 99
        write (iw,'(2/10X,''HEAT OF FORMATION ='',F15.6,'' KCALS/MOLE'')') escf
        time2 = seconds(1)
        tscf = time2 - time1
        call compfg (xparam, .TRUE., escf, .FALSE., grad, .TRUE.)
        time3 = seconds(1)
        tder = time3 - time2
        if (prnt) then
          if (ts) then
          else
            write (iw, '(2/10X,''CARTESIAN COORDINATE DERIVATIVES'',2/3X,&
        &     ''NUMBER  ATOM'',9X,''X'',12X,''Y'',12x,''Z'',/)')
            l = 0
            iu = 0
            do i = 1, natoms
              if (labels(i) == 99) cycle
              l = l + 1
              il = iu + 1
              iu = il + 2
              if (labels(i) == 107) iu = il
              write (iw, '(I6,6X,A2,F14.6,2F13.6)') l, elemnt(labels(i)), (grad(j),j=il,iu)
            end do
            write(iw,*)
          end if
        end if
        gnorm = dsqrt(ddot(nvar,grad,1,grad,1))
        if (.not. mozyme .and. index(keywrd, " AUX") /= 0) then
          call chrge (p, q)
          do i = 1, numat
            l = nat(i)
            q(i) = tore(l) - q(i)
          end do
          sum = dipole(p,coord,fmatrx,1)
          call symtrz (ca, eigs, 1, .true.)
          if (uhf) call symtrz (cb, eigb, 1, .TRUE.)
          call to_screen("To_file: Normal output")
        end if
        write (iw, '(2/10X,''GRADIENT NORM ='',F10.5)') gnorm
        sum = max(10.d0, sqrt(float(numat))*2.d0)
        if (gnorm > sum .and. .not. ts) then
          if (index(keywrd,' LET ') /= 0) then
             write (iw, &
        '(3/1X,''** GRADIENT IS VERY LARGE, BUT SINCE "LET" IS USED, CALCULATION WILL CONTINUE'')')
          else
            write (iw, &
        '(3/1X,"** GRADIENT IS TOO LARGE TO ALLOW THE FORCE MATRIX TO BE CALCULATED, (LIMIT =",f5.1,") **",2/)') sum
            write (iw, "(/,A)") " THE GEOMETRY IS NOT AT A STATIONARY POINT ON THE POTENTIAL ENERGY SURFACE"
            write (iw, "(/,A)") " EITHER ADD 'LET' OR OPTIMIZE THE GEOMETRY SO THAT THE GEOMETRY IS AT A STATIONARY POINT."
            write (iw, "(A)") " (FOR A TRANSITION STATE, USE 'TS', FOR A GROUND STATE DO A NORMAL GEOMETRY OPTIMIZATION)"
            call to_screen("To_file: ERROR: GRADIENT IS TOO LARGE TO ALLOW FORCE MATRIX TO BE CALCULATED")
            call mopend ("Gradient in FORCE is too large, geometry is not at a stationary point on the P.E.S.")
            goto 99
          end if
        end if
      end if
      sum = max(1.d0, sqrt(float(numat))*0.2d0)
      if (index(keywrd,'THERMO') /= 0 .and. gnorm > sum) then
        write (iw, &
      '(2/30X,''**** WARNING ****'',2/10X,'' GRADIENT IS VERY LARGE FOR A THERMO CALCULATION'', &
      & /10X,'' RESULTS ARE LIKELY TO BE INACCURATE IF THERE ARE'')')
        write (iw, &
          '(10X,'' ANY LOW-LYING VIBRATIONS (LESS THAN ABOUT 400CM-1)'')')
        write (iw, &
      '(10X," GRADIENT NORM SHOULD BE LESS THAN ABOUT",f4.1," FOR THERMO",&
      &/10X,'' TO GIVE ACCURATE RESULTS'')') sum
      end if
      if ( .not. mozyme .and. .not. restrt) call mullik()
      if (tscf > 0.01D0) then
        write (iw, '(2/10X,''TIME FOR SCF CALCULATION ='',F8.2)') tscf
        write (iw, '( /10X,''TIME FOR DERIVATIVES     ='',F8.2)') tder
      end if
      if (ndeold > 0) write (iw, &
      '(2/10X,''SYMMETRY WAS SPECIFIED, BUT CANNOT BE USED HERE'')')
      c = 1.d0
      if ( .not. ts) call axis (a, b, c, rot)
      allocate (store_coord(3,numat))
      store_coord(:,:numat) = coord(:,:numat)
      geo(:,:numat) = coord(:,:numat)
      if (rot(1,1) > 2.d0) goto 99 ! dummy use of rot
      nvib = 3*numat - 6
      if (ts) then
        nvib = nvar
      else
        if (abs(c) < 1.D-20) nvib = nvib + 1
        if (id /= 0) nvib = 3*numat - 3
      end if
      if (prnt .and. prt_orientation) then
        write (iw, &
      '(/10X,''ORIENTATION OF MOLECULE IN FORCE CALCULATION'')')
        write (iw, &
      '(/,4X,''NO.'',7X,''ATOM'',9X,''X'',9X,''Y'',9X,''Z'',/)')
      end if
      l = 0
      if (.not. (prnt .and. prt_orientation)) then
        l = l + count(labels(:natoms)/=99)
      else
        do i = 1, natoms
          if (labels(i) == 99) cycle
          l = l + 1
          write (iw, '(I6,9X,A2,3X,3F10.4)') l, elemnt(labels(i)), (coord(j,l),&
            j=1,3)
        end do
      end if
      call symtrz (cnorml, cnorml, 2, .FALSE.)
      call fmat (fmatrx, nvib, tscf, tder, deldip, escf, ff, ts)
      if (moperr) goto 99
      na(1) = 0
      na(:natoms) = nar(:natoms)
      nb(:natoms) = nbr(:natoms)
      nc(:natoms) = ncr(:natoms)
      geo(:,:natoms) = geoa(:,:natoms)
      if (nvib < 0) then
        ndep = ndeold
        nvar = 0
        iflepo = -1
        goto 99
      end if
!
!   THE FORCE MATRIX IS PRINTED AS AN ATOM-ATOM MATRIX RATHER THAN
!   AS A 3N*3N MATRIX, AS THE 3N MATRIX IS VERY CONFUSING!
!
      ij = 0
      iu = 0
      do i = 1, nvar/3
        il = iu + 1
        iu = il + 2
        im1 = i - 1
        ju = 0
        do j = 1, im1
          jl = ju + 1
          ju = jl + 2
          sum = 0.D0
          do ii = il, iu
            do jj = jl, ju
              sum = sum + fmatrx((ii*(ii-1))/2+jj)**2
            end do
          end do
          ij = ij + 1
          store(ij) = sqrt(sum)
        end do
        ij = ij + 1
        store(ij) = sqrt(fmatrx(((il+0)*(il+1))/2)**2+fmatrx(((il+1)*(il+2))/2)&
          **2+fmatrx(((il+2)*(il+3))/2)**2+2.D0*(fmatrx(((il+1)*(il+2))/2-1)**2&
          +fmatrx(((il+2)*(il+3))/2-2)**2+fmatrx(((il+2)*(il+3))/2-1)**2))
      end do
      if (debug) then
        write (iw, '(2/10X,'' FULL FORCE MATRIX, INVOKED BY "DFORCE"'')')
        if (index(keywrd, " NOREOR") == 0) then
          write(iw,'(/10x,a)')" Caution: NOREOR is NOT present, therefore system will be oriented"
          write(iw,'(10x,a)')" so that the moments of inertia are along the Cartesian axes."
        end if
        i = -nvar
        call vecprt (fmatrx, i)
      end if
      if (prnt .and. prt_force) then
        sum = 1.d-21*fpc_10*a0**2/(4.184*ev*fpc_9)
        write(iw,'(/,a,f7.4,a)')"(To convert to Hartree/Bohr^2, multiply by (10^(-21) x N x a0^2)/(4.184 x 627.51) = ", &
          sum,")"
        write (iw, '(2/10X,'' FORCE MATRIX IN MILLIDYNES/ANGSTROM'')')
        if (ts) then
          i = nvar/3
        else
          i = numat
        end if
        call vecprt (store, i)
      end if
      l = (nvar*(nvar + 1))/2
      store(:l) = fmatrx(:l)
      if (prnt) call axis (a, b, c, rot)
!
!  A molecule is linear IF one or two of a, b, and c are much smaller than the largest of a, b, and c.
!
      linear = abs((a/(a + b + c))*(b/(a + b + c))*(c/(a + b + c))) < 1.D-10
      if (prnt) write (iw, &
        '(2/10X,''HEAT OF FORMATION ='',F15.6,'' KCALS/MOLE'')') escf
      coord(:,:numat) = store_coord(:,:numat)
      if (large) then
        store_numat = numat
        numat = nvar/3
        call frame (store, numat, 0)
        numat = store_numat
        call rsp (store, nvar, freq, cnorml)
        call phase_lock(cnorml, nvar)
        do i = nvib + 1, nvar
          j = int((freq(i)+50.D0)*0.01D0)
          freq(i) = freq(i) - j*100
        end do
        if (prnt) then
          write (iw, '(2/10X,''TRIVIAL VIBRATIONS, SHOULD BE ZERO'')')
          write (iw, &
      '(/, F9.4,''=TX'',F9.4,''=TY'',F9.4,''=TZ'',F9.4,''=RX'',F9.4,''=RY'',F9.4,''=RZ'')') &
      & (freq(i),i=nvib + 1,nvar)
          call symtrz (cnorml, freq, 2, .TRUE.)
          write (iw, '(2/''      MOLECULAR POINT GROUP   :   '',A4)') name
          write (iw, '(2/10X,'' EIGENVECTORS  '')')
          call matou1 (cnorml, freq, nvib, nvar, nvib, 5)
          write (iw, &
      '(2/10X,''FORCE CONSTANTS IN MILLIDYNES/ANGSTROM'' ,'' (= 10**5 DYNES/CM)'',/)')
          write (iw, '(8F10.5)') (freq(i),i=1,nvib)
! CONVERT TO WEIGHTED FMAT
          write (iw, '(2/10X,'' ASSOCIATED EIGENVECTORS'')')
          i = -nvar
          call matout (cnorml, freq, nvib, i, nvar)
        end if
      end if
      n_trivial = nvar - nvib
      call freqcy (fmatrx, freq, travel, force_const, .TRUE., deldip, ff, oldf, ts)
!
!  CALCULATE ZERO POINT ENERGY
!
!.
!   N AVOGADRO'S NUMBER
!   H PLANCK'S CONSTANT
!   C SPEED OF LIGHT
!   CONST=0.5*N*H*C/(1.D10*4.184)
      const = 0.5D0*fpc_10*fpc_6*fpc_8/(1.D10*4.184D0)
      sum = 0.D0
      ni = 0
      do i = 1, nvib
        if (freq(i) > 0) then
          sum = sum + freq(i)
        else
          ni = ni + 1
        end if
      end do
      zpe = sum*const
      if (prnt) then
        write (iw, '(/9X,'' ZERO POINT ENERGY'', F12.3,'' KCAL/MOL'')') zpe
        if (ni /= 0) then
          write (iw, '(2(/9X,A))') &
            'NOTE: SYSTEM IS NOT A GROUND STATE, THEREFORE ZERO POINT', &
            'ENERGY IS NOT MEANINGFULL. ZERO POINT ENERGY PRINTED'
          write (iw, '(9X,A,I3,A)') 'DOES NOT INCLUDE THE', ni, &
            ' IMAGINARY FREQUENCIES'
        end if
      end if
      summ = 0.D0
      do i = 1, nvar
        sum1 = 1.D-20
        do j = 1, nvar
          sum1 = sum1 + cnorml(j+(i-1)*nvar)**2
        end do
        sum1 = 1.D0/sqrt(sum1)
        grad = 0.D0
        do k = 1, 3
          sum = 0.D0
          do j = 1, nvar
            sum = sum + cnorml(j+(i-1)*nvar)*deldip(k,j)
          end do
          summ = summ + abs(sum)
          trdip(k,i) = sum*sum1
        end do
        dipt(i) = sqrt(trdip(1,i)**2+trdip(2,i)**2+trdip(3,i)**2)
      end do
      if (prnt .and. large) then
        write (iw, &
      '(2/10X,'' FREQUENCIES, REDUCED MASSES AND VIBRATIONAL DIPOLES''/)')
        nto6 = nvar/6
        nrem6 = nvar - nto6*6
        iinc1 = -5
        if (nto6 >= 1) then
          do i = 1, nto6
            write (iw, '(/)')
            iinc1 = iinc1 + 6
            iinc2 = iinc1 + 5
            write (iw, '(3X,''I'',10I10)') (j,j=iinc1,iinc2)
            write (iw, '('' FREQ(I)'',6F10.4,/)') (freq(j),j=iinc1,iinc2)
            write (iw, '('' MASS(I)'',6F10.5,/)') (redmas(j, 1),j=iinc1,iinc2)
            write (iw, '('' DIPX(I)'',6F10.5)') (trdip(1,j),j=iinc1,iinc2)
            write (iw, '('' DIPY(I)'',6F10.5)') (trdip(2,j),j=iinc1,iinc2)
            write (iw, '('' DIPZ(I)'',6F10.5,/)') (trdip(3,j),j=iinc1,iinc2)
            write (iw, '('' DIPT(I)'',6F10.5)') (dipt(j),j=iinc1,iinc2)
          end do
        end if
        if (nrem6 >= 1) then
          write (iw, '(/)')
          iinc1 = iinc1 + 6
          iinc2 = iinc1 + (nrem6 - 1)
          write (iw, '(3X,''I'',10I10)') (j,j=iinc1,iinc2)
          write (iw, '('' FREQ(I)'',6F10.4)') (freq(j),j=iinc1,iinc2)
          write (iw, '(/,'' MASS(I)'',6F10.5)') (redmas(j, 1),j=iinc1,iinc2)
          write (iw, '(/,'' DIPX(I)'',6F10.5)') (trdip(1,j),j=iinc1,iinc2)
          write (iw, '('' DIPY(I)'',6F10.5)') (trdip(2,j),j=iinc1,iinc2)
          write (iw, '('' DIPZ(I)'',6F10.5)') (trdip(3,j),j=iinc1,iinc2)
          write (iw, '(/,'' DIPT(I)'',6F10.5)') (dipt(j),j=iinc1,iinc2)
        end if
      end if
      if (ts) then
        write(iw,'(//25x,a,/)')"Atoms used in the FORCETS calculation"
        do i = 1, nvar, 3
          j = loc(1,i)
          if (maxtxt == 26) then
            write(iw,'(i4,4x,a2,"(",a26,")",3(F13.8," +1"))')(i + 2)/3, elemnt(nat(j)),txtatm(j), geo(:,j)
          else
            write(iw,'(10x,i4,4x,a2,3(F13.8," +1"))')(i + 2)/3, elemnt(nat(j)), geo(:,j)
          end if
        end do
      end if
      if (prnt .and. prt_normal_coords) then
        write (iw, '(2/10X,'' NORMAL COORDINATE ANALYSIS (Total motion = 1 Angstrom)'')')
        j = nvar
        i = -nvar
        call matou1 (cnorml, freq, nvib, j, nvib, 5)
      end if
!
!   CARRY OUT IRC IF REQUESTED.
!
      if (index(keywrd,'IRC') + index(keywrd,'DRC') /= 0) then
        if (index(keywrd, " HTML") /= 0) then
           if (index(keywrd,' DIPOLE') /= 0) call write_path_html(2)
           call write_path_html(1)
        end if
        loc(1,:nvar) = 0
        loc(2,:nvar) = 0
        nvar = nvaold
        loc(1,:nvar) = locold(1,:nvar)
        loc(2,:nvar) = locold(2,:nvar)
        call xyzint (coord, numat, na, nb, nc, 1.D0, geo)
        last = 1
        store_coord = coord
        geo(:,:numat) = store_coord(:,:numat)
        na = 0
        jloop = 0
        this_point = 0
        if (ts) then
          allocate(velocity(3*numat))
          if (index(keywrd, " FORCETS") /= 0) then
            k = 1
          else
            k = 0
          end if
!
!  Generate a complete normal mode
!
          velocity = 0.d0
          j = 3
          do i = 1, numat
            if (lopt(1,i) == k) then
              velocity(i*3 - 2) = cnorml(j - 2)
              velocity(i*3 - 1) = cnorml(j - 1)
              velocity(i*3 - 0) = cnorml(j - 0)
              j = j + 3
            end if
          end do
          na_store = 0
          l_normal_html = .false.
          call drc (velocity, freq)
        else
          if (.not. allocated(velocity)) &
            allocate(velocity(3*numat))
          velocity(:3*numat) = cnorml(:3*numat)
          l_normal_html = .false.
          call drc (cnorml, freq)
          cnorml(:3*numat) = velocity(:3*numat)
        end if
        if (moperr) return
        i = Index (keywrd, " IRC")
        if (i /= 0) then
          if  (Index (keywrd(i+1 + Index(keywrd(i+1:), " ") - 2:), "*") /= 0) then
!
! a double sided IRC
!
!   First, reverse the reaction path already written.
!
            If (index(keywrd,' DIPOLE') /= 0) call reverse_trajectory(2)
            call reverse_trajectory(1)
            call reverse_aux
            call l_control("REVERSE", len("REVERSE"), 1)
            last = 1
            geo(:,:numat) = store_coord(:,:numat)
            na = 0
            jloop = 0
            numcal = numcal + 1
            loc(1,:nvar) = 0
            loc(2,:nvar) = 0
            nvar = nvaold
            loc(1,:nvar) = locold(1,:nvar)
            loc(2,:nvar) = locold(2,:nvar)
            if (ts) then
              velocity = 0.d0
              j = 3
              do i = 1, numat
                if (lopt(1,i) == k) then
                  velocity(i*3 - 2) = -cnorml(j - 2)
                  velocity(i*3 - 1) = -cnorml(j - 1)
                  velocity(i*3 - 0) = -cnorml(j - 0)
                  j = j + 3
                end if
              end do
               call drc (velocity, freq)
            else
              cnorml = -cnorml
              call drc (cnorml, freq)
            end if
          end if
          ndep = 0
          lopt(:,:numat) = 1
        else
          ndep = ndeold
          nvar = 0
          geo(:,:natoms) = geoa(:,:natoms)
        end if
        if (index(keywrd, " PDBOUT") /= 0) close(ixyz, status = 'delete', iostat=i)
        goto 99
      end if
      call freqcy (fmatrx, freq, deldip, force_const, .FALSE., deldip, ff, oldf, ts)
      if (prt_normal_coords) then
      write (iw, '(2/10X,'' MASS-WEIGHTED COORDINATE ANALYSIS (NORMAL COORDINATES)'')')
        i = -nvar
        j = nvar
        call matou1 (cnorml, freq, nvib, j, nvib, 5)
      end if
      if (.not. ts) call anavib (freq, dipt, nvar, cnorml, store, nvib, fmatrx, ff)
  98 continue
      if (index(keywrd,'THERMO') /= 0) then
        call gmetry (geo, coord)
        sym = irot(igroup)
        write (iw, '(/5x,A,A,A,I3)') ' SYMMETRY NUMBER FOR POINT-GROUP ', &
        & name, '=', irot(igroup)
        i = index(keywrd,' TRANS')
!
!   "I" IS GOING TO MARK THE BEGINNING OF THE GENUINE VIBRATIONS.
!
        if (i /= 0) then
          i = index(keywrd,' TRANS=')
          if (i /= 0) then
            i = nint(1 + reada(keywrd,i))
            j = nvib - i + 1
            if (i - 1 > 0) &
              write (iw, '(2/1X,''THE LOWEST'',I3,'' VIBRATIONS ARE NOT'',/, &
      & '' TO BE USED IN THE THERMO CALCULATION'')') i - 1
          else
            write (iw, '(2/10X,''SYSTEM IS A TRANSITION STATE'')')
            i = 2
            j = nvib - 1
          end if
        else
          if (numat > 1) write (iw, '(2/10X,''SYSTEM IS A GROUND STATE'')')
          i = 1
          j = nvib
        end if
        if (.not. allocated(freq)) then
          allocate (freq(1))
          freq = 0.d0
        end if
        if (.not. allocated(store)) allocate(store(1))
        store(:nvib) = freq(:nvib)
        call thermo (a, b, c, linear, sym, store(i), j, escf)
      end if
      if (allocated(store_coord)) coord(:,:numat) = store_coord(:,:numat)
      call to_screen("To_file: Force output")
      if (numat > 1) then
        fmatrx = oldf*1.d-5
        if ( .not. ts .and. numat > 1) call intfc (oldf, xparam, geoa, nar, nbr, ncr)
      end if
      na(1) = 0
      nvar = 0
      ndep = ndeold

  99  if (allocated(dipt))    deallocate (dipt)
      if (allocated(travel))  deallocate (travel)
      if (allocated(freq))    deallocate (freq)
      if (allocated(redmas))  deallocate (redmas)
      if (store_natoms /= natoms + id) then
        natoms = -30; numat = -30
      end if
      return
      end subroutine force
      subroutine phase_lock(vecs,n)
      implicit none
   !
      integer :: n
      double precision, dimension (n,n) :: vecs
!  Local
      integer :: i, j
      double precision :: asum, ssum
!
!  The convention used in MOPAC is that the phase of an eigenvector is
!  defined such that the largest coefficient is positive.
!
      do i = 1, n
        asum = 0.d0
        ssum = 0.d0
        do j = 1, n
          if (Abs(vecs(j,i)) > asum) then
            asum = Abs(vecs(j,i))
            ssum = vecs(j,i)
          end if
        end do
        if (ssum < 0.d0) then
!
!  The largest coefficient was negative - therefore reverse the phase of the eigenvector
!
          vecs(1:n,i) = -vecs(1:n,i)
        end if
      end do
      return
      end subroutine phase_lock
