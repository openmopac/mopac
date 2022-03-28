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

      double precision function meci ()
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use meci_C, only : rjkaa, rjkab, microa, microb, nmeci, nelec, nstate, &
      & lab, labsiz, cdiagi, cdiag, occa, maxci, nalmat, conf, nmos, &
      & vectci, xy, dijkl, ispin, eig, ispqr,  k, deltap, root_requested
      use symmetry_C, only : jndex, namo, state_spin, state_Irred_Rep, &
      state_QN
      USE common_arrays_C, only : nfirst, nlast, eigs, c
      USE molkst_C, only : nopen, nclose, fract, nelecs, norbs, numcal, &
      & numat, keywrd, last, lm61, rjkab1, msdel, line
      USE chanel_C, only : ir, iw
      use cosmo_C, only : iseps
      USE mndod_C, only : fx
      implicit none
!
      integer , dimension(:), allocatable :: nfa
      integer :: icalcn, ndoubl, j, l, mdim, lroot, ne, lima, &
        limb, i1, j1, ntot, izero, limci, nupp, &
        ndown, loc, i2, j2, iroot, i, ii, ji, m, maxvec, &
        iuj, iofset, il, iu, irep(maxci), smult, jroot
      integer, dimension (20, 20) :: qn_temp
      integer, dimension(:,:), allocatable :: nperma, npermb
      double precision, dimension(:), allocatable :: eiga, cimat, spin, cij, ckl, wcij, &
        delta
      double precision, dimension(:), allocatable :: diag
      double precision, dimension(nmeci, nmeci) :: deltapp
      double precision, dimension(:,:), allocatable :: oscil
      double precision, dimension(3) :: work
      double precision :: xx, x, gse, y, sum, summ
!      double precision, allocatable :: overlap(:)
      double precision, external :: diagi, reada
      logical :: debug, large, prnt, lspin, lspin1, peci, first1, bigprt, sing, &
        doub, trip, quar, quin, sext, prnt2, geook, getmic, sept, octe, none, &
        cis, cisd, cisdt
      character, dimension(11) :: tspin*8
      character :: srep(maxci)*15, root_ir*4, tex_root*1, num1*1, num2*1

      save spin, debug, large, lspin1, first1, sing, doub, trip, quar, quin, &
        sext, prnt2, geook, getmic, sept, tspin, icalcn, ndoubl, j, l, delta, &
        mdim, lroot, smult, ne, xx, lima, limb, eiga, octe, none, peci, cis, cisd, cisdt, &
        oscil, nperma, npermb, nfa, limci, nupp, ndown, root_ir
!-----------------------------------------------
!**********************************************************************
!
!                 PROGRAM MECI
!
!   A MULTI-ELECTRON CONFIGURATION INTERACTION CALCULATION
!
!   WRITTEN BY JAMES J. P. STEWART, AT THE
!              FRANK J. SEILER RESEARCH LABORATORY
!              USAFA, COLORADO SPRINGS, CO 80840
!
!              1985
!
!**********************************************************************
!
!
!   MATRICES FOR PERMUTATION WORK
!
!
!   MATRICES FOR ONE AND TWO ELECTRON INTEGRALS
!
!
!   SPIN MATRICES
!
!
!   MATRICES FOR SEC.DET., VECTORS, AND EIGENVALUES.
!
      data icalcn/ 0/
      data tspin/ 'SINGLET ', 'DOUBLET ', 'TRIPLET ', 'QUARTET ', 'QUINTET ', &
      'SEXTET  ', 'SEPTET  ', 'OCTET   ', 'NONET   ', 'DECET   ', '??????? '/
      meci = 0.d0
      if (numat == 0) then
        if (allocated(spin))   deallocate(spin)
        if (allocated(nalmat)) deallocate(nalmat)
        if (allocated(eig))    deallocate(eig)
        if (allocated(conf))   deallocate(conf)
        if (allocated(vectci)) deallocate(vectci)
        if (allocated(ispin))  deallocate(ispin)
        if (allocated(oscil))  deallocate(oscil)
        if (allocated(ispqr))  deallocate(ispqr)
        if (allocated(occa))   deallocate(occa)
        if (allocated(deltap)) deallocate(deltap)
        if (allocated(nperma)) deallocate(nperma)
        if (allocated(npermb)) deallocate(npermb)
        if (allocated(nfa))    deallocate(nfa)
        if (allocated(delta))  deallocate(delta)
        if (allocated(eiga))   deallocate(eiga, rjkaa, xy, rjkab, dijkl)
        meci = 0.d0
        return
      end if
      if (icalcn /= numcal) then
        icalcn = numcal
        first1 = .TRUE.
        limci = 0
        mdim = maxci
        geook = index(keywrd,'GEO-OK') /= 0
        lspin1 = index(keywrd,'ESR') /= 0
        debug = index(keywrd,'DEBUG') /= 0
        prnt2 = index(keywrd,'MECI') /= 0
        debug = debug .and. prnt2
        large = index(keywrd,'LARGE') /= 0
        peci  = index(keywrd,' PECI') /= 0
        cis   = index(keywrd,' CIS ') /= 0
        cisd  = index(keywrd,' CISD ') /= 0
        cisdt = index(keywrd,' CISDT ') /= 0
        root_ir = "XXXX"
        ndoubl = 99
        i = index(keywrd,'C.I.=(')
        if (i /= 0) then
          j = index(keywrd(i:i+10),',') + i - 1
          ndoubl = nint(reada(keywrd,j))
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=(') + 5))
          ndoubl = min(ndoubl,nclose)
        else if (index(keywrd,'C.I.=') /= 0) then
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=') + 5))
        else
          nmos = nopen - nclose
        end if
        nmos = min(nmos,norbs)
        if (allocated(occa))    deallocate(occa)
        if (allocated(deltap))  deallocate(deltap)
        allocate(occa(nmos), deltap(nmos, nmos))
        lroot = 1
        if (index(keywrd,'EXCI') /= 0) lroot = 2
        i = index(keywrd," ROOT")
        if (i /= 0) then
          do j = i + 6, i + 9
            if (keywrd(j:j) > "9" .or. keywrd(j:j) < "0") exit
          end do
          lroot = Nint (reada (keywrd(:j), i + 6))
          do j = i + 7, i + 14
            if (keywrd(j:j) == " ") exit
            if (keywrd(j:j) >= "A" .and. keywrd(j:j) <= "Z") then
              root_ir = keywrd(j:)
              i = Index(root_ir," ")
              if (i /= 0) root_ir(i:) = " "
              do k = 2,4
                if (root_ir(k:k) >="A" .and. root_ir <= "Z") &
                & root_ir(k:k) = Char(Ichar(root_ir(k:k)) + Ichar('a') - Ichar('A'))
              end do
            end if
            if (root_ir /= " ") exit
          end do
        end if
        if (ndoubl == 99) then
          j = max(min((nclose + nopen + 1)/2 - (nmos - 1)/2,norbs - nmos + 1),1)
        else
          j = nclose - ndoubl + 1
          if (fract > 1.99D0) j = j + 1
        end if
        l = 0
        if (nclose - j + 1 > 0) then
          occa(:nclose-j+1) = 1
          l = nclose - j + 1
        end if
        j = max(j,nclose + 1)
        if (nopen - j + 1 > 0) then
          occa(l+1:nopen-j+1+l) = fract*0.5D0
          l = nopen - j + 1 + l
        end if
        do i = 1, nmeci
          l = l + 1
          if (l > nmos) exit
          occa(l) = 0.D0
        end do
        sing = index(keywrd,' SING') + index(keywrd,' EXCI') + index(keywrd,&
          ' BIRAD') /= 0
        doub = index(keywrd,' DOUB') /= 0
        trip = index(keywrd,' TRIP') /= 0
        quar = index(keywrd,' QUAR') /= 0
        quin = index(keywrd,' QUIN') /= 0
        sext = index(keywrd,' SEXT') /= 0
        sept = index(keywrd,' SEPT') /= 0
        octe = index(keywrd,' OCTE') /= 0
        none = index(keywrd,' NONE') /= 0
        smult = -1
        if (sing) smult = 1
        if (doub) smult = 2
        if (trip) smult = 3
        if (quar) smult = 4
        if (quin) smult = 5
        if (sext) smult = 6
        if (sept) smult = 7
        x = 0.D0
        do j = 1, nmos
          x = x + occa(j)
        end do
        xx = x + x
        ne = nint(xx)
        i = (nelecs - ne + 1)/2
        nelec = (nelecs - ne + 1)/2
      end if
      root_requested = 0
      prnt = debug .or. last==3 .and. prnt2
      bigprt = prnt .and. large
!
!    TEST TO SEE IF THE SET OF ENERGY LEVELS USED IN MECI IS COMPLETE,
!    I.E., ALL COMPONENTS OF DEGENERATE IRREDUCIBLE REPRESENTATIONS
!    ARE USED.  IF NOT, THEN RESULTS WILL BE NONSENSE.  GIVE USERS A
!    CHANCE TO REALLY FOUL THINGS UP BY ALLOWING JOB TO CONTINUE IF
!    'GEO-OK' IS SPECIFIED.
!
      if (nelec < 0) then
        meci = 0.D0
        call mopend (&
           ' NUMBER OF ELECTRONS IN M.O.S BELOW ACTIVE SPACE IS LESS THAN ZERO')
        return
      end if
      if (nelec + nmos > norbs) then
        meci = 0.D0
        call mopend (&
       ' UPPER BOUND OF ACTIVE SPACE IS GREATER THAN THE NUMBER OF ORBITALS!')
        return
      end if
      if (allocated(eiga)) then
        deallocate(eiga, rjkaa, xy, rjkab, dijkl, delta)
      end if
      allocate(eiga(nmos), rjkaa(nmos,nmos), xy(nmos, nmos, nmos, nmos), &
      & rjkab(nmos, nmos), dijkl(norbs, nmos,(nmos*(nmos + 1))/2))
      allocate(cij(lm61), ckl(lm61), wcij(lm61), delta(nmeci*norbs))
      eiga(:nmos) = eigs(1+nelec:nmos+nelec)
      lspin = lspin1 .and. last==3
      if (bigprt) then
        write (iw, '(''  INITIAL EIGENVALUES'')')
        write (iw, '(5F12.6)') (eiga(i),i=1,nmos)
        write (iw, '(2/10X,''NUMBER OF ELECTRONS IN C.I. ='',F5.1)') xx
      end if
      i = nmos + nelec + 1
      if (.not.geook .and. nelec>0) then
        if (abs(eigs(nelec+1)-eigs(nelec))<2.D-2 .or. nelec+1+nmos .lt. norbs &
        .and. abs(eigs(nelec+1+nmos)-eigs(nelec+nmos))<2.D-2) then
          write (iw, '(3/10X,A)') 'DEGENERATE ENERGY LEVELS DETECTED IN MECI'
          write (iw, '(10X,A)') &
            'SOME OF THESE LEVELS WOULD BE TREATED BY MECI,'
          write (iw, '(10X,A)') 'WHILE OTHERS WOULD NOT.  THIS WOULD RESULT IN'
          write (iw, '(10X,A)') 'NON-REPRODUCIBLE ELECTRONIC ENERGIES.'
          write (iw, '(/2X,A)') 'EIGENVALUES OF ACTIVE SPACE'
          write (iw, '(8F10.4)') (eigs(i),i=nelec + 1,nelec + nmos)
          write (iw, '(/2X,A)') 'EIGENVALUES BELOW ACTIVE SPACE'
          write (iw, '(8F10.4)') (eigs(i),i=max(1,nelec - 5),nelec)
          write (iw, '(/2X,A)') 'EIGENVALUES ABOVE ACTIVE SPACE'
          write (iw, '(8F10.4)') (eigs(i),i=nelec + nmos + 1,min(norbs,nelec + &
            nmos + 5))
          write (iw, *)
          meci = 0.D0
          call mopend ('JOB STOPPED. TO CONTINUE, SPECIFY "GEO-OK".')
          return
        end if
      end if
!      allocate(overlap((norbs*(norbs + 1))/2))
!      call fill_overlap_matrix(overlap)
      if (bigprt) then
        write (iw, '(2/10X,''EIGENVECTORS'',/)')
        do i = 1, norbs
          write (iw, '(6F12.6)') (c(i,j+nelec),j=1,nmos)
        end do
      end if
      call ijkl (c(1,nelec+1), c, nelec, nmos, dijkl, cij, ckl, wcij, xy)
      do i = 1, nmos
        do j = 1, nmos
          rjkaa(i,j) = xy(i,i,j,j) - xy(i,j,i,j)
          rjkab(i,j) = xy(i,i,j,j)
        end do
      end do
      rjkab1 = rjkab(1,1)
!
!  Divide M.O. eigenvalue by the overlap integral of the M.O.
!
!     if (allocated(overlap)) deallocate(overlap)
!      write (iw, '(/,5x,a,/)') 'EIGENVALUES BEFORE DIVISION BY M.O. OVERLAP'
!      write (iw, '(8F10.5)') (eiga(i),i=1,nmos)
!      allocate(overlap((norbs*(norbs + 1))/2))
!      call fill_overlap_matrix(overlap)
!      write(iw,*)" M.O.  overlap"
!      do j = 1, nmos
!        sum = 0.d0
!        l = 0
!        do i = 1,norbs
!          do k = 1, i - 1
!            l = l + 1
!            sum = sum + c(i, j + nelec)*overlap(l)*c(k, j + nelec)
!          end do
!            l = l + 1
!            sum = sum + c(i, j + nelec)*0.5d0*c(k, j + nelec)
!        end do
!        sum = sum*2.d0
!        write(iw,'(i4, f10.4)')j + nelec, sum
!        eiga(j) = eiga(j)/sum
!      end do
!      write (iw, '(/,5x,a,/)') 'EIGENVALUES AFTER DIVISION BY M.O. OVERLAP'
!      write (iw, '(8F10.5)') (eiga(i),i=1,nmos)
      do i = 1, nmos
        x = 0.D0
        do j = 1, nmos
          x = x + (rjkaa(i,j)+rjkab(i,j))*occa(j)
        end do
        eiga(i) = eiga(i) - x
      end do
      if (bigprt) then
        write (iw, 120)
  120   format(/,5x,'EIGENVALUES AFTER REMOVAL OF INTER-ELECTRONIC',&
          ' INTERACTIONS',/)
        write (iw, '(8F10.5)') (eiga(i),i=1,nmos)
        write (iw, '(3/10X,''TWO-ELECTRON J-INTEGRALS'',/)')
        do i1 = 1, nmos
          write (iw, '(10F8.4)') (rjkab(i1,j1),j1=1,min(10,nmos))
        end do
        write (iw, '(3/10X,''TWO-ELECTRON K-INTEGRALS'',/)')
        do i1 = 1, nmos
          write (iw, '(10F8.4)') (rjkab(i1,j1) - rjkaa(i1,j1),j1=1,min(nmos,10)&
            )
        end do
      end if
      i = 1
      if (nmos > 0) then
        j = 1
        rjkaa(:nmos,:nmos) = rjkaa(:nmos,:nmos)*0.5D0
        j = nmos + 1
        i = nmos + 1
      end if
      if (first1) then
        getmic = index(keywrd,'MICROS') /= 0
        if (getmic) then
          k = nint(reada(keywrd,index(keywrd,'MICROS')))
          lab = k
          if (prnt) write (iw, '(''    MICROSTATES READ IN'')')
          ntot = nint(xx)
          rewind ir
          do i = 1, 10000
            read (ir, '(A)') line
            call upcase (line, 80)
            if (index(line,'MICRO') == 0) cycle
            exit
          end do
          do i = 1, 1000
            read (ir, '(A)', end=190, err=190) line
            call upcase (line, 80)
            if (index(line,'MICRO') /= 0) go to 200
          end do
  190     continue
          meci = 0.D0
          call mopend (&
             'MICROSTATES SPECIFIED BY KEYWORDS BUT MISSING FROM DATA')
          return
  200     continue
          if (allocated(microa)) deallocate(microa)
          if (allocated(microb)) deallocate(microb)
          allocate(microa(nmos, lab), microb(nmos,lab))
          do i = 1, lab
            read (ir, '(A)') line
            izero = max(0,min(index(line,'0'),index(line,'1')) - 1)
            do j = 1, nmos
              if (line(j+izero:j+izero) /= '1') line(j+izero:j+izero) = '0'
              if (line(j + nmos + izero:j + nmos + izero) /= '1') &
              line(j + nmos + izero:j + nmos + izero) = '0'
              microa(j,i) = ichar(line(j+izero:j+izero)) - ichar('0')
              microb(j,i) = ichar(line(j+nmos+izero:j+nmos+izero)) - ichar('0')
            end do
            if (prnt) write (iw, '(20I3)') (microa(j,i),j=1,nmos), (microb(j,i),j=1,nmos)
            k = 0
            do j = 1, nmos
              k = k + microa(j,i) + microb(j,i)
            end do
            if (k == ntot) cycle
            ntot = k
            xx = k
            write (iw, &
      '(/10X,''NUMBER OF ELECTRONS IN C.I. REDEFINED TO:'',I4,/)') k
          end do
          if (getmic) go to 240
          limci = 0
        end if
        if (msdel == 0 .and. mod(ne,2) == 1) msdel = 1
        nupp = (ne + msdel)/2
        ndown = ne - nupp
        if (nupp*ndown < 0) then
          write (iw, '(/10X,''IMPOSSIBLE VALUE OF DELTA S'')')
          meci = 0.D0
          call mopend ('IMPOSSIBLE VALUE OF DELTA S')
          return
        end if
        if (limci == 0) then
          lima = nint(fx(nmos + 1)/(fx(nupp + 1)*fx(nmos - nupp + 1)))
          limb = nint(fx(nmos + 1)/(fx(ndown + 1)*fx(nmos - ndown + 1)))
          if (allocated(nperma)) deallocate(nperma)
          if (allocated(npermb)) deallocate(npermb)
          allocate(nperma(nmos, lima), npermb(nmos, limb), stat = i)
          if (i /= 0) then
            meci = 0.d0
            call mopend("A problem occurred during memory assignment. The number of configurations is too large. ")
            return
          end if
        end if
        lab = lima*limb

        if (peci) limci = 4
        if (cis .or. cisd .or. cisdt) then
          if (cis)   limci = 2
          if (cisd)  limci = 4
          if (cisdt) limci = 6
        end if
        call perm (nperma, nupp, nmos, lima, limci)
        call perm (npermb, ndown, nmos, limb, limci)
        lab = lima*limb
        if (allocated(microa)) deallocate(microa)
        if (allocated(microb)) deallocate(microb)
        allocate(microa(nmos, lab), microb(nmos,lab), stat = i)
        if (i /= 0) then
          meci = 0.d0
          call mopend("A problem occurred during memory assignment. The number of configurations is too large. ")
          return
        end if
      else
        if (.not.getmic) lab = lima*limb
      end if
  240 continue
      if (allocated(spin))   deallocate(spin)
      if (allocated(nalmat)) deallocate(nalmat)
      if (allocated(eig))    deallocate(eig)
      if (allocated(conf))   deallocate(conf)
      if (allocated(vectci)) deallocate(vectci)
      if (allocated(ispin))  deallocate(ispin)
      if (allocated(oscil))  deallocate(oscil)
      if (allocated(ispqr))  deallocate(ispqr)
      allocate(spin(lab), nalmat(lab), eig(lab + 1), vectci(30*lab), &
      & ispin(lab), oscil(3,lab + 4), ispqr(lab,nmeci + 1), stat = i)
      if (i /= 0) then
          meci = 0.d0
          call mopend("A problem occurred during memory assignment. The number of configurations is too large. ")
          return
        end if
      ispqr = -99999
      if (prnt) write (iw, 250) (nupp - ndown)*0.5D0
  250 format(10x,'COMPONENT OF SPIN  = ',f4.1)
      gse = 0.0D0
      do i = 1, nmos
        gse = gse + eiga(i)*occa(i)*2.D0
        gse = gse + xy(i,i,i,i)*occa(i)*occa(i)
        do j = i + 1, nmos
          gse = gse + 2.D0*(2.D0*xy(i,i,j,j)-xy(i,j,i,j))*occa(i)*occa(j)
        end do
      end do
      if (limci /= 0) then
!
! REMOVE ALL UNWANTED EXCITATIONS
!
        loc = 0
        do i = 1, lima
          do j = 1, limb
!
!        DETERMINE THE TYPE OF EXCITATION
!
            k = 0
            do l = 1, nmos
              k = k + abs(nperma(l,i)-nperma(l,1)) + abs(npermb(l,j)-npermb(l,1))
            end do
            if (peci) then
              l = 1
              if (nmos > 0) then
                do l = 1, nmos
                  if (.not.k<=2 .and. .not.nperma(l,i)==npermb(l,j)) k = 1000
                end do
                l = nmos + 1
              end if
            end if
!
!        COPY DET AND INCREMENT LOC
!
            if (k > limci) cycle
            loc = loc + 1
            if (loc <= 0) cycle
            l = 1
            if (nmos > 0) then
              microa(:nmos,loc) = nperma(:nmos,i)
              microb(:nmos,loc) = npermb(:nmos,j)
              l = nmos + 1
            end if
          end do
        end do
        lab = loc
      end if
      num1 = char(ichar("2") + int(log10(lab*1.0001)))
      if (prnt) write (iw, "(/,/,10x,' NO OF CONFIGURATIONS CONSIDERED =',i"//num1//")") lab
      if (allocated(cdiag))  deallocate(cdiag)
      if (lab > maxci) then
        meci = 0.d0
        call mopend("Too many configurations requested")
        write(iw,"(a,i7)")"  Number requested:", lab
        return
      end if
      allocate(conf(lab**2), cimat((lab*(lab + 1))/2 + 9), diag(lab), cdiag(lab), stat = i)
      if (i /= 0) then
        meci = 0.d0
        call mopend("A problem occurred during memory assignment. The number of configurations is too large. ")
        return
      end if
      if (getmic .or. limci/=0) then
        do i = 1, lab
          diag(i) = diagi(microa(1,i),microb(1,i),eiga,xy,nmos) - gse
          cdiag(i) = cdiagi
        end do
        go to 400
      end if
      j = 0
      i = 0
      do i1 = 1, lima
        do i2 = 1, limb
          i = i + 1
          cimat(i) = diagi(nperma(1,i1),npermb(1,i2),eiga,xy,nmos) - gse
        end do
      end do
      if (lab <= mdim) then
        k = 0
        do i = 1, lima
          do j = 1, limb
            k = k + 1
            diag(k) = cimat(k)
            l = 1
            if (nmos > 0) then
              microa(:nmos,k) = nperma(:nmos,i)
              microb(:nmos,k) = npermb(:nmos,j)
              l = nmos + 1
            end if
          end do
        end do
      else
!
!   SELECT THE MDIM LOWEST MICROSTATES
!
        do lab = 1, mdim
          x = 1000.D0
          i = 0
          j1 = 0
          j2 = 0
          do i1 = 1, lima
            do i2 = 1, limb
              i = i + 1
              if (cimat(i) >= x) cycle
              x = cimat(i)
              j = i
              j1 = i1
              j2 = i2
            end do
          end do
!
!   MICROSTATE J IS THE LOWEST IN THE REMAINING SET
!
          k = 1
          if (nmos > 0) then
            microa(:nmos,lab) = nperma(:nmos,j1)
            microb(:nmos,lab) = npermb(:nmos,j2)
            k = nmos + 1
          end if
          diag(lab) = cimat(j)
          cimat(j) = 1.D8
        end do
        lab = mdim
      end if
  400 continue
      if (lab > maxci) then
        meci = 0.d0
        call mopend("Too many configurations requested")
        write(iw,"(a,i7)")"  Number requested:", lab
        write(iw,"(a,i7)")"  Max. no. allowed:", maxci
        return
      end if
      do i = 1, lab
        k = 0
        x = 0.D0
        do j = 1, nmos
          x = x + microa(j,i)*microb(j,i)
          k = k + microa(j,i)
        end do
        nalmat(i) = k
        spin(i) = 4.D0*x - (xx - 2*nalmat(i))**2
      end do
!
!   BEFORE STARTING, CHECK THAT THE ROOT WANTED CAN EXIST
!
      if (lab < lroot) then
        call mopend (&
       'C.I. IS OF SIZE LESS THAN ROOT SPECIFIED. MODIFY SIZE OF C.I. OR ROOT NUMBER.')
        write (iw, '(10X,''MODIFY SIZE OF C.I. OR ROOT NUMBER'')')
        meci = 0.D0
        return
      end if
      if (prnt) then
        write (iw, &
      '(/,'' CONFIGURATIONS CONSIDERED IN C.I.      '',/,'' M.O. NUMBER :      '',20I4)') &
       (i, i = nelec + 1, nelec + min(20, nmos))
       if (nmos > 20) write (iw, '(20X,20I4)') (i,i=nelec + 21, nelec + min(40,nmos))
        write (iw, '(''          ENERGY'')')
        num1 = char(ichar("2") + int(log10(lab*1.0001)))
        num2 = char(ichar("6") - int(log10(lab*1.0001)))
        do i = 1, lab
          write (iw, '(/7x,'//num2//'x,I'//num1//',5X,20I4)') i, (microa(k,i),k=1,min(20,nmos))
          if (nmos > 20) write (iw, '(20X,20I4)') (microa(k,i),k=21,min(40,nmos))
          write (iw, '(6X,F10.4,4X,20I4)') diag(i), (microb(k,i),k=1,min(20,nmos))
          if (nmos > 20) write (iw, '(20X,20I4)') (microb(k,i),k=21,min(40,nmos))
        end do
      end if
      call mecih (diag, cimat, nmos, lab, xy)
      if (bigprt) then
        write (iw, '(2/,'' C.I. MATRIX'')')
        call vecprt (cimat, (-lab))
      end if
!
!  Sometimes exact degeneracies cause problems with RSP, so perturb the C.I. matrix
!  to destroy exact degeneracy.
!
      do i = 1, lab
        cimat((i*(i + 1))/2) = cimat((i*(i + 1))/2) + 1.d-10*i*(-1)**i
      end do
      labsiz = min(lab,lroot + 10)
      if (last == 3 .and. prnt2 .or. root_requested /= 1 .and. first1 .or. root_ir /= "XXXX") then
        labsiz = lab
      end if
      call rsp (cimat, lab, eig, conf)
      if (bigprt .or. last > 0 .or. root_ir /= "XXXX") call symtrz (c(1,nelec+1), eig, 3, .TRUE.)
      if (bigprt) then
        write (iw, '(2/20X,''STATE VECTORS'',2/)')
        if (last == 3) then
          call matou1 (conf, eig, lab, lab, lab, 4)
        else
          call matout (conf, eig, lab, lab, lab)
        end if
      end if
      if (prnt) then
        write (iw, &
      '(/,''  STATE       ENERGY (EV)        Q.N.  SPIN   SYMMETRY              POLARIZATION'')')
        write (iw, &
      '(''         ABSOLUTE     RELATIVE'',28X,''    X           Y           Z'',/)')
      end if
      iroot = 0
      jroot = 1
      cimat = 0.1D0
      ispin(:lab) = 0
      do i = 1, labsiz
        if (last == 0 .and. root_ir == "XXXX") then
          namo(i) = ' '
          jndex(i) = i
        end if
        x = 0.5D0*xx
        ii = (i - 1)*lab
        do j = 1, lab
          ji = j + ii
          x = x - conf(ji)*conf(ji)*spin(j)*0.25D0
          k = ispqr(j,1)
          if (k == 1) cycle
          l = 2
          if (k - 1 > 0) then
            do l = 1, k - 1
              x = x + conf(ji)*conf(ispqr(j,l+1)+ii)*2.D0
            end do
            l = k + 1
          end if
        end do
        y = (-1.D0 + sqrt(1.D0 + 4.D0*x))*0.5D0
        ispin(i) = nint(y*2.D0 + 1)
        if (Abs(ispin(i) - y*2.d0 - 1.d0) > 0.2d0) then
          ispin(i) = 11
        else
          ispin(i) = Min(ispin(i), 11)
        end if
        cimat(j) = cimat(j) + 1
      end do
!
!   Reset jndex to allow for spin.
!
    qn_temp = 0
    srep = " "
    irep = 0
    if (allocated(nfa))     deallocate(nfa)
    i = Max(maxci, nmos + 3, lab)
    allocate(nfa(i))
    nfa(1:lab) = jndex(1:lab)
    j = 0
    l = 0
    do i = 1, lab
      if (i > 1) then
!
!  If part of a degenerate manifold, use previous quantum number
!
        if (nfa(i) == nfa(i-1) .and. namo(i) == namo(i-1)) then
          jndex(i) = jndex(i-1)
          if (smult < 0 .or. ispin(i) == smult) iroot = iroot + 1
          cycle
        end if
      end if
!
!  Find Irreducible Representation
!
      do k = 1, j
        if (srep(k) == namo(i)) exit
      end do
      if (k > j) then
        j = j + 1
        srep(j) = namo(i)
      end if
      do m = 1, l
        if (irep(m) == ispin(i)) exit
      end do
      if (m > l) then
        l = l + 1
        irep(l) = ispin(i)
      end if
      if (k>20 .or. m > 20) then
        write(iw,*)" An error has ocurred in MECI, here are some diagnostics"
        write(iw,"(a,i5)")"Number of states (lab):              ", lab
        write(iw,"(a,i5)")"Number of states calculated (labsiz):", labsiz
        write(iw,"(a,i5)")"Error ocurred on state number:       ", i
        write(iw,"(a,2i5)")"Value of 'k' and 'm':", k,m
        if (m > 20) then
          write(iw,*)" Values of ispin"
          write(iw,"(10i5)")ispin(:lab)
          write(iw,*)" Values of irep"
          write(iw,"(10i5)")irep(:m)
        end if
        stop 'Dimension error in MECI.F90'
      end if
      qn_temp(k, m) = qn_temp(k, m) + 1
      jndex(i) = qn_temp(k, m)
      if (root_requested < 1 .and. (smult < 0 .or. ispin(i) == smult)) then
        iroot = iroot + 1
        if (iroot >= lroot .and. jroot <= lroot .and. root_ir == "XXXX" .or. &
          &   lroot == jndex(i) .and. root_ir == namo(i))  then
!
!    Select spin root "smult" if spin is defined, otherwise use all spins
!    Select state with lroot if no I.R. is given, otherwise select the "lroot"'th
!    state of names I.R.
!
!   The root requested has been found
!
          root_requested = i
          meci = eig(i)
          m = 0
          do k = i, min(lab,i + 22)
            if (abs(eig(k)-eig(i)) > 1.D-2) go to 510
            if (lab > 0) then
              vectci(m+1:lab+m) = conf(1+lab*(k-1):lab*k)
              m = lab + m
            end if
          end do
          k = min(lab,i + 4) + 1
  510     continue
          nstate = k - i
        end if
        jroot = iroot
      end if
    end do
      if (root_requested < 1) then
        call mopend ("ROOT REQUESTED DOES NOT EXIST IN C.I.")
        write (iw,'(a,i5,a)')"  Root requested:",lroot,root_ir
        if (sing) write(iw,'(a)')"  State requested must be a Singlet."
        if (doub) write(iw,'(a)')"  State requested must be a Doublet."
        if (trip) write(iw,'(a)')"  State requested must be a Triplet."
        if (quar) write(iw,'(a)')"  State requested must be a Quartet."
        if (quin) write(iw,'(a)')"  State requested must be a Quintet."
        if (sext) write(iw,'(a)')"  State requested must be a Sextet."
        if (octe) write(iw,'(a)')"  State requested must be a Octet."
        if (none) write(iw,'(a)')"  State requested must be a Nonet."
        if (smult < -0.1d0) write(iw,'(a)')"  State requested can have any spin."
        write (iw,'(/,a,/)')"      Roots available:"
        do i = 1, labsiz
          k = min(labsiz, i + 14)
          do j = i + 1, k
            if (.not.(eig(j)-eig(j-1)>1.D-4 .or. jndex(j)/=jndex(j-1) &
             .or. namo(j)/=namo(j-1))) cycle
            exit
          end do
          if (i > 1) then
            if (eig(i) - eig(i-1)>1.D-4 .or. jndex(i)/=jndex(i-1) .or. &
            namo(i) /= namo(i-1)) &
            write (iw, "(i6,2x,a8,1x,a4)")jndex(i), tspin(ispin(i)), namo(i)
          else
            write (iw, "(i6,2x,a8,1x,a4)") jndex(i), tspin(ispin(i)), namo(i)
          end if
        end do
        meci = 0.d0
        return
      end if
      state_Irred_Rep = namo(root_requested)
      state_spin = tspin(ispin(root_requested))
      state_QN = jndex(root_requested)

      if (prnt) then
        call ciosci (c(1,nelec+1), root_requested, oscil, conf)
        deltapp = 0.d0
        delta = 0.d0
        if (iseps) call dmecip (c, deltapp, delta, eig, vectci, conf)
        do i = 1, labsiz
          tex_root = " "
          k = min(labsiz, i + 14)
          do j = i + 1, k
            if (.not.(eig(j)-eig(j-1)>1.D-4 .or. jndex(j)/=jndex(j-1) &
             .or. namo(j)/=namo(j-1))) cycle
            exit
          end do
          work = 0.D0
          do k = i, j - 1
            work = work + oscil(:,k)**2
            l = 4
          end do
          do k = i, j - 1
            if(root_requested == k) tex_root = "+"
          end do
          if (work(1) + work(2) + work(3) > 1.D-9) then
            k = 3
          else
            k = 0
          end if
          if (i > 1) then
            if (eig(i) - eig(i-1)>1.D-4 .or. jndex(i)/=jndex(i-1) .or. namo(i)&
              /=namo(i-1)) write (iw, 580) i, tex_root, eig(i), eig(i) - eig(1), &
              jndex(i), tex_root, tspin(ispin(i)), namo(i), (sqrt(work(j)),j=1,k)
          else
            write (iw, 580) i, tex_root, eig(i), eig(i) - eig(1), jndex(i), &
            tex_root, tspin(ispin(i)), namo(i), (sqrt(work(j)),j=1,k)
          end if
        end do
      end if
      if (prnt .and. msdel /= 0) write (iw, '(/10x, A, /10x, a)') &
        ' ''RELATIVE ENERGY'' Is relative to the lowest state calculated.', &
        ' This may or may not be the ground state.'
      if (prnt) write (iw, '(10x, A)')  ' The "+" symbol indicates the root used.'
  580 format(i5,a,2f12.6,i6,a1,2x,a8,4x,a4,3f12.4)
      if (root_requested == 0) then
        call mopend (&
       'THE STATE REQUIRED IS NOT PRESENT IN THE SET OF CONFIGURATIONS AVAILABLE')
        write (iw, &
      '(/ 4X,''NUMBER OF STATES ACCESSIBLE USING CURRENT KEY-WORDS'',/)')
        do i = 1, 7
          if (cimat(i) <= 0.5D0) cycle
          write (iw, '((24X,A8,I4))') tspin(i), nint(cimat(i))
        end do
        meci = 0.D0
        return
      end if
      maxvec = 0
      if (lspin) maxvec = root_requested + min(3,lab - root_requested)
      if (lspin .and. (ne/2)*2==ne) write (iw, &
        '(''   ESR SPECIFIED FOR AN EVEN-ELECTRON SYSTEM'')')
      do iuj = root_requested, maxvec
        iofset = (iuj - 1)*lab
        write (iw, &
      '(2/,''      MICROSTATE CONTRIBUTIONS TO STATE EIGENFUNCTION'',I3)') iuj
        write (iw, '(5F13.6)') (conf(i+iofset),i=1,lab)
        do i = 1, lab
          conf(i) = conf(i+iofset)**2
        end do
!                                             SECOND VECTOR!
        do i = 1, nmos
          sum = 0.D0
          do j = 1, lab
            sum = sum + (microa(i,j)-microb(i,j))*conf(j)
          end do
          eiga(i) = sum
        end do
        write (iw, &
          '(/,''    SPIN DENSITIES FROM EACH M.O., ENERGY:''    ,F7.3)') eig(&
          iuj)
        write (iw, '(5F12.6)') (eiga(i),i=1,nmos)
        write (iw, *)
        j = 0
        do i = 1, numat
          if (nlast(i) - nfirst(i) > 5) j = 1
        end do
        write (iw, *) '     SPIN DENSITIES FROM EACH ATOMIC ORBITAL'
        if (j == 1) then
          write (iw, '(a,/)') ' Atom   s        p(x)      p(y)      p(z)     d(x2-y2)   '//&
          &'d(xz)     d(z2)     d(yz)     d(xy)     TOTAL'
        else
          write (iw, '(a,/)') ' Atom   s        p(x)      p(y)      p(z)      TOTAL'
        end if
        do i = 1, numat
          il = nfirst(i)
          iu = nlast(i)
          l = 0
          summ = 0.D0
          do k = il, iu
            l = l + 1
            sum = 0.D0
            do j = 1, nmos
              sum = sum + c(k,j+nelec)**2*eiga(j)
            end do
            summ = summ + sum
            eigs(l) = sum
          end do
          if (l == 9) then
            write (iw, '(i4,10f10.6)') i, (eigs(k),k=1,9), summ
          else if (l == 4) then
            write (iw, '(i4,4f10.6,50x,f10.6)') i, (eigs(k),k=1,4), summ
          else if (l == 1) then
            write (iw,'(i4,f10.6,80X,f10.6)') i, eigs(1), summ
          end if
        end do
      end do
      first1 = .FALSE.
      return
      end function meci
