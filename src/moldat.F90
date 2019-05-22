      subroutine moldat(mode) 
!-----------------------------------------------
!
!   If mode == 1 run silently
!   If mode /= 1 generate normal MOPAC output
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
      use common_arrays_C, only : nfirst, nlast, nat, uspd, txtatm, &
      & pdiag, labels, coord, geo, atmass, na, nbonds, ibonds
!
      USE molmec_C, only : nnhco
!
      USE molkst_C, only : natoms, norbs, nalpha, nbeta, nclose, nopen, &
      & nelecs, fract, numat, mpack, keywrd, n2elec, lm61, moperr, line, &
      & uhf, id, msdel, mol_weight, method_PM6, method_PM7, &
      is_PARAM, formula, ispd, mozyme, nvar, rhf, old_chrge, jobnam, &
      N_3_present, Si_O_H_present, nalpha_open, nbeta_open, pdb_label, &
      method_rm1, use_ref_geo
!
      USE parameters_C, only : natorb, uss, upp, udd, tore, &
      dorbs, zd, zs, zp
!
      USE symmetry_C, only : name
!
      use meci_C, only: nmos
!
      USE chanel_C, only : iw, log, ilog
!
      use elemts_C, only : elemnt
!
!***********************************************************************
!
!       MOLDAT works out essential molecular data, such as orbital
!              counters (NFIRST, NLAST), number of electrons,
!              starting populations, etc.
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.4G  18:38:35  03/15/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reada_I 
      use mopend_I 
      use refer_I 
      use gmetry_I 
      use symtrz_I 
      use vecprt_I 
      use to_screen_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mode 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kharge, i, ndorbs, ia, ib, nheavy, nnull, ii, &
        k, k1, j, ielec, ndoubl, ne, nupp, &
        ndown, l, iminr, jminr, icount, &
        ireal, jreal, ni, n1, n4, n9
      real(double) :: elecs,  c(1), yy, w, sum, rmin   ! c(1) is for a dummy call
      logical :: debug, exci, sing, doub, trip, quar, quin, sext, sept, &
        octe, none, birad, halfe, odd
      character :: num1*1
      real(double), dimension (:), allocatable :: rxyz
!-----------------------------------------------
      debug = index(keywrd,'MOLDAT') /= 0 
!
!   SPECIAL MODIFIERS FOR LIMITATIONS IN AVAILABLE PARAMETERS
!
      i = index(keywrd,' CHARGE=') 
      if (i /= 0) then
        kharge = nint(reada(keywrd,i)) 
        old_chrge = kharge
      else
        kharge = 0
      end if
      elecs = -kharge 
      ndorbs = 0 
      if (uss(1) > (-1.D0)) then 
      call mopend (&
           'THE HAMILTONIAN REQUESTED IS NOT AVAILABLE IN THIS PROGRAM') 
        return  
      endif
!
!  Determine the number of atomic orbitals on each element
!
       do i = 1,100
        dorbs(i) = (zd(i) > 1.d-8) 
        if (dorbs(i)) then
          natorb(i) = 9                   ! Element has "d"-orbitals
        else if (zp(i) > 1.d-20) then
          natorb(i) = 4                   ! Element has "p"-orbitals
        else if (zs(i) > 1.d-20) then
          natorb(i) = 1                   ! Element has one "s"-orbital
        else
          natorb(i) = 0                   ! Element has no orbitals
        end if
      end do
      numat = 0 
      ia = 1 
      ib = 0 
      nheavy = 0 
      nnull = 0 
      do ii = 1, natoms 
        if (labels(ii)/=99 .and. labels(ii)/=107) then 
          numat = numat + 1 
          nat(numat) = labels(ii) 
          nfirst(numat) = ia 
          ni = nat(numat) 
          elecs = elecs + tore(ni)
          ib = ia + natorb(ni) - 1 
          if (natorb(ni) == 9) ndorbs = ndorbs + 5 
          nlast(numat) = ib 
          uspd(ia) = uss(ni) 
          if (ia /= ib) then 
            k = ia + 1 
            k1 = ia + 3 
            do j = k, k1 
              uspd(j) = upp(ni) 
            end do 
            if (ib > ia) then 
              nheavy = nheavy + 1 
            else 
              nnull = nnull + 1 
            endif 
            if (k1 /= ib) then 
              k = k1 + 1 
              uspd(k:ib) = udd(ni) 
            endif 
          endif 
        endif 
        ia = ib + 1 
      end do 
      if (numat == 1) then 
        if (index(keywrd,'FORCE') /= 0 .and. index(keywrd,' THERMO') == 0) then 
          call mopend ('A SINGLE ATOM HAS NO VIBRATIONAL MODES') 
          return  
        endif 
      endif 
      if (mode /= 1) call refer 
      call gmetry (geo, coord) 
      if (mode /= 1)then
        call empiri()
        write(iw,"(/,a,/)")formula(:len_trim(formula))  
        if (moperr) return      
      end if
      allocate(rxyz((numat*(numat + 1))/2), stat = i)
      if (i /= 0) then
        if (mode /= 1) write(iw,*)" The interatomic distance array could not be assigned"
        if (mode /= 1) call mopend("The interatomic distance array could not be assigned")
!
!  Do not attempt to check interatomic distances, or anything related to the geometry
!
        return
      end if
      if (pdb_label .and. index(keywrd, " GEO-OK") + index(keywrd, " 0SCF") == 0) then
!
!  Sanity check - does the atom name match its label?
!
        l = 0
        if (index(keywrd, "RESIDUES") == 0) then
          do i = 1, numat
            j = labels(i)
            if (j /= 1 .and. (j < 6 .or. j > 8)) cycle
!
!  Check H, C, N, and O, only
!
            if (txtatm(i)(14:14) >= 'a' .and. txtatm(i)(14:14) <= 'z') &
              txtatm(i)(14:14) = char(ichar(txtatm(i)(14:14)) + ichar("A") - ichar("a"))
            if(elemnt(j)(2:2) /= txtatm(i)(14:14) .and. elemnt(j)(2:2) /= txtatm(i)(13:13)) then
              num1 = char(Int(log10(i + 0.5)) + ichar("2")) 
              write(line,'(a,i'//num1//',a,a)')"Atom name for atom", i, &
              " ("//elemnt(j)(2:2)//") does not match its label """//txtatm(i)//""""
              write(iw,'(10x,a)')trim(line)
              l = l + 1
              if (l == 20) then
                write(iw,'(/10x,a)') "Remaining errors not printed"
                exit
              end if
            end if
          end do
        end if
        if (l > 0) then
          call mopend("Atom name does not match atom label. To suppress this error, add GEO-OK")
          return
        end if
      end if
!
!   WRITE OUT THE INTERATOMIC DISTANCES
!
      rmin = 100.D0 
      l = 0 
      do i = 1, numat 
        do j = 1, i 
          l = l + 1 
          rxyz(l) = sqrt((coord(1,i) - coord(1,j))**2 + &
                         (coord(2,i) - coord(2,j))**2 + &
                         (coord(3,i) - coord(3,j))**2) 
          if (.not.(rmin > rxyz(l) .and. i /= j .and. (nat(i) < 103 .or. nat(j) < 103))) cycle  
          iminr = i 
          jminr = j 
          rmin = rxyz(l) 
        end do 
      end do 
      call setcup
      if (moperr) return
!
!  The "or" here is needed because symmetry theory uses connectivity
!  (If "or" is replaced by "and" a really weird bug is introduced. DON'T DO IT!
!
      if (index(keywrd, " ADD-H") == 0 .or. index(keywrd, " PDBOUT") /= 0 ) then
        call set_up_dentate
        call check_cvs(.false.)
        call check_H(i)
      end if
      if (use_ref_geo) call big_swap(0,1)
      if (id == 0) then
        if ((index(keywrd, " ADD-H") == 0 .and. .not. mozyme) .or. numat < 200) then
          call symtrz (c, pdiag, 1, .FALSE.) 
          if (moperr) return  
          if (mode /= 1) write (iw, '(2/''      MOLECULAR POINT GROUP   :   '',A4)') name 
        end if
      end if
      mol_weight = 0.D0 
      do i = 1, numat 
        mol_weight = mol_weight + atmass(i) 
      end do 
      n9 = 0
      n4 = 0
      n1 = 0
      do i = 1, natoms
        j = labels(i)
        if (j > 0 .and. j /= 99 .and. j /= 107) then
          k = natorb(j)
          if (k == 1) then
            n1 = n1 + 1
          else if (k == 4) then
            n4 = n4 + 1
          else if (k == 9) then
            n9 = n9 + 1
          end if
        end if
      end do
      ispd = n9
      if (id == 0) then
         n2elec = 2025*n9 + 100*n4 + n1 + 2025*(n9*(n9 - 1))/2 + 450*n9*n4 + 45*n9*n1 + &
     & 100*(n4*(n4 - 1))/2 + 10*n4*n1 + (n1*(n1 - 1))/2 + 10       
      else
         n2elec = 2025*n9 + 100*n4 + n1 + 2025*(n9*(n9 + 1))/2 + 450*n9*n4 + 45*n9*n1 + &
     & 100*(n4*(n4 + 1))/2 + 10*n4*n1 + (n1*(n1 + 1))/2 + 10
      end if
      if (.not. mozyme .and. n2elec < 0) then
        if (Index (keywrd, " RESEQ") + Index (keywrd, " SITE=") + index(keywrd, " ADD-H") &
          + index(keywrd, " 0SCF") == 0) then
          call mopend(" Data set '"//trim(jobnam)//"' exists, but is too large to run.")
          write(iw,'(10x,a)')"(Maximum number of two-electron integrals allowed: 2147483647)"
          write(iw,'(10x,a, i10,a)')"(Number of two-electron integrals exceeded this by:", &
            2147483647 + n2elec,")"
          return
        end if
      end if
      if (ispd > 0 .and. index(keywrd, " ESP") /= 0) then
        line = " The ESP method does not work with 'd' orbitals"
        write(iw,'(//,10x,a)')trim(line)
        call mopend(trim(line))
        if (method_pm6) then
          write(iw,'(10x,a)')" (If possible, use MNDO, AM1, or PM3, with a 1SCF to prevent geometry changes)"
        end if
        return
      end if
!
!  Set all solid-state variables at this point
!
      lm61 = 45*n9 + 10*n4 + n1
      norbs = nlast(numat) 
      if (id > 0) then
        if (norbs < 18 .and. n9 > 0 .or. norbs < 8 .and. n4 > 0) then
          if (n9 > 0) then
            line = "Minimum number of orbitals allowed in this system is 18"
          else
            line = "Minimum number of orbitals allowed in this system is 8"
          end if
          write(iw,'(//,10x,a)')trim(line)
          call mopend(trim(line))
          return
        end if
      end if
      if (index(keywrd, " GRAPH") /= 0) then
        if (norbs > 9500) then
          write(line,'(a)')" The system is too large for 'GRAPH' to be used."
          write(iw,'(//10x,a,/)')trim(line)
          call mopend(trim(line))
          return
        end if
      end if
      mpack = (norbs*(norbs+1))/2
      if (index(keywrd," 0SCF") /= 0) mpack = 1
!
!   NOW TO CALCULATE THE NUMBER OF LEVELS OCCUPIED
!
      sing = index(keywrd,' SING') /= 0 
      doub = index(keywrd,' DOUB') /= 0 
      trip = index(keywrd,' TRIP') /= 0 
      quar = index(keywrd,' QUAR') /= 0 
      quin = index(keywrd,' QUIN') /= 0 
      sext = index(keywrd,' SEXT') /= 0 
      sept = index(keywrd,' SEPT') /= 0 
      octe = index(keywrd,' OCTE') /= 0 
      none = index(keywrd,' NONE') /= 0 
      exci = index(keywrd,' EXCI') /= 0 
      birad = exci .or. index(keywrd,'BIRAD')/=0 
      if (index(keywrd,'C.I.')/=0 .and. uhf) then 
        write (iw, '(2/10X,''C.I. NOT ALLOWED WITH UHF '')')
        call mopend ('C.I. NOT ALLOWED WITH UHF')  
        return  
      endif 
!
! NOW TO WORK OUT HOW MANY ELECTRONS ARE IN EACH TYPE OF SHELL
!
      nalpha = 0 
      nbeta = 0 
!
!      PROTECT DUMB USERS FROM DUMB ERRORS!
!
      nelecs = nint(max(elecs,0.D0)) 
      nelecs = min(2*norbs,nelecs) 
      if (.not. rhf .and. .not. uhf) then
        if (mod(nelecs,2) == 0) then
          rhf = .true.
        else
          uhf = .true.          
        end if
      end if
      if (uhf .and. index(keywrd," 0SCF") == 0) then
        if (index(keywrd," POLAR") /= 0) then
            write(iw,'(//10x,a)')" Keyword POLAR used with an odd-electron system."
            write(iw,'(10x,a)')" By default, UHF is used for odd-electron systems."
            write(iw,'(10x,a)')" POLAR does not work with UHF."
            write(iw,'(10x,a)')" To correct this fault, add RHF to the keyword line, and re-run."
            call mopend("POLAR used with a UHF system.  Add keyword RHF and re-run")
            return
          end if
          if (index(keywrd," ESP") /= 0) then
            write(iw,'(//10x,a)')" Keyword ESP used with a UHF calculation."
            write(iw,'(10x,a)')" By default, UHF is used for odd-electron systems."
            write(iw,'(10x,a)')" ESP does not work with UHF."
            write(iw,'(10x,a)')" To correct this fault, add RHF to the keyword line, and re-run."
            call mopend("ESP used with a UHF system.  Add keyword RHF and re-run")
            return
          end if
      end if
      odd = (mod(nelecs,2) /= 0)
!
! Explicit (named) declaration of spin state
!
      msdel = 100000 ! Set to a very high value, so if it is explicitly set, that will be obvious
      if (sing) then 
        if (odd) then 
         call mopend (&
             'SINGLET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' SINGLET STATE CALCULATION'')') 
          msdel = 0
        endif 
      else if (doub) then 
        if (.not. odd) then 
         call mopend (&
             'DOUBLET SPECIFIED WITH EVEN NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' DOUBLET STATE CALCULATION'')') 
          msdel = 1
        endif 
      else if (trip) then 
        if (odd) then 
         call mopend (&
             'TRIPLET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' TRIPLET STATE CALCULATION'')') 
          msdel = 2
        endif 
      else if (quar) then 
        if (.not. odd) then 
          call mopend (&
             'QUARTET SPECIFIED WITH EVEN NUMBER OF ELECTRONS, CORRECT FAULT')
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' QUARTET STATE CALCULATION'')') 
          msdel = 3
        endif 
      else if (quin) then 
        if (odd) then 
          call mopend (&
             'QUINTET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' QUINTET STATE CALCULATION'')') 
          msdel = 4
        endif 
      else if (sext) then 
        if (.not. odd) then 
          call mopend (&
             'SEXTET SPECIFIED WITH EVEN NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' SEXTET STATE CALCULATION'')') 
          msdel = 5
        endif 
      else if (sept) then 
        if (odd) then 
          call mopend (&
             'SEPTET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT') 
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' SEPTET STATE CALCULATION'')') 
          msdel = 6
        endif 
      else if (octe) then 
        if (.not. odd) then 
          call mopend (&
             'OCTET SPECIFIED WITH EVEN NUMBER OF ELECTRONS, CORRECT FAULT')  
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' OCTET STATE CALCULATION'')') 
          msdel = 7
        endif 
      else if (NONE) then 
        if (odd) then 
          call mopend (&
             'NONET SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT') 
          return  
        else 
          if (mode /= 1) write (iw, '(2/'' NONET STATE CALCULATION'')') 
          msdel = 8
        endif 
      endif 
!
!  Generic declaration of spin state
!
      i = Index (keywrd, " MS")
      if (i /= 0) then
        if (msdel < 100000) then
          write(iw,'(a)')" MS cannot be used with other keywords that define spin"
          write(iw,'(a)')" Specify either MS or other keyword, but not both."
          call mopend("MS cannot be used with other keywords that define spin")
          call to_screen(" MS cannot be used with other keywords that define spin")
          call to_screen(" Specify either MS or other keyword, but not both.")
          return
        end if
        i = Nint (2*reada (keywrd, Index (keywrd, " MS")))
        if((odd .eqv. mod(i,2) == 0) .and. index(keywrd,' 0SCF') == 0) then
          write(iw,"(10x,a,i5)")" Number of electrons in system =", nelecs
          write(iw,"(10x,a,f5.1)")" Value of MS supplied =", i*0.5d0
          write(iw,"(10x,a)")"Correct the error and re-submit"
          call mopend ("Value of MS not consistent with number of electrons")
          return
        end if
        i = index(keywrd,' MS') 
        if (i /= 0   .and. index(keywrd,' 0SCF') == 0) then 
          msdel = nint(2.d0*reada(keywrd,index(keywrd,' MS')))
          if (Mod(nelecs+msdel,2) == 1) then
            write (iw, '(//10 x, "Impossible value of MS ")')
            call mopend ("Impossible value of MS")
            return
          end if
        endif 
      end if
      if (msdel > 99999) msdel = 0 ! Set default value of msdel
      
!
! At this point, msdel is known. Now work out UHF and RHF quantities
!
      if (uhf) then
        nbeta = (nelecs - msdel)/2 
        nalpha = nelecs - nbeta
        nopen = 0
        nclose = 0
        fract = 0.d0
        nalpha_open = nalpha
        nbeta_open = nbeta
        if (mode /= 1 .and. .not. mozyme) write (iw, &
    '(2/10X,''UHF CALCULATION, NO. OF ALPHA ELECTRONS ='',I5,/27X, &
    & ''NO. OF BETA  ELECTRONS ='',I5)') nalpha, nbeta 
        i = index(keywrd,'OPEN(') 
        if (i /= 0) then 
          j = index(keywrd(i:i+10),',') + i - 1 
          nmos = nint(reada(keywrd,j)) 
          ielec = nint(reada(keywrd,index(keywrd,'OPEN(') + 5))  
          fract = (ielec*1.d0)/nmos 
          if (nclose < 0   .and. index(keywrd,' 0SCF') == 0) then 
            write (iw, '(A)') ' IMPOSSIBLE NUMBER OF FILLED SHELLS' 
            call mopend ('IMPOSSIBLE NUMBER OF FILLED SHELLS') 
            return   
          end if
        else
          ielec = 0
          nmos = 0
        endif 
         nalpha = nalpha - ielec
         nalpha_open = nalpha + nmos
         nbeta_open = nbeta
      else  !  The RHF option
!
!   NOW TO DETERMINE OPEN AND CLOSED SHELLS
!
        ielec = 0 
        nmos = 0 
        fract = 0
        if (exci .or. birad) then 
          if (mode /= 1) then 
            if (birad) write (iw, '(2/'' SYSTEM IS A BIRADICAL'')') 
            if (exci) write (iw, '(2/'' EXCITED STATE CALCULATION'')') 
          endif 
          if (odd) then 
            write (iw, &
      '(2/10X,''SYSTEM SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT '')') 
            call mopend (&
               'SYSTEM SPECIFIED WITH ODD NUMBER OF ELECTRONS, CORRECT FAULT') 
            return  
          endif          
          ielec = 2 
          nmos = 2 
        else if ((nelecs/2)*2 /= nelecs) then 
          ielec = 1 
          nmos = 1 
        endif 
       
        i = index(keywrd,'OPEN(') 
        if (i /= 0) then 
          j = index(keywrd(i:i+10),',') + i - 1 
          nmos = nint(reada(keywrd,j)) 
          ielec = nint(reada(keywrd,index(keywrd,'OPEN(') + 5)) 
        endif 
        nclose = nelecs/2 
        nopen = nelecs - nclose*2 
        if (ielec /= 0 .and. norbs > 0) then 
          if (((nelecs/2)*2==nelecs .neqv. (ielec/2)*2==ielec)  .and. index(keywrd,' 0SCF') == 0) then 
            write (iw, &
              '('' IMPOSSIBLE NUMBER OF OPEN SHELL ELECTRONS'')') 
            write (iw, '(A,I5)') ' NUMBER OF ELECTRONS IN SYSTEM:    ', nelecs 
            write (iw, '(A,I5)') ' NUMBER OF ELECTRONS IN OPEN SHELL:', ielec 
            call mopend ('IMPOSSIBLE NUMBER OF OPEN SHELL ELECTRONS') 
            return  
          endif 
          nclose = nclose - ielec/2 
          nopen = nmos 
          if (nclose + nopen > norbs) then 
            call mopend (&
       'NUMBER OF DOUBLY FILLED PLUS PARTLY FILLED LEVELS GREATER THAN TOTAL NUMBER OF ORBITALS') 
            return  
          endif 
          if (nmos == 0) then 
            write (iw, *) ' NUMBER OF M.O.s IN OPEN(M,N) IS ZERO!'
            call mopend (' NUMBER OF ELECTRONS IN OPEN(M,N) IS ZERO!')   
            return  
          endif 
          fract = ielec*1.D0/nmos 
          if (nclose < 0   .and. index(keywrd,' 0SCF') == 0) then 
            write (iw, '(A)') ' IMPOSSIBLE NUMBER OF CLOSED SHELLS' 
            call mopend ('IMPOSSIBLE NUMBER OF CLOSED SHELLS') 
            return  
          endif 
          if (mode /= 1) write (iw, '(/6X,''THERE ARE'',I5,'' DOUBLY FILLED LEVELS'')') nclose 
        endif 
        if (mode /= 1 .and. .not. mozyme) then
          num1 = char(Int(log10(nclose + 0.5)) + ichar("2")) 
          write (iw, '(2/6X,''RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS ='',I'//num1//',/)') nclose 
          if (nopen/=0 .and. abs(fract-1.D0)<1.D-4) &
            write (iw, '(/23X,''NO. OF SINGLY OCCUPIED LEVELS ='',I5)') nopen 
          if (nopen/=0 .and. abs(fract-1.D0)>1.D-4) &
            write (iw, '(/23X,''NO. OF LEVELS WITH OCCUPANCY'',F6.3,''  ='',I3)') fract, nopen 
        end if
        i = index(keywrd,'C.I.=(') 
        if (i /= 0) then 
          j = index(keywrd(i:i+10),',') + i - 1 
          ndoubl = nint(reada(keywrd,j)) 
          if (ndoubl > nclose) then 
            write (iw, '(7x,A,I3)') 'NUMBER OF DOUBLY FILLED LEVELS IN C.I. RESET TO', nclose 
            ndoubl = nclose 
          endif 
          i = nint(reada(keywrd,index(keywrd,'C.I.=(') + 5)) 
          i = i - ndoubl 
          if (nopen > i) then 
            write (iw, '(2/,'' NUMBER OF OPEN-SHELLS ALLOWED IN C.I. IS LESS '',/,    &
      &''    THAN THAT SPECIFIED BY OTHER KEYWORDS'')') 
            call mopend ('NUMBER OF OPEN-SHELLS ALLOWED IN C.I. IS LESS THAN THAT SPECIFIED BY OTHER KEYWORDS') 
            return  
          endif 
          if (i + nclose > norbs) then 
            write (iw, '(A,/,A,I3)') &
      ' NUMBER OF M.O.s REQUESTED IN C.I. IS GREATER THAN THE NUMBER OF ORBITALS' 
            call mopend (&
       'NUMBER OF M.O.s REQUESTED IN C.I. IS GREATER THAN THE NUMBER OF ORBITALS') 
            return  
          endif 
        endif 
        nopen = nopen + nclose 
      endif 
      if (index(keywrd, " CIS") /= 0) then
        msdel = 0
        else
!
!  WORK OUT IF DEFINED SPIN-STATE ALLOWED
!
        
      end if
      if (msdel/=0 .and. .not. uhf) then 
!
!   MSDEL = NUMBER OF ALPHA ELECTRONS - NUMBER OF BETA ELECTRONS
!
        ndoubl = 99 
        i = index(keywrd,'C.I.=(') 
        if (i /= 0) then 
          j = index(keywrd(i:i+10),',') + i - 1 
          ndoubl = nint(reada(keywrd,j)) 
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=(') + 5)) 
        else if (index(keywrd,'C.I.=') /= 0) then 
          nmos = nint(reada(keywrd,index(keywrd,'C.I.=') + 5)) 
        else 
          nmos = nopen - nclose 
          nmos = min0(norbs,nmos) 
        endif 
        if (ndoubl == 99) then 
          j = max(min((nclose + nopen + 1)/2 - (nmos - 1)/2, norbs - nmos + 1),1) 
        else 
          j = nclose - ndoubl + 1 
        endif 
        ne = int(max(0.D0,nclose - j + 1.D0)*2.D0 + &
          max(0.D0,(nopen - nclose)*fract) + 0.5D0) 
        nupp = (ne + msdel)/2
        ndown = ne - nupp 
        if (nmos == 0) then
          nupp = 0
          ndown = 0
        end if
!
!  NUPP  = NUMBER OF ALPHA ELECTRONS IN ACTIVE SPACE
!  NDOWN = NUMBER OF BETA  ELECTRONS IN ACTIVE SPACE
!
        if (nupp*ndown < 0 .or. nupp > nmos .or. ndown > nmos .or. nmos == 0) then 
          write (iw, '(A)')' SPECIFIED SPIN COMPONENT NOT SPANNED BY ACTIVE SPACE' 
          write (iw, '(a,i3)')' Number of alpha electrons in active space:', nupp 
          write (iw, '(a,i3)')'  Number of beta electrons in active space:', ndown 
          write (iw, '(a,i3)')'                      Size of active space:', nmos 
          if (Index (keywrd, " MS") /= 0) then
            call web_message(iw,"ms.html")
          else
            call web_message(iw,"active_space.html")
          end if
          call mopend ('SPECIFIED SPIN COMPONENT NOT SPANNED BY ACTIVE SPACE') 
          return  
        endif 
      endif 
      halfe = (nopen > nclose .and. Abs(fract -2.D0) > 1.d-20 .and. Abs(fract) > 1.d-20 &
      & .or. index(keywrd,'C.I.')/=0) 
      if (halfe) halfe = (.not. (index(keywrd,'EXCI') /= 0 .or. &
                                 index(keywrd,'ROOT') /= 0 .and. index(keywrd,'ROOT=1') == 0)) 
      if (halfe .and. id /= 0 .and. index(keywrd,' NOANCI') + index(keywrd,' 0SCF') == 0) then 
        write (iw, *) 
        write (iw, *) ' ''NOANCI'' MUST BE USED FOR RHF OPEN-SHELL SYSTEMS' 
        write (iw, *) ' THAT INVOLVE TRANSLATION VECTORS'  
        call mopend (&
       '"NOANCI" MUST BE USED FOR RHF OPEN-SHELL SYSTEMS THAT INVOLVE TRANSLATION VECTORS') 
        return  
      endif 
      yy = dble(kharge)/(norbs + 1.D-10) 
      pdiag(:norbs) = 0.D0 
      do i = 1, numat 
        ni = nat(i) 
        if (nlast(i) - nfirst(i) == (-1)) cycle  
        l = nfirst(i) - 1 
        if (nlast(i) - nfirst(i) == 0) then 
!
!    Hydrogen
!
          l = l + 1 
          pdiag(l) = tore(ni) - yy 
        else if (nlast(i) - nfirst(i) == 3) then 
!
!    Normal heavy atom
!
          w = tore(ni)*0.25D0 - yy 
          pdiag(l+1:4+l) = w 
        else 
!
!   This atom has a 'd' shell
!
          if (ni<21 .or. ni>30 .and. ni<39 .or. ni>48 .and. ni<57) then 
!
!   Main Group Element:  The "d" shell is formally empty.
!
            w = tore(ni)*0.25D0 - yy 
            pdiag(l+1:4+l) = w 
            l = 4 + l 
            pdiag(l+1:5+l) = -yy 
          else if (ni < 99) then 
!
!   Transition metal
!
            sum = tore(ni) - 9*yy 
!   First, put 2 electrons in the 's' shell
            l = l + 1 
            pdiag(l) = max(0.D0,min(sum,2.D0)) 
            sum = sum - 2.D0 
            if (sum > 0.D0) then 
!
!   Now put as many electrons as possible into the 'd' shell
!
              l = l + 3 
              do j = 1, 5 
                l = l + 1 
                pdiag(l) = max(0.D0,min(sum*0.2D0,2.D0)) 
              end do 
              sum = sum - 10.D0 
              if (sum > 0) then 
!
!   Put the remaining electrons in the 'p' shell
!
                l = l - 8 
                pdiag(l+1:3+l) = sum/3.D0 
              endif 
            endif 
          endif 
        endif 
      end do 
      call setup_nhco(ii)
!
!  Does the system contain a N bonded to three atoms, at least
!  two of which are not hydrogen atoms?
!
      do i = 1, numat
        if (nat(i) == 7 .and. nbonds(i) == 3) then
        j = 0
        if (nat(ibonds(1,i)) == 1) j = 1
        if (nat(ibonds(2,i)) == 1) j = j + 1
        if (nat(ibonds(3,i)) == 1) j = j + 1
        if ( j < 2)  N_3_present = .true.
        end if
        if (N_3_present) exit
      end do
!
!  Does the system contain a Si-O-H group?
!
      do i = 1, numat
        if (nat(i) == 14) then
          do j = 1, nbonds(i)
            if (nat(ibonds(j,i)) == 8) then
              l = ibonds(j,i)
              do k = 1, nbonds(l)
                if (nat(ibonds(k,l)) == 1) Si_O_H_present = .true.
              end do
            end if
          end do
        end if
        if (Si_O_H_present) exit
      end do
  !    Si_O_H_present = .false.
      if (nvar == 0) then
        if (index(keywrd," GRAD") /= 0 .and. index(keywrd, " 0SCF") == 0) then
          if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE") == 0) then
            line = " Keyword GRADIENTS used, but geometry has no variables."
            call mopend(trim(line))
            return
          end if
        end if
      end if
      if (mod(nelecs, 2) == 1 .and. .not. uhf .and. nnhco > 0 .and. mode /= 1 &
      .and. (index(keywrd," 1SCF") == 0 .or. index(keywrd," GRAD") /= 0) .and. &
      nvar > 10 .and. index(keywrd," GEO-OK") == 0 .and. abs(kharge) < 400) then
      write(iw,"(/)")
       write(iw,"(10x,a)") &
       " This system has an odd number of electrons and also", &
       " contains peptide linkages.  The CPU time required for", &
       " calculating the gradients is likely to be large. If ""UHF"" is ", &
       " added, then the gradient calculation will run much faster. ", &
       " Either add ""GEO-OK"" or ""UHF"" to the keyword list.", &
       " Also, check that the system should be a radical. A quick " , &
       " way to do this is to add keyword ""LEWIS"" and run it again.", &
       " This will print out the Lewis structure."
       call mopend("radical run with RHF.  Either add 'GEO-OK' or 'UHF'.")
       return
      end if
      if (nnhco > 3 .and. &
        index(keywrd," GEO-OK") == 0 .and. &
        index(keywrd," RESEQ") == 0 .and. &
        index(keywrd," 0SCF") == 0 .and. &
        index(keywrd," CHARGES") == 0 .and. &
        index(keywrd," 1SCF") == 0) then
        j = 1
        do i = 1, natoms
          if (na(i) > 0) j = j + 1
        end do
        if (nvar > 0 .and. j == natoms .and. .not. is_PARAM) then
          if (index(keywrd," 0SCF") + index(keywrd," ADD-H") + index(keywrd," SITE=") == 0) then
            write(iw,"(5x,a)")"The system contains peptide linkages, but all atoms are in internal coordinates.", &
            "This is likely to result in problems in geometry optimization.", &
            "To correct this, add ""XYZ"" to the keyword line.", &
            "(If internal coordinates should be used, add ""GEO-OK"" to the keyword line.)"
            call mopend("Peptides should be run using Cartesian coordinates")
            return 
          end if
        end if
      end if
      if (mode /= 1 .and. ii > 1) then 
        write (iw, '(A,I4,2A)') ' THERE ARE', ii/2, ' PEPTIDE LINKAGES IDENTIFIED IN THIS SYSTEM' 
        write (iw, '(A)') ' Keyword "NOMM" has been used, therefore a Molecular Mechanics correction will not be used' 
      else if (mode /= 1 .and. nnhco /= 0 .and. (method_RM1 .or. method_pm6 .or. method_PM7) &
      & .and. (index(keywrd, "MMOK") == 0)) then 
        if (index(keywrd, " LEWIS") + index(keywrd, " 0SCF") + &
          index(keywrd, " RESEQ") + index(keywrd, " CHARGES")== 0) then
          if (method_RM1) line = "RM1"
          if (method_PM6) line = "PM6"
          if (method_PM7) line = "PM7"
      !      write (iw, '(A)') ' When peptide bonds are present, and '//line(:3)//&
      !      &' is used, then keyword "MMOK" or "NOMM" must be used'
      !      call mopend("Keyword MMOK or NOMM must be used with "//line(:3)//" for this system")
      !      return
          end if
      endif            
      if (mode /= 1 .and. index(keywrd,'PRTINT') /= 0) then  
        write (iw, '(2/10X,''  INTERATOMIC DISTANCES'')') 
        call vecprt (rxyz, numat) 
      endif 
      if (rmin < 0.9D0 .and. index(keywrd,'CHECK') /= 0 .or. &
          rmin < 0.2D0 .and. index(keywrd,'GEO-OK') == 0 .or. &
          rmin < 1.d-4) then 
        icount = 0 
        ireal = iminr 
        jreal = jminr 
        do i = 1, natoms 
          if (labels(i)==99 .or. labels(i)==107) cycle  
          icount = icount + 1 
          if (icount == iminr) ireal = i 
          if (icount /= jminr) cycle  
          jreal = i 
        end do
        if (index(keywrd,' 0SCF') == 0) then
          write(line,"('ATOMS',i6,' AND',i6,' ARE SEPARATED BY',f7.4,' ANGSTROMS.')")ireal, jreal, rmin
          if (index(keywrd,'CHECK') /= 0) then
            call mopend(trim(line))
            write(iw,'(/10x,a)')'FAULT DETECTED BY KEYWORD "CHECK'
            if (log) write (ilog, '(//10x,a,//10x,a)')trim(line), 'FAULT DETECTED BY KEYWORD "CHECK'          
          else
            call mopend(trim(line))
            write(iw,'(/10x,a)')'TO CONTINUE CALCULATION SPECIFY "GEO-OK"'
            num1 = char(Int(log10(max(ireal,jreal)     + 1.0)) + ichar("1")) 
            if (pdb_label) then
              write(iw,"(/10x,a,i"//num1//",a)")"(Label for atom ", ireal, ": """//txtatm(ireal)//""")"   
              write(iw,"( 10x,a,i"//num1//",a)")"(Label for atom ", jreal, ": """//txtatm(jreal)//""")"  
            end if
            if (log) write (ilog, '(//10x,a,//10x,a)')trim(line), 'TO CONTINUE CALCULATION SPECIFY "GEO-OK"'   
          end if         
          write (iw, '(/,a)')'   NOTE THAT THE ATOM NUMBERS CORRESPOND TO THE'//&
            ' ''CARTESIAN COORDINATES'' ATOM LIST.'
          if (index(keywrd,' 0SCF') == 0) then 
            if (index(keywrd,' ADD-H') /= 0) then 
              write(line,'(a,i6,a,i6)')'SEVERE ERROR IN GEOMETRY, FIX FAULT BEFORE CONTINUING. Atoms:',ireal, ' and', jreal
              call mopend (trim(line))   
              numat = 0
              natoms = 0
              return
            end if
            if (rmin < 1.d-4) then 
              write(line,'(a,i6,a,i6)')'GEOMETRY IN ERROR, FIX FAULT BEFORE CONTINUING. Atoms:',ireal, ' and', jreal
              call mopend (trim(line))
            else
              call mopend ("GEOMETRY IN ERROR.  TO CONTINUE CALCULATION SPECIFY 'GEO-OK'.")
            end if
            return  
          endif 
        end if
      endif 
      if (.not.debug) return  
      write (iw, 290) numat, norbs, ndorbs, natoms 
  290 format('   NUMBER OF REAL ATOMS:',i4,/,'   NUMBER OF ORBITALS:  ',i4,/,&
        '   NUMBER OF D ORBITALS:',i4,/,'   TOTAL NO. OF ATOMS:  ',i4) 
      write (iw, 300) (uspd(i),i=1,norbs) 
  300 format('   ONE-ELECTRON DIAGONAL TERMS',/,10(/,10f8.3)) 
      write (iw, 310) (pdiag(i),i=1,norbs) 
  310 format('   INITIAL P FOR ALL ATOMIC ORBITALS',/,10(/,10f8.3)) 
      return  
      end subroutine moldat 
      
subroutine setcup 
   !***********************************************************************
   !
   !   SETCUP determines CUTOFP and the number of unit cells in each
   !          direction that will be needed in order to allow CUTOFP
   !   to be correctly used.
   !
   !***********************************************************************
    use molkst_C, only: cutofp, id, keywrd, l1u, l2u, l3u, l123, l11, l21, l31, line
    use common_arrays_C, only : tvec
    use chanel_C, only: iw
    use reada_I
    use volume_I
    implicit none
    integer :: i
    double precision :: area12, area13, area23, r1, r12, r13, r2, &
         & r23, r3, tv1, tv2, tv3, vol, sum
    cutofp = 1.d10
    l1u = 0
    l2u = 0
    l3u = 0
    l123 = 1
    if (id == 0) return
    i = Index (keywrd, " CUTOFP")
    if (i /= 0) then
      cutofp = reada (keywrd, i+7)
    else
      cutofp = 30.d0
    end if
   !
   !   CALCULATE L1U, L2U, and L3U
   !
    if (id == 1) then  !   Polymer case
      !
      !   TV1 =  length of polymer.
      !
      tv1 = Sqrt (tvec(1, 1)**2+tvec(2, 1)**2+tvec(3, 1)**2)
      if (tv1 < 1.d0) then
        line = "  Length of translation vector is too small."
        write(iw,'(a)')trim(line)
        call mopend (trim(line))
        return
      end if
      l1u = Int (cutofp*4.d0/3.d0/tv1) + 1
    else if (id == 2) then !  Layer system
      !
      !   TV1, TV2 = Distances across unit cell.
      !
      r1 = Sqrt (tvec(1, 1)**2 + tvec(2, 1)**2 + tvec(3, 1)**2)
      r2 = Sqrt (tvec(1, 2)**2 + tvec(2, 2)**2 + tvec(3, 2)**2)
      if (r1 < 1.d0) then
        line = "  Length of first translation vector is too small."
        write(iw,'(a)')trim(line)
        call mopend (trim(line))
        return
      end if
      if (r1 < 1.d0) then
         line = "  Length of second translation vector is too small."
        write(iw,'(a)')trim(line)
        call mopend (trim(line))
        return
      end if
      r12 = Sqrt ((tvec(1, 2)-tvec(1, 1))**2 + (tvec(2, 2)-tvec(2, 1))**2 &
           & + (tvec(3, 2)-tvec(3, 1))**2)
      tv1 = r1 * Sin (Acos((r1**2 + r2**2 - r12**2)/(2*r1*r2)))
      tv2 = r2 * Sin (Acos((r1**2 + r2**2 - r12**2)/(2*r1*r2)))
      l1u = Int (cutofp*4.d0/3.d0/tv1) + 1
      l2u = Int (cutofp*4.d0/3.d0/tv2) + 1
    else !  Solid-state (three-dimensional crystal)
      !
      !   TV1,TV2,TV3 = Distances between faces of unit cell
      !
      r1 = Sqrt (tvec(1, 1)**2 + tvec(2, 1)**2 + tvec(3, 1)**2)
      r2 = Sqrt (tvec(1, 2)**2 + tvec(2, 2)**2 + tvec(3, 2)**2)
      r3 = Sqrt (tvec(1, 3)**2 + tvec(2, 3)**2 + tvec(3, 3)**2)
      if (r1 < 1.d0 .or. r2 < 1.d0 .or. r3 < 1.d0) then
        if (r1 < 1.d0) then
          line = "  Length of first translation vector is too small."
          write(iw,'(a)')trim(line)
          call mopend (trim(line))
          return
        else if (r2 < 1.d0) then
          line = "  Length of second translation vector is too small."
          write(iw,'(a)')trim(line)
          call mopend (trim(line))
          return
        else 
          line = "  Length of third translation vector is too small."
          write(iw,'(a)')trim(line)
          call mopend (trim(line))
          return
        end if
      end if
      r12 = Sqrt ((tvec(1, 2)-tvec(1, 1))**2 + (tvec(2, 2)-tvec(2, 1))**2 &
           & + (tvec(3, 2)-tvec(3, 1))**2)
      r13 = Sqrt ((tvec(1, 3)-tvec(1, 1))**2 + (tvec(2, 3)-tvec(2, 1))**2 &
           & + (tvec(3, 3)-tvec(3, 1))**2)
      r23 = Sqrt ((tvec(1, 3)-tvec(1, 2))**2 + (tvec(2, 3)-tvec(2, 2))**2 &
           & + (tvec(3, 3)-tvec(3, 2))**2)
      area12 = r1 * r2 * Sin (Acos((r1**2 + r2**2 - r12**2)/(2*r1*r2)))
      area13 = r1 * r3 * Sin (Acos((r1**2 + r3**2 - r13**2)/(2*r1*r3)))
      area23 = r2 * r3 * Sin (Acos((r2**2 + r3**2 - r23**2)/(2*r2*r3)))
      vol = volume (tvec, id)
      if (vol < 1.d0) then
        write (iw,'(//10X,A,F9.6,A)') &
             & "Volume of unit cell unreasonably small:", vol, &
             & " Cubic Angstroms"
        call mopend("Volume of unit cell unreasonably small")
        return
      end if
      tv1 = vol / area23
      tv2 = vol / area13
      tv3 = vol / area12
      sum = 4.0d0
      if ((tv1 < sum .or. tv2 < sum .or. tv3 < sum) .and. Index (keywrd, " GEO-OK") == 0) then
        write(iw,"(//10x,a,f9.3,a   )") "Translation vector length 1:", r1," Angstroms"
        write(iw,"(  10x,a,f9.3,a   )") "Translation vector length 2:", r2," Angstroms"
        write(iw,"(  10x,a,f9.3,a,/ )") "Translation vector length 3:", r3," Angstroms"  
        if (tv1 < sum) write (iw,'(10X,A,F6.3,A,f6.3)') &
             & "Distance between faces 2 and 3 is unreasonably small:", tv1," Angstroms, min:", sum
        if (tv2 < sum) write (iw,'(10X,A,F6.3,A,f6.3)') &
             & "Distance between faces 1 and 3 is unreasonably small:", tv2," Angstroms, min:", sum
        if (tv3 < sum) write (iw,'(10X,A,F6.3,A,f6.3)') &
             & "Distance between faces 1 and 2 is unreasonably small:", tv3," Angstroms, min:", sum
        if (tv1 < sum .or. tv2 < sum .or. tv3 < sum) then
          call mopend("One or more translation vectors are unreasonably small")
        write(iw,'(10x,a)')"To over-ride this safety check, add keyword 'GEO-OK'"
        return
        end if
      end if
      l1u = Int ((cutofp*4.d0/3.d0)/tv1) + 1
      l2u = Int ((cutofp*4.d0/3.d0)/tv2) + 1
      l3u = Int ((cutofp*4.d0/3.d0)/tv3) + 1    
    end if
    l123 = (2*l1u + 1)*(2*l2u + 1)*(2*l3u + 1)
    l11 = min(l1u,1)
    l21 = min(l2u,1)
    l31 = min(l3u,1)
end subroutine setcup
subroutine write_cell(iprt)
  use molkst_C, only: mol_weight, escf, numat, keywrd, gnorm, &
   line, mers, gui
  use funcon_C, only: a0, ev, fpc_9, fpc_10
  use common_arrays_C, only: nat, tvec
  use ef_C, only : nstep
  use reada_I
  use volume_I
  implicit none
  integer, intent(in) :: iprt
  integer :: i, j, k, l, z, m, old_nstep = -1
  double precision :: ta, tb, tc, tab, tbc, tac, talpha, tbeta, tgamma, vol
  integer, dimension (100) :: nel
  save :: old_nstep
    if (iprt < 0 .or. gui) return
  if (iprt == 0) then
    if(old_nstep == nstep) return
    old_nstep = nstep
  end if

!
!  Write out unit cell lengths and angles
!
  ta = Sqrt(tvec(1,1)**2 + tvec(2,1)**2 + tvec(3,1)**2)
  tb = Sqrt(tvec(1,2)**2 + tvec(2,2)**2 + tvec(3,2)**2)
  tc = Sqrt(tvec(1,3)**2 + tvec(2,3)**2 + tvec(3,3)**2)
  tab = Sqrt((tvec(1,1)-tvec(1,2))**2 + (tvec(2,1)-tvec(2,2))**2 + (tvec(3,1)-tvec(3,2))**2)
  tac = Sqrt((tvec(1,1)-tvec(1,3))**2 + (tvec(2,1)-tvec(2,3))**2 + (tvec(3,1)-tvec(3,3))**2)
  tbc = Sqrt((tvec(1,3)-tvec(1,2))**2 + (tvec(2,3)-tvec(2,2))**2 + (tvec(3,3)-tvec(3,2))**2)
  talpha = 57.295779513d0*Acos((tb**2 + tc**2 - tbc**2)/(2*tc*tb))
  tbeta  = 57.295779513d0*Acos((ta**2 + tc**2 - tac**2)/(2*ta*tc))
  tgamma = 57.295779513d0*Acos((ta**2 + tb**2 - tab**2)/(2*ta*tb))
  tab = 1.d0
  if (index(keywrd," BCC") /= 0) tab = 2.d0
!
!  Work out the empirical formula
!
  i = Index(keywrd, " Z=")
  if (i /= 0) then
    i = Nint(reada(keywrd,i))
    z = mers(1)*mers(2)*mers(3)*i
    if (Index (keywrd, " BCC") /= 0) z = z/2
  else
    nel = 0
    do i = 1, numat
      nel(nat(i)) = nel(nat(i)) + 1
    end do
    j = 0
    do i = 1, 100
      if (nel(i) > 0) then
        j = j + 1
        nel(j) = nel(i)
       end if
    end do
    k = 10000
    do i = 1, j
      if(nel(i) < k) k = nel(i)
    end do
!
!  k is the smallest number of atoms of any element in the formula
!
    do i = 1, 20
      m = 0
      do l = 1, j
        if (Abs((i*nel(l))/k - (i*1.d0*nel(l))/k) > 1.d-5) m = 1 
      end do
      if (m == 0) exit
    end do
!
!  Number of empirical units  = k/i
!
    z = k/i
  end if
  vol = volume (tvec, 3)
  if (mers(1) > 0 .and. mers(2) > 0 .and. mers(3) > 0) then
    write (line,'(i4,a,3f7.3,3f7.2,a,f8.2,a,f6.3,a,f10.3, a, f7.2)')nstep + 1," a, b, c, alpha, beta, gamma:", &
    & tab*ta/mers(1), tab*tb/mers(2), tab*tc/mers(3), talpha, tbeta, tgamma, &
    " Vol:",vol/(mers(1)*mers(2)*mers(3))," Density:",mol_weight * 1.d24 / fpc_10 / vol, &
    " HoF:",escf/z," Grad:",gnorm/sqrt(z*1.d0)
    write(iprt,"(a)")line(:len_trim(line))
  end if       
  return
end subroutine write_cell
subroutine write_unit_cell_HOF(iprt)
  use molkst_C, only : keywrd, escf, numat, line, mers
  use common_arrays_C, only : nat
  use reada_I
  use to_screen_I
  implicit none
  integer, intent(in) :: iprt
  integer :: i, j, k, l, z, nel(107), m
    if (iprt < 0) return
!
!  Work out the empirical formula
!
        i = Index(keywrd, " Z=")
        if (i /= 0) then
          i = Nint(reada(keywrd,i))
          z = mers(1)*mers(2)*mers(3)*i
          if (Index (keywrd, " BCC") /= 0) z = z/2
        else
          nel = 0
          do i = 1, numat
            nel(nat(i)) = nel(nat(i)) + 1
          end do
          j = 0
          do i = 1, 100
            if (nel(i) > 0) then
              j = j + 1
              nel(j) = nel(i)
            end if
          end do
          k = 1000
          do i = 1, j
            if(nel(i) < k) k = nel(i)
          end do
!
!  k is the smallest number of atoms of any element in the formula
!
          do i = 1, 10
            m = 0
            do l = 1, j
              if (Abs((i*nel(l))/k - (i*1.d0*nel(l))/k) > 1.d-5) m = 1 
            end do
            if (m == 0) exit
          end do
!
!  Number of empirical units  = k/i
!
         z = k/i
       end if
        
       write (line, "(10x,'H.o.F. per unit cell    =',f17.5,' KCAL, for',i4,' unit cells')") &
         & escf/z, z 
       if (iprt == 0) then
         call to_screen(line)
       else
         write(iprt,"(a)")line(:len_trim(line))
       end if
      end subroutine write_unit_cell_HOF
      subroutine write_pressure(iprt)
      use common_arrays_C, only: loc, tvec, na
      use molkst_C, only: nvar, pressure, line, press
      use funcon_C, only: fpc_10
      use common_arrays_C, only : grad, xparam, labels
      use volume_I
      use to_screen_I
      implicit none
      integer, intent(in) :: iprt
      integer :: m, i1,i2, k, l, i, ndim
      double precision :: xi
      double precision, dimension (nvar) :: dsum, dsum1
      double precision, external :: ddot
      if (iprt < 0) return  
        m = 0
        i1 = 0
        i2 = 0
        ndim = 0
        do i = 1, nvar
          k = loc(1, i)
          l = labels(k)
          xi = xparam(i)
          if (l == 107 .and. (m == 0 .or. k == i1)) then
!
!  Atom is a Tv
!
            i1 = k
            m = m + 1
            dsum(m) = grad(i)
            dsum1(m) = xparam(i)
            if (m == 3) then
              if (na(k) /= 0) then
                if (iprt == 0) then
                  call to_screen("The pressure required to constrain translation vectors")
                  call to_screen("can only be calculated if Cartesian coordinates are used.")
                else
                  if (Abs(pressure) > 0.01d0) then
                    write(iprt,*)"The pressure required to constrain translation vectors"
                    write(iprt,*)"can only be calculated if Cartesian coordinates are used."
                  end if
                end if                
                return
              end if
!
!  Determine the scalar of the component of the gradient vector in the 
!  direction of the translation vector
!
              xi = ddot(3,dsum, 1,dsum1,1)
!
!  Convert this into a pressure = force per unit area
!
              xi = -(4184.d0*10.d0**30)/fpc_10 * xi/volume (tvec, 3)
              if (Abs(xi) < 1.d-20) cycle ! suppress printing if gradients are zero
              if (i2 == 0) then
                write(line,'(a)') "          Pressure required to constrain translation vectors"
                if (iprt == 0) then
                  call to_screen(trim(line))
                else
                  write(iprt,*)trim(line)
                end if   
                i2 = 1
              end if
              xi = (xi - pressure * (4184.d0*10.d0**30)/ fpc_10)*1.d-9
              ndim = ndim + 1
              press(ndim) = xi
              write(line,'(10x,a,i4,a,f7.2,a)')"Tv(", k,")  Pressure:",xi," GPa"
               if (iprt == 0) then
                  call to_screen(trim(line))
                else
                  write(iprt,*)trim(line)
                end if   
              m = 0
            end if
          end if
        end do   
  end subroutine write_pressure
  subroutine setup_nhco(ii)
  USE molmec_C, only : nnhco, nhco, htype
  USE molkst_C, only : numat,  method_am1, method_pm3, method_mndo, method_pm6, method_PM7, method_rm1, &
   keywrd
  use common_arrays_C, only : nat
  implicit none
  integer, intent (out) :: ii
!
!  Local
!
  integer :: j, i, k, l, m, jj
  double precision, external :: distance

  nnhco = 0 
  !
  !   SET UP MOLECULAR-MECHANICS CORRECTION TO -(C=O)-(NH)- LINKAGE
  !   THIS WILL BE USED IF NOMM HAS NOT BEEN SPECIFIED.
  !
                      htype = 0.d0
    if (method_mndo) htype = 6.1737D0 
    if (method_am1)  htype = 3.3191D0 
    if (method_pm3)  htype = 7.1853D0 
    if (method_rm1)  htype = 2.4127D0 
    if (method_pm7)  htype = 3.1595D0 
    if (method_pm6)  htype = 2.5000D0 
    ii = 0 
    if (index(keywrd,'NOMM') /= 0) ii = 1 
  !
  !   IDENTIFY O=C-N-H SYSTEMS VIA THE INTERATOMIC DISTANCES MATRIX
  !
  l230: do j = 1, numat 
      if (nat(j) /= 6) cycle  l230 
      do i = 1, numat 
        if (nat(i) /= 8) cycle  
        if (distance(i, j) > 1.3D0) cycle  
        do k = 1, numat 
          if (nat(k) /= 7) cycle  
          if (distance(k, j) > 1.6D0) cycle  
          do l = 1, numat 
            if (nat(l) /= 1) cycle  
            if (distance(k, l) > 1.3D0) cycle  
  !
  !   WE HAVE A H-N-C=O SYSTEM.  THE ATOM NUMBERS ARE L-K-J-I
  !   NOW SEARCH OUT ATOM ATTACHED TO NITROGEN, THIS SPECIFIES
  !   THE SYSTEM X-N-C=O
  !
            l190: do m = 1, numat 
              if (m==k .or. m==l .or. m==j) cycle  l190 
              if (distance(m, k) > 1.7D0) cycle  l190 
              do jj = 1, nnhco, 2 
                if (nhco(3,jj) /= k) cycle  
                cycle  l190 
              end do 
              nnhco = nnhco + 1 
              nhco(1,nnhco) = i 
              nhco(2,nnhco) = j 
              nhco(3,nnhco) = k 
              nhco(4,nnhco) = m 
              nnhco = nnhco + 1 
              nhco(1,nnhco) = i 
              nhco(2,nnhco) = j 
              nhco(3,nnhco) = k 
              nhco(4,nnhco) = l 
              if (ii /= 0) then 
                ii = ii + 2 
                nnhco = nnhco - 2 
              endif 
              cycle  l230 
            end do l190 
          end do 
        end do 
      end do 
  end do l230 
    end subroutine setup_nhco




