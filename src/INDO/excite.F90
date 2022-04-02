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

! ____________________________________________________________________________
! _____________________ ARBITRARY - EXCITATION DET GENERATOR ___________________
! ____________________________________________________________________________

      subroutine excite (istsym, e, ci)

!     **** This generates spin - adapted determinates as SCF excitations	****
!     **** Uses fast code for single excitations from SCF reference det	****
!     **** or uses completely general code for Complete Active Space, 	****
!     **** single or double excitations from any number of reference	****
!     **** determinates.						****
!     **** The routine makes two passes: first, do determine the CI	****
!     **** state identity, second to generate the MP2 perturbation	****
!     **** corrections to all of the CI matrix elements.		****
      use reimers_C, only:  n, nb2, matind, nbt, nel, nshell, norb, &
          norbl, norbh, vnn, ixprd, nmrep, stwt, nptg, nr, ns, nslwr, &
          nsupr, istr, nsym, mofrag, icifrag, nconf, multci, nci, &
          ncore, nol, noh, nvl, nvh, nov, fastci, emaxci, nese, nciout, &
          icisym, nmp2, edef, ecore, eec, ee2, &
          ifst, afst, aor1, bor1, isc, tr1, nfact, vv, ifbinc, krefd, &
          aor, bor, aod, bod, ao1, bo1, aos, bos, cc0, gamma, beta, tot, &
          nbtmo, aocc, bocc, spintr, evalmo, nspn, nex, dbls, mspn, nold, &
          nvhd, nolc, nvhc, mrci

      USE molkst_C, only : norbs, nopen, nclose
      USE chanel_C, only : iw
      use cosmo_C, only  : useps

      implicit none
!     **** MP2 storage in excite and spcls ****
      double precision  e(nconf), ci(matind(nconf) + nconf + 1), &
                        ener, tel, xj12, xk12
      integer           vvs(4), vvd(4), a1, a, b, b1, maxunp, &
                        istsym(nconf), &
                        i1, icifrag1, iirefd, irefd, ismin, ispn, &
                        isynca, isyncb, &
                        j1, jconf, jconfm, &
                        jmin, maxoc, &
                        nalpha, nbeta, ndif, nelopn, nend, norbopn, nov2, nov4, &
                        nperm,  nv, nvs, nvd
      double precision, allocatable :: aii(:)
      logical           exok, maxmp2
      logical*1, allocatable :: aoc(:), boc(:)
      character*8       param
      integer           i, j, k
      double precision  ediel
      logical*1, external :: permca
      param = " "
      if (param == "xxxxxxxx") return
      ecore = 0.d0

!     **** determine n! ****
      n = norbs
      nb2 = n*(n + 1)/2
      icisym(1) = nconf
      emaxci = 1.d6

      nfact(0) = 1
      do i = 1, 10
        nfact(i) = nfact(i - 1) * i
      end do

      maxunp = multci
      if (dbls) maxunp = maxunp + 1
      if (mrci) maxunp = maxunp + min(noh - nolc + 1, nvhc - nvl + 1)
      mspn = nfact(maxunp*2) / nfact(maxunp)**2

!     **** determine basis - vector type of each MO ****
!     **** if type j = s, px, py, etc is present, then nbtmo bit j is set to 1 ****
      if(allocated(nbtmo)) deallocate(nbtmo)
      allocate(nbtmo(n))
      do i = 1, n
        do j = 0, 8
          ifbinc(j) = 0
        end do
        do k = 1, n
          if (abs(cc0(i, k)) > 1.d-5) ifbinc(nbt(k)) = 1
        end do
        nbtmo(i) = 0
        do j = 8, 0, - 1
          nbtmo(i) = nbtmo(i)*2 + ifbinc(j)
        end do
      end do
!     **** apply symmetry weights to MO coeffs ****

      if (allocated(stwt))  deallocate(stwt)
      if (allocated(isc))   deallocate(isc)
      if (allocated(istr))  deallocate(istr)
      if (allocated(nsym))  deallocate(nsym)
      allocate(stwt(n))
      allocate(isc(n))
      allocate(istr(8, n))
      allocate(nsym(n))

!     ************** debug: use this to turn off symmetry info ***********
      nptg = 1
      nr = 1
      ns(1) = 1
      nslwr(1) = 1
      nsupr(1) = n
      do i = 1, n
        stwt(i) = 1.D0
        nsym(i) = 1
        istr(1, i) = i
      end do

      do i1 = 1, n
!	**** stwt is set to nstr**(1/4) ****
        stwt(i1) = sqrt (1.D0/stwt(i1))
        isc(i1) = 0
      end do
      do i1 = 1, n
        i = istr(1, i1)
        if (isc(i).eq.0) then
          isc(i) = 1
          do j = 1, n
            cc0(j, i) = cc0(j, i) * stwt(i1)
          end do
        end if
      end do

!     ********** determine if it's a default or a specified CI **********

      fastci = .not. dbls .and. .not. mrci .and. (nopen.eq.nclose)

!     ************ initial setup of integral storage ***********
      nov = nvh - nol + 1
      nov2 = nov*(nov + 1)/2
      nov4 = nov2*(nov2 + 1)/2

      if(allocated(eec)) deallocate(eec)
      if(allocated(ee2)) deallocate(ee2)
      allocate(eec(nov2))
      allocate(ee2(nov4))

      if (fastci) then
        ifst = 0
      else

!       ********** calc energy of core of electrons (with nucleii) ************
!       **** code replaced in SCF routine by evaluating Fock matrix given *****
!       **** the density from the CI core only ... much much faster! **********

        tot = tot + vnn
        ecore = vnn
        do i = 1, ncore
          call inttr (i, i, i, i, cc0, gamma, nbtmo, xj12, xk12, 0)
          ecore = ecore + 2.D0*beta(matind(i) + i) + xj12
          do j = 1, i - 1
            call inttr (i, i, j, j, cc0, gamma, nbtmo, xj12, xk12, 0)
            ecore = ecore + 4.D0*xj12 - 2.D0*xk12
          end do
        end do
        ecore = ecore - tot

!       ******* set to default value storage for all possible integrals ******
!       **** eec(i, a) = sum_j - 2(i, a|j, j) + (i, j|a, j) with j in core and	****
!       **** i, a in x - space; matrix stored in symmetric mode.		****
!       **** eec(i, j, a, b) = (i, j|a, b) with i, j, a, b in x - space.		****

        edef = 1.D20
        j = nov*(nov + 1)/2
        do i = 1, j
          eec(i) = sngl(edef)
        end do
        j = j*(j + 1)/2
        do i = 1, j
          ee2(i) = sngl(edef)
        end do

      end if

      jconf = 0
      nmp2 = 0
      jconfm = 0
      do k = 1, nr
        nci(k) = 0
      end do
!     ******************** form determinants ***************************
!     **** loop over reference determinates; first read alpha occupancy ****z
      if(allocated(vv))  deallocate(vv)
      if(allocated(aor)) deallocate(aor)
      if(allocated(bor)) deallocate(bor)
      if(allocated(aor1)) deallocate(aor1)
      if(allocated(bor1)) deallocate(bor1)
      if(allocated(aos)) deallocate(aos)
      if(allocated(bos)) deallocate(bos)
      if(allocated(aod)) deallocate(aod)
      if(allocated(bod)) deallocate(bod)
      if(allocated(ao1)) deallocate(ao1)
      if(allocated(bo1)) deallocate(bo1)
      allocate(vv(n))
      allocate(aor(nov))
      allocate(bor(nov))
      allocate(aor1(nov))
      allocate(bor1(nov))
      allocate(aoc(nov))
      allocate(boc(nov))
      allocate(aos(nov))
      allocate(bos(nov))
      allocate(aod(nov))
      allocate(bod(nov))
      allocate(ao1(nov, mspn))
      allocate(bo1(nov, mspn))

      irefd = 0

      do while (irefd.eq.0)
        irefd = irefd + 1
        maxmp2=.false.

!       **** open shell contributions ****
        tel = 0.D0
        do i = 2, nshell
          tel = tel + nel(i)
        end do
        nelopn = nint(tel)
        nalpha = (nelopn + multci - 1) / 2
        nbeta = nelopn - nalpha

        if (nbeta < 0) then
          nbeta = 0
          nalpha = nelopn
        end if

        nv = 1
        vv(1) = norb(1) + nalpha

        maxoc = vv(nv)

!	**** alpha set occupied/vacant flags ****
        do i = nol, vv(1)
          aor(i - ncore) = .true.
        end do
        do i = vv(1) + 1, nvh
          aor(i - ncore) = .false.
        end do
        do j = 2, nv
          aor(vv(j) - ncore) = .true.
        end do
        write (iw, "(//, ' Reference determinate nber', i4, ' is ', a5, ' 1 -', 25i4)") &
          irefd, 'ALPHA', (vv(i), i=1, nv)

!	**** determine beta occupancy ****
!	**** auto include ground state reference det ****
        vv(1) = norb(1) + nbeta
        maxoc = max (maxoc, vv(nv))

!	**** beta set occupied/vacant flags ****
        do i = nol, vv(1)
          bor(i - ncore) = .true.
        end do
        do i = vv(1) + 1, nvh
          bor(i - ncore) = .false.
        end do
        do j = 2, nv
          bor(vv(j) - ncore) = .true.
        end do
! Subtract dielectric energy off of ecore
        if (useps) then
          call staticsolv(aor, bor, ediel)
          ecore = ecore - ediel
        end if

        if (irefd.eq.1) then
!	  **** save first reference determinate for excitation naming ****
          do i = 1, nov
            aor1(i) = aor(i)
            bor1(i) = bor(i)
          end do
        end if

        if(allocated(ao1))    deallocate(ao1)
        if(allocated(bo1))    deallocate(bo1)
        if(allocated(aocc))   deallocate(aocc)
        if(allocated(bocc))   deallocate(bocc)
        if(allocated(nspn))   deallocate(nspn)
        if(allocated(krefd))  deallocate(krefd)
        if(allocated(spintr)) deallocate(spintr)
        if(allocated(aii))    deallocate(aii)
        if(allocated(tr1))    deallocate(tr1)
        allocate(ao1(nov, mspn))
        allocate(bo1(nov, mspn))
        allocate(aocc(nov, mspn, nex))
        allocate(bocc(nov, mspn, nex))
        allocate(nspn(nconf))
        allocate(krefd(nconf))
        allocate(spintr(mspn, nconf))
        allocate(aii(n))
        allocate(tr1(mspn))

        do i=1, n
          aii(i) = evalmo(i)
        end do

        do i=1, nov
          do j=1, mspn
            do k=1, nex
              aocc(i, j, k) = .False.
              bocc(i, j, k) = .False.
            end do
          end do
        end do

!	**** spin adapt and include in config list ****

        call spclas(istsym, jconf, aor, bor, ao1, bo1, e, &
     &        aocc, bocc, nspn, spintr, beta, cc0, gamma, nbtmo, aii, nfact, &
     &        irefd, krefd, .true.)
!        write (6, *) '1   ', jconf, e(jconf), aor, '   ', bor

!	********* read in excitation control parameters ********
        nvs = 4
        vvs(1) = nol
        vvs(2) = noh
        vvs(3) = nvl
        vvs(4) = nvh

!	**** set up values for double excitations if no DOUBLE card read ****
        nvd = 4
        if (dbls) then
          vvd(1) = nold
          vvd(2) = noh
          vvd(3) = nvl
          vvd(4) = nvhd
        else
          vvd(1) = 1
          vvd(2) = 0
          vvd(3) = 1
          vvd(4) = 0
        end if
!	**** set up values for CAS if no CAS card read ****
!	**** CAS amongst all open shell orbitals ****
        norbopn = norbh(nshell) - norbh(1)
        if (norbopn.eq.0) then
!         **** RHF, turn CAS off ****
          nv = 1
          vv(1) = maxoc + 1
        else
          nv = norbopn
          j = 0
          do i = norbl(2), norbh(nshell)
            j = j + 1
            vv(j) = i
          end do
        end if

        if (mrci) then
          fastci = .false.
          nv = 0
          do i = nolc, nvhc
            nv = nv + 1
            vv(nv) = i
          end do

!         **** ensure CAS orbitals are in ascending order for sort speed ***
          do i = 1, nv - 1
            a = n
            jmin = 0
            do j = i, nv
              if (vv(j) < a) then
                a = vv(j)
                jmin = j
              end if
            end do
            vv(jmin) = vv(i)
            vv(i) = a
          end do
          if (vv(1).le.ncore .or. vv(nv) > nvh) &
              stop 'in EXCITE: CAS RANGE'
        end if

!       **** check and print excitation generators ****

        if (nv > 1) write (iw, " (' CAS excitations amongst orbitals:', 25i4)") (vv(i), i=1, nv)
        if (vvs(2) > 0) then
          write (iw, "(1x, a6, ' excitations FROM orbs', i4, ' to', i4, ' INTO orbs', i4, ' to', i4)") &
            'SINGLE', (vvs(i), i=1, nvs)
          call chkdeg (vvs, nvs, aii)
        end if
        if (vvd(2) > 0) then
          write (iw, "(1x, a6, ' excitations FROM orbs', i4, ' to', i4, ' INTO orbs', i4, ' to', i4)") &
            'DOUBLE', (vvd(i), i=1, nvd)
          call chkdeg (vvd, nvd, aii)
        end if
        write (iw, *)

!       **** first, CAS excitations amongst a finite nber of orbitals ****

!       **** determine nber of alpha and beta els in space of nv orbs	****
!       **** and hence the nber of permutations possible		****
        nalpha = 0
        nbeta = 0
        do i = 1, nv
          if (aor(vv(i) - ncore)) nalpha = nalpha + 1
          if (bor(vv(i) - ncore)) nbeta = nbeta + 1
        end do
        nend = nint (2.D0**nv)
        isynca = nend

        do while (permca (aor, aoc, nalpha, vv, nv, isynca))
          isyncb = nend
          do while (permca (bor, boc, nbeta, vv, nv, isyncb))

!         **** spin adapt and include in config list ****
            call spclas(istsym, jconf, aoc, boc, ao1, bo1, e, &
     &          aocc, bocc, nspn, spintr, beta, cc0, gamma, nbtmo, aii, nfact, &
     &          irefd, krefd, .false.)
!            write (6, *) 'CAS ', jconf, e(jconf), aoc, '   ', boc

!	    ***** generate single excitations from all CAS dets *****

            do i1 = vvs(2), vvs(1), - 1
              ifst = i1
              i = i1 - ncore
              do a1 = vvs(3), vvs(4)
                afst = a1
                a = a1 - ncore
!               **** this is the i - >a excitation ****

                if (aoc(i) .and. .not. aoc(a) .or.&
     &              boc(i) .and. .not. boc(a) .or. i.eq.a) then
!                 **** fill in occupancy from CAS det ****
                  do k = 1, nov
                    aos(k) = aoc(k)
                    bos(k) = boc(k)
                  end do
!                 **** select alpha spin if available, else beta spin ****
                  if (aos(i) .and. .not. aos(a)) then
                    aos(i) = .false.
                    aos(a) = .true.
                  else if (i.ne.a) then
                    bos(i) = .false.
                    bos(a) = .true.
                  end if
!	          **** spin adapt and include in config list ****
                  call spclas(istsym, jconf, aos, bos, ao1, bo1, e, &
     &              aocc, bocc, nspn, spintr, beta, cc0, gamma, nbtmo, aii, nfact, &
     &              irefd, krefd, .false.)
!                  write (6, *) 'CAS2', jconf, e(jconf), aos, '   ', bos

!	          ******** double excitations *********

                  do j1 = vvd(2), vvd(1), - 1
                    j = j1 - ncore
                    do b1 = vvd(3), vvd(4)
                      b = b1 - ncore
!                     **** this is the i - >a & j - >b double excitation ****

!                     **** eliminate double counting so that J>=I and B>=A ****
                      exok = .true.
                      if (vvs(1).le.j .and. j.le.vvs(2) .and.&
     &                    vvs(3).le.b .and. b.le.vvs(4) )&
     &                    exok = j.ge.i .and. a.ge.b
                      if (exok .and. j.eq.a) then
!                       **** eliminate possible double counting ****
!                       **** proceed only if no single excitation i - >b ***
                        exok = (.not. aoc(i) .or. aoc(b)) .and.&
     &                        (.not. boc(i) .or. boc(b) .or.&
     &                        i < vvs(1) .or. i > vvs(2)  .or.&
     &                        b < vvs(3) .or. b > vvs(4) )
                      end if
                      if (exok .and. i.eq.b)&
!		        **** proceed only if no single excitation j - >a ***
     &                  exok = (.not. aoc(j) .or. aoc(a)) .and.&
     &                        (.not. boc(j) .or. boc(a) .or.&
     &                        j < vvs(1) .or. j > vvs(2) .or.&
     &                        a < vvs(3) .or. a > vvs(4) )

                      if (exok .and. ( aos(j) .and. .not. aos(b) .or.&
     &                                 bos(j) .and. .not. bos(b) ) ) then
!		        **** fill in occupancy from single excited det ****
                        do k = 1, nov
                          aod(k) = aos(k)
                          bod(k) = bos(k)
                        end do

!                       **** do only one of alpha or beta promotion ****
                        if (aod(j) .and. .not.aod(b)) then
                          aod(j) = .false.
                          aod(b) = .true.
                        else
                          bod(j) = .false.
                          bod(b) = .true.
                        end if

!                       **** spin adapt and include in config list ****
                        call spclas(istsym, jconf, aod, bod, ao1, bo1, &
     &                      e, aocc, bocc, nspn, spintr, beta, &
     &                      cc0, gamma, nbtmo, aii, nfact, irefd, krefd, .false.)
!                        write (6, *) 'CAS3', jconf, e(jconf), aod, '   ', bod

!                       **** end of double excitation ****
                      end if
                    end do
                  end do

!                 **** end of single excitation ****
                end if
              end do
            end do

!	    **** end of CAS excitation ****
          end do
        end do
!       **** end of processing of excitation cards ****
!	**** end of reference det, position for next ALPHA read ****
!        if (param.eq.'   ALPHA') backspace ix
        if (maxmp2) then
          jconfm = jconfm + 1
          do while (jconfm.le.jconf .and. krefd(jconfm).ne. - 1)
            jconfm = jconfm + 1
          end do
        end if
      end do

!     **************** reorganize basis set order **********************

      nconf =  min (nconf, jconf)
      nciout = min (nciout, nconf)
      nese =   min (nese,  nciout)
      jconf =  nconf
      write (iw, "(/' CI excitations=', i5, ':', 8(5x, a3, '=', i3))") nconf, (nmrep(k, nptg), nci(k), k=1, nr)

      if (allocated(icifrag))  deallocate(icifrag)
      if (allocated(mofrag))   deallocate(mofrag)
      allocate(icifrag(nconf))
      allocate(mofrag(n))
      do i = 1, n
        mofrag(i) = 1
      end do

!     **** determine symmetry of states ****
!     **** determine if state is in a diagonal fragment block ****

      do i = 1, nconf
        istsym(i) = 1
        icifrag(i) = 0
        do j = 1, nov
          if (aocc(j, 1, i) .neqv. bocc(j, 1, i)) then
            istsym(i) = ixprd (istsym(i), nsym(j + ncore))
            if (icifrag(i).eq.0) then
              icifrag(i) = mofrag(j + ncore)
            else
              if (icifrag(i).ne.mofrag(j + ncore)) icifrag(i) = -1
            end if
          end if
        end do
        ndif = 0
        do j = 1, nov
          if (aocc(j, 1, i) .neqv. bocc(j, 1, i)) then
            ndif = ndif + 1
          end if
        end do
      end do

!     **** sort states by symmetry ****

      do k = 1, nconf
        jmin = 0
        ismin = 10000
        do j = k, nconf
          if (istsym(j) < ismin) then
            jmin = j
            ismin = istsym(j)
          end if
        end do
!	**** swap states, maintaining energy ordering ****
        ismin = istsym(jmin)
        ener = e(jmin)
        nperm = nspn(jmin)
        icifrag1 = icifrag(jmin)
        iirefd = krefd(jmin)
        do ispn = 1, nperm
          tr1(ispn) = spintr(ispn, jmin)
          do i = 1, nov
            ao1(i, ispn) = aocc(i, ispn, jmin)
            bo1(i, ispn) = bocc(i, ispn, jmin)
          end do
        end do

        do j = jmin, k + 1, - 1
          j1 = j - 1
          istsym(j) = istsym(j1)
          e(j) = e(j1)
          nspn(j) = nspn(j1)
          icifrag(j) = icifrag(j1)
          krefd(j) = krefd(j1)
          do ispn = 1, nspn(j1)
            spintr(ispn, j) = spintr(ispn, j1)
            do i = 1, nov
              aocc(i, ispn, j) = aocc(i, ispn, j1)
              bocc(i, ispn, j) = bocc(i, ispn, j1)
            end do
          end do
        end do

        istsym(k) = ismin
        e(k) = ener
        nspn(k) = nperm
        icifrag(k) = icifrag1
        krefd(k) = iirefd
        do ispn = 1, nperm
          spintr(ispn, k) = tr1(ispn)
          do i = 1, nov
            aocc(i, ispn, k) = ao1(i, ispn)
            bocc(i, ispn, k) = bo1(i, ispn)
          end do
        end do
      end do

!     **** count the number of CI states of each symmetry ****

      do k = 1, nr
        nci(k) = 0
      end do
      do i = 1, nconf
        k = istsym(i)
        nci(k) = nci(k) + 1
      end do

!     **** zero the CI matrix ****
      do i = 1, matind(nconf) + nconf
        ci(i) = 0.D0
      end do

      return

      end subroutine excite

!     **********************************************************************

      subroutine chkdeg (vv, nv, evalmo)

!     **** this checks an exc list for degenerate boundary MO eigenvalues ****
      use reimers_C, only: n
      USE chanel_C, only : iw

      implicit none
      integer :: nv, vv(nv)
      double precision   evalmo(n), diff, tol
      integer            i, j, k, idir, ifn

      data tol  /1.d-1/
      do j = 1, 4
        i = vv(j)
        idir = 0
        if ( (j.eq.1 .or. j.eq.3) .and. i > 1) then
          diff = abs (evalmo(i) - evalmo(i - 1))
          if (diff  <  tol) then
            idir = 1
            ifn = n
          end if
        else if ( (j.eq.2 .or. j.eq.4) .and. i < n .and. i.ge.1) then
          diff = abs (evalmo(i) - evalmo(i + 1))
          if (diff  <  tol) then
            idir = -1
            ifn = 1
          end if
        end if

        if (idir.ne.0) then
          if (diff > tol/5.) then
            write (iw, "(/' CI WARNING: MO', i4, ' EVAL DIFF <', f7.4)") i, tol
          else
            write (iw, "(/' CI WARNING: MO', i4, ' EVAL DIFF <', f7.4)") i, tol/5.
            do k = i, ifn, idir
              if (abs(evalmo(k) - evalmo(i)) < tol) vv(j) = vv(j) + idir
            end do
          end if
        end if
      end do

      return
      end subroutine chkdeg

!     **********************************************************************

      subroutine spclas (istsym, jconf, ao, bo, ao1, bo1, e, &
     &        aocc, bocc, nspn, spintr, beta, c, gamma, nbtmo, aii, nfact, &
     &        irefd, krefd, isrefd)

!     **** this takes a determinate AO, BO and finds all other DETs	****
!     **** belonging to its CLASS (ie, dets coupled by S**2).  It calcs	****
!     **** their interaction energies, diagonalizes S**2 to deduce the	****
!     **** linear combinations which are eigenvalues of spin, and hence	****
!     **** transforms the energy matrix to deduce the spin - adapted	****
!     **** CONFiguration energies E.  Those of spin MULTCI are selected	****
!     **** and the lowest NCONF stored if their energy is < EMAXCI.	****
!     **** AOCC and BOCC store the dets if at least one energy included	****
!     **** ISNG stores locn of all singly occ orbs in det		****
      use reimers_C, only: n, nb2, matind, nconf, multci, ncore, nol, &
          nvh, nov, fastci, icisym, ecore, ifst, afst, ixprd, ndiff, nex, &
          mspn, nsym
      USE chanel_C, only : iw
      use cosmo_C, only : useps

      implicit none
!     **** MP2 storage in excite and spcls ****
      integer            isng(mspn), nbtmo(n), nel(nov), nspn(nconf), &
                         istsym(nconf), krefd(nconf), nfact(0:mspn)
      double precision   s2(mspn, mspn), eval(mspn), wk(mspn), &
                         h(mspn, mspn), s2eig(mspn, mspn), s2eval(mspn), &
                         gamma(n, n), c(n, n), e(nconf), beta(nb2), &
                         aii(n), spintr(mspn, nov)
      integer            i, j, k, i0, i1, idet, ier, iisym, ih, il, ip, &
                         irefd, isync, j1, jconf, jdet, &
                         k1, k2, nalpha, nbeta, ndeg, &
                         nelj, nend, nperm, nsng
      double precision   ciener, e1, eec1, j12, k12, s1, tol

      logical*1          ao(nov), bo(nov), ao1(nov, mspn), bo1(nov, mspn), &
                         aocc(nov, mspn, nex), bocc(nov, mspn, nex)
      logical            isrefd
      double precision       ediel

!     *********** check that its of the correct symmetry ************
      iisym = 1
      do j = 1, nov
        if (ao(j) .neqv. bo(j)) iisym = ixprd (iisym, nsym(j + ncore))
      end do
      if (.not. isrefd .and. icisym(iisym).eq.0) then
        return
      end if

!     *********** check that this det isnt already included in CI list ********

!     **** calculate the total number of electrons in each orb ****
      do i = 1, nov
        nel(i) = 0
        if (ao(i)) nel(i) = 1
        if (bo(i)) nel(i) = nel(i) + 1
      end do

      do 200 i = 1, jconf
        do j = 1, nov
          nelj = 0
          if (aocc(j, 1, i)) nelj = 1
          if (bocc(j, 1, i)) nelj = nelj + 1
          if (nelj.ne.nel(j)) goto 200
        end do

!	**** can only get here if two dets are in the same spin class ****
        return
200   continue

!     ***** if its closed shell with singles excitation; fast code ****

      if (fastci) then
        if (ifst > 0) then
         call inttr (ifst, ifst, afst, afst, c, gamma, nbtmo, j12, k12, 0)

          eval(1) = aii(afst) - aii(ifst) - j12
          if (multci.eq.1) eval(1) = eval(1) + 2.D0*k12
        else
!	  **** ground state energy ****
          eval(1) = 0.D0
        end if

        nperm = 1
        h(1, 1) = 1.D0
        wk(1) = multci
        do i = 1, nov
          ao1(i, 1) = ao(i)
          bo1(i, 1) = bo(i)
        end do

        goto 500
      end if

!     **** determine nber of sngl occ alpha and beta orbitals ****

      nsng = 0
      nalpha = 0
      nbeta = 0
      do i = 1, nov
        if (nel(i).eq.1) then
          nsng = nsng + 1
          if (nsng > mspn) stop 'NSNG in SPCLAS'
          isng(nsng) = i + ncore
          if (ao(i)) then
            nalpha = nalpha + 1
          else
            nbeta = nbeta + 1
          end if
        end if
      end do
      i = abs (nalpha - nbeta) + 1

!     **** generate all possible perturbations ****

      nperm = nfact(nsng) / nfact(nalpha) / nfact(nbeta)
      if (nperm > mspn) then
        write (iw, *) 'ERROR: max nber of spin components=', mspn, &
     &     ' found=', nperm, ' for det:'
        write (iw, '(1x, 131l1)') (ao(i), i=1, nov)
        write (iw, '(1x, 131l1)') (bo(i), i=1, nov)
        stop 'in SPCLAS'
      end if

      isync = 0
      nend = nint (2.D0**nsng)
      idet = 1

!     **** PERTS generates next perturbation, is false if all done ****

      do idet = 1, nperm
        call perms (ao, bo, ao1(1, idet), &
     &              bo1(1, idet), isng, nalpha, nbeta, isync, nend)
      end do

!     **** calculate the energy matrix for all dets in this class ****

      do idet = 1, nperm
        e1 = 0.D0

        do i1 = nol, nvh
          i = i1 - ncore
          if (nel(i) > 0) then
            e1 = e1 + beta(matind(i1) + i1) * nel(i)

!	    **** matrix element <i, i_bar || i, i_bar> ****
            if (nel(i).eq.2) then
              call inttr (i1, i1, i1, i1, c, gamma, nbtmo, j12, k12, 0)
              e1 = e1 + j12
            end if

!	    **** add core terms sum_j (i, i|j, j) ****
            e1 = e1 + eec1 (i1, i1, c, gamma, nbtmo) * nel(i)

!	    **** inner loop over all non - core electrons up to current one ****
            do j1 = nol, i1 - 1
              j = j1 - ncore
              if (nel(j) > 0) then
!		**** matrix elements <ij||ij> via J = [ii||jj], K = [ij|ij] ****
                call inttr (i1, i1, j1, j1, c, gamma, nbtmo, j12, k12, 0)

                if (nel(i).eq.1) then
                  if (nel(j).eq.1) then
!		    **** only one electron in each orbital ****
                    e1 = e1 + j12
                    if (ao1(i, idet) .eqv. ao1(j, idet))&
     &                e1 = e1 - k12
                  else
!		    **** one double occ, one single ****
                    e1 = e1 + 2.D0*j12 - k12
                  end if
                else
                  if (nel(j).eq.1) then
!		    **** one double occ, one single ****
                    e1 = e1 + 2.D0*j12 - k12
                  else
!		    **** both double occ ****
                    e1 = e1 + 4.D0*j12 - 2.D0*k12
                  end if
                end if
              end if
            end do

          end if
        end do
! RMG - add solvent correction
        ediel = 0.D0
        if (useps) then
          call staticsolv(ao, bo, ediel)
        end if
        h(idet, idet) = e1 + ecore + ediel
      end do

!     ***** determine interaction energies between all dets in class ****

      do idet = 2, nperm
        do jdet = 1, idet - 1
          h(idet, jdet) = ciener (ao1(1, idet), bo1(1, idet), &
     &      ao1(1, jdet), bo1(1, jdet), c, gamma, nbtmo, beta)
          h(jdet, idet) = h(idet, jdet)
        end do
      end do

!     **** diagonalize energy block ****
      call tred2e (mspn, nperm, h, eval, wk, h)
      call tql2e  (mspn, nperm,  eval, wk, h, ier)

!     **** form S**2 matrix in det basis ****

      do i = 1, nperm
        s2(i, i) = (nalpha + nbeta)*0.5D0 + (nalpha - nbeta)**2*0.25D0
        do j = 1, i - 1
!	  **** sign of off diag elem is given by det ordering convention ****
!	  **** alpha beta alpha beta .... as used also in ALIGN function ****
!	  **** first, ensure that only 2 alpha orbs are different ****
          ndiff = 0
          do k = 1, nov
            if (ao1(k, i) .neqv. ao1(k, j)) ndiff = ndiff + 1
          end do
          s2(i, j) = 0.D0
          if (ndiff.eq.2) s2(i, j) = 1.D0
          s2(j, i) = s2(i, j)
        end do
      end do

!     **** if degenerate energy eigenvalues, then form S**2 in eigvec basis ****

      tol = 1.d-5
      il = 1
      do while (il < nperm)
        ih = il + 1
        do while (ih.le.nperm .and. abs(eval(il) - eval(ih))  <  tol)
          ih = ih + 1
        end do
        if (ih > nperm .or. abs(eval(il) - eval(ih))  >  tol) ih = ih - 1
        ndeg = ih - il + 1
        if (ndeg > 1) then
!	  **** found NDEG degenerate eigenvalues, form S**2 in evec basis ****
          i0 = il - 1
          do i = 1, ndeg
            do j = 1, i
              s2eig(i, j) = 0.D0
              do k1 = 1, nperm
                s1 = 0.D0
                do k2 = 1, nperm
                  s1 = s1 + s2(k1, k2) * h(k2, i0 + j)
                end do
                s2eig(i, j) = s2eig(i, j) + h(k1, i0 + i) * s1
              end do
            end do
          end do
!	  **** diagonalize S**2 in eigvec basis ****
          call tred2e (mspn, ndeg, s2eig, s2eval, wk, s2eig)
          call tql2e  (mspn, ndeg,      s2eval, wk, s2eig, ier)
!	  **** transform energy eigvecs to be eigvecs of S**2 also ****
          do j = 1, nperm
            do i = 1, ndeg
              wk(i) = 0.D0
              do k = 1, ndeg
                wk(i) = wk(i) + s2eig(k, i) * h(j, i0 + k)
              end do
            end do
            do i = 1, ndeg
              h(j, i0 + i) = wk(i)
            end do
          end do
        end if

        il = il + ndeg
      end do

!     **** determine diag element of transformed S**2 matrix ****

      tol = 1.d-5
      do i = 1, nperm
        wk(i) = 0.D0
        do j = 1, nperm
          s1 = 0.D0
          do k = 1, nperm
            s1 = s1 + s2(k, j) * h(k, i)
          end do
          wk(i) = wk(i) + h(j, i) * s1
        end do
!	**** use S**2 = s(s + 1) and spin mult = 2s + 1 to get multiplicity ****
        wk(i) = sqrt (1.D0 + 4.D0*wk(i))
      end do

!     **** check for non - integral multiplicity ****

      do i = 1, nperm
        if (abs(mod(wk(i) + tol, 1.D0))  >  2*tol) then
          write (iw, *) 'ERROR in S^2 eigenvalue', (sngl(wk(j)), j=1, nperm)
          do k = 1, nperm
            write (iw, *) k, (ao1(j, k), j=1, nov), '  ', (bo1(j, k), j=1, nov)
          end do
          stop 'spin in SPCLAS'
        end if
      end do

500   continue

!     **** select only eigvecs of the correct multiplicity		****
!     **** CI: store in mastor store, truncating by max nber and energy	****
!     **** PM2: evaluate interactions with chosen states		****
      do ip = 1, nperm
        if (nint(wk(ip)) .eq. multci) then

!	    *********** store the state in CI master index ***********

            call sortst (jconf, eval(ip), e, ao1, bo1, aocc, bocc, nperm, &
     &        nspn, h(1, ip), spintr, irefd, krefd, iisym, istsym, isrefd)

        end if
      end do
      return
      end subroutine spclas

!     **********************************************************************

      subroutine perms (aoi, boi, ao, bo, isng, nalpha, nbeta, isync, imax)

!     **** this forms all possible permutations of NALPHA and NBETA ****
!     **** electrons, one per orbital, in orbs ISNG		    ****
!     **** ISYNC is increased one per call; IMAX is terminating value **
      use reimers_C, only: ncore, nov

      implicit none
      integer           isng(nov), i, j, &
                        ialpha, ibeta, imax, io, isync, &
                        jsync, kk, nalpha, nbeta, ntot
      logical*1         aoi(nov), boi(nov), ao(nov), bo(nov)

!     **** copy accross initial determinate ****
      do i = 1, nov
        ao(i) = aoi(i)
        bo(i) = boi(i)
      end do

!     **** loop increasing isync looking for allowed pattern ****

100   continue

      if (isync.ge.imax) stop 'PERMS: could not find permutation'
      isync = isync + 1

      ntot = nalpha + nbeta
      ialpha = 0
      ibeta = 0
      jsync = isync
      do i = 1, ntot
        io = isng(i) - ncore
        j = jsync / 2
        kk = jsync - 2*j
        jsync = j
        ao(io) = kk.eq.0
        if (ao(io)) then
          ialpha = ialpha + 1
          if (ialpha > nalpha) goto 100
        else
          ibeta = ibeta + 1
          if (ibeta > nbeta) goto 100
        end if
      end do

!     **** it can only get here if pattern is allowed ****
      do i = 1, ntot
        io = isng(i) - ncore
        bo(io) = .not. ao(io)
      end do

      return
      end subroutine perms

!     **********************************************************************

      logical*1 function permca (aoi, ao, nalpha, icas, ntot, isync)

!     **** this forms all possible permutations of NALPHA electrons in	***
!     **** in the NTOT orbitals specified by ISNG			***
!     **** ISYNC is decreased one per call; 0 is terminating value	***
      use reimers_C, only: ncore, nov

      implicit none
      integer :: nalpha, icas, ntot, isync
      logical*1         aoi(nov), ao(nov)
      dimension         icas(nov)
      integer :: i, j, kk, io, jsync, ialpha

!     **** copy accross initial determinate ****
      do i = 1, nov
        ao(i) = aoi(i)
      end do

      permca = .false.

!     **** loop increasing isync looking for allowed pattern ****

100   continue

      isync = isync - 1
      if (isync < 0) return

      ialpha = 0
      jsync = isync
      do i = 1, ntot
        io = icas(i) - ncore
        j = jsync / 2
        kk = jsync - 2*j
        jsync = j
        ao(io) = kk.eq.0
        if (ao(io)) then
          ialpha = ialpha + 1
          if (ialpha > nalpha) goto 100
        end if
      end do
      if (ialpha.ne.nalpha) goto 100

!     **** it can only get here if pattern is allowed ****
      permca = .true.

      return
      end function permca

!     **********************************************************************

      subroutine sortst (k, ener, e, ao1, bo1, aocc, bocc, nperm, nspn, &
     &                tr, spintr, irefd, krefd, iisym, istsym, isrefd)

!     **** this keeps track of the lowest - energy NCONF determinates	****
!     **** a new state is inserted into the list in energy order	****
!     **** a new state is inserted into the list in energy order	****
      use reimers_C, only: nconf, nci, icisym, emaxci, nov, fastci, &
          nex, mspn

      implicit none
      double precision  e(nconf), tr(mspn), spintr(mspn, nconf), &
                        ener
      integer           nspn(nconf), krefd(nconf), istsym(nconf), &
                        k, nperm, irefd, iisym, &
                        i, j, iii, ispn, j1, jo
      logical           isrefd
      logical*1         aocc(nov, mspn, nex), bocc(nov, mspn, nex), &
     &                  ao1(nov, mspn), bo1(nov, mspn)

!     **** first, check the energy max cuttoff criterion ****
      if (ener > emaxci) then
        return
      end if

!     **** first check to see if max states of this symmetry exceeded ****
      if (.not. isrefd .and. nci(iisym).eq.icisym(iisym)) then
!	**** find locn of state ****
        i = k
        do while (istsym(i).ne.iisym)
          i = i - 1
        end do
        if (ener > e(i)) then
          return
        end if
!	**** remove this state ****
        k = k - 1
        nci(iisym) = nci(iisym) - 1
        do j1 = i, k
          j = j1 + 1
          e(j1) = e(j)
          istsym(j1) = istsym(j)
          nspn(j1) = nspn(j)
          krefd(j1) = krefd(j)
          do ispn = 1, nspn(j)
            spintr(ispn, j1) = spintr(ispn, j)
            do jo = 1, nov
              aocc(jo, ispn, j1) = aocc(jo, ispn, j)
              bocc(jo, ispn, j1) = bocc(jo, ispn, j)
            end do
          end do
        end do
      end if

!     **** if full, spot check to see if its energy is too large ****
      if (k.eq.nconf) then
        if (abs(ener - e(k))  <  1.d-5) then
!	  **** degenerate with last state, remove them both ****
          nci(istsym(k)) = nci(istsym(k)) - 1
          k = k - 1
          return
        else
          if (ener > e(k)) then
            return
          end if
        end if
!	**** insert new into filled list, pop top one off ****
        nci(istsym(k)) = nci(istsym(k)) - 1
        k = k - 1
      end if

!     **** see where new state fits in energy order (crap code to fool Intel debugger) ****
      i = k
      iii = 1
      do while (i.ge.1 .and. iii.eq.1)
        if (ener < e(i)) then
          i = i - 1
        else
          iii = 0
        end if
      end do
      i = i + 1

! RMG - add to avoid moving 1st config to solve triplet fastci issues
      if (fastci .and. i.eq.1 .and. k.ge.1) i = i + 1

!     **** move those above it up one spot ****
      do j = k, i, - 1
        j1 = j + 1
        e(j1) = e(j)
        istsym(j1) = istsym(j)
        nspn(j1) = nspn(j)
        krefd(j1) = krefd(j)
        do ispn = 1, nspn(j)
          spintr(ispn, j1) = spintr(ispn, j)
          do jo = 1, nov
            aocc(jo, ispn, j1) = aocc(jo, ispn, j)
            bocc(jo, ispn, j1) = bocc(jo, ispn, j)
          end do
        end do
      end do

!     **** insert new state ****
      e(i) = ener
      istsym(i) = iisym
      nspn(i) = nperm
      krefd(i) = irefd

      do ispn = 1, nperm
        spintr(ispn, i) = tr(ispn)
        do jo = 1, nov
          aocc(jo, ispn, i) = ao1(jo, ispn)
          bocc(jo, ispn, i) = bo1(jo, ispn)
        end do
      end do

      nci(iisym) = nci(iisym) + 1
      k = k + 1

      return
      end subroutine sortst
