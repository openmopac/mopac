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

subroutine chklew (mb, numbon, l, large, debug)
!***********************************************************************
!
!   CHKLEW  identifies all the obvious Lewis structural elements such
!           as bonds, lone pairs, and ions.
!
!   On exit, NUMBON(1):  Number of sigma bonds.
!            NUMBON(2):  Number of lone pairs.
!            NUMBON(3):  Number of pi-bonds.
!
!***********************************************************************
    use molkst_C, only: natoms, numat, keywrd, norbs, line, maxtxt
    use chanel_C, only: iw
    use elemts_C, only: elemnt
    use parameters_C, only: tore, ndelec, natorb
    use common_arrays_C, only : txtatm, labels, nat, nbonds, ibonds, &
      nfirst, nlast, pibonds
    use MOZYME_C, only: icharges, ions, Lewis_elem, Lewis_tot, Lewis_max, ib, iz
    implicit none
    logical, intent (in) :: debug
    integer, intent (in) :: large
    integer, intent (out) :: l
    integer, dimension (3), intent (inout) :: numbon
    integer, dimension (numat), intent (inout) :: mb
!
    logical, external :: arom, arom2
    logical :: big, graphi, lnext, first_pi
    integer :: i, i1, i2, ii, iset, j, jj, k, loop, m, nbii, ni_loc, &
         & npi, alloc_stat, not_aromatic_pi
    character :: plus_seven*1, num1*1, num2*1
    integer, dimension(:), allocatable :: ipi, ir5, mpii
    intrinsic Index, Min, Nint
    double precision, external :: distance

!
!   Set LET to .TRUE.  At some later date, have LET defined by a keyword
!   If LET is true, then the job will not stop if there is a very large
!   number of ions.
!
    if (allocated(Lewis_elem)) deallocate(Lewis_elem)
    Lewis_max = max(norbs,1)
    allocate(Lewis_elem(2,Lewis_max))
    Lewis_tot = 0
    big = (Index (keywrd, " LARGE") /= 0)
    graphi = .false.
    first_pi = .true.
    allocate (ipi(6*numat), ir5(numat), mpii(numat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("chklew")
      goto 1100
    end if
    do i = 1, numat
!
!  Do not allow more than 4 bonds per atom.  Later on,
!  more can be added, if necessary.
!
      ib(i) = nlast(i) - nfirst(i) + 1
      iz(i) = Nint (tore(nat(i)))
    end do
!   Construction of the Lewis structure is based on a set of "rules",
!   These are, in order of usage:
!
! *  Every atom that is bonded to another atom forms a sigma bond to
!    that atom.
! *  The number of lone pairs on an atom is equal to the difference
!    between the number of unused valence electrons and unused orbitals.
! *  For unsaturated systems:
!    Open-ended pi-bonds are created when an atom pi-bonds to exactly
!    one other atom.
!    Pi bonds attached to open-ended pi-bonds are created next.
!    If no open-ended pi-bonds can be created, a delocalized pi-bond is
!    assigned.  Only atoms that pi-bond to exactly two other atoms are
!    involved here. Pi-bonds in five membered rings are NOT allowed, if
!    the ring is attached to one or more six-membered rings.
!
!   If any atoms have explicit charges, assign them first.
!
    j = 0
    do i = 1, natoms
      if (labels(i) /= 99) then
        j = j + 1
      end if
      plus_seven = txtatm(i)(7:7)
      txtatm(i)(7:7) = " "
      if (Index (txtatm(i)(:6), "+") /= 0) then
        icharges = icharges + 1
        iz(j) = iz(j) - 1
        ions(j) =  1
      else if (Index (txtatm(i)(:12), "-") /= 0) then
        icharges = icharges + 1
        iz(j) = iz(j) + 1
        ions(j) = -1
      end if
      txtatm(i)(7:7) = plus_seven
    end do
!
!   Fill "d" shell if necessary
!
    do i = 1, numat
      if (ib(i) /= nlast(i) - nfirst(i) + 1) then
        j = nat(i)
        if (ndelec(j) /= 0) then
!
!   NDELEC(I) is the number of "d" electrons NOT involved in the
!   Lewis structure.
!
          do k = 1, ndelec(j) / 2
            call add_Lewis_element(i,0,0, numbon(2))
          end do
        end if
      end if
    end do
!
    if (debug) then
      write (iw,*) &
           & " Arrays before making Sigma Framework (ATOM CORE ORBS BONDS)"
      j = Min (numat, large)
      do k = 1, j, 20
        l = Min (k+19, j)
        write (iw, "(20I4)") (i, i=k, l)
        write (iw, "(20I4)") (iz(i), i=k, l)
        write (iw, "(20I4)") (ib(i), i=k, l)
        write (iw, "(20I4)") (nbonds(i), i=k, l)
      end do
    end if
    numbon = 0
!
!   First set: The SIGMA bond framework
!
    do ii = 1, numat
      mb(ii) = 0
      nbii = nbonds(ii)
      do i = 1, nbii
        jj = ibonds(i, ii)
        if (jj >= ii .and. ib(ii) > 0 .and. ib(jj) > 0) then
!
! A sigma bond exists if:
!  (a) The atoms are bonded together.
!  (b) The second atom has a higher atom number than the first atom
!      (This prevents "double counting")
!  (c) Both atoms have unused atomic orbitals available
!
          call add_Lewis_element(ii,jj,0, numbon(1))
        end if
      end do
    end do
    do i = 1, numat
      if (ib(i) /= 0) go to 1000
    end do
    l = 0
    goto 1070
1000 continue
    if (debug) then
      write (iw,*) &
           & " Core charges after making Sigma Framework (ATOM CORE ORBS)"
      j = Min (numat, large)
      do k = 1, j, 20
        l = Min (k+19, j)
        write (iw, "(20I4)") (i, i=k, l)
        write (iw, "(20I4)") (iz(i), i=k, l)
        write (iw, "(20I4)") (ib(i), i=k, l)
      end do
    end if
    if (index(keywrd, " SETPI") /= 0) then
!
!  The user has supplied some pi bonds.
!
      j = 0
      do loop = 1, numat
        ii = pibonds(loop,1)
        if (ii == 0) exit
        jj = pibonds(loop,2)
        do i = 1, nbonds(ii)
          if (ibonds(i, ii) == jj) exit
        end do
        if (i > nbonds(ii)) then
          if (j == 0) write(iw, '(/10x, a)') &
            "Using ""SETPI"" a pi-bond has been defined between two atoms that are not bonded together"
          if (maxtxt == 26) then
            if (j == 0) write(iw,'(/27x, a, 27x, a)')"Atom-1", "Atom-2"
            write(iw,'(17x, a, 5x, a)')'"'//txtatm(ii)(:26)//'"', '"'//txtatm(jj)(:26)//'"'
          else
            if (j == 0) write(iw,'(/15x, a, 15x, a)')"Atom-1", "Atom-2"
            write(iw,'(14x,a2,i5, 14x,a2,i5)')elemnt(nat(ii)), ii, elemnt(nat(jj)), jj
          end if
          j = 1
        end if
        call add_Lewis_element(ii,jj,0, numbon(3))
        if (ib(ii) < 0 .or. ib(jj) < 0) then
         if (Index (keywrd, " RESEQ")+Index (keywrd, " SITE=")+Index (keywrd, " 0SCF") == 0) then
           num1 = char(ichar("2") + int(log10(ii*1.01)))
           num2 = char(ichar("2") + int(log10(jj*1.01)))
            write(line,'(a,i'//num1//',a,i'//num2//')') &
            " Impossible user-supplied pi bond between "//elemnt(nat(ii)), ii, &
            " and "//elemnt(nat(jj)), jj
            call mopend(trim(line))
            return
          end if
        end if
      end do
    end if
    if (j == 1) then
      call mopend("AN ERROR IN ASSIGNING PI-BONDS USING ""SETPI"" HAS BEEN DETECTED")
      return
    end if
!
!   Second set: Assign all lone-pairs
!
!
    l = 1
    do i = 1, numat
      if (iz(i) > ib(i)) then
        if (ib(i) > 0) then
          j = Min (ib(i), iz(i)-ib(i))
          do k = 1, j
           call add_Lewis_element(i,0,0, numbon(2))
          end do
        end if
      end if
      if (iz(i) < ib(i)) then
  !
  ! This is a "virtual" lone pair - empty orbital on an atom
  !
        if (ib(i) > 0) then
          j = ib(i) - iz(i)
          do k = 1, j
            call add_Lewis_element(0,i,0, ii) ! ii is not used here.
          end do
        end if
      end if
    end do
!
    do i = 1, numat
      if (ib(i) /= 0) go to 1010
    end do
    l = 0
    goto 1070
1010 continue
    if (debug) then
      write (iw,*) " Core charges after adding Lone Pairs (ATOM CORE ORBS)"
      j = Min (numat, large)
      do k = 1, j, 20
        l = Min (k+19, j)
        write (iw, "(20I4)") (i, i=k, l)
        write (iw, "(20I4)") (iz(i), i=k, l)
        write (iw, "(20I4)") (ib(i), i=k, l)
      end do
    end if
!
!  Third set: pi-bonds attached to five-membered rings.
!
    npi = 0
    not_aromatic_pi = 0
    do ii = 1, numat
      mpii(ii) = 0
      ir5(ii) = 0
      l = 0
      if (ib(ii) /= 0) then
        nbii = nbonds(ii)
        do i = 1, nbii
          jj = ibonds(i, ii)
          if (ib(jj) /= 0) then
            l = l + 1
          end if
        end do
      end if
      mb(ii) = l
    end do
    do i = 1, numat
      if (mb(i) > 2 .and. ir5(i) == 0) then
        call ring5 (i, mb, ir5)
      end if
    end do
    do ii = 1, numat
!
!  If NOT an atom on a five-membered ring, skip loop
!
      if (ir5(ii) /= 0) then
        nbii = nbonds(ii)
        do i = 1, nbii
          jj = ibonds(i, ii)
!
!  If the other end of the pi bond is IN the same five membered
!  ring, skip loop
!
          if (ir5(jj) /= ir5(ii)) then
            if (mb(jj) /= 0 .and. ib(jj) /= 0 .and. ib(ii) > 0) then
              iset = 1
              mb(ii) = mb(ii) - 1
              mb(jj) = mb(jj) - 1
              mpii(ii) = jj
              mpii(jj) = ii
              call add_Lewis_element(ii,jj,0, numbon(3))
              if (debug .and. big) then
                write (iw, "(A,I5,A,I5)") " PI BOND BETWEEN ATOMS", &
               & ii, " AND", jj
              end if
!
!  Store atom numbers of attached atoms that might
!  pi-bond - these have priority.
!
              npi = npi + 1
              ipi(npi) = ii
              npi = npi + 1
              ipi(npi) = jj
            end if
          end if
        end do
      end if
    end do
    do
!
!  Fourth set: pi-bonds (non-ring systems)
!
      iset = 0
      do ii = 1, numat
        l = 0
        if (ib(ii) /= 0) then
          nbii = nbonds(ii)
          do i = 1, nbii
            jj = ibonds(i, ii)
            if (ib(jj) /= 0) then
!
!  Special case: C-S-X treat S as di-valent
!
              if (nat(ii) == 6 .and. nat(jj) == 16 .and. nbonds(jj) > 1 .or. &
                nat(ii) == 16 .and. nbonds(ii) > 1.and. nat(jj) == 6) cycle
              l = l + 1
            end if
          end do
        end if
        mb(ii) = l
      end do
      if (debug) then
        write (iw,*) " Number of atoms pi bonding to each atom (MB)"
        l = Min (numat, large)
        do j = 1, l, 20
          k = Min (j+19, numat)
          write (iw, "(20I4)") (i, i=j, k)
          write (iw, "(20I4)") (mb(i), i=j, k)
        end do
      end if
!
!  If in a graphitic lattice, use rules for poly-conjugated system only.
!
      if ( .not. graphi) then
        do loop = 1, 2
!
!  Do this twice - first time start from carbon end, second time
!  start from any end.  This distinguishes C-C-C-N from N-C-C-C
!
          do i1 = 1, numat
!
!   Identify the end of a ene or polyene system
!
            if (mb(i1) == 1) then
!
!  If heteroatom at end of polyene, start from carbon end.
!
              if (loop /= 1 .or. Nint (tore(nat(i1))-4.d0) == 0) then
                ii = i1
                nbii = nbonds(ii)
                do i = 1, nbii
                  jj = ibonds(i, ii)
                  if (natorb(nat(ii)) == 9 .or. natorb(nat(jj)) == 9  ) cycle
                  if (nat(ii) == 8 .or. nat(jj) == 8) then
                    if (nat(ii) == 6 .or. nat(jj) == 6) then
                      if (distance(ii,jj) > 1.315d0) cycle  !  Exclude C-O, but keep C=O
                    end if
                  end if
                  do while (mb(jj) /=  0 .and. ib(jj) /= 0 .and. ib(ii) > 0)
                    iset = 1
                    not_aromatic_pi = not_aromatic_pi + 1
                    mb(ii) = mb(ii) - 1
                    mb(jj) = mb(jj) - 1
                    mpii(ii) = jj
                    mpii(jj) = ii
                    call add_Lewis_element(ii,jj,0, numbon(3))
                    if (debug .and. big) then
                      write (iw, "(A,I5,A,I5)") " PI BOND BETWEEN ATOMS", ii, " AND", jj
                    end if
                    npi = npi + 1
                    ipi(npi) = ii
                    npi = npi + 1
                    ipi(npi) = jj
!
!  YNE groups take precedence over ENE groups, therefore:
!
                    if (ib(ii) <= 0 .or. ib(jj) <= 0) exit
                  end do
                end do
              end if
            end if
          end do
        end do
      end if
!
!   If any pi bonds are formed, go back and see if any open-ended
!   pi systems have been uncovered.  This detects conjugated ene-yne
!   systems.
!
      if (iset /= 1) then
        do i = 1, numat
          if (ib(i) /= 0) go to 1020
        end do
        go to 1060
1020    if (debug) then
          write (iw,*) " After removing open-ended pi systems"
          write (iw,*) " Number of bonds to each atom (MB, IB)"
          l = Min (numat, large)
          do j = 1, l, 20
            k = Min (j+19, numat)
            write (iw, "(20I4)") (i, i=j, k)
            write (iw, "(20I4)") (mb(i), i=j, k)
            write (iw, "(20I4)") (ib(i), i=j, k)
          end do
          write (iw,*) " Core charges (ATOM CORE ORBS)"
          j = Min (numat, large)
          do k = 1, j, 20
            l = Min (k+19, j)
            write (iw, "(20I4)") (i, i=k, l)
            write (iw, "(20I4)") (iz(i), i=k, l)
            write (iw, "(20I4)") (ib(i), i=k, l)
          end do
        end if
!
!   Remaining Pi system: pi-systems that were attached to open-ended
!   pi-bonds.
!
        do ii = 1, numat
          l = 0
          if (ib(ii) /= 0) then
            nbii = nbonds(ii)
            do i = 1, nbii
              jj = ibonds(i, ii)
              if (ib(jj) /= 0) then
                l = l + 1
              end if
            end do
          end if
          mb(ii) = l
        end do
!
!  If there are no open-ended pi-systems left, then a ring system
!  must be opened.  The first ring to open is a ring that was attached
!  to an open-ended pi-system, and which is bonded to the least number
!  of other pi-atoms.  This should allow graphitic lattices to be
!  calculated
!
        if (first_pi) then
          not_aromatic_pi = 0
          first_pi = .false.
        end if
        l = npi
        do i = 1, npi
          ni_loc = nbonds(ipi(i))
          loop1: do j = 1, ni_loc
            k = ibonds(j, ipi(i))
            if (ib(k) /= 0) then
              do m = 1, l
                if (ipi(m) == k) cycle loop1
              end do
              l = l + 1
              ipi(l) = k
            end if
          end do loop1
        end do
        npi = l
        l = 0
        loop2: do i = 1, npi
          do j = 1, l
            if (ipi(j) == ipi(i)) cycle loop2
          end do
          if (ib(ipi(i)) /= 0) then
            l = l + 1
            ipi(l) = ipi(i)
          end if
        end do loop2
        npi = l
!
!   Now to start building graphite lattice.
!
        loop = 1
        do
          iset = 0
          do i1 = npi, 1, -1
            ii = ipi(i1)
            if (mb(ii) /= 0) then
              nbii = nbonds(ii)
              do i = 1, nbii
                jj = ibonds(i, ii)
                if (natorb(nat(ii)) == 9 .or. natorb(nat(jj)) == 9  ) cycle
                if (mb(jj) /= 0) then
                  if (ib(ii) > 0 .and. ib(jj) > 0) then

!
!   Rules of priority:  LOOP determines the order of
!   priority. LOOP=1 is Rule 1, LOOP=2 is Rule 2, etc.
!   For all bonds, the order of application of the rules
!   is 1, then 2 then 3, etc.
!
!   Rule 1.  If a six membered ring has two Pi bonds
!            already, then add a third Pi bond to make
!            it an aromatic ring.
!   Rule 2.  If a six membered ring has one Pi bond,
!            then add another Pi bond to the same ring,
!            so that Rule 1 can be used.
!   Rule 3.  Any Pi bond that can be made, is made.
!
                    if (loop == 1 .and. &
                         & arom (ii, jj, mpii) &
                         & .or. loop == 2 .and. &
                         & arom2 (ii, jj, mpii) &
                         & .or. loop == 3) then
                      iset = 1
                      mpii(ii) = jj
                      mpii(jj) = ii
                      call add_Lewis_element(ii, jj, 0, numbon(3))
!
!  Check for acetylenic bond
!
                      if (ib(ii) == 1 .and. ib(jj) == 1) then
                        call add_Lewis_element(ii, jj, 0, numbon(3))
                      end if
!
!   Atom JJ was pi-bonding.  Add JJ to list of
!   potential pi bonding atoms, and add on all
!   atoms that might pi bond to II or JJ.
!
                      npi = npi + 1
                      ipi(npi) = jj
                      graphi = .true.
                      l = npi
                      ni_loc = nbonds(ii)
                      do j = 1, ni_loc
                        k = ibonds(j, ii)
                        if (ib(k) /= 0) then
                          l = l + 1
                          ipi(l) = k
                        end if
                      end do
                      ni_loc = nbonds(jj)
                      do j = 1, ni_loc
                        k = ibonds(j, jj)
                        if (ib(k) /= 0) then
                          l = l + 1
                          ipi(l) = k
                        end if
                      end do
                      npi = l
!
! Eliminate duplicates and atoms that have already
! pi bonded
!
                      l = 0
                      loop3: do i2 = 1, npi
                        do j = 1, l
                          if (ipi(j) == ipi(i2)) cycle loop3
                        end do
                        if (ib(ipi(i2)) /= 0) then
                          l = l + 1
                          ipi(l) = ipi(i2)
                        end if
                      end do loop3
                      npi = l
                    end if
                  end if
                end if
              end do
              if (iset == 1 .and. loop /= 1) go to 1030
            end if
          end do
          loop = loop + 1
          if (loop /= 4) then
            cycle
          else
            exit
          end if
1030      loop = 1
        end do
        graphi = .false.
        npi = 0
!
!   Remaining Pi system: rings
!
!   Order of creation of pi bonds:
!   1:  Pi bonds involving atoms bonded to 2 pi bonds and which are
!       adjacent to atoms that form 3 pi bonds (pi bonds in
!       aromatic rings or pi bonds attached to aromatic rings)
!   2:  Pi bonds involving atoms bonded to 2 pi atoms and which are
!       not necessarily adjacent to 3 pi bonds.
!   3:  Any other remaining pi system.
!
        do loop = 1, 3
          do ii = 1, numat
            l = 0
            if (ib(ii) /= 0) then
              nbii = nbonds(ii)
              do i = 1, nbii
                jj = ibonds(i, ii)
                if (ib(jj) /= 0) then
                  l = l + 1
                end if
              end do
            end if
            mb(ii) = l
          end do
          do ii = 1, numat
            if (mb(ii) /= 0 .and. (mb(ii) /= 3 .or. loop >= 3)) then
              nbii = nbonds(ii)
              lnext = .false.
              do i = 1, nbii
                if (mb(ibonds(i, ii)) == 3) then
                  lnext = .true.
                end if
              end do
              lnext = (lnext .or. loop > 2)
!
!  LNEXT is TRUE if atom II is next to an atom bonded to
!  3 other atoms (bonded to a PI system?)
!
              do i = 1, nbii
                jj = ibonds(i, ii)
                if (lnext .and. mb(jj) /= 0 .and. (mb(jj) /= 3 .or. loop > 2)) then
                  if (natorb(nat(ii)) == 9 .or. natorb(nat(jj)) == 9  ) cycle
                  if (ib(ii) > 0 .and. ib(jj) > 0) go to 1040
                end if
              end do
            end if
          end do
        end do
        exit
1040    mpii(ii) = jj
        mpii(jj) = ii
        call add_Lewis_element(ii,jj,0, numbon(3))
        if (debug .and. big) then
          write (iw, "(A,I5,A,I5)") " PI BOND BETWEEN ATOMS", ii, &
               & " AND", jj
        end if
        npi = 1
        ipi(npi) = ii
        npi = npi + 1
        ipi(npi) = jj
      end if
    end do
!
!  Work out the number of aromatic carbon atoms
!
    do i = 1, numat
      if (ib(i) > 0) go to 1050
    end do
    l = 0
    goto 1070
1050 if (debug) then
      write (iw,*) " After removing closed pi systems"
      write (iw,*) " Number of bonds to each atom"
      k = Min (20, numat)
      write (iw, "(20I4)") (i, i=1, k)
      write (iw, "(20I4)") (mb(i), i=1, k)
      write (iw,*) " Core charges (ATOM CORE ORBS)"
      j = Min (numat, large)
      do k = 1, j, 20
        l = Min (k+19, j)
        write (iw, "(20I4)") (i, i=k, l)
        write (iw, "(20I4)") (iz(i), i=k, l)
        write (iw, "(20I4)") (ib(i), i=k, l)
      end do
    end if
    l = 1
    goto 1070
1060 continue
!
!  Work out the number of aromatic carbon atoms
!
    i = 0
    do j = 1, numat
      if (mpii(j) /= 0) i = i + 1
    end do
    l = 0
1070 continue
    deallocate (ipi, ir5, mpii)
1100 continue
    return
end subroutine chklew

subroutine ring5 (i, mb, ir5)
    use common_arrays_C, only : nbonds, ibonds
    use molkst_C, only : numat
    implicit none
    integer, intent (in) :: i
    integer, dimension (numat), intent (in) :: mb
    integer, dimension (numat), intent (inout) :: ir5
!
    integer :: j, j1, k, k1, l, l1, m, m1, n1, ni, nj, nk, nm
    integer :: nrings = 1
    ni = nbonds(i)
    do j1 = 1, ni
      j = ibonds(j1, i)
      if (mb(j) >= 3) then
        nj = nbonds(j)
        do k1 = j1 + 1, ni
          k = ibonds(k1, i)
          if (mb(j) >= 3) then
            nk = nbonds(k)
!
!   Atoms J and K are attached to atom I.
!   Find atoms, other than I, that are attached to J and K.
!
            do l1 = 1, nj
              l = ibonds(l1, j)
              if (l /= i) then
                if (mb(l) >= 3) then
                  do m1 = 1, nk
                    m = ibonds(m1, k)
                    if (m /= i) then
                      if (mb(m) >= 3) then
!
!  Atoms L and M are attached to J and K respectively.
!  The connectivity is thus: M-K-I-J-L.  Now, check:
! is M attached to L?
!
                        nm = nbonds(m)
                        do n1 = 1, nm
                          if (ibonds(n1, m) == l) go to 1000
                        end do
                      end if
                    end if
                  end do
                end if
              end if
            end do
          end if
        end do
      end if
    end do
    return
!
! Flag atoms I, J, K, L, and M as being members of a five-membered ring.
!
1000 nrings = nrings + 1
    ir5(i) = nrings
    ir5(j) = nrings
    ir5(k) = nrings
    ir5(l) = nrings
    ir5(m) = nrings
end subroutine ring5
logical function arom (ii, jj, mpii)
    use common_arrays_C, only : nbonds, ibonds
    use molkst_C, only : numat
    implicit none
    integer, intent (in) :: ii, jj
    integer, dimension (numat), intent (in) :: mpii
!
    integer :: i, iia, iib, j, j2, jja, jjb, jjc, ni, nj, njb
!***********************************************************************
!
!  AROM determines if the PI bond between II and JJ is part of an
!       aromatic ring.  It does this by checking to see if two other pi
!       bonds already exist in the ring.  If they do, then II and JJ
!       form the third and final pi bond.
!
!***********************************************************************
    ni = nbonds(ii)
    nj = nbonds(jj)
    do i = 1, ni
      iia = ibonds(i, ii)
      if (iia /= jj) then
        if (mpii(iia) /= 0) then
!
!  IIA is attached to II and forms a pi bond
!  IIB is a pi-bonded atom attached to IIA
!
          iib = mpii(iia)
          do j = 1, nj
            jja = ibonds(j, jj)
            if (jja /= ii) then
              if (mpii(jja) /= 0) then
!
!  JJA is attached to JJ and forms a pi bond
!  JJB is a pi-bonded atom attached to JJA
!
                jjb = mpii(jja)
                njb = nbonds(jjb)
                do j2 = 1, njb
                  jjc = ibonds(j2, jjb)
!
! Is any atom attached to JJB the same as IIB?
!
                  if (jjc == iib) then
                    arom = .true.
                    return
                  end if
                end do
              end if
            end if
          end do
        end if
      end if
    end do
    arom = .false.
end function arom
logical function arom2 (ii, jj, mpii)
    use common_arrays_C, only : nbonds, ibonds
    use molkst_C, only : numat
    implicit none
    integer, intent (in) :: ii, jj
    integer, dimension (numat), intent (in) :: mpii
    integer :: i, i1, iia, iib, iic, j, j1, jja, jjb, jjc, ni, nia, nib, nj, &
   & nja, njb
!***********************************************************************
!
!  AROM2 determines if the PI bond between II and JJ is part of an
!        aromatic ring.  It does this by detecting a second PI bond
!        in the ring, and then identifying the final PI bond.
!        The atom numbers of the final PI bond are put into II and JJ
!
!***********************************************************************
    ni = nbonds(ii)
    nj = nbonds(jj)
!
!  Was II attached to a PI bonded system?
!
    do i = 1, ni
      iia = ibonds(i, ii)
      if (iia /= jj) then
        nia = nbonds(iia)
        if (mpii(iia) /= 0) then
          iib = mpii(iia)
!
!   IIA is attached to to II
!   IIB is pi bonded to IIA
!
          nib = nbonds(iib)
          do j = 1, nj
            jja = ibonds(j, jj)
            if (jja /= ii) then
              nja = nbonds(jja)
!
!   JJA is attached to JJ
!
              do i1 = 1, nib
                iic = ibonds(i1, iib)
                do j1 = 1, nja
!
!   Is any atom attached to IIB also attached to JJA?
!
                  if (iic == ibonds(j1, jja)) then
                    arom2 = .true.
                    return
                  end if
                end do
              end do
            end if
          end do
        end if
      end if
    end do
!
!  Was JJ attached to a PI bonded system?
!
    do j = 1, nj
      jja = ibonds(j, jj)
      if (jja /= ii) then
        nja = nbonds(jja)
        if (mpii(jja) /= 0) then
          jjb = mpii(jja)
!
!   JJA is attached to to JJ
!   JJB is pi bonded to JJA
!
          njb = nbonds(jjb)
          do i = 1, ni
            iia = ibonds(i, ii)
            if (iia /= jj) then
              nia = nbonds(iia)
!
!  IIA is attached to II
!
              do j1 = 1, njb
                jjc = ibonds(j1, jjb)
                do i1 = 1, nia
!
!   Is any atom attached to JJB also attached to IIA?
!
                  if (jjc == ibonds(i1, iia)) then
                    arom2 = .true.
                    return
                  end if
                end do
              end do
            end if
          end do
        end if
      end if
    end do
    arom2 = .false.
end function arom2
subroutine add_Lewis_element(atom_i, atom_j, charge, element_type)
!
!  All the elements of the Lewis structure are stored in Lewis_elem
!  These are used in constructing the starting LMOs
!
!  atom_i : If atom_j is non-zero, atom_i is the atom number of one of the two atoms in a bond.
!  atom_i : If atom_j is zero, atom_i is the atom number of an occupied (normal) lone-pair
!           If atom_i is negative, -atom_i is the atom number of either an anion or a cation,
!           the choice is made later.
!  atom_j : If atom_i is non-zero, atom_j is the atom number of the other of the two atoms in a bond.
!  atom_j : If atom_i is zero, atom_j is the atom number of a virtual lone-pair
!  charge : The charge implied by this Lewis element (If Na virtual lone pair, then +1)
!  Element_type : Type of Lewis structural element (sigma, lone pair, pi bond, etc.)
!
  use MOZYME_C, only: ions, Lewis_elem, Lewis_tot, Lewis_max, iz, ib
  integer, intent (in) :: atom_i, atom_j, charge
  integer, intent (inout) :: element_type
  integer, dimension(:,:), allocatable :: Lewis_save
    Lewis_tot = Lewis_tot + 1
! Expand Lewis_elem buffer if maximum size is exceeded
    if (Lewis_tot > Lewis_max) then
      Lewis_max = 2*Lewis_max
      allocate(Lewis_save(2,Lewis_max))
      Lewis_save(:,:Lewis_max) = Lewis_elem(:,:Lewis_max)
      deallocate(Lewis_elem)
      allocate(Lewis_elem(2,2*Lewis_max))
      Lewis_elem(:,:Lewis_max) = Lewis_save(:,:Lewis_max)
      deallocate(Lewis_save)
      Lewis_max = 2*Lewis_max
    endif
    Lewis_elem(1,Lewis_tot) = atom_i
    Lewis_elem(2,Lewis_tot) = atom_j
    if (atom_i > 0 .and. atom_j > 0) then
!
! Diatomic bond (a line in the Lewis structure)
! Decrement both atoms by one electron and one atomic orbital
!
      iz(atom_i) = iz(atom_i) - 1
      iz(atom_j) = iz(atom_j) - 1
      ib(atom_i) = ib(atom_i) - 1
      ib(atom_j) = ib(atom_j) - 1
      element_type = element_type + 1
    else if (atom_i /= 0) then
      if (atom_i > 0) then
        if (charge == -1) then
!
! A lone pair (two dots in the Lewis structure) arising from a unit negative charge
! Decrement the atom by one electron and one atomic orbital.  The other electron
! comes from the negative charge
!
          iz(atom_i) = iz(atom_i) - 1
        else if (charge == 0) then
!
! A normal lone pair (two dots in the Lewis structure)
! Decrement the atom by two electrons and one atomic orbital.
!
          iz(atom_i) = iz(atom_i) - 2
        end if
        ib(atom_i) = ib(atom_i) - 1
        element_type = element_type + 1  ! Lone pair
      else
!
! At this point, the type is unknown, so simply decrement ib,
! and don't change iz.
!
        ib(-atom_i) = ib(-atom_i) - 1
      end if
    else
      if (charge == 2) then
!
! A virtual lone pair arising from a di-cation
! Decrement the atom by two electrons and one atomic orbital.
! (Think of Ca(++) - this is calcium after removal of two electrons
!
        iz(atom_j) = iz(atom_j) - 2
      else if (charge == 1) then
!
! A virtual lone pair arising from a cation
! Decrement the atom by one electrons and one atomic orbital.
! (Think of Na(+) - this is sodium after removal of one electron
!
        iz(atom_j) = iz(atom_j) - 1
      end if
      ib(atom_j) = ib(atom_j) - 1
    end if
    if (charge == 0) return
    ions(atom_i + atom_j) = ions(atom_i + atom_j) + charge
end subroutine add_Lewis_element
