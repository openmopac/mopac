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

subroutine geochk ()
!***********************************************************************
!                                                                      *
!  GEOCHK DOES SEVERAL VERY DIFFERENT THINGS:
!
!   (A) IT CHECKS THE GEOMETRY TO MAKE SURE IT CONFORMS TO THE LEWIS
!       STRUCTURE CONVENTIONS.
!
!   (B) IT IDENTIFIES ALL IONIZED ATOMS. LATER ON, MAKVEC WILL NEED TO
!       KNOW WHICH ATOMS ARE IONIZED.
!
!   (C) IT CALCULATES THE CHARGE ON THE SYSTEM.
!
!   (D) (OPTIONAL) IT RESEQUENCES THE ATOMS SO THAT ALL ATOMS OF EACH
!       RESIDUE ARE CONTIGUOUS.  IF THIS IS DONE, THE JOB CANNOT BE
!       CONTINUED, SO THE NEW GEOMETRY IS PRINTED OUT AND THE RUN
!       STOPPED.
!
!                                                                      *
!***********************************************************************
    use common_arrays_C, only : geo, coord, nfirst, nlast, breaks, &
       & labels, nat, na, nb, nc, nbonds, ibonds, txtatm, atmass, pdiag, &
         txtatm1, chains, all_comments, coorda, loc, lopt, l_atom, break_coords
    use MOZYME_C, only : ions, icharges, angles, allres, iz, ib, tyr, allr, &
    iopt, nres, at_res, Lewis_tot, Lewis_elem, noccupied, nvirtual, &
    odd_h, tyres, start_res, lstart_res, uni_res
!    
    use molkst_C, only: natoms, numat, nelecs, keywrd, moperr, maxtxt, mozyme, &
      maxatoms, line, nalpha, nbeta, uhf, nclose, nopen, norbs, numcal, id, &
      ncomments, numat_old, nvar, prt_coords, prt_topo, allkey, txtmax, pdb_label
    use chanel_C, only: iw, iarc, archive_fn, log, ilog
    use atomradii_C, only: atom_radius_covalent
    use parameters_C, only : ams, natorb, tore, main_group
    use elemts_C, only : elemnt
    implicit none
!
    integer, parameter :: max_sites = 400
    character :: padding*40, txtatm_1*27, txtatm_2*27, tmp*130
    character, allocatable :: Lewis_formatted(:,:)*20, temp_txtatm(:)*27
    character (len=1), dimension (:), allocatable :: atom_charge
    character :: ion_names(-6:6)*12, charge(max_sites,3)*1, num
    double precision, dimension(:), allocatable :: radius
    logical, save :: debug, let, lres, lreseq, times, opend, charges, l_protein, &
      done, neutral(100), lsite, ter, residues, lbreaks, l_use_old_labels, &
    l_names(max_sites), first, first_prt, header, l_rama, l_salt
    logical, dimension (:), allocatable :: ioptl
    integer :: i, ibad, ichrge, irefq, ires, ii, jj, m, nfrag, io, kk, kkk, near_ions_store(2), &
   & j, jbad, k, l, large, n1, new, alloc_stat, mres, near_ions(2, 100), atomic_charges(8), &
     maxtxt_store, nn1, n_new, new_res(max_sites), j2, mbreaks, icalcn = -10, nnumat, delta_res, max_frag
    integer, dimension(:), allocatable ::  mb    
    integer, save :: numbon(3), num_ions(-6:6)
    integer, dimension(:), allocatable :: nnbonds
    integer, dimension(:,:), allocatable :: iibonds
    double precision :: sum, r_ions(100)
    character :: new_name(max_sites)*3, old_name(max_sites)*3, new_chain(max_sites)*1, &
      new_alt(max_sites)*1
    double precision, external :: reada
    data numbon / 3 * 0 /
    data atomic_charges /1, 4, 3, 2, -1, -2, -3, -4/
    data ion_names / &
      "Less than -5", &
      "Penta-anion ", &
      "Tetra-anion ", &
      "Tri-anion   ", &
      "Di-anion    ", &
      "Anion       ", &
      "(not used)  ", &
      "Cation      ", &
      "Di-cation   ", &
      "Tri-cation  ", &
      "Tetra-cation", &
      "Penta-cation", &
      "More than +5" /
!
    intrinsic Abs, Index, Min, Nint, Allocated
    if (icalcn /= numcal) then
      first_prt = .true.
      icalcn = numcal
    end if
    l_rama = (index(keywrd, " writemo_rama") /= 0)
!
    if (allocated(ions))    deallocate (ions)
    if (allocated(iopt))    deallocate(iopt)
    if (allocated(at_res))  deallocate(at_res)
    if (Allocated (iz))     deallocate (iz)
    if (Allocated (ib))     deallocate (ib)
    allocate (ions(maxatoms), iz(maxatoms), ib(maxatoms), mb(maxatoms), at_res(maxatoms), &
            & atom_charge(maxatoms), ioptl(maxatoms), iopt(maxatoms), radius(maxatoms), &
            stat=alloc_stat)
    if (alloc_stat /= 0) then
      call mopend("Failed to allocate arrays in GEOCHK")
      go to 1100
    end if
    l_protein = .false.
    ib(:) = 0
    at_res = 0
    mres = 0
    ions = 0
    j = 0
    lbreaks = (breaks(1) /= -300)
    if (.not. lbreaks) mbreaks = 0
    do i = 1, ncomments
      if (index(all_comments(i), "REMARK   2") == 0) then
        j = j + 1
        all_comments(j) = all_comments(i)
      end if
    end do
    ncomments = j
!
    done = .false.
    n_new = 0
    line = trim(keywrd)
    i = Index (line, " XENO")
    if (i /= 0) then
      do
        if (line(i:i) == "(") exit
        line = line(2:)
      end do
      j = index(line(i:), ") ") 
      if (j == 0) then
        call mopend("Closing parenthesis for XENO not found")
        return
      end if      
      line = ","//line(i + 1: j + i)
      do
        if (line == " ") exit
        do j = 1, 100
          if (line(1:1) == "," .or. line(1:1) == ";" .or. line(1:1) == ")") exit
          line = line(2:)
        end do
        if (line(1:1) == ")") exit
        line = line(2:)
        k = index(line, "=")
        if (k == 0) then
          call mopend("Equals sign (""="") expected in XENO but not found")
          return
        end if  
        l = k + 2
        j = k + 4
        n_new = n_new + 1
        new_chain(n_new) = line(1:1)
        if (new_chain(n_new) >= "0" .and. new_chain(n_new) <= "9") new_chain(n_new) = "A"
        new_res(n_new) = nint(reada(line, 1))
        if (line(l:l) == "," .or. line(l:l) == ";" .or. line(l:l) == ")") then
!
!   A one-letter residue name has been detected.
!
          new_name(n_new) = line(k + 1:k + 1)
        else if (line(j:j) == "," .or. line(j:j) == ";" .or. line(j:j) == ")") then
!
!   A three-letter residue name has been detected.
!
          new_name(n_new) = line(k + 1:k + 3)
        else  
          write(iw,'(/10x,a)')"XENO keyword used: ("//trim(line)
          call mopend("RESIDUE NAMES IN KEYWORD XENO MUST BE ONE OR THREE CHARACTERS LONG")
          return
        end if
      end do
    end if      
    do
      i = Index (keywrd, " xeno")
      if (i == 0) exit
      keywrd(i:i+6) = " XENO=("
    end do
    do i = 1, n_new
      do j = i + 1, n_new
        if (new_chain(i) == new_chain(j) .and. new_res(i) == new_res(j) .and. &
          new_name(i) /= new_name(j)) then
          write(line,'(a,a,i4,a)')" XENO residue ",new_chain(i), new_res(i),&
          " occurs more than once.  Remove extra definition."
          call mopend(trim(line))
          return
        end if
      end do
      if (new_name(i)(2:3) == "  ") then
        do j = 1, 20
          if (tyr(j) == new_name(i)(1:1)) then
            new_name(i) = tyres(j)
            exit
          end if
        end do
      end if
    end do
!
!   Add or remove hydrogen atoms, as necessary.
!
    i =  index(keywrd," SITE=(IONIZE)")
    if (i > 0) then
      line = "SITE=(COO,NH3,ARG(+),SO4,PO4)"
      keywrd = keywrd(:i)//trim(line)//keywrd(i + 14:)
    end if      
    i = index(keywrd," SITE=")
    lsite = .false.
!
!   Identify sites that either should be ionized or should not be ionized.
!
!   In array "neutral", these are:
!
!   1:   -COOH
!   2:   -COO(-)
!   3:   -NH3(+)
!   4:   -NH2
!   5:   -Arg(+)-
!   6:   -Arg-
!   7:   -His(+)-
!   8:   -His-
!   9:    SO4(=)
!  10:    PO4(=)
!
    
    i = index(keywrd," SITE=(")
    if (i /= 0 .and. index (keywrd, " ADD-H") == 0) then
!
!  Quick check that obvious structures are present
!
      j = index(keywrd(i:), ") ") + i
      line = keywrd(i + 7:j)
      do 
        i = index(line,'"')
        if (i == 0) exit
        j = index(line(i + 1:), '"') + i + 1
        if (line(j:j) /= "(") then
         call mopend("ERROR IN SITE KEYWORD AT THE END OF '"//line(i:j)//"'")
          write(iw,'(10x,a)')"The last three characters must be one of '(+)', '(0)', or '(-)', e.g,"//line(i:j-1)//"(+)"
          return
        end if
        j = j + 1
        if (line(j:j) /= "+" .and. line(j:j) /= "0" .and. line(j:j) /= "-") then
          call mopend("ERROR IN SITE KEYWORD AT THE END OF '"//line(i:j)//"'")
          write(iw,'(10x,a)')"The last character must be a '+', a '0', or a '-'"
          return
        end if
        j = j + 1
        if (line(j:j) == "+" .or. line(j:j) == "0" .or. line(j:j) == "-") j = j + 1
        if (line(j:j) /= ")" .and. line(j:j) /= "2" .and. line(j:j) /= "3" ) then
          call mopend("ERROR IN SITE KEYWORD AT THE END OF '"//line(i:j)//"'")
          write(iw,'(10x,a)')"The last character must be a '2', a '3', or a ')'"
          return
        end if
        j = j + 1
        if (line(j:j) /= "," .and. line(j:j) /= ")") then
          call mopend("ERROR IN SITE KEYWORD AT THE END OF '"//line(i:j)//"'")
          write(iw,'(10x,a)')"The last character must be a comma ',' or a close parenthesis ')'"
          return
        end if
        line = line(j:)
      end do  
      call lewis(.true.)
      do
        i = index(keywrd," SITE=")
        if (keywrd(i + 7:i + 7) == ")") exit
        if (i /= 0 .and. index (keywrd, " ADD-H") == 0) then
          if (moperr) return
          j = index(keywrd(i:), ") ") + i
          call upcase(keywrd(i:j), j - i + 1)
          allkey = keywrd(i:j)
          neutral = .false.
          if (index(keywrd(i:j),"""") == 0) then
            if (index(keywrd(i:j),"COOH") /= 0) then
              neutral(1) = .true.
            else if (index(keywrd(i:j),"COO") /= 0) then
              neutral(2) = .true.
            end if
            if (index(keywrd(i:j),"NH3") /= 0) then
              neutral(3) = .true.
            else if (index(keywrd(i:j),"NH2") /= 0) then
              neutral(4) = .true.
            end if
            if (index(keywrd(i:j),"ARG(+)") /= 0) then
              neutral(5) = .true.
            else if (index(keywrd(i:j),"ARG") /= 0) then
              neutral(6) = .true.
            end if
            if (index(keywrd(i:j),"HIS(+)") /= 0) then
              neutral(7) = .true.
            else if (index(keywrd(i:j),"HIS") /= 0) then
              neutral(8) = .true.
            end if
            if (index(keywrd(i:j),"SO4") /= 0) then
               neutral(9) = .true.
            end if
            if (index(keywrd(i:j),"PO4") /= 0) then
               neutral(10) = .true.
            end if        
          end if
          do k = 1, 10
            if (neutral(k)) exit
          end do
          if (k > 10) then
            if (index(keywrd(i:j), "SALT") /= 0) then
              call find_salt_bridges(numbon, numbon, 0, 0)
              if (moperr) return
              if (index(keywrd(i:j), "SALT") /= 0) then
                i = index(keywrd(i:j), "SALT") + i - 2
                do j = i + 1, len_trim(keywrd)
                  if (keywrd(j:j + 1) == ") " .or. keywrd(j:j + 1) == ",") exit
                end do
                keywrd = keywrd(:i)//trim(keywrd(j:))
                cycle
              end if
            end if
    !
    !  Now check for specific residues
    !
            i = index(keywrd," SITE=")
            j = index(keywrd(i:), ") ") + i
            do i = i, len_trim(keywrd)
              if (keywrd(i:i) == "(") exit
            end do
            j = j - 1
            m = 0
            charge = "*"
            do
              i = i + 1
              do
                if (keywrd(i:i) /= """") exit
                do i = i + 1, len_trim(keywrd)
                  if (keywrd(i:i) == ")") exit
                end do
                i = i + 2
              end do
              if (i >= j) exit
              m = m + 1
              new_chain(m) = keywrd(i:i)
              num = char(Int(log10(m*1.0)) + ichar("1") + 1) 
              if (new_chain(m) < "A" .or. new_chain(m) > "Z") then             
                write(iw,'(//10x,a,/10x,a,i'//num//',a)')" There is a fault in the SITE keyword", &
                " Chain letter for residue", m, " is not in the range 'A' to 'Z'"
              end if          
              i = i + 1
              k = ichar(keywrd(i:i)) - ichar("0")
              if ((k < 0 .or. k > 9) .and. keywrd(i:i) /= "-") then
                write(iw,'(//10x,a,/10x,a,i'//num//',a)')" There is a fault in the SITE keyword", &
                " The entry for site", m, " is faulty"
              end if  
              if (keywrd(i:i) == "-") then
                ii = -1
                k = 0
              else
                ii = 1
              end if            
              do
              i = i + 1
              if (keywrd(i:i) < "0" .or. keywrd(i:i) > "9") exit
                k = k*10 + ichar(keywrd(i:i)) - ichar("0")  
              end do
              new_res(m) = k*ii
              if (txtmax == 27 .and. keywrd(i:i) /= "(") then
                new_alt(m) = keywrd(i:i)
                if ((new_alt(m) < "A" .or. new_alt(m) > "Z") .and. new_alt(m) /= " ") &             
                  write(iw,'(//10x,a,/10x,a,i'//num//',a)')" There is a fault in the SITE keyword", &
                  " Alternate residue for residue", m, " is not in the range 'A' to 'Z'"
                i = i + 1
              else
                new_alt(m) = " "
              end if                   
              i = i + 1
              charge(m,1) = keywrd(i:i)
              if (charge(m,1) /= "+" .and. charge(m,1) /= "-" .and. charge(m,1) /= "0") then
                write(iw,'(//10x,a,/10x,a,i4,a)')" There is a fault in the SITE keyword", &
                " Charge on site", m, " is not '+', '0', or '-'"   
              end if
              i = i + 1
              charge(m,2) = keywrd(i:i)
              if (charge(m,2) == "+" .or. charge(m,2) == "-" .or. charge(m,2) == "0") then
                i = i + 1
              else
                charge(m,2) = "*"
              end if
              i = i + 1
              continue
              if (keywrd(i:i) == ")") i = i + 1
            end do
          else
            m = 0
          end if      
          lsite = .true.
          call update_txtatm(.true., .true.)
          do j = 1, id
            nat(numat + j) = 107
          end do
          i = numat
          if (index(keywrd, " NOSITE") == 0 .or. index(keywrd, " CVB") == 0) then
            call store_and_restore_Tv("STORE")
            call site(neutral, new_chain, new_res, new_alt, charge, m, max_sites, allkey)
            call store_and_restore_Tv("RESTORE")
          end if
          if (moperr) return
          l_atom(i:natoms) = .true.
          call update_txtatm(.true., .true.)
          if (index(keywrd, " LET") /= 0) moperr = .false.
          if (moperr) return
        end if
        i = index(keywrd," SITE=")
        if (i == 0) exit
        do j = i + 1, len_trim(keywrd)
          if (keywrd(j:j + 1) == ") ") exit
        end do
        keywrd = keywrd(:i)//trim(keywrd(j + 1:))
      end do
    end if
!
!  Store charge, if present
!
    do i = 1, natoms
      atom_charge(i) = txtatm(i)(2:2)
!
!  Prevent atom number being mis-read as a charge.
!
      if (txtatm(i)(2:2) /= "+" .and. txtatm(i)(2:2) /= "-" .and. txtatm(i)(2:2) /= "0") &
      & atom_charge(i) = " "
    end do
!
!    ASSIGN LOGICALS USING KEYWRD
!
    lres = (Index (keywrd, " RESI") + Index (keywrd, " RESEQ") /= 0)
    if (.not. lres) lres = (Index (keywrd, " PDBOUT") /= 0 .and. maxtxt /= txtmax)
    if (.not. lres) lres = (Index (keywrd, " ADD-H") /= 0 .and. Index (keywrd, " NORESEQ") == 0)
    lreseq = (Index (keywrd, " NORESEQ") == 0 .and. Index (keywrd, " RESEQ") /= 0)
    let = (Index (keywrd, " 0SCF")+Index (keywrd, " LET")+ &
   & Index (keywrd, " RESEQ")+Index (keywrd, " GEO-OK") /= 0)
    times = (Index (keywrd, " TIMES") /= 0)
    if (times) then
      call timer (" START OF GEOCHK")
    end if
    debug = (Index (keywrd, " GEOCHK") /= 0)
    if (Index (keywrd, " CHARGE=") /= 0) then
      irefq = Nint (reada (keywrd, Index (keywrd, " CHARGE=")))
    else
      irefq = 0
    end if
    call extvdw_for_MOZYME (radius, atom_radius_covalent)
    if (moperr) return
!
    if (Index (keywrd, " LARGE") /= 0) then
      large = 1000000
    else
      large = 20
    end if
    nnumat = numat
!
!   WORK OUT WHAT ATOMS ARE BONDED TO EACH OTHER.
!
    call lewis (.true.)
    if (moperr) return   
    allocate (nnbonds(numat), iibonds(15,numat), stat=alloc_stat)
    if (alloc_stat /= 0) then
      call memory_error ("geochn")
      go to 1100
    end if
!
!  Store nbonds and ibonds in case they are modified within this subroutine
!
    nnbonds = nbonds
    iibonds = ibonds
    if (moperr) then
      if (Index (keywrd, " GEO-OK") == 0) then
        write (iw,*) " GEOMETRY CONTAINS FAULTS. TO CONTINUE CALCULATION SPECIFY ""GEO-OK"""
        go to 1100
      else
        moperr = .false.
      end if
    end if
!
!   ZERO OUT IONS
!
!
!  FIND THE NITROGEN ATOM OF THE N END OF THE PROTEIN.
!
    do i = 1, numat
      ioptl(i) = .false.
    end do
    call findn1 (n1, ioptl, io, delta_res)
    if (lreseq) then
!
!   RESEQUENCE THE ATOMS.  WHEN RESEQ IS CALLED, THE SEQUENCE OF
!   ATOMS IS CHANGED.  A CONSEQUENCE OF THIS IS THAT THE SCF CANNOT
!   BE RUN.
!
! First, delete all bonds between ATOMs and HETATMs
!
      do i = 1, numat
        if (txtatm(i)(:4) /= "ATOM") cycle
        do j = 1, nbonds(i)
          k = ibonds(j,i)
          if (txtatm(k)(:4) /= "ATOM") ibonds(j,i) = 0
        end do
        l = 0
        do j = 1, nbonds(i)
          k = ibonds(j,i)
          if (k > 0) then
            l = l + 1
            ibonds(l,i) = k
          end if
        end do
        nbonds(i) = l
      end do
      if (index(keywrd, "RESID") /= 0) txtatm1(:numat) = " "
      new = 0
      iz = -1000
      do
        if (n1 /= 0) call reseq (ioptl, iz, n1, new, io)
          if (n1 == -100) goto 90
        if (moperr) go to 1100
        call findn1 (n1, ioptl, io, delta_res)
        if (n1 == 0) exit
      end do
      if (new /= numat) then
        do i = 1, numat
          if (nat(i) /= 1 .and. .not. ioptl(i)) then
!
!   IDENTIFY ALL NON-PROTEIN MOLECULES IN THE SYSTEM
!
            call moiety (ioptl, iz, i, new)
          end if
        end do
        do i = 1, numat
          if (nat(i) == 1 .and. .not. ioptl(i)) then
!
!   IDENTIFY ALL HYDROGENS ATTACHED TO RESIDUE-LIKE SPECIES
!
            new = new + 1
            iz(new) = i
          end if
        end do
        if (new /= numat) then
          write (iw,*) " THERE IS A FAULT IN RESEQ"
          write (iw,'(a,i5)') "  Number of atoms found in data-set:  ", numat
          write (iw,'(a,i5)') "  Number of atoms after re-sequencing:", new
          if (new < numat) then
            write(iw,'(a)') "Atoms missing (Use original numbering system)"
            ioptl(:numat) = .true.
            do i = 1, new
              ioptl(iz(i)) = .false.              
            end do
            do i = 1, numat
              if (ioptl(i)) write(iw,'(i5)')i
            end do         
          else
          end if
          call mopend("THERE IS A FAULT IN RESEQ")
          go to 1100
        end if
!
!  Unconditionally, convert geometry into Cartesian coordinates
!
        geo(:,:numat) = coord(:,:numat)
        do i = 1, numat
          if (na(i) /= 0) exit
        end do
        if (i <= numat) then
          call mopend("SOME COORDINATES WERE IN INTERNAL. THESE HAVE BEEN CHANGED TO CARTESIAN")
          moperr = .false.
        end if
        na = 0
      end if
      l = 1
      do i = 1, numat
        nfirst(i) = l
        j = iz(i)
        mb(j) = i
        do k = 1, 3
          geo(k, i) = coord(k, j)
        end do
        labels(i) = nat(j)
        nlast(i) = nfirst(i) + natorb(labels(i)) - 1
        l = nlast(i) + 1 
      end do
      nat(:numat) = labels(:numat)
      coord(:,:numat) = geo(:,:numat)
      done = .true.
!
!   DUMMY ATOMS ARE EXCLUDED, THEREFORE
!
      natoms = numat
!
!
!   REARRANGE ATOMS TO SUIT THE NEW NUMBERING SYSTEM
!
      do i = 1, numat
        j = iz(i)
        ib(i) = Min (nbonds(j), 4)
        l = Min (nbonds(i), 4)
        do k = 1, l
          ibonds(k+4, i) = mb(ibonds(k, i))
        end do
      end do
      do i = 1, numat
        nbonds(i) = ib(i)
        do k = 1, nbonds(i)
          ibonds(k, i) = ibonds(k+4, iz(i))
        end do
      end do
    end if
    noccupied = 0  
    if (Index (keywrd, " RESEQ") + Index (keywrd, " SITE=") + index(keywrd, " ADD-H") &
      + index(keywrd, " 0SCF") == 0) then
!
!  EXAMINE THE GEOMETRY - IDENTIFY THE LEWIS ELEMENTS (SIGMA BONDS,
!  LONE PAIRS, CATIONS, PI BONDS, ANIONS, OTHER CHARGES)
!
      numbon = 0
      call chklew (mb, numbon, l, large, debug)
      if (moperr) return
      l = 0
      do i = 1, numat
        l = l + Abs (iz(i))
      end do
      if (l /= 0) then         
!
!  THERE ARE IONS.  IDENTIFY THEM.
!
        call chkion (mb, numbon(2), atom_charge)
        if (lreseq) moperr = .false.
        if (moperr) return
      end if
      do i = 1, numat
        ions(i) = nint(tore(nat(i)))
      end do
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) then
          noccupied = noccupied + 1
          j = Lewis_elem(1,i)
          if (Lewis_elem(2,i) > 0) then
            k = Lewis_elem(2,i)
            ions(k) = ions(k) - 1! one electron from a bond
            ions(j) = ions(j) - 1! one electron from a bond
          else
            ions(j) = ions(j) - 2! two electrons from a lone pair
          end if
        end if
      end do
      nvirtual = 0
      do i = 1, Lewis_tot
        if (Lewis_elem(2,i) > 0) nvirtual = nvirtual + 1
      end do
    else
      nvar = 0
      do i = 1, numat
        do j = 1, numat_old          
          if (abs(coord(1,i) - coorda(1,j)) > 0.1d0) cycle
          if (abs(coord(2,i) - coorda(2,j)) > 0.1d0) cycle
          if (abs(coord(3,i) - coorda(3,j)) > 0.1d0) cycle
          exit
        end do 
        if (j <= numat_old) then
!
!  Atom "i" in coord has the same coordinates as atom "j" in coorda
!
          do k = 1, 3
            if (lopt(k,j) == 1) then
              nvar = nvar + 1
              loc(1,nvar) = i
              loc(2,nvar) = k
            end if
          end do
        else
          do k = 1, 3
            nvar = nvar + 1
            loc(1,nvar) = i
            loc(2,nvar) = k
          end do          
        end if
      end do
      do j = numat + 1, numat + id
        do k = 1,3
          if (lopt(k,j) == 1) then
            nvar = nvar + 1
            loc(1,nvar) = j
            loc(2,nvar) = k
          end if
        end do
      end do       
    end if 
!
!   Use the following block to debug the construction of the Lewis structure
!
    if ( .false.) then
!
!  Sanity check - are all atomic orbitals accounted for?
!
      k = abs(norbs - noccupied - nvirtual)
      do i = 1, numat
        if (iz(i) /= 0) then
          k = k + 1
        end if
        if (ib(i) /= 0) then
          k = k + 1
        end if
      end do
      iz = 0
      do i = 1, Lewis_tot
        j = Lewis_elem(1,i)
        if (j > 0) iz(j) = iz(j) + 1
        j = Lewis_elem(2,i)
        if (j > 0) iz(j) = iz(j) + 1
      end do
      do i = 1, numat
        if (iz(i) - natorb(nat(i)) /= 0) then
          k = k + 1
        end if
      end do
      if (k /= 0) then
        write(iw,*)" An error has been detected."
      end if
    end if
!
!
    if (lres) then      
      txtatm(:) = " "
      if (maxtxt == 0) then
        maxtxt = 26
        txtmax = 26
      end if
      angles = 0.d0
      allres = " "
      if (done) then        
!
!   THE EARLIER CALL TO RESEQ MEANS THAT N1 MIGHT HAVE MOVED.
!   SO FIND N1 AGAIN.
!
        ioptl(:numat) = .false.
        call findn1 (n1, ioptl, io, delta_res)
      end if
      l_protein = (n1 /= 0)
      ib(:numat) = -100000
      nfrag = 0
      ires = 0
      uni_res = 0
      odd_h = .true.
!
!  Break all intra-chain bonds, so that the residues can easily be
!  identified.
!
      call lyse !
      allr = " "
      do max_frag = 1, 10000
        if (start_res(max_frag) == -200) exit
      end do
      max_frag = max_frag - 1  
      lstart_res = .false.
      nfrag = 0
      do      
        if (max_frag > 0) then
          do nfrag = 1, max_frag
            if (.not. lstart_res(start_res(nfrag) + 1)) exit
          end do
        else
          nfrag = nfrag + 1
        end if
        if (max_frag > 0) then
          if (nfrag > max_frag) exit
          ires = start_res(nfrag) - delta_res
        end if
!
! General bug-trap: if the number of fragments is unreasonably large,
! assume the system is unrecognizable and exit
!
        if (nfrag > 999) then
          call mopend("STRUCTURE UNRECOGNIZABLE")
          inquire(unit=iarc, opened=opend) 
          if (opend) close(iarc, status = "DELETE")          
          go to 1100
        endif
!
        call names (ioptl, ib, n1, ires, nfrag, io, uni_res, mres)
        if (moperr) return
!    
        nn1 = n1
        call findn1 (n1, ioptl, io, delta_res)
!
! n1 is the start of the next fragment
! delta_res is the distance back along the chain to the start of the next fragment.
! 
        if (.not. l_protein) nfrag = 0
        if (n1 == 0) exit  
        if (n1 == nn1) ioptl(n1) = .true.   
        if (.not. lbreaks) then
          mbreaks = mbreaks + 1
          breaks(mbreaks) = n1
!
! Find the last atom that has been defined.  Use this for defining the break
!
          do i = 1, numat
            if (.not. ioptl(i)) exit
          end do
          if (i - 1 > 0) break_coords(:,mbreaks) = coord(:,i - 1)          
        end if             
      end do
!
!  Re-evaluate all residues
!
      j = 1
      allres(j) = txtatm(1)(18:20)
      do i = 2, natoms
        if (txtatm(i) == " ")exit
        if (nat(i) /= 1 .and. txtatm(i)(23:27) /= txtatm(i - 1)(23:27)) then
        j = j + 1
        allres(j) = txtatm(i)(18:20)
        end if
      end do
      ires = j      
      iopt(:natoms) = ib(:natoms)
!
!   LABEL THE ATOMS IN ANY NON-PROTEIN MOLECULES IN THE SYSTEM
!
      nfrag = nfrag + 1  
      if (start_res(max(1,nfrag)) == -200) then
        if (.not. l_protein) ires = 0       
      else
        ires = start_res(nfrag)         
      end if
      call lewis(.true.)
      call ligand (ires, start_res, nfrag)
!
!  If ligands are present, set a break at the end of the protein, to separate protein from ligands
!
      if (.not. lbreaks) then
        mbreaks = mbreaks + 1
        breaks(mbreaks) = n1
!
! Find the last atom that has been defined.  Use this for defining the break
!
        do i = 1, numat
          if (.not. ioptl(i)) exit
        end do
        if (i - 1 > 0) break_coords(:,mbreaks) = coord(:,i - 1)          
      end if         
      nres = uni_res
      maxtxt = txtmax
      if (.not. pdb_label) pdb_label = (maxtxt > 25)
!
!  Add chain letters
!
      call reset_breaks()
      i = index(keywrd, " RESI")
      if (i > 0) then
        j = index(keywrd(i + 1:), " ") + i
        j = index(keywrd(i + 1: j), "0")
      end if
      if (i > 0 .and. j > 0) then
        do i = 1, numat
          txtatm(i)(13:16) = txtatm1(i)(13:16)
        end do
      end if
      if (moperr) return
      mbreaks = 1
      line = txtatm(1)(23:txtmax)
      do i = 1, numat
        txtatm(i) = txtatm(i)(:21)//chains(mbreaks)//txtatm(i)(23:)
        if (i == breaks(mbreaks)) then
          mbreaks = mbreaks + 1
        else  
          if (txtatm(i)(23:txtmax) == "    ") cycle          
          if (mbreaks > 1) then
!
!  This is a workaround to avoid having a break that already has been used.
!  (This can occur when a heterogroup that contains a peptide bond is positioned
!   after other heterogroups that do not contain a peptide bond.)
!
            if (i - 1 == breaks(mbreaks - 1)) then
              line = txtatm(i)(23:txtmax)
              cycle
            end if
          end if
          if (txtatm(i)(:6) == "HETATM") then
            if (txtatm(i)(23:26) /= line(1:4)) then
              if (txtatm(i - 1)(:4) /= "ATOM") then
                if (txtatm(i - 1)(18:20) /= "HOH" .or. txtatm(i)(18:20) /= "HOH") then
                  mbreaks = min(26, mbreaks + 1)
                end if
              end if
            end if
            txtatm(i) = txtatm(i)(:21)//chains(mbreaks)//txtatm(i)(23:)            
          end if
        end if
        line = txtatm(i)(23:txtmax)
      end do 
!
!   Check for unknowns
!
      do i = 1, numat
        do j = 1, 20
          if (txtatm(i)(18:20) == tyres(j)) exit
        end do
        if (j == 21) then
          if (txtatm(i)(15:16) /= "  ") cycle
          j = 1
          do k = i + 1, numat
            if (txtatm(k)(14:txtmax) == txtatm(i)(14:txtmax)) then
              j = j + 1
              if (j < 10) then
                write(txtatm(k)(15:15),'(i1)')j
              else
                 write(txtatm(k)(15:16),'(i2)')min(j, 99)
              end if                
            end if
          end do
          if (j > 1) write(txtatm(i)(15:15),'(i1)')1
        end if 
      end do
      do i = 1, natoms 
        if (txtatm(i)(1:6) /= "HETATM" .and. txtatm(i)(1:6) /= "ATOM  ") then
!
!  Catch all atoms that have still not been identified
!
          write(line, '(a,i5, a, a)')"HETATM", i, " "//elemnt(labels(i)), "   HET"
          txtatm(i) = trim(line)
        end if  
      end do
      if (n_new /= 0) then
!
!  Re-name residues to use the XENO name
!
        call update_txtatm(.true., .true.)
        mbreaks = 1
        old_name(:n_new) = "---"
        l_names = .false.
        do i = 1, natoms
          j = nint(reada(txtatm(i)(23:),1))
          do k = 1, n_new
            if (new_chain(k) == chains(mbreaks) .and. new_res(k) == j) then
              if (old_name(k) == "---") old_name(k) = txtatm(i)(18:20)
              txtatm(i) = txtatm(i)(:17)//new_name(k)//txtatm(i)(21:)
              if (txtatm(i)(:6) /= "HETATM") l_names(k) = .true.
              exit
            end if
          end do
          ter = (i == breaks(mbreaks))
          if (ter) mbreaks = mbreaks + 1
        end do
        first = .true.
        do i = 1, n_new
          if (l_names(i)) then
            if (first) then
              write(iw,'(/,a)') "      Residue names that have been changed"
              write(iw,'(/,a)') "      Residue No.  Calculated name   XENO name"
              first = .false.
            end if
            if (old_name(i) ==  new_name(i)) then
              write(iw,'(i3,i9,3x,a1,8x,a3,13x,a3, a)')i, new_res(i), new_chain(i), old_name(i), new_name(i), &
                "  Name not changed!"
            else
              write(iw,'(i3,i9,3x,a1,8x,a3,13x,a3)')i, new_res(i), new_chain(i), old_name(i), new_name(i)
            end if 
          end if
        end do
        if (log) then
          write(ilog,'(/,a)') "      Residue names that have been changed"
          write(ilog,'(/,a)') "      Residue No.  Calculated name   XENO name"
          do i = 1, n_new
            if (old_name(i) ==  new_name(i)) then
              write(ilog,'(i3,i9,3x,a1,8x,a3,13x,a3, a)')i, new_res(i), new_chain(i), old_name(i), new_name(i), &
              "  Name not changed!"
            else
              write(ilog,'(i3,i9,3x,a1,8x,a3,13x,a3)')i, new_res(i), new_chain(i), old_name(i), new_name(i)
            end if        
          end do
        end if
      end if  
      l_use_old_labels = (index(keywrd," SITE=") /= 0 .and. index(keywrd, " ADD-H") == 0)
      l_use_old_labels = .true.
      call update_txtatm(l_use_old_labels, .true.)
      call write_sequence
    end if
90   if (index(keywrd, " PDBOUT") /= 0) then
      allocate(temp_txtatm(natoms))
!
! Assign atom numbers
!
      temp_txtatm = txtatm
      do i = 1, natoms
        write (txtatm(i), "(a,i5,a)")temp_txtatm(i)(1:6),i,temp_txtatm(i)(12:maxtxt)         
      end do
      deallocate(temp_txtatm)
      
!
!  Identify atoms where chain breaks occur
!
      if (Index (keywrd, " RESEQ") + index(keywrd, " ADD-H") /= 0) call reset_breaks()
      if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") == 0) then
        if (allocated(txtatm1) .and. index(keywrd, " RESID") /= 0) then
          if (txtatm1(1) /= " ") call compare_sequence(n_new)
        end if
      end if
    end if
!
!  Edit keywords to remove text that would not be used in the next calculation. 
!
    call delete_ref_key("SITE", len_trim("SITE"), ') ', 2)
    i = index(keywrd, "SITE=()")
    if (i /= 0) keywrd(i:i + 6) = " "
    if (index(keywrd, " RESEQ") /= 0) call delete_ref_key("RESEQ", len_trim("RESEQ"), ' ', 1)
    numat = nnumat
    natoms = max(natoms, numat)
    if (index(keywrd, " ADD-H") /= 0) return
    if (index(keywrd, "CHARGES") == 0 .and. index(keywrd, "CONTROL_no_MOZYME") /= 0 .or. index(keywrd, " RESEQ") /= 0) then
      if ( index(keywrd," PDBOUT") /= 0) then
        line = archive_fn(:len_trim(archive_fn) - 3)//"pdb"
        i = iarc
        inquire(unit=i, opened=opend) 
        if (opend) close(i, status = 'keep', iostat=j)  
        open(unit=i, file=line, status='UNKNOWN', position='asis') 
        rewind i 
        call pdbout(i)
        close (i)
      end if
      if (index(keywrd, " RESEQ") /= 0) then
        line = archive_fn(:len_trim(archive_fn) - 3)//"arc"
        inquire(unit=iarc, opened=opend) 
        if (opend) close(iarc, status = 'keep', iostat=i)  
        open(unit=iarc, file=line, status='UNKNOWN', position='asis') 
        rewind iarc 
        call geout (iarc)
      end if
      if (index(keywrd, " RESEQ") == 0) return
    end if
    ibad = 0
    if (.not. mozyme) goto 99
!
!  MODIFY IONS SO THAT IT REFERS TO ALL ATOMS (REAL AND DUMMY)
!
    j = 0
    iz = ions
    ions = 0
    do i = 1, natoms
      if (labels(i) == 99) then
        ions(i) = 0       
      else
         j = j + 1
        ions(i) = iz(j)
      end if
    end do
    if (lres .and. .not. lreseq) then
!
!   CHECK ALL IONS TO SEE IF ANY RESIDUE IS A DI-ION
!
      outer_loop: do i = 1, numat
        if(ions(i) /= 1 .and. ions(i) /= -1) cycle! WARNING
      end do outer_loop
      if (ibad /= 0) then
        write (iw,*)
      end if
      jbad = ibad
      ibad = 0
      ibad = ibad + jbad
    end if
!
!
    charges = (lsite .or. index(keywrd, "CHARGES") /= 0)
    ichrge = -noccupied*2
    do i = 1, numat
      ichrge = ichrge + nint(tore(nat(i)))
    end do
    line = " "
    if (noccupied /= 0 .and. (index(keywrd," LEWIS") > 0 .or. noccupied*2 /= nelecs)) then
      maxtxt = 0
      do i = 1, numat
        maxtxt = max(maxtxt, len_trim(txtatm(i)))
      end do
      if (maxtxt == 0) then
        j = 1
      else
        j = maxtxt/2 + 2
      end if
      if (lreseq) call update_txtatm(.true., .false.)
      if (prt_topo) then
          line = " "
        write (iw, "(/,A,/)") "   TOPOGRAPHY OF SYSTEM"
        write (iw,*) "  ATOM No. "//line(:j)//"  LABEL  "//line(:j)//"Atoms connected to this atom"
        if (j == 0) then
          do i = 1, numat
            write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // "  ", (ibonds(j, i), j=1, nbonds(i))
          end do
        else
          if (maxtxt > 2) then
            do i = 1, numat
              write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)) // " (" // txtatm(i) (:maxtxt) // ") ", &
                (ibonds(j, i), j=1, nbonds(i))
            end do
          else
            do i = 1, numat
              write (iw, "(I7,9X,A,9I7)") i, elemnt (nat(i)), (ibonds(j, i), j=1, nbonds(i))
            end do
          end if        
        end if 
      end if
    end if
    if (noccupied /= 0 .and. index(keywrd," LEWIS") > 0) then
      write (iw, "(/37x,A,/)") "   Lewis Structure"
      if (index(keywrd, " LARGE") /= 0)  &
        write (iw,"(23x,a,/)") "  ATOMS IN OCCUPIED LOCALIZED MOLECULAR ORBITALS"
      l = 4
      allocate(Lewis_formatted(Lewis_tot/l + 5, l))
      Lewis_formatted = " "
      k = 0
      j = 0
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) j = j + 1
      end do
      m = j/l + 1    
      write(iw,"(4('     LMO  Atom  Atom    '))")  
      ii = 0
      jj = 1  
      do i = 1, Lewis_tot
        if (Lewis_elem(1,i) > 0) then
          k = k + 1
          ii = ii + 1
          if (ii > m) then
            ii = 1
            jj = jj + 1
          end if
          if (Lewis_elem(2,i) > 0) then
            write (Lewis_formatted(ii,jj), "(I8,I6,I6)") k, (Lewis_elem(j,i),j=1,2)
          else
            write (Lewis_formatted(ii,jj), "(I8,I6)") k, Lewis_elem(1,i)
          end if
        end if
      end do
      do i = 1, Lewis_tot
        if (Lewis_formatted(i,1) == " ") exit
        write(iw,"(10a)")(Lewis_formatted(i,j)//"    ", j = 1, l)," "
      end do
      if (index(keywrd, " LARGE") /= 0) then
        write (iw,"(/23x,a,/)") "  ATOMS IN UNOCCUPIED LOCALIZED MOLECULAR ORBITALS"
        Lewis_formatted = " "
        j = 0
        do i = 1, Lewis_tot
          if (Lewis_elem(2,i) > 0) j = j + 1
        end do
        m = j/l + 1
        k = 0
        ii = 0
        jj = 1  
        do i = 1, Lewis_tot
          if (Lewis_elem(2,i) > 0) then
            k = k + 1
            ii = ii + 1
            if (ii > m) then
              ii = 1
              jj = jj + 1
            end if
            if (Lewis_elem(1,i) > 0) then
              write (Lewis_formatted(ii,jj), "(I8,I6,I6)") k, (Lewis_elem(j,i),j=1,2)
            else
              write (Lewis_formatted(ii,jj), "(I8,I6)") k, Lewis_elem(2,i)
            end if
          end if
        end do
        do i = 1, Lewis_tot
          if (Lewis_formatted(i,1) == " ") exit
          write(iw,"(10a)")(Lewis_formatted(i,j)//"    ", j = 1, l)," "
        end do
      end if
      deallocate(Lewis_formatted)
      write (iw,*)
      if (noccupied > 0) then
        write (iw,"(a,/)") "          Type          Number of Lewis structural elements identified"
        if (log) write (ilog,"(/,a,/)") "          Type          Number of Lewis structural elements identified"
      end if
      if (numbon(1) /= 0) then
        write (iw, "(A,I6)") "         SIGMA BONDS   ", numbon (1)
        if (log) write (ilog, "(A,I6)") "         SIGMA BONDS   ", numbon (1)
      end if
      if (numbon(2) /= 0) then
        write (iw, "(A,I6)") "         LONE PAIRS    ", numbon (2)
        if (log) write (ilog, "(A,I6)") "         LONE PAIRS    ", numbon (2)
      end if
      if (numbon(3) /= 0) then
        write (iw, "(A,I6)") "         PI BONDS      ", numbon (3)
        if (log) write (ilog, "(A,I6)") "         PI BONDS      ", numbon (3)
      end if
      if (index(keywrd," LEWIS") > 0 .or. (noccupied*2 /= nelecs .and. noccupied /= 0)) then
        write(iw,"(/,a,i6)")" Number of filled levels from atoms and charge:", nelecs/2
        write(iw,"(a,i6)")" Number of filled levels from Lewis structure: ", noccupied
        if (log) write (ilog,"(/,a,i6)")" Number of filled levels from atoms and charge:", nelecs/2
        if (log) write (ilog,"(a,i6)")" Number of filled levels from Lewis structure: ", noccupied
      end if
      l = 0
      m = 0
      num = char(Int(log10(numat + 1.0)) + ichar("1") + 1) 
      do i = 1, numat
        if (.not. main_group(nat(i))) then
!
!  Element is a transition metal.  Work out its formal oxidation state
!
          k = ions(i)
          do j = 1, Lewis_tot
            if (Lewis_elem(1,j) /= 0 .and. Lewis_elem(2,j) /= 0) then
              if (Lewis_elem(1,j) == i) k = k + 1
              if (Lewis_elem(2,j) == i) k = k + 1
            end if
          end do
          if (m == 0) write(iw,*)
          m = 1
          write(iw,"(10x,a,i"//num//",a,sp,i3)")" Formal oxidation state of atom",i,", a "//trim(elemnt(nat(i)))//", is", k
          if (k < 0) l = 1
          if (k > 3) l = 1
        end if
      end do
      if (l == 1) call web_message(iw,"Lewis_structures.html")
    end if    
!
! Check for sulfate and phosphate
!
      do i = 1, numat
        if (nat(i) == 16 .and. txtatm(i)(18:20) == "SO4") then
          k = 2
          ions(i) = 0
          do j = 1, nbonds(i)
            l = ibonds(j,i)
            if (ions(l) == -1) then
              ions(l) = 0
              k = k - 1
              if (k == 0) exit
            end if
          end do
        end if
        if (nat(i) == 15 .and. txtatm(i)(18:20) == "PO4") then
          k = 1
          ions(i) = 0
          do j = 1, nbonds(i)
            l = ibonds(j,i)
            if (ions(l) == -1) then
              ions(l) = 0
              k = k - 1
              if (k == 0) exit
            end if
          end do
        end if
      end do
    num_ions = 0
    do i = 1, numat
      j = (Min(6, Max(-6, ions(i))))
      num_ions(j) = num_ions(j) + 1
    end do
    i = 0
    do j = 1,6
      i = i + num_ions(j) + num_ions(-j)
    end do
    if (i == 0) then
      write(iw,'(/10x, a)')"NO CHARGES FOUND."
    else      
      if (index(keywrd," LEWIS") > 0) then
        write (iw,*)
        write (iw,"(a,/)") "          Type           Number of charged sites identified"
        if (log) write(ilog,"(/,a,/)") "          Type           Number of charged sites identified"
        do i = 1,6
        if (num_ions(i)  > 0) then
          write (iw, "(9x,a,2x,i6)") ion_names(i) , num_ions(i)
          if (log) write(ilog, "(9x,a,2x,i6)") ion_names(i) , num_ions(i)
        end if
        if (num_ions(-i) > 0) then
          write (iw, "(9x,a,2x,i6)") ion_names(-i), num_ions(-i)
          if (log) write(ilog, "(9x,a,2x,i6)") ion_names(-i), num_ions(-i)
        end if
        end do
        i = 0
        do j = 1,6
          i = i + num_ions(j) * j
        end do
        write(iw,"(SP/,a,i5)")" SUM OF POSITIVE CHARGES", i
        if (log) write(ilog,"(SP/,a,i5)")" SUM OF POSITIVE CHARGES", i
        i = 0
        do j = 1,6
          i = i + num_ions(-j) * j
        end do
        write(iw,"(a,i5)")" SUM OF NEGATIVE CHARGES", -i
        if (log) write(ilog,"(a,i5)")" SUM OF NEGATIVE CHARGES", -i
      end if
      padding = " "
      do i = 1, numat
        if (ions(i) /= 0) then
          if (nat(i) == 15 .or. nat(i) == 16) then
            kk = 0
            do jj = 1, nbonds(i)
              ii = ibonds(jj,i)
              if (nat(ii) == 8) then
                if (ions(ii) == -1) then
!
!  Found a PO4 or SO4.  Neutralize the P-O or S-O Zwitterion.
!
                  ions(i) = ions(i) - 1
                  ions(ii) = ions(ii) + 1
                  kk = 1
                  if (ions(i) == 0) exit
                end if
              end if
            end do
            if (kk == 1) cycle
          end if
        end if
      end do
      do i = 1, numat
        if (ions(i) /= 0) exit
      end do
      maxtxt_store = maxtxt
      if (maxtxt < 0) maxtxt = 14
      l_salt = .false.
      if (first_prt .and. .not. l_rama) then
        if (i <= numat) then
          if (maxtxt > 1) then        
             l = max(1,(17 - maxtxt/2))
             residues = (index(keywrd, " RESID") /= 0)
          end if
        else
          write(iw,'(/18x,a)') "NO CHARGED ATOMS FOUND."        
        end if
        ioptl = .true.
        m = 0
        header = .true.
        do n1 = 1, 8
          j = atomic_charges(n1)
          first = .true.
          k = 0
          kkk = 0
          do i = 1, numat
            if (ions(i) == j) then
              if (first) then
                write(iw,*)
                first = .false.
                end if
              line = " "
              jj = 0
              if (j == 1) then
                line = " "
                do ii = 1, numat
                  if (ions(ii) < 0 .and. ioptl(ii)) then
                    call identify_nearby_counterions(i, ii, jj, near_ions, r_ions, line)   
                    if (line == "Salt bridge!") then
                      ioptl(i) = .false.
                      ioptl(ii) = .false.
                      exit
                    end if
                  end if
                end do
              else if (j == -1) then
                if (.not. ioptl(i)) cycle
                line = " "
                do ii = 1, numat
                  if (ions(ii) > 0 .and. ioptl(ii)) then
                    call identify_nearby_counterions(i, ii, jj, near_ions, r_ions, line)   
                    if (line == "Salt bridge!") exit
                  end if
                end do
              end if
              if (jj > 1) then
  !
  !  Sort near ions into increasing distance
  !
                jj = min(jj,100)
                do ii = 1, jj
                  do kk = ii + 1, jj
                    if (r_ions(kk) < r_ions(ii)) then
                      sum = r_ions(kk)
                      r_ions(kk) = r_ions(ii)
                      r_ions(ii) = sum
                      near_ions_store = near_ions(:,kk)
                      near_ions(:,kk) = near_ions(:,ii)
                      near_ions(:,ii) = near_ions_store
                    end if
                  end do 
                end do
                do ii = 2, jj
                  if (r_ions(ii) > 5.d0) then
                    jj = ii - 1
                    exit
                  end if
                end do
              end if
              if (maxtxt > 1) then
                if (line == "Salt bridge!") then
                  if (j == 1) then
!
! There is a counterion, therefore this is a salt-bridge
!
                    if (residues) then
                      txtatm_1 = txtatm(near_ions(1,1))  
                      txtatm_2 = txtatm(near_ions(2,1))
                    else
                      txtatm_1 = txtatm1(near_ions(1,1))
                      txtatm_2 = txtatm1(near_ions(2,1))
                    end if
                    num = "5"
                    if (elemnt(nat(near_ions(2,1)))(1:1) /= " ") num = "5"
                      if (txtatm_2(maxtxt:maxtxt) == ")") then
                        ii = maxtxt + 1
                        txtatm_2 = "("//trim(txtatm_2)
                        kk = maxtxt - 1
                      else
                        ii = maxtxt
                        kk = maxtxt
                      end if
                      if ( .not. l_salt) then
                        if (maxtxt > 1) then        
                           write(line,"(a)")"   Ion     <--------PDB Label------->   Charge     "// &
                             &"Distance        PDB label for Salt Bridge   Charge"
                           write(iw,'(/,a,/)') trim(line)
                           if (log) write(ilog,"(/,a)") trim(line)
                        else
                           write(iw,"(/,a)")"     Ion Atom No.  Type    Charge"
                           if (log) write(ilog,"(/,a,/)")"    Ion Atom No.  Type    Charge"
                        end if
                      end if
                      k = k + 1
                      write(tmp,"(i5, 2x, 3x, a, SP,i5,S, f13.2, 9x, a, a, a)") &
                      & k, padding(:l-5)//"("//txtatm_1(1:kk)//")"//padding(:l-15), ions(i), &
                      & r_ions(1), "("//txtatm_2(:ii)//")", "   -1"
                      write(iw,'(a)') trim(tmp)
                      l_salt = .true.
                      if (log) write(ilog,"(a)") trim(line)
                  end if
                else
!
! No nearby counterions, therefore this is not a salt-bridge
!
                   if (residues) then
                    txtatm_1 = txtatm(i)
                  else
                    txtatm_1 = txtatm1(i)
                  end if
                  kk = len_trim(txtatm_1)
                  if (txtatm_1(kk:kk) == ")") then
                    kk = maxtxt - 1
                    ii = l - 6
                    j2 = l - 1
                  else
                    kk = maxtxt
                    ii = l - 5
                    j2 = l - 5
                  end if
                  kkk = kkk + 1
                  write(line,"(i5, 2x, 3x, a,SP,i5,S,a)") &
                    kkk, padding(:ii)//"("//txtatm_1(1:kk)//")"//padding(:j2), ions(i)
                  if (m == 0) then
!
! Find the first free unit, and use that number for the scratch file
!
                    do m = 10, 100
                      inquire(unit=m, opened = opend) 
                      if (.not. opend) then
                      open(unit=m, status='SCRATCH') 
                      exit
                      end if
                    end do
                  end if
                  if (j /= 1) then
                    if (header) then
                      if (maxtxt < 2) then
                        write(iw,'(10x, a, /)') "All other ions"
                        if (maxtxt < 2)write(iw,"(a)")"   Ion Atom No.  Type    Charge"
                      else
                        write(iw,'(/17x, a)') "All other ions"
                      end if
                      header = .false.
                    end if
                    write(iw,'(a)')trim(line)
                    else
!
! Store text of ions with a charge of +1
!
                  write(m,'(a)')trim(line)
                  end if
                  if (log) write(ilog,'(a)')trim(line)
                end if
              else
                k = k + 1
                write(line,"(i5,3x,i5,5x,a2,SP,i9,S)")k, i, elemnt(nat(i)),ions(i)
                if (header) then
                  if (maxtxt < 2) then
                    write(iw,'(10x, a, /)') "All other ions"
                    if (maxtxt < 2)write(iw,"(a)")"   Ion Atom No.  Type    Charge"
                  else
                    write(iw,'(/17x, a)') "All other ions"
                  end if
                  header = .false.
                end if
                write(iw,'(a)') trim(line)
                if (log) write(ilog,"(a)")trim(line)
              end if            
            end if
          end do
          if (j == 2 .and. m /= 0) then
!
!  Write out ions that have a charge of -1
!
            rewind (m)
            if (header) then
              if (maxtxt < 2) then
                    write(iw,'(10x, a, /)') "All other ions"
                    if (maxtxt < 2)write(iw,"(a)")"   Ion Atom No.  Type    Charge"
                  else
                    write(iw,'(/17x, a)') "All other ions"
                  end if
              header = .false.
            end if
            do k = 1, 10000
              read(m,'(a)', iostat = ii)line
              if (ii /= 0) exit
              if (k == 1) write(iw,*)
              write(iw, '(a)')trim(line)
            end do
            close (m)
          end if
        end do
      end if
      maxtxt = maxtxt_store
    end if
    if (first_prt .and. noccupied /= 0 .and. .not. l_rama) then
      num = char(ichar("3") + max(int(log10(abs(ichrge) + 0.05)),0))
      if (index(keywrd," CHARGE=") /= 0) then
        if (ichrge == irefq) then
          write (iw, "(SP/10x,A,I"//num//", a)") "COMPUTED CHARGE ON SYSTEM:", ichrge, &
         ", THIS IS THE SAME AS THE CHARGE DEFINED IN THE DATA-SET."
        else
          i = index(keywrd," CHARGE=") + 1
          j = i -2 + index(keywrd(i:), " ")
          write (iw, "(SP/2x,A,I"//num//", a)") "COMPUTED CHARGE ON SYSTEM:", ichrge, &
         ", THIS IS NOT THE SAME AS THE CHARGE DEFINED IN THE DATA-SET: """//keywrd(i:j)//"""."
        end if
      else
        if (.not. charges) write (iw, "(SP/10x,A,I"//num//")") "COMPUTED CHARGE ON SYSTEM:", ichrge
      end if      
      if (log) write (ilog, "(SP/10x,A,I"//num//")") "COMPUTED CHARGE ON SYSTEM:", ichrge
      end if
      first_prt = .false.      
    if (index(keywrd, " LEWIS") /= 0 .and. index(keywrd," 0SCF") == 0) &
      call mopend ("RUN STOPPED BECAUSE KEYWORD LEWIS WAS USED.")
!
99  continue  
    atmass(1:numat) = ams(nat(1:numat))
    if (done .and. .not. lreseq) then
      call xyzint (coord, numat, na, nb, nc, 1.d0, geo)
    end if
    if (lreseq) then
      geo(:,:numat) = coord(:,:numat)
      na = 0
    end if   
    if (lreseq .or. lsite) then    
      archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"arc"
      inquire(unit=iarc, opened=opend) 
      if (opend) close(iarc, status = 'keep', iostat=i)  
      open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis') 
      rewind iarc 
      if (index(keywrd, " PDBOUT") /= 0) call delete_ref_key("PDBOUT", len_trim("PDBOUT"), ' ', 1)
      call geout (iarc)
      if (index(keywrd, " PDBOUT") /= 0) then        
        line = archive_fn(:len_trim(archive_fn) - 3)//"pdb" 
        i = iarc
        inquire(unit=i, opened=opend) 
        if (opend) close(i)
        open(unit=i, file=line, status='UNKNOWN', position='asis') 
        rewind i 
        coord(:,:numat) = geo(:,:numat)
        call update_txtatm(.true., .true.) 
        if (index(keywrd, " NORES") == 0) call rectify_sequence()
        call reset_breaks()
        call pdbout(i)
        if (allocated(coorda)) deallocate (coorda)
        allocate(coorda(3,numat))
        txtatm1(:numat) = txtatm(:numat)
        coorda(:,:numat) = coord(:,:numat)
        close (i)
      end if
    end if
    if (.not. mozyme) goto 98
    if (irefq /= ichrge .and. .not. lreseq .or. charges .or. index(keywrd, " Move") /= 0) then
!
!  THE CALCULATED CHARGE DOES NOT MATCH THAT DEFINED BY CHARGE=N.
!  THEREFORE, THE USER HAS MADE A MISTAKE.  WRITE OUT CHARGES
!  FOUND HERE.
!
!
      if (icharges /= 0) write (iw,*)   
      line = " "
      if (index(keywrd, " 0SCF") == 0) line = "JOB STOPPED BECAUSE"
      i = len_trim(line)
      if (i > 0) i = i + 1
      j = max(abs(ichrge), abs(irefq))
      k = min(ichrge, irefq)
      num = char(ichar("2") + max(int(log10(j + 0.05)),0))
      if (k < 0) num = char(ichar(num) + 1)
      if (charges) then
        if (index(keywrd," CHARGES") +index(keywrd," CHARGE=") == 0) then
          call mopend (line(:i)//"CHARGES MODIFIED BY SITE COMMAND")
        else if (noccupied > 0) then
          call write_sequence
          if (irefq /= ichrge) then
            if (index(keywrd," CHARGE=") /= 0) then
              write (iw, "(10x,A,SP,I"//num//",A)") "CHARGE SPECIFIED IN DATA SET: ", irefq," IS INCORRECT." 
            else
              write (iw, "(/10x,A,SP,I"//num//",a)") "COMPUTED CHARGE ON SYSTEM: ", ichrge, &
                "  (THIS DOES NOT AGREE WITH THE DEFAULT CHARGE OF ZERO)"
            end if
          else
            if (.not. charges) write (iw, "(SP/3x,A,I"//num//", a)") "COMPUTED CHARGE ON SYSTEM = ", ichrge, &
            ". THIS AGREES WITH THE CHARGE IN THE DATA-SET."             
          end if
          archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"arc"
          inquire(unit=iarc, opened=opend) 
          if (opend) close(iarc, status = 'keep', iostat=i)  
          open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis') 
          rewind iarc 
          if (index(keywrd, " PDBOUT") /= 0) call delete_ref_key("PDBOUT", len_trim("PDBOUT"), ' ', 1)
          call geout (iarc)
          call mopend ("JOB STOPPED HERE BECAUSE KEYWORD ""CHARGES"" WAS USED")
        end if
        call delete_ref_key("SITE", len_trim("SITE"), ') ', 2)
      end if    
      if (lreseq .or. index(keywrd, " Move") /= 0) then
!
!  Check to see if any atoms are in internal coordinates
!
        j = 0
        do i = 1, numat
          if (na(i) /= 0) j = j + 1
        end do
        if (j /= 0) then
          if (index(keywrd, " Move") /= 0) then
            call mopend("WHEN HYDROGEN ATOMS ARE ADDED OR DELETED, ALL ATOMS MUST BE IN CARTESIAN COORDINATES")
            return
          else if (index(keywrd, "GEO-OK") == 0) then
            call mopend("WHEN SYSTEM IS RESEQENCED, EITHER ADD ""GEO-OK"" OR ALL ATOMS MUST BE IN CARTESIAN COORDINATES")
          end if
        end if
        call store_and_restore_Tv("STORE")
        call move_hydrogen_atoms
        call store_and_restore_Tv("RESTORE")
        call lewis(.false.)
        if (index(keywrd, " ADD-H") > 0) then
          call mopend("ADD-H: SYSTEM HAS BEEN HYDROGENATED")
        else if (index(keywrd, " Move") > 0) then
          call mopend("HYDROGEN ATOMS ADDED OR DELETED")
        else
          call mopend("GEOMETRY RESEQUENCED")
        end if
        moperr = .false.
        if (prt_coords) call geout (iw)
      end if
      if (charges .or. lreseq) return
      if (index(keywrd, " LEWIS") /= 0) then
        if (index(keywrd, " 0SCF") == 0) then
          if (.not. moperr) call mopend (line(:i)//"KEYWORD LEWIS USED")
          return
        end if
      end if 
!
!
      if (Index (keywrd, " 0SCF") + Index (keywrd, " RESEQ") + Index (keywrd, " LEWIS") == 0) then    
        write(iw,*)" "
        if (index(keywrd," CHARGE=") /= 0) then
          if (index(keywrd," GEO-OK") /= 0) then
            write (iw, "(10x,A)") "KEYWORD 'GEO-OK' WAS PRESENT, SO THE CHARGE HAS BEEN RESET."
            write (iw, "(10x,A)") "IF THE NEW CHARGE IS INCORRECT, EITHER MODIFY THE STRUCTURE OR "//&
            "USE KEYWORD 'SETPI' TO CORRECT THE LEWIS STRUCTURE." 
          else
            write (iw, "(10x,A,SP,I"//num//",A)") "In the data-set supplied, the charge specified (", irefq,") is incorrect."  
          end if
        else
          num = char(ichar("3") + max(int(log10(abs(ichrge) + 0.05)),0))
          write (iw, "(10x,A,SP,I"//num//",a)") "BECAUSE KEYWORD ""CHARGE"" WAS NOT USED, THE COMPUTED CHARGE OF", &
          ichrge, " WILL BE USED IN THIS RUN."
          sum = -dfloat(ichrge)/norbs
          pdiag(:norbs) = pdiag(:norbs) + sum
        end if             
        nelecs = nelecs + irefq - ichrge     
        nclose = nelecs/2
        nopen = nclose
        nalpha = 0
        nbeta = 0
        uhf = .false.     
        if (irefq /= ichrge) then
          if (mod(irefq - ichrge, 2) == 0) then
            write(iw,"(/10x,a)")"If the Lewis structure used by MOZYME is incorrect, "// &
            &"use keywords such as CVB or SETPI to correct it"
            call web_message(iw,"setpi.html")
          else
            if (index(keywrd," CHARGE=") /= 0) &
              write(iw,"(/10x,a)")"The charge keyword, Lewis structure, or the chemical formula is faulty"
          end if            
          if (Index(keywrd," GEO-OK") == 0 .and. index(keywrd," CHARGE=") /= 0 ) then  
            call mopend("CHARGE SPECIFIED IS INCORRECT. CORRECT THE ERROR BEFORE CONTINUING")
            write(iw,"(/10x,a)")"If that is done, then the correct charge will be used."
            return  
          else if (id > 0) then
            write(iw,"(/10x,a)")"INFINITE SYSTEMS MUST HAVE A ZERO CHARGE ON THE UNIT CELL"
            call mopend("Unit cell has a charge. Correct fault and re-submit ")
            return
          else 
            call fix_charges(ichrge)  
          end if
        end if
      end if
    end if
98  if (done) then
      if (prt_coords) write (iw, "(//10X,A,//)") " GEOMETRY AFTER RE-SEQUENCING"
      call update_txtatm(.true., .true.) 
      if (prt_coords) call geout (iw)
      if (index(keywrd, "0SCF") + index(keywrd, " RESEQ") == 0 .or. &
         index(keywrd, " PDBOUT") == 0) then
        inquire(unit=iarc, opened=opend) 
        if ( .not. opend) open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
        rewind(iarc)
        call geout (iarc)
      end if
      if ( index(keywrd," PDBOUT") /= 0) then
        if (prt_coords) call pdbout(1)
      end if
      line = "GEOMETRY RESEQUENCED"
      call mopend(trim(line))
      go to 1100
    end if
    if (ibad /= 0 .and. .not. let) then
      call mopend("ERROR")
      go to 1100
    end if
    if (times) then
      call timer (" END OF GEOCHK")
    end if
    if (lreseq .or. lsite) then     
      if (lsite) then
        call mopend ("RUN STOPPED BECAUSE KEYWORD ""SITE"" USED")
        call delete_ref_key("SITE", len_trim("SITE"), ') ', len_trim(') '))
      else if (lreseq) then
        if (index(keywrd, " ADD-H") > 0) then
          call mopend("ADD-H: SYSTEM HAS BEEN HYDROGENATED")
        else if (index(keywrd, " Move") > 0) then
          call mopend("HYDROGEN ATOMS ADDED OR DELETED")
        else
          call mopend("GEOMETRY RESEQUENCED")
        end if
      end if   
      moperr = .false.
      if (prt_coords) call geout (iw)
      return
    end if
!
!  MODIFY IONS SO THAT IT REFERS TO REAL ATOMS ONLY
!
    j = 0
    do i = 1, natoms
      if (labels(i) /= 99) then
        j = j + 1
        iz(i) = j
      end if
    end do
!
!  Restore charges, if present
!
    if (maxtxt /= 27) then
      do i = 1, natoms
        if(atom_charge(i) /= " ") txtatm(i)(2:2) = atom_charge(i)
      end do
    end if
!
!  Restore nbonds and ibonds in case they are modified within this subroutine
!
    nbonds(:numat) = nnbonds(:numat)
    ibonds(:,:numat) = iibonds(:,:numat)
1100 continue
    if (Allocated (nnbonds))      deallocate (nnbonds, iibonds)
    if (Allocated (iz))           deallocate (iz)
    if (Allocated (ib))           deallocate (ib)
    if (Allocated (mb))           deallocate (mb)
    if (Allocated (atom_charge))  deallocate (atom_charge)
    if (Allocated (ioptl))        deallocate (ioptl)
    return
  end subroutine geochk
  subroutine identify_nearby_counterions(i_in, j_in, n_jj, near_ions, r_ions, line)
!
! Identify charged atoms 
! If the charge is +1 or -1, then test for salt bridge
!
  use MOZYME_C, only : ions
  use common_arrays_C, only : txtatm, nbonds, ibonds, nat
  implicit none
  integer, intent(in) :: i_in, j_in
  integer, intent (inout) :: n_jj, near_ions(2, 100)
  double precision, intent (inout) :: r_ions(100)
  character, intent (out) :: line*100
!
  double precision :: sum, sum_min
  integer :: i, j, k, l, m, atom_i, atom_j, i_set(3), ni, j_set(3), nj
  character :: atm*3
  logical :: l_salt_i, l_salt_j
  double precision, external :: distance
  if (ions(i_in) == 1) then
    atom_i = i_in
    atom_j = j_in 
  else
    atom_j = i_in
    atom_i = j_in
  end if
  ni = 1
  i_set(1) = atom_i
  nj = 1
  j_set(1) = atom_j
  if (ions(atom_i)*ions(atom_j) == -1) then
!
!  This might be a salt-bridge.  Search for the shortest contact distance
!
! Look for cations
!
    if (txtatm(atom_i)(18:20) == "HIS") then
!
! Search for ND1 and NE2 in the imidazolium ring
!
      if (txtatm(atom_i)(14:16) == "CE1" .or. txtatm(atom_i)(14:16) == "ND1" .or. txtatm(atom_i)(14:16) == "NE2") then
        if (txtatm(atom_i)(14:16) == "CE1") then
!
!  Use the nitrogen atoms on each side of the current atom
!
          ni = 0
          do i = 1, nbonds(atom_i)
            j = ibonds(i,atom_i)
            if (nat(j) == 7) then
              ni = ni + 1
              i_set(ni) = j
            end if
          end do     
        else                                                  ! The following block has not been tested.
          if (txtatm(i_in)(14:16) == "ND1") then
            atm = "NE2"
          else
            atm = "ND1"
          end if          
!
!  Use the current nitrogen atom and the nitrogen atom on the other side of "CE1"
!
          do i = 1, nbonds(atom_i)
            j = ibonds(i,atom_i)
            if (txtatm(j)(14:16) == "CE1") then
              do k = 1, nbonds(j)
                l = ibonds(k,j)
                if (txtatm(l)(14:16) == atm) then    
                  ni = ni + 1
                  i_set(ni) = l
                end if          
              end do
            end if
          end do     
        end if
      end if
    else if (txtatm(atom_i)(18:20) == "ARG") then 
      if (txtatm(atom_i)(14:15) == "CZ" .or. txtatm(atom_i)(14:16) == "NH1" .or. txtatm(atom_i)(14:16) == "NH2") then
        if (txtatm(atom_i)(14:15) == "CZ") then
!
!  Use the terminal nitrogen atoms on each side of the current atom
!
          ni = 0
          do i = 1, nbonds(atom_i)
            j = ibonds(i,atom_i)
            ni = ni + 1
            i_set(ni) = j
          end do     
        end if
      end if
    end if
!
! Look for anionic oxygen atoms
!
   if (txtatm(atom_i)(24:26) == "107" .and. txtatm(atom_j)(24:26) == "245") then
    continue
   end if
    if (txtatm(atom_j)(14:14) == "O") then
      do i = 1, nbonds(atom_j)
        j = ibonds(i, atom_j)
        if (txtatm(j)(14:14) == "C" .and. nbonds(j) == 3) then
!
! Search for both oxygen atoms of a -COO 
!
          nj = 0
          do l = 1, nbonds(j)
            m = ibonds(l,j)
            if (txtatm(m)(14:14) == "O") then    
              nj = nj + 1
              j_set(nj) = m
            end if          
          end do
        end if
      end do
    end if  
  end if
  
!
! Now find the shortest distance between the two ions
!
  sum_min = 1.d10
  k = 0
  l = 0
  do i = 1, ni
    atom_i = i_set(i)
    do j = 1, nj
      atom_j = j_set(j)
      sum = distance(atom_i, atom_j)
      if (sum < sum_min) then
        sum_min = sum
        k = atom_i
        l = atom_j
      end if
    end do
  end do
  atom_i = k
  atom_j = l
  if (sum_min < 99.d0) then
    n_jj = n_jj + 1
    if (n_jj > 100) return
    if(ions(i_in) == 1) then
      near_ions(1,n_jj) = atom_i
      near_ions(2,n_jj) = atom_j
    else
      near_ions(1,n_jj) = atom_j
      near_ions(2,n_jj) = atom_i
    end if      
    r_ions(n_jj) = sum_min
    l_salt_i = .false.
    l_salt_j = .false.        
    if (txtatm(atom_i)(18:20) == "HIS") then
      if (txtatm(atom_i)(14:16) == "ND1" .or. txtatm(atom_i)(14:16) == "NE2") l_salt_i = .true.
    else if (txtatm(atom_i)(18:20) == "ARG") then
       if (txtatm(atom_i)(14:16) == "NH1" .or. txtatm(atom_i)(14:16) == "NH2" &
         .or. txtatm(atom_i)(14:15) == "NE") l_salt_i = .true.
    else if (txtatm(atom_i)(14:14) == "N") then
       if (nbonds(atom_i) == 4) l_salt_i = .true.
    end if
    if (txtatm(atom_j)(14:14) == "O") l_salt_j = .true.
    if (l_salt_i .and. l_salt_j .and. sum_min < 4.d0) line = "Salt bridge!"
  end if 
  return
  end subroutine identify_nearby_counterions
 
