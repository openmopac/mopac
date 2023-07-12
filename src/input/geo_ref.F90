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

    subroutine geo_ref
!
! geo_ref reads in the geometry reference file, and makes it as similar
! as possible to the current geometry.  It does this by:
!
! Requiring files to include PDB-type information
!
! The atoms of the reference geometry are put into the same sequence as in the current geometry.
!
! If any atoms cannot be mapped from reference to current geometry, then put those
! atoms at the end of the data sets.
!
! Rotate and translate the geometries to put them into maximum coincidence,
! i.e., minimize the RMS difference between the two geometries.
!
! During this process, swap pairs of atoms, if that will lower the RMS difference.
!
!
      use chanel_C, only : ir, iw, output_fn, job_fn, iarc, archive_fn
!
      use molkst_C, only : numat, keywrd, nvar, id, natoms, moperr, line, refkey, density, &
        maxtxt, numat_old, koment, title, geo_ref_name, geo_dat_name, arc_hof_2, arc_hof_1, &
        keywrd_txt, pdb_label, ncomments, refkey_ref, backslash, formula, keywrd_quoted
!
      use parameters_C, only : ams
!
      use common_arrays_C, only : geo, coord, xparam, na, nb, nc, nat, c, labels, txtatm, lopt, &
        txtatm1, geoa, atmass, loc, lopt_store
!
      use elemts_C, only: atom_names
      implicit none
!
      integer :: i, j, k, l, ii, jj, i4, j4, k4, iquit, numat_dat, numat_ref, ub
      integer, allocatable :: map_atoms_A(:), atom_no(:)
      integer, external :: quoted
      character, allocatable :: tmp_txt(:)*27, diffs(:)*80
      double precision, allocatable :: tmp_geoa(:,:)
      double precision :: dum1, dum2, sum, rms, rms_min, sum1, sum2, sum3, &
        toler, xmin, sum4
      double precision, external :: reada
      logical :: intern = .true., exists, bug, any_bug, swap, first, let, l_0SCF_HTML, opend
      logical, allocatable :: same(:), ok(:)
      character :: line_1*1000, line_2*1000, num*2, geo_dat*7
!
!   For Geo-Ref to work, some very specific conditions must be satisfied.
!   So before attempting a GEO_REF calculation, check that the data are okay
!
      do j = 2, natoms
        if (na(j) /= 0) exit
      end do
      l_0SCF_HTML = (index(keywrd," 0SCF") /= 0 .and. index(keywrd, " HTML ") /= 0)
!
! Force geometry to be Cartesian
!
      call gmetry (geo, coord)
      geo(:,:numat) = coord(:,:numat)
      na(:natoms) = 0
      j = 0
      do i = 1, natoms - id
        if (labels(i) /= 99) then
          j = j + 1
          labels(j) = labels(i)
          txtatm(j) = txtatm(i)
        end if
      end do
      txtatm1(:numat) = txtatm(:numat)
      id = 0
      if (moperr) return
      allocate(geoa(3,natoms + 300), c(3,natoms + 300)) ! Generous safety factor for second geometry
      i = quoted('GEO_REF=')
      if (i < -10) stop ! dummy use of "i" to prevent FORCHECK from flagging a possible error
      j = len_trim(line)
      if (line(j:j) == '"') line(j:j) = " "
      line_1 = trim(line)
      call upcase(line_1, len_trim(line_1))
      geo_ref_name = trim(line)
      i = index(keywrd_quoted," GEO_DAT")
      if (i > 0) then
        i = index(keywrd_quoted(i:), '"') + i
        j = index(keywrd_quoted(i + 2:),'" ') + i
        geo_dat_name = keywrd_quoted(i:j)
      else
        geo_dat_name = trim(job_fn)
      end if
      line = trim(geo_ref_name)
      inquire (file=trim(line), exist = exists)
      if (.not. exists) then
        call add_path(line)
        inquire (file=trim(line), exist = exists)
        if ( .not. exists) then
          inquire (file=trim(line)//".mop", exist = exists)
          if (exists) line = trim(line)//".mop"
        end if
        if ( .not. exists) then
          inquire (file=trim(line)//".txt", exist = exists)
          if (exists) line = trim(line)//".txt"
        end if
        if (.not. exists) then
          call mopend ("GEO_REF file '"//trim(line)//"' does not exist.")
          if (index(keywrd, " 0SCF") /= 0) then
            moperr = .false.
            call l_control("GEO_REF", len("GEO_REF"), -1)
            keywrd_txt = trim(keywrd)
          end if
          return
        end if
      end if
      open(unit = 99, file = trim(line), iostat = i)
      if (i /= 0) then
        call mopend ("Problem opening file "//trim(line))
        return
      end if
      if (l_0SCF_HTML) then
        if (index(geo_dat_name, "/") + index(geo_dat_name, backslash) == 0) then
          line = trim(geo_dat_name)
          call upcase(line, len_trim(line))
          line = line(len_trim(line) - 3:)
        else
          do i = len_trim(geo_dat_name), 1, -1
            if (geo_dat_name(i:i) == "/" .or. geo_dat_name(i:i) == backslash) exit
          end do
          geo_dat_name = geo_dat_name(i + 1:)
        end if
        line = geo_dat_name(:len_trim(geo_dat_name) - 3)//"pdb"
        call add_path(line)
        inquire (file=trim(line), exist = exists)
        if (exists .and. index(keywrd, " LET") == 0) then
           call mopend("A PDB FILE WITH THE SAME NAME AS THE GEO_DAT FILE ALREADY EXISTS")
           write(iw,'(10x,a)')"GEO_DAT FILE: """//trim(line)//""""
           write(iw,'(/10x,a)')"To continue, over-writing this file, add keyword ""LET"""
           return
        end if
        if (index(geo_ref_name, "/") + index(geo_ref_name, backslash) == 0) then
          line = trim(geo_ref_name)
          call upcase(line, len_trim(line))
          line = line(len_trim(line) - 3:)
        else
          do i = len_trim(geo_ref_name), 1, -1
            if (geo_ref_name(i:i) == "/" .or. geo_ref_name(i:i) == backslash) exit
          end do
          geo_ref_name = geo_ref_name(i + 1:)
        end if
        line = geo_ref_name(:len_trim(geo_ref_name) - 3)//"pdb"
        call add_path(line)
        inquire (file=trim(line), exist = exists)
        if (exists .and. index(keywrd, " LET") == 0) then
           call mopend("A PDB FILE WITH THE SAME NAME AS THE GEO_REF FILE ALREADY EXISTS")
           write(iw,'(10x,a)')"GEO_REF FILE: """//trim(line)//""""
           write(iw,'(/10x,a)')"To continue, over-writing this file, add keyword ""LET"""
        end if
        line = geo_dat_name
        i = len_trim(line) - 4
        call upcase(line, i)
        i = max(i, len_trim(line_1) - 4)
        if (line(:i) == line_1(:i)) then
          call mopend("The GEO_DAT and GEO_REF filenames are too similar.")
          write(iw,'(10x,a)')"(They would both write to the same PDB file.)"
        end if
        if (moperr) return
      end if
      nat(:natoms) = labels(:natoms)
      i = natoms
      rewind 99
      do i = 1, 100000
        read(99,"(a)", iostat = j) num
        if (j /= 0) exit
      end do
      rewind 99
      allocate(map_atoms_A(numat), tmp_txt(numat),  atom_no(numat))
      tmp_txt(:numat) = txtatm(:numat)
      if (i > numat) then
        deallocate (txtatm, txtatm1, labels, atmass)
        deallocate (geoa, c, coord, lopt, na, nb, nc)
        allocate(txtatm(i), txtatm1(i), labels(i), atmass(i))
        allocate(geoa(3,i), c(3,i), coord(3,i),lopt(3,i), na(i), nb(i), nc(i))
      end if
       line = trim(geo_ref_name)
      call upcase(line, len_trim(line))
      if (index(line, '.ARC') /= 0) then
        do
          read(99,"(a)", iostat=ii)line_1
          if (ii /= 0) exit
          if (index(line_1, "HEAT OF FORMATION") > 0) arc_hof_2 = reada(line_1,20)
          if (index(line_1," FINAL GEOMETRY OBTAINED") /= 0) exit
        end do
        if (ii /= 0) rewind(99)
      end if
      line_1 = trim(keywrd)
      if (index(keywrd_quoted,"GEO_DAT") /= 0) then
        if (geo_ref_name == job_fn) then
          i = index(keywrd_quoted," GEO_REF") + 11
          do
            if (keywrd_quoted(i:i) == '"' .or. keywrd_quoted(i:i) == "'") exit
            i = i + 1
          end do
          density = reada(keywrd_quoted, i + 1)
          write(iw,'(/10x,a,f8.3,a)')"A restraining force of",density," kcal/mol/A^2 will be used"
          geoa(:,:numat) = geo(:,:numat)
          ii = numat
          goto 96
        end if
      end if
      i = 0
      do
        i = i + 1
        read(99,"(a)", err = 99, iostat = ii)refkey_ref(i)  !  Dummy read over first three lines
        if (refkey_ref(i)(4:4) == "(" .and. refkey_ref(i)(31:31) == ")" .and. refkey_ref(i)(36:36) == ".") then
          rewind(99)
          exit
        end if
        if (ii /= 0) then
          write(iw,'(10x,a)')"Line 1 of data file defined by GEO_REF: """//trim(refkey_ref(1))//""""
          line = trim(refkey_ref(1))
          call upcase(line, len_trim(line))
          if (index(keywrd,"GEO_DAT") /= 0) then
            call mopend("A data file defined by GEO_REF cannot use GEO_DAT to point to another file.")
            return
          else
            write(line,'(a,i1,a)')"The data file defined by GEO_REF contains only ", i, " lines."
            call mopend(trim(line))
            return
          end if
        end if
        if (refkey_ref(i)(1:1) == "*") i = i - 1
        if (i == 1) then
          keywrd = refkey_ref(1)
          call upcase (keywrd, len_trim(keywrd))
          call l_control("GEO_DAT", len_trim("GEO_DAT"), -1)
          call l_control("GEO_REF", len_trim("GEO_REF"), -1)
          line = trim(keywrd)
          if (index(line, " +") + index(line, "++") /= 0) i = i - 1
        end if
        if (i == 3) exit
      end do
      keywrd = trim(line_1)
      refkey_ref(4:6) = "    NULL"
      goto 97
99    write(iw,*)" File' "//trim(line)//"' is faulty"
      return
97    i = quoted('GEO_REF=')
      do
        if (keywrd_quoted(i:i) == '"' .or. keywrd_quoted(i:i) == "'") exit
        i = i + 1
      end do
      if (keywrd_quoted(i + 1: i + 1) == " ") then
        if (.not. l_0SCF_HTML .and. index(keywrd, "LOCATE-TS") + index(keywrd, "SADDLE") == 0) &
          write(iw,'(/10x,a)')"By default, no restraining force will be used"
        density = 0.d0
      else
        density = reada(keywrd_quoted, i + 1)
        write(iw,'(/10x,a,f8.3,a)')"A restraining force of",density," kcal/mol/A^2 will be used"
      end if
      if (id == 3 .and. abs(density) > 1.d-10 .and. index(keywrd, " 0SCF") == 0) then
        write(line,'(a)')"For solids, the restraining force in GEO_REF must be 0.0."
        call mopend(trim(line))
        return
      end if
      numat_dat = numat
      nat(numat + 1:) = 0
      atom_no(:numat) = nat(:numat)
!
!  Put atoms in the reference data set into the same order as those in the current data set
!
      ii = numat
      j = maxtxt
      call getgeo (99, labels, geoa, c, lopt, na, nb, nc, intern)
      if (natoms == 0) then
        call mopend("NO ATOMS DETECTED IN GEO_REF FILE")
        return
      end if
      if (numat > ii) then
!
! Increase the size of arrays
!
        deallocate (nat)
        allocate (nat(numat + 10))
        nat(:ii) = atom_no(:ii)
      end if
      if (index(keywrd, " OPT ") /= 0) lopt = 1
      deallocate(c)
      maxtxt = j
      if (moperr .and. natoms == 0) then
        call mopend("An error was detected in the reference data set.")
        return
      end if
      if (natoms == -2) then
        j = ir
        ir = 99
        rewind(99)
!
!  Do NOT remove CHAINS - no point in adding atoms that cannot be used
!
        call getpdb(geoa)
        ir = j
        numat_old = numat
      end if
!
! Force geometry to be Cartesian, so that both systems are "clean"
!
      call gmetry (geoa, coord)
      geoa(:,:numat) = coord(:,:numat)
      na(:natoms) = 0
      j = 0
      do i = 1, natoms - id
        if (labels(i) /= 99) then
          j = j + 1
          labels(j) = labels(i)
          txtatm(j) = txtatm(i)
        end if
      end do
      numat = numat - id
      numat_ref = numat
      id = 0
      if (index(keywrd, " 0SCF") == 0 .and. ii /= numat) then
        num = char(ichar("1") +int(log10(ii + 0.05)))
        write(iw,'(/10x,a,i'//num//')')"Number of atoms in """//trim(geo_dat_name)//""" = ", ii
        num = char(ichar("1") +int(log10(numat + 0.05)))
        write(iw,'(/10x,a,i'//num//')')"Number of atoms in """//trim(geo_ref_name)//""" = ", numat
        call mopend("Number of atoms in both systems must be the same, unless keyword ""0SCF"" is present")
        return
      end if
96    txtatm1(:ii) = tmp_txt(:ii)
      i = max(ii, numat, 26)
      allocate (tmp_geoa(3,i), same(i), ok(i), diffs(i))
      diffs = " "
      same = .false.
      ok = .false.
      call l_control("GEO_DAT", len_trim("GEO_DAT"), -1)
      call l_control("GEO_REF", len_trim("GEO_REF"), -1)
      if (index(keywrd," 0SCF") + index(keywrd, " LOCATE-TS") + index(keywrd, " SADDLE") /= 0) then
        if (maxtxt == 0 .or. txtatm(1) == " " .or. txtatm1(1) == " ") then
          if (index(keywrd, " GEO-OK") /= 0) then
            if (txtatm(1) == " " .and. txtatm1(1) == " ") call geochk()
            if (txtatm(1) == " ") then
              txtatm(:numat) = txtatm1(:numat)
            else
              txtatm1(:numat) = txtatm(:numat)
            end if
            maxtxt = 26
          else
            call mopend("In this job, atoms in both data-sets must have PDB-type labels")
            write(iw,'(a)')"Or add ""GEO-OK"" if one of the two data sets has PDB-type "// &
              "labels and the other does not have PDB information."
            return
          end if
        end if
!
!  Before starting, check that the chains are the same.  If they're not, then
!  delete chain information.
!
        do j = 1, ii
          k = ichar(txtatm1(j)(22:22)) - ichar("A") + 1
          if (k > 0) ok(k) = .true.
        end do
        do j = 1, numat
          k = ichar(txtatm(j)(22:22)) - ichar("A") + 1
          if (k > 0) same(k) = .true.
        end do
        do i = 1, 26
          if (same(i) .and. ok(i)) exit
        end do
        if (i == 27) then
          write(iw,'(/10x,a)')" The two systems have different chain letters, so all chain letters changed to ""X"""
          txtatm1(:ii)(22:22) = "X"
          txtatm(:numat)(22:22) = "X"
        end if
        ok = .false.
        same = .false.
!
! Put atoms in the GEO_REF file into the same order as those in the data-set
!
!    GEO_REF arrays:   geoa, txtatm, labels
!    dataset arrays:   geo, txtatm1, nat
!
!  Align atoms in three passes.  In the first pass, all atoms with the same label are aligned.
!  In the second pass, hydrogen atoms that might be mis-labeled within a residue are assigned.
!  In the third pass, any remaining hydrogen atoms are assigned.
!
! If deuterium atoms are used, the PDB name will have "D" instead of "H"
! so convert the "D"'s to "H"'s for the purpose of ordering the atoms.
!
        do i = 1, numat
          if (labels(i) == 1 .and. (index(txtatm(i)(:17), " D") /= 0)) then
            j = index(txtatm(i), " D")
            if (txtatm(i)(j + 1:j + 2) /= "DH") exit
          end if
        end do
        do j = 1, numat
          if (labels(j) == 1 .and. (index(txtatm1(j)(:17), " D") /= 0)) then
            k = index(txtatm1(j), " D")
            if (txtatm1(j)(k + 1:k + 2) /= "DH") exit
          end if
        end do
        if ((i > numat .and. j <= numat) .or. (i <= numat .and. j > numat)) then
          do i = 1, numat
            if (labels(i) == 1) then
              j = index(txtatm(i)(:17), " D")
              if (j /= 0) then
                txtatm(i)(j + 1:j + 1) = "H"
              end if
              j = index(txtatm1(i)(:17), " D")
              if (j /= 0) then
                txtatm1(i)(j + 1:j + 1) = "H"
              end if
            end if
          end do
        end if
!
!   Delete all atoms except those in the backbone
!
        if (index(keywrd, " COMPARE(BACKBONE)") /= 0) then
          j = 0
          do i = 1, numat
            continue
            if (txtatm(i)(13:15) == " N " .or. txtatm(i)(13:15) == " CA" .or. txtatm(i)(13:15) == " C ") then
              j = j + 1
              txtatm(j) = txtatm(i)
              geoa(:,j) = geoa(:,i)
              labels(j) = labels(i)
            end if
          end do
          j = 0
          do i = 1, numat
            if (txtatm1(i)(13:15) == " N " .or. txtatm1(i)(13:15) == " CA" .or. txtatm1(i)(13:15) == " C ") then
              j = j + 1
              txtatm1(j) = txtatm1(i)
              geo(:,j) = geo(:,i)
              nat(j) = nat(i)
            end if
          end do
          numat = j
          ii = numat
!
!   End of delete atoms
!
        end if
        swap = (index(keywrd, " NOSWAP") == 0)
        k = 0
        l = 0
        do j = 1, ii
          j4 = 0
!
!  First pass
!
          do i = 1, numat
            if (same(i)) cycle
            if (txtatm1(j)(12:) == txtatm(i)(12:)) then
              k = k + 1
              tmp_txt(k) = trim(txtatm(i))
              coord(:,k) = geo(:,j)
              tmp_geoa(:,k) = geoa(:,i)
              labels(k) = nat(j)
              same(i) = .true.
              ok(j) = .true.
              j4 = 1
              exit
            end if
          end do
        end do
        let = (index(keywrd, " LOCATE-TS") + index(keywrd, " SADDLE")  /= 0 .or. &
          (index(keywrd," GEO-OK+") /= 0 .and. index(keywrd, " COMPARE") /= 0))
        do j = 1, ii
!
!  Second pass
!
          if(ok(j)) cycle
          do i4 = 1, numat
            if (same(i4) .or. ok(i4)) cycle
            if (txtatm1(j)(14:14) == "H" .and. txtatm(i4)(14:14) == "H" .and. &
              txtatm1(j)(18:) == txtatm(i4)(18:)) then
              if (txtatm1(j)(15:15) == txtatm(i4)(15:15)) then
              k = k + 1
              tmp_txt(k) = txtatm(i4)
              coord(:,k) = geo(:,j)
              tmp_geoa(:,k) = geoa(:,i4)
              labels(k) = nat(j)
              same(i4) = .true.
              ok(j) = .true.
              exit
              end if
            end if
          end do
          if (.not. ok(j) .and. .not. let) then
            l = l + 1
            write(diffs(l)(:40),'(3x,a)')txtatm1(j)
          end if
        end do
        if (let) then
          first = .true.
          do j = 1, ii
!
!  Third pass
!
            if(ok(j)) cycle
            do i4 = 1, numat
              if (same(i4)) cycle
              if (txtatm1(j)(14:14) == "H" .and. txtatm(i4)(14:14) == "H") then
                if (first) then
                  write(iw,'(/16x,a)')"Hydrogen atoms that are on different residues"
                  if (trim(job_fn) == trim(geo_dat_name)) then
                    line = "dataset"
                  else
                    line = "GEO_DAT"
                  end if
                  write(iw,'(/6x,a,14x,a,/)')"Atoms in "//trim(line)//" only","     Atoms in GEO_REF only"
                  first = .false.
                end if
                line = txtatm1(j)
                line(41:) = txtatm(i4)
                write(iw,'(4x,a)')trim(line)
                k = k + 1
                tmp_txt(k) = txtatm(i4)
                coord(:,k) = geo(:,j)
                tmp_geoa(:,k) = geoa(:,i4)
                labels(k) = nat(j)
                same(i4) = .true.
                ok(j) = .true.
                exit
              end if
            end do
            if (.not. ok(j)) then
              l = l + 1
              write(diffs(l)(:40),'(3x,a)')txtatm1(j)
            end if
          end do
          if (.not. first) write(iw,'(/18x,a)')"(Atom labels from GEO_REF will be used)"
        end if
        geo(:,:k) = coord(:,:k)
        ii = k
        jj = 0
        do i = 1, numat
          if (.not. same(i)) then
            jj = jj + 1
            write(diffs(jj)(41:),'(3x,a)')txtatm(i)
          end if
        end do
        nat(:k) = labels(:k)
        do i = 1, k
          atmass(i) = ams(nat(i))
        end do
        txtatm(:k) = tmp_txt(:k)
        if (index(keywrd, " GEO-OK") == 0) then
          geo_dat = "dataset"
        else
          geo_dat = "GEO_DAT"
        end if
        l = max(l, jj)
        if (l > 0) then
          write(iw,'(/28x,a)')"Differences in atoms sets"
          write(iw,'(/10x,a,19x,a,/)')"Atoms in "//geo_dat//" only",    "Atoms in GEO_REF only"
          do i = 1, l
            write(iw,'(i4,a)')i, trim(diffs(i))
          end do
          if (.not. let .and. index(keywrd, " COMPARE") == 0) &
            write(iw,'(/10x,a)')"If ""LET"" is used, then some or all of these differences will be ignored"
        end if
        geoa(:,:k) = tmp_geoa(:,:k)
        tmp_txt(:k) = txtatm(:k)
        if (numat_dat /= numat_ref) then
          if (index(keywrd, " GEO-OK") == 0 .and. .not. l_0SCF_HTML) then
            write(line,'(a, i5, a, i5)')"Number of atoms in data-set:",numat_dat, ", in GEO_REF:", numat_ref
            call mopend(trim(line))
            i = index(keywrd_txt," GEO_DAT") + 9
            if (i > 9) then
              write(iw,'(/10x, a,/10x,a, /)')"Docking can only be done when the number of atoms in ", &
              "the files defined by GEO_DAT and GEO_REF are the same."
            else
             write(iw,'(/5x, a, /)')"Docking can only be done when the number of atoms in the data-set and GEO_REF are the same."
            end if
            i = (k*100)/max(numat,ii)
            write(iw,'(10x,a,i3,a)')"(The two systems have",i, &
            "% of atoms in common. To continue, but using","only those atoms that are common to both systems, add ""GEO-OK"")"
            i = index(keywrd_txt," GEO_DAT") + 9
            if (i > 9) then
              write(iw,'(/10x,a)')"GEO_DAT name: """//trim(geo_dat_name)//""""
            else
              write(iw,'(/10x,a)')"Data-set name: "//trim(geo_dat_name)
            end if
            write(iw,'(/10x,a)')"GEO_REF name: """//trim(geo_ref_name)//""""
            return
          else
            write(iw,'(/10x, a)')"Empirical formulae of data-set and GEO_REF are different"
          end if
!
! First do GEO_REF
! Store the nat for GEO_DAT in nc
! Load the labels for geo_ref into nat
!
          nc(:numat_dat) = nat(:numat_dat)
          nat(:numat) = labels(:numat)
          call empiri()
          write(iw,'(/10x,a,i5)')   "Empirical formula of system in GEO_REF:"//trim(formula(30:))
!
! Now do GEO_DAT
! First store the data for GEO_REF
!
          i = numat
          numat = numat_dat
          nat(:numat) = atom_no(:numat)
          call empiri()
!
! Restore GEO_REF
!
          numat = i
          nat(:numat) = nc(:numat)
          write(iw,'(10x,a,i5)')"Empirical formula of system in "//geo_dat//":"//trim(formula(30:))
          write(iw,'(/10x,a,i5,/)') "Number of atoms common to both systems:",ii
        end if
        natoms = k
        numat = k
        if (pdb_label) txtatm1(:numat) = txtatm(:numat)
        txtatm(:numat) = tmp_txt(:numat)
        coord(:,:numat) = geoa(:,:numat)
      end if
      if (numat /= numat_dat .and. index(keywrd, " LOCATE-TS") /= 0) then
        if (maxtxt == 0 .or. txtatm(1) == " ") then
          if (index(keywrd, " GEO-OK") /= 0) then
            if (maxtxt == 0) then
              maxtxt = 26
              txtatm1(:numat) = txtatm(:numat)
            else
              txtatm(:numat) = txtatm1(:numat)
            end if
          else
            call mopend("In this job, atoms in both data-sets must have PDB-type labels")
            write(iw,'(a)')"Or add ""GEO-OK"" if one of the two data sets has "// &
              "PDB-type labels and the other does not have PDB information."
            return
          end if
        else
          call mopend("ERRORS DETECTED IN PDB LABELS DURING DOCKING MUST BE CORRECTED BEFORE LOCATE-TS CAN BE RUN.")
          write(iw,'(a)')"Or add ""GEO-OK"" if one of the two data sets has"// &
            " PDB-type labels and the other does not have PDB information."
          return
        end if
      end if
      if (moperr) then
        i = index(keywrd_txt," GEO_REF")
        j = index(keywrd(i + 10:),'"') + i + 8
        line = keywrd(i + 10:j)
        write(line,'(a)')"Fault detected in GEO_REF data set: '"//trim(line)//"'"
        write(iw,'(a)')trim(line)
        call mopend(trim(line))
        return
      end if
      call gmetry(geo, coord)
      if (natoms == 0) return
      geo(:,:natoms) = coord(:,:natoms)
      call gmetry(geoa, coord)
      geoa(:,:natoms) = coord(:,:natoms)
      na(:natoms) = 0
      if (index(keywrd, " NOREOR") /= 0) then
        call geo_diff(sum, rms, .true.)
        sum1 = sum/numat
        sum2 = sqrt(rms/numat)
        goto 98
      end if
      if (pdb_label .and. txtatm(1)(:4) == "ATOM" .and. txtatm1(1)(:4) == "ATOM") then
!
!  Put atoms in the reference data set into the same order as those in the current data set
!
        bug = .false.
        any_bug = .false.
!
!  Check for faults in the data-set labels
!
        call compare_txtatm(bug, any_bug)
        bug = .false.
!
!  Check for faults in the geo_ref labels
!
        do i = 1, numat
          if (labels(i) == 1) then
            l = 1
            do j = i + 1, numat
              if (labels(j) == 1) then
                if (txtatm(i)(12:) == txtatm(j)(12:)) then
                  do k = 13, 16
                    if (txtatm(j)(k:k) == " ") then
                      l = l + 1
                      txtatm(j)(k:k) = char(l + ichar("0"))
                      txtatm1(j)(k:k) = char(l + ichar("0"))
                      exit
                    end if
                  end do
                end if
              end if
              if (l == 9) exit
            end do
            if (l > 1) then
              do k = 13, 16
                if (txtatm(i)(k:k) == " ") then
                  txtatm(i)(k:k) = "1"
                  txtatm1(i)(k:k) = "1"
                  exit
                end if
              end do
            end if
          end if
        end do
!
!  At this point, txtatm1 holds the data from the input-dataset.
!                 txtatm holds the data from the geo_ref dataset
!
        k = 0
        do i = 1, numat
          do j = i + 1, numat
            if (txtatm(i)(12:) == txtatm(j)(12:)) exit
          end do
          if (j <= numat) then
            if ( .not. bug) then
              ii = index(keywrd_txt," GEO_REF")
              jj = index(keywrd(ii + 10:),'"') + ii + 8
              line = keywrd(ii + 10:jj)
              write(iw,'(/10x,a,/)')"Atoms in the GEO_REF file '"// &
              trim(geo_ref_name)//"' with the same labels"
            end if
            write(iw,'(10x,a,i6,a,i6,a)')"Atoms", i, " and", j, &
              ";  Labels: ("//txtatm(i)(:maxtxt)//") and ("//txtatm(j)(:maxtxt)//")"
              bug = .true.
              k = k + 1
          end if
        end do
        if (bug) then
          if (index(keywrd, " GEO-OK") == 0) &
            write(iw,'(/10x,a,/)')"Edit GEO_REF data set and re-run."
          any_bug = .true.
          bug = .false.
        end if
!
!  Check for labels in the data-set that are not in geo_ref
!
        do i = 1, numat
          do j = 1, numat
            if (txtatm1(i)(12:) == txtatm(j)(12:)) exit
          end do
          if (j > numat) then
            if ( .not. bug) then
              ii = index(keywrd_txt, "GEO_DAT=")
              if (ii > 0) then

                write(iw,'(/10x,a,/)')"Atoms in the GEO_DAT file '"// &
                trim(geo_dat_name)//"' with no equivalent in the geo_ref file"
              else
                write(iw,'(/10x,a,/)')"Atoms in the data-set file '"// &
                trim(geo_dat_name)//"' with no equivalent in the geo_ref file"
              end if
            end if
            write(iw,'(10x,a,i6,a,i6,a)')"Atom", i, ";  Label: ("//txtatm1(i)(:maxtxt)//")"
              bug = .true.
              any_bug = .true.
              k = k + 1
          end if
        end do
        bug = .false.
        do i = 1, numat
          do j = 1, numat
            if (txtatm(i)(12:) == txtatm1(j)(12:)) exit
          end do
          if (j > numat) then
            if ( .not. bug) then
              write(iw,'(/10x,a,/)')"Atoms in the GEO_REF file '"// &
                & trim(geo_ref_name)//"' with no equivalent in the data-set file"
            end if
            write(iw,'(10x,a,i6,a,i6,a)')"Atom", i, ";  Label: ("//txtatm(i)(:maxtxt)//")"
              bug = .true.
              any_bug = .true.
              k = k + 1
          end if
        end do
        if (any_bug) then
          if (index(keywrd, "GEO-OK") == 0) then
            any_bug = .false.
            if (k == 1) then
              call mopend("Fault detected in atom labels in a GEO_REF run.")
            else
              call mopend("Faults detected in atom labels in a GEO_REF run.")
            end if
            call mopend("(To continue with the current data set, use 'GEO-OK')")
            call web_message(iw,"geo_ref.html")
            return
          else
            if (k == 1) then
              write(iw,'(/10x,a)')"Fault detected in atom labels in a GEO_REF run,"
            else
              write(iw,'(/10x,a)')"Faults detected in atom labels in a GEO_REF run,"
            end if
            if (index(keywrd, " COMPARE") /= 0) then
              write(iw,'(10x,a)')"but because COMPARE was used, corrective action will be taken."
            else
              write(iw,'(10x,a)')"but because GEO-OK was used, corrective action will be taken."
            end if
          end if
        end if
!
!   Search all atoms in the geo_dat_name set to find which atom maps up with the
!   corresponding atom in GEO_REF file.
!
!   (Atoms in the GEO_REF file are never re-arranged or moved.)
!
!  First pass:  Map all atoms that are an exect match
!
        map_atoms_A = -1
        do i = 1, numat
          do j = 1, numat
            if (txtatm(j)(12:) == txtatm1(i)(12:)) exit
          end do
          if (j <= numat) then
            map_atoms_A(j) = i
          end if
        end do
        continue
!
! Second pass: If any atoms are not mapped, then
! find atoms that are NOT mapped, and put them into the map.
!
!   (This code should do nothing unless there is a bug higher up.)
!
        do i = 1, numat
          if (map_atoms_A(i) == -1) then
            do k = 1, numat
              do j = 1, numat
                if (map_atoms_A(j) == k) exit
              end do
              if (j > numat) then
                map_atoms_A(i) = k
                exit
              end if
            end do
            continue
          end if
        end do
        continue
!
!  At this point, the atom "map_atoms_A(i)" in GEO_DAT is replaced with atom i
!
        tmp_txt(:numat) = txtatm1(:numat)
        coord(:,:numat) = geo(:,:numat)
        do i = 1, numat
          j = map_atoms_A(i)
          geo(:,i) = coord(:,j)
          txtatm1(i) = tmp_txt(j)
          atom_no(i) = j
        end do
      else
        i = size(atom_no)
        do i = 1, numat
          atom_no(i) = i
        end do
      end if
      nat(:numat) = labels(:numat)
      coord(:,:numat) = geoa(:,:numat)
!                                                            Add atom numbers at this point
      do i = 1, numat
        write(txtatm1(i)(7:11),'(i5)')i
        write(txtatm (i)(7:11),'(i5)')i
      end do
      rms_min = 1.d6
      if (.not. allocated(same)) allocate(same(numat))
      call geo_diff(sum, rms, .false.)
      toler = 3.d0
      i4 = 1
      j4 = numat
      k4  = 1
      iquit = 0
      xmin = 1.d8
      exists = .true.
      if (index(keywrd, "NOREOR") /= 0) then
        ub = -1
        call geo_diff(sum, rms, .false.)
        sum1 = sum/numat
        sum2 = sqrt(rms/numat)
      else
        ub = 30
      end if
      do ii = 1, ub
        call dock (geo, geoa, sum)
        sum = 0.d0
        do i = 1, numat
          j = j + 1
          sum = sum + sqrt((geo(1,i) - geoa(1,i))**2 + &
                            (geo(2,i) - geoa(2,i))**2 + &
                            (geo(3,i) - geoa(3,i))**2)
        end do
        if (swap .and. ii > 1) then
!
!  Check that atoms are in maximum coincidence.
!  If they are not, then re-arrange geoa
!
          do i = 1, numat
            sum = (geo(1,i) - geoa(1,i))**2 + &
                  (geo(2,i) - geoa(2,i))**2 + &
                  (geo(3,i) - geoa(3,i))**2
            same(i) = (sum < toler)
          end do
          jj = 0
          do i = i4, j4, k4
            if (txtatm1(i)(18:20) == "TRP") cycle ! Add more of these tests, as needed
            sum = (geo(1,i) - geoa(1,i))**2 + &
                  (geo(2,i) - geoa(2,i))**2 + &
                  (geo(3,i) - geoa(3,i))**2
            if (sum > toler) then
              sum = 1.d8
              k = i
              do j = i4, j4, k4
                if (txtatm1(j)(18:20) == "TRP") cycle ! Add more of these tests, as needed
                if ((txtatm1(i)(:6) == txtatm1(j)(:6)    .and.  &
                  txtatm1(i)(22:)   == txtatm1(j)(22:)   .and.  &
                  txtatm1(i)(14:14) == txtatm1(j)(14:14) .or. &
!
!  Always allow water molecules to swap
!
                (txtatm1(i)(18:20) == "HOH" .and. txtatm1(j)(18:20) == "HOH")) .and.  &
                  labels(i) == labels(j) .and. .not. (same(i) .and. same(j))) then
                  num = txtatm1(i)(15:15)
                  if (num == " " .or. (num >= "0" .and. num <= "9")) then
                    num = txtatm1(j)(15:15)
                    if (num == " " .or. (num >= "0" .and. num <= "9")) num = "+"
                  end if
                  if (txtatm1(i)(15:15) /= txtatm1(j)(15:15) .and. num /= "+") cycle
                  sum2 = ((geo(1,j) - geoa(1,j))**2 + &
                          (geo(2,j) - geoa(2,j))**2 + &
                          (geo(3,j) - geoa(3,j))**2)
                  sum3 = ((geo(1,i) - geoa(1,i))**2 + &
                          (geo(2,i) - geoa(2,i))**2 + &
                          (geo(3,i) - geoa(3,i))**2)
                  dum1 = ((geo(1,i) - geoa(1,j))**2 + &
                          (geo(2,i) - geoa(2,j))**2 + &
                          (geo(3,i) - geoa(3,j))**2)
                  dum2 = ((geoa(1,i) - geo(1,j))**2 + &
                          (geoa(2,i) - geo(2,j))**2 + &
                          (geoa(3,i) - geo(3,j))**2)
!
!  If swapping the atoms around will improve the overlap, then do so.
!
!  dum1 = distance contribution from atom i if atom i and j were swapped around.
!  dum2 = distance contribution from atom j if atom i and j were swapped around.
!  sum2 = distance contribution from atom j if the atoms were not swapped around
!  sum3 = distance contribution from atom i if the atoms were not swapped around
!  sum  = smallest distance contribution, i.e., best choice for swapping
!
                  if (dum1 + dum2 < sum - 0.1d0 .and. dum1 + dum2 < sum2 + sum3 - 0.1d0) then
                    sum = dum1 + dum2
                    sum4 = sum2 + sum3
                    k = j
                  end if
                end if
              end do
              if (i /= k) then
                jj = jj + 1
                if (.false.) then
                if (exists) then
                  exists = .false.
                  write(iw,'(" ",20("----"))')
                  if (pdb_label) then
                    write(iw,'(/4x,a,//8x,a,/)')"  List of atoms in GEO_REF that were swapped in order to maximize overlap", &
                    "      Atom Label          and         Atom Label         Difference to GEO_REF"
                  else
                    write(iw,'(/4x,a,//8x,a,/)')"  List of atoms that were swapped in order to maximize overlap", &
                    "  Atom name   Atom No. and Atom No.  Difference to GEO_REF"
                  end if
                end if
                if (pdb_label) then
                  write(iw,'(i5, a, 2x, a, f12.2)')jj,  " ("//txtatm1(i)(:maxtxt)//")",  &
                  " ("//txtatm1(k)(:maxtxt)//")", sqrt(sum) - sqrt(sum4)
                else
                    write(iw,'(i5,1x, a, i9,  i13, 1x, f12.2)')jj, atom_names(nat(i)), atom_no(i),  atom_no(k), &
                  sqrt(sum) - sqrt(sum4)
                end if
                end if
                do j = 1,3
                  sum = geoa(j,i)
                  geoa(j,i) = geoa(j,k)
                  geoa(j,k) = sum
                end do
                line_1 = txtatm1(i)
                txtatm1(i) = txtatm1(k)
                txtatm1(k) = trim(line_1)
                j = atom_no(i)
                atom_no(i) = atom_no(k)
                atom_no(k) = j
              end if
              same(i) = .true.
            end if
          end do
        end if
        if (.not. exists) write(iw,*)" "
        k = 0
        do i = 1, numat
          do l = 1,3
            k = k + 1
            xparam(k) = geo(l,i)
          end do
        end do
          call geo_diff(sum, rms, .false.)
          sum1 = sum/numat
          sum2 = sqrt(rms/numat)
        if (abs(rms_min - rms) < 1.d-4 .and. ii > 2) exit
        if (xmin > rms) then
          xmin = rms
          iquit = 0
        else
          iquit = iquit + 1
          if (iquit > 3) exit
        end if
        rms_min = rms
        toler = max(1.d0, 0.7d0*toler)
        k4 = i4
        i4 = j4
        j4 = k4
        if (i4 == 1) then
          k4 = 1
        else
          k4 = -1
        end if
      end do
!
!  Docking complete. Now un-swap GEO_REF
!
      do i = 1, numat
        do k = i, numat
          if (txtatm1(i) == txtatm(k)) then
            if (i == k) exit
            do j = 1,3
              sum = geoa(j,i)
              geoa(j,i) = geoa(j,k)
              geoa(j,k) = sum
              sum = geo(j,i)
              geo(j,i) = geo(j,k)
              geo(j,k) = sum
            end do
            txtatm1(k) = txtatm1(i)
          end if
        end do
      end do
      txtatm1(:numat) = txtatm(:numat)
      if (sum1 < 1.d-10 .and. index(keywrd, " LOCATE-TS") /= 0) then
        write(line,'(a)')" The data-set and GEO_REF geometries are the same.  Correct the fault and re-run."
        call mopend (trim(line))
        return
      end if
      write(iw,'(/24x,a)')"After docking"
      call geo_diff(sum, rms, .true.)
98    if (trim(job_fn) == trim(geo_dat_name)) then
        line = "dataset"
      else
        line = "GEO_DAT"
      end if
      i = ichar("5") + max(0, int(log10(sum1*numat + 1.d-20)))
      if (i > 57) then
        num(1:1) = "1"
        num(2:2) = char(i - 10)
      else
        num(1:1) = char(i)
        num(2:2) = " "
      end if
      write(iw,'(/3x,a,f'//trim(num)//'.3,a,f8.4,a,f8.4,a)') &
        "Difference between "//trim(line)//" and GEO_REF: ", sum1*numat, &
        " = total,", sum1, " = Average,", sum2," = RMS, in Angstroms"
      if (index(keywrd," 0SCF") /= 0 .and. sum1 < 1.d-6) then
        call mopend("The two systems are identical, so the job will stop here.")
        return
      end if
      close(99)
      if (index(keywrd," NOCOM") /= 0) ncomments = 0
      if (index(keywrd, " COMPARE") /= 0) then
!
! Write out the reference geometry file needed for comparing the two geometries using JSmol.
!
        line_1 = geo_ref_name
        call upcase(line_1, len_trim(line_1))
        i = len_trim(line_1)
        j = max(1,i - 7)
        if (line_1(i - 3: i) == ".TXT" .and. line_1(j:j) == ".") &
        geo_ref_name = geo_ref_name(:len_trim(geo_ref_name) - 4)
        if (index(keywrd, " COMPARE") /= 0) then
          if (index(trim(geo_ref_name), "'") /= 0) then
            call mopend("When keyword ""COMPARE"" is used the GEO_REF filename"// &
             " must not contain a single quotation mark, i.e., a ""'"" mark")
            return
          end if
          if (index(trim(geo_dat_name), "'") /= 0) then
            if (index(keywrd, "GEO_DAT") /= 0) then
              call mopend("When keyword ""COMPARE"" is used the GEO_REF filename"// &
               " must not contain a single quotation mark, i.e., a ""'"" mark")
            else
              call mopend("When keyword ""COMPARE"" is used the data-set filename"// &
               " must not contain a single quotation mark, i.e., a ""'"" mark")
            end if
            return
          end if
          do i = len_trim(job_fn), 1, -1
            if (job_fn(i:i) == "/" .or. job_fn(i:i) == backslash) exit
          end do
          do j = len_trim(job_fn) - 4, 2, -1
            if (job_fn(j:j) /= " ") exit
          end do
          line_1 =  job_fn(i + 1:j)//"_"
        else
          line_1 = " "
        end if
        line = trim(line_1)//geo_ref_name(:len_trim(geo_ref_name) - 3)//"pdb"
        call add_path(line)
        line_2 = trim(line)
        open(unit = 99, file = trim(line), iostat = i)
        call l_control("HTML", len_trim("HTML"), -1)
        coord(:,:numat) = geo(:,:numat)
        geo(:,:numat) = geoa(:,:numat)
        call pdbout(99)
        geo(:,:numat) = coord(:,:numat)
        call l_control("HTML", len_trim("HTML"), 1)
        close (99)
      end if
      coord(:,:natoms) = geo(:,:natoms)   ! Store input geometry
      if (index(keywrd, " LOCATE-TS") == 0) then
        do i = 1, nvar
          xparam(i) = geoa(loc(2,i), loc(1,i))
        end do
      else
        if (.not. allocated(lopt_store)) then
!
! Store lopt for input data-set or GEO_DAT in case it will be needed in LOCATE_TS
!
          allocate(lopt_store(3,natoms))
          lopt_store = lopt(:3,:natoms)
        end if
        nvar = 0
        do i = 1, natoms
          do j = 1,3
            if (lopt_store(j,i) == 1) then
              nvar = nvar + 1
              xparam(nvar) = geoa(j,i)
            end if
          end do
        end do
      end if
      geo(:,:natoms) = geoa(:,:natoms)
!
! The following "big_swap" saves the geometry and xparam for the rotated input geometry
!
      call big_swap(0,2)
      geo(:,:natoms) = coord(:,:natoms)   ! Restore input geometry
      do i = 1, nvar
        xparam(i) = geo(loc(2,i), loc(1,i))
      end do
      if (index(keywrd, " COMPARE") /= 0) then
        call wrt_diffs()
        if (abs(arc_hof_1) + abs(arc_hof_2) > 1.d-4) write(iw,*)
        if (trim(job_fn) == trim(geo_dat_name)) then
          line = "dataset"
        else
          line = "GEO_DAT"
        end if
        if (abs(arc_hof_1) > 1.d-4) &
          write(iw,"(10x,a,f12.3,a)")  "Heat of formation of "//trim(line)//" system:",arc_hof_1," Kcal/mol"
        if (abs(arc_hof_2) > 1.d-4) &
          write(iw,"(10x,a,f12.3,a)")  "Heat of formation of GEO_REF system:",arc_hof_2," Kcal/mol"
        if (abs(arc_hof_1) > 1.d-4 .and. abs(arc_hof_2) > 1.d-4) &
          write(iw,"(40x,a,f12.3,a)")  "Diff.:", arc_hof_1 - arc_hof_2," Kcal/mol"
        call analyze_h_bonds()
      end if
      do ii = 1, 6
        line = " "//trim(refkey(ii))
        call upcase(line, len_trim(line))
        i = index(line," GEO_REF")
        if (i /= 0) exit
      end do
      inquire(unit = 99, opened = opend, name = line)
      if (index(keywrd, " TS ") /= 0 .and. index(keywrd, " LOCATE-TS") == 0) then
        write(iw,'(/,a)')"    The average of the supplied and reference geometry will be written to:"
        write(iw,'(a)')"'"//trim(line)//"'"
        geo(:,:natoms) = 0.5d0*(geoa(:,:natoms) + geo(:,:natoms))
        call geout(99)
      else
        if (index(keywrd, " COMPARE") /= 0) then
          do i = len_trim(output_fn), 1, -1
            if (output_fn(i:i) == "/" .or. output_fn(i:i) == backslash) exit
          end do
          if (i > 0) then
            write(iw,'(//,a)')"    The following files will be written to """//output_fn(:i)//""":"
          else
            write(iw,'(//,a)')"    The following files will be written:"
          end if
          line_1 = job_fn
          call upcase(line_1, len_trim(line_1))
          j = len_trim(line_1)
          if (j > 7) then
            if (line_1(j - 3: j) == ".TXT" .and. line_1(j - 7:j - 7) == ".") &
            job_fn = job_fn(:len_trim(job_fn) - 4)
          end if
          line = job_fn(i + 1:len_trim(job_fn) - 3)//"html"
          write(iw,'(14x, a)')"'"//trim(line)//"'"
!
          line_1 = geo_dat_name
          call upcase(line_1, len_trim(line_1))
          i = len_trim(line_1)
          j = max(1,i - 7)
          if (line_1(i - 3: i) == ".TXT" .and. line_1(j:j) == ".") &
          geo_dat_name = geo_dat_name(:len_trim(geo_dat_name) - 4)
          line = geo_dat_name(:len_trim(geo_dat_name) - 3)//"pdb"
          if (index(keywrd, " COMPARE") /= 0) then
            write(iw,'(a)')"From GEO_REF: '"//trim(line_2)//"'"
          else
            write(iw,'(a)')"'"//trim(line_2)//"'"
          end if
!
          line = geo_ref_name(:len_trim(geo_ref_name) - 3)//"pdb"
        end if
        line = refkey_ref(1)
        do i = 1,3
          line = refkey(i)
          refkey(i) = refkey_ref(i)
          refkey_ref(i) = trim(line)
        end do
        refkey_ref(4) = koment
        refkey_ref(5) = title
        if (l_0SCF_HTML) then
!
! Write out the input geometry file needed for comparing the two geometries using JSmol.
!
            close (99)
            if (sum2 > 9.9995d0) then
              num = "6"
            else
              num = "5"
            end if
            write(line,'(a,f'//num//'.3,a)')"RMS_DIFF=",sum2
            call l_control(trim(line), len_trim(line), 1)
            write(line,'(a,f12.3,a)')"DIFF=",sum1*numat
            call l_control(trim(line), len_trim(line), 1)
            if (index(keywrd, " COMPARE") /= 0) then
              do i = len_trim(job_fn), 1, -1
                if (job_fn(i:i) == "/" .or. job_fn(i:i) == backslash) exit
              end do
              do j = len_trim(job_fn) - 4, 2, -1
                if (job_fn(j:j) /= " ") exit
              end do
              line_1 =  job_fn(i + 1:j)//"_"
            else
              line_1 = " "
            end if
            if (index(keywrd, " COMPARE") /= 0 .and. &
              geo_dat_name(:len_trim(geo_dat_name) -3) == geo_ref_name(:len_trim(geo_ref_name) -3)) then
              line = trim(line_1)//geo_dat_name(:len_trim(geo_dat_name) - 4)//"_a.pdb"
            else
              line = trim(line_1)//geo_dat_name(:len_trim(geo_dat_name) - 4)//".pdb"
            end if
            coord(:,:numat) = geo(:,:numat)
            call add_path(line)
            open(unit = 99, file = trim(line), iostat = i)
            if (index(keywrd, " COMPARE") /= 0) then
              write(iw,'(a)')"From GEO_DAT: '"//trim(line)//"'"
            else
              write(iw,'(a)')"'"//trim(line)//"'"
            end if
            call pdbout(99)
        end if
        do i = 1,3
          line = refkey_ref(i)
          refkey_ref(i) = refkey(i)
          refkey(i) = trim(line)
        end do
        koment = trim(refkey_ref(4))
        title = trim(refkey_ref(5))
        close(99)
      end if
      deallocate(same)
      if (index(keywrd, " LOCATE-TS") == 0) then
        if (index(keywrd," 0SCF") + index(keywrd, " TS ") /= 0) then
!
! If a RESTART, generate an ARC file
!
          if (index(keywrd, " RESTART") /= 0) then
            inquire(unit=iarc, opened=opend)
            if (opend) close(iarc)
            archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"arc"
            open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
            rewind iarc
            call geout (iarc)
          end if
          call mopend("GEO_REF with 0SCF: Job complete")
          return
        end if
      end if
      geo(:,:natoms) = coord(:,:natoms)
      coord(:,:natoms) = geoa(:,:natoms)
      return
  end subroutine geo_ref
  subroutine compare_txtatm(bug, any_bug)
    use molkst_C, only : numat, maxtxt, keywrd, line, keywrd_txt
    use common_arrays_C, only : txtatm, txtatm1, nat
    use chanel_C, only : iw, job_fn
    implicit none
    logical, intent (inout) :: bug, any_bug
    integer :: i, j, ii,jj
    if (index(keywrd, "GEO-OK") /= 0) return
    if (index(keywrd, "GEO_REF") == 0) return
    bug = .false.
    if (maxtxt /= 27) return
    do i = 1, numat
      do j = 1, 2
        if (txtatm(i)(20:20) /= " ") exit
        txtatm(i)(18:20) = " "//txtatm(i)(18:19)
      end do
      do j = 1, 2
        if (txtatm1(i)(20:20) /= " ") exit
        txtatm1(i)(18:20) = " "//txtatm1(i)(18:19)
      end do
    end do
    do i = 1, numat
      do j = i + 1, numat
        if (txtatm1(i)(12:) == txtatm1(j)(12:)) exit
      end do
      if (j <= numat .and. nat(i) /= 1) then
        if ( .not. bug) then
          ii = index(keywrd_txt, "GEO_DAT=")
         if (ii > 0) then
            jj = index(keywrd(ii + 9:),'"') + ii + 7
            line = keywrd(ii + 9:jj)
            write(iw,'(/10x,a,/)')"Atoms in the GEO_DAT file '"// &
            trim(line)//"' with the same labels"
          else
            write(iw,'(/10x,a,/)')"Atoms in the data-set file '"// &
            trim(job_fn)//"' with the same labels"
          end if
          write(iw,'(10x,a,i6,a,i6,a)')"Atoms", i, " and", j, &
          ";  Labels: ("//txtatm1(i)//") and ("//txtatm1(j)//")"
        end if
        bug = .true.
      end if
    end do
    if (bug) then
      call mopend("Error in data detected while using GEO_REF")
      write(iw,"(5x,a)")"(To continue with the current data set, use 'GEO-OK')"
      any_bug = .true.
    end if
    return
  end subroutine compare_txtatm
