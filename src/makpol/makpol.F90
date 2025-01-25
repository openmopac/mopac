! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

program makpol
!***********************************************************************
!
!   Program MAKPOL- A utility for MOPAC.  MAKPOL can read in a PDB file
!                   with suffix ".ent" for example Urea.ent, or a MOPAC
!                   <Filename>.DAT file which contains the MERS keyword.
!                   It will produce a MOPAC data-file named <filename>.mop
!                   that is suitable for running solid-state systems.
!
!                   If keyword "TV" is not present, then a second file
!                   named Make_<filename>.dat will also be created.  
!                   This file is the same as <filename>.mop.  If its
!                   keyword MERS is edited, and MAKPOL run using
!                   Make_<filename>.dat, then a new <filename>.mop will
!                   be generated. 
!
!***********************************************************************
    use common_systm, only : iw, l_int, ir, ifiles_1, ndep, natoms, id, numat, line, l_use_tv
    use common_symult, only : depmul
    use common_keywrd, only : keywrd
    use common_jobnam, only : jobnam
    use common_geosym, only : locpar, idepfn, locdep
    use common_elemts, only : elemnt
    use common_titles, only : koment, title
    use common_ctvec, only : tvec
    use common_sizes, only : numatm
    implicit none
    integer :: i, j, k, l, ll, m, n, nmax, mers(3), store_mers(3), loc(3*numatm), &
      labels(numatm), na(numatm), nb(numatm), nc(numatm), N_Tv(3), M_Tv(3), Tv_atoms(3), &
      store_labels(numatm)
    double precision :: r, rmers(3), coord(3,numatm), geo(3,numatm), xtoc(3,3), &
      vec(3), Tv_n(3), sum, angle, zero(3) = 0.d0, store_tvec(3,3), store_coord(3,numatm)    
    character*1 :: n1, n2, n3
    logical :: exists, debug, pdb_file, l_cell, bug
    double precision, external :: dist, snapth, reada
    call getdat
!
!  Subroutines GETDAT, GETTXT, and GETGEO are from MOPAC. Only minor edits have been made
!  so that they will work in MAKPOL
!
    call gettxt
    if (index(keywrd, "HEADER") /= 0 .and. index(keywrd, " CSD") /= 0)  keywrd = " "
    line = "keywords for MAKPOL.txt"
    call add_path(line)
    inquire (file=trim(line), exist = exists)
    if (exists) then
      open(7, file = trim(line))
      read(7,'(a)', iostat = i) keywrd(len_trim(keywrd) + 2:)
      call upcase(keywrd, len_trim(keywrd))
      debug = (index(keywrd, "DEBUG") /= 0)
      if (debug) then
        write(iw,'(/10x,a)')"Keywords supplied by data-set and ""keywords for MAKPOL.txt"""
        write(iw,'(/10x,a)')trim(keywrd)
      end if
      line = trim(keywrd)
!
!  Delete duplicate keywords
!
      i = 0
      do 
        i = i + 1
        if (i >= len_trim(keywrd)) exit
        if (keywrd(i:i) /= " ") then
          j = i
          do 
            if (keywrd(j:j + 2) == "TV=" .or. keywrd(j:j + 3) == "CELL") then ! jump over these keywords
              do
                j = j + 1
                if (keywrd(j:j) == ")") exit
              end do
            end if
            j = j + 1
            if ((keywrd(j:j) < "A" .or. keywrd(j:j) > "Z") .and. &
                (keywrd(j:j) < "0" .or. keywrd(j:j) > "9") .and. &
                (keywrd(j:j) /= "-" .and. keywrd(j:j) /= "_") .and. keywrd(j:j) /= "." ) exit
          end do
          do
            k = index(keywrd(i + 2:), keywrd(i:j - 1))
            if (k == 0) exit
            do m = k + j + 1, 100
              if (keywrd(m:m) == " ") exit
            end do
            keywrd(k + i + 1: m) = " "
          end do
          do i = j, 100
            if (keywrd(i:i) == " ") exit            
          end do
        end if
      end do
  !  else
  !    write(iw,'("*",/,a,/,a,/,"*")')"* This message will not be printed if the file ""keywords for MAKPOL.txt""", &
  !    "* is present and is blank or contains keywords"
    end if
    if (keywrd(1:1) /= " ") keywrd = " "//trim(keywrd)
    if (koment == " ") then
      koment = jobnam
      do i = 1, len_trim(koment)
        if (koment(i:i) == "_") koment(i:i) = " "
      end do
    end if
!
    i = Index (" "//keywrd, " INT") 
    l_int = (i > 0)
    if (l_int) then
      keywrd (i:i + 3) = " "
    end if
    i = index(keywrd, " TV ")
    if (i > 0) then
      l_use_tv = .false.
      keywrd = keywrd(:i)//trim(keywrd(i + 4:))
    end if    
    call l_control("XYZ", len("XYZ"), -1)
    l_cell = (index(keywrd, " CELL") /= 0)
!
!  Read in the geometry in either MOPAC or PDB format
!
    call getgeo (ir, labels, geo, coord, loc, na, nb, nc)    
    if (natoms == -2) then
      rewind (ir)
      call getpdb(ir, geo, coord, na, nb, nc, labels, iw)
      rewind (ir)
      if (debug) then
       write(iw,'(/10x,a,/)')"Geometry read in from PDB file (atom number, label, coordinates)"      
        do i = 1, numat
          write (iw, "(2I6,3F18.12,3I4)") i, labels (i), (coord(k, i), k = 1, 3)
        end do
      end if
!
!  Assume geometry is in PDB format - look for translation vectors
!  and put them in keywrd as a keyword "TV=(a b c)"
!
      do
        read(ir,'(a)', iostat = i) line
        if (i /= 0) exit
        if (line(:5) == "CRYST") then
          i = 7
          do 
            if (line(i:i) /= " ") exit
            i = i + 1
          end do
          keywrd = trim(keywrd)//" CELL("//line(i:55)//")"
         if (l_cell) then
           call cell(coord, tvec)
           labels(numat + 1:numat + 3) = 99
           numat = natoms + 3
           id = 3
         end if
          l_cell = .true.
          exit
        end if
      end do
      koment = " "
      title = " "
      pdb_file = .true.
    else
!
! Geometry is in MOPAC format - look for atomic numbers 2, 93, or 107.  
! Use the coordinates of TV atoms as translation vectors, this is done in gmetrn,
! and the translation vectors are stored in the array tvec.
!
       if (index(keywrd, " XYZ") /= 0) natoms = natoms + 3
       numat = natoms
       call gmetrn (geo, coord, labels, na, nb, nc)
       j = 0
       do i = 1, numat   ! Delete all dummy atoms
!
!  Atoms that were translation vectors are given atomic number 98 by gmetrn,
!  at this point, these atoms are converted into atoms of the same type as atom 1
!  This is so that the atom types from a MOPAC data set and from PDB files are the same.
!
         if (labels(i) == 98) labels(i) = labels(1)
         if (labels(i) /= 99) then
           j = j + 1
           labels(j) = labels(i)
         end if
       end do
       if (l_cell) then
         call cell(coord, tvec)
         labels(numat + 1:numat + 3) = 99
         numat = natoms + 3
         id = 3
       end if
       geo(:,:numat) = coord(:,:numat)
!
!  The geometry is now in Cartesian coordinates, in geo and coord,
!  and the translation vectors are in tvec
!
       na = 0
       l = 0
       if (id > 0) then
         do i = 1, id
           sum = sqrt(tvec(1,i)**2 + tvec(2,i)**2 + tvec(3,i)**2)
           l = l + 1
           if (l == 1) write(keywrd(len_trim(keywrd) + 2:),'(a)')"TV=("
           write(keywrd(len_trim(keywrd) + 1:),'(f8.4)')sum
         end do
       else          
         do i = 2, natoms
          if (labels(i) == 93 .or. labels(i) == 2) labels(i) = 107
          if (labels(i) == 107) then
            sum = dist(coord, 1, i, tvec, xtoc, vec)
            l = l + 1
            if (l == 1) write(keywrd(len_trim(keywrd) + 2:),'(a)')"TV=("
            write(keywrd(len_trim(keywrd) + 1:),'(f8.4)')sum
          end if
         end do
       end if
       if (l /= 0) write(keywrd(len_trim(keywrd) + 1:),'(a)')")"
       pdb_file = .false.
    end if
!
!                                                                         All data now read in
!
! Geometry has been read in, and the translation vectors are known 
! and currently stored in the keyword.
!
!  Clean up keywrd (some keywords might have been added by subroutine "getgeo")
!
    call l_control("CHAINS", len("CHAINS"), -1)
    call l_control("START_RES", len("START_RES"), -1)
    if (labels(1) == 107) then
      write(iw,*)" Atom 1 must not be a translation vector"
      stop
    end if
    if (l_cell) then
         call cell(coord, tvec)
         do id = 1, 3
           Tv_n(id) = sqrt(tvec(1,id)**2 + tvec(2,id)**2 + tvec(3,id)**2)
           N_Tv(id) = 1 + int(11.d0/Tv_n(id))
           M_Tv(id) = 1
         end do
         numat = natoms + 3
         id = 3
         k = 3
    end if
    store_tvec = tvec
    store_coord(:,:numat) = coord(:,:numat)
    store_labels(:numat) = labels(:numat)
    i = index(keywrd, " TV=(")
    if (i /= 0 ) then
!
! Translation distances have been found in keywrd.
!
      Tv_n = 0.d0
      j = index(keywrd(i:), ")") + i
      l = 0
      do k = i, j
        if (keywrd(k:k) == ".") l = l + 1
      end do
      id = l
      k = index(keywrd(i:), ".") + i - 4
      do m = 1, id
        Tv_n(m) = reada(keywrd, k)
        if (m == id) exit
        k = index(keywrd(k + 5:), ".") + k + 1
      end do  
!
! Use Tv_n to work out tvec
!
      zero(:) = geo(:,1)
      do i = 1, natoms
        geo(:,i) = geo(:,i) - zero
      end do
      mers = 0
      do i = 1, id
        do j = 1, numat
          if (labels(j) /= labels(1)) cycle
          do k = 1, i - 1
            if (mers(k) == j) exit
          end do
          if (k < i) cycle
          sum = dist(geo, 1, j, tvec, xtoc, vec) 
          if (abs(sum - Tv_n(i)) < 0.01d0) then
            store_tvec(:,i) = geo(:,j)
            mers(i) = j
            exit
          end if
        end do
      end do
      tvec = store_tvec 
      if (natoms  == numat) then
        do i = 1, id
          geo(:,i + numat) = tvec(:,i)
        end do
        numat = numat + id
      end if

      if (debug) then
        write(iw,'(/5x,a)')"Data at the start of building the unit cell"
        write(iw,'(/5x,a, i4)')"Number of atoms:",numat
        write(iw,'(5x,a, i4)')"Number of Tv:   ",id
        write(iw,'(5x,a)')"KEYWORD: """//trim(keywrd)//""""
        write(iw,'(5x,a,3f9.4)')"Tv_n:           ",Tv_n(:id)
        write(iw,'(/15x,a,9F9.4)')"TVEC           "
        write(iw,'(2x,3f9.4)')tvec(:,:id)
        write(iw,'(/10x,a,/)')"Geometry (atom number, label, coordinates"      
        do l = 1, numat
          write (iw, "(2I6,3F18.12,3I4)") l, labels (l), (geo(k, l), k = 1, 3)
        end do
      end if
      call l_control("TV=", len("TV="), -1)
      if (l_cell) then
        do k = 1, 3
          N_Tv(k) = 1 + int(11.d0/Tv_n(k))
        end do
        k = 3
      else      
!
!  Now work out the size of the smallest translation vector lengths.
!  An earlier run might have made a multiple unit cell.
!  Store the new translation vectors in the array "tvec"
!
! Unconditionally set atom 1 at the origin
!
        zero(:) = geo(:,1)
        do i = 1, natoms
          geo(:,i) = geo(:,i) - zero
        end do
        k = 0
        loop: do i = 2, numat
          if (labels(i) /= labels(1) .and. (pdb_file .or. labels(i) /= 107)) cycle
          sum = dist(geo, 1, i, tvec, xtoc, vec) 
          do j = 1, id
!
! By default, the primitive unit cell is worked out, 
! if the Tv keyword is NOT supplied in the keywords in the data set 
!
            n = int(Tv_n(j)/sum) + 2
            if (.not. l_use_tv) n = 1
            do ll = n, 1, -1
              if (abs(sum - Tv_n(j)/ll) < 0.002d0) exit
            end do          
            if (ll > 1) then
!
!   Atom "i" is at a distance Tv_n(j)/ll from atom 1.
!
              if (ll > 1) then
!   Check that there is a second atom at the same fractional distance
!
                do l = 2, numat
                  if (l /= i .and. labels(l) /= 107) exit
                end do
                do m = 1, numat
                  if (labels(m) /= labels(l) .and. labels(m) /= 107) cycle              
                  sum = dist(geo, l, m, tvec, xtoc, vec)  
                  if (abs(sum - Tv_n(j)/ll) < 0.002d0) exit
                end do     
                if (m > numat) cycle
              end if
              do m = 2, numat
                if (m == i .or. (labels(m) /= labels(1) .and. labels(m) /= 107)) cycle
                sum = dist(geo, 1, m, tvec, xtoc, vec)  
                n = int(Tv_n(j)/sum) + 2
                if (.not. l_use_tv) n = 1
                do l = n, 1, -1
                  if (abs(sum - Tv_n(j)/l) < 0.002d0) then
                    call bangle (geo, i, 1, m, angle)
                    if (angle > 3.1d0) angle = angle - 3.1415926d0
                    if (angle > 0.001d0) then
                      continue
                      cycle
                    end if 
!
!   Atom "m" is at a distance Tv_n(j)/l from atom 1 and makes a straight line with atom "i".
!
                    if (ll > l) then
                      labels(i) = 107
                      labels(m) = 99
                      if (debug) write(iw,'(5x,a, i4, 5x, a, i4)')"Old Tv atom:", m, "New Tv atom:", i
                    else if (l > ll) then
                      labels(m) = 107
                      labels(i) = 99
                      if (debug) write(iw,'(5x,a, i4, 5x, a, i4)')"Old Tv atom:", i, "New Tv atom:", m
                    else
                      labels(m) = 107  
                      if (debug) write(iw,'(5x,a, i4, 5x, a, i4)')"Atom", m, "identified as a Tv"
                    end if
                    k = k + 1
!
! N_Tv holds the number of translation vectors that will be needed for a normal
! calculation.  If keyword MERS was provided, then information in N_Tv will be ignored.
!
                    N_Tv(k) = 1 + int(11.d0/Tv_n(j))
!
!  Tv_n(j) represents "l" or "ll" fundamental translation vector distances
!  Scale it back to exactly one fundamental translation vector distance.
!
                    Tv_n(j) = -Tv_n(j)/max(l, ll)
                    M_Tv(k) = max(l, ll)
                    cycle loop
                  end if
                end do
              end do
              if (Tv_n(j) > 0.1d0) then
!
! Atom "i" is at a distance of a Tv, but there is no other atom to verify that it is the
! correct atom, however the possibility of a mistake is low, so use it.
!
! If this possibility is a cause for concern, extend this test to check the angles or 
! check that there are no other atoms that are at a distance of Tv_n(j) from atom 1.
! Or, best option, use a larger PDB file, by selecting a 3 by 3 by 3 packing.
!
                labels(i) = 107
                k = k + 1
                N_Tv(k) = 1 + int(11.d0/Tv_n(j))
                M_Tv(k) = max(l, ll)
                Tv_atoms(k) = i
                if (debug) write(iw,'(5x,a, i4, 5x, a, i4)')"Atom", i, "identified as a Tv"
!
!  Tv_n(j) represents "l" fundamental translation vectors
!
                Tv_n(j) = -Tv_n(j)/max(l, ll)
                 cycle loop
              end if
            end if
          end do 
        end do loop
      end if
      i = 0
      if (id == 3) then
        i = 0
        do j = 1, 3
          if (Tv_n(j) < -0.001d0) i = i + 1
        end do
      end if
      if (.not. l_cell .and. id == 3 .and. i == 0) then
!
!  No Tv were found, so,
!  and reset the original Tv labels back to 107.
!
        do i = natoms + 1, natoms + 3
          labels(i) = 107
        end do
        natoms = natoms + 3
        do j = 1, id
          do k = 1, 3
            N_Tv(k) = 1 + int(11.d0/Tv_n(j))
            M_Tv(k) = 1
          end do
        end do
        k = 3
      end if      
      Tv_n = abs(Tv_n)
    else if (id /= 3) then
      write(iw,'(/10x,a)')"Translation vector distances not found"
      stop
    end if
    if (debug) then
        write(iw,'(/5x,a)')"Data after constructing Tv_n and N_Tv"
        write(iw,'(/5x,a,3f9.4)')"Tv_n:           ",Tv_n(:id)
        write(iw,'(/5x,a,3i9)')"N_Tv:           ",N_Tv(:id)
    end if  

!
!   Tv and the number of unit cells in each direction, N_Tv, are known.
!   tvec is still not known.
!
!  Write out the number of unit cells in each direction as a keyword MERS=(n,n,n)
!
    if (index(keywrd, " MERS") == 0) then
      if (k == 1) then
        write(keywrd(len_trim(keywrd) + 1:),'(a,i1,")")')" MERS=(", N_Tv(1)
      elseif (k == 2) then
        write(keywrd(len_trim(keywrd) + 1:),'(a,(i1,","),i1,")")')" MERS=(", N_Tv(1), N_Tv(2)
      elseif (k == 3) then
        write(keywrd(len_trim(keywrd) + 1:),'(a,2(i1,","),i1,")")')" MERS=(", N_Tv
      end if
    end if  
    if (Index (keywrd, " SYMM") + Index (keywrd, " SYM ") /= 0) then
      ifiles_1 = 0
      call getsym (locpar, idepfn, locdep, depmul, na)
      if (ndep /= 0) then
        call symtnn (geo, na)
      end if
    end if
    natoms = numat
    if (natoms == 0) then
      write (iw,*) " The system appears to have no atoms!"
      stop
    end if
    close (ir, status = "DELETE")
    numat = natoms
    if (Index (keywrd, " SNAP") /= 0) then
!
!   If any angles are near to important angles (such as 109.47...)
!   snap the angle to the exact angle
!
      do i = 1, numat
        if (na(i) /= 0) then
          geo(2, i) = snapth (geo(2, i))
          geo(3, i) = snapth (geo(3, i))
        end if
      end do
    end if
!
!  Protect any dummy atoms
!
    do i = 1, numat
      if (labels(i) == 99) then
        labels(i) = 98
      end if
    end do
!
!   Convert to Cartesian coordinates
!
    xtoc = tvec
    k = numat - id
    i = id
    call gmetrn (geo, coord, labels, na, nb, nc)
!
    if (l_cell) then
      call cell(coord, tvec)
      labels(numat - 2:numat) = 99
    end if
    id = i
!
!  xtoc: Translation vectors of the input data
!  tvec: Smallest valid translation vectors (primitive vectors)
!  Tv_n: Smallest valid translation distances
!  M_Tv: Number of primitive vectors needed to conver tvec into xtoc
!  n_tv: Number of primitive vectors requested.
!
    if (abs(Tv_n(1) - Tv_N(2)) + abs(Tv_n(1) - Tv_n(3)) < 0.004) then
!
!  All three primitive translation distances are equal, so
!  use the order of Tv found in the input data.
!
      do i = 1, 3
        tvec(:,i) = xtoc(:,i)/max(1,M_Tv(i))
      end do
      if (index(keywrd, " BCC") /= 0 .and. M_Tv(1)*M_Tv(2)*M_Tv(3) /= 1) tvec = tvec*0.5d0
    else
!
! Sort tvec into its original sequence.  The primitive Tv distances were worked out
! from the atoms of the input data.
!
      xtoc = tvec
      do j = 1, id
        do i = 1, id
          sum = sqrt(xtoc(1,i)**2 + xtoc(2,i)**2 + xtoc(3,i)**2)
          if (abs(sum - tv_n(j)) < 0.002d0) then
            tvec(:,j) = xtoc(:,i)
            coord(:,j + k) = tvec(:,j)
            xtoc(:,i) = -1.d3
            tv_n(j) = -tv_n(j)
          end if
        end do
      end do
      tv_n = -tv_n
    end if
!
!  Check to see if a severe error has been made
!
    bug = .false.
    if (id == 3) then
      bug = (bug .or. ((tvec(1,1) - tvec(1,2))**2 + (tvec(2,1) - tvec(2,2))**2 + (tvec(3,1) - tvec(3,2))**2 < 1.d0))
      bug = (bug .or. ((tvec(1,1) - tvec(1,3))**2 + (tvec(2,1) - tvec(2,3))**2 + (tvec(3,1) - tvec(3,3))**2 < 1.d0))
      bug = (bug .or. ((tvec(1,2) - tvec(1,3))**2 + (tvec(2,2) - tvec(2,3))**2 + (tvec(3,2) - tvec(3,3))**2 < 1.d0))
    end if
    if (bug) then
      if (l_cell) then
        write(iw,'(a)')" A severe fault is present in keyword ""CELL"" "// &
          "Two translation vectors are almost the same."
      else
        if (index(keywrd, " GEO-OK") /= 0) then
          tvec = store_tvec
          coord(:,:numat) = store_coord(:,:numat)
          labels(:numat) = store_labels(:numat)
        else
          write(iw,'(a)')" A severe fault was found in trying to calculate the translation vectors.", &
          " Try using keyword ""CELL"" or add the keyword ""GEO-OK"".", &
          " If ""GEO-OK"" is used, then the supplied translation vectors will be used,", &
          " and no attempt will be made to calculate them.", &
          " See http://openmopac.net/manual/makpol.html"
        end if
      end if
      i =index(keywrd, " GEO-OK")
      if (i /= 0) then
        line = trim(keywrd)
        keywrd = line(:i - 1)//trim(line(i + 7:))
      else
        stop
      end if
    end if
    if (debug) then
      write(iw,'(/10x,a,i4)')"Number of atoms in input data:", numat - id
    end if
    na(:numat) = 0
    if (id == 0) id = i
!
!   Safety check - make sure that no two atoms are at the same position
!
    do i = 1, numat
      if (labels(i) == 98) then
        labels(i) = 99
      end if
    end do
    k = 0
    do i = 1, numat
      l = 0
      do j = 1, i - 1
        r = (coord(1, i)-coord(1, j))**2 + (coord(2, i)-coord(2, j))**2 + &
       &  (coord(3, i)-coord(3, j))**2
        if (r < 0.1d0 .and. labels(i) < 99 .and. labels(j) < 99 .and. labels(i) /= labels(j)) then
          write (iw, "(A,I4,A,I4,A)") " In the data set supplied" // ", atoms",&
         &  i, " and", j, " are at the same position"
          write (iw, "(A)") " (These atoms are of type " // elemnt (labels(i)) &
         & // " and " // elemnt (labels(j)) // " respectively)"
           stop
         end if
         if (r < 0.1d0 .and. labels(i) < 99 .and. labels(j) < 99 ) l = 1
      end do
      if (l == 0) then
        k = k + 1
        coord(1,k) = coord(1,i)
        coord(2,k) = coord(2,i) 
        coord(3,k) = coord(3,i)
        labels(k) = labels(i)
      end if
    end do
    numat = k
    do i = 1, numat
      if (labels(i) == 99) then
        labels(i) = 98
      end if
    end do
!
!   Remove dummy atoms from LABELS
!
    do i = numat, 1, -1
      if (labels(i) < 98) exit
    end do
    numat = i
    natoms = numat
!
!   The translation vectors are in TVEC
!
!  If keyword MERS is NOT present, write in a default value of MERS
!
    j = Index (keywrd, "MERS=")
    if (j == 0) then
      j = index(keywrd, "              ") + 1
      if (id == 1) then
         keywrd(j:j + 13) = "MERS=(1) "
      else if (id == 2) then
         keywrd(j:j + 13) = "MERS=(1,1) "
      else
         keywrd(j:j + 13) = "MERS=(1,1,1) "
      end if
    end if
    i = Index (keywrd(j:), " ") + j
    k = 0
    do l = 1, 3
      j = j + k
      if (l > 1 .and. k == 0) exit
      rmers(l) = reada (keywrd(j:), 1)
      k = Index (keywrd(j:i), ",")
    end do
    if( l-1 /= id) then
      j = Index (keywrd, "MERS=")
      write(iw,"(10x,3a,i3)")  "Keyword MERS is: "&
    &,'"'//keywrd(j:i - 2)//'",', " number of Tv in system:",id
      write(iw,"(10x,a)")"Correct error and re-submit"
      stop
    end if
    do i = 1, id
      if (Abs(rmers(i) - Nint(rmers(i)))< 0.001d0) then
        mers(i) = nint(rmers(i))
      else 
!
!  Special case:  Number of mers is not integer
!
        do j = 1, 100
          if (Abs(rmers(i)*j - Nint(rmers(i)*j)) < 0.01d0) exit
        end do
!
!   Multiply by j to get integer
!
        mers(i) = Nint(rmers(i)*j)
        do k = 1,3
          tvec(k,i) = tvec(k,i)/j
        end do
!
!  Do a quick check: is there an atom at the first translation distance?
!
        l = 0
        do k = 2, natoms
          do j = -1,1,2
            if (Abs(coord(1,k) - tvec(1,i)) + Abs(coord(2,k) - tvec(2,i)) &
            & + Abs(coord(3,k) - tvec(3,i)) < 0.1d0)  l = 1
          end do
          if (l == 1) exit
        end do
        if (l == 0) then
          write(iw,"(a,i2,3a)") "For translation vector",i," there are no atoms at", &
          & " the correct distance - check MERS"
          write(iw,"(2a)")"Keyword = ",keywrd(1:len_trim(keywrd))
          stop
        end if
      end if
    end do
    if (index(keywrd, " BCC") /= 0) then
      if (mers(1)*mers(2)*mers(3) /= 1) then
        do i = 1, 3
          if (mod(mers(i),2) /= 0) then
            write(iw,'(//10x,a)')" When BCC is used, all MERS must be even"
            stop
          end if
        end do  
      end if
    end if
!
! Write out keyword MERS again
!
    j = Index (keywrd, "MERS=")
    i = Index (keywrd(j:), " ") + j - 1
    keywrd(j:i) = " "
    n1 = char(ichar("1") + int(log10(mers(1)*1.01)))
    if (id > 1) n2 = char(ichar("1") + int(log10(mers(2)*1.01)))
    if (id > 2) n3 = char(ichar("1") + int(log10(mers(3)*1.01)))
    if (id == 1) then
        write(keywrd(j:i),'("MERS=(",i'//n1//'")")')mers(1)     
    else if (id == 2) then
      write(keywrd(j:i),'("MERS=(",i'//n1//',",",i'//n2//',")")')mers(1),mers(2)
    else 
      write(keywrd(j:i),'("MERS=(",i'//n1//',",",i'//n2//',",",i'//n3//'")")')mers
    endif
    
    mers(2) = Max (1, mers(2))
    mers(3) = Max (1, mers(3))
    store_mers = mers
!
! First pass - check that all atoms are real, and that there are no duplicates.
!
    mers = 2
    if (id < 3) mers(3) = 1
    if (id < 2) mers(2) = 1
!
!   If "LET" or "SORT" is present, then do NOT put all atoms into the CUC.
!   Use LET or SORT if band-structure work is to be done. 
!   If setcuc is called here, then the positions of atoms might get changed -
!   assume the user knows what is wanted.
!
    if (Index (keywrd, " LET") + Index (keywrd, " SORT") == 0) then
!
! Put all atoms into the central unit cell.
!
      call setcuc (coord, natoms, tvec, id, xtoc)
!
! Detect redundancies (atoms within 0.5A of each other are considered to be the same atom)
!
      do i = 2, numat
        do j = 1, i - 1
          if (labels(j) == 200) cycle
          r = dist(coord, i, j, tvec, xtoc, vec)
          if (r < 0.5d0 .and. labels(i) < 99 .and. labels(j) < 99) then
            labels(i) = 200
            exit
          end if
        end do
      end do
!
! Remove redundancies 
!
      k = 0
      do i = 1, numat      
        if (labels(i) /= 200) then
          k = k + 1
          coord(1,k) = coord(1,i)
          coord(2,k) = coord(2,i) 
          coord(3,k) = coord(3,i)
          labels(k) = labels(i)
        end if
      end do
      numat = k
    end if 
!
!   If "LET" or "SORT" is present, then do NOT join fragments together.
!   Use LET or SORT if band-structure work is to be done.
!
    if (Index (keywrd, " LET") + Index (keywrd, " SORT") == 0)  &
!
!   Join together fragments, so that discrete molecules can be recognized.
!
    call join_bits_together(coord, labels, numat, tvec, xtoc)   
    mers = store_mers
    nmax = numat * mers(1) * mers(2) * mers(3) + 2 * id
    i = index(keywrd, " CELL")
    if (i /= 0) then
        j = index(keywrd(i:), ")") + i 
        keywrd(i:) = trim(keywrd(j:))
    end if
    i = Index (keywrd, " SORT")
    if (i /= 0) keywrd(i:i + 4) = " "
    i = 0
    call maklay (labels, coord, geo, na, nb, nc, mers, loc, nmax) 
end program makpol
double precision function snapth (theta)
!
!.. Implicit Declarations ..
    implicit none
!
!.. Formal Arguments ..
    double precision, intent (in) :: theta
!
!.. Local Scalars ..
    integer :: i, j, k
    double precision :: angle, const, phase, sum
!
!.. Intrinsic Functions ..
    intrinsic Abs, Acos, Asin, Cos, Int, Mod, Nint, Sign, Sqrt
!
! ... Executable Statements ...
!
    angle = Cos (theta)
    phase = Sign (1.d0, theta)
    if (Abs (angle) < 1.d-4) then
!
!   Cos(Theta) is zero - this is 90 or 270 degrees
!
      const = 2 * Asin (1.d0)
      if (Abs (theta) < const) then
        snapth = phase * Acos (0.d0)
      else
        snapth = phase * Acos (0.d0) + const
      end if
    else
      sum = Abs (1/angle)**2
      do i = 1, 7
        j = Nint (sum*i)
        if (Abs (j - sum*i) < 1.d-4) go to 1000
      end do
      snapth = theta
      return
1000  sum = Sqrt ((1.d0*i)/j)
      const = 2 * Asin (1.d0)
      k = Int (Abs (theta)/const)
      if (Mod(k, 2) == 0) then
!
!   Theta is in domain -180 - 0 - +180 degrees
!
        snapth = phase * Acos (sign(sum, angle))
      else
!
!   Theta is in domain -180 - 360 - +180 degrees
!
        snapth = phase * (2*const-Acos (sign(sum, angle)))
      end if
    end if
end function snapth
double precision function dist(coord, i, j, tvec, xtoc, vec)
  use common_systm, only :  id
  double precision :: coord(3,*), tvec(3,3), xtoc(3,3), vec(3)
  integer :: i, j
!
  double precision, dimension (3) :: xyz
  integer :: k, l, jj
  do k = 1,3
    vec(k) = coord(k,i) - coord(k,j)
  end do
      xyz = 0.d0
      do k = 1, 3        
        do l = 1, id
          xyz(l) = xyz(l) + vec(k) * xtoc(l, k)
        end do
      end do
!
!   XYZ holds the crystal fractional coordinates
!
      do jj = 1, id
        l = Int (xyz(jj)+2000.5d0) - 2000
        do k = 1, 3
          vec(k) = vec(k) - l * tvec(k, jj)
        end do
      end do
      dist = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
      return
end function dist
subroutine join_bits_together(coord, labels, numat, tvec, xtoc)  
  use common_systm, only : id
  use mod_atomradii, only: atom_radius_covalent, atom_radius_vdw
  implicit none
  integer :: numat, labels(numat)
  double precision, dimension (3,numat) :: coord
  double precision, dimension (3,3) :: tvec, xtoc
!
!  Local
!
  logical, allocatable, dimension (:,:) :: bonded
  logical, dimension (numat) :: used
  double precision :: coord_local(3,numat), xyz(3), vec(3)
  integer, dimension (numat) :: fragment, bonded_to
!
  double precision :: dist, Rij, sum, all_rij(107, 107)
  integer :: i, j, l, ii, jj, k, infrag
  external dist 
  allocate (bonded(numat,numat))
  atom_radius_covalent =  atom_radius_vdw*0.6d0
  do i = 1, 83
    do j = 1,83
      all_rij(i,j) = atom_radius_covalent(i) + atom_radius_covalent(j)
    end do
  end do
!
!  Add specific RIJ here
!
   all_rij(8,1)   = 1.4d0! Oxygen - hydrogen
   all_rij(1,8)   = 1.4d0! Oxygen - hydrogen
   all_rij(35,5)  = 2.3d0! Boron - bromine
   all_rij(5,35)  = 2.3d0! Boron - bromine
   all_rij(8,11)  = 2.38d0! Sodium - oxygen
   all_rij(11,8)  = 2.38d0! Sodium - oxygen
   all_rij(7,12)  = 2.3d0! Magnesium - nitrogen
   all_rij(12,7)  = 2.3d0! Magnesium - nitrogen
   all_rij(8,12)  = 2.2d0! Magnesium - oxygen
   all_rij(12,8)  = 2.2d0! Magnesium - oxygen
   all_rij(35,12) = 2.8d0! Magnesium - bromine
   all_rij(12,35) = 2.8d0! Magnesium - bromine
   all_rij(53,12) = 2.8d0! Magnesium - iodine
   all_rij(12,53) = 2.8d0! Magnesium - iodine
   all_rij(13,53) = 2.7d0! Aluminum - iodine
   all_rij(53,13) = 2.7d0! Aluminum - iodine
   all_rij(8,15)  = 1.7d0! Phosphorus - oxygen
   all_rij(15,8)  = 1.7d0! Phosphorus - oxygen
   all_rij(19,7)  = 2.85d0! Potassium - nitrogen
   all_rij(7,19)  = 2.85d0! Potassium - nitrogen
   all_rij(53,15) = 2.6d0! Phosphorus - iodine
   all_rij(15,53) = 2.6d0! Phosphorus - iodine 
   all_rij(19,8)  = 2.5d0! Potassium - oxygen
   all_rij(8,19)  = 2.5d0! Potassium - oxygen
   all_rij(8,20)  = 2.3d0! Calcium - oxygen
   all_rij(20,8)  = 2.3d0! Calcium - oxygen
   all_rij(53,20) = 3.3d0! Calcium - iodine
   all_rij(20,53) = 3.3d0! Calcium - iodine
   all_rij(21,8)  = 2.6d0! Scandium - oxygen
   all_rij(8,21)  = 2.6d0! Scandium - oxygen
   all_rij(22,7)  = 2.3d0! Titanium - nitrogen
   all_rij(7,22)  = 2.3d0! Titanium - nitrogen
   all_rij(22,8)  = 2.6d0! Titanium - oxygen
   all_rij(8,22)  = 2.6d0! Titanium - oxygen
   all_rij(22,17) = 2.6d0! Titanium - chlorine
   all_rij(17,22) = 2.6d0! Titanium - chlorine 
   all_rij(22,35) = 2.8d0! Titanium - bromine
   all_rij(35,22) = 2.8d0! Titanium - bromine 
   all_rij(22,53) = 2.9d0! Titanium - iodine
   all_rij(53,22) = 2.9d0! Titanium - iodine 
   all_rij(23,8)  = 2.4d0! Vanadium - oxygen  
   all_rij(8,23)  = 2.4d0! Vanadium - oxygen
   all_rij(23,35) = 2.5d0! Vanadium - bromine  
   all_rij(35,23) = 2.5d0! Vanadium - bromine
   all_rij(24,17) = 2.55d0! Chromium - chlorine
   all_rij(17,24) = 2.55d0! Chromium - chlorine
   all_rij(25,6)  = 1.9d0! Manganese - carbon  
   all_rij(6,25)  = 1.9d0! Manganese - carbon
   all_rij(25,8)  = 2.4d0! Manganese - oxygen  
   all_rij(8,25)  = 2.4d0! Manganese - oxygen
   all_rij(25,17) = 2.7d0! Manganese - chlorine  (MnCl6(4-))
   all_rij(17,25) = 2.7d0! Manganese - chlorine
   all_rij(25,35) = 2.8d0! Manganese - bromine  (MnBr6(4-))
   all_rij(35,25) = 2.8d0! Manganese - bromine
   all_rij(26,17) = 2.4d0! Iron - chlorine
   all_rij(17,26) = 2.4d0! Iron - chlorine
   all_rij(27,16) = 2.4d0! Cobalt - sulfur
   all_rij(16,27) = 2.4d0! Cobalt - sulfur
   all_rij(27,35) = 2.6d0! Cobalt - bromine
   all_rij(35,27) = 2.6d0! Cobalt - bromine
   all_rij(27,53) = 2.7d0! Cobalt - iodine
   all_rij(53,27) = 2.7d0! Cobalt - iodine 
   all_rij(28,6)  = 2.2d0! Nickel - carbon
   all_rij(6,28)  = 2.2d0! Nickel - carbon
   all_rij(28,7)  = 2.3d0! Nickel - nitrogen 
   all_rij(7,28)  = 2.3d0! Nickel - nitrogen
   all_rij(28,8)  = 2.6d0! Nickel - oxygen  (NiO6)
   all_rij(8,28)  = 2.6d0! Nickel - oxygen
   all_rij(28,16) = 2.4d0! Nickel - sulfur
   all_rij(16,28) = 2.4d0! Nickel - sulfur
   all_rij(28,17) = 2.4d0! Nickel - chlorine
   all_rij(17,28) = 2.4d0! Nickel - chlorine
   all_rij(29,7)  = 2.2d0! Copper - nitrogen 
   all_rij(7,29)  = 2.2d0! Copper - nitrogen
   all_rij(29,8)  = 2.6d0! Copper - oxygen (CuO6)
   all_rij(8,29)  = 2.6d0! Copper - oxygen
   all_rij(29,17) = 2.6d0! Copper - chlorine
   all_rij(17,29) = 2.6d0! Copper - chlorine
   all_rij(29,35) = 3.1d0! Copper - bromine
   all_rij(35,29) = 3.1d0! Copper - bromine
   all_rij(29,53) = 3.1d0! Copper - iodine
   all_rij(53,29) = 3.1d0! Copper - iodine
   all_rij(30,35) = 2.6d0! Zinc - bromine
   all_rij(35,30) = 2.6d0! Zinc - bromine 
   all_rij(31,8)  = 2.2d0! Gallium - oxygen
   all_rij(8,31)  = 2.2d0! Gallium - oxygen 
   all_rij(31,16) = 2.5d0! Gallium - sulfur
   all_rij(16,31) = 2.5d0! Gallium - sulfur 
   all_rij(31,17) = 2.5d0! Gallium - chlorine
   all_rij(17,31) = 2.5d0! Gallium - chlorine 
   all_rij(31,35) = 2.5d0! Gallium - bromine
   all_rij(35,31) = 2.5d0! Gallium - bromine 
   all_rij(31,53) = 2.6d0! Gallium - iodine
   all_rij(53,31) = 2.6d0! Gallium - iodine 
   all_rij(32,35) = 2.6d0! Germanium - bromine
   all_rij(35,32) = 2.6d0! Germanium - bromine 
   all_rij(32,53) = 2.6d0! Germanium - iodine
   all_rij(53,32) = 2.6d0! Germanium - iodine 
   all_rij(33,16) = 2.4d0! Arsenic - sulfur
   all_rij(16,33) = 2.4d0! Arsenic - sulfur 
   all_rij(33,35) = 2.5d0! Arsenic - bromine
   all_rij(35,33) = 2.5d0! Arsenic - bromine 
   all_rij(33,53) = 2.8d0! Arsenic - iodine
   all_rij(53,33) = 2.8d0! Arsenic - iodine 
   all_rij(47,6)  = 2.6d0! Silver - carbon
   all_rij(6,47)  = 2.6d0! Silver - carbon 
   all_rij(48,8)  = 2.9d0! Cadmium - oxygen
   all_rij(8,48)  = 2.9d0! Cadmium - oxygen 
   all_rij(34,17) = 2.4d0! Selenium - chlorine
   all_rij(17,34) = 2.4d0! Selenium - chlorine 
   all_rij(34,30) = 2.6d0! Selenium - zinc
   all_rij(30,34) = 2.6d0! Selenium - zinc 
   all_rij(34,34) = 2.5d0! Selenium - selenium
   all_rij(34,34) = 2.5d0! Selenium - selenium 
   all_rij(34,35) = 2.7d0! Selenium - bromine
   all_rij(35,34) = 2.7d0! Selenium - bromine 
   all_rij(34,48) = 2.8d0! Selenium - cadmium
   all_rij(48,34) = 2.8d0! Selenium - cadmium 
   all_rij(34,53) = 2.7d0! Selenium - iodine
   all_rij(53,34) = 2.7d0! Selenium - iodine 
   all_rij(34,80) = 2.8d0! Selenium - mercury
   all_rij(80,34) = 2.8d0! Selenium - mercury 
   all_rij(38,7)  = 3.0d0! Strontium - nitrogen
   all_rij(7,38)  = 3.0d0! Strontium - nitrogen 
   all_rij(38,53) = 3.6d0! Strontium - iodine
   all_rij(53,38) = 3.6d0! Strontium - iodine 
   all_rij(39,35) = 3.0d0! Yttrium - bromine  
   all_rij(35,39) = 3.0d0! Yttrium - bromine 
   all_rij(39,53) = 3.1d0! Yttrium - iodine  
   all_rij(53,39) = 3.1d0! Yttrium - iodine
   all_rij(42,8)  = 2.6d0! Molybdenum - oxygen
   all_rij(8,42)  = 2.6d0! Molybdenum - oxygen
   all_rij(42,17) = 2.7d0! Molybdenum - chlorine
   all_rij(17,42) = 2.7d0! Molybdenum - chlorine 
   all_rij(42,35) = 2.7d0! Molybdenum - bromine 
   all_rij(35,42) = 2.7d0! Molybdenum - bromine 
   all_rij(43,15) = 2.6d0! Technetium - phosphorus  
   all_rij(15,43) = 2.6d0! Technetium - phosphorus
   all_rij(43,16) = 2.6d0! Technetium - sulfur  
   all_rij(16,43) = 2.6d0! Technetium - sulfur
   all_rij(43,35) = 2.6d0! Technetium - bromine
   all_rij(35,43) = 2.6d0! Technetium - bromine
   all_rij(43,53) = 2.6d0! Technetium - iodine
   all_rij(53,43) = 2.6d0! Technetium - iodine
   all_rij(44,16) = 2.5d0! Ruthenium - sulfur  
   all_rij(16,44) = 2.5d0! Ruthenium - sulfur
   all_rij(44,17) = 2.5d0! Ruthenium - chlorine  
   all_rij(17,44) = 2.5d0! Ruthenium - chlorine
   all_rij(44,35) = 2.8d0! Ruthenium - bromine  
   all_rij(35,44) = 2.8d0! Ruthenium - bromine
   all_rij(44,53) = 2.9d0! Ruthenium - iodine  
   all_rij(53,44) = 2.9d0! Ruthenium - iodine
   all_rij(45,35) = 2.6d0! Rhodium - bromine  
   all_rij(35,45) = 2.6d0! Rhodium - bromine
   all_rij(45,53) = 2.8d0! Rhodium - iodine  
   all_rij(53,45) = 2.8d0! Rhodium - iodine
   all_rij(46,6)  = 2.3d0! Palladium - carbon
   all_rij(6,46)  = 2.3d0! Palladium - carbon 
   all_rij(46,7)  = 2.4d0! Palladium - nitrogen
   all_rij(7,46)  = 2.4d0! Palladium - nitrogen
   all_rij(46,16) = 2.5d0! Palladium - sulfur  
   all_rij(16,46) = 2.5d0! Palladium - sulfur
   all_rij(46,17) = 2.5d0! Palladium - chlorine  
   all_rij(17,46) = 2.5d0! Palladium - chlorine
   all_rij(46,35) = 2.5d0! Palladium - bromine  
   all_rij(35,46) = 2.5d0! Palladium - bromine
   all_rij(46,53) = 2.8d0! Palladium - iodine  
   all_rij(53,46) = 2.8d0! Palladium - iodine
   all_rij(49,6)  = 2.6d0! Indium - carbon  
   all_rij(6,49)  = 2.6d0! Indium - carbon
   all_rij(49,7)  = 2.5d0! Indium - nitrogen   
   all_rij(7,49)  = 2.5d0! Indium - nitrogen
   all_rij(49,15) = 2.6d0! Indium - phosphorus  
   all_rij(15,49) = 2.6d0! Indium - phosphorus
   all_rij(49,16) = 2.7d0! Indium - sulfur  
   all_rij(16,49) = 2.7d0! Indium - sulfur
   all_rij(49,17) = 2.7d0! Indium - chlorine  
   all_rij(17,49) = 2.7d0! Indium - chlorine
   all_rij(49,35) = 2.8d0! Indium - bromine  
   all_rij(35,49) = 2.8d0! Indium - bromine
   all_rij(49,52) = 2.9d0! Indium - tellurium  
   all_rij(52,49) = 2.9d0! Indium - tellurium
   all_rij(49,53) = 2.9d0! Indium - iodine  
   all_rij(53,49) = 2.9d0! Indium - iodine
   all_rij(50,6)  = 2.6d0! Tin - carbon  
   all_rij(6,50)  = 2.6d0! Tin - carbon
   all_rij(50,16) = 2.7d0! Tin - sulfur  
   all_rij(16,50) = 2.7d0! Tin - sulfur
   all_rij(50,17) = 2.5d0! Tin - chlorine  
   all_rij(17,50) = 2.5d0! Tin - chlorine
   all_rij(50,35) = 2.5d0! Tin - bromine  
   all_rij(35,50) = 2.5d0! Tin - bromine
   all_rij(51,6)  = 2.6d0! Antimony - carbon
   all_rij(6,51)  = 2.6d0! Antimony - carbon
   all_rij(51,8)  = 2.3d0! Antimony - oxygen
   all_rij(8,51)  = 2.3d0! Antimony - oxygen
   all_rij(51,16) = 2.9d0! Antimony - sulfur
   all_rij(16,51) = 2.9d0! Antimony - sulfur
   all_rij(51,17) = 3.0d0! Antimony - chlorine
   all_rij(17,51) = 3.0d0! Antimony - chlorine
   all_rij(51,35) = 3.0d0! Antimony - bromine
   all_rij(35,51) = 3.0d0! Antimony - bromine
   all_rij(51,51) = 2.9d0! Antimony - antimony
   all_rij(51,51) = 2.9d0! Antimony - antimony
   all_rij(51,53) = 3.5d0! Antimony - iodine
   all_rij(53,51) = 3.5d0! Antimony - iodine
   all_rij(52,8)  = 2.7d0! Tellurium - oxygen
   all_rij(8,52)  = 2.7d0! Tellurium - oxygen
   all_rij(52,16) = 2.9d0! Tellurium - sulfur
   all_rij(16,52) = 2.9d0! Tellurium - sulfur
   all_rij(52,17) = 2.7d0! Tellurium - chlorine
   all_rij(17,52) = 2.7d0! Tellurium - chlorine
   all_rij(52,35) = 2.8d0! Tellurium - bromine
   all_rij(35,52) = 2.8d0! Tellurium - bromine
   all_rij(52,52) = 3.4d0! Tellurium - tellurium
   all_rij(52,52) = 3.4d0! Tellurium - tellurium
   all_rij(52,53) = 3.1d0! Tellurium - iodine
   all_rij(53,52) = 3.1d0! Tellurium - iodine
   all_rij(55,7)  = 3.0d0! Cesium - nitrogen
   all_rij(7,55)  = 3.0d0! Cesium - nitrogen
   all_rij(56,8)  = 2.6d0! Barium - oxygen
   all_rij(8,56)  = 2.6d0! Barium - oxygen
   all_rij(56,35) = 3.4d0! Barium - bromine
   all_rij(35,56) = 3.4d0! Barium - bromine
! all_rij(56,56) = 2.4d0! Barium - barium
! all_rij(56,56) = 2.6d0! Barium - barium
   all_rij(57,7)  = 2.9d0! Lanthanum - nitrogen  
   all_rij(7,57)  = 2.9d0! Lanthanum - nitrogen
   all_rij(57,8)  = 2.7d0! Lanthanum - oxygen  
   all_rij(8,57)  = 2.7d0! Lanthanum - oxygen
   all_rij(57,17) = 3.0d0! Lanthanum - chlorine  
   all_rij(17,57) = 3.0d0! Lanthanum - chlorine
   all_rij(57,35) = 3.2d0! Lanthanum - bromine  
   all_rij(35,57) = 3.2d0! Lanthanum - bromine
   all_rij(57,53) = 3.4d0! Lanthanum - iodine  
   all_rij(53,57) = 3.4d0! Lanthanum - iodine
   all_rij(71,6)  = 2.9d0! Lutetium - carbon   
   all_rij(6,71)  = 2.9d0! Lutetium - carbon
   all_rij(71,7)  = 2.6d0! Lutetium - nitrogen   
   all_rij(7,71)  = 2.6d0! Lutetium - nitrogen
!   all_rij(62,8)  = 2.43d0! Samarium - oxygen   
!   all_rij(8,62)  = 2.43d0! Samarium - oxygen
   all_rij(71,8)  = 2.7d0! Lutetium - oxygen   
   all_rij(8,71)  = 2.7d0! Lutetium - oxygen
   all_rij(71,17) = 2.6d0! Lutetium - chlorine
   all_rij(17,71) = 2.6d0! Lutetium - chlorine
   all_rij(71,35) = 2.6d0! Lutetium - bromine
   all_rij(35,71) = 2.6d0! Lutetium - bromine
   all_rij(71,53) = 3.0d0! Lutetium - iodine
   all_rij(53,71) = 3.0d0! Lutetium - iodine
   all_rij(72,6)  = 2.7d0! Hafnium - carbon   
   all_rij(6,72)  = 2.7d0! Hafnium - carbon
   all_rij(72,17) = 2.6d0! Hafnium - chlorine   
   all_rij(17,72) = 2.6d0! Hafnium - chlorine
   all_rij(73,6)  = 2.6d0! Tantalum - carbon   
   all_rij(6,73)  = 2.6d0! Tantalum - carbon
   all_rij(73,16) = 2.6d0! Tantalum - sulfur   
   all_rij(16,73) = 2.6d0! Tantalum - sulfur
   all_rij(74,7)  = 2.4d0! Tungsten - nitrogen 
   all_rij(7,74)  = 2.4d0! Tungsten - nitrogen
   all_rij(74,8)  = 2.6d0! Tungsten - oxygen 
   all_rij(8,74)  = 2.6d0! Tungsten - oxygen
   all_rij(74,35) = 2.7d0! Tungsten - bromine
   all_rij(35,74) = 2.7d0! Tungsten - bromine
   all_rij(74,53) = 3.0d0! Tungsten - iodine
   all_rij(53,74) = 3.0d0! Tungsten - iodine
   all_rij(75,8)  = 2.4d0! Rhenium - oxygen
   all_rij(8,75)  = 2.4d0! Rhenium - oxygen
   all_rij(75,35) = 2.7d0! Rhenium - bromine
   all_rij(35,75) = 2.7d0! Rhenium - bromine
   all_rij(75,53) = 2.9d0! Rhenium - iodine
   all_rij(53,75) = 2.9d0! Rhenium - iodine
   all_rij(76,35) = 2.9d0! Osmium - bromine  
   all_rij(35,76) = 2.9d0! Osmium - bromine
   all_rij(76,53) = 3.0d0! Osmium - iodine  
   all_rij(53,76) = 3.0d0! Osmium - iodine
   all_rij(78,6)  = 2.3d0! Platinum - carbon
   all_rij(6,78)  = 2.3d0! Platinum - carbon
   all_rij(78,7)  = 2.4d0! Platinum - nitrogen
   all_rij(7,78)  = 2.4d0! Platinum - nitrogen 
   all_rij(78,17) = 2.5d0! Platinum - chlorine   
   all_rij(17,78) = 2.5d0! Platinum - chlorine
   all_rij(78,35) = 2.6d0! Platinum - bromine  
   all_rij(35,78) = 2.6d0! Platinum - bromine
   all_rij(78,53) = 2.9d0! Platinum - iodine  
   all_rij(53,78) = 2.9d0! Platinum - iodine
   all_rij(79,6)  = 2.3d0! Gold - carbon   
   all_rij(6,79)  = 2.3d0! Gold - carbon
   all_rij(79,15) = 2.6d0! Gold - phosphorus   
   all_rij(15,79) = 2.6d0! Gold - phosphorus
   all_rij(79,16) = 2.4d0! Gold - sulfur   
   all_rij(16,79) = 2.4d0! Gold - sulfur
   all_rij(79,17) = 2.5d0! Gold - chlorine   
   all_rij(17,79) = 2.5d0! Gold - chlorine
   all_rij(79,35) = 2.7d0! Gold - bromine   
   all_rij(35,79) = 2.7d0! Gold - bromine
   all_rij(79,53) = 2.8d0! Gold - iodine   
   all_rij(53,79) = 2.8d0! Gold - iodine
   all_rij(80,6)  = 2.4d0! Mercury - carbon 
   all_rij(6,80)  = 2.4d0! Mercury - carbon   
   all_rij(80,8)  = 3.0d0! Mercury - oxygen 
   all_rij(8,80)  = 3.0d0! Mercury - oxygen
   all_rij(80,17) = 2.5d0! Mercury - chlorine  
   all_rij(17,80) = 2.5d0! Mercury - chlorine
   all_rij(80,35) = 2.5d0! Mercury - bromine   
   all_rij(35,80) = 2.5d0! Mercury - bromine
   all_rij(80,53) = 2.5d0! Mercury - iodine 
   all_rij(53,80) = 2.5d0! Mercury - iodine
   all_rij(82,8)  = 2.5d0! Thallium - oxygen 
   all_rij(8,81)  = 2.5d0! Thallium - oxygen  
   all_rij(81,17) = 3.1d0! Thallium - chlorine
   all_rij(17,81) = 3.1d0! Thallium - chlorine
   all_rij(81,35) = 3.1d0! Thallium - bromine 
   all_rij(35,81) = 3.1d0! Thallium - bromine
   all_rij(81,53) = 3.1d0! Thallium - iodine
   all_rij(53,81) = 3.1d0! Thallium - iodine
   all_rij(82,6)  = 2.7d0! Lead - carbon
   all_rij(6,82)  = 2.7d0! Lead - carbon
   all_rij(82,8)  = 3.0d0! Lead - oxygen
   all_rij(8,82)  = 3.0d0! Lead - oxygen
   all_rij(82,17) = 3.1d0! Lead - chlorine
   all_rij(17,82) = 3.1d0! Lead - chlorine
   all_rij(82,35) = 3.1d0! Lead - bromine
   all_rij(35,82) = 3.1d0! Lead - bromine
   all_rij(82,35) = 3.3d0! Lead - iodine
   all_rij(53,82) = 3.3d0! Lead - iodine
   all_rij(83,8)  = 2.7d0! Bismuth - oxygen
   all_rij(8,83)  = 2.7d0! Bismuth - oxygen
   all_rij(83,17) = 3.0d0! Bismuth - chlorine
   all_rij(17,83) = 3.0d0! Bismuth - chlorine
   all_rij(83,35) = 3.2d0! Bismuth - bromine
   all_rij(35,83) = 3.2d0! Bismuth - bromine
   all_rij(83,53) = 3.2d0! Bismuth - iodine
   all_rij(53,83) = 3.2d0! Bismuth - iodine
!
!   Special functions
!
!  Cd_present = .false.
!  do i = 1, numat
!    if (labels(i) == 48) Cd_present = .true.
!  end do
!  if (Cd_present) then
!    all_rij(19,7)  = 2.5d0! Potassium - nitrogen
!    all_rij(7,19)  = 2.5d0! Potassium - nitrogen
!    all_rij(19,8)  = 2.5d0! Potassium - oxygen
!    all_rij(8,19)  = 2.5d0! Potassium - oxygen
!  end if
   


   do i = 1, numat
    do j = 1, i - 1
      Rij = all_rij(labels(i),labels(j))      
      sum = dist(coord, i, j, tvec, xtoc, vec)
      if ((labels(i) == 35)) then
        Rij = Rij
      end if
      bonded(i,j) = (sum < Rij .and. (labels(i) /= 1 .or. labels(j) /= 1))
      bonded(j,i) = bonded(i,j)     
    end do
    bonded(i,i) = .true.
  end do
!
! Avoid bridging hydrogens (short hydrogen bonds, etc.)
!
  do i = 1, numat
    if (labels(i) == 1) then
      sum = 1.d10
      do j = 1, numat
        if (.not. bonded(i,j)) cycle
        if (i == j) cycle
        Rij = dist(coord, i, j, tvec, xtoc, vec)
        if (Rij < sum) then
          sum = Rij
          k = j
        end if
      end do
      do j = 1, numat
        bonded(j,i) = (j == k .or. j == i)
      end do
    end if
  end do
!
!  Start making up fragments
!
  used = .false.  
  do i = 1, numat
    if ( .not. used(i)) then 
      fragment = 0
      infrag = 1
      fragment(infrag) = i
      do j = 1, numat
        jj = fragment(j)
        if (jj == 0) exit
        do k = 1, numat
          if (bonded(jj,k) .and. .not. used(k)) then
!
!  atom ii (in the fragment) is bonded to atom k.  Check if k has already been 
!  added to the fragment.
!
            do ii = 1, infrag
              if (fragment(ii) == k) exit
            end do
            if (ii > infrag) then
              infrag = infrag + 1
              fragment(infrag) = k
              bonded_to(infrag) = jj
            end if
          end if
        end do
      end do
!
!  At this point, a complete molecule exists in "fragment".
!  Mark all the atoms as used.
!
      do j = 1, infrag
        used(fragment(j)) = .true.
        do k = 1,3
          coord_local(k,j) = coord(k, fragment(j))
        end do  
      end do
!
!  Join all the atoms together.  Take each atom in turn and find where it's neighbor is.
!
      do ii = 2, infrag
        Rij = dist(coord, fragment(ii), bonded_to(ii), tvec, xtoc, vec)
        do j = 1, ii
          if (fragment(j) == bonded_to(ii)) exit
        end do
        do k = 1, 3
          coord_local(k, ii) = coord_local(k, j) + vec(k)
        end do      
      end do
!
!  Put center of fragment inside the unit cell (stops molecules from being mainly outside the cell)
!
      vec = 0.d0        
      do ii = 1, infrag
        do k = 1,3
          vec(k) = vec(k) + coord_local(k,ii)
        end do       
      end do
      vec = vec / infrag!  vec now holds the center of the molecule
      xyz = 0.d0
      do k = 1, 3        
        do j = 1, id
          xyz(j) = xyz(j) + vec(k) * xtoc(j, k)
        end do
      end do
!
!  If atom 1 is involved, and it's at the origin, do not move fragment
!
      do j = 1, infrag
        if (fragment(j) == 1) exit
      end do
      if (j < infrag + 1) then
        if (coord_local(1,j)**2 + coord_local(2,j)**2 + coord_local(2,j)**2  < 0.1d0) xyz = 0.d0
      end if
!
!   XYZ holds the crystal fractional coordinates of the center
!
      do ii = 1, infrag
        do j = 1, id
          l = Int (xyz(j) + 20.7d0) - 20! was 20.5d0, then 20.7
          do k = 1, 3
            coord_local(k, ii) = coord_local(k, ii) - l * tvec(k, j)
          end do
        end do
      end do
      do ii = 1, infrag
        do k = 1, 3
          coord(k,fragment(ii)) = coord_local(k,ii)
        end do
      end do
    end if
  end do
end subroutine join_bits_together
