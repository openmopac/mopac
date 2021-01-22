subroutine getpdb (geo) 
    use molkst_C, only: maxtxt, natoms, numat, na1, keywrd, line, &
      ncomments, nbreaks, title, txtmax, pdb_label
    use common_arrays_C, only : txtatm, atmass, na, nb, nc, break_coords, &
      labels, lopt, coord, all_comments, breaks, txtatm1
    use chanel_C, only: iw, ir
    use parameters_C, only: ams
    use elemts_C, only : elemnt
    use reada_I
    implicit none
 !
 !***********************************************************************
 !
 !  GETPDB READS IN THE GEOMETRY, IN BROOKHAVEN PROTEIN DATA BASE FORMAT.
 !
 !  ON INPUT   IR  = CHANNEL NUMBER FOR READ, NORMALLY 5
 !
 ! ON OUTPUT LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
 !           GEO    = CARTESIAN COORDINATES, IN ANGSTROMS
 !           LOPT   = INTEGER ARRAY FILLED WITH '1'S.
 !                    (PDB FILES DO NOT SUPPORT OPTIMIZATION FLAGS)
 !           NA     = INTEGER ARRAY FILLED WITH '0'S.
 !           NB     = INTEGER ARRAY FILLED WITH '0'S.
 !           NC     = INTEGER ARRAY FILLED WITH '0'S.
 !***********************************************************************\
    double precision :: geo(3,natoms)
!
    integer, parameter :: maxel = 115
    character, save :: comma, space
    character :: typea, typer, ch
    character :: ele*3, letter*26
    logical :: leadsp, lxyz, lchain, last_atom = .true., l_pdb, first = .true., &
      l_letter(26)
    integer :: i, icomma, ii, j, k, khar, label, nline, npdb, nvalue, n_water = 0, &
      old_natoms = 0, numerr
    integer :: new_elements, defined_elements, previous_res, current_res
    double precision :: degree
    character :: commas(20), new_key_1*200, new_key_2*200 = " ", &
      line1*400
    character (len=2), dimension (maxel), save :: element
    character (len=4), dimension (20) :: txtpdb
    integer, dimension (20) :: ntxt_loc
    integer, dimension (maxel) :: ielem
    character, allocatable :: tmp_comments(:)*81
    intrinsic Abs, Ichar, Index, Nint
    data commas / 20 * "," /
    data ntxt_loc / 20 * 2 /
    data comma, space, txtpdb / ",", 21 * " " /
 !
 ! ... Executable Statements ...
 !
    do i = 1, 93
      if (elemnt(i)(1:1) == " ") then
        element(i) = elemnt(i)(2:2)//elemnt(i)(1:1)
      else
        element(i) = elemnt(i)(1:1)//char(ichar(elemnt(i)(2:2)) + ichar("A") - ichar("a"))
      end if
      ielem(i) = i
    end do
    element(i) = "X "
    ielem(i)   = 99
    i = i + 1
    element(i) = "L "
    ielem(i)   = 0
    i = i + 1
    element(i) = "D "
    ielem(i)   = 1
    defined_elements = i
    numerr = 0
    new_elements = 0 
    line = trim(keywrd)
    j = 0
    do i = 1, len_trim(keywrd)
      if (keywrd(i:i) == '"') j = 1 - j
      if (j == 1) keywrd(i:i) = " "
    end do
    i = Index (keywrd, " CHAINS")
    if (i > 0) then
      i = index(keywrd(i:), "(") + i
      j = index(keywrd(i:), ") ") + i - 2
      letter = keywrd(i:j)
    else
      letter = " "
    end if
    i = Index (keywrd, " PDB(")
    l_pdb = (Index (keywrd, " PDB ") /= 0 .or. i /= 0)
    keywrd = trim(line)
    if (i /= 0) then
!
!  THE USER IS GOING TO DEFINE VARIOUS SYMBOLS AND ATOMIC NUMBERS
!
      i = i + 5
      j = Index (keywrd(i:241), ")") + i
      if (j /= i) then
        k = defined_elements + 1
        new_elements = 1
        do
          element(k) = keywrd(i:j)
          ii = Index (element(k), ":")
          if (ii /= 0) then
            element (k) (ii:2) = " "
          end if
          ii = Index (keywrd(i:j), ":") + i
          if (ii == i) exit
          ielem(k) = Nint (reada (keywrd, ii))
          i = Index (keywrd(ii:j), ",") + ii
          if (i /= ii) then
            k = k + 1
            new_elements = new_elements +1
          else
            go to 1000
          end if
        end do
      end if
      call mopend ("KEYWORD PDB IS INCORRECT")
    end if
1000 i = Index (keywrd, " ALT_A=")
    if (i /= 0) then
      typea = keywrd(i+7:i+7)
      if (typea == "*") typea = " "
    else
      typea = "A"
    end if
    i = Index (keywrd, " ALT_R=")
    if (i /= 0) then
      typer = keywrd(i+7:i+7)
      if (typer == "*") typer = " "
    else
      typer = "A"
    end if
    maxtxt = 0
    natoms = 0
    numat = 0
    npdb = 0
    maxtxt = 26
    nline = 0
    previous_res = 2000
    ncomments = 0
    nbreaks = 0
    new_key_1 = " "
    lchain = .true.
    if (.not. allocated(tmp_comments)) allocate(tmp_comments(10000))
!
    outer_loop: do
!
      read (ir, "(A)", end=1020, err=1020) line
!
      nline = nline + 1
      if (natoms > 0 .and. line == " ") exit
      if( .not. ( line(1:4) == "ATOM" .or. line(1:6) == "HETATM" ) ) then
        if (index(line,"ATOM  ") + index(line,"HETATM") + index(line,"TITLE ") + index(line,"HEADER") + &
                                   index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
            index(line,"HELIX ") + index(line,"SHEET ") + index(line,"REMARK") + index(line,"USER  ") + &
            index(line,"EXPDTA") + index(line,"AUTHOR") + index(line,"REVDAT") + index(line,"JRNL  ") + &
            index(line,"DBREF ") + index(line,"SEQRES") + index(line,"HET   ") + index(line,"HETNAM") + &
            index(line,"LINK  ") + index(line,"CRYST1") + index(line,"SCALE" ) + index(line,"ORIGX" ) + &
            index(line,"FORMUL") + index(line,"SEQRES") + index(line,"CONECT") /= 0) then
          if (line(1:6) /= "CONECT") then
            ncomments = ncomments + 1
            tmp_comments(ncomments) = "*"//line(:80)
          end if
          cycle
        end if
        if (line(1:3) == "TER" .and. natoms > 0) then
          nbreaks = nbreaks + 1
          if (nbreaks > 400) then
            call mopend("Maximum number (400) of breaks (lines starting with TER) exceeded")
            return
          end if 
          line1 = txtatm(natoms)(22:)
          do i = natoms, 1, -1
            j = i
            if (labels(i) /= 1) exit
            if (txtatm(i)(22:) /= line1(:5)) then
              j = j + 1
              exit
            end if
          end do
          break_coords(:,nbreaks) = geo(:,i)
          breaks(nbreaks) = i
          lchain = .true.
        end if
        if (line(1:3) == "END") exit
        if (l_pdb) then
          if (first) then
            first = .false.
            write(iw,'(//10x,a,/)')" Keyword PDB present, but the following lines are not in PDB format"
          end if
          write(iw,'(a)')trim(line)
        end if
        cycle
      endif
      if (letter /= " " .and. index(letter, line(22:22)) == 0) cycle
      if (line(17:17) /= " ") then
        if (typea == " ") then
          write (iw,*) nline, " lines have been read in"
          write (iw,*) " Faulty line:", line
          call mopend ("The data set contains alternative location" // &
         & " indicators.  Keyword ALT_A must be used")
           call web_message(iw,"ALT_A.html")
         return
        end if
        if (line(17:17) /= typea) cycle
      end if
      natoms = natoms + 1
!   CLEAN THE INPUT DATA
      call upcase (line, 80)
      icomma = Ichar (comma)
      do i = 1, 80
        khar = Ichar (line(i:i))
        if (khar == icomma) then
          line(i:i) = space
        end if
      end do
      if (line(14:14) == "O" .or. line(14:14) == "H") n_water = n_water + 1
      if (line(14:14) /= "H") then
!
! If HETATM and chain is zero, then use the current chain letter
!
        if (line(:6) == "HETATM" .and. line(22:22) == "0") line(22:22) = "A"
        current_res = nint(reada(line,23))
        if (line(27:27) /= " ") txtmax = 27
        if (lchain .and. (line(:6) /= "HETATM" .or. line(22:22) /= " " .or. ch /= " " &
          .or. natoms < 100 .or. last_atom)) then
          lchain = .false.
          if (new_key_1 /= " ") then
             do i = 23, 25
                if (line(i:i) /= " ") exit
              end do
            new_key_1 = trim(new_key_1)//" "//line(i:26)//line(22:22)
            new_key_2 = trim(new_key_2)//line(22:22)
            previous_res = 2000
          else
             do i = 23, 25
                if (line(i:i) /= " ") exit
              end do
            new_key_1 = " START_RES=("//line(i:26)//line(22:22)
            new_key_2 = " CHAINS=("//line(22:22)
          end if
        else         
          if (current_res - previous_res > 1 .or. (natoms > 2 .and. line(22:22) /= ch)) then
            if (line(:6) /= "HETATM" .or. line(22:22) /= " ") then
              do i = 23, 25
                if (line(i:i) /= " ") exit
              end do
              if (line(22:22) /= ch) then
!
!  Letter changed, so make a break.  There was a fault(?) in the original PDB file
!
                if (n_water /= natoms - old_natoms) new_key_1 = trim(new_key_1)//" "//line(i:26)//line(22:22)
                n_water = 0
                old_natoms = natoms 
                new_key_2 = trim(new_key_2)//line(22:22)
              else
                if (n_water /= natoms - old_natoms) new_key_1 = trim(new_key_1)//"-"//line(i:26)
                n_water = 0
                old_natoms = natoms 
              end if
            end if
          end if          
        end if
        previous_res = current_res   
        ch = line(22:22)      
      end if
      last_atom = (line(:4) == "ATOM")
!
!  Format of txtatm: "( nnnnn rrr mmmm)abc"
!  where: nnnnn = atom number (integer)
!         rrr   = residue three-letter abbreviation (text)
!         mmmm  = residue number (integer)
!         a     = A space or an integer, e.g., 1, 2 or 3, etc.
!         b     = A space or an integer, e.g., 1, 2 or 3, etc.
!         c     = A space or a Greek letter, e.g. a: alpha, g: gamma, etc.
!
!  Check that the residue number is a legal quantity
!
       if (line(18:20) == "***") line(18:20) = "UNK"
       do i = 23, 26
         if (line(i:i) /= " ") then
           if (line(i:i) > "9" .or. line(i:i) < "0") then
             if (line(i:i) /= "-") then
               if (numerr < 10) then
                 if (numerr == 0) write(iw,'(10x,a,/)')"THE RESIDUE SEQUENCE NUMBER IN PDB LINE(S): "
                 write(iw,'(10x,a)')""""//trim(line)//""""
                 numerr = numerr + 1
                 exit
               else if (numerr < 30) then
                 write(iw,'(/10x,a,/)')"Remaining errors suppressed"
               else
                 numerr = 100
               end if
                write(iw,'(10x,a)')"CONTAINS ONE OR MORE NON-NUMERIC CHARACTERS"
                call mopend("ERROR DETECTED IN RESIDUE SEQUENCE NUMBER IN A PDB LINE")
               return
             end if
           end if
         end if
       end do
       txtatm(natoms) = line(1:27)
!
!  If the user wants to keep the original text, use txtatm1
!
      txtatm1(natoms) = line(1:27)
!
      leadsp = .true.
      nvalue = 0
      do i = 1, 80
        if (leadsp .and. line(i:i) /= space) then
          nvalue = nvalue + 1
        end if
        leadsp = (line(i:i) == space)
      end do

!
! ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
!

      if( line(1:4) == "ATOM" .or. line(1:6) == "HETATM" ) then
!
! Is there an element field present?
!
        ele = line(77:78)
        if (ele(1:1) /= " " .and. ele(2:2) == " ") ele = " "
        if( scan( " 0123456789.", ele(1:1) ) /= 0 ) then
          ele(1:1) = ele(2:2)
          ele(2:2) = " "
          if( scan( " 0123456789.", ele(1:1) ) /= 0 ) then ! no
!
! Use 13-14 field
!
            ele = line(13:14)
!
! Special case - old insight format
!
            if ( line(1:4) == "ATOM" .and. ele(1:1) == "H" ) then
            ele(2:2) = " "
            else if( scan( " 0123456789", ele(1:1) ) /= 0 ) then
            ele(1:1) = ele(2:2)
            ele(2:2) = " "
            endif
            
            if( scan( " 0123456789", ele(2:2) ) /= 0 ) then
              ele(2:2) = " "
            endif
          endif
        endif
      else
        natoms = natoms -1
        cycle
      endif

      do i = defined_elements+1, defined_elements + new_elements
        if (ele == element(i)) go to 1010
      end do

      do i = 1, defined_elements
        if (ele == element(i)) go to 1010
      end do

      natoms = natoms - 1

      do i = 1, npdb
        if (ele == txtpdb(i)) cycle outer_loop
      end do
      npdb = npdb + 1
      txtpdb(npdb) = ele
      if (ele(2:2) == " ") then
        ntxt_loc(npdb) = 1
      end if
      write (iw, "('  UNRECOGNIZED SPECIES: (',A,')' ,A,' ON LINE',I6)") ele &
     & (1:ntxt_loc(npdb)), ele(ntxt_loc(npdb)+1:3), nline
      cycle
1010  label = ielem(i)
      if (label == 0 .or. label == 99) then
        natoms = natoms - 1
      else
 !
 ! ALL O.K.
 !
        labels(natoms) = label
        if (label /= 99) then
          numat = numat + 1
          atmass(numat) = ams(label)
        end if
        geo(1, natoms) = reada (line(31:38), 1)
        geo(2, natoms) = reada (line(39:46), 1)
        geo(3, natoms) = reada (line(47:54), 1)
        lopt(1, natoms) = 1
        lopt(2, natoms) = 1
        lopt(3, natoms) = 1
        na(natoms) = 0
        nb(natoms) = 0
        nc(natoms) = 0
      end if
    end do outer_loop
!
!  Ensure that all atoms have a chain-letter
!
    l_letter = .false.
    do i = 1, natoms
      if (txtatm(i)(22:22) /= " ") l_letter(ichar(txtatm(i)(22:22)) - ichar("A") + 1) = .true.
    end do
    do i = 1, 26
      if (.not. l_letter(i)) exit
    end do
    if (i <= 26) then
      j = i - 1
      do i = 1, natoms
        if (txtatm(i)(22:22) == " ") txtatm(i)(22:22) = char(ichar("A") + j)
      end do
    end if        
1020 if (npdb /= 0) then
      write (iw,*) " THE SPECIES THAT WERE NOT RECOGNIZED CAN BE"
      write (iw,*) " RECOGNIZED BY USING THE FOLLOWING KEYWORD"
      write (iw,*)
      commas (npdb) = ")"
      write (iw,*) " PDB(", (txtpdb(i) (1:ntxt_loc(i)), ":Z", commas(i), i=1, &
     & npdb)
      write (iw,*)
      write (iw,*) " WHERE 'Z' IS/ARE ATOMIC NUMBERS"

      call mopend( "UNRECOGNIZED ELEMENT ")
    end if
    txtatm(:numat)(17:17)  = " " !  Delete structural disorder
    txtatm1(:numat)(17:17) = " " !  Delete structural disorder
    if (title == '    NULL  ')  title = " "
    if (allocated(all_comments)) deallocate(all_comments)
    allocate(all_comments(ncomments + 100))
    all_comments(:ncomments) = tmp_comments(:ncomments)
    deallocate(tmp_comments)
    if (index(keywrd,"START_RES") == 0 .and. new_key_1(12:) /= "(" ) then
      j = 0
      line1 = new_key_1
      new_key_1 = " "
      ii = min(399, len_trim(line1))
      do i = 2, ii 
        if (line1(i:i + 1) == "  ") cycle
        if (line1(i - 1:i) == "- ") cycle
        j = j + 1
        new_key_1(j:j) = line1(i:i)
      end do
      if (new_key_1(12:12) == " ") new_key_1 = new_key_1(:11)//new_key_1(13:j)
      keywrd = " "//trim(new_key_1)//") "//trim(keywrd)
    end if
    if (index(keywrd,"CHAINS") == 0 .and. new_key_2(9:) /= "(") keywrd = " "//trim(new_key_2)//") "//trim(keywrd)
    if (line /= " ") then
!
! Dummy read to end of PDB file
!
      do
        read(ir,'(a)',iostat = i) line
        if (i /= 0) exit
        if (line == " ") exit
      end do
    end if
    maxtxt = txtmax
    pdb_label = .true.
 !
 !  SWITCH:   Coordinates are in Cartesian.  Should they be internal?
 !
    lxyz = (Index (keywrd, " INT") == 0)
    do i = 1, natoms
      do j = 1, 3
        coord(j, i) = geo(j, i)
      end do
    end do
    if (lxyz) then
      na1 = 99
      return
    end if
    degree = 57.29577951308232d0
    call xyzint (coord, natoms, na, nb, nc, degree, geo)
 !
 !  UNCONDITIONALLY SET FLAGS FOR INTERNAL COORDINATES
 !
    do i = 1, 3
      do j = i, 3
        lopt(j, i) = 0
      end do
    end do
    if (natoms > 2) then
      if (Abs (geo(2, 3)-180.d0) < 1.d-4 .or. Abs (geo(2, 3)) < 1.d-4) then
        write (iw, "(A)") " DUE TO PROGRAM BUG, THE FIRST THREE " // &
                        & "ATOMS MUST NOT LIE IN A STRAIGHT LINE."
        call mopend ("FIRST 3 ATOMS MUST NOT LIE IN A LINE.")
        return
      end if
    end if
 !
 !  COORDINATES ARE INTERNAL, BUT IN DEGREES.  CONVERT TO RADIANS
 !
    degree = 1.7453292519943d-02
    do j = 2, 3
      do i = 1, natoms
        geo(j, i) = geo(j, i) * degree
      end do
    end do
    na1 = 0
    if (index(keywrd,"START_RES") == 0) then
      j = 0
      do i = 1, len_trim(new_key_1)
        if (new_key_1(i:i + 1) == " ") cycle
      end do
      keywrd = " "//trim(new_key_1)//") "//trim(keywrd)
    end if
    if (index(keywrd,"CHAINS") == 0) then
    end if
    return
end subroutine getpdb
