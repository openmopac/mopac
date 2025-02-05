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

      subroutine getgeo(iread, labels, geo, xyz, lopt, na, nb, nc, int)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
      use parameters_C, only : ams
!
      use molkst_C, only : natoms, keywrd, numat, maxtxt, line, moperr, &
        numcal, numcal0, id, units, Angstroms, arc_hof_1, arc_hof_2, keywrd_txt, pdb_label
!
      use chanel_C, only : iw, ir, input_fn, end_fn, iend
!
      use common_arrays_C, only :atmass, simbol, txtatm, na_store, nat, l_atom, p
!
      USE maps_C, ONLY: react
!
      use funcon_C, only : a0
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iread
      integer , intent(out) :: labels(*)
      integer , intent(out) :: lopt(3,*)
      integer, intent(out)  :: na(*)
      integer, intent(out)  :: nb(*)
      integer, intent(out)  :: nc(*)
      double precision, intent(out)  :: geo(3,*), xyz(3,*)
      logical :: int, lmop, solid, exists, get_q

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(40) :: istart
      integer :: i, icapa, icapz, iserr, k, icomma, khar, nvalue, label, j, ndmy, &
      jj, ltl, max_atoms, ii, ios
      double precision :: weight, real, sum
      logical :: lxyz, velo, leadsp, ircdrc, saddle, mini, l_gaussian
      character , dimension(107) :: elemnt*2
      character :: space, nine, zero, comma, string*200, ele*2, no
      double precision, external :: reada
      save elemnt, space, nine, zero, comma
!-----------------------------------------------
!***********************************************************************
!
!   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
!          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
!
!  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
!             AMS    = DEFAULT ATOMIC MASSES.
!
! ON OUTPUT LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
!           GEO    = INTERNAL COORDINATES, IN ANGSTROMS, AND DEGREES.
!                    OR CARTESIAN COORDINATES, DEPENDING ON WHETHER
!                    KEYWORD ' XYZ' IS PRESENT
!           LOPT   = INTEGER ARRAY, A '1' MEANS OPTIMIZE THIS PARAMETER,
!                    '0' MEANS DO NOT OPTIMIZE, AND A '-1' LABELS THE
!                    REACTION COORDINATE.
!           NA     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           NB     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           NC     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           ATMASS = ATOMIC MASSES OF ATOMS.
!***********************************************************************
      data (elemnt(i),i=1,107)/ 'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F'&
        , 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', &
        'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', &
        'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', &
        'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', &
        'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM'&
        , 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL'&
        , 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP'&
        , 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', '+T', '-T', 'CB', '++', '+', '--'&
        , '-', 'TV'/
      data comma, space, nine, zero/ ',', ' ', '9', '0'/
      ircdrc = index(keywrd,' IRC') + index(keywrd,' DRC') /= 0
      if (index(keywrd," FORCETS") /= 0) ircdrc = .false.
      icapa = ichar('A')
      icapz = ichar('Z')
      get_q = (index(keywrd, " 0SCF") > 0 .and. index(keywrd, " HTML") > 0 .and. index(keywrd_txt, " GEO_REF") == 0)
      if (get_q) then
        if (allocated(p)) deallocate(p)
        i = size(atmass)
        allocate(p(i))
        p = -1.d3
      end if
      saddle = (index(keywrd, "SADDLE") > 0)
      lxyz = (index(keywrd, " XYZ") > 0 .or. saddle .or. (index(keywrd, " LOCATE-TS") /= 0))
      int = (index(keywrd, " INT ") > 0)
      velo = (index(keywrd,' VELO') > 0)
      lmop = (Index (keywrd, " MOPAC") /= 0)
      if (lmop) then
!
! Don't set lmop if " MOPAC" is inside a quoted string
! 
        i = 0
        j = Index (keywrd, " MOPAC") 
        do k = j, len_trim(keywrd)
          if (keywrd(k:k) == '"') i = i + 1
        end do
        lmop = (mod(i, 2) == 0)
      end if
      max_atoms = size(txtatm)
      maxtxt = 0
      istart(1:1) = shape(simbol)
      if (istart(1) < natoms*3) then
        natoms = 0
        numat = 0
        return
      end if
      simbol(:natoms*3) = '---------'
      natoms = 0
      numat = 0
      nvalue = 0
      label = 0
      l_gaussian = .false.
      weight = 0.d0
     ! call dbreak()
      iserr = 0
      mini = (index(keywrd, " MINI") /= 0)
!
!  Set distance units undefined, but with a default of Angstroms
!
      units = .false.
      Angstroms = .true.
      if (index(keywrd, " A0 ") /= 0) then
        units = .true.
        Angstroms = .false.
      else  if (index(keywrd, " Ang") /= 0) then
        units = .true.
      end if
      ii = 0
write(*,*) "before main loop of getgeo"
   20 continue
      read (iread, '(A241)', iostat=ios, end=120, err=210) line
      if (line == '$coord') go to 20
      if (line == '$end') go to 20
      if (line(1:1) == '*') go to 20
      if (line == ' ') then
        if(natoms == 0 .and. numcal == 1+numcal0) then
!
!  Check:  Is this an ARC file?
!
          numcal = 2+numcal0
          rewind (iread)
          sum = 0.d0
          do i = 1, 10000
            read (iread, '(A)', iostat=ios, end=120, err=210) line
            if (index(line, "HEAT OF FORMATION") > 0) sum = reada(line,20)
            if (index(line, "FINAL GEOMETRY OBTAINED") > 0) exit
            if (index(line, "GEOMETRY IN CARTESIAN COORDINATE") > 0) exit
            if (index(line, "GEOMETRY IN MOPAC Z-MATRIX") > 0) exit
          end do
!
!  Yes, it's an arc file
!
          if (abs(sum) > 1.d-4) then
            if (abs(arc_hof_1) < 1.d-4) then
              arc_hof_1 = sum
            else
              arc_hof_2 = sum
            end if
          end if
          ii = 0
          if (index(line, "FINAL GEOMETRY OBTAINED") > 0) then
!
! Read in keywords, title and comment
!
            call gettxt
            saddle = (index(keywrd, "SADDLE") > 0)
            lxyz = (index(keywrd, " XYZ") > 0 .or. saddle .or. (index(keywrd, " LOCATE-TS") /= 0))
            int = (index(keywrd, " INT ") > 0)
            velo = (index(keywrd,' VELO') > 0)
            lmop = (Index (keywrd, " MOPAC") /= 0)
            read (iread, '(A)', iostat=ios, end=120, err=210) line
            ii = 3
          else
            natoms = -3
            return
          end if
        end if
        if (ii == 0) goto 120
        ii = 0
      end if
!
!   If two quotation marks are side-by-side, force a space in between them
!
      do
        i = index(line, '""')
        if (i /= 0) then
          string = trim(line)
          line = string(:i)//" "//trim(string(i + 1:))
        else
          exit
        end if
      end do
      ltl = len_trim(line)
      icomma = ichar(comma)
      do i = 1, ltl
        khar = ichar(line(i:i))
        if (khar == icomma .or. khar == 9) line(i:i) = space
      end do
      if (natoms == max_atoms) then
        write(iw,"(//10x,a)")" Maximum number of atoms exceeded"
        write(iw,"(10x,a,i6)")" Maximum allowed:", max_atoms
        call mopend("Maximum number of atoms exceeded")
        return
      end if
      if (natoms < 10) then
!
!  Check for anything that might suggest that the data-set is a PDB
!
        if (index(line(:4),"ATOM") + index(line(:6),"HETATM") + index(line,"TITLE ") + index(line,"HEADER") + &
                                         index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
                  index(line,"HELIX ") + index(line,"SHEET ") + index(line,"REMARK") + index(line,"USER  ") + &
                  index(line,"EXPDTA") + index(line,"AUTHOR") + index(line,"REVDAT") + index(line,"JRNL  ") + &
                  index(line,"DBREF ") + index(line,"SEQRES") + index(line,"HET   ") + index(line,"HETNAM") + &
                  index(line,"LINK  ") + index(line,"CRYST1") + index(line,"SCALE" ) + index(line,"ORIGX" ) + &
                  index(line,"FORMUL") + index(line,"SEQRES") + index(line,"CONECT") /= 0 .and. &
                  index(line(:6), "(") == 0) then
            natoms = -2
            return
        end if
      end if
!
!   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
!
      i = index(line,'(')
      if (i /= 0) then
!
!  YES, ELEMENT IS LABELLED.
!
        k = index(line,')')
        if (k == 0) then
          write(iw,"(a,i5,a)")" Atom",natoms," has an opening parenthesis but no closing parenthesis"
          write(iw,"(a)")" Line :'"//line(:len_trim(line))//"'"
          write (iw, '(/,A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS'
          atmass(1) = -1.D0
          nat(1) = 0
          l_atom = .true.
          call geout (iw)
          call mopend ('GEOMETRY IS FAULTY')
          return
        end if
        if (line(k:k) == ")") then
          jj = k - 1
        else
          jj = k
        end if
        txtatm(natoms + 1) = line(i + 1:jj)
        if (maxtxt == 27 .and. k - i - 1 /= 27) then
          if (index(keywrd, " GEO-OK") == 0) then
            no = char(Nint(log10(natoms + 1.1)) + ichar("1"))
            write(iw,'(//10x,a,i'//no//',a,i3,a)') "Atom label length of 27 detected, but atom ", natoms + 1, " has a label of", &
            k - i - 1, " characters."
            write(iw,'(/,a)')"Faulty line = '"//trim(line)//"'"
            write (iw,'(/10x,a)') "If this is intended, add keyword ""GEO-OK"""
            write(line,'(a,i'//no//')')"Fault detected in atom number ", natoms + 1
            call mopend (trim(line))
            stop
          end if
        end if
        maxtxt = max(maxtxt,k - i - 1)
        if (maxtxt > 38) then
          if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + &
            index(line,"HEADER") + index(line,"ANISOU") + index(line,"COMPND") + &
            index(line,"SOURCE") + index(line,"KEYWDS") + index(line,"USER ")  + &
            index(line,"HELIX") + index(line,"SHEET") /= 0)  goto 70
          write(iw,"(a)")" Atom labels must not exceed 38 characters"
          string = " "
          write(iw,"(a)")string(:i)//"                1         2         3         4"
          write(iw,"(a)")string(:i)//"       1234567890123456789012345678901234567890"
          write(iw,"(a)")"Line :"""//line(:len_trim(line))//""""
          write (iw, '(/,A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS'
          atmass(1) = -1.D0
          nat(1) = 0
          call geout (iw)
          call mopend ('GEOMETRY IS FAULTY')
          return
        end if
        string = line(1:i-1)//line(k+1:)
        line = string
      else
        txtatm(natoms + 1) = ' '
      end if
!   CLEAN THE INPUT DATA
      call upcase (line, ltl)
      if (index(line, " NEXT") /= 0) goto 120
!
!   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      istart(:10) = len_trim(line) + 1
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER AND STORE IN ISTART
      leadsp = .TRUE.
      nvalue = 0
      i = 0
      do jj = 1, ltl
        i = i + 1
        if (leadsp .and. line(i:i) /= space) then
          nvalue = nvalue + 1
          istart(nvalue) = i
          if (i > 2) then
            if (line(i - 1:i - 1) ==  "-" .and. line(i - 2:i - 2) /=  " ") istart(nvalue) = i - 1
          end if
          if (line(i:i) == '"') then
!
!  Run to other end of quoted text
!
            do j = 1, 40
              i = i + 1
              if (line(i:i) == '"') exit
            end do
          end if
        end if
!
!  set leadsp true if a space is detected, or if a "-" sign is found and it's part of a new number
!
        leadsp = (line(i:i) == space .or. (i /= istart(max(nvalue,1)) .and. line(i:i) ==  "-"))
        if (i > 1) then
          if ((line(i-1:i) == "D-" .or. line(i-1:i) == "E-")) leadsp = .false.
        end if
      end do
      l_gaussian = (nvalue == 7)
      if (nvalue == 4) then !  Cartesian coordinates without optimization flags
!
! Check: is the fourth datum only 1 or 2 characters long?
!
        i = index(line(istart(4):), " ") + istart(4) - 2
        k = 0
        if (i > istart(4) + 1) k = 1
!
! Is is alphabetic?
!
        do j = istart(4), i
          if (line(j:j) < "A" .or. line(j:j) > "Z") k = 1
        end do
!
! Is the first datum definitely not pure alphabetic?
!
        i = index(line(istart(1):), " ") + istart(1) - 2
        do j = istart(1), i
          if (line(j:j) >= "A" .and. line(j:j) <= "Z") exit
        end do
        if (j > i .and. k == 0) then
          natoms = natoms + 1
          geo(1,natoms) = reada(line,istart(1))
          geo(2,natoms) = reada(line,istart(2))
          geo(3,natoms) = reada(line,istart(3))
          lopt(1,natoms) = 1
          lopt(2,natoms) = 1
          lopt(3,natoms) = 1
          na(natoms) = 0
          nb(natoms) = 0
          nc(natoms) = 0
          do i = 1, 107
            if (line(istart(4): istart(4) + 1) /= elemnt(i)) cycle
            labels(natoms) = i
            atmass(natoms) = ams(i)
            exit
          end do
          if (.not. units) then
!
!  distance units were not defined, so define them
!
            Angstroms = .false.
            units = .true.
          end if
          goto 20
        end if
      end if
!
! ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
!
      weight = 0.D0
      string = line(istart(1):istart(2)-1)
      if (string == "+3") string = "+T"
      if (string == "-3") string = "-T"
      if (string(1:1) >= zero .and. string(1:1) <= nine) then
!  ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
        label = nint(reada(string,1))
        if (label == 0) go to 120
        if (label < 0 .or. label > 107) then
          write (iw, '(''  ILLEGAL ATOMIC NUMBER'')')
          go to 210
        end if
        go to 70
      end if
!  ATOMIC SYMBOL USED
      real = abs(reada(string,1))
      if (real < 1.D-15) then
!   NO ISOTOPE
        ele = string(1:2)
      else
        weight = real
        if (string(2:2) >= zero .and. string(2:2) <= nine) then
          ele = string(1:1)
        else
          ele = string(1:2)
        end if
      end if
!   CHECK FOR ERROR IN ATOMIC SYMBOL
      if (ele(1:1)=='-' .and. ele(2:2)/='-') ele(2:2) = ' '
      do i = 1, 107
        if (ele /= elemnt(i)) cycle
        label = i
        go to 70
      end do
      if (ele(1:1) == 'X') then
        label = 99
        go to 70
      end if
      if (ele == "D ") then
        label = 1
        weight = 2.014d0
        go to 70
      else if (ele == "T ") then
        label = 1
        weight = 3.016d0
        go to 70
      else
         if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + index(line,"HEADER") + &
        index(line,"ANISOU") + index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
        index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") + index(line,"USER ")  + &
        index(line, "SEQRES") + index(line,"ENDMDL") /= 0) goto 70
        if (trim(line) == " *                    *") return
        write (iw, '(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')') ele
        write(iw,'(/,"  Faulty line: """,a,"""",/)')trim(line)
        go to 210
      end if
!
! ALL O.K.
!
 70   continue
      natoms = natoms + 1
      nb(natoms) = 0
      nc(natoms) = 0
      if (label /= 99) numat = numat + 1
      if (weight /= 0.D0) then
        atmass(numat) = weight
      else
        if (label /= 99) atmass(numat) = ams(label)
      end if
      labels(natoms) = label
      if (nvalue == 4) then !  Cartesian coordinates without optimization flags
        geo(1,natoms) = reada(line,istart(2))
        geo(2,natoms) = reada(line,istart(3))
        geo(3,natoms) = reada(line,istart(4))
        lopt(1,natoms) = 1
        lopt(2,natoms) = 1
        lopt(3,natoms) = 1
      else
        geo(1,natoms) = reada(line,istart(2))
        geo(2,natoms) = reada(line,istart(4))
        geo(3,natoms) = reada(line,istart(6))
        lopt(1,natoms) = nint(reada(line,istart(3)))
        lopt(2,natoms) = nint(reada(line,istart(5)))
        lopt(3,natoms) = nint(reada(line,istart(7)))
        if (nvalue == 8 .and. get_q) p(natoms) = reada(line,istart(8))
        if (nvalue == 11 .and. get_q) p(natoms) = reada(line,istart(11))
        do i = 3, 7, 2
          if (.not.(ichar(line(istart(i):istart(i)))>=icapa .and. ichar(line(&
            istart(i):istart(i)))<=icapz .and. natoms>1)) cycle
          iserr = 1
        end do
      end if
      pdb_label = (maxtxt > 25)
      if (line(istart(10):istart(10)) == '"' .or. line(istart(9):istart(9)) == '"' .or. &
        line(istart(8):istart(8)) == '"') then
        call  l_control("CONTROL_NABC_in_PDB", len_trim("CONTROL_NABC_in_PDB"), 1)
        end if
      string = trim(line)
      if (line(istart(10):istart(10)) == '"') then
        call txt_to_atom_no(line, istart(10), .false.)
      end if
      if (line(istart(9):istart(9)) == '"')   call txt_to_atom_no(line, istart(9), .false.)
      if (line(istart(8):istart(8)) == '"')   call txt_to_atom_no(line, istart(8), .false.)
      if (moperr) then
        write(iw,'(/10x,a, i5)')" Error detected in definition of atom no.:", natoms
        write(iw,'(/,a)')" Text of faulty atom: '"//trim(elemnt(labels(natoms)))//&
        & "("//txtatm(natoms)//") "//trim(string(istart(2):))//"'"
        return
      end if
      j = len_trim(line)
      if (ltl /= j) then
        i = index(line(istart(8):), " ") + istart(8)
        nvalue = 8
        leadsp = .true.
        do i = i, j
          if (leadsp .and. line(i:i) /= space) then
            nvalue = nvalue + 1
            istart(nvalue) = i
          end if
          leadsp = (line(i:i) == " ")
        end do
      end if
      sum = reada(line,istart(8))
      i = index(line(istart(8):), " ") + istart(8) ! Find end of 8'th datum
      if (index(line(istart(8):i), ".") /= 0) sum = 0.D0 ! if 8'th datum contains a decimal, then it's not a connectivity
      na(natoms) = nint(sum)
      if (natoms == 1) na(1) = 0
      if (lmop .and. natoms == 2) then
        na(2) = 1
      else if (lmop .and. natoms == 3 .and. na(3) == 0) then
        na(3) = 2
        nb(3) = 1
        geo(2,3) = geo(2,3)*1.7453292519943D-02
        geo(3,3) = 0.D0
        lopt(3,3) = 0
      else if (na(natoms) > 0) then
        nb(natoms) = nint(reada(line,istart(9)))
        nc(natoms) = nint(reada(line,istart(10)))
        geo(2:3,natoms) = geo(2:3,natoms)*1.7453292519943D-02
      end if

!
!  SPECIAL CASE OF USERS FORGETTING TO ADD DIHEDRAL DATA FOR ATOM 3
!
      if (natoms == 3) then
        if (lopt(1,3) /= 2 .and. lopt(2,3) /= 2 .and. lopt(3,3) == 2) then
          na(3) = 1
          nb(3) = 2
          geo(3,3) = 0.D0
          lopt(3,3) = 0
        else if (lopt(3,3)==1 .and. line(istart(6):istart(6) + 1) == "2 ") then
          na(3) = 2
          nb(3) = 1
          geo(3,3) = 0.D0
          geo(2,3) = geo(2,3)*1.7453292519943D-02
          lopt(3,3) = 0
        end if
      end if
      if ( .not. mini .and. .not. l_gaussian) then
        if ((lopt(1,natoms) > 1 .or. lopt(2,natoms) > 1 .or. lopt(3,natoms) > 1) .and. natoms > 1) then
          if (nvalue == 5 .and. lopt(2,natoms) > 1) &
            l_gaussian = (line(istart(5):istart(5)) >= "A" .and. line(istart(5):istart(5)) <= "Z")
          if (nvalue == 3 .and. lopt(1,natoms) > 1) &
            l_gaussian = (line(istart(3):istart(3)) >= "A" .and. line(istart(3):istart(3)) <= "Z")
          if (l_gaussian) then
!
!  Geometry contins indications that it is in Gaussian format, so return and try reading it again.
!
            natoms  = -1
            return
          end if
          iserr = 1
          write(iw,'(/10x,a,i5)')"Faulty atom:", natoms
          write(iw,'(/10x,a)')"Faulty line: """//trim(line)//""""
          call mopend("Unless MINI is used, optimization flags must be 1, 0, or -1")
          numcal = 2+numcal0
          if ((lopt(1,natoms) > 10 .or. lopt(2,natoms) > 10 .or. lopt(3,natoms) > 10) .and. natoms > 1) &
            write(iw,'(/10x,a)')" If the geometry is in Gaussian format, add keyword ""AIGIN"" and re-run"
          return
        end if
      else
         i = max(lopt(1,natoms),lopt(2,natoms),lopt(3,natoms))
         j = min(lopt(1,natoms),lopt(2,natoms),lopt(3,natoms))
        if (max(i, -j) > 1) then
          if (i /= -j) then
            do k = 1, 3
              if (lopt(k,natoms) /= 0) lopt(k,natoms) = sign(2,lopt(k,natoms))
            end do
          end if
        end if
      end if
      if (iserr == 1) then
      if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + index(line,"HEADER") + &
        index(line,"ANISOU") + index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
        index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") + index(line,"USER ")  + &
        index(line, "SEQRES") /= 0) then
!
!  Geometry is definitely PDB
!
        natoms = -2
        return
      end if
!
!  MUST BE GAUSSIAN GEOMETRY INPUT
!
        do i = 2, natoms
          do k = 1, 3
            j = nint(geo(k,i))
            if (abs(geo(k,i)-j) <= 1.D-5) cycle
!
!   GEOMETRY CANNOT BE GAUSSIAN
!
           write (iw, '(A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS'
            atmass(1) = -1.D0
            nat(1) = 0
            call geout (iw)
            call mopend ('GEOMETRY IS FAULTY')
            return
          end do
        end do
        natoms = -1
        return
      end if
write(*,*) "loop iter", natoms
if(natoms > 5583) stop
      go to 20
!***********************************************************************
! ALL DATA READ IN, CLEAN UP AND RETURN
!***********************************************************************
  120 continue
write(*,*) "end of loop"
stop
      if (natoms == 0) then
        if (numcal == 1+numcal0) call mopend (' Error detected while reading geometry')
        return
      end if
      if ( .not. Angstroms) then
!
!  Convert from A0 to Angstroms
!
        do i = 1, natoms
          geo(1,i) = geo(1,i)*a0
          if (na(i) == 0) then
            geo(2,i) = geo(2,i)*a0
            geo(3,i) = geo(3,i)*a0
          end if
        end do
      end if
      if (maxtxt > 0) then
        do i = 1, natoms
          j = max(1, len_trim(txtatm(i)))
          if (j /= maxtxt) then           !  Pad out text with blanks
            txtatm(i)(j + 1:) = " "       !  (Removes rubbish like "null" characters)
          end if
        end do
      end if
      na(1) = 0
      nb(1) = 0
      nc(1) = 0
      nb(2) = 0
      nc(2) = 0
      nc(3) = 0
!
!     READ IN VELOCITY VECTOR, IF PRESENT
!
      if (velo) then
        j = 0
        do i = 1, natoms
          if (na(i) /= 0) j = 1
        end do
        if ((j /= 0 .or. int .and. index(keywrd,' LET') == 0) .and. index(keywrd,' 0SCF') == 0) then
          call mopend ('COORDINATES MUST BE CARTESIAN WHEN VELOCITY VECTOR IS USED.')
          return
        end if
        if (allocated(react)) deallocate(react)
        allocate(react(3*natoms))
        do i = 1, natoms
          read (iread, '(A)') line
          call nuchar (line, len_trim(line), react((i-1)*3+1), ndmy)
          if (ndmy == 3) cycle
          call mopend ('THERE MUST BE EXACTLY THREE VELOCITY DATA PER LINE.')
          return
        end do
      end if
      if (ircdrc .and. index(keywrd,' 0SCF') == 0) then
        if (numat /= natoms) then
          call mopend ('Only real atoms are allowed in IRC and DRC calculations.')
          return
        end if
        call gmetry(geo, xyz)
        call xyzint (xyz, numat, na, nb, nc, 1.D0, geo)
        geo(:,1:numat) = xyz(:,1:numat)
        na_store = na(:numat) ! Store na - it will be used in printing the DRC
        na(:numat) = 0
      end if
      solid = .false.
      do i = 1, natoms
        if (labels(i) <= 0) then
          write (iw, '('' ATOMIC NUMBER OF '',I3,'' ?'')') labels(i)
          if (i == 1) then
            write (iw, '(A)') ' THIS WAS THE FIRST ATOM'
          else
            write (iw, '(A)') &
              '    GEOMETRY UP TO, BUT NOT INCLUDING, THE FAULTY ATOM'
            natoms = i - 1
            call geout (iw)
          end if
          call mopend ('Error in READMO')
          return
        end if
        if (labels(i) == 107) solid = .true.
        if (solid .and. labels(i) < 99) then
          call mopend ('TRANSLATION VECTORS MUST BE AT THE END OF THE DATA SET')
          return
        end if
        if (i == 1 .or. na(i) == 0) cycle
        j = 0
        if (na(i) == nb(i).and. i > 1 .or. &
          (na(i) == nc(i) .or. nb(i) == nc(i)) .and. i > 2 .or. &
          nb(i)*nc(i) == 0 .and. i > 3) j = 1  !  Error condition
        if (na(i) >= i .or. nb(i) >= i .or. nc(i) >= i) then
!
! An atom is being defined using a connectivity involving one or more atoms
! whose positions have not yet been defined.  This is likely to be an error.
!
          j = j + 1
!
!  Test: The atom is defined in internal coordinates, so check that
!  the dependent atoms are defined in Cartesian coordinates
!
          if ((na(i) < i .or. na(max(1,na(i))) == 0) .and. &
              (nb(i) < i .or. na(max(1,nb(i))) == 0) .and. &
              (nc(i) < i .or. na(max(1,nc(i))) == 0) ) j = j - 1
        end if
        if (j == 0) cycle
        j = max(i, na(i), nb(i), nc(i))
        no = char(Nint(log10(j + 1.1)) + ichar("1"))
        write (line, '(" ATOM NUMBER ",I'//no//'," IS ILL-DEFINED")') i
        write (iw, '(//10x, a)') trim(line)
        write(iw,'(/10x," Connectivity of atom ",i'//no//',": NA=",i'//no//',", NB=",i'//no//',", NC=",i'//no//')') &
        i, na(i), nb(i), nc(i)
        if (na(i) > i) then
          if (na(na(i)) /= 0) write(iw,'(/10x,a)')" NA is defined using internal coordinates"
        end if
        if (nb(i) > i) then
          if (na(nb(i)) /= 0) write(iw,'(/10x,a)')" NB is defined using internal coordinates"
        end if
        if (nc(i) > i) then
          if(na(nc(i)) /= 0) write(iw,'(/10x,a)')" NC is defined using internal coordinates"
        end if
          call web_message(iw,"geometry_specification.html")
        if (i == 1) then
          return
        end if
        write (0, '(//10x, a)') trim(line)
        write (iw, '(/,''  GEOMETRY READ IN'',/)')
        nat = 0
        call geout (iw)
        call mopend (trim(line) )
        inquire (file = input_fn, exist = exists)
        if (exists) close(ir, status = 'delete', err = 99)
        inquire (file = end_fn, exist = exists)
        if (exists) then
          open(unit=iend, file=end_fn, iostat=j)
          close(iend, status = 'delete', iostat=j)
        end if
99      stop
      end do
      if (natoms > 0) then
        call gmetry (geo, xyz)
        if (get_q) then
          j = 0
          k = 1
          do i = 1, natoms
            if (labels(i) /= 99) then
              j = j + 1
              p(j) = p(i)
              if (p(j) < -10.d0) k = 0
            end if
          end do
          if (k == 0) deallocate(p)
        end if
      else
        call mopend ("No atoms!")
      end if
      if (moperr) return
!
!  Switch for converting between coordinate systems.
!
      if (lxyz .or. int) then
!
!    Coordinates should all be Cartesian or all be internal.
!    First, unconditionally convert to Cartesian.
!
        k = 0
        do i = 1, natoms
          do j = 1, 3
            k = k + lopt(j,i)
          end do
        end do
        j = na(3)
!
!  Get rid of dummy atoms
!
        numat = 0
        do i = 1, natoms
          if (labels(i) /= 99) then
            numat = numat + 1
            labels(numat) = labels(i)
            txtatm(numat) = txtatm(i)
            lopt(:,numat) = lopt(:,i)
          end if
          geo(:,i) = xyz(:,i)
        end do
        na(:natoms) = 0
!
!   If everything is marked for optimization then unconditionally mark the first
!   three atoms for optimization
!
        if (k >= 3*numat - 6 .and. j /= 0) lopt(:,:min(3, numat)) = 1
        natoms = numat
        if (saddle .or. (index(keywrd, " LOCATE-TS") /= 0)) then
          if (index(keywrd, " LOCATE-TS") == 0) lopt(:,:numat) = 1  ! In a saddle calculation, all parameters must be optimizable.
          na(:natoms) = 0
        end if
      end if
      if (int) then
!
!  Coordinates should be internal. Unconditionally convert to internal
!
        if (id > 0) then
!
!  System is a solid, so move translation vectors to suit atom 1
!
          do i = numat - id + 1, numat
            xyz(:,i) = xyz(:,i) + xyz(:,1)
          end do
        end if
        call xyzint (xyz, numat, na, nb, nc, 1.d0, geo)
!
!  UNCONDITIONALLY SET FLAGS FOR INTERNAL COORDINATES
!
        do i = 1, 3
          lopt(i:3,i) = 0
        end do
      end if
      return
! ERROR CONDITIONS
  210 continue
! gfortran flags all EOF reads past the first one as a 5001 error, and MOPAC doesn't avoid this behavior at the moment
      if (ios == 5001) goto 120
      j = natoms - 1
      write (iw, '('' DATA CURRENTLY READ IN ARE: '',/)')
      do k = 1, j
        if (na(k) > 0) geo(2:3,k) = geo(2:3,k)/1.7453292519943D-02
        write (iw, '(3x,a,2x,3(f10.5,2x,i2,2x),3(i2,1x))') &
        elemnt(labels(k)), (geo(jj,k),lopt(jj,k),jj=1,3), &
        na(k), nb(k), nc(k)
      end do
      natoms = 0
      return
      end subroutine getgeo
      subroutine web_message(channel, txt)
        implicit none
        character (*) ::txt
        integer :: channel
        write(channel,'(/10x,a,/)')"For more information, see: HTTP://OpenMOPAC.net/Manual/"//trim(txt)
        return
      end subroutine web_message
