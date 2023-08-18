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

subroutine pdbout (mode1)
    use molkst_C, only: numat, natoms, ncomments, verson, line, nbreaks, id, &
      maxtxt, keywrd, nelecs, escf, stress, backslash, l_normal_html
    use chanel_C, only: iw, input_fn
    use elemts_C, only: elemnt
    use common_arrays_C, only: txtatm, nat, all_comments, p, labels, &
      breaks, txtatm1, geo, l_atom
    USE parameters_C, only : tore
    implicit none
    integer, intent (in) :: mode1
!
    character :: ele_pdb*2, idate*24, num*1, x*1
    integer :: i, i1, i2, iprt, k, nline
    logical :: ter, html, ter_ok, l_irc_drc, charge
    double precision :: sum
    intrinsic Abs, Char, Ichar
    double precision, allocatable :: q2(:), coord(:,:)
!
    html = (index(keywrd, " HTML") /= 0)
    ter_ok = (index(keywrd, " NOTER") == 0)
    charge = (index(keywrd, " PRTCHAR") /= 0)
    l_irc_drc = (index(keywrd, " IRC") + index(keywrd, " DRC") /= 0)
    ter = .false.
    allocate (q2(numat), coord(3, natoms))
    if (allocated(p)) then
      if (index(keywrd, " 0SCF") /= 0) then
        sum = 0.d0
        i2 = 0
        do i = 1, numat
          sum = sum + p(i)
          if (txtatm(i) /= " ") then
            i2 = i2 + 1
            txtatm(i2) = txtatm(i)
          end if
        end do
        if (sum > 0.d0) then
          nelecs = nelecs - nint(sum)
          q2(:numat) = p(:numat)
        else
          q2(:numat) = 0.d0
        end if
      else
        call chrge (p, q2)
        q2(:numat) = tore(nat(:numat)) - q2(:numat)
      end if
    else
      q2 = 0.d0
    end if
!
! Convert geometry into Cartesian coordinates, preserving dummy atoms
!
    do i = 1, numat
      if (labels(i) == 99) labels(i) = 199
    end do
    call gmetry(geo, coord)
    do i = 1, numat
      if (labels(i) == 199) labels(i) = 99
    end do
    nline = 0
    if (mode1 == 1) then
      iprt = iw
    else
      iprt = Abs (mode1)
    end if
    call fdate (idate)
    i = len_trim(input_fn)
    if (ncomments > 0) then
      if (index(all_comments(1),"HEADER") == 0) then
        line = "HEADER  data-set: "//input_fn(:i - 5)
        if(len_trim(line) > 80) then
          do
            k = index(line, backslash)
            if (k == 0) exit
            i2 = index(line,"data-set:") + 9
            line = line(:i2)//line(k + 1:)
          end do
        end if
        line(81:) = " "
        write (iprt, "(A)") trim(line)
        line = "REMARK  MOPAC, Version: "//verson//" Date: "//idate(5:11)//idate(21:)//idate(11:16)
        write (iprt, "(A)") trim(line)
      end if
      do i = 1, ncomments
        line = all_comments(i)(:7)
        if (index(line,"ATOM  ") + index(line,"HETATM") + index(line,"TITLE ") + index(line,"HEADER") + &
            index(line,"ANISOU") + index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
            index(line,"HELIX ") + index(line,"SHEET ") + index(line,"REMARK") + index(line,"USER  ") + &
            index(line,"EXPDTA") + index(line,"AUTHOR") + index(line,"REVDAT") + index(line,"JRNL  ") + &
            index(line,"DBREF ") + index(line,"SEQRES") + index(line,"HET   ") + index(line,"HETNAM") + &
            index(line,"LINK  ") + index(line,"CRYST1") + index(line,"SCALE" ) + index(line,"ORIGX" ) + &
            index(line,"FORMUL") + index(line,"SEQRES") == 0) cycle
        write(iprt,"(a)")all_comments(i)(2:len_trim(all_comments(i)))
      end do
    else
      line = "HEADER  data-set: "//input_fn(:i - 5)
      if(len_trim(line) > 80) then
        do
          k = index(line, backslash)
          if (k == 0) exit
          i2 = index(line,"data-set:") + 9
          line = line(:i2)//line(k + 1:)
        end do
      end if
      line(81:) = " "
      write (iprt, "(A)") trim(line)
      if (ncomments == 0) then
        line = "REMARK  MOPAC, Version: "//verson
        write (iprt, "(A)") trim(line)
        line = "REMARK  Date: "//idate(5:11)//idate(21:)//idate(11:16)
        write (iprt, "(A)") trim(line)
        if (.not. l_irc_drc) line = "REMARK  Heat of Formation ="
        sum = escf - stress
        if (abs(sum) > 4.99999d-4) then
          i = max(int(log10(abs(sum))), 0)
          if ( sum < 0.d0) i = i + 1
          if (i < 4) then
            num = char(ichar("6") + i)
            write (iprt, "(A, f"//num//".3, a)") trim(line), sum, " Kcal/mol"
          else
            num = char(ichar("6") + i - 10)
            write (iprt, "(A, f1"//num//".3, a)") trim(line), sum, " Kcal/mol"
          end if
        end if
      end if
    end if
    if (maxtxt == 0 .and. txtatm(1) /= " ") then
      txtatm1(:natoms) = txtatm(:natoms)
    end if
    i1 = 0
    i2 = 0
    nbreaks = 1
    do i = 1, natoms - id
      if (labels(i) /= 99) i1 = i1 + 1
      i2 = i2 + 1
!
!   i:  all atoms, real and dummy
!   i1: real atom only
!   i2: all atoms, real and dummy, plus TER.  I.e., PDB record entry number
!

      nline = nline + 1
      if (ter_ok) ter = (i == breaks(nbreaks))
      if (ter) nbreaks = nbreaks + 1
      if (i1 == 0 .and. i /= 1) cycle
      if ( .not. l_atom(max(1,i1))) cycle
      if (elemnt(labels(i)) (1:1) == " " .or. labels(i) == 99) then
        ele_pdb(1:1) = " "
        ele_pdb(2:2) = elemnt(labels(i)) (2:2)
      else
        ele_pdb(1:1) = elemnt(labels(i)) (1:1)
        ele_pdb (2:2) = Char (Ichar ("A")-Ichar ("a")+Ichar (elemnt(labels(i)) (2:2)))
      end if
      if (maxtxt == 27) then
        num = "1"
      else
        num = "2"
      end if
      if (txtatm(i)(15:16) == "**") txtatm(i)(15:16) = "99"
      if (txtatm(i)(23:26) == "****") txtatm(i)(23:26) = "9999"
      x = txtatm(i)(13:13)
      if (x == "X") txtatm(i)(13:13) = " "
      if (txtatm(i)(14:14) /= "X") then
        if (charge) then
          write (iprt, "(a,i5,a,f1"//num//".3,f8.3,f8.3,a,f7.2,f10.3, a2,a)") txtatm(i)(1:6),i2,txtatm(i)(12:maxtxt), &
          & (coord(k, i), k=1, 3), "  1.0", 0.d0, q2(i1), ele_pdb, " "
        else
          write (iprt, "(a,i5,a,f1"//num//".3,f8.3,f8.3,a,f7.2,a, a2,a)") txtatm(i)(1:6),i2,txtatm(i)(12:maxtxt), &
          & (coord(k, i), k=1, 3), "  1.0",q2(i1)*10.d0,"      PROT", ele_pdb, " "
        end if
      else
        write (iprt, "(a,i5,a,f1"//num//".3,f8.3,f8.3,a,f7.2,a, a2,a)") txtatm(i)(1:6),i2,txtatm(i)(12:maxtxt), &
        & (coord(k, i), k=1, 3), "  1.0 ",0.d0,"      PROT", ele_pdb, " "
      end if
      txtatm(i)(13:13) = x
      if (ter) then
        i2 = i2 + 1
        write(iprt, "(a,i8,a)")"TER",i2,"      "//txtatm(i)(18:26)  !  Write out the TERMINAL marker
      end if
    end do
    write(iprt, "(a)")"END"
    if (html) then
      if (l_normal_html) then
        l_normal_html = .false.
        call write_html
      end if
    end if
  end subroutine pdbout
  subroutine write_html
    use chanel_C, only: input_fn
    use molkst_C, only : line, koment, title, keywrd, numcal, geo_ref_name, geo_dat_name, &
      maxtxt, backslash, natoms, id
    use common_arrays_C, only : txtatm, p
    use mozyme_C, only : tyres, tyr
    implicit none
    integer :: nres, i, j, k, nprt, ncol, biggest_res, iprt, icalcn = -1, it
    character :: res_txt(4000)*10, l_res*14, n_res*14, wrt_res*14, num*1, line_1*400
    integer, parameter :: limres = 260
    logical :: exists, l_prt_res, l_geo_ref, l_compare, l_het_only
    double precision, external :: reada
    save :: icalcn
    it = 0
    iprt = 27
    if (icalcn == numcal .and. index(keywrd, " Write_Escf") == 0) return
    icalcn = numcal
    l_compare = (index(keywrd, " COMPARE") /= 0)
    l_prt_res = (index(keywrd," NORJSMOL") == 0)
    line = input_fn(:len_trim(input_fn) - 4)//"html"
    open(unit=iprt, file=trim(line))
    nres = 0
    if (l_prt_res) then
!
!   Identify all residues
!

      l_het_only = .true.
      do i = 1, natoms - id
        if (txtatm(i)(1:4) == "ATOM") l_het_only = .false.
        if (txtatm(i)(18:20) == "HOH") cycle
        if (txtatm(i)(18:20) == "WAT") cycle
        if (txtatm(i)(18:20) == "DOD") cycle
        if (txtatm(i)(18:20) == "SO4") cycle
        j = nint(reada(txtatm(i), 23))
        if (j < 0) then
           write(l_res,"(a3,i4.3,a)")txtatm(i)(18:20),j,":"//txtatm(i)(22:22)
        else
           write(l_res,"(a3,i4.4,a)")txtatm(i)(18:20),j,":"//txtatm(i)(22:22)
        end if
        if (maxtxt == 27) then
          l_res(10:10) = txtatm(i)(27:27)
        else
          l_res(10:10) = " "
        end if
!
! If the first character is a number, change it to "Q"
!
        do k = 1, 6
          if (l_res(k:k) /= " ") exit
        end do
        if (l_res(k:k) >= "0" .and. l_res(k:k) <= "9") l_res(k:k) = "Q"
!
! Change any spaces to "J"
!
        do k = k + 1, 9
          if (l_res(k:k) == " ") then
            if (txtatm(i)(:6) /= "HETATM") l_res(k:k) = "J"
          end if
        end do
          do k = 1, nres
          if (res_txt(k) == l_res .or. l_res(1:7) == "   Q000") exit
        end do
        if (k > nres) then
          nres = nres + 1
          res_txt(nres) = trim(l_res)
          do k = 1, len_trim(l_res)
            if (res_txt(nres)(k:k) == " ") res_txt(nres)(k:k) ="Q"
          end do
        end if
      end do
!
! Remove redundancies
!
      j = 0
      do i = 1, nres
        do k = 1, j
          if (res_txt(k) == res_txt(i)) exit
        end do
        if (k > j) then
          j = j + 1
          res_txt(j) = res_txt(i)
        end if
      end do
      nres = j
    end if
!
!  Heading
!
    if (len_trim(koment) == 0 .or. (index(koment(:8), " NULL") /= 0)) then
      write(iprt,"(a)") "<HTML><HEAD><TITLE>"//input_fn(:len_trim(input_fn) - 5)//"</TITLE></HEAD>"
    else
      do it = 1, len_trim(koment)
        if (koment(it:it) /= " ") exit
      end do
      write(iprt,"(a)") "<HTML><HEAD><TITLE>"//trim(koment(it:))//"</TITLE></HEAD>"
    end if
    write(iprt,"(a)")"<style type=""text/css"">"
    write(iprt,"(a)")".auto-style4 {"
    write(iprt,"(a)")"margin-top: 0px;"
    write(iprt,"(a)")"text-align: center;"
    write(iprt,"(a)")"line-height: 85%;"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")".auto-style5 {"
    write(iprt,"(a)")"text-decoration: none;"
    write(iprt,"(a)")"}"
    write(iprt,"(a)")"</style>"
    write(iprt,"(a)") "<!--"," ","   Start of JSmol script"," ", "-->"
    write(iprt,"(a)") "<meta charset=""utf-8""> <script type=""text/javascript"" src=""../jsmol/JSmol.min.js""></script>  "
    write(iprt,"(a)") "<script type=""text/javascript"">"
    write(iprt,"(/a)") "$(document).ready(function() {Info = {"
    write(iprt,"(10x,a)") "width: 1500,", "height: 1000,", "color: ""0xB0B0B0"",", &
    "disableInitialConsole: true, ", "addSelectionOptions: false,", "j2sPath: ""../jsmol/j2s"",", &
    "jarPath: ""../jsmol/java"",", "use: ""HTML5"", script: ", " "
    write(iprt,"(a)")"// Data set to be loaded", " "
    if (index(keywrd, " GRAPHF") /= 0) then
      line = input_fn(:len_trim(input_fn) - 4)//"mgf"
    else
      line = input_fn(:len_trim(input_fn) - 4)//"pdb"
    end if
    do i = len_trim(line), 1, -1
      if (line(i:i) == "/" .or. line(i:i) == backslash) exit
    end do
    write(iprt,"(a)") """load "//backslash
    l_geo_ref = (index(keywrd, " COMPARE") /= 0)
    if (l_geo_ref) then
      if (l_compare) then
        line_1 = line(i + 1:len_trim(line)-4)//"_"
      else
        line_1 = " "
      end if
      if (geo_dat_name(:len_trim(geo_dat_name) -3) == geo_ref_name(:len_trim(geo_ref_name) -3)) then
        line = trim(line_1)//geo_dat_name(:len_trim(geo_dat_name) - 4)//"_a.pdb"//"' '"// &
        trim(line_1)//geo_ref_name(:len_trim(geo_ref_name) - 3)//"pdb"//"'; "//backslash
      else
        line = trim(line_1)//geo_dat_name(:len_trim(geo_dat_name) -3)//"pdb"//"' '"// &
        trim(line_1)//geo_ref_name(:len_trim(geo_ref_name) - 3)//"pdb"//"'; "//backslash
      end if
      write(iprt,"(10x,a)")"FILES '"//trim(line)
    else
      write(iprt,"(10x,a)")"'"//trim(line(i + 1:))//"'; "//backslash
    end if
    write(iprt,"(10x,a)")"set measurementUnits ANGSTROMS; "//backslash
    if (l_geo_ref) then
      write(iprt,"(10x,a)")"set bondRadiusMilliAngstroms (25); "//backslash, "spacefill 10%; "//backslash
    else
      write(iprt,"(10x,a)")"set bondRadiusMilliAngstroms (50); "//backslash, "spacefill 15%; "//backslash
    end if
    write(iprt,"(10x,a)")"set display selected; "//backslash, "hBonds calculate; "//backslash
    write(iprt,"(10x,a)")"set defaultDistanceLabel '%0.3VALUE %UNITS'; "//backslash
    if (index(keywrd, " COMPARE") /= 0) then
      write(iprt,"(10x,a)")"select */2.1; color bonds green;  select none;"//backslash
    else
      write(iprt,"(10x,a)")"select none; "//backslash
    end if
    write(iprt,"(10x,a)")"set perspectivedepth off; "//backslash, &
    "connect 0.8  1.5 (hydrogen) (phosphorus) create; "//backslash
    if (l_geo_ref) then
      write(iprt,"(10x,a)")"set zoomLarge false; frame 0;"" "
    else
      write(iprt,"(10x,a)")"set zoomLarge false;"" "
    end if
    write(iprt,"(a)") "} "
    write(iprt,"(a)") "$(""#mydiv"").html(Jmol.getAppletHtml(""jmolApplet0"",Info))}); "
    write(iprt,"(a)") "</script>"
    write(iprt,"(a)") "<!--"," ","   End of JSmol script"," ", "-->"
    if (len_trim(koment) == 0 .or. (index(koment(:8), " NULL") /= 0) ) then
      write(iprt,"(a)") "<h1 align=""center"">"//input_fn(:len_trim(input_fn) - 5)//"</h1>"
    else
      line_1 = trim(koment(it:))
      line = " "
      do i = 1, len_trim(line_1)
        j = len_trim(line) + 1
        if (line_1(i:i) == " ") then
          line(j:) = "&nbsp;"
        else
          line(j:j) = line_1(i:i)
        end if
      end do
    write(iprt,"(a)") "<h1 align=""center"">"//trim(line)//"</h1>"
    end if
    if (len_trim(title) /= 0 .and. (index(title(:8), " NULL") == 0) ) &
      write(iprt,"(a)") '<h2 align="center">'//trim(title)//'</h2>'
    if (index(keywrd, " COMPARE") /= 0) then
      write(iprt,"(a)")"<h2 align=""center"">Compare """//trim(geo_dat_name)// &
      & """ and ""<span style=""color:green"">"//trim(geo_ref_name)//"</span>""</h2>"
    end if
    write(iprt,"(a)") "<BODY  BGCOLOR=""#ffffff"">"
!
!  Start of table.  The table consists of two side-by-side cells
!  The first entry in the first cell is a table of size three rows, two columns
!
    write(iprt,"(a)")"<TABLE>", "<TD>", "<TABLE>", "<TR>"
    if (nres < limres) then
!
!   Element(1,1) and (2,1)
!
      write(iprt,"(a)") "<TD colspan=""2"">"
      call write_data_to_html(iprt)
      write(iprt,"(a)") "</TD></TR><TR>"
    end if
!
!   Element(1,2)
!
    if (index(keywrd, " COMPARE") /= 0) then
      write(iprt,"(a)") "<TD>Toggle display<br><a href=""javascript:Jmol.script(jmolApplet0,'"
    else
      write(iprt,"(a)") "<TD><a href=""javascript:Jmol.script(jmolApplet0,'"
    end if
    write(iprt,"(a)") "if (isOK1);  Display *; zoom 0; isOK1 = FALSE; else hide *; "
    line = " "
    biggest_res = 0
    do i = 1, nres
      n_res = "["//res_txt(i)(:3)//"]"
      j = nint(reada(res_txt(i),4))
      biggest_res = max(biggest_res, j)
      if (j < 0) then
        num = char(Int(log10(-j + 1.0)) + ichar("1"))
        write(n_res(6:),'(a1,i'//num//',a)')"_", -j, res_txt(i)(8:9)
      else
        num = char(Int(log10(j + 1.0)) + ichar("2"))
        write(n_res(6:),'(i'//num//',a)')j, res_txt(i)(8:9)
      end if
      wrt_res = n_res(2:4)//n_res(6:)
      wrt_res(2:2) = char(ichar(wrt_res(2:2)) + ichar("a") - ichar("A"))
      wrt_res(3:3) = char(ichar(wrt_res(3:3)) + ichar("a") - ichar("A"))
      j = len_trim(line)
      if (j > 120) then
        j = 0
        write(iprt,"(a)")trim(line)
        line = " "
      end if
      l_res = res_txt(i)(:7)//res_txt(i)(9:9)
      if (l_res(1:1) >= "0" .and. l_res(1:1) <= "9") &
        l_res(1:1) = char(ichar("A") + ichar(l_res(1:1)) - ichar("0"))
      k = index(l_res, "-")
      if (k > 0) l_res(k:k) = "_"
      write(line(j + 2:),"(a)") l_res//" = FALSE;"
    end do
    write(iprt,"(a)")trim(line)
    write(iprt,"(a)") "isOK1 = TRUE; isOK2 = FALSE; lzoom = TRUE; lcenter = TRUE;"
    if (index(keywrd, " COMPARE") /= 0) then
      write(iprt,"(a)") "endif;')""> 1&<span style=""color:green"">2</span></a>"
      write(iprt,"(a)")"<a href=""javascript:Jmol.script(jmolApplet0,'if (isOKone);  Display add */1.1; zoom 0; isOKone = FALSE;", &
        "else hide add */1.1; zoom 0; isOKone = TRUE;  lzoom = TRUE; lcenter = TRUE; endif;')""> 1</a>"
      write(iprt,"(a)")"<a href=""javascript:Jmol.script(jmolApplet0,'if (isOKtwo);  Display add */2.1; zoom 0; isOKtwo = FALSE;", &
        "else hide add */2.1; zoom 0; isOKtwo = TRUE;  lzoom = TRUE; "// &
      & "lcenter = TRUE; endif;')""><span style=""color:green"">2</span></a>"
    else
      write(iprt,"(a)") "endif;')"">Toggle display all</a>"
    end if
    write(iprt,"(a)")"</TD>"
!
!   Element(2,2)
!
    write(iprt,"(a)")"<TD><a href=""javascript:Jmol.script(jmolApplet0,'"
    write(iprt,"(a)")"if (lcenter);  lcenter = FALSE; else lcenter = TRUE; center {visible}; endif;"
    write(iprt,"(a)")"')"">Toggle center picture</a> </TD>","</TR> <TR>"
!
!   Element(1,3)
!
    write(iprt,"(a)")"<TD><a href=""javascript:Jmol.script(jmolApplet0,'console;')"">Console</a> </TD> "
    write(iprt,"(a)") ""
!
!   Element(2,3)
!
!
!   Element(1,4)
!
    write(iprt,"(a)")"<TD> <a href=""javascript:Jmol.script(jmolApplet0,'if (lzoom); zoom 0; "// &
      "select */2.1; color bonds green; select none; lzoom = FALSE; else lzoom = TRUE; "// &
    "end if ')"">Fit to screen</a> </TD></TR><TR>"
!
!   Element(2,4)
!
    write(iprt,"(a)")"<TD><a href=""javascript:Jmol.script(jmolApplet0,'display within(3,visible);"// &
      "select */2.1; color bonds green; select none; if (lzoom); zoom 0; end if ')"">Near Neighbors</a> </TD> "
     write(iprt,"(a)")"<TD><a href=""javascript:Jmol.script(jmolApplet0,'"// &
     "connect (all) (all) delete; connect; display add connected(visible); select *; hbonds calculate; select "// &
      " */2.1; color bonds green; select none; if (lzoom); delay 0.001; zoom 0; end if ')"">Connected Neighbors</a> </TD> "
    if (allocated(p)) then
      if (p(1) > -900.d0) then
        write(iprt,"(a)")"<TR><TD>"
        write(iprt,"(a)")"<a href=""javascript:Jmol.script(jmolApplet0,'if (!lcharge_x); "
        write(iprt,"(a)")"frame 0; var use = {visible}; frame 1; select off; var sel = use;"
        write(iprt,"(a)")"var z = 0; for (var i IN @sel){z = 3}"
        write(iprt,"(a)")"if (z = 3); use = sel; endif;"
        write(iprt,"(a)")"for (var x IN @use){select @x; var txt =  (x.temperature > 0 ? "//backslash// &
          "'+"//backslash//"':"//backslash//"'"//backslash//"')"// &
          "+format("//backslash//"'%1.3f"//backslash//"',x.temperature*0.1 ); label @txt; color label black;"
        write(iprt,"(a)")"set labelOffset 0 0;}  select @sel; lcharge_x= TRUE;"
        write(iprt,"(a)")"else lcharge_x= FALSE; var use = {visible}; var sel = {selected};"
        write(iprt,"(a)")"var z = 0; for (var i IN @sel){z = 3}"
        write(iprt,"(a)")"if (z = 3); use = sel; endif;"
        write(iprt,"(a)")"select @use; label OFF; select @sel; end if ')"">Charges as Nos.</a>"

        write(iprt,"(a)")"</TD><TD>"
        write(iprt,"(a)")"<a href=""javascript:Jmol.script(jmolApplet0,'if (!lcharge_s); "
        write(iprt,"(a)")"frame 0; var use = {visible}; frame 1; select off; var sel = use;"
        write(iprt,"(a)")"var z = 0; for (var i IN @sel){z = 3}"
        write(iprt,"(a)")"if (z = 3); use = sel; endif;"
        write(iprt,"(a)")"for (var x IN @use)"
        write(iprt,"(a)")"{select @x; var txt =  @x.temperature*0.05;"
        write(iprt,"(a)")"if (@txt > 0){spacefill @txt; color atom deepskyblue;}"
        write(iprt,"(a)")"if (!@txt > 0){txt = -txt; spacefill @txt; color atom deeppink;}}"
        write(iprt,"(a)")"select @sel; lcharge_s= TRUE;"
        write(iprt,"(a)")"else lcharge_s= FALSE; var use = {visible}; var sel = {selected};"
        write(iprt,"(a)")"var z = 0; for (var i IN @sel){z = 3}"
        write(iprt,"(a)")"if (z = 3); use = sel; endif;"
        write(iprt,"(a)")"select @use; spacefill 15%; color cpk; select @sel; end if ')"">Charges as Sizes</a>"
        write(iprt,"(a)")"</TD></TR>"
      end if
    end if
    if (index(keywrd,  " IRC") + index(keywrd,  " STEP=") /= 0) then
!
!  Animation instructions
!
      write(iprt,"(a)")"<TR><TD> <a href=""javascript:Jmol.script(jmolApplet0,'model first;')"">first</a> </TD>"
      write(iprt,"(a)")"<TD> <a href=""javascript:Jmol.script(jmolApplet0,'model last;')"">last</a> </TD></TR>"
      write(iprt,"(a)")"<TR><TD> <a href=""javascript:Jmol.script(jmolApplet0,'model prev;')"">previous</a></TD> "
      write(iprt,"(a)")"<TD> <a href=""javascript:Jmol.script(jmolApplet0,'model next;')"">next</a> </TD></TR>"
      write(iprt,"(a)")"<TR><TD> <a href=""javascript:Jmol.script(jmolApplet0,'animation mode loop 0 0;animation play;')"">"// &
      "loop</a>&nbsp;&nbsp;"
      write(iprt,"(a)")"<a href=""javascript:Jmol.script(jmolApplet0,'animation off;')"">off</a> </TD> "
      write(iprt,"(a)")"<TD> <a href=""javascript:Jmol.script(jmolApplet0,'animation mode palindrome 0 0;animation play;')"">"// &
      "palindrome</a> </TD></TR>"
!
! End of animation instructions
!
    end if
    write(iprt,"(a)") "<TR>"
!
!   Element(1,5)
!
    line = input_fn(:len_trim(input_fn) - 4)//"txt"
    line_1 = trim(line)
    call upcase(line_1, len_trim(line_1))
    i = 0
    do j = 1, len_trim(line) - 4
      if (line(j:j) == "/" .or. line(j:j) == backslash) then
        i = j
      end if
    end do
    if (i /= 0) line = line(i + 1:)
    write(iprt,"(a)")"<TD colspan=""2""><a href=""javascript:Jmol.script(jmolApplet0,'script common.txt;')"">"// &
    "<strong style=""font-size:20px"">Common Script</strong>"
    write(iprt,"(a)")"</a>&nbsp; (Read file from:<br> ""<a href=""common.txt""  target=""_blank"">common.txt</a>"")</TD>"
    write(iprt,"(a)")"</TR> <TR>"
    write(iprt,"(a)")"<TD colspan=""2""><a href=""javascript:Jmol.script(jmolApplet0,'script "//backslash// &
    "'"//trim(line)//backslash//"';')"">"//"<strong style=""font-size:20px"">Specific Script</strong>"
    write(iprt,"(a)")"</a>&nbsp; (Read file from:<br> ""<a href="""//trim(line)// &
    """  target=""_blank"">"//trim(line)//"</a>"")</TD>"
    write(iprt,"(a)")"</TR></TABLE>"
!
!  The second entry in the first cell is a table of size (nres/5) rows, ncol columns.
!  Write all residues in table form
!
    if (nres > 0) then
      write(iprt,"(a)")"<p align=""center"">Toggle Individual Residues</p>"
      write(iprt,"(a)") "<TABLE>"
      write(iprt,"(a)") "<TR>"
!
! Number of rows is a maximum of 36.
! If the number of residues is small, put 8 residues on a line
!
      ncol = max(7, nres/16) + 1
      nprt = 1
      do i = 1, nres
        n_res = "["//res_txt(i)(:3)//"]"
        j = nint(reada(res_txt(i),4))
        if (j < 0) then
          num = char(Int(log10(-j + 1.0)) + ichar("2"))
          write(n_res(6:),'(i'//num//',a)')j, res_txt(i)(8:9)
        else
          num = char(Int(log10(j + 1.0)) + ichar("1"))
          write(n_res(6:),'(i'//num//',a)')j, res_txt(i)(8:9)
        end if
        if (maxtxt == 27 .and. res_txt(i)(10:10) /= " ") then
          j = index(n_res,":")
          n_res = n_res(:j - 1)//"^"//res_txt(i)(10:10)//trim(n_res(j:))
        end if
        wrt_res = n_res(2:4)//n_res(6:)
        l_res = res_txt(i)(:7)//res_txt(i)(9:10)
  !
  !  Residue names must not start with a number.  If they do, convert them into a letter
  !  0 => A, 1 => B, etc.
  !
        if (l_res(1:1) >= "0" .and. l_res(1:1) <= "9") &
          l_res(1:1) = char(ichar("A") + ichar(l_res(1:1)) - ichar("0"))
        do j = 1, 23
          if (n_res(2:4) == tyres(j)) exit
        end do
        if (j > 23) then
          n_res(:5) = " "
          wrt_res(1:1) = "X"
        else
          wrt_res(1:1) = tyr(j)
          if (wrt_res(1:1) == "?") wrt_res(1:1) = "X"
        end if
        j = index(l_res, "-")
        if (j > 0) l_res(j:j) = "_"
        if (n_res(2:4) == "UNK") then
          if (res_txt(i)(1:3) == "UNK") then
            if (l_het_only) then
              n_res = res_txt(i)(4:7)
            else
              n_res = "[UNK]"//res_txt(i)(4:7)
            end if
          else
            n_res = res_txt(i)(4:7)
          end if
        end if
        write(iprt,"(2a)") "<TD> <a href=""javascript:Jmol.script(jmolApplet0,'if (!"//l_res, &
        ");   display ADD "//trim(n_res)//";  "//l_res//" = TRUE; else "
        write(iprt,"(2a)") " hide ADD "//trim(n_res)//";  "//l_res//" = FALSE; end if; ", &
        "if (lcenter); center {visible}; endif; if (lzoom); zoom 0; endif;')"" class=""auto-style5"">  "
        j = index(wrt_res, "^")
        if (j > 0) wrt_res = wrt_res(:j - 1)//wrt_res(j + 1:)
        if (index(wrt_res,":") /= 0) then
          j = index(wrt_res,":") + 1
          write(iprt,"(2a)")"<p class=""auto-style4"">"//wrt_res(1:1)//"<br>"//wrt_res(4:j - 2)//"</p> </a></TD>"
        else
          write(iprt,"(2a)")"<p class=""auto-style4"">"//wrt_res(1:1)//"<br>"//wrt_res(4:)//"</p> </a></TD>"
        end if
        if (nprt == ncol) then
          nprt = 1
          write(iprt,'(a)')"</TR> <TR>"
        else
          nprt = nprt + 1
        end if
      end do
      write(iprt,"(a)") "</TR>"
      write(iprt,"(a)") "</TABLE><br> "
    end if
!
!  End of table of residues, and end of the first cell, now to write the JSmol picture (in the cell (1,2))
!
    write(iprt,"(a)") " </TD><TD>"
    if (nres >= limres) then
!
!   Element(1,1) and (2,1)
!
      call write_data_to_html(iprt)
    end if
    write(iprt,"(a)") "<span id=mydiv></span><a href=""javascript:Jmol.script(jmolApplet0)""></a></TD></TABLE>"
    write(iprt,"(a)") "</BODY>"
    write(iprt,"(a)") "</HTML>"
!
!  Write out a simple script file
!
    close (iprt)
    line = input_fn(:len_trim(input_fn) - 4)
    call add_path(line)
    inquire (file=trim(line)//"txt", exist = exists)
    if (.not. exists) then
      open(unit=iprt, file=trim(line)//"txt")
      i = 0
      do j = 1, len_trim(line)
        if (line(j:j) == "/" .or. line(j:j) == backslash) i = j
      end do
      if (i /= 0) line = line(i + 1:)
      write(iprt,"(a)")"#","# Script for use with the HTML file """//trim(line)//"html""","#"
    end if
    close (iprt)
  end subroutine write_html
  subroutine write_data_to_html(iprt)
    use molkst_C, only : numat, formula, escf, nelecs, keywrd, arc_hof_1, arc_hof_2, &
      density, id, mol_weight, stress
    USE parameters_C, only : tore
    use common_arrays_C, only: nat, tvec
    use chanel_C, only: log
    use funcon_C, only : fpc_10
    implicit none
    integer :: iprt
    integer :: i, j, k, l, eqls, colon
    character :: idate*24, line*120
    logical :: store_log
    double precision :: sum
    integer, external :: quoted
    double precision, external :: volume
!
!  Print out information on the system: formula, number of atoms, heat of formation, date, etc.
!
    write(iprt,"(a)")"<TABLE>"
    call fdate (idate)
    write(iprt,"(a)")  "<TR><TD> Date:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>"// &
      idate(5:11)//idate(21:)//idate(11:16)//"</TD></TR>"
    write(iprt,"(a,i5,a)")  "<TR><TD> No. atoms:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>", &
      numat,  "</TD></TR>"
    store_log = log
    log = .false.
    call empiri()
    log = store_log
    colon = index(formula, ":") + 1
    eqls = index(formula, "=") - 1
    line = " "
    do i = colon, eqls
      j = i + 1
      if ((formula(i:i) < "0" .or. formula(i:i) > "9") .and. formula(j:j) >= "0" .and. formula(j:j) <= "9") then
        line = trim(line)//formula(i:i)
        line = trim(line)//"<sub>"
      else if ((formula(j:j) < "0" .or. formula(j:j) > "9") .and. formula(i:i) >= "0" .and. formula(i:i) <= "9") then
        line = trim(line)//formula(i:i)
        line = trim(line)//"</sub>"
      else
        line = trim(line)//formula(i:i)
      end if
    end do
    write(iprt,"(a)")  "<TR><TD> Formula:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>"//trim(line)//"</TD></TR>"
    if (index(keywrd, " 0SCF") == 0 .or. quoted('GEO_REF="') == 0) then
       if (id == 3) then
        sum = volume(tvec,3)
        density = mol_weight*1.D24/fpc_10/sum
        if (density > 1.d-1) write(iprt,"(a,f7.3,a)")  &
        "<TR><TD> Density:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>",density,"</TD></TR>"
      end if
    else
      if (abs(arc_hof_1) > 1.d-4) then
        if ( quoted('GEO_DAT="') /= 0) then
          write(iprt,"(a,f12.3,a)")  "<TR><TD> GEO_DAT:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>", &
          arc_hof_1," kcal/mol</TD></TR>"
        else
          write(iprt,"(a,f12.3,a)")  "<TR><TD> Dataset:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>", &
          arc_hof_1," kcal/mol</TD></TR>"
        end if
      end if
      if (abs(arc_hof_2) > 1.d-4) &
      write(iprt,"(a,f12.3,a)")  "<TR><TD> <span style=""color:green""> GEO_REF:</span></TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>", &
        arc_hof_2," kcal/mol</TD></TR>"
      i = index(keywrd, " RMS_DIFF")
      if (i > 0) then
        i = index(keywrd(i:), "F=") + i + 1
        j = index(keywrd(i:), " ") + i - 2
        k = index(keywrd, " DIFF=") + 6
        l = k
        do l = l, len_trim(keywrd)
          if (keywrd(l:l) >= "0" .and. keywrd(l:l) <= "9") exit
        end do
        l = index(keywrd(l:)," ") + l - 2
        write(iprt,"(a)")  "<TR><TD> RMSD: "//keywrd(i:j)//"&Aring;</TD><TD> </TD><TD>Diff: "//keywrd(k:l)//"&Aring;</TD></TR>"
      end if
    end if
    if (index(line, "H") /= 0 .and. nelecs > 0) then
!
!  Net charge has to be worked out the hard way
!
      sum = -nelecs
      do i = 1, numat
        sum = sum + tore(nat(i))
      end do
      i = nint(sum)
      if (i /= 0) then
        write(iprt,"(a,SP,i6,SS,a)")"<TR><TD> Net charge:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>",i,"</TD></TR>"
      else
          write(iprt,"(a)")"<TR><TD> Net charge:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>Zero</TD></TR>"
      end if
    end if
    if (abs(escf) > 1.d-10) then
      if (abs(escf) > 1.d-10) &
        write(iprt,"(a, f12.3)")"<TR><TD> Heat of Formation:</TD><TD> &nbsp;&nbsp; &nbsp;</TD><TD>", &
        escf - stress," Kcal/mol</TD></TR>"
    end if
    write(iprt,"(a)")  "</TABLE>"
  end subroutine write_data_to_html
