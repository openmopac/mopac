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

      subroutine datin(iw)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : partyp, n_partyp, n_partyp_alpb, v_par, t_par
      use Common_arrays_C, only : ijpars, parsij
      use molkst_C, only : keywrd, lpars, line, backslash
      use chanel_C, only : iext
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
      integer, intent(in) :: iw
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer :: i, j, l, iparam, k, ielmnt, jelmnt, nref, &
      loop, mpar

      double precision :: param
      logical :: exists, lv_par(60) = .false.
      character (len=500), dimension (10) :: file
      character :: text*80,  ch*1
      character, dimension(107) :: elemnt*2
      character, external :: get_a_name*300
      integer, external :: end_of_keyword
      double precision, external :: reada

      save elemnt
!-----------------------------------------------
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/
    if (.not. allocated(ijpars))  allocate(ijpars(5,5000), parsij(5000))
    i = Index(keywrd, "EXTERNAL=") + Index(keywrd, "PARAMS=")
    t_par     = "Add a description of this parameter near line 50 in datin.F90  "
    t_par(1)  = "Used in ccrep for scalar correction of C-C triple bonds."
    t_par(2)  = "Used in ccrep for exponent correction of C-C triple bonds."
    t_par(3)  = "Used in ccrep for scalar correction of O-H term."
    t_par(4)  = "Used in ccrep for exponent correction of O-H term."
    t_par(5)  = "Used in ccrep for offset correction of O-H term."
    t_par(7)  = "Used in dftd3 to set ""s6""  in D3H4"
    t_par(8)  = "Used in dftd3 to set ""alp"" in D3H4"
    t_par(9)  = "Used in dftd3 to set ""rs6"" in D3H4"
    t_par(10) = "Used in dftd3 to set ""s18"" in D3H4."
    t_par(11) = "Used in ccrep for exponent correction of C-H term."
    t_par(12) = "Used in ccrep for offset correction of C-H term."
    t_par(13) = "Used in ccrep for scalar correction of C-C term."
    t_par(14) = "Used in ccrep for exponent correction of C-C term."
    t_par(15) = "Used in ccrep for offset correction of C-C term."
    t_par(16) = "Used in ccrep for scalar correction of H-H term."
    t_par(17) = "Used in ccrep for exponent correction of H-H term."
    t_par(18) = "Used in ccrep for offset correction of H-H term."
    t_par(19) = "Used in ccrep for scalar correction of C-H term."
    t_par(20) = "Used in ccrep for scalar correction of C-O term."
    t_par(21) = "Used in ccrep for exponent correction of C-O term."
    t_par(22) = "Used in ccrep for offset correction of C-O term."
    t_par(23) = "Used in ccrep for scalar correction of S-O term."
    t_par(24) = "Used in ccrep for exponent correction of S-O term."
    t_par(25) = "Used in ccrep for offset correction of S-O term."
    t_par(26) = "Used in ccrep for scalar correction of O-N term."
    t_par(27) = "Used in ccrep for exponent correction of O-N term."
    t_par(28) = "Used in ccrep for offset correction of O-N term."
    t_par(29) = "Used in ccrep for scalar correction of F-H term."
    t_par(30) = "Used in ccrep for exponent correction of F-H term."
    t_par(31) = "Used in ccrep for offset correction of F-H term."
    nref = 0
    k = Index (keywrd(i:), "=") + i
    j = end_of_keyword(keywrd, len_trim(keywrd), k)
!
! k = start of reference data directory list
! j = end of list.
! in between are the names of the reference directories, separated by ";"
!
    do l = 1, 20
      i = Index(keywrd(k:j),";")
      if (i /= 0) then
        nref = nref + 1
        file(nref) = trim(get_a_name(keywrd(k:j), len_trim(keywrd(k:j))))
        k = k + i
      end if
      if (i == 0) then
!
!  Last entry
!
        nref = nref + 1
        file(nref) = keywrd(k:j)
        exit
      end if
    end do
    if (file(nref)(1:1) == '"') file(nref) = file(nref)(2:)
    i = len_trim(file(nref))
    if (file(nref)(i:i) == '"') file(nref) = file(nref)(:i - 1)

  !
  !   Read in parameters from a previous run - these will overwrite
  !   the default values of the parameters.
  !
    lpars = 0
    do 10 loop = 1, nref
      mpar = 0
      call add_path(file(loop))
      inquire (file=trim(file(loop)), exist = exists)
      if (.not. exists) then
        if (index(keywrd,' 0SCF') == 0) call mopend("EXTERNAL file: '"//trim(file(loop))//"' does not exist!")
        line = trim(file(loop))
        do i = len_trim(line), 1, -1
          if (line(i:i) == backslash .or. line(i:i) == "/") then
!  using inquire to check for the existence of a directory is a non-standard Intel extension to Fortran 90
!  there is no standard, system-independent way to do this as directories do not exist in Fortran language specifications
!  I'm just commenting this out for now, since it is just some extra error messages [JEM 2019.05.24]
!            inquire (directory = line(:i), exist = exists)
!            if (exists) then
!              write(iw,"(10x,a)")"  (but folder: '"//line(:i)//"' does exist.)"
!            else
!              if (verson(7:7) == "W") write(iw,"(5x,a)")" (Note: the folder: '"//line(:i)//"'  also does not exist.)"
!            end if
            exit
          end if
        end do
        return
      end if
      open (unit=iext, form="FORMATTED", status="OLD", file=trim(file(loop)), action="READ", iostat = i)
      if (i /= 0) then
        if (lpars > 0) exit
        if (loop == 1) then
          write(line,'(a)')" EXTERNAL file """//trim(file(loop))//""" could not be opened"
          write(iw,'(/,a)')trim(line)
          if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 ) then
            call mopend(trim(line))
            inquire (file=trim(file(loop)), exist = exists)
            if (exists) then
              write(line,'(a)')" (The EXTERNAL file exists, but could not be read)"
              write(iw,'(a)')trim(line)
              call mopend(trim(line))
            else
              write(line,'(a)')" (The EXTERNAL file does not exist)"
              write(iw,'(a)')trim(line)
              call mopend(trim(line))
            end if
          end if
        end if
      end if
      rewind (iext,err = 10)
      do
        read (iext, "(A60)", err=11, end=11) text
        call upcase (text, 80)
        if (Index (text, "END") /= 0 .or. text == " ")  exit
        if (text(1:1) == "*") cycle
        if (index(text, "PAR") /= 0) then
          i = index(text, "PAR") + 3
          j = ichar(text(i:i)) - ichar("0")
          i = i + 1
          line = text (1:80)
          if (text(i:i) /= " ")then
            j = j*10 + ichar(text(i:i)) - ichar("0")
            i = i + 1
          end if
          v_par(j) = reada(text, i)
          lv_par(j) = .true.
          lpars = lpars + 1
          ijpars(1, lpars) = j
          ijpars(2, lpars) = 41
          parsij(lpars) =  reada(text, i)
          mpar = 1
          cycle
        end if
        line = text (1:80)
        do
          if (line(1:1) /= " ") exit
          line = line(2:)
        end do
        j = len_trim(line)
        do i = j, 1, -1
          ch = line(i:i)
          if (ichar(ch) > 126) exit
          if (ichar(ch) < 32) exit
        end do
        if (i /= 0) then
          call mopend(" Non-standard characters detected in the EXTERNAL file")
          l = 0
          write(iw,'(/, "Faulty line")')
          do
            l = l + 20
            k = min(l,len_trim(line))
            if (k < l - 19) exit
            write(iw,'(a,i3,30i4)')"Position:", (j, j = l - 19, k)
            write(iw,'(a,i3,30i4)')"   ASCII:", (ichar(line(j:j)), j = l - 19, k)
            write(iw,'(a,a3,29a4)')"    Char:", (line(j:j), j = l - 19, k)
            write(iw,*)
            if (k == i) exit
          end do
          write(iw,'(//,a,//)')" (Edit the EXTERNAL to remove non-standard characters using a primitive editor such as vi)"
          return
        end if
!
! Clean up line - delete anything after the third set of spaces
!
        do
          j = 0
          do i = 1, len_trim(line)
            if (line(i:i + 1) == "  ") then
              j = 1
              line(i:) = line(i + 1:)
            end if
          end do
          if (j == 0) exit
        end do
        i = Index (line, " ")
        i = index(line(i + 1:)," ") + i + 1
        i = index(line(i + 1:)," ") + i
        text = line(1:i)
!
!  Force in spaces needed for parsing
!
        i = index(text," ")
        text(i:) = "  "//text(i:)
        i = index(text(i + 3:)," ") + i + 3
        text(i:) = "  "//text(i:)
        line = text
        do j = 1, n_partyp
          if (Index(" "//text, " "//partyp(j)) /= 0) go to 1000
        end do
        write(iw,"(3a)")" EXTERNAL parameter type: '",trim(text),"' unrecognized"
        close(iext)
        goto 99
1000    iparam = j
        jelmnt=0
        if (iparam == n_partyp_alpb .or. iparam == n_partyp_alpb + 1) then
!
!  This is a di-atomic parameter - read in the other element number
!
          i = Index(text, partyp(j))+5
          do j = 1, 99
            if (Index (" "//text(i:i+2), " "//elemnt(j)) /= 0) exit
          end do
          jelmnt = j
          end if
        i = Index (line, " ")
        text = line(i:)
        line = text
        do i = 1, 88
          if (line(1:1) /= " ") exit
          text = line(2:60)
          line = text
        end do
        text = " " // line (1:79)
        do j = 1, 100
          if (Index (text, " "//elemnt(j)) /= 0) go to 1100
        end do
        write(iw,"(3a)")" EXTERNAL element type: '",trim(text),"' unrecognized"
        close(iext)
        goto 99
1100    param = reada (text, Index (text, " "//elemnt(j)))
        if (j > jelmnt) then
          ielmnt = j+200*jelmnt
        else
          ielmnt = jelmnt+200*j
        end if
        do i = 1, lpars
          if (ijpars(1, i) == ielmnt .and. ijpars(2, i) == iparam) go to 1200
        end do
        lpars = lpars + 1
        i = lpars
1200    ijpars(1, i) = ielmnt
        ijpars(2, i) = iparam
        if (Abs(param) < 1.d-7) param = 1.d-7  ! Don't allow new parameters to be exactly zero
        parsij(i) = param
        text = " "
        mpar = 1
      end do
11  continue
    if(mpar == 1) then
      write(iw,'(/,3a)')" Parameters read in from file: """, file(loop)(:len_trim(file(loop))),""""
    else
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 ) then
        call mopend("No parameters read in from '"//file(loop)(:len_trim(file(loop)))//"'")
      end if
    end if
10  continue
    close (iext, status="KEEP")
    call write_params(iw, lv_par)
    do i = 1, lpars
      call update(ijpars(2, i), ijpars(1, i), parsij(i), 0.d0)
    end do
    close(iext)
 99   return
    end subroutine datin
!
!
!
    ! copy of datin for processing INTERNAL parameters
      subroutine datex(ir, iw)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameters_C, only : partyp, n_partyp, n_partyp_alpb, v_par, t_par
      use Common_arrays_C, only : ijpars, parsij
      use molkst_C, only : keywrd, lpars, line, backslash
      use chanel_C, only : iext
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
      integer, intent(in) :: iw, ir
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      integer :: i, j, l, iparam, k, ielmnt, jelmnt, nref, &
      loop, mpar

      double precision :: param
      logical :: exists, lv_par(60) = .false.
      character (len=500), dimension (10) :: file
      character :: text*80,  ch*1
      character, dimension(107) :: elemnt*2
      character, external :: get_a_name*300
      integer, external :: end_of_keyword
      double precision, external :: reada

      character :: infile * 100
      logical :: first = .true.

      save elemnt
!-----------------------------------------------
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/
    if (.not. allocated(ijpars))  allocate(ijpars(5,5000), parsij(5000))
    i = Index(keywrd, "INTERNAL=") + Index(keywrd, "PARAMS=")
    t_par     = "Add a description of this parameter near line 50 in datin.F90  "
    t_par(1)  = "Used in ccrep for scalar correction of C-C triple bonds."
    t_par(2)  = "Used in ccrep for exponent correction of C-C triple bonds."
    t_par(3)  = "Used in ccrep for scalar correction of O-H term."
    t_par(4)  = "Used in ccrep for exponent correction of O-H term."
    t_par(5)  = "Used in ccrep for offset correction of O-H term."
    t_par(7)  = "Used in dftd3 to set ""s6""  in D3H4"
    t_par(8)  = "Used in dftd3 to set ""alp"" in D3H4"
    t_par(9)  = "Used in dftd3 to set ""rs6"" in D3H4"
    t_par(10) = "Used in dftd3 to set ""s18"" in D3H4."
    t_par(11) = "Used in ccrep for exponent correction of C-H term."
    t_par(12) = "Used in ccrep for offset correction of C-H term."
    t_par(13) = "Used in ccrep for scalar correction of C-C term."
    t_par(14) = "Used in ccrep for exponent correction of C-C term."
    t_par(15) = "Used in ccrep for offset correction of C-C term."
    t_par(16) = "Used in ccrep for scalar correction of H-H term."
    t_par(17) = "Used in ccrep for exponent correction of H-H term."
    t_par(18) = "Used in ccrep for offset correction of H-H term."
    t_par(19) = "Used in ccrep for scalar correction of C-H term."
    t_par(20) = "Used in ccrep for scalar correction of C-O term."
    t_par(21) = "Used in ccrep for exponent correction of C-O term."
    t_par(22) = "Used in ccrep for offset correction of C-O term."
    t_par(23) = "Used in ccrep for scalar correction of S-O term."
    t_par(24) = "Used in ccrep for exponent correction of S-O term."
    t_par(25) = "Used in ccrep for offset correction of S-O term."
    t_par(26) = "Used in ccrep for scalar correction of O-N term."
    t_par(27) = "Used in ccrep for exponent correction of O-N term."
    t_par(28) = "Used in ccrep for offset correction of O-N term."
    t_par(29) = "Used in ccrep for scalar correction of F-H term."
    t_par(30) = "Used in ccrep for exponent correction of F-H term."
    t_par(31) = "Used in ccrep for offset correction of F-H term."
    nref = 0
    k = Index (keywrd(i:), "=") + i
    j = end_of_keyword(keywrd, len_trim(keywrd), k)

    ! internal parameter specification must take place after a line containing
    ! INTERNAL and continues until the end of the file
    rewind(ir)
    do 2
       read(ir, '(A)', end=3) infile
       if (first) then
          first = .false.
          goto 2
       end if
       call upcase(infile, 80)
       i = Index(infile, "INTERNAL")
       if (i /= 0) then
          ! found the parameters, move on
          goto 3
       end if
2      continue
3      continue
!
! k = start of reference data directory list
! j = end of list.
! in between are the names of the reference directories, separated by ";"
!
    do l = 1, 20
      i = Index(keywrd(k:j),";")
      if (i /= 0) then
        nref = nref + 1
        file(nref) = trim(get_a_name(keywrd(k:j), len_trim(keywrd(k:j))))
        k = k + i
      end if
      if (i == 0) then
!
!  Last entry
!
        nref = nref + 1
        file(nref) = keywrd(k:j)
        exit
      end if
    end do
    if (file(nref)(1:1) == '"') file(nref) = file(nref)(2:)
    i = len_trim(file(nref))
    if (file(nref)(i:i) == '"') file(nref) = file(nref)(:i - 1)

  !
  !   Read in parameters from a previous run - these will overwrite
  !   the default values of the parameters.
  !
    lpars = 0
    do 10 loop = 1, nref
      mpar = 0
      do
        read (ir, "(A60)", err=11, end=11) text
        call upcase (text, 80)
        if (Index (text, "END") /= 0 .or. text == " ")  exit
        if (text(1:1) == "*") cycle
        if (index(text, "PAR") /= 0) then
          i = index(text, "PAR") + 3
          j = ichar(text(i:i)) - ichar("0")
          i = i + 1
          line = text (1:80)
          if (text(i:i) /= " ")then
            j = j*10 + ichar(text(i:i)) - ichar("0")
            i = i + 1
          end if
          v_par(j) = reada(text, i)
          lv_par(j) = .true.
          lpars = lpars + 1
          ijpars(1, lpars) = j
          ijpars(2, lpars) = 41
          parsij(lpars) =  reada(text, i)
          mpar = 1
          cycle
        end if
        line = text (1:80)
        do
          if (line(1:1) /= " ") exit
          line = line(2:)
        end do
        j = len_trim(line)
        do i = j, 1, -1
          ch = line(i:i)
          if (ichar(ch) > 126) exit
          if (ichar(ch) < 32) exit
        end do
        if (i /= 0) then
          call mopend(" Non-standard characters detected in the INTERNAL file")
          l = 0
          write(iw,'(/, "Faulty line")')
          do
            l = l + 20
            k = min(l,len_trim(line))
            if (k < l - 19) exit
            write(iw,'(a,i3,30i4)')"Position:", (j, j = l - 19, k)
            write(iw,'(a,i3,30i4)')"   ASCII:", (ichar(line(j:j)), j = l - 19, k)
            write(iw,'(a,a3,29a4)')"    Char:", (line(j:j), j = l - 19, k)
            write(iw,*)
            if (k == i) exit
          end do
          write(iw,'(//,a,//)')" (Edit the INTERNAL to remove non-standard characters using a primitive editor such as vi)"
          return
        end if
!
! Clean up line - delete anything after the third set of spaces
!
        do
          j = 0
          do i = 1, len_trim(line)
            if (line(i:i + 1) == "  ") then
              j = 1
              line(i:) = line(i + 1:)
            end if
          end do
          if (j == 0) exit
        end do
        i = Index (line, " ")
        i = index(line(i + 1:)," ") + i + 1
        i = index(line(i + 1:)," ") + i
        text = line(1:i)
!
!  Force in spaces needed for parsing
!
        i = index(text," ")
        text(i:) = "  "//text(i:)
        i = index(text(i + 3:)," ") + i + 3
        text(i:) = "  "//text(i:)
        line = text
        do j = 1, n_partyp
          if (Index(" "//text, " "//partyp(j)) /= 0) go to 1000
        end do
        write(iw,"(3a)")" INTERNAL parameter type: '",trim(text),"' unrecognized"
        goto 99
1000    iparam = j
        jelmnt=0
        if (iparam == n_partyp_alpb .or. iparam == n_partyp_alpb + 1) then
!
!  This is a di-atomic parameter - read in the other element number
!
          i = Index(text, partyp(j))+5
          do j = 1, 99
            if (Index (" "//text(i:i+2), " "//elemnt(j)) /= 0) exit
          end do
          jelmnt = j
          end if
        i = Index (line, " ")
        text = line(i:)
        line = text
        do i = 1, 88
          if (line(1:1) /= " ") exit
          text = line(2:60)
          line = text
        end do
        text = " " // line (1:79)
        do j = 1, 100
          if (Index (text, " "//elemnt(j)) /= 0) go to 1100
        end do
        write(iw,"(3a)")" INTERNAL element type: '",trim(text),"' unrecognized"
        goto 99
1100    param = reada (text, Index (text, " "//elemnt(j)))
        if (j > jelmnt) then
          ielmnt = j+200*jelmnt
        else
          ielmnt = jelmnt+200*j
        end if
        do i = 1, lpars
          if (ijpars(1, i) == ielmnt .and. ijpars(2, i) == iparam) go to 1200
        end do
        lpars = lpars + 1
        i = lpars
1200    ijpars(1, i) = ielmnt
        ijpars(2, i) = iparam
        if (Abs(param) < 1.d-7) param = 1.d-7  ! Don't allow new parameters to be exactly zero
        parsij(i) = param
        text = " "
        mpar = 1
      end do
11  continue
    if(mpar == 1) then
      write(iw,'(/,3a)')" Parameters read in from main input file"
    else
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 ) then
        call mopend("No parameters read in from '"//file(loop)(:len_trim(file(loop)))//"'")
      end if
    end if
10  continue
    call write_params(iw, lv_par)
    do i = 1, lpars
      call update(ijpars(2, i), ijpars(1, i), parsij(i), 0.d0)
    end do
 99   return
    end subroutine datex
    !
    !
    !
    subroutine write_params(iw, lv_par)
      USE parameters_C, only : partyp, n_partyp_alpb, n_partyp, v_par, t_par
      use Common_arrays_C, only : ijpars, parsij
      use molkst_C, only : lpars, line
      implicit none
      integer, intent(in) :: iw
      logical, intent (in) :: lv_par(60)
      logical :: lold = .true.
!
!  Local
!
      integer :: i, j, k, l, il, iu, ii, jj, iparam, ielmnt, jelmnt
      character :: elemnt(107)*2, elemnt2*2
      save elemnt, lold
!-----------------------------------------------
      data (elemnt(i),i=1,107)/ 'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O '&
        , 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', &
        'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA'&
        , 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', &
        'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE'&
        , 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', &
        'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR'&
        , 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', &
        'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', 'FM'&
        , 'MD', 'CB', '++', '+', '--', '-', 'TV'/
      do j = 1, 107
        do k = 1, n_partyp
          if (k == n_partyp_alpb + 1) cycle
          if (k == n_partyp_alpb) then
           il = 1
           iu = 98
          else
           il = 0
           iu = 0
          end if
          do ii = il, iu
            do jj = 0, il
              do i = 1, lpars
                iparam = ijpars(2, i)
                ielmnt = ijpars(1, i)
                jelmnt = mod(ielmnt,200)
                l = ielmnt/200
                if (iparam == k + jj .and. jelmnt == j .and. l == ii) then
                  if (lold) then
                    write (iw, "(//,8X,A)") " Parameters read in"
                    write (iw, "(/5X, ' Parameter Type  Element    Parameter')")
                    lold = .false.
                  end if
                  if(l /= 0) then
                    elemnt2 = elemnt(l)
                    if(elemnt2(2:2) == " ")elemnt2 = elemnt2(1:1)
                  else
                    elemnt2 = " "
                  end if
                  write (line, "(12X,A7,7X,A2,F16.8)") partyp (iparam)//elemnt2, elemnt(jelmnt), parsij(i)
                  if (iparam /= 41 .and. iparam < n_partyp) then
                    write (iw, "(a)") trim(line)
                  else
                   continue
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
!
! Write out the global parameters that have been read in
!
      do i = 1, 9
        if (lv_par(i)) write(iw,"(12x,'PAR', i1, f28.8, a)") i, v_par(i), "    "//trim(t_par(i))
      end do
       do i = 10, 60
        if (lv_par(i)) write(iw,"(12x,'PAR', i2, f27.8, a)") i, v_par(i), "    "// trim(t_par(i))
      end do
      return
  end subroutine write_params
    integer function end_of_keyword(input, len, start)
!
! Given a keyword, input(start:), find the end of the keyword
!
  integer, intent (in) :: len, start
  character, intent (in) :: input*(len)
!
! Local
!
  integer :: j, k
  logical :: quotation_mark
    quotation_mark = .false.
    k = start
    do
      if (input(k:k) /= " ") exit
      k = k + 1
    end do
    j = k - 1
    do
      j = j + 1
      if (j > len) exit
      if (input(j:j) == ' ') exit
      if (input(j:j) == '"')  quotation_mark = (.not.  quotation_mark)
      if (quotation_mark) then
        do
          j = j + 1
          if (j > len) exit
          if (input(j:j) == '"')  quotation_mark = (.not.  quotation_mark)
          if (.not. quotation_mark) exit
        end do
      end if
    end do
    end_of_keyword = j
    return
  end function end_of_keyword
  character(LEN=300) function get_a_name(input, len)
!
! Given a string of the type "aaa;bbb;ccc" separate
! off the first string, here "aaa"
!
! Note: In input, the first character MUST be the first character of the name or
!       a quotation mark
!
  integer, intent (in) :: len
  character, intent (in) :: input*(len)
!
! Local variables
!
  integer :: istart, istop
  logical :: quotation_marks
  istart = 1
  do 
    if (input(istart:istart) /= " ") exit
    istart = istart + 1
  end do
  quotation_marks = (input(istart:istart) == '"')
  if (quotation_marks) then
!
!    name starts and ends with quotation marks
!
    istart = istart + 1
    istop = index(input(istart:), '"') + istart - 2
  else
    do istop = istart + 1, len
      if (input(istop:istop)== " " .or. input(istop:istop) == ";") exit
    end do
    istop = istop - 1
  end if
  get_a_name = input(istart:istop)
  return
  end function get_a_name
