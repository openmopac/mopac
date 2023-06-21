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

      subroutine gettxt
      use chanel_C, only: ir, iw, isetup, input_fn
      use molkst_C, only: keywrd, keywrd_quoted, koment, title, refkey, gui, numcal, line, &
        moperr, allkey, backslash
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, k, l, m, ipath
      character :: filen*300, oldkey*3000, line1*3000, path*240, ch*1
      logical :: aux, exists, setup_present, zero_scf
      character (len = 300), external :: get_text
      integer, external :: quoted
!-----------------------------------------------
      koment = " "
      title = " "
      refkey = '    NULL  '
      keywrd_quoted = " "
      oldkey = " "
      aux = (index(keywrd, "AUX") /= 0)
      read (ir, '(A2000)', end=100, err=100) refkey(1)
      keywrd = refkey(1)
      do
        j = index(keywrd, "++")
        if (j /= 0) then
          do
            j = index(keywrd, "++")
            if (j== 0) exit
            i = j
            if (keywrd(i - 1:i - 1) == " ")then
              line = keywrd(:i - 1)//" "//trim(keywrd(i + 2:))
            else
              line = keywrd(:i - 1)//" "//trim(keywrd(i + 2:))
            end if
            keywrd = trim(line)
          end do
          read (ir, '(A2000)', end=100, err=100) line
          if (keywrd(i - 1:i - 1) == " ")then
              keywrd = trim(keywrd)//" "//trim(line)
            else
              keywrd = trim(keywrd)//trim(line)
            end if
        else
          refkey(1) = trim(keywrd)
          exit
        end if
      end do
      oldkey = trim(keywrd)
      call upcase (keywrd, len_trim(keywrd))
      zero_scf = (index(keywrd, "0SCF") /= 0)
      do i = len_trim(input_fn), 2, -1
        if (input_fn(i:i) == backslash .or. input_fn(i:i) == "/") exit
      end do
      ipath = i
      if  (ipath > 2) then
        path = input_fn(:ipath)
        filen = trim(path)//'SETUP'
      else
        filen = 'SETUP'
      end if
      inquire (file=filen, exist = exists)
      i = len_trim(keywrd)
      allkey = keywrd(:i)
      setup_present = (index(keywrd,'SETUP') /= 0)
      if (setup_present) then
        if (index(allkey, " + ") /= 0) then
          call mopend(" Keywords ""SETUP"" and ""+""cannot be used together")
!
!  Dummy read to go to the end of the data set
!
          do
            read(ir,*,err=98, iostat = i)filen
            if (i /= 0) exit
          end do
  98      return
        end if
        i = index(keywrd,'SETUP=')
        if (i /= 0) then
          filen = trim(get_text(oldkey, i + 6, 1))
        else
          i = index(keywrd,'SETUP')
          j = index(keywrd(i:), " ") + i - 1
          filen = keywrd(i:j)
        end if
        call add_path(filen)
        inquire (file=filen, exist = exists)
        if (.not. exists) then
          inquire (file=trim(filen)//".txt", exist = exists)
          if (exists) filen = trim(filen)//".txt"
        end if
        if (.not. exists) then
          if (setup_present .and. .not. zero_scf) then
            write (line, '(A)') "SETUP FILE """//trim(filen)//""" MISSING."
            numcal = 2
            if (.not. gui )write(0,'(//30x,a)')' SETUP FILE "'//trim(filen)//'" MISSING'
            call mopend (trim(line))
            return
          end if
        else
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED',position='REWIND', iostat=i)
          if (i /= 0) then
            if (.not. zero_scf) then
              call mopend ('COULD NOT OPEN SETUP FILE: '//trim(filen))
              if (zero_scf) moperr = .false.
              return
            end if
          end if
          rewind isetup
          refkey(2) = " "
          i = -1
          do
            read (isetup, '(A)', end=61, err=50) line1
61          if (line1(1:1) == "*") cycle
            j = len_trim(line1)
            if (j == 3000) line1 = " "
            if (i == -1) then
              refkey(2) = trim(refkey(2))//trim(line1)
            else
              if (refkey(2)(i - 1:i - 1) == " ")then
                refkey(2) = trim(refkey(2))//" "//trim(line1)
              else
                refkey(2) = trim(refkey(2))//trim(line1)
              end if
            end if
            j = index(refkey(2), "++")
            if (j /= 0) then
              do
                j = index(refkey(2), "++")
                if (j == 0) exit
                i = j
                if (refkey(2)(i - 1:i - 1) == " ")then
                  line = refkey(2)(:i - 1)//" "//trim(refkey(2)(i + 2:))
                else
                  line = refkey(2)(:i - 1)//" "//trim(refkey(2)(i + 2:))
                end if
                refkey(2) = trim(line)
                line = refkey(2)(:i - 1)//trim(refkey(2)(i + 2:))
              end do
            else
              oldkey = trim(oldkey)//" "//trim(refkey(2))
              exit
            end if
          end do
          close (isetup)
          call upcase (refkey(2), len(refkey(2)) )
!
!  Check for " -" signs in setup file
!
          i = ichar(refkey(2)(1:1))
          if (i == 0) then
            refkey(2) = " "
          else if (refkey(2)(1:1) /= " ") then
            refkey(2) = " "//refkey(2)(:len_trim(refkey(2)))
          end if
          do
            i = index(refkey(2), " -")
            if (i == 0) exit
!
! Is the minus sign inside a quoted text
!
            j = 0
            do k = i, len_trim(refkey(2))
              if (refkey(2)(k:k) == '"') j = j + 1
            end do
            if (mod(j,2) == 1) then
!
!  Yes, it's inside a quoted string, so protect it and carry on
!
              refkey(2)(i:i) = "*"
              cycle
            end if
            j = index(refkey(2)(i + 2:), " ") + i + 1
            do
              k = index(" "//keywrd, " "//refkey(2)(i + 2: j - 1))
              if (k == 0) exit
              l = index(keywrd(k + 1:)," ") + k + 1
              keywrd(k:) = keywrd(l:)
              if (keywrd == " ") exit
            end do
            refkey(2)(i:) = refkey(2)(j:)
          end do
!
!   Check for keywords in SETUP that are present in the keyword line, and delete them
!   ( Keywords on keyword line take precedence.)
!
!  First, move keywords that contain quoted text into keywrd_quoted so that they
!  will not affect keyword analyses
!
          i = 1
          do
            i = i + 1
            if (i > len_trim(refkey(2))) exit
            if (refkey(2)(i - 1:i - 1) == " " .and. refkey(2)(i:i) /= " ") then
!
!  Found a keyword. Now look for the end of the keyword
!
              do j = i + 1, len_trim(refkey(2))
                if (refkey(2)(j:j) == " ") exit
              end do
              line = refkey(2)(i:j)
              do k = 1,100
                if ((line(k:k) < "A" .or. line(k:k) > "Z") .and. &
                    (line(k:k) < "0" .or. line(k:k) > "9") .and. &
                    (line(k:k) /= "_")) exit
              end do
              k = k - 1
              if (index(allkey, " "//line(:k)) > 0 .and. line(:5) /= "SETUP") then
                k = index(keywrd, " "//line(:k)) + k + 1
                if (keywrd(k:k) < "A" .or. keywrd(k:k) > "Z") then
                  refkey(2) = refkey(2)(:i - 1)//refkey(2)(j + 1:)
                  i = i - 1
                end if
              end if
            end if
          end do
!
!  Check for " -" signs in keywrd line
!
          m = 0
          do
            i = index(keywrd(m + 1:), " -") + m
            if (i == m) exit
!
! Is the minus sign inside a quoted text
!
            j = 0
            do k = i, len_trim(keywrd)
              if (keywrd(k:k) == '"') j = j + 1
            end do
            if (mod(j,2) == 1) then
!
!  Yes, it's inside a quoted string, so protect it and carry on
!
              keywrd(i:i) = char(0)
              m = i
              cycle
            end if
            j = index(keywrd(i + 2:), " ") + i + 1
            do
              k = index(" "//refkey(2), " "//keywrd(i + 2: j - 1))
              if (k == 0) exit
              l = index(refkey(2)(k + 1:)," ") + k + 1
              refkey(2)(k:) = refkey(2)(l:)
            end do
            i = index(keywrd, " -")
            j = index(keywrd(i + 2:), " ") + i + 1
            keywrd(i:) = keywrd(j:)
            m = i
          end do
          do
            i = index(keywrd, char(0))
            if (i == 0) exit
            keywrd(i:i) = " "
          end do
          i = len_trim(keywrd)
          keywrd(i + 1:) = refkey(2)(:3000 - 1 - i)
          refkey(1) = trim(keywrd)
          refkey(2) = refkey(3)
!
! Delete SETUP keyword
!
          i = index(refkey(1)," SETUP")
          ch = " "
          do j = i + 6, len_trim(refkey(1))
            if (refkey(1)(j:j) == ch) exit
            if (refkey(1)(j:j) == '"') then
              if (ch == '"') then
                ch = " "
              else
                ch = '"'
              end if
            end if
          end do
          refkey(1) = refkey(1)(:i)//trim(refkey(1)(j + 1:))
        end if
        read (ir, '(A)', end=100, err=100) koment, title
      else if (index(allkey(1:i),' +') /= 0) then
!
!  READ SECOND KEYWORD LINE
!
        read (ir, '(A)', end=100, err=100) refkey(2)
        oldkey(len_trim(oldkey) + 2:) = trim(refkey(2))
        i = index(allkey(1:i),' +')
        keywrd(i:i + 1) = " "
        i = len_trim(keywrd)
        keywrd(i + 2:) = trim(refkey(2))
        call upcase (keywrd(i + 1:), len(keywrd) - i)
        k = len_trim(keywrd)
        allkey = keywrd(:k)
        if (index(keywrd,'SETUP') /= 0) then
          i = index(keywrd,'SETUP=')
          if (i /= 0) then
            j = index(keywrd(i:),' ')
            filen = oldkey(i+6:i+j-2)
            keywrd(i: i + j - 1) = " "
          else
            filen = 'SETUP'
            i = index(keywrd,'SETUP')
            keywrd(i: i + 5) = " "
          end if
          keywrd(i:i+6) = " "
          call add_path(filen)
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED', &
            position='REWIND')
          rewind isetup
          read (isetup, '(A)', end=30, err=30) refkey(2)
          close(isetup)
          i = len_trim(keywrd) + 1
          keywrd(i:) = refkey(2)(:1001 - i)
          call upcase (keywrd, len_trim(keywrd))
   30     continue
        else if (index(allkey(i + 1:),' +') /= 0) then
!
!  READ THIRD KEYWORD LINE
!
          read (ir, '(A)', end=100, err=100) refkey(3)
          allkey = refkey(3)
          i = index(allkey, " +")
          if (i /= 0) then
            write(iw,"(a)")" A maximum of three lines of keywords are allowed."
            write(iw,"(a)")" On the third line of keywords is a '+' sign, implying more lines of keywords."
            write(iw,"(a)")" Remove the '+' sign from the third line of keywords, and re-run."
            call web_message(iw,"plus.html")
            call mopend("A maximum of three lines of keywords are allowed.")
            return
          end if
          i = index(keywrd(1:len_trim(keywrd)),' +')
          keywrd(i:i + 1) = "  "
          i = len_trim(keywrd)
          keywrd(i + 2:) = refkey(3)(:1001 - i)
          i = len_trim(oldkey)
          oldkey(i:i) = " "
          oldkey(i+1:) = trim(refkey(3))
          call upcase (keywrd, len_trim(keywrd))
        end if
!
!  READ TITLE LINE
!
        read (ir, '(A)', end=100, err=100) koment, title
      else if (index(keywrd,'&') /= 0) then
        i = index(keywrd,'&')
        keywrd(i:i) = ' '
        oldkey(i:i) = ' '
        i = len_trim(keywrd)
        read (ir, '(A)', end=100, err=100) refkey(2)
        keywrd(i + 1:) = " "//refkey(2)(:1001 - i)
        oldkey(len_trim(oldkey) + 2:) = refkey(2)(:1001 - i)
        call upcase (keywrd, len_trim(keywrd))
        i = len_trim(keywrd)
        if (index(keywrd,'SETUP') /= 0) then
          i = index(keywrd,'SETUP=')
          if (i /= 0) then
            j = index(keywrd(i:),' ')
            filen = keywrd(i + 6:i + j)
            keywrd(i: i + j) = " "
          else
            filen = 'SETUP'
            i = index(keywrd,'SETUP')
            keywrd(i:i + 6) = " "
          end if
          call add_path(filen)
          open(unit=isetup, file=filen, status='UNKNOWN', form='FORMATTED', &
            position='REWIND')
          rewind isetup
          read (isetup, '(A)', end=39, err=40) keywrd(len_trim(keywrd) + 2:)
39        close (isetup)
          call upcase (keywrd, len_trim(keywrd))
          read (ir, '(A)', end=100, err=100) title
   40     continue
        else if (index(keywrd(i:),'&') /= 0) then
          j = index(keywrd,'&')
          keywrd(j:j) = ' '
          oldkey(j:j) = ' '
          read (ir, '(A)', end=100, err=100) refkey(3)
          read(refkey(3), '(a)') keywrd(i + 1:)
          read(refkey(3), '(a)') oldkey(i + 1:)
          call upcase (keywrd, len_trim(keywrd))
        else
          read (ir, '(A)', end=100, err=100) title
        end if
      else
        read (ir, '(A)', end=100, err=100) koment, title
      end if
      call split_keywords(oldkey)
      go to 60
50    continue
      if (zero_scf) go to 60
      numcal = 2
      call mopend ('SETUP FILE "'//trim(filen)//'" MISSING')
      write(iw,'(a)') " (Setup file name: '"//trim(filen)//"')"
      return
60    continue
      line = " "
      goto 99
100   continue
      call split_keywords(oldkey)
      
      if (numcal > 1) then
        if (index(keywrd,"OLDGEO") /= 0) return ! User forgot to add extra lines for title and comment
        if (aux) keywrd = " AUX"
        line = "JOB ENDED NORMALLY"
        if (quoted('GEO_DAT"') == 0) then
          line = ' ERROR IN READ OF FIRST THREE LINES' 
        else
          line = " "
        end if
      end if
99    continue
!
!  Look for non-standard characters.  If Apple's "text editor" is used,
!  convert the fancy '"' (four characters) into a normal '"'.
!
      exists = .false.
      j = len_trim(keywrd)
      do i = j, 3, -1
        if (ichar(keywrd(i:i)) /= 157) cycle
        if (ichar(keywrd(i - 1:i - 1)) == 128 .and. &
            ichar(keywrd(i - 2:i - 2)) == 226 .and. &
            ichar(keywrd(i - 3:i - 3)) == 156) then
          if (ichar(keywrd(i - 4:i - 4)) /= 128) then
            keywrd = keywrd(:i - 3)//'"'//trim(keywrd(i + 1:))
          else
            keywrd = keywrd(:i - 6)//'"'//trim(keywrd(i + 1:))
          end if
          exists = .true.
        end if
      end do
      do k = 1, 6
        if (index(refkey(k), " NULL") /= 0) exit
        j = len_trim(refkey(k))
        do
          i = index(refkey(k), char(160))
          if (i > 0) then
            refkey(k)(i:i) = " "
          else
            exit
          end if
        end do
      end do
      if (exists) then
        do k = 1, 6
          if (index(refkey(k), " NULL") /= 0) exit
          j = len_trim(refkey(k))
          do i = j, 3, -1
            if (ichar(refkey(k)(i:i)) /= 157) cycle
            if (ichar(refkey(k)(i - 1:i - 1)) == 128 .and. &
                ichar(refkey(k)(i - 2:i - 2)) == 226 .and. &
                ichar(refkey(k)(i - 3:i - 3)) == 156) then
              if (ichar(refkey(k)(i - 4:i - 4)) /= 128) then
                refkey(k) = refkey(k)(:i - 3)//'"'//trim(refkey(k)(i + 1:))
              else
                refkey(k) = refkey(k)(:i - 6)//'"'//trim(refkey(k)(i + 1:))
              end if
            end if
          end do
        end do
      end if
!
! Convert all bad slashes into good slashes
!
      do
        i = index(keywrd, backslash)
        if (i == 0) exit
        keywrd(i:i) = "/"
      end do
!
! Convert all fancy quotation marks into normal ASCII quotation marks
!
      do
        i = index(keywrd, char(148))
        if (i == 0) exit
        keywrd(i:i) = '"'
      end do
      return
  end subroutine gettxt
  character (len = 300) function get_text(line, i_start, zero)
    implicit none
    character :: line*(*)
    integer :: i_start, zero
!
!  Return text between character i_start and the next space.
!  If character i_start is '"' or ''', return text between character i_start + 1 and the closing character.
!
    integer, parameter :: num_lim = 2
    character :: limit(num_lim)*1, ch*1
    integer :: i, j
    data limit /"""","'"/
    do i = 1, num_lim
      if (line(i_start:i_start) == limit(i)) exit
    end do
    j = i_start
    if (i > num_lim) then
      ch = " "
      i = i_start
    else
      j = j + 1
      ch = limit(i)
      i = i_start + 1
    end if
    do
      if (line(i + 1:i + 1) == ch) exit
      i = i + 1
    end do
    get_text = line(j:i)
    if (zero == 0) line(i_start:i + 1) = " "
    return
  end function get_text
!
!
!
  subroutine split_keywords(oldkey)
  use molkst_C, only: keywrd, keywrd_quoted, line
  implicit none
  character, intent (in) :: oldkey*3000
  integer :: i, j, k, loop
  integer, parameter :: n_quoted = 4
  character :: quoted_keywords(n_quoted)*20, old*3000, quotation*1
  logical :: first
!
!  Split the keyword string "keywrd" into two strings: "keywrd" and "keywrd_quoted".
!  On exit:
!  First string:  keywrd        - holds all the keywords that must not contain quotation marks.
!  Second string: keywrd_quoted - holds all the keywords that can contain quotation marks.
!
!  If a keyword that can contain quotation marks does not have quotation marks, 
!  then quotation marks will be added.
!
!  If a keyword occurs more than once, only the first occurance will be kept in "keywrd_quoted", 
!  all others will be deleted.
!
  quoted_keywords(1) = ' GEO_DAT='
  quoted_keywords(2) = ' GEO_REF='
  quoted_keywords(3) = ' EXTERNAL='
  quoted_keywords(4) = ' SETUP='
  old = trim(oldkey)
  keywrd_quoted = " "
  do
    i = index(old, " + ")
    if (i /= 0) then
      line = trim(old)
      old = line(:i)//trim(line(i + 3:))
    else
      exit
    end if
  end do
  keywrd = trim(old)
  call upcase (keywrd, len_trim(keywrd))
  if (keywrd(1:1) /= " ") then
    line = trim(keywrd)
    keywrd = " "//trim(line)
    line = trim(old)
    old = " "//trim(line)
  end if
  do loop = 1, n_quoted
    first = .true.
    do
      i = index(keywrd, trim(quoted_keywords(loop)))
      if (i /= 0) then
!
! Found a keyword that should be quoted.  Make sure that it has a quotation mark
!
        k = len_trim(quoted_keywords(loop))
        j = i + k
        if (keywrd(j:j) /= '"') then
          quotation = '"'
          j = index(keywrd(j:), ' ') + j - 1
        else
          quotation = ' '
          j = index(keywrd(j + 1:), '"') + j + 1
          j = index(keywrd(j:), ' ') + j - 1
        end if
!
! Copy the keyword that contains quotation marks to keywrd_quoted
!    
        if (first) then
          line = trim(keywrd_quoted)
          keywrd_quoted = trim(line)//keywrd(i:i + k -1)//trim(quotation)//old(i + k:j - 1)//trim(quotation)  
          first = .false.
        end if
!
! Delete the keyword from keywrd that contains the quotation marks
!  
        line = trim(keywrd)
        keywrd = line(1:i - 1)//trim(line(j:))
        line = trim(old)
        old = line(1:i - 1)//trim(line(j:))
      else
        exit
      end if
    end do
  end do
  do
    i = index(keywrd_quoted, "\")
    if (i == 0) exit
    keywrd_quoted(i:i) = "/"
  end do
  return
  end subroutine split_keywords
  
  
  integer function quoted(key)
!
! If the test keyword "key" is present in "keywrd_quoted", then the function is set to
! the location in the string, and the string variable "line" is set to the text.
! If the test keyword is not present, the function is set to zero and "line" is empty.  
!
  use molkst_C, only: keywrd_quoted, line
  implicit none
  character, intent (in) :: key*(*)
  integer :: i, j, k
  integer, external :: end_of_keyword
!
! Search for a keyword that contains a string that contains text that could be mistaken
! for keywords.
!
    line = " "
    i = index(keywrd_quoted, trim(key))
    if (i /= 0) then
!
! Keyword contains quotation marks, so:
!
      k = Index (keywrd_quoted(i:), "=") + i
      j = end_of_keyword(keywrd_quoted, len_trim(keywrd_quoted), k) - 2
      i = i + 1 + len_trim(key) 
      if (trim(key) == "GEO_REF=") then
        if (index(keywrd_quoted(i:j), '"') /= 0) then
!
! Special treatment for "GEO_REF="
!
          do j = j, 1, -1
            if (keywrd_quoted(j:j) == '"') exit
          end do
        end if
      end if
      line = keywrd_quoted(i:j)
      end if
    quoted = i
    return
  end function quoted
