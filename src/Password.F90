subroutine password()
  use molkst_C, only: ijulian, line, academic, site_no, verson
  use reada_I
  use chanel_C, only : iw0
  implicit none
  integer :: date_of_creation
  integer :: i, j, jj, k, pass, p_new = 0, screen, ii(2)
  integer, parameter :: n_ac = 56, n_co = 77
  character :: mopac_pw*300, julian*8, academic_text(n_ac)*90, commercial_text(n_co)*90
  integer, external :: getenvqq
  logical :: exists, new, site_ok, opend
  integer, external :: iargc 
  character (len=8), external :: jdate
  data academic_text &
  
  /"Academic End User Software License Agreement", &
  "Important:  ", &
  "This License Agreement (""Agreement"") is a legal agreement between you, ", &
  "the end user (either an individual or an entity), and Stewart Computational ", &
  "Chemistry (""SCC"").  By typing the word ""Yes"" below, you signify that you ", &
  "have read the SCC License Agreement and accept its terms.  If you do not agree ", & 
  "to the Agreement's terms, please delete all copies of this Software.", &
  "Grant of License. ", & 
  "SCC grants, and you hereby accept, a non-exclusive license to install and use", &
  "the enclosed software product (""Software"") in accordance with the", &
  "terms of this Agreement.  This licensed copy of the Software may only be used ", &
  "at a single site. You may not: (a) electronically transfer the Software from", &
  "one site to another, (b) distribute copies of the Software to other sites, or ", &
  "(c) modify or translate the Software without the prior written consent of SCC.  ", &
  "The Software may be placed on a file or disk server connected to a network. ", &
  "You may make only those copies of the Software which are necessary ", &
  "to install and use it as permitted by this Agreement, or are for purposes of ", &
  "backup and archival records; all copies shall bear SCC's copyright and ", &
  "proprietary notices.  You may not make copies of any accompanying materials ", &
  "without prior, written notice from SCC.", &
  " ", &
  "Ownership.  ", &
  "The Software is and at all times shall remain the sole property of SCC or ", &
  "its Licensors.  This ownership is protected by the copyright laws of the ", &
  "United States and by international treaty provisions. You may not modify, ", &
  "decompile, reverse engineer, or disassemble the Software.", &
  " ", &
  "Assignment Restrictions.  ", &
  "You shall not rent, lease, or otherwise sublet the Software or any part thereof.  ", &
  " ", &
  "No Warranty.", &
  "SCC offers no warranties whatsoever.", &
  " ", &
  "SCC's Liability.  ", &
  "In no event shall SCC or its agents be liable for any indirect, special, or", &
  "consequential damages, such as, but not limited to, loss of anticipated profits", &
  "or other economic loss in connection with or arising out of the use of the ", &
  "software by you or the services provided for in this Agreement, even if SCC ", &
  "has been advised of the possibility of such damages.  ", &
  " ", &
  "No Other Warranties.  ", &
  "SCC and its agents disclaims other implied warranties, including, but not  ", &
  "limited to,implied warranties of merchantability or fitness for any purpose,  ", &
  "and implied warrantiesarising by usage of trade, course of dealing, or   ", &
  "course of performance. ", &
  " ", &
  "Governing Law.  ", &
  "This Agreement shall be construed according to the Laws of the State of Colorado.", &
  " Colorado.", &
  " ", &
  "----------------------------------------------------------------------------------", &
  "Stewart Computational Chemistry", &
  "15210 Paddington Circle, Colorado Springs, Colorado  80921-2512, USA.", &
  "URL: http://openmopac.net/                          mail: MrMOPAC@ATT.net", &
  "----------------------------------------------------------------------------------", &
  "++++" /
  data commercial_text /&
'Commercial End User Software License Agreement ', &
' ', &
'Important:  ', &
'This License Agreement ("Agreement") is a legal agreement between you, the end ', &
'user (either an individual or an entity), and Stewart Computational Chemistry ("SCC").  ', &
'By typing "Yes" below, you signify that you have read the SCC License Agreement ', &
'and accept its terms.  If you do not agree to the Agreement''s terms, please ', &
'delete all copies of this Software.', &
' ', &
'Grant of License.  ', &
'SCC grants, and you hereby accept, a non-exclusive license to install and use ', &
'one copy of the enclosed software product ("Software") in accordance with the ', &
'terms of this Agreement.  This licensed copy of the Software may only be used ', &
'on a single computer. You may not: (a) electronically transfer the Software from ', &
'one computer to another, (b) distribute copies of the Software to others, or ', &
'(c) modify or translate the Software without the prior written consent of SCC.  ', &
'The Software may be placed on a file or disk server connected to a network, ', &
'provided that a license has been purchased or legally acquired and registered ', &
'for every computer with access to that server or some license management software ', &
'is present and in use to allow only licensed users to access the Software.  ', &
'You may make only those copies of the Software which are necessary to install ', &
'and use it as permitted by this Agreement, or are for purposes of backup and ', &
'archival records; all copies shall bear SCC''s copyright and proprietary notices.  ', &
'You may not make copies of any accompanying materials without prior, written ', &
'notice from SCC.', &
' ', &
'Ownership.  ', &
'The Software is and at all times shall remain the sole property of SCC or its ', &
'Licensors.  This ownership is protected by the copyright laws of the United States ', &
'and by international treaty provisions.  Upon expiration or termination of this ', &
'Agreement, you shall promptly delete all copies of the Software and accompanying ', &
'materials.  You may not modify, decompile, reverse engineer, or disassemble the ', &
'Software.', &
' ', &
'Assignment Restrictions.  ', &
'You shall not rent, lease, or otherwise sublet the Software or any part thereof.  ', &
'You may transfer on a permanent basis the rights granted under this license provided ', &
'you transfer this Agreement and all copies of the Software, including prior versions, ', &
'and all accompanying materials (included in the product download).  The recipient ', &
'must agree to the terms of this Agreement in full and register this transfer with SCC.', &
' ', &
'SCC Limited Warranty.', &
'SCC''s sole warranty with respect to the Software is that it shall be free of errors ', &
'in program logic or documentation attributable to SCC and of errors that prevent the ', &
'performance of the principal computing functions of the Software.  SCC warrants this ', &
'for a period of ninety (90) days from the date of receipt of the Software.', &
' ', &
'SCC''s Liability.  ', &
!123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
'In no event shall SCC or its agents be liable for any indirect, special, or consequential', &
'damages, such as, but not limited to, loss of anticipated profits or other economic loss', &
'in connection with or arising out of the use of the software by you or the services ', &
'provided for in this Agreement, even if SCC has been advised of the possibility of ', &
'such damages.  SCC''s entire liability and your exclusive remedy shall be, at SCC''s ', &
'discretion, to return the Software and proof of purchase to SCC for either (a) return ', &
'of any license fee, or (b) correction or replacement of Software that does not meet ', &
'the terms of this limited warranty.', &
' ', &
'No Other Warranties.  ', &
'SCC and its agents disclaims other implied warranties, including, but not limited to,  ', &
'implied warranties of merchantability or fitness for any purpose, and implied warranties ', &
'arising by usage of trade, course of dealing, or course of performance.  Some states ', &
'do not allow the limitation of the duration or liability of implied warranties, so ', &
'the above restrictions might not apply to you.', &
' ', &
'Governing Law.  ', &
'This Agreement shall be construed according to the Laws of the State of Colorado.', &
' ', &
'----------------------------------------------------------------------------------', &
'Stewart Computational Chemistry', &
'15210 Paddington Circle, Colorado Springs, Colorado  80921-2512, USA.', &
'URL: http://openmopac.net/                          Email: MrMOPAC@ATT.net', &
'----------------------------------------------------------------------------------', &
"*******************************************************************************", &
"** Cite this work as: MOPAC2016, James J. P. Stewart, Stewart Computational  **", & 
"** Chemistry, Colorado, USA      web: HTTP://OpenMOPAC.net                   **", &
"*******************************************************************************", &
'++++'/
  if (site_no == 999) then 
    ijulian = 1000
    verson(7:7) = "W"
    return                                   !
  end if
  if (iw0 == 0) then
    screen = 6
  else
    screen = 0
  end if
  ii(1) = 0
  ii(2) = screen
  i = getenvqq("MOPAC_LICENSE",mopac_pw) 
  if (i /= 0) then
    i = len_trim(mopac_pw)
    if (mopac_pw(2:2) ==":" .or. index(mopac_pw, "$") /= 0) then  !  Is a Windows system
      if (mopac_pw(i:i) /= "\") then !  Add a terminal slash if it is not already present
        mopac_pw(i + 1:i + 1) ="\"
      end if
      verson(7:7) = "W"
      inquire (directory = trim(mopac_pw), exist = exists)
      if ( .not. exists) then
        write(screen,*)" Environmental variable ""MOPAC_LICENSE"" exists, but its value" 
        write(screen,*)" does not represent a folder on this computer"
        write(screen,*)" (An attempt was made to detect the folder '"//trim(mopac_pw)//"')"
        stop
      end if
      inquire (file = trim(mopac_pw)//"MOPAC2016.exe", exist = exists)
      if (exists) then
        mopac_pw= trim(mopac_pw)//"password for MOPAC2016"
      else
        write(screen,*)" Environmental variable ""MOPAC_LICENSE"" exists, but its value"
        write(screen,*)" does not represent a folder on this computer that contains MOPAC2016.exe"
        write(screen,"(a)")"(An attempt was made to detect the file '"//trim(mopac_pw)//"MOPAC2016.exe')"
        stop
      end if
    else if (index(mopac_pw, "/") /= 0) then ! It is a Linux system
      if (mopac_pw(i:i) /= "/") then !  Add a terminal slash if it is not already present
        mopac_pw(i + 1:i + 1) ="/"
      end if
      inquire (file = trim(mopac_pw)//"mopac2016.exe", exist = exists)
      if (exists) then
        verson(7:7) = "M"
      else
        inquire (file = trim(mopac_pw)//"MOPAC2016.exe", exist = exists)
        verson(7:7) = "L"
      end if
      if (exists) then
        mopac_pw= trim(mopac_pw)//"password_for_mopac2016"
      else
        write(screen,*)" Environmental variable ""MOPAC_LICENSE"" exists, and the directory it points to exists,"
        write(screen,*)" but the directory does not contain the MOPAC2016 executable"
        write(screen,*)" (MOPAC_LICENSE points to directory '"//trim(mopac_pw)//"')"
        stop
      end if   
    end if
  else
!
!  Environmental variable not set.  Test for defaults
!
    inquire (directory = "C:\", exist = exists)
    if (exists) then
!
! Platform is Windows
!  
      inquire (directory = "C:\program files\", exist = exists)
      if (exists) then
        inquire (directory = "C:\program files\mopac", exist = exists)
        if ( .not. exists) then
          write(screen,*) 
          write(screen,*) '  The MOPAC2016.exe file MUST be put in folder C:\Program Files\mopac'
          write(screen,*) '  unless the Environmental Variable MOPAC_LICENSE is set.'
          write(screen,*) ' '
          write(screen,*) '  If a different location is used, set the location in the Environmental Variable "MOPAC_LICENSE".'
          write(screen,*) '  For example, the Variable: MOPAC_LICENSE could be given the Value: M:\Utility.'
          call web_message(screen,"password_inaccessible.html")
          write(screen,*) " Correct this fault before continuing"
          read(5,"(a)")line
          stop
        end if
        mopac_pw   = 'C:\program files\mopac\password for MOPAC2016'
      else
!
!  The Operating System is WINDOWS, but the folder "program files" does not exist.
!  So put the password in the root folder, "C:\".
!
        mopac_pw   = 'C:\password for MOPAC2016'
      end if
      verson(7:7) = "W"
    else
!
! Check for Linux or Macintosh.  This is a poor test: Macintosh platforms are case-insensitive
! and Linux platforms are case-sensitive.  This is used in deciding whether to use "W" or "M".
! The test is poor because there is no guarantee that the case sensitivity will always be the same.
! but since the only use of the case is in the output, this is not very important.
!
      inquire (file = "/opt/mopac/mopac2016.exe", exist = exists)
      if (exists) then
        verson(7:7) = "M"
        mopac_pw = "/opt/mopac/password_for_mopac2016"
      else
        inquire (file = "/opt/mopac/MOPAC2016.exe", exist = exists)
        verson(7:7) = "L"
        mopac_pw = "/opt/mopac/password_for_mopac2016"
      end if
      if (.not. exists) then
!
! Platform is neither Linux nor Macintosh
!
        write(screen,'(a)')' The MOPAC executable must be put into a directory named "/opt/mopac"'
        write(screen,'(a)')' (If it cannot be put there, use environmental variable MOPAC_LICENSE to re-define the directory)'
        call web_message(screen,"trouble_shooting.html#default location")
        stop
      end if
    end if
  end if
  j = 0
 98 continue
  j = j + 1
  inquire (file = trim(mopac_pw), exist = exists)
  if (.not. exists) then
!
! Password for MOPAC2016 does not exist, so check for Password for MOPAC2012
!
    do k = len_trim(mopac_pw), 1, -2
      if (mopac_pw(k:k) == "/" .or. mopac_pw(k:k) == "\") exit
    end do
    i = index(mopac_pw(k:),"MOPAC20") + index(mopac_pw(k:),"mopac20") + 5 + k
    mopac_pw = mopac_pw(:i)//"12"
    inquire (file = trim(mopac_pw), exist = exists)
    if (exists) then
      write(screen,*)"An old license file for MOPAC2012 was found." 
      write(screen,*)"An attempt will be made to write a license for MOPAC2016."
      p_new = 76
      open(p_new, file = mopac_pw(:i)//"16", form="unformatted", iostat = j)
      if (j /= 0) then
        write(screen,'(/,a)')" The attempt to write the license file for MOPAC2016 failed."
        if (j == 9) then
          write(screen,*)"Permission to write to the file:", """"//mopac_pw(:i)//"16"" denied"
          write(screen,*)"For information on setting permissions, see:"
          write(screen,*)"HTTP://OpenMOPAC.net/Windows_set_permissions.html"
          write(screen,*)"(Delete this window when it is no longer wanted.)"
        end if
        do
          inquire(unit=0, opened=opend) 
          if (opend) call sleep(1)
        end do
        stop
      end if
    else
      mopac_pw = mopac_pw(:i)//"16"
    end if
  end if
  if (exists .and. j < 3) then    
    do jj = 1, 3
      open(unit=7, file=trim(mopac_pw), form="unformatted", iostat = i)
      if (i == 0) exit
      call sleep(1)
    end do
    if (i /= 0) then
      write(screen,*) 'Password file "' //trim(mopac_pw) //'" exists, ','but is currently inaccessible'
      call web_message(screen,"password_inaccessible.html")
      write(screen,*) " Correct this fault before continuing"
      read(5,"(a)")line
      stop
    end if  
    read(7, iostat=i)pass, site_no
    if (i /= 0) then
      if (i == -1) then
        close (7, STATUS = "DELETE",iostat=i)        
        if (i /= 0) then
          write(line,'(2a)')" Unable to delete corrupt password file in '", trim(mopac_pw)//"'"
          write(screen,'(//10x,a,/)')trim(line)
          write(screen,'(10x,a)')" Press 'Enter' to finish"
          call mopend(trim(line))
          read(*,*)
          stop
        end if
        goto 98
      else
        write(line,'(2a)')" Unable to read password file in '", trim(mopac_pw)//"'"
      end if
      write(screen,'(//10x,a,/)')trim(line)
      write(screen,'(10x,a)')" (To correct this fault, delete the file '"//trim(mopac_pw)//"' and try again.)"
      call mopend(trim(line))
      stop
    end if     
    new = .false.
    if  (p_new /= 0) then
      write(p_new, iostat=i)pass, site_no
      close (p_new)
    end if
  else
!
!  "exists" is true only if a password exists in the appropriate place.
!   
    new = .true.    
    line = " "
!
!
!                    Read in the password
!
!
!
    i = iargc()
    if (i == 0) then
    !  i = from_MAC_address()
      write(screen,'(//10x,a,/)')" MOPAC is a semiempirical quantum chemistry program"
      write(screen,'(10x,a)')" To install the MOPAC license, use the command"
      write(screen,'(10x,a)')" 'MOPAC2016.exe <license-key>'"
      write(screen,'(10x,a)')" e.g., 'MOPAC2016.exe 1234567a1234567'"
      call web_message(screen,"running_MOPAC.html")
      write(screen,'(10x,a)')" Press (return) to continue"
      read(5,*) 
      stop
    end if
    call getarg (1, line)
    do i = 1, 60 
      if(line(1:1) == " ")then
        line = line(2:)
      else if (index(line,"\") /= 0 .or. index(line,"/") /= 0) then
        j = index(line,"\")
        if (j > 0) line = line(j + 1:)
        j = index(line,"/")
        if (j > 0) line = line(j + 1:)
      else
        exit
      end if
    end do
    pass = nint(reada(line,1)) 
    if (pass > 99999999 .or. pass < 999999) then
      do i = 1, 2
        if (i == 2 .and. ii(i) == 0) exit
        write(ii(i),'(a)')" The licence key has not been read in,"
        write(ii(i),'(a)')" instead, an attempt was made to run MOPAC with the file """//trim(line)//"""" 
        write(ii(i),'(/,a)')" To activate MOPAC, run it using your personal license key."    
        write(ii(i),'(/,a)')" The command should be of the type ""MOPAC2016.exe 12345678a12345678""" 
        write(ii(i),'(a)')" or ""MOPAC2016.exe 12345678a12345678.mop"" (this file can be empty)" 
        write(ii(i),'(/,a)')" (Replace the key ""12345678a12345678"" with the license key you were given.)"
      end do
      call sleep(1000)  
      stop
    end if
    i = index(line,"a") + index(line,"A")
    if (i /= 0) then
      site_no = nint(reada(line, i))
      if (site_no > 99999999 .or. site_no < 999999) then
        write(screen,'(/,a)')" There must be 7 or 8 digits in the second set of numbers in the licence key" 
         if (screen == 6) then
        write(0,'(a)')" Licence key: """//trim(line)//""" read in" 
        write(0,'(a)')" (The licence key must be in the format ""12345678a12345678"")"    
        write(0,'(/,a)')" There must be 7 or 8 digits in the second set of numbers in the licence key" 
      end if
     call sleep (1000)
        stop
      end if
      i = Mod(pass,1523)
      j = Mod(site_no,1511)
    end if
    site_ok = .false.
    if (i == 0 .and. pass > 0 .and. j == 0 .and. site_no > 0)then
      if (screen == 6) &
      write(0,'(/,a)')" Minimize this window so you can see the license agreement." 
!
!  This is the Academic license
!
      do i = 1, n_ac
        if (Index(academic_text(i),"++++") /= 0) exit
        write(screen,"(a)")academic_text(i)(:len_trim(academic_text(i)))
        if (mod(i,52) == 0) then
          write(screen,*)" Press (return) to exit"
          read(5,*) 
        end if
      end do
      site_ok = .true.
    end if
    i = Mod(pass,1511)
    j = Mod(site_no,1523)      
    if (i == 0 .and. pass > 0 .and. j == 0 .and. site_no > 0)then
      if (screen == 6) &
      write(0,'(/,a)')" Minimize this window so you can see the license agreement." 
!
!  This is the Commercial license
!
      ijulian = ijulian + pass/1523
      if (ijulian < -59) then
        write(screen,"(///10x,a)")" This version of MOPAC is now out-of-date"
        close(7, status="delete")
        stop
      end if
      do i = 1, n_co
        if (Index(Commercial_text(i),"++++") /= 0) exit
        write(screen,"(a)")Commercial_text(i)(:len_trim(Commercial_text(i)))
        if (mod(i,52) == 0) then
          write(screen,*)" Press (return) to exit"
          read(5,*) 
        end if
      end do
      site_ok = .true.
    end if
    if ( .not. site_ok) then
      write(screen,"(//,a)")" A password file was supplied, but the eight-digit numerical codes were not"
      write(screen,"(a)")" recognized.  Check the digits of the numerical codes and try again"
      write(screen,"(a)")" "
      write(screen,"(a)")" Reminder: To install the password, run MOPAC with the file name equal to "
      write(screen,"(a)")" the password.  For example, if the password is 1234567a12345678, then run "
      write(screen,"(/,a)")" ""MOPAC2016.exe 1234567a12345678"""
      close(7, status="delete", iostat = i)
      write(screen,"(//,a)") " Password invalid"
      open(unit=7, file=trim(mopac_pw), form="unformatted", iostat = i)
      close(7, status="delete", iostat = i)
      if (screen == 6) &
      write(0,'(/,a)')" An error ocurred. Minimize this window so you can see the error message." 
      call sleep(1000)    
      stop
    end if
    write(screen,"(a,//,a)")"Scroll down to see the next part", &
      "Are these conditions are acceptable?"
    do
      write(screen,"(a)")" Type a ""Yes"" or ""No"" on the following line and press (return) to continue"
      read(5,"(a)")line
      call upcase(line,len_trim(line))
      if (index(line,"ES") /= 0) then
        exit
      else if (index(line,"O") /= 0) then
        write(screen,*)" Please send a message to MrMOPAC@ATT.net giving details", &
        " why the conditions are unacceptable."
        close(7, status="delete")
        stop
      end if
    end do
    open(unit=7, file=trim(mopac_pw), form="unformatted", iostat=i)
    if (i == 9) then
      write(screen,'(/)')
      write(screen,'(a)')" The file: ", " """//trim(mopac_pw)//""""," cannot be created - check permissions"
      if (verson(7:7) == "W") &
        write(screen,'(/,a,//)')" For instructions, see: http://openmopac.net/Windows_set_permissions.html"  
      write(screen,'(/,a)')" (Delete this window when it is no longer wanted.)"
      do
        inquire(unit=0, opened=opend) 
        if (opend) call sleep(1)
      end do
      stop     
    end if
    write(7, iostat=i)pass, site_no
    if (i /= 0) then
      write(screen,'(a)')" Cannot write password to file '"//trim(mopac_pw)//"'"
      do i = 1,100
        write(screen,'(a)')' '
      end do
      stop
    end if
    close(7, status="keep")
    write(screen,"(//a,////)")"         Password for MOPAC2016 successfully installed. Enjoy!"
    call sleep(3)
    stop
  end if
!
!   A valid password is now in the file mopac_pw
!        
!  On entry to password, ijulian is set to minus the number of days that have
!  passed since 1 January 2009
!
!  The password should be translate at the number of days since 1 January 2007 
!  plus the number of days the password is good for.
!  Academic sites get 366 days
!  Commercial evaluation sites get 31 days
!  All other commercial sites get 366 or more days
!
!
  do k = 1, 3
    open(unit=7, file=trim(mopac_pw), form="unformatted", iostat=i)
    if (i == 0) exit
    call sleep(1)
  end do
  if (i /= 0) then
    write(line,'(2a)')" Password file: '"//trim(mopac_pw)//"' exists, but it cannot be opened."
    write(screen,'(//10x,a,//)')trim(line)
    call mopend(trim(line))
    stop
  end if
  read(7, iostat=i)pass, site_no
  if (i < -1 .or. i > 0) then
    write(line,'(2a)')" Password file: '"//trim(mopac_pw)//"' exists, but it cannot be read."
    write(screen,'(//10x,a,//)')trim(line)
    call mopend(trim(line))
    stop
  end if
  call getdatestamp(line, julian)
  julian(3:) = julian(4:)
!
!  Set date_of_creation to the date of creation of the executable
!
  date_of_creation = ((ichar(julian(2:2)) - ichar('0'))*365 + &
                      (ichar(julian(3:3)) - ichar('0'))*100 + &
                      (ichar(julian(4:4)) - ichar('0'))*10  + &
                      (ichar(julian(5:5)) - ichar('0')))    + 1096
  julian = jdate() 
!
!  Set ijulian to minus the number of days lapsed since 1-Jan-2009
!
  ijulian = - (((ichar(julian(2:2)) - ichar('0'))*365 + &
                (ichar(julian(3:3)) - ichar('0'))*100 + &
                (ichar(julian(4:4)) - ichar('0'))*10 + &
                (ichar(julian(5:5)) - ichar('0'))) + 1096)
!
!  CHECK TO SEE IF THE LICENSE IS ACADEMIC
!
  i = Mod(pass,1523)
  j = Mod(site_no,1511)
  site_ok = .false.
  if (i == 0 .and. pass > 0 .and. j == 0 .and. site_no > 0)then
!
!  Prevent MOPAC running if other specific programs are present
!
! do i = 2009, 2015
!   write(line,'(a,i4,a)')"C:\Program files\Cambridgesoft\ChemOffice",i,"\Chem3D"
!   inquire (file=trim(line), exist = exists)
!   if (exists) exit
! end do
! if (i < 2016) then
!   write(screen,"(/10x,a)")" MOPAC requires a commercial license key to work with Chem3D,"
!   write(screen,"( 10x,a)")" and this key was not found. The free Academic license key "
!   write(screen,"( 10x,a)")" for stand-alone MOPAC is not valid for Chem3D. "
!   write(screen,"( 10x,a)")" The commercial license key is available at a discount to "
!   write(screen,"( 10x,a)")" Academics from: http://cacheresearch.com/mopac.html"
!   write(screen,"(/10x,a)")" Press press 'enter' to continue in Demo mode only (accepts  "
!   write(screen,"( 10x,a)")" only molecules with 1-12, 50-60, or 110-120 atoms)."
!   read(5,*)
!   ijulian = ijulian + 10000 !  
!   site_no = -1
!   return
! end if
!
!  This is the Academic license.  Set the executable to expire 365 days after the executable is constructed.
!
!
!  ijulian           = -(current number of days since the end of 2006)
!  date_of_creation  = number of days since the end of 2006 until the executable was made.
!
    ijulian = ijulian + date_of_creation + 365 
    academic = .true.
    site_no = site_no/1511 - 10000
    if (ijulian < 0) then
      write(screen,"(///,10(10x,a,/))") &
        " Your MOPAC executable, Version: "//verson//", has expired.", &
        " Please go to web-site: http://openmopac.net/Download_MOPAC_Executable_Step2.html to get a new version of MOPAC.", &
        " Do NOT request a new license key.  Academic license keys do not include time limits - the limit is in the program.", &
        " ", &
        " Press (enter) to continue."
      ijulian = 0
      read(5,'(i5)', iostat = i)k
    end if
    site_ok = .true.
  end if
! j = from_MAC_address()
!
!  CHECK TO SEE IF THE LICENSE IS COMMERCIAL
!
  i = Mod(pass,1511)
  if (j /= i) then
    j = i
  end if
  j = Mod(site_no,1523)
  if (i == 0 .and. pass > 0 .and. j == 0 .and. site_no > 0)then
!
!  This is the commercial license
!
    ijulian = ijulian + pass/1511 
    academic = .false.
    site_no = site_no/1523 - 10000
    if (ijulian < 0) then
      write(screen,"(///,10(10x,a,/))") &
 " The license for this copy of MOPAC has expired", &
 " Please go to web-site: HTTP://OpenMOPAC.net for instructions on", &
 " renewing the license", &
 " ", &
 " Press (return) to continue."
      read(5,'(i5)', iostat = i)k
    end if
    site_ok = .true.
    if (ijulian < -60) then
      write(screen,"(a)")
      write(screen,"(a)")" The time allowed by the password has expired."
      open(unit=7, file=trim(mopac_pw), form="unformatted")
      close(7, status="delete")
      stop
    end if      
  end if
  if (.not. site_ok) then
    write(screen,"(//,a)")" Before MOPAC2016 can be used, a password must be supplied"
    write(screen,"(/,a)")" To get a password, send an E-mail request to MrMOPAC@ATT.net"
    write(screen,"(a)")" "
    write(screen,"(a)")" To install the password, run MOPAC with the file name equal to the password"
    write(screen,"(a)")" For example, if the password is 123456789, then run ""MOPAC2016.exe 123456789"""
    write(screen,"(//,a)") "Password invalid"
    open(unit=7, file=trim(mopac_pw), form="unformatted")
    close(7, status="delete")
    stop
  end if
  if (new) then
    write(screen,"(//a,////)")"         Password for MOPAC2016 successfully installed. Enjoy!"
    stop
  end if 
end subroutine password
!function from_MAC_address()
!  use dfport, only : system
!  use molkst_C, only: verson
!  implicit none
!  integer :: from_MAC_address
!  character :: line*120, location*30, command*60
!  integer :: i, j, mac_no
!  location = "get_MAC_address.txt"
!
!  if (verson(7:7) == "W") then
!!
!!  Windows-specific code
!!
!    command = 'ipconfig /all  > '//trim(location)
!    j = SYSTEM(trim(command))
!    if (j /= 0) then
!      j = 0
!    end if
!    open(unit=8, file=trim(location), iostat = j)
!    if (j /= 0) then
!      j = 0
!    end if
!    do 
!      read(8,'(a)', iostat = j)line
!      if (j /= 0) then
!        j = 0
!      end if
!      if (index(line, "Physical Address") /= 0) exit
!    end do
!    close(8)
!    command = 'del /F /q '//trim(location)  
!    j = SYSTEM(trim(command))
!    if (j /= 0) then
!      j = 0
!    end if
!    line = line(index(line,":") + 2:)         
!    mac_no = 0
!    do i = 1, 6
!      mac_no = 10*mac_no + 100*ichar(line(i*3 - 2: i*3 - 2)) + 10*ichar(line(i*3 - 1:i*3 - 1))
!    end do
!    from_MAC_address = mod(mac_no,71) + 100*mod(mac_no,73) + 10000*mod(mac_no,79) + 1000000*mod(mac_no, 83)
!    return
!  end if
!        
!end function from_MAC_address
