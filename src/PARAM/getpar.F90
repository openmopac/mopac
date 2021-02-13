 subroutine getpar
!**********************************************************************
!
!   GETPAR 1: Reads in any reference parameters from a previous run.
!          2: Reads in the list of parameters to be optimized.
!
!**********************************************************************
    use param_global_C, only : valvar, valold, toplim, botlim,  &
    numvar, ifiles_8, locked, contrl, &
    maxpms, locvar
    use chanel_C, only: iw, ir
    use parameters_C, only : udd, zd, f0sd, g2sd, main_group, f0sd_store, g2sd_store, &
    dorbs, partyp, defmin, defmax, n_partyp_fn, n_partyp_alpb, n_partyp 
    implicit none
    character :: text*240, text2*240,koment*120, elemnt2*2
    logical :: buggy = .false., use_list
    integer :: i, j, k, l, jelmnt, io_stat, &
    & npar, mpar, nele, ipar, iele, n_ele, nlines
    character (len=2), dimension (107) :: elemnt
    character (len=7), dimension (2,1000) :: parlist
    external upcase 
    double precision, external :: reada
    intrinsic Index
    data (elemnt(i), i=1, 107) / "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", &
   & "F ", "NE", "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", &
   & "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", &
   & "AS", "SE", "BR", "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", &
   & "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", &
   & "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", &
   & "YB", "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", &
   & "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U ", "NP", "PU", &
   & "AM", "CM", "BK", "MI", "XX", "+3", "-3", "CB", "++", "+", "--", "-", "TV" /
    do i = 1,90
    if (main_group(i)) f0sd(i) = 0.d0
    if (main_group(i)) g2sd(i) = 0.d0
    if (f0sd_store(i) < 1.d-5) f0sd_store(i) = f0sd(i)
    if (g2sd_store(i) < 1.d-5) g2sd_store(i) = g2sd(i)
    dorbs(i) = (Abs(udd(i)) > 1.d-8 .and. zd(i) > 1.d-1) 
    end do
!
!  Read in parameters to be optimized
!
    numvar = 0
    locked = .false.
    nlines = 0
    if (index(contrl, " NOLIM") /= 0) then
      defmax = 400.d0
      defmin = -400.d0
    end if
!
!  Read in each parameter to be optimized, one parameter per line
!
    do
      read (ir, "(A240)", iostat=io_stat) text
      if(io_stat /= 0 .or. text == " ") goto 97
      call upcase (text, len_trim(text))
      nlines = nlines + 1

      if (Index (text, "END") /= 0) then
        close(unit=iw,status="DELETE")
        if (buggy) then
          write(ifiles_8,"(a)")" One or more parameters have an exactly zero value"
          write(ifiles_8,"(a)")" Either set value to non-zero or remove parameter"
          call finish
        end if
        if (numvar == 0 .and. nlines > 1) then
          write(ifiles_8,"(a)")" None of the parameters read in have non-zero values"
          write(ifiles_8,"(a)")" Either set one or more values to non-zero or remove parameters"
          call finish
        end if
        i = numvar
        do numvar = 1, i
           call print_par
        end do
        numvar = i
        return
      end if
!
!   Check for par symbols
!
      i = index(text, "PAR")
      if (i /= 0) then
        numvar = numvar + 1
        if (numvar > maxpms) then
          write (ifiles_8, "(A,I5)") " Too many parameters to be optimized. Max:", maxpms
          stop
        end if
        locvar(1, numvar) = 41 ! Parameter type is "PAR"
        locvar(2, numvar) = nint(reada(text, i))
        if (mpar == 1) locked(numvar) = .true.
        botlim(numvar) = -1.d3
        toplim(numvar) =  1.d3
        i = index(text(i:)," ") + i
        do i = i, 67
!
!  Check for bottom limit
!
          if (text(i:i) /= " ") then
            botlim(numvar) = reada (text, i)
            exit
          end if
        end do
!
!  Check for top limit
!
        i = Index (text(i:), " ") + i
        do i = i, 67
          if (text(i:i) /= " ") then
            toplim(numvar) = reada (text, i)
            exit
          end if
        end do
        call extract_parameter(locvar(1, numvar), locvar(2, numvar), valvar(numvar))
        valold(numvar) = valvar(numvar)
        goto 96
      end if
      if (index(text, "PAR") /= 0) then
        i = index(text,"PAR") + 4
        text2 = text
        j = ichar(text(i-1:i-1)) - ichar("0") 
        if (text(i:i) /= " ") j = j*10 + ichar(text(i:i)) - ichar("0")
        if (j < 13) then       
          k = mod(j - 1,3) + 1
          j = (j - 1)/3 + 1
          write(text,"('FN',2i1,' XX ',a)") k, j, trim(text2(i + 1:))
        else
          j = j - 12
          k = mod(j - 1,3) + 1
          j = (j - 1)/3 + 1
           write(text,"('FN',2i1,' MI ',a)") k, j, trim(text2(i + 1:))
         end if
      end if
      use_list = .false.
!
!  If "text" involves parentheses, parse the line to extract the different 
!  parameters to be optimized
!
      k = Index(text,"(") + 1
      if (k /= 1) then
        use_list = .true.
        do
          i = index(text," )") 
          if (i /= 0) then
            text(i:) = text(i + 1:)
          else
            exit
          end if
        end do
        l = Index(text,")")
        do n_ele = 1,84
          nele = n_ele
          do
            if (text(k:k) /= " ") exit
            k = k + 1
          end do
          j = Min(Index(text(k:)," ") + k - 1 , l-1)
          if (j == k - 1) j = l - 1
!
! "k" marks the start of the element symbol
! "j" marks the end of the element symbol
!
          if(j - k > 2 .or. j - k > 1 .and. text(j:j) /= " ")then
            if (text(k:j) == "ALL") then
              nele = 0
              do i = 1,57 
                nele = nele + 1
                parlist(1,nele) = elemnt(i)                 
              end do
              do i = 71, 83
                nele = nele + 1
                parlist(1,nele) = elemnt(i)
              end do
              goto 23
            else
!
!   Trap instances where chemical symbols are not separated by spaces
!
              write(ifiles_8,"(3a)")" Element symbol more than two characters, symbol: '",text(k:j),"'"
              call finish
            end if
          end if
          parlist(1,n_ele) = text(k:j)
          k = j + 1
          if (k >=l) exit
        end do
  23    text = text(l+1:)
        k = Index(text,"(") + 1
        l = Index(text,")")
        do npar = 1,50
          do
            if (text(k:k) /= " ") exit
            k = k + 1
          end do
          j = Min(Index(text(k:)," ") + k - 1 , l-1)
          if (j == k - 1) j = l - 1
!
! "k" marks the start of the parameter symbol
! "j" marks the end of the parameter symbol
!
          parlist(2,npar) = text(k:j)
          k = j + 1
          if (k >=l) exit
        end do   
      else
        npar = 1
        nele = 1
      end if
      text2 = text
      do iele = 1, nele ! iele loop
        if (Index(text2,"ALL") /= 0) then
            do j = 1, 99
            k = Index (" "//parlist(1,iele), " "//elemnt(j))
            if (k /= 0) then
              call what_parameters(j,parlist,mpar)
              exit
            end if
          end do
        else
          mpar = npar
        end if                  
        mpar_loop: do ipar = 1,mpar
          if (use_list) text = parlist(2,ipar)//" "//parlist(1,iele)
          koment = "    " // text (1:76)
          do i = 1, 80
            if (koment(1:1) /= " ") exit
            text = koment(2:80)
            koment = text(:120)
          end do
          i = Index (koment, " ")
          text = koment(1:i)
          i = len_trim(text)
          if (text(i:i) == ")") text(i:i) = " "
          numvar = numvar + 1
          if (numvar > maxpms) then
            write (ifiles_8, "(A,I5)") " Too many parameters to be optimized. Max:", maxpms
            stop
          end if
          do j = 1, n_partyp
            if (text(:5) == partyp(j)) go to 1500
          end do
          if (text == " ") then
            numvar = numvar - 1
            cycle
          end if
          go to 2000
    1500  locvar(1, numvar) = j
          if (mpar == 1) locked(numvar) = .true.
          jelmnt=0
          if (j == n_partyp_alpb .or. j == n_partyp_alpb + 1) then
!
!  This is a di-atomic parameter - read in the other element number
!
            i = Index(text, partyp(j))+5
            do j = 1, 99
              if (Index (" "//text(i:i+2), " "//elemnt(j)) /= 0) exit
            end do
            jelmnt = j
          end if
          i = Index (koment, " ")
          text = koment(i:)
          koment = text(:120)
          do i = 1, 50
            if (koment(1:1) /= " ") exit
            text = koment(2:50)
            koment = text(:120)
          end do
          text = " " // koment (1:49)
          botlim(numvar) = defmin(locvar(1, numvar))
          toplim(numvar) = defmax(locvar(1, numvar))
!Gallo                  
! Instead of setting the general upper and lower bounds to -1000 and 1000,
! respectively, they are defined in parameters_C.F90 also for par == T
!          if (par) then
!            botlim(numvar) = -1000.d0
!            toplim(numvar) =  1000.d0
!          end if
!Gallo                  
!
!  Read in lower and upper limits
!
          do i = 4, 67
!
!  Check for bottom limit
!
            if (text(i:i) /= " ") then
              botlim(numvar) = reada (text, i)
              exit
            end if
          end do
!
!  Check for top limit
!
          i = Index (text(i:), " ") + i
          do i = i, 67
            if (text(i:i) /= " ") then
              toplim(numvar) = reada (text, i)
              exit
            end if
          end do
          do j = 1, 100
            k = Index (text, " "//elemnt(j))
            if (k /= 0) go to 1900
          end do
          if (J == 101 .or. text /= " ") then
            write (ifiles_8, "(A)") " Element not recognized:"//text(1:len_trim(text)+1)
            stop
          end if
    1900  locvar(2, numvar) = j + jelmnt*200
!
!  Check: Is this parameter already marked for optimization?
!
          do i = 1, numvar - 1
              if (locvar(1,i) == locvar(1,numvar)) then
                k = locvar(2, i)/200
                if(k /= 0) then
                  l = locvar(2,i) - k*200
                  if ((k == j .and. l == jelmnt) .or. (k == jelmnt .and. l == j)) then
                    numvar = numvar -1
                    cycle mpar_loop
                  end if
                else
                  if (locvar(2,i) == j) then
                    numvar = numvar -1
                    cycle mpar_loop
                  end if
                end if
              end if
            end do
          call extract_parameter(locvar(1, numvar), locvar(2, numvar), valvar(numvar))
!
!  Is the parameter one that cannot be optimized, i.e. a transition metal with Gss - Hsp?
!
          if (locvar(2,numvar) < 200) then
            if (.not. main_group(locvar(2,numvar))) then
              if (locvar(1,numvar) > 9 .and. locvar(1,numvar) < 15) then
                numvar = numvar -1
                cycle
              end if
            end if
          end if
          if (Abs(valvar(numvar)) < 1.d-20) then
            numvar = numvar - 1
            cycle
          end if
          if (.not. buggy) buggy = (Abs(valvar(numvar)) < 1.d-6)
          i = locvar(2,numvar)
          if(i > 200) then
            j = i/200
            i = i - j*200
            elemnt2 = elemnt(j)
            if(elemnt2(2:2) == " ")elemnt2 = elemnt2(1:1)
          else
            elemnt2 = " "
          end if
!    call lockit(valvar(numvar), numvar)
!
!   Crude way to lock various parameters
!
          if (locvar(1,numvar) == 19 .and. locvar(2,numvar) == 17) then
            if (Abs(botlim(numvar) -0.5d0) < 1.d-3) botlim(numvar) = 2.d0! Lock Cl ZSN
          end if
          j = locvar(1,numvar) - n_partyp_fn + 1
          text(2:2) = " "
          if (i == 98) j = j + 12
          if (j < 10) then
            text(1:1) = char(j + ichar("0"))
          else if (j < 20) then
            text(2:2) = char(j - 10 + ichar("0"))
            text(1:1) = "1"
          else 
            text(2:2) = char(j - 20 + ichar("0"))
            text(1:1) = "2"
          end if
          valold(numvar) = valvar(numvar)
        end do mpar_loop
      end do ! iele loop
!
!  Check on atomic and diatomic parameters complete.
!
96 continue
    end do ! End of loop over lines read in
2000 write (ifiles_8, "(A)") text
    write (ifiles_8, "('   NAME NOT FOUND')")
    write (ifiles_8,*)" PARAM data set:"
    rewind (ir)
    do
     read (ir, "(A80)", end=199, err=199) text
     write (ifiles_8,'(a)')text(1:len_trim(text))
    end do!  All parameters loop
199 stop
97  return    
end subroutine getpar
subroutine what_parameters(z,parlist,mpar)
  use parameters_C, only : alpb
  implicit none
  integer :: z, mpar
  double precision alp
  character (len=7), dimension (2,1000) :: parlist
  integer :: i
  character (len=2), dimension (99) :: elemnt
  data elemnt /  "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE", &
   & "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", "SC", "TI", &
   & "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR&
   &",     "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", &
   & "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PD", &
   & "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF&
   &",     "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI"&
   & , "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U ", "NP", "PU", "AM", &
   & "CM", "BK", "CF", "XX" /
  parlist(2,1) = "USS"
  parlist(2,2) = "UPP"
  parlist(2,3) = "UDD"
  parlist(2,4) = "ZSN"
  parlist(2,5) = "ZPN"
  parlist(2,6) = "ZDN"
  parlist(2,7) = "ZS"
  parlist(2,8) = "ZP"
  parlist(2,9) = "ZD"
  parlist(2,10) = "BETAS"
  parlist(2,11) = "BETAP"
  parlist(2,12) = "BETAD"
  parlist(2,13) = "GSS"
  parlist(2,14) = "GSP"
  parlist(2,15) = "GPP"
  parlist(2,16) = "GP2"
  parlist(2,17) = "HSP"
  parlist(2,18) = "F0SD"
  parlist(2,19) = "G2SD"
  parlist(2,20) = "ALP" 
    parlist(2,21) = "FN11"
    parlist(2,22) = "FN21"
    parlist(2,23) = "FN31"
    parlist(2,24) = "FN12"
    parlist(2,25) = "FN22"
    parlist(2,26) = "FN32"
    parlist(2,27) = "FN13"
    parlist(2,28) = "FN23"
    parlist(2,29) = "FN33"
    parlist(2,30) = "FN14"
    parlist(2,31) = "FN24"
    parlist(2,32) = "FN34"
    mpar=33
  parlist(2,mpar) = "POC"
  do i =1,83
    alp = alpb(i,z)
    if (alp > 0.001d0) then
      mpar = mpar + 1
      parlist(2,mpar) = "ALPB_"//elemnt(i)
      mpar = mpar + 1
      parlist(2,mpar) = "XFAC_"//elemnt(i)
    end if
  end do

end subroutine what_parameters
subroutine lockit(value, numvar)
  implicit none
  double precision :: value
  integer :: numvar
! double precision :: multiplier = 1.d7
    value = value
    numvar = numvar
       return
!
!  Truncate value so that the last digits are all zero's
!
!     i = int(value)
!     value = i + nint((value - i)*multiplier)/multiplier      
!    if (.not. locked(numvar)) then
!
!  Modify value of parameter so that the last digit printed is not zero
!
!    value =  value + 0.2d0/multiplier
!  end if
  end subroutine lockit
  subroutine print_par
    use param_global_C, only : valvar, toplim, botlim, numvar, ifiles_8, locvar, &
      penalty, contrl
    use parameters_C, only : partyp 
!
!  Local
!
    double precision :: penalty_fn
    integer :: i, j
    logical :: lfirst = .true., l_prt
    character :: typtxt, text*5, elemnt2*2, num*1
       character (len=2), dimension (107) :: elemnt
    data (elemnt(i), i=1, 107) / "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", &
   & "F ", "NE", "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", &
   & "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", &
   & "AS", "SE", "BR", "KR", "RB", "SR", "Y ", "ZR", "NB", "MO", "TC", "RU", "RH", &
   & "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I ", "XE", "CS", "BA", "LA", &
   & "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", &
   & "YB", "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", &
   & "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", "PA", "U ", "NP", "PU", &
   & "AM", "CM", "BK", "MI", "XX", "+3", "-3", "CB", "++", "+", "--", "-", "TV" /
    save
    l_prt = (index(contrl, " SURVEY") == 0) 
    if (lfirst) then
      if (l_prt) then
        write (ifiles_8, "(//,10X,A)") "    PARAMETERS TO BE OPTIMIZED"
        write (ifiles_8, "(//5X, ' PARAMETER TYPE  ELEMENT    PARAMETER      LOWER LIMIT',&
            &'   UPPER LIMIT')")
      end if
      lfirst = .false.
    end if
    i = locvar(2,numvar)
    if (i > 200) then
      j = i/200
      i = i - j*200
      elemnt2 = elemnt(j)
      if(elemnt2(2:2) == " ")elemnt2 = elemnt2(1:1)
    else
      elemnt2 = " "
    end if
    typtxt = " "
    penalty_fn = penalty*((Max(0.d0, valvar(numvar)-toplim(numvar))+ &
    & Min(0.d0, valvar(numvar)-botlim(numvar))))**2
    if (penalty_fn > 1.d-6) typtxt = "*"
    if (l_prt) then
      if (penalty_fn > 1.d-6) then
        if (locvar(1, numvar) == 41) then
           write (ifiles_8, "(12X,A,11X,F16.8,A1,2F16.2)") &
          & text(:5),  valvar(numvar), "*", &
          & botlim(numvar), toplim(numvar)
        else
          write (ifiles_8, "(12X,A,7X,A,F16.8,A1,3F16.2)") &
          & partyp(locvar(1, numvar))//elemnt2, elemnt(i), valvar(numvar), typtxt, &
          & botlim(numvar), toplim(numvar), penalty_fn
        end if
      else
        if (locvar(1, numvar) == 41) then
          text = "PAR"
          if (i < 10) then
            num = "1"
          else
            num = "2"
          end if
          write(text(4:),'(i'//num//')')i
          write (ifiles_8, "(12X,A,11X,F16.8,A1,2F16.2)") &
          & text(:5),  valvar(numvar)
        else
          write (ifiles_8, "(12X,A,7X,A,F16.8,A1,2F16.2)") &
          & partyp (locvar(1,numvar))//elemnt2, elemnt(i), valvar(numvar), typtxt, &
          & botlim(numvar), toplim(numvar)
        end if
      end if
    end if
    return
    end subroutine print_par
