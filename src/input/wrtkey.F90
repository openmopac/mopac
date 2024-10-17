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

subroutine wrtkey
  use molkst_C, only : moperr, allkey, keywrd
  implicit none
  integer :: i, j, k, l, m
  integer, parameter :: n_protected_keywords = 14
  character :: protected_keywords(n_protected_keywords)*10
  data protected_keywords /"SITE=", "C.I.=", "I.D.=", "METAL=", "POLAR=", &
       "POLAR(", "PDB(", "OPEN(", "OPT(", "LARGE=", "EXTERNAL=", &
       "C.A.S.=", "C.I.D.=", "C.I.="/
!**********************************************************************
!
!  WRTKEY CHECKS ALL KEY-WORDS AND PRINTS THOSE IT RECOGNIZES.  IF IT
!  FINDS A WORD IT DOES NOT RECOGNIZE THE PROGRAM WILL BE STOPPED.
!
!**********************************************************************
!
!   Write out the control keywords. First, tidy up: delete equals signs and quotation marks,
!   and make sure there is a space at the start of the line.
!
!  Do not tidy up allkey earlier in the job, instead fill allkey here from keywrd,
!  and do all the tidying up at this one point.  The old style of tidying up allkey
!  as the job progressed was very error-prone and hard to debug.
!
  allkey = trim(keywrd)
  j = 1
  do
!
! Delete equals sign and everything after it that belongs to it
!
    i = index(allkey(j:), "=") + j
    if (i /= j) then
      if (allkey(i:i) == "(") then
        j = index(allkey(i + 1:), ') ') + i
      else if (allkey(i:i) == '"') then
        j = index(allkey(i + 1:), '"') + i
      else
        j = j + 1
        cycle
      end if
      allkey(i - 1:j) = " "
      j = i
    else
      do
        i = index(allkey, "(")
        if (i == 0) exit
        if (i > 6) then
          if (allkey(i-6:i) == "OUTPUT(") exit
        end if
        j = index(allkey(i + 1:), ')') + i
        allkey(i:j) = " "
      end do
      exit
    end if
  end do
  j=1
  do
    i = index(allkey(j:), '"') + j
    if (i /= j) then
      j = index(allkey(i + 1:), '"') + i
      allkey(i - 1:j) = " "
      j = j + 1
    else
      exit
    end if
  end do
!
!  Special treatment for keywords that must not be modified because they contain data that is
!  printed in the output.
!
  do k = 1, n_protected_keywords
    j = 1
    do
      i = index(keywrd(j:), " "//trim(protected_keywords(k))) + j
      if (i == j) exit
      m = i + len_trim(protected_keywords(k)) 
      if (keywrd(m:m) == "(" .or. keywrd(m - 1:m - 1) == "(") then
!
! Keyword starts with "(" so search for closing ")"
!
        j = 1
        do
          m = m + 1
          if (keywrd(m:m) == "(") j = j + 1
          if (keywrd(m:m) == ")") j = j - 1
          if (j == 0) exit
        end do
        j = m   
        l = 1000
      else
        l = index(keywrd(i:), ' ') 
        l = index(keywrd(i:l), '" ')    
        if (l == 0) exit
      endif
      if (j > 0 .and. l > 0) then
        j = min(j,l)
      else if (l > 0) then
        j = l
      end if
      allkey(i - 1:j) = keywrd(i - 1:j)
    end do
  end do
!
! Now print out everything to do with keywords, starting with control keywords
!
  call wrtcon (allkey)
  if (moperr) return
!
!   Write out the keywords controlling the working
!
  call wrtwor (allkey)
!
!   Write out the keywords controlling the output
!
  call wrtout (allkey)
!
!   Check the keywords to see if two or more conflict
!
  call wrtchk (allkey)
end subroutine wrtkey



subroutine wrtchk (allkey)
  use molkst_C, only: keywrd, is_PARAM, uhf, method_mndo, method_am1, &
  method_pm3, method_mndod, method_pm6, method_rm1, rhf, mozyme, line, &
  method_pm7, method_pm6_org, method_pm8, method_indo, koment, title, refkey
  use chanel_C, only: iw, input_fn
  implicit none
  logical :: birad, exci, ci, trip
  integer :: i, j, k, l, nmos
  character :: ch
  character (len=3000), intent (inout) :: allkey
  logical, external :: myword
  double precision, external :: reada
  birad = (Index (keywrd, " BIRAD") /= 0)
  exci  = (Index (keywrd, " EXCI") /= 0)
  ci    = (Index (keywrd, " C.I.") /= 0)
  trip  = (Index (keywrd, " TRIP") /= 0)
  uhf   = (Index (keywrd, " UHF") /= 0)
  rhf   = (Index (keywrd, " RHF") + index(keywrd, " MECI") /= 0  &
    .or. ci .or. (.not. uhf .and. Index(keywrd, " OPEN") /= 0))
  if (.not. (birad .or. exci .or. ci .or. (Index (keywrd, " INDO") /= 0) .or. (index(keywrd, " OPEN") /= 0))) then
    if ((index(keywrd, " MECI") /= 0) .or. (index(keywrd, " CIS") /= 0)) then
      if (index(keywrd, " CIS ") /= 0) call mopend("Keyword ""CIS"" used, but no C.I. defined")
      if (index(keywrd, " CISD ") /= 0) call mopend("Keyword ""CISD"" used, but no C.I. defined")
      if (index(keywrd, " CISDT ") /= 0) call mopend("Keyword ""CISDT"" used, but no C.I. defined")
      if ((index(keywrd, " MECI") /= 0)) call mopend("Keyword ""MECI"" used, but no C.I. defined")
      write(iw,'(/10x,a)')"Keyword ""C.I."" or ""OPEN"" must be present for the other C.I. keywords to be used"
      return
    end if
  end if
  if (.not. (mozyme .or. (index(keywrd," PDBOUT") + index(keywrd," RESID") + index(keywrd," ADD-H")/= 0)) &
    .and. index(keywrd, " 0SCF") == 0) then
!
! Check for keywords that require MOZYME to be present
!
    if (index(keywrd," PDBOUT") /= 0) then
      call mopend("Keyword PDBOUT only works when MOZYME or 0SCF is also present")
      return
    end if
    if (index(keywrd," RESID") /= 0) then
      call mopend("Keyword RESIDUES only works when MOZYME or 0SCF is also present")
      return
    end if
    if (index(keywrd," CVB") /= 0) then
      call mopend("Keyword CVB only works with MOZYME")
      return
    end if
    if (index(keywrd," SETPI") /= 0) then
      call mopend("Keyword SETPI only works with MOZYME ")
      return
    end if
  end if
  if (mozyme) then
    if (index(keywrd, " UHF") /= 0) then
      call mopend("Keyword UHF cannot be used with MOZYME")
    end if
    if (index(keywrd, " ENPART") /= 0)  then
      call mopend("Keyword ENPART is not available with MOZYME")
    end if
    if (index(keywrd, " LOCAL") /= 0)  then
      call mopend("Keyword LOCAL is not available with MOZYME")
    end if
    if (index(keywrd," RE-LOC") /= 0) then
      i = index(keywrd," RE-LOC")
      j = index(keywrd(i + 7:), " ") + i + 7
      if (index(keywrd(i:j), "=")  + index(keywrd," DENOUT") + &
        index(keywrd," VEC") + index(keywrd," ALLVE") == 0) then
        call mopend("Keyword RE-LOC only has meaning if DENOUT or VECTORS or ALLVEC is present")
      end if
      if (index(keywrd, "BANANA") /= 0) &
        call mopend("Keyword BANANA only works with keyword LOCAL, i.e. NOT with MOZYME")
      if (index(keywrd, "RABBIT") /= 0) &
        call mopend("Keyword RABBIT only works with keyword LOCAL, i.e. NOT with MOZYME")
    end if
    if (index(keywrd," 1SCF") /= 0 .and. index(keywrd," RAPID") /= 0) then
      call mopend("RAPID cannot be used with 1SCF")
      return
    end if
  end if
!
!   MULLIK does not work with UHF
!
  if (Index (keywrd, " MULLIK") /= 0 ) then
    if (uhf) then
      call mopend ("MULLIKEN POPULATION NOT AVAILABLE WITH UHF")
      return
    end if
  end if
!
!   Check UHF with nonallowed keywords
!
  if (uhf) then
    if (rhf) then
      call mopend("UHF and RHF cannot both be used")
    end if
    if (birad .or. exci .or. ci) then
      write (iw, "(//10X, ' UHF USED WITH EITHER BIRAD, EXCITED OR C.I. ')")
      write (iw, '(/ / 10 x, " IMPOSSIBLE OPTION REQUESTED,")')
      go to 1020
!
!  POLAR and UHF do not work
!
    else if (Index (keywrd, " POLAR") /= 0) then
      call mopend("POLAR does not work with UHF")
      go to 1020
    end if
  else if (exci .and. trip) then
    write (iw, "(//10X,' EXCITED USED WITH TRIPLET')")
    write (iw, 10000)
    go to 1020
  end if
!
!   PMEP only works with AM1
!
  if (Index (keywrd, " PMEP") /= 0 ) then
    if (.not. method_am1) then
      write (iw, "(A)") " PMEP only works with AM1"
      call mopend ("PMEP only works with AM1")
      return
    end if
  end if
!
!  XYZ and INT cannot simultaneously be present
!
  if (Index (keywrd, " INT ") /= 0 .and. Index (keywrd, " XYZ") /= 0) then
    call mopend ("INT cannot be used with XYZ")
  end if
!
!   Check T-PRIO MUST have DRC
!
  if (Index (keywrd, " T-PRIO") /= 0 .and. Index (keywrd, " DRC") == 0) then
    write (iw, "(//10X,'T-PRIO AND NO DRC')")
    write (iw, 10000)
    go to 1020
  end if
!
!   Check that only one method is used
!
    i = 0
    if (method_am1)     i = i + 1
    if (method_pm3)     i = i + 1
    if (method_pm6)     i = i + 1
    if (method_pm7)     i = i + 1
    if (method_PM8)     i = i + 1
    if (method_mndo)    i = i + 1
    if (method_mndod)   i = i + 1
    if (method_rm1)     i = i + 1
    if (method_indo)    i = i + 1
    if (i > 1) then
      write (iw, '(//10 x, " ONLY ONE OF MNDO, MNDOD, AM1, PM3, RM1, PM6, PM7, PM8, AND INDO ALLOWED")')
10000 format(/ / 10 x, "IMPOSSIBLE OPTION REQUESTED")
      call mopend("IMPOSSIBLE OPTION REQUESTED")
      return
    end if
      !
      !   Check that only one geometry option is used
      !
    i = 0
    line = trim(keywrd)
    j = 0
    l = len_trim(line)
    do
      j = j + 1
      if (j > l) exit
      if (line(j:j) == '"') then
        k = j
        do
          j = j + 1
          if (j > l) exit
          if (line(j:j) == '"') exit
        end do
        line(k:j) = " "
      end if
      if (j >= l) exit
    end do
    if (Index (line, " BFGS") /= 0)    i = i + 1
    if (Index (line, " LBFGS") /= 0)   i = i + 1
    if (Index (line, " EF") /= 0)      i = i + 1
    if (Index (line, " TS") /= 0)      i = i + 1
    if (Index (line, " SIGMA") /= 0)   i = i + 1
    if (Index (line, " NLLSQ") /= 0)   i = i + 1
    if (Index (line, " FORCE") + Index (line, " IRC") + Index (line, " DRC") /= 0)   i = i + 1
    if (i > 1) then
      call mopend ("MORE THAN ONE GEOMETRY OPTION HAS BEEN" // &
      & " SPECIFIED. CONFLICT MUST BE RESOLVED BEFORE JOB WILL RUN.")
      return
    end if
    if (Index (keywrd, " HESSIAN") /= 0) then
      if (Index (keywrd, " EF") == 0) then
        call mopend(" Keyword EF must be present if HESSIAN is used")
        return
      end if
    end if
    if (Index (keywrd, " PKA") /= 0 .and. .not. method_pm6) then
      write(iw,"(a)")" *"," *  The pKa option only works with PM6."
      write(iw,"(a)")" *  Either remove keyword pKa or change method to PM6."," *"
      call mopend("The pKa option only works with PM6. Calculation abandoned.")
        return
    end if
    if (Index(keywrd, " C.I.=") /= 0) then
      j = Index (keywrd, " C.I.=(")
      if (j /= 0) then
        nmos = Nint (reada (keywrd, Index (keywrd, "C.I.=(")+5))
      else
        nmos = Int (reada (keywrd, Index (keywrd, "C.I.")+4))
      end if
      if (nmos > 29 .and. .not. method_indo) then
        write (iw, "(' *',/,a)") " *      Maximum size of open space = 29 M.O.s"
        write (iw, "(a, i3,/,' *')") " *      Size requested:", nmos
        call mopend("Maximum size of open space = 29 M.O.s")
        return
      else if (nmos < 1) then
        write (iw, "(' *',/,a)") " *      Minimum size of open space = 1 M.O.s"
        write (iw, "(a, i3,/,' *')") " *      Size requested:", nmos
        call mopend("Minimum size of open space = 1 M.O.s")
        return
      end if
    end if
! INDO keywords incompatible with non-INDO methods
    if (.not. method_indo) then
      if (Index (keywrd, " C.I.D.") /= 0) then
        write (6,*) "C.I.D. is only compatible with INDO"
        call mopend("C.I.D. is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " C.A.S.") /= 0) then
        write (6,*) "C.A.S. is only compatible with INDO"
        call mopend("C.A.S. is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " MRCI") /= 0) then
        write (6,*) "MRCI is only compatible with INDO"
        call mopend("MRCI is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " MAXCI") /= 0) then
        write (6,*) "MAXCI is only compatible with INDO"
        call mopend("MAXCI is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " WRTCI") /= 0) then
        write (6,*) "WRTCI is only compatible with INDO"
        call mopend("WRTCI is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " TDIP") /= 0) then
        write (6,*) "TDIP is only compatible with INDO"
        call mopend("TDIP is only compatible with INDO")
        return
      end if
      if (Index (keywrd, " WRTCONF") /= 0) then
        write (6,*) "WRTCONF is only compatible with INDO"
        call mopend("WRTCONF is only compatible with INDO")
        return
      end if
! INDO incompatible with many standard MOPAC keywords
    else
! Incompatible CI options
      if (Index (keywrd, " CIS ") /= 0 .and. Index (keywrd, " CISD ") /= 0) then
        write (6,*) "CIS and CISD are incompatible keywords"
        call mopend("CIS and CISD are incompatible keywords")
        return
      end if
      if (Index (keywrd, " CIS ") /= 0 .and. Index (keywrd, " C.I.D.") /= 0) then
        write (6,*) "CIS and C.I.D. are incompatible keywords"
        call mopend("CIS and C.I.D. are incompatible keywords")
        return
      end if
! Too high spin
      if (Index (keywrd, " QUIN") /= 0 .or. Index (keywrd, " SEXT") /= 0 .or. &
          Index (keywrd, " SEPT") /= 0 .or. Index (keywrd, " OCTE") /= 0 .or. &
          Index (keywrd, " NONE") /= 0) then
        write (6,*) "Spin is incompatible with INDO"
        call mopend("Spin is incompatible with INDO")
        return
      end if

! Keywords in wrtcon - exit if a bad combination specified
      if (Index (keywrd, " FORCE") /= 0) then
        write (6,*) "FORCE is incompatible with INDO"
        call mopend("FORCE is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " LOCATE_TS") /= 0) then
        write (6,*) "LOCATE_TS is incompatible with INDO"
        call mopend("LOCATE_TS is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " MOZ") /= 0) then
        write (6,*) "MOZYME is incompatible with INDO"
        call mopend("MOZYME is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " LEWIS") /= 0) then
        write (6,*) "LEWIS is incompatible with INDO"
        call mopend("LEWIS is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " RESEQ") /= 0) then
        write (6,*) "RESEQ is incompatible with INDO"
        call mopend("RESEQ is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " CHARGES") /= 0) then
        write (6,*) "CHARGES is incompatible with INDO"
        call mopend("CHARGES is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " RAPID") /= 0) then
        write (6,*) "RAPID is incompatible with INDO"
        call mopend("RAPID is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " SITE") /= 0) then
        write (6,*) "SITE is incompatible with INDO"
        call mopend("SITE is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " PDBOUT") /= 0) then
        write (6,*) "PDBOUT is incompatible with INDO"
        call mopend("PDBOUT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " QMMM") /= 0) then
        write (6,*) "QMMM is incompatible with INDO"
        call mopend("QMMM is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " COMPAR") /= 0) then
        write (6,*) "COMPAR is incompatible with INDO"
        call mopend("COMPAR is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " BZ") /= 0) then
        write (6,*) "BZ is incompatible with INDO"
        call mopend("BZ is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " EXCI") /= 0) then
        write (6,*) "EXCI is incompatible with INDO"
        call mopend("EXCI is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " VELO") /= 0) then
        write (6,*) "VELO is incompatible with INDO"
        call mopend("VELO is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ESR") /= 0) then
        write (6,*) "ESR is incompatible with INDO"
        call mopend("ESR is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " NOMM") /= 0) then
        write (6,*) "NOMM is incompatible with INDO"
        call mopend("NOMM is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " MMOK") /= 0) then
        write (6,*) "MMOK is incompatible with INDO"
        call mopend("MMOK is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " CISDT") /= 0) then
        write (6,*) "CISDT is incompatible with INDO"
        call mopend("CISDT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " COSCCH") /= 0) then
        write (6,*) "COSCCH is incompatible with INDO"
        call mopend("COSCCH is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " INVERT") /= 0) then
        write (6,*) "INVERT is incompatible with INDO"
        call mopend("INVERT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ESP") /= 0) then
        write (6,*) "ESP is incompatible with INDO"
        call mopend("ESP is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " NSURF") /= 0) then
        write (6,*) "NSURF is incompatible with INDO"
        call mopend("NSURF is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " SCINCR") /= 0) then
        write (6,*) "SCINCR is incompatible with INDO"
        call mopend("SCINCR is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " SLOPE") /= 0) then
        write (6,*) "SLOPE is incompatible with INDO"
        call mopend("SLOPE is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " PDBOUT") /= 0) then
        write (6,*) "PDBOUT is incompatible with INDO"
        call mopend("PDBOUT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " NOTER") /= 0) then
        write (6,*) "NOTER is incompatible with INDO"
        call mopend("NOTER is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " XENO") /= 0) then
        write (6,*) "XENO is incompatible with INDO"
        call mopend("XENO is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " CHAINS") /= 0) then
        write (6,*) "CHAINS is incompatible with INDO"
        call mopend("CHAINS is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " SETPI") /= 0) then
        write (6,*) "SETPI is incompatible with INDO"
        call mopend("SETPI is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " PINOUT") /= 0) then
        write (6,*) "PINOUT is incompatible with INDO"
        call mopend("PINOUT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " REORTH") /= 0) then
        write (6,*) "REORTH is incompatible with INDO"
        call mopend("REORTH is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " RE-LOCAL") /= 0) then
        write (6,*) "RE-LOCAL is incompatible with INDO"
        call mopend("RE-LOCAL is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ADD_H") /= 0) then
        write (6,*) "ADD_H is incompatible with INDO"
        call mopend("ADD_H is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " CONNOL") /= 0) then
        write (6,*) "CONNOL is incompatible with INDO"
        call mopend("CONNOL is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ESPRST") /= 0) then
        write (6,*) "ESPRST is incompatible with INDO"
        call mopend("ESPRST is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " UHF") /= 0) then
        write (6,*) "UHF is incompatible with INDO"
        call mopend("UHF is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " STATIC") /= 0) then
        write (6,*) "STATIC is incompatible with INDO"
        call mopend("STATIC is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " POLAR") /= 0) then
        write (6,*) "POLAR is incompatible with INDO"
        call mopend("POLAR is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " RESI") /= 0) then
        write (6,*) "RESI is incompatible with INDO"
        call mopend("RESI is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " NOOPT") /= 0) then
        write (6,*) "NOOPT is incompatible with INDO"
        call mopend("NOOPT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " OPT") /= 0) then
        write (6,*) "OPT is incompatible with INDO"
        call mopend("OPT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " P=") /= 0) then
        write (6,*) "P= is incompatible with INDO"
        call mopend("P= is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " CUBE") /= 0) then
        write (6,*) "CUBE is incompatible with INDO"
        call mopend("CUBE is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ESPGRID") /= 0) then
        write (6,*) "ESPGRID is incompatible with INDO"
        call mopend("ESPGRID is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " MICROS") /= 0) then
        write (6,*) "MICROS is incompatible with INDO"
        call mopend("MICROS is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " OPEN") /= 0) then
        write (6,*) "OPEN is incompatible with INDO"
        call mopend("OPEN is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " MS=") /= 0) then
        write (6,*) "MS= is incompatible with INDO"
        call mopend("MS= is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " ROOT") /= 0) then
        write (6,*) "ROOT is incompatible with INDO"
        call mopend("ROOT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " STEP") /= 0) then
        write (6,*) "STEP is incompatible with INDO"
        call mopend("STEP is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " MERS") /= 0) then
        write (6,*) "MERS is incompatible with INDO"
        call mopend("MERS is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " POINT") /= 0) then
        write (6,*) "POINT is incompatible with INDO"
        call mopend("POINT is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " IRC") /= 0) then
        write (6,*) "IRC is incompatible with INDO"
        call mopend("IRC is incompatible with INDO")
        return
      end if
      if (Index (keywrd, " DRC") /= 0) then
        write (6,*) "DRC is incompatible with INDO"
        call mopend("DRC is incompatible with INDO")
        return
      end if

! Keywords in wrtout - print a warning but let the calculation proceed
      i = 0
      if (Index (keywrd, " AUX")      /= 0) i = i + 1
      if (Index (keywrd, " OUTPUT")   /= 0) i = i + 1
      if (Index (keywrd, " NOCOM")    /= 0) i = i + 1
      if (Index (keywrd, " ISOTOPE")  /= 0) i = i + 1
      if (Index (keywrd, " DENOUT")   /= 0) i = i + 1
      if (Index (keywrd, " OLDENS")   /= 0) i = i + 1
      if (Index (keywrd, " CIOSCI")   /= 0) i = i + 1
      if (Index (keywrd, " ENPART")   /= 0) i = i + 1
      if (Index (keywrd, " GRAD")     /= 0) i = i + 1
      if (Index (keywrd, " MINI")     /= 0) i = i + 1
      if (Index (keywrd, " VEC")      /= 0) i = i + 1
      if (Index (keywrd, " EIGEN")    /= 0) i = i + 1
      if (Index (keywrd, " CARTAB")   /= 0) i = i + 1
      if (Index (keywrd, " CHARST")   /= 0) i = i + 1
      if (Index (keywrd, " DCART")    /= 0) i = i + 1
      if (Index (keywrd, " DERI1")    /= 0) i = i + 1
      if (Index (keywrd, " DERI2")    /= 0) i = i + 1
      if (Index (keywrd, " DERITR")   /= 0) i = i + 1
      if (Index (keywrd, " DERIV")    /= 0) i = i + 1
      if (Index (keywrd, " DERIVO")   /= 0) i = i + 1
      if (Index (keywrd, " MAKVEC")   /= 0) i = i + 1
      if (Index (keywrd, " DIIS")     /= 0) i = i + 1
      if (Index (keywrd, " FLEPO")    /= 0) i = i + 1
      if (Index (keywrd, " FMAT")     /= 0) i = i + 1
      if (Index (keywrd, " DFORCE")   /= 0) i = i + 1
      if (Index (keywrd, " HCORE")    /= 0) i = i + 1
      if (Index (keywrd, " FREQCY")   /= 0) i = i + 1
      if (Index (keywrd, " LINMIN")   /= 0) i = i + 1
      if (Index (keywrd, " RAMA")     /= 0) i = i + 1
      if (Index (keywrd, " SYMOIR")   /= 0) i = i + 1
      if (Index (keywrd, " GROUP")    /= 0) i = i + 1
      if (Index (keywrd, " SYMTRZ")   /= 0) i = i + 1
      if (Index (keywrd, " DENS")     /= 0) i = i + 1
      if (Index (keywrd, " SPIN")     /= 0) i = i + 1
      if (Index (keywrd, " PRNT")     /= 0) i = i + 1
      if (Index (keywrd, " DISP")     /= 0) i = i + 1
      if (Index (keywrd, " DEP")      /= 0) i = i + 1
      if (Index (keywrd, " BONDS")    /= 0) i = i + 1
      if (Index (keywrd, " ALLBO")    /= 0) i = i + 1
      if (Index (keywrd, " SYBYL")    /= 0) i = i + 1
      if (Index (keywrd, " FOCK")     /= 0) i = i + 1
      if (Index (keywrd, " LARGE")    /= 0) i = i + 1
      if (Index (keywrd, " HTML")     /= 0) i = i + 1
      if (Index (keywrd, " PDBOUT")   /= 0) i = i + 1
      if (Index (keywrd, " AIGOUT")   /= 0) i = i + 1
      if (Index (keywrd, " GRAP")     /= 0) i = i + 1
      if (Index (keywrd, " 1ELEC")    /= 0) i = i + 1
      if (Index (keywrd, " SUPER")    /= 0) i = i + 1
      if (Index (keywrd, " ALLVEC")   /= 0) i = i + 1
      if (Index (keywrd, " H-PRIO")   /= 0) i = i + 1
      if (Index (keywrd, " X-PRIO")   /= 0) i = i + 1
      if (Index (keywrd, " T-PRIO")   /= 0) i = i + 1
      if (Index (keywrd, " MECI")     /= 0) i = i + 1
      if (Index (keywrd, " HESSIAN")  /= 0) i = i + 1
      if (Index (keywrd, " LOCAL")    /= 0) i = i + 1
      if (Index (keywrd, " MULLIK")   /= 0) i = i + 1
      if (Index (keywrd, " PI")       /= 0) i = i + 1
      if (Index (keywrd, " POTWRT")   /= 0) i = i + 1
      if (Index (keywrd, " PMEP")     /= 0) i = i + 1
      if (Index (keywrd, " PRTMEP")   /= 0) i = i + 1
      if (Index (keywrd, " MINMEP")   /= 0) i = i + 1
      if (Index (keywrd, " PL")       /= 0) i = i + 1
      if (Index (keywrd, " PKA")      /= 0) i = i + 1
      if (Index (keywrd, " POPS")     /= 0) i = i + 1
      if (Index (keywrd, " BCC")      /= 0) i = i + 1
      if (Index (keywrd, " Z=")       /= 0) i = i + 1
      if (Index (keywrd, " HYPERF")   /= 0) i = i + 1
      if (Index (keywrd, " THERMO")   /= 0) i = i + 1
      if (Index (keywrd, " K=")       /= 0) i = i + 1
      if (i > 0) then
        write (6,*) "WARNING: Specified print options may not behave as expected with INDO"
        write (iw,*) "**********"
        write (iw,*) "*  WARNING: Specified print options may not behave as expected with INDO"
        write (iw,*) "**********"
      end if
    end if
       !******************************************************************
       !
       !  Check to see that all keywords have been recognized
       !
       !*******************************************************************
       !
       !    DUMMY IF STATEMENT TO REMOVE AMPERSAND, PLUS SIGNS AND OBSOLETE KEYWORDS, IF PRESENT
       !
      if (myword(allkey, " SETUP")    ) i = 1
      if (myword(allkey, " &")        ) write (iw,'(" *  &          - THIS IS A DEPRECATED KEYWORD, USE "" ++ "" INSTEAD")')
      if (myword(allkey, " +")        ) write (iw,'(" *  +          - THIS IS A DEPRECATED KEYWORD, USE "" ++ "" INSTEAD")')
      if (myword(allkey, " CONTROL")  ) i = 2
      if (myword(allkey, " DIIS")     ) write (iw,'(" *  DIIS       - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " NODIIS")   ) write (iw,'(" *  NODIIS     - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " ROT")      ) write (iw,'(" *  ROT        - THIS IS AN OBSOLETE KEYWORD, IT WILL BE IGNORED")')
      if (myword(allkey, " HEADER") .or. myword(allkey, " USER"))  then
        write (iw,'(" *  HEADER     - DATA SET IS IN PROTEIN DATA BANK FORMAT")')
        i = index(keywrd, " HEADER")
        keywrd = " ADD-H PDBOUT "//keywrd(:i)
        refkey(1) = trim(keywrd)
        refkey(2) = " NULL "
        koment = " From PDB file: '"//input_fn(1:len_trim(input_fn)   - 5)//"'"
        title = "   "
        return
      end if

      if (allkey == " ") return
!
!  Compress unrecognized key-words
!
   if (allkey /= ' ' .and. .not. is_PARAM) then
        j = 0
        do i = 1, 3000 - 1
          if (allkey(i:i+1) == '  ') cycle
          j = j + 1
          ch = allkey(i:i)
          allkey(j:j) = ch
        end do
        j = max(1,j)
        l = index(keywrd,' DEBUG')
        if (l /= 0) then
          write (iw, '('' *  DEBUG KEYWORDS USED:  '',A)') allkey(2:j)
        else
          write (line, "('UNRECOGNIZED KEY-WORDS: (',A,')')") allkey (2:j)
          call to_screen(trim(line))
          call mopend (trim(line))
          call mopend ("IF THESE ARE DEBUG KEYWORDS, ADD THE KEYWORD ""DEBUG"".")
        end if
      return
    end if
1020 if (is_PARAM) return
  call mopend ("CALCULATION ABANDONED, SORRY!")
end subroutine wrtchk


subroutine wrtcon (allkey)
  use molkst_C, only: keywrd, keywrd_quoted, numat, pressure, id, mozyme, mers, natoms, maxtxt, txtmax, &
    line, old_chrge
  use cosmo_C, only : fepsi, nspa
  use chanel_C, only: iw, job_fn
  use meci_C, only : nmos
  use conref_C, only : fpcref
  use common_arrays_C, only :  lopt, na, txtatm
  implicit none
  character (len=3000), intent (out) :: allkey
  double precision :: epsi, sum
  integer :: i, ielec, ilevel, j, k, l
  logical :: l_add_H = .false., l_temp
  character :: num*1, num1*1
  character(len=300), external :: get_a_name
  integer, external :: quoted
  logical, external :: myword
  double precision, external :: reada
  save :: l_add_H
!
!  COSMO does not work with polymers or other infinite systems
!
  if (id > 0 .and. Index(keywrd, "EPS=") /= 0) call l_control("EPS=", len_trim("EPS="), -1)
!
  i = index(keywrd, " LOCATE_TS")
  if (i /= 0) keywrd(i:i + 9) = " LOCATE-TS"
  l_temp = (index(keywrd," MOZ") + index(keywrd," LEWIS") + index(keywrd," LOCATE-TS") + &
    index(keywrd, " RESEQ") + index(keywrd, " CHARGES") + &
    index(keywrd, " RAPID") + index(keywrd, " SITE=(") /= 0)
  mozyme = (l_temp .or. index(keywrd, " PDBOUT") /= 0)
  if (index(keywrd, "CHARGE=") == 0 .and. old_chrge /= 0 .and. index(keywrd, " OLDGEO") /= 0) then
!
!  Use the charge from the previous calculation
!
    i = index(keywrd,"           ")
    if (old_chrge > 99) then
      write(keywrd(i:i+11),'(" CHARGE=",i3)')old_chrge
    else if (old_chrge > 9) then
      write(keywrd(i:i+11),'(" CHARGE=",i2)')old_chrge
    else if (old_chrge > 0) then
      write(keywrd(i:i+11),'(" CHARGE=",i1)')old_chrge
    else if (old_chrge > -10) then
      write(keywrd(i:i+11),'(" CHARGE=",i2)')old_chrge
    else
      write(keywrd(i:i+11),'(" CHARGE=",i3)')old_chrge
    end if
   allkey(i:i+11) = keywrd(i:i+11)
  else
    old_chrge=0
  end if
  if (index(keywrd, " COMPAR") /= 0) &
    call l_control("0SCF HTML GEO-OK LET NOCOM", len_trim("0SCF HTML GEO-OK LET NOCOM"), -1)
  if (myword(allkey, " NEXT"))   write (iw, '(" *  NEXT       - DO NOT USE A BLANK LINE AFTER THE PREVIOUS JOB")')
  if (myword(allkey, " CCDC "))  i = 0 ! Dummy assignment   - to clear CCDC
  if (myword(allkey, " MNDO "))  write (iw, '(" *  MNDO       - The MNDO Hamiltonian to be used")')
  if (myword(allkey, " AM1 "))   write (iw, '(" *  AM1        - The AM1 Hamiltonian to be used")')
  if (myword(allkey, " PM3 "))   write (iw, '(" *  PM3        - The PM3 Hamiltonian to be used")')
  if (myword(allkey, " MNDOD"))  write (iw, '(" *  MNDOD      - The MNDOD Hamiltonian to be used")')
  if (myword(allkey, " PM6 "))   write (iw, '(" *  PM6        - The PM6 Hamiltonian to be used")')
  if (myword(allkey, " PM7-TS")) write (iw, '(" *  PM7-TS     - Calculate barrier height using PM7-TS")')
  if (myword(allkey, " PM7"))    write (iw, '(" *  PM7        - The PM7 Hamiltonian to be used")')
  if (myword(allkey, " PM6-ORG"))write (iw, '(" *  PM6-ORG    - The PM6-ORG Hamiltonian to be used")')
  if (myword(allkey, " PM8"))    write (iw, '(" *  PM8        - The PM8 Hamiltonian to be used (IN DEVELOPMENT)")')
  if (myword(allkey, " SPARKL")) write (iw, '(" *  SPARKLE    - Use SPARKLES when they exist.")')
  if (myword(allkey, " RM1 "))   write (iw, '(" *  RM1        - The RM1 Hamiltonian to be used")')
  if (myword(allkey, " PM5 "))   write (iw, '(" *  PM5        - METHOD NOT SUPPORTED. DEFAULT METHOD USED INSTEAD (See above)")')
  if (myword(allkey, " INDO "))  write (iw, '(" *  INDO       - The INDO Hamiltonian to be used")')
  if (myword(allkey,' MRCI'))    write (iw, '(" *  MRCI       -  Use multi-reference CI (specify ref dets using C.A.S.)")')
  if (myword(allkey,' MAXCI='))  write (iw, '(" *  MAXCI =",i5," maximum CI states computed")') &
    nint(reada(keywrd,index(keywrd, ' MAXCI')))
  if (myword(allkey,' WRTCI='))  write (iw, '(" *  WRTCI =",i5," maximum CI states printed")') &
    nint(reada(keywrd,index(keywrd, ' WRTCI')))
  if (myword(allkey,' TDIP'))    write (iw, '(" *  TDIP       -  Print transition dipoles between excited states")')
  if (myword(allkey,'WRTCONF=')) &
    write (iw, '(" *  WRTCONF =",f8.4," - Print configurations with coefficients greater than cutoff")') &
    reada(keywrd,index(keywrd, ' WRTCONF'))
!
!  The lack of space before QMMM on the next line is deliberate   - at allows MOL_QMMM as an option
!
  if (myword(allkey, "QMMM "))   write (iw, '(" *  QMMM       - Generate energies and gradients for use in MM codes")')
  if (myword(allkey, " COMPAR")) write (iw, '(" *  COMPARE    - Compare two geometries")')
  if (myword(allkey, " BZ"))     write (iw, '(" *  BZ         - Write file <name>.brz for use by program BZ")')
  if (myword(allkey, " BIRAD"))  write (iw, '(" *  BIRADICAL  - SYSTEM HAS TWO UNPAIRED ELECTRONS")')
  if (myword(allkey, " EXCI"))   write (iw, '(" *  EXCITED    - FIRST EXCITED STATE IS TO BE OPTIMIZED")')
  if (myword(allkey, " VELO"))   write (iw, '(" *  VELOCITY   - INPUT STARTING VELOCITIES FOR DRC")')
  if (myword(allkey, " DIPO"))   write (iw, '(" *  DIPOLE     - PRINT DIPOLE INSTEAD OF ENERGY IN TRAJECTORIES")')
  if (myword(allkey, " GEO-OK")) write (iw, '(" *  GEO-OK     - OVERRIDE INTERATOMIC DISTANCE AND OTHER SAFETY CHECKS")')
  if (myword(allkey, " CHECK"))  write (iw, '(" *  CHECK      - RUN EXTRA INTERATOMIC DISTANCE CHECKS")')
  if (myword(allkey, " PM6-D")) then
    if (index(keywrd, "PM6-DH+") /= 0) then
      write (iw, '(" *  PM6-DH+    - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING PM6-DH+")')
    elseif (index(keywrd, "PM6-DH2X") /= 0) then
      write (iw, '(" *  PM6-DH2X   - CORRECT DISPERSION, HALOGEN AND HYDROGEN BOND TERMS USING PM6-DH2X")')
    elseif (index(keywrd, "PM6-D2X") /= 0) then
      write (iw, '(" *  PM6-D2X    - CORRECT DISPERSION AND HALOGEN BOND TERMS USING PM6-D2X")')
    elseif (index(keywrd, "PM6-D3 ") /= 0) then
      write (iw, '(" *  PM6-D3     - CORRECT DISPERSION USING GRIMME''s D3 METHOD")')
    elseif (index(keywrd, "PM6-D3H4") /= 0) then
      write (iw, '(" *  PM6-D3H4   - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING THE D3H4 METHOD")')
    elseif (index(keywrd, "PM6-D3(H4)") /= 0) then
      write (iw, '(" *  PM6-D3(H4) - CORRECT DISPERSION USING THE D3H4 METHOD")')
    elseif (index(keywrd, "PM6-DH2") /= 0) then
      write (iw, '(" *  PM6-DH2    - CORRECT DISPERSION AND HYDROGEN BOND TERMS USING PM6-DH2")')
    else
      write (iw, '(" *  PM6-D      - CORRECT DISPERSION TERMS USING PM6-D")')
    end if
  else if (myword(allkey, " PM6-H"))  then
    write (iw, '(" *  PM6-H      - CORRECT HYDROGEN BOND TERMS USING PM6-H")')
  end if
  if (myword(allkey, " JCTC")) i = 6
  if (myword(allkey, " AIGIN"))  write (iw, '(" *  AIGIN      - GEOMETRY MUST BE IN GAUSSIAN FORMAT")')
  if (myword(allkey, " ESR"))    write (iw, '(" *  ESR        - RHF SPIN DENSITY CALCULATION REQUESTED")')
  if (myword(allkey, " NOMM"))   write (iw, '(" *  NOMM       - DO NOT MAKE MM CORRECTION TO CONH BARRIER")')
  if (myword(allkey, " MMOK"))   write (iw, '(" *  MMOK       - APPLY MM CORRECTION TO CONH BARRIER")')
  if (myword(allkey, " CIS "))   write (iw, '(" *  CIS        - C.I. USES 1 ELECTRON EXCITATIONS ONLY")')
  if (myword(allkey, " CISD "))  write (iw, '(" *  CISD       - C.I. USES 1 AND 2 ELECTRON EXCITATIONS")')
  if (myword(allkey, " CISDT ")) write (iw, '(" *  CISDT      - C.I. USES 1, 2 AND 3 ELECTRON EXCITATIONS")')
  if (myword(allkey, " SING"))   write (iw, '(" *  SINGLET    - SPIN STATE DEFINED AS A SINGLET")')
  if (myword(allkey, " DOUB"))   write (iw, '(" *  DOUBLET    - SPIN STATE DEFINED AS A DOUBLET")')
  if (myword(allkey, " TRIP"))   write (iw, '(" *  TRIPLET    - SPIN STATE DEFINED AS A TRIPLET")')
  if (myword(allkey, " QUAR"))   write (iw, '(" *  QUARTET    - SPIN STATE DEFINED AS A QUARTET")')
  if (myword(allkey, " QUIN"))   write (iw, '(" *  QUINTET    - SPIN STATE DEFINED AS A QUINTET")')
  if (myword(allkey, " SEXT"))   write (iw, '(" *  SEXTET     - SPIN STATE DEFINED AS A SEXTET")')
  if (myword(allkey, " SEPT"))   write (iw, '(" *  SEPTET     - SPIN STATE DEFINED AS A SEPTET")')
  if (myword(allkey, " OCTE"))   write (iw, '(" *  OCTET      - SPIN STATE DEFINED AS A OCTET")')
  if (myword(allkey, " NONE"))   write (iw, '(" *  NONET      - SPIN STATE DEFINED AS A NONET")')
  if (myword(allkey, " COSCCH")) write (iw, '(" *  COSCCH     - ADD IN COSMO CHARGE CORRECTIONS")')
  if (myword(allkey, " FIELD"))  write (iw, '(" *  FIELD      - AN EXTERNAL ELECTRIC FIELD IS TO BE USED")')
  if (myword(allkey, " NOREOR")) write (iw, '(" *  NOREOR     - DO NOT ALLOW THE SYSTEM TO BE REORIENTATED")')
  if (myword(allkey, " INVERT")) write (iw, '(" *  INVERT     - REVERSE ALL OPTIMIZATION FLAGS")')
  if (myword(allkey, " ESP "))   write (iw, '(" *  ESP        - ELECTROSTATIC POTENTIAL CALCULATION")')
  if (myword(allkey, " NSURF"))  write (iw, '(" *  NSURF      - NUMBER OF LAYERS")')
  if (myword(allkey, " NOGPU"))  write (iw, '(" *  NOGPU      - DO NOT USE GPU ACCELERATION")')
  if (myword(allkey, " SCALE"))  write (iw, '(" *  SCALE      - SCALING FACTOR FOR VAN DER WAALS DISTANCE")')
  if (myword(allkey, " SCINCR")) write (iw, '(" *  SCINCR     - INCREMENT BETWEEN LAYERS")')
  if (myword(allkey, " SLOPE"))  write (iw, '(" *  SLOPE      - SLOPE   - USED TO SCALE MNDO ESP CHARGES")')
  if (myword(allkey, " PDBOUT")) write (iw, '(" *  PDBOUT     - PRINT GEOMETRY IN PDB FORMAT")')
  if (myword(allkey, " NOTER"))  write (iw, '(" *  NOTER      - DO NOT PUT ""TER""S IN PDB FILE")')
  if (myword(allkey, " CHARGES"))write (iw, '(" *  CHARGES    - IDENTIFY AND PRINT CHARGED ATOMS")')
  if ( .not. l_add_H) l_add_H = (index(keywrd, " ADD-H") /= 0)
  if (myword(allkey, " NORES"))  write (iw, '(" *  NORESEQ    - DO NOT RESEQUENCE Z-MATRIX INTO NORMAL PDB FORMAT")')
  if (myword(allkey, " MACRO"))  write (iw, '(" *  MACRO      - MODIFY METHODS TO SMOOTH NDDO - POINT-CHARGE TRANSITION")')
  if (myword(allkey, " XENO"))   then
                                 write (iw, '(" *  XENO       - RENAME PROTEIN FRAGMENTS AND SMALL MOLECULES")')
    if ( .not. (Index (keywrd, " RESI") /= 0 )) then
      call mopend("XENO requires keyword RESIDUES to also be used")
      return
    end if
  end if
  if (myword(allkey, " RESEQ"))  then
                                 write (iw, '(" *  RESEQ      - RESEQUENCE Z-MATRIX INTO NORMAL PDB FORMAT")')
    if (index(keywrd, "RESID") == 0 .and. maxtxt /= txtmax) then
      call mopend("RESEQ only works when the atom labels are in PDB format")
      write(iw,'(a,/)')"(Before using RESEQ, run a job using keyword RESIDUES to add PDB atom labels.)"
      return
    end if
  end if
  i = index(keywrd," START_RES") + 1
  if (i > 1) then
    j = index(keywrd(i:),")")
    if (j == 0) then
      write (iw, '(" ***START_RES must be followed by open and close parentheses, e.g., START_RES=(1 10)")')
      call web_message(iw,"start_res.html")
      call mopend("Keyword START_RES not properly defined")
      return
    end if
    j = j + i
    allkey(i:j) = " "
    write (iw, '(" *  START_RES  - STARTING RESIDUE NUMBERS DEFINED")')
    write (iw, '(" *  Keyword:     ",a)')keywrd(i:j)
    if(myword(allkey, "START_RES")) then
      call mopend("Only one keyword START_RES allowed")
      return
    end if
  end if

  if (quoted('GEO_DAT=')  > 0) then
    i = index(keywrd_quoted," GEO_DAT")
    j = index(keywrd_quoted(i + 10:),'"') + i + 9 
    write (iw, '(" *  GEO_DAT    - DATA SET GEOMETRY IS IN FILE """,a,"""")') &
      keywrd_quoted(i + 10:j - 1)
  end if
  if (quoted('GEO_REF=')  > 0) then
    i = index(keywrd_quoted," GEO_REF")
    j = index(keywrd_quoted(i + 10:),'"') + i + 9 
    k = index(keywrd_quoted(i:j), "SELF")
    if (index(keywrd_quoted(i:j), "SELF") /= 0) then
      line = " *  GEO_REF=""SELF""    - USE MOPAC DATA SET """//keywrd_quoted(i + 10:j - 5) &
        //trim(job_fn)//""" AS REFERENCE GEOMETRY"
       write (iw, '(a)')trim(line)
    else
      write (iw, '(" *  GEO_REF    - REFERENCE GEOMETRY IS IN FILE """,a,"""")')keywrd_quoted(i + 10:j - 1)
    end if
    if (keywrd_quoted(j + 1:j + 1) == " ") then
      write(iw,'(a)') " *               (NO BIAS TOWARDS THE REFERENCE GEOMETRY WILL BE APPLIED)"
    else
      sum = reada(keywrd_quoted(j:), 1)
      if (sum < 1.d-15) then
        write(iw,'(a)') " *               (NO BIAS TOWARDS THE REFERENCE GEOMETRY WILL BE APPLIED)"
      else
        i = nint(log10(sum))
        if (abs(i) > 5) then
          write(iw,'(a, g11.4, a)') " *               (A BIAS OF",sum, &
          " KCAL/MOL/ANGSTROM^2 TOWARDS THE REFERENCE GEOMETRY WILL BE APPLIED)"
        else         
          num = char(ichar("4") + abs(i))
          if (i < 0) then
            num1 = char(ichar("1") - i)
          else
            num1 = "1"
          end if 
          write(iw,'(a, f'//num//'.'//num1//', a)') " *               (A BIAS OF",sum, &
          " KCAL/MOL/ANGSTROM^2 TOWARDS THE REFERENCE GEOMETRY WILL BE APPLIED)"
        end if
      end if
    end if
  end if
  i = index(allkey," CHAIN") + 1
  if (i > 1) then
    j = index(keywrd(i:),")") + i
    if (j == i) then
      write (iw, '(" ***CHAINS must be followed by open and close parentheses, e.g., CHAINS=(ABCD)")')
      call web_message(iw,"chains.html")
      call mopend("Keyword CHAINS not properly defined")
      return
    end if
    if (.not. myword(allkey, " CHAIN")) return
                                 write (iw, '(" *  CHAINS     - PDB CHAIN LETTERS EXPLICITLY DEFINED")')
                                 write (iw, '(" *  Keyword:     ",a)')keywrd(i:j)
  end if
  if (myword(allkey, " NEWPDB")) write (iw, '(" *  NEWPDB     - CONVERT PDB ATOM FORMAT INTO THE MODERN PDB VERSION, VERSION-3")')
  if (myword(allkey, " GEOCHK")) write (iw, '(" *  GEOCHK     - PRINT WORKING IN SUBROUTINE GEOCHK")')
  if (myword(allkey, " LEWIS"))  write (iw, '(" *  LEWIS      - PRINT OUT LEWIS STRUCTURE, THEN STOP")')
  if (myword(allkey, " SETPI"))  write (iw, '(" *  SETPI      - SOME OR ALL PI BONDS EXPLICITLY SET BY USER")')
  if (myword(allkey, " PDB "))   write (iw, '(" *  PDB        - INPUT GEOMETRY IS IN PDB FORMAT")')
  if (myword(allkey, " PDB("))   write (iw, '(" *  PDB(txt)   - SYMBOLS IN PDB FILE ARE DEFINED BY USER")')
  if (myword(allkey, " MOZ"))    write (iw, '(" *  MOZYME     - USE LOCALIZED M.O.s IN SOLVING THE SCF EQUATIONS")')
  if (myword(allkey, " RAPID"))  write (iw, '(" *  RAPID      - IN MOZYME SCF, USE ATOMS BEING OPTIMIZED ONLY")')
  if (myword(allkey, " PINOUT")) write (iw, '(" *  PINOUT     - WRITE OUT THE LMOs WHEN READING OR WRITING A ''.DEN'' FILE")')
  if (myword(allkey, " REORTH")) write (iw, '(" *  REORTH     - RE-ORTHOGONALIZE LMOs EACH 10 SCF CALCULATIONS")')
  if (myword(allkey, " RE-LOCAL")) write(iw,'(" *  RE-LOCAL   - RE-LOCALIZE LMOs AT END OF CALCULATION")')
  if (myword(allkey, " SWAP"))   write(iw,  '(" *  SWAP       - THIS KEYWORD IS NOW REDUNDANT. SEE KEYWORD NOSWAP")')
  if (myword(allkey, " NOSWAP")) write(iw,  '(" *  NOSWAP     - DO NOT SWAP ATOMS EVEN IF IT WILL IMPROVE OVERLAP IN GEO_REF")')
  if (myword(allkey, " A0 "))    write (iw, '(" *  A0         - INPUT GEOMETRY IS IN ATOMIC UNITS (A0)")')
  if (myword(allkey, " ANG"))    write (iw, '(" *  ANGSTROMS- INPUT GEOMETRY IS IN ANGSTROMS")')
  if (myword(allkey, " ADD-H"))  write (iw, '(" *  ADD-H      - ADD HYDROGEN ATOMS TO SATISFY VALENCE")')
  if (myword(allkey, " CONNOL")) write (iw, '(" *  CONNOLLY   - USE CONNOLLY SURFACE")')
  if (myword(allkey, " ESPRST")) write (iw, '(" *  ESPRST     - RESTART OF ELECTRIC POTENTIAL CALCULATION")')
  if (myword(allkey, " UHF"))    write (iw, '(" *  UHF        - UNRESTRICTED HARTREE-FOCK CALCULATION")')
  if (myword(allkey, " RHF"))    write (iw, '(" *  RHF        - RESTRICTED HARTREE-FOCK CALCULATION")')
  if (myword(allkey, " STATIC")) write (iw, '(" *  STATIC     - CALCULATE STATIC FIELD POLARIZABILITIES")')
  i = quoted(' SETUP=')
  if (i /= 0) then
    if (len_trim(line) == 0) then
      allkey(i:i + 7) = " " 
      line = "SETUP or SETUP.TXT"
    end if
    i = len_trim(line)
    if (i < 26) then
      write (iw, '(" *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE """, a, """")') trim(line)
    else
      write (iw, '(" *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE:")')
      write (iw, '("                 """, a, """")') trim(line)
    end if
  else
    i = index(allkey," SETUP")
    if (i /= 0) then
      j = index(allkey(i + 6:), " ") + i + 4
      write (iw, '(" *  SETUP      - EXTRA KEYWORDS TO BE READ FROM FILE """,a,"""")')allkey(i + 1:j)
      allkey(i + 1:j) = " "
    end if
  end if
  if (myword(allkey, " MAX"))    write (iw, '(" *  MAX        - GRID SIZE 23*23 ")')
  if (myword(allkey, " COSWRT")) write (iw, '(" *  COSWRT")')
  if (myword(allkey, " OLDCAV")) write (iw, '(" *  OLDCAV")')
  if (myword(allkey, " SYM ") .or. myword(allkey, " SYMM")) &
    write (iw, '(" *  SYMMETRY   - SYMMETRY CONDITIONS TO BE IMPOSED")')
  if (myword(allkey, " POLAR")) then
    write (iw, '(" *  POLAR      - CALCULATE FIRST, SECOND AND THIRD-ORDER POLARIZABILITIES")')
    if (id /= 0) then
      call mopend("KEYWORD ""POLAR"" CANNOT BE USED WITH PERIODIC BOUNDARY CONDITIONS")
      write(iw,'(2x,a)')"(TV in the geometry implies PBC, but PBC cannot allow an electric field gradient.) "
    end if
  end if
  if (myword(allkey, " RESI")) then
    write (iw, '(" *  RESIDUES   - DETERMINE THE SEQUENCE OF RESIDUES")')
    if (index(keywrd," RESIDUES0") /= 0) then
      if (index(keywrd," ADD-H") + index(keywrd," RESEQ") + index(keywrd, " SITE=") + index(keywrd, " SITE(") /= 0) then
        call mopend("Keyword RESIDUES0 cannot be used when a keyword that changes the geometry is present")
        if (index(keywrd," ADD-H") /= 0) write(iw,'(10x,a)')"Keyword ADD-H changes the geometry."
        if (index(keywrd," RESEQ") /= 0) write(iw,'(10x,a)')"Keyword RESEQ changes the geometry."
        if (index(keywrd," SITE=(") /= 0) write(iw,'(10x,a)')"Keyword SITE changes the geometry."
        call web_message(iw,"Residues.html")
      end if
    end if
    if(myword(allkey, " NOZERO")) &
      write (iw, '(" *  NOZERO     - THE NUMBER ZERO IN THE RESIDUE SEQUENCE IS DELETED")')
    if(myword(allkey, " ZERO")) &
      write (iw, '(" *  ZERO       - A RESIDUE IS ALLOWED TO HAVE THE NUMBER ZERO")')
  end if
  if (myword(allkey, " NORES"))  &
    write (iw, '(" *  NORES      - THIS IS THE DEFAULT.  USE ""RESIDUES"" IF RESIDUES ARE TO BE CALCULATED")')
  if (myword(allkey, " SITE"))   then
    write (iw, '(" *  SITE       - SET IONIZATION LEVELS OF IONIZABLE RESIDUES ")')
    k = 0
    do i = 1, natoms
      if (na(i) /= 0) then
        k = k + 1
        if (k == 1) then
          call mopend("KEYWORD ""SITE"" REQUIRES ALL ATOMS TO BE IN CARTESIAN COORDINATES")
          write(iw,'(/10x,a)') "Atoms in internal coordinates:"
        end if
        write(iw,'(12x, a)') txtatm(i)
      end if
    end do
    if (k > 0) write(iw,'(/10x,a)') "(Keyword ""XYZ"" converts internal coordinates to Cartesian coordinates.)"
    if (k > 5) write(iw,'(/10x,a)') "Remaining errors suppressed"
    if (index(keywrd, " RESIDU") /= 0) then
      call mopend("THE PRESENCE OF KEYWORD ""RESIDUES"" PREVENTS KEYWORD ""SITE"" FROM WORKING CORRECTLY")
      write(iw,*)
      return
    end if
    i = index(keywrd, " SITE=(")
    do j = i, len_trim(keywrd)
      if (keywrd(j:j + 1) == ") ") exit
    end do
    allkey(i:j) = " "
  end if
  if (myword(allkey, " NOOPT"))   then
    line = trim(keywrd)
    do
      i = index(line, " NOOPT")
      if (i == 0) exit
      j = index(line(i + 5:)," ") + i + 3
      if (j   - i == 8) then
        line(j:j) = char( ichar(line(j:j)) + ichar("a")   - ichar("A"))
        write (iw, '(" *  NOOPT-",a2,"   - DO NOT OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a2)')line(j - 1:j), line(j   - 1:j)
      else if (j   - i == 7) then
        write (iw, '(" *  NOOPT-",a1,"    - DO NOT OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a1)')line(j:j), line(j:j)
      else if (j   - i == 5) then
        write (iw, '(" *  NOOPT      - DO NOT OPTIMIZE ANY COORDINATES")')
        if (index(keywrd, " OPT ") /= 0) then
          write (iw, '(" *  OPT        - OPTIMIZE COORDINATES OF ALL ATOMS")')
          write (iw, '(" *             - USE ""NOOPT"" OR ""OPT"" BUT NOT BOTH")')
          call mopend("OPT AND NOOPT CANNOT BOTH BE PRESENT")
          return
        end if
      else
        write (iw, '(" ***NOOPT      - IMPROPER USE OF NOOPT KEYWORD, KEYWORD USED = ''", a,"''")')line(i + 1:j)
        call mopend("IMPROPER USE OF OPT KEYWORD")
        return
      end if
      line(i:j) = " "
    end do
  end if
  line = trim(allkey)
  i = 0
  if (myword(line, " OPT-")) i = i + 1
  if (myword(line, " OPT ")) i = i + 1
  if (myword(line, " OPT(")) i = i + 1
  if (myword(line, " OPT=")) i = i + 1
  if (i > 1) then
    call mopend("MORE THAN ONE TYPE OF ""OPT"" REQUESTED. THIS IS NOT ALLOWED")
    write(iw,'(/10x,a,//,a,/)')"KEYWORDS SUPPLIED:", trim(keywrd)
    return
  end if
  line = trim(allkey)
  if (myword(allkey, " OPT-") .or. myword(allkey, " OPT ") .or. myword(allkey, " OPT(") &
    .or. myword(allkey, " OPT="))   then
    i = index(line, " OPT=")
    if (i /= 0) line(i + 4:) = line(i + 5:)
    l = 0
    do
      i = index(line, " OPT-") + index(line, " OPT ") + index(line, " OPT(") + index(line, " OPT=(")
      if (i == 0) exit
      j = index(line(i + 4:)," ") + i + 2
      if (index(line, " OPT-") /= 0) then
        if (j - i == 6) then
          line(j:j) = char( ichar(line(j:j)) + ichar("a")   - ichar("A"))
          write (iw, '(" *  OPT-",a2,"     - OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a2)')line(j - 1:j), line(j - 1:j)
        else if (j - i == 5) then
          write (iw, '(" *  OPT-",a1,"      - OPTIMIZE COORDINATES OF ATOMS OF TYPE ",a1)')line(j:j), line(j:j)
        else
          l = 1
        end if
      else if (index(line(i:), " OPT(") /= 0) then
        j = index(line(i + 1:), ") ")
        if (j /= 0) then
          i = i + 5
          k = index(line(i:j + i - 5), "=")
          if (k /= 0) then
            write(iw,'(" *  OPT        - OPTIMIZE COORDINATES OF ALL ATOMS WITHIN ")')
            j = i + j
            do
              sum = reada(line, k + i)
              write(iw,'(" *",14x,f5.2, " ANGSTROMS OF ATOMS WITH PDB LABEL ",a)') sum, line(i: k + i -2)
              l = index(line(i:j), ",")
              if (l /= 0) then
                i = i + l
                k = index(line(i:j), "=")
              else
                exit
              end if
            end do
          else
            write(iw,'(" *  OPT        - OPTIMIZE COORDINATES OF ALL ATOMS IN RESIDUES:",/," *",15x,a, &
            & a)') line(k + i:j + i - 6), line(i: k + i -2)
          end if
        else
          l = 1
        end if
        exit
      else if (index(line(i:), " OPT ") /= 0) then
        write (iw, '(" *  OPT        - OPTIMIZE COORDINATES OF ALL ATOMS")')
      else
        l = 1
      end if
      line(i:j) = " "
    end do
    if (l == 1) then
      i = index(line, " OPT")
      j = index(line(i + 1:), " ")
      write (iw, '(" ***OPT        - IMPROPER USE OF OPT KEYWORD, KEYWORD USED = ''", a,"''")')line(i + 1:j)
      call mopend("IMPROPER USE OF OPT KEYWORD")
      return
    end if
  end if
!
!  Keywords involving equals sign
!
!**********************************************************************
!
    !  Tension or pressure in a solid-state calculation
    !
    pressure = 0.d0
    if (myword(allkey, " P=")) then
      pressure = reada (keywrd, Index (keywrd, " P="))
      if (id == 1) then
        write (iw, '(" *  P=         - TENSION IN POLYMER=", g13.6, " NEWTONS PER MOLE")') pressure
        pressure = pressure * 10.d0 ** (-13) / 4.184d0
      else if (id == 2) then
      else if (id == 3) then
        i = Index (keywrd, " P=")
        j = Index (keywrd(i + 3:), " ") + i + 3
        if (index(keywrd(i + 3: j),"GP") /= 0) then
         write (iw, '(" *  P=         - PRESSURE ON SOLID=", f7.3, &
             & " GIGAPASCALS")') pressure
          pressure = pressure*1.d9
        else
          write (iw, '(" *  P=         - PRESSURE ON SOLID=", g13.6, &
             & " NEWTONS PER SQUARE METER")') pressure
        end if
!
!  Multiply by N, Avogadro's Number, to convert from J/M**3 per molecule
!  to J/M**3/mol.
!  Divide by 4184 to convert from J/M**3/mol to Kcal/M**3/mol
!  Divide by 10**30 to convert from Kcal/M**3/mol to Kcal/Angstrom**3/mol
!
        pressure = (fpcref(1,10)*pressure) / (4184.d0*10.d0**30)
      else
        write (iw, *) " Keyword 'P=n.nn' is not allowed here"
        call mopend("Keyword 'P=n.nn' is not allowed here")
        return
      end if
    end if
!
!                       ESP grid
!
    i = -100
    if (myword(allkey, " CUBE")) then
      write (iw,'(" *  CUBE       - USE GAUSSIAN CUBE FILE FOR ELECTROSTATICS")')
    else if (myword(allkey, " ESPGRID")) then
      i = Nint (reada (keywrd, Index (keywrd, " ESPGRID")))
      write (iw,'(" *  ESPGRID    - GENERATE ELECTROSTATICS CUBE FILE  =", i5," POINTS ON SIDE")') i
    end if
    if (i /= -100 .and. i < 2 .or. i > 100) then
      if (i < 2) write (iw,'(//," *   Value of number of points is too small",//)')
      if (i > 100) write (iw,'(//," *   Value of number of points is too large",//)')
      call mopend("Number of points in electrostatics is outside limits")
      return
    end if

!
!
!               Electronic quantities (mainly C.I.)
!
!                       CHARGE
!
  if (myword(allkey, " CHARGE=")) then
    i = Nint (reada (keywrd, Index (keywrd, " CHARGE=")))
    if (i == 0) then
      write (iw,'(3(" *", /), " *", 15 x, "  CHARGE ON SYSTEM = 0", 3 (/, " *"))')
    else
      num = char(ichar("3") + max(int(log10(abs(i) + 0.05)),0))
      write (iw,'(3(" *", /), " *", 15 x, "  CHARGE ON SYSTEM =", SP,i'//num//', 3 (/, " *"))') i
    end if
    if (id /= 0 .and. i /= 0) then
      write(iw,"(/10x,a)")"INFINITE SYSTEMS MUST HAVE A ZERO CHARGE ON THE UNIT CELL"
      call mopend("Unit cell has a charge. Correct fault and re-submit ")
    end if
  end if
  if (Index(allkey, " C.I.D.") /= 0) then
    j = Index (keywrd, " C.I.D.=(")
    if (j /= 0) then
      j = Index (keywrd(j:j+10), ",") + j   - 1
      nmos = Nint (reada (keywrd, Index (keywrd, "C.I.D.=(")+7))
      num = char(ichar("2") + int(log10(nmos*1.0)))
      l = Int (reada (keywrd, j))
      num1 = char(ichar("2") + int(log10(l*1.0)))
      write (iw,'(" * C.I.D.=(N,M)-", i'//num1// &
        ', " DOUBLY FILLED LEVELS USED IN A CI DOUBLES INVOLVING", i'//num//', " M.O.''S")') l, nmos
     if (.not. myword(allkey, " C.I.D.=(")) return ! Impossible option used to delete keyword
    else if (myword(allkey, " C.I.D.=")) then
      nmos = Int (reada (keywrd, Index (keywrd, "C.I.D.")+7))
      write (iw,' (" *  C.I.D.=N   -", i2, " M.O.S TO BE USED IN CI DOUBLES")') nmos
    else
      write (iw, "(' *',/,a)") " *      C.I.D. keyword must be of form 'C.I.D.=n' or 'C.I.D.=(n1,n2)' (See manual)"
      call mopend("C.I.D. keyword must be of form 'C.I.D.=n' or 'C.I.D.=(n1,n2)'")
      return
    end if
  end if
!
!                       C.I.=(n1,n2)
!
  if (Index(allkey, " C.I.") /= 0) then
    j = Index (keywrd, " C.I.=(")
    if (j /= 0) then
      j = Index (keywrd(j:j+10), ",") + j   - 1
      nmos = Nint (reada (keywrd, Index (keywrd, "C.I.=(")+5))
      num = char(ichar("2") + int(log10(nmos*1.0)))
      l = Int (reada (keywrd, j))
      num1 = char(ichar("2") + int(log10(l*1.0)))
      write (iw,'(" *  C.I.=(N,M) -", i'//num1//', " DOUBLY FILLED LEVELS USED IN A C.I. INVOLVING", i' &
        //num//', " M.O.''S")') l, nmos
     if (.not. myword(allkey, " C.I.=(")) return ! Impossible option used to delete keyword
    else if (myword(allkey, " C.I.=")) then
      nmos = Int (reada (keywrd, Index (keywrd, "C.I.")+5))
      i = Index (keywrd, "C.I.=")
      j = Index(keywrd(i + 1:), " ") + i
      line = " *  "//keywrd(i:j)
      num = char(ichar("2") + int(log10(nmos*1.0)))
      write (iw, "(a, i"//num//", a)")line(1:15)//"-", nmos, " M.O.S TO BE USED IN C.I."
    else
      write (iw, "(' *',/,a)") " *      C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)' (See manual)"
      call mopend("C.I. keyword must be of form 'C.I.=n' or 'C.I.=(n1,n2)'")
      call web_message(iw,"CI=nm.html")
      call web_message(iw,"CI=n.html")
      return
    end if
  end if
  if (Index(allkey, " C.A.S.") /= 0) then
    j = Index (keywrd, " C.A.S.=(")
    if (j /= 0) then
      j = Index (keywrd(j:j+10), ",") + j   - 1
      nmos = Nint (reada (keywrd, Index (keywrd, "C.A.S.=(")+7))
      write (iw,'(" * C.A.S.=(N,M)-", i3, " DOUBLY FILLED LEVELS USED IN A",/, &
     & " *             CAS REF DET GENERATION INVOLVING ", i2, " M.O.''S")')  &
     Int (reada (keywrd, j)), nmos
     if (.not. myword(allkey, " C.A.S.=(")) return ! Impossible option used to delete keyword
    else if (myword(allkey, " C.A.S.=")) then
      nmos = Int (reada (keywrd, Index (keywrd, "C.A.S.")+7))
      write (iw,' (" *  C.A.S.=N   -", i2, " M.O.S TO BE USED IN CAS REF DET GENERATION")') nmos
    else
      write (iw, "(' *',/,a)") " *      C.A.S. keyword must be of form 'C.A.S.=n' or 'C.A.S.=(n1,n2)' (See manual)"
      call mopend("C.A.S. keyword must be of form 'C.A.S.=n' or 'C.A.S.=(n1,n2)'")
      return
    end if
  end if
!
!                       MECI Microstates read in
!
  if (myword(allkey, " MICROS")) &
    write (iw,'(" *  MICROS=N -", i4, " MICROSTATES TO BE SUPPLIED FOR C.I.")') &
    Int (reada (keywrd, Index (keywrd, " MICROS")))
!
!                       Fractionally occupied degenerate Open shell
!
  if (myword(allkey, " OPEN")) then
    i = Index (keywrd, " OPEN")
    j = Index (keywrd(i:i+10), ",") + i   - 1
    ilevel = Nint (reada (keywrd, j))
    ielec = Nint (reada (keywrd, Index (keywrd, " OPEN")+6))
    write (iw,'(" *  OPEN(M,N)  - THERE ARE", i2, " ELECTRONS IN", i2, " LEVELS")') &
     ielec, ilevel
  end if
!
!                       Magnetic component of spin
!
  if (myword(allkey, " MS=")) &
    write (iw,'(" *  MS=        - IN MECI, MAGNETIC COMPONENT OF SPIN =", f5.1)') &
    reada (keywrd, Index (keywrd, " MS="))
!
!   Select root of C.I. matrix
!
  if (myword(allkey, " ROOT")) then
    i = Index (keywrd, " ROOT") + 6
    j = index(keywrd(i:), " ") + i - 2
    line = keywrd(i:j)
    do k = 3,5
      if (line(k:k) >="A" .and. line(k:k) <= "Z") &
      & line(k:k) = Char(Ichar(line(k:k)) + Ichar('a') - Ichar('A'))
    end do
    write (iw,'(" *  ROOT       - IN A C.I. CALCULATION, ROOT """, a, """ TO BE OPTIMIZED.")') trim(line)
    if (keywrd(i:i) > "9" .or. keywrd(i:i) < "0") then
      call mopend("THE VALUE OF ""ROOT"" MUST START WITH AN INTEGER")
    end if
  end if
!**********************************************************************
!
!   SOLVATION KEYWORDS
!

  if (myword(allkey, " EPS=") .or. (index(keywrd," PKA") .ne. 0)) then
    if (Index(keywrd," EPS=CRSDEF") /= 0) then
      fepsi = 1.d0
    else
      if (index(keywrd," PKA") .ne. 0) then
        epsi = 78.4d0
        line = "EPS=78.4 "
      else
        i = Index (keywrd, "EPS=")
        j = Index(keywrd(i + 1:), " ") + i
        line = " *  "//keywrd(i:j)
        epsi = reada (keywrd, i)
      end if
      write (iw, "(a, a)")line(1:15),"- USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL"
      fepsi = (epsi-1.d0) / (epsi+0.5d0)
    end if
    if (epsi < -100.d0) then
      write(iw,"(a)") " *  "
      write(iw,"(a)") " *  The lower bound of a dielectric is 1.00 = vacuum"
      write(iw,"(2a)")" *  Keywords: ",keywrd(:len_trim(keywrd))
      write(iw,"(a)") " *  "
      call mopend("Nonsense value of dielectric constant supplied")
      return
    end if
  else if (myword(allkey, " COSMO")) then
    epsi = 999.d0
    write (iw, "(a, f6.2, a)")" *  EPS=", epsi," - USE ANDREAS KLAMT'S COSMO IMPLICIT SOLVATION MODEL"
    fepsi = (epsi-1.d0) / (epsi+0.5d0)
  else if (numat < 100) then  !  Limit surface area calculation to 100 atoms   - saves time
    epsi = 78.4d0   !  For surface area and volume only
    fepsi = (epsi-1.d0) / (epsi+0.5d0)
  end if
  if (myword(allkey, " DIPL"))  write (iw,'(" *  DIPL=", f7.3)') reada (keywrd, Index (keywrd, " DIPL"))
  if (myword(allkey, " RSOLV")) write (iw,'(" *  RSOLV=", f7.3)') reada (keywrd, Index (keywrd, " RSOLV"))
  if (myword(allkey, " DELSC")) write (iw,'(" *  DELSC=", f7.3)') reada (keywrd, Index (keywrd, " DELSC"))
  if (myword(allkey, " DISEX")) write (iw,'(" *  DISEX=", f7.3)') reada (keywrd, Index (keywrd, " DISEX"))
  i = Index(keywrd, "N**2")
  if (i /= 0) then
    j = Index(keywrd(i + 1:), " ") + i
    line = " *  "//keywrd(i:j)
    write (iw, "(a, a)")line(1:14), " - USED IN COSMO-CI"
    if (.not. myword(allkey, "N**2")) return
  end if
  if (myword(allkey, " NSPA")) then
    i = Nint (reada (keywrd, Index (keywrd, " NSPA")))
    write (iw,'(" *  NSPA=", i5)') i
    nspa = i
  else
    nspa = 42
  end if
  if (myword(allkey, " ROTX"))  write (iw,'(" *  ROTX=", i5)') Nint (reada (keywrd, Index (keywrd, " ROTX")))
!**********************************************************************
!
!                       Geometric quantities
!
!                       Grids, Reaction paths, etc.
!
 i = 0
 if (myword(allkey, " STEP1"))  then
   write (iw,'(" *  STEP1      - FIRST STEP-SIZE IN GRID =", f7.2)') &
   reada (keywrd, Index (keywrd, "STEP1")+6)
   if (index(keywrd, " POINT1") == 0) then
    call l_control("POINT1=11", len_trim("POINT1=11"), 1)
   end if
   i = 2
 end if
 if (myword(allkey, " STEP2"))  then
   write (iw,'(" *  STEP2      - SECOND STEP-SIZE IN GRID =", f7.2)') &
   reada (keywrd, Index (keywrd, "STEP2")+6)
   if (index(keywrd, " POINT2") == 0) then
    call l_control("POINT2=11", len_trim("POINT2=11"), 1)
   end if
 end if
 if (myword(allkey, " STEP="))  then
   write (iw,'(" *  STEP       - STEP-SIZE IN PATH =", f8.3)') &
   reada (keywrd, Index (keywrd, " STEP") + 4)
   if (index(keywrd, " POINT") == 0) call mopend ("KEYWORD POINT MISSING")
   if (index(keywrd, " 1SCF") /= 0) call mopend ("KEYWORD ""1SCF"" CANNOT BE USED WITH KEYWORD ""STEP""")
   i = 1
 end if
 if (i > 0) then
!
! Check that the appropriate flags are set
!
   k = 0
   do j = 1, natoms
     do l = 1, 3
       if (lopt(l,j) < 0) k = k + 1
     end do
   end do
   if (i /= k) then
     if (index(keywrd," 0SCF") == 0) then
       write(iw,'(10x, a,i1)')"Number of optimization flags requested by keywords: ", i
       write(iw,'(10x, a,i1)')"Number of optimization flags found in data-set:     ", k
       call mopend("The number of optimization flags set to '-1' does not match the keywords used")
     end if
   end if
 end if
 if (myword(allkey, " MERS"))then
   j = index(keywrd," MERS")
   mers = 0
   k = 0
   i = Index (keywrd(j + 1:), " ") + j
   do l = 1, 3
     j = j + k
     if (l > 1 .and. k == 0) exit
     mers(l) = Nint (reada (keywrd(j:), 1))
     k = Index (keywrd(j:i), ",")
   end do
   i = mers(1)
   do l = 2, 3
     if (mers(l) == 0) exit
     i = i*mers(l)
   end do
   write (iw,'(" *  MERS=N     - NUMBER OF FUNDAMENTAL UNIT CELLS USED:", i3)') i
 else
   mers = 0
 end if
 if (myword(allkey, " POINT1")) then
   write (iw,'(" *  POINT1     - NUMBER OF ROWS IN GRID =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT1")+7))
   if (index(keywrd, " STEP1") == 0) then
     write (iw,'("*",/," *             - **** KEYWORD STEP1 MISSING ****",/,"*")')
     call mopend("KEYWORD STEP1 MISSING")
   end if
 end if
 if (myword(allkey, " POINT2")) then
   write (iw,'(" *  POINT2     - NUMBER OF COLUMNS IN GRID =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT2")+7))
   if (index(keywrd, " STEP2") == 0) then
     write (iw,'("*",/," *             - **** KEYWORD STEP2 MISSING ****",/,"*")')
     call mopend("KEYWORD STEP2 MISSING")
   end if
 end if
 if (myword(allkey, " POINT")) then
   write (iw,'(" *  POINT      - NUMBER OF POINTS IN PATH =", i3)') &
   Nint (reada (keywrd, Index (keywrd, "POINT")+6))
   if (index(keywrd, " STEP") == 0) then
     write (iw,'("*",/," *             - **** KEYWORD STEP MISSING ****",/,"*")')
     call mopend("KEYWORD STEP MISSING")
   end if
 end if
 if (myword(allkey, " KINE"))   write (iw,'(" *  KINETIC=   - ", f9.3, " KCAL KINETIC ENERGY ADDED TO DRC")') &
    reada (keywrd, Index (keywrd, " KINE"))
!
!                       Specific intrinsic reaction coordinate
!
 if (myword(allkey, " IRC="))  then
    write (iw,'(" *  IRC=N      - INTRINSIC REACTION COORDINATE", i3, " DEFINED")') &
    Nint (reada (keywrd, Index (keywrd, " IRC=")))
  else if (myword(allkey, " IRC")) then
    write (iw,'(" *  IRC        - INTRINSIC REACTION COORDINATE CALCULATION")')
  end if
!
!                       DRC (may have half-life)
!
  if (myword(allkey, " DRC=")) then
    write (iw,'(" *  DRC=       - HALF-LIFE FOR KINETIC ENERGY LOSS =", f9.2, " *10**(-15) SECONDS")') &
     reada (keywrd, Index (keywrd, " DRC="))
  else if (myword(allkey, " DRC")) then
    write (iw,'(" *  DRC        - DYNAMIC REACTION COORDINATE CALCULATION")')
  end if
   mozyme = l_temp
!
!                       External parameters read from file
!
  if (myword(allkey, " EXTERNAL")) write (iw, '(" *  EXTERNAL   - DEFAULT PARAMETERS RESET USING DATA IN INPUT FILE")')
  if (quoted('EXTERNAL=')  > 0) then
    i = index(keywrd_quoted," EXTERNAL=")
    i = i + 10
    line = get_a_name(keywrd_quoted(i:), len_trim(keywrd_quoted(i:)))
    write (iw, '(" *",/," *  EXTERNAL=n -  DEFAULT PARAMETERS RESET USING DATA IN FILES: ",/," *",17x, a)') '"'//trim(line)//'"'
    do
      j = index(keywrd_quoted(i:), ";")
      if (j /= 0) then
        i = i + j
        line = get_a_name(keywrd_quoted(i:), len_trim(keywrd_quoted(i:)))
        write (iw, '(" *", 10x, a)')'   and "'//trim(line)//'"'
      else
        exit
      end if
    end do
  end if
  return
end subroutine wrtcon
subroutine wrtout (allkey)
  use molkst_C, only : keywrd, mozyme, maxtxt, line, prt_coords, prt_gradients, prt_cart, prt_charges, prt_pops, &
    prt_topo, prt_force, prt_normal_coords, prt_orientation, prt_velocity, backslash, is_PARAM
  use chanel_C, only: iw0, iw, log, input_fn
  implicit none
  character (len=1000), intent (inout) :: allkey
  integer :: i, j
  logical, external :: myword
  double precision, external :: reada
  intrinsic Index, Nint
  if (myword(allkey, " PRTINT"))  write (iw,'(" *  PRTINT     - INTERATOMIC DISTANCES TO BE PRINTED")')
  if (myword(allkey, " PRTCHAR")) write (iw,'(" *  PRTCHARGE  - PRINT CHARGES IN ARC FILE AND PDB OUTPUT")')
  if (index(allkey, " AUX") > 0) then
    i = index(allkey, " AUX(") + 1
    if (i > 1) then
      j = index(allkey(i + 4:), ") ") + i + 3
      allkey(i:j) = " "
    else
      i = index(allkey, " AUX")
      allkey(i:i + 3) = " "
    end if
    write (iw,'(" *  AUX        - OUTPUT AUXILIARY INFORMATION")')
    if (myword(allkey, " AUX")) i = i + 1
  end if
  if (index(allkey, " OUTPUT") > 0) then
    i = index(allkey, " OUTPUT(") + 1
    if (i > 1) then
      j = index(allkey(i + 4:), ") ") + i + 3
      line = allkey(i + 6:j - 1)
    else
      line = " "
    end if
    prt_coords        = (index(line, "C") /= 0)
    prt_gradients     = (index(line, "G") /= 0)
    prt_cart          = (index(line, "X") /= 0)
    prt_charges       = (index(line, "Q") /= 0)
    prt_pops          = (index(line, "P") /= 0)
    prt_topo          = (index(line, "T") /= 0)
    prt_force         = (index(line, "F") /= 0)
    prt_normal_coords = (index(line, "N") /= 0)
    prt_orientation   = (index(line, "O") /= 0)
    prt_velocity      = (index(line, "V") /= 0)
    if (prt_gradients .or. i > 1) then
      write (iw,'(" *  OUTPUT     - REDUCE OUTPUT, BUT PRINT:")')
    else
      write (iw,'(" *  OUTPUT     - MINIMIZE OUTPUT")')
    end if
    if (prt_coords)        write (iw,'(" *           C - COORDINATES")')
    if (prt_gradients)     write (iw,'(" *           G - GRADIENTS")')
    if (prt_cart)          write (iw,'(" *           X - CARTESIAN COORDINATES")')
    if (prt_charges)       write (iw,'(" *           Q - FRACTIONAL ATOMIC CHARGES")')
    if (prt_pops)          write (iw,'(" *           P - ATOMIC ORBITAL POPULATIONS")')
    if (prt_topo)          write (iw,'(" *           T - TOPOGRAPHY (ATOM CONNECTIVITY)")')
    if (prt_force)         write (iw,'(" *           F - FORCE MATRIX")')
    if (prt_normal_coords) write (iw,'(" *           N - NORMAL COORDINATES")')
    if (prt_orientation)   write (iw,'(" *           O - ORIENTATION)")')
    if (prt_velocity)      write (iw,'(" *           V - VELOCITY)")')
    if (.not. myword(allkey, " OUTPUT")) return ! An impossible option
    if (index(keywrd, " FORCE ") /= 0 .and. .not. prt_force ) then
      if (index(keywrd, " GEO-OK") + index(keywrd, " ISOTOPE") == 0) then
        call mopend("KEYWORD ""FORCE"" PRESENT BUT KEYWORD ""OUTPUT"" SUPPRESSES RESULTS OF ""FORCE"" CALCULATION")
        write (iw,'(10x, "(Either add ""GEO-OK"" or ""ISOTOPE"", remove ""OUTPUT"", or use ""OUTPUT(F)"")")')
      end if
    end if
  else
    prt_coords        = .TRUE.
    prt_gradients     = .TRUE.
    prt_cart          = .TRUE.
    prt_charges       = .TRUE.
    prt_pops          = .TRUE.
    prt_topo          = .TRUE.
    prt_force         = .TRUE.
    prt_normal_coords = .TRUE.
    prt_orientation   = .TRUE.
    prt_velocity      = .TRUE.
  end if
  if (myword(allkey, " MOPAC"))   write (iw,'(" *  MOPAC      - USE OLD MOPAC CONVENTION FOR FIRST THREE ATOMS")')
  if (myword(allkey, " NOINT"))   write (iw,'(" *  NOINTER    - THIS KEYWORD HAS BEEN REPLACED BY PRTINT")')
  if (myword(allkey, " NOCOM"))   write (iw,'(" *  NOCOMMENTS - SUPPRESS PDB COMMENTS")')
  if (myword(allkey, " ISOTOPE")) write (iw,'(" *  ISOTOPE    - SAVE THE FORCE MATRIX FOR USE LATER")')
  if (myword(allkey, " DENOUT"))  write (iw,'(" *  DENOUT     - SAVE DENSITY MATRIX FOR USE LATER")')
  if (myword(allkey, " OLDENS"))  write (iw,'(" *  OLDENS     - INITIAL DENSITY MATRIX READ OFF DISK")')
  if (myword(allkey, " CIOSCI"))  write (iw,'(" *  CIOSCI     - PRINT WORKING IN SUBROUTINE CIOSCI")')
  if (myword(allkey, " ENPART"))  write (iw,'(" *  ENPART     - ENERGY TO BE PARTITIONED INTO COMPONENTS")')
  if (myword(allkey, " NOXYZ"))   write (iw,'(" *  NOXYZ      - CARTESIAN COORDINATES NOT TO BE PRINTED")')
  if (myword(allkey, " PRTXYZ"))  write (iw,'(" *  PRTXYZ     - PRINT CARTESIAN COORDINATES")')
  if (myword(allkey, " NOTXT"))   write (iw,'(" *  NOTXT      - DO NOT PRINT TEXT ASSOCIATED WITH AN ATOM")')
  if (myword(allkey, " GRAD"))    write (iw,'(" *  GRADIENTS  - ALL GRADIENTS TO BE PRINTED")')
  if (myword(allkey, " MINI"))    write (iw,'(" *  MINI       - REDUCE OUTPUT BY ONLY PRINTING FLAGGED ATOMS")')
  if (mozyme) then
    if (myword(allkey, " VEC"))     write (iw,'(" *  VECTORS    - FINAL MOLECULAR ORBITALS TO BE PRINTED")')
  else
    if (myword(allkey, " VEC"))     write (iw,'(" *  VECTORS    - FINAL EIGENVECTORS TO BE PRINTED")')
  end if
  if (myword(allkey, " EIGEN"))   write (iw,'(" *  EIGEN      - CONVERT LMO''s INTO EIGENVECTORS")')
  if (myword(allkey, " CARTAB"))  write (iw,'(" *  CARTAB     - PRINT ALL THE CHARACTER TABLES")')
  if (myword(allkey, " CHARST"))  write (iw,'(" *  CHARST     - PRINT DETAILS OF SUBROUTINE CHARST")')
  if (myword(allkey, " DCART"))   write (iw,'(" *  DCART      - PRINT DETAILS OF SUBROUTINE DCART")')
  if (myword(allkey, " DERI1"))   write (iw,'(" *  DERI1      - PRINT DETAILS OF SUBROUTINE DERI1")')
  if (myword(allkey, " DERI2"))   write (iw,'(" *  DERI2      - PRINT DETAILS OF SUBROUTINE DERI2 ")')
  if (myword(allkey, " DERITR"))  write (iw,'(" *  DERITR     - PRINT DETAILS OF SUBROUTINE DERITR")')
  if (myword(allkey, " DERIV"))   write (iw,'(" *  DERIV      - PRINT DETAILS OF SUBROUTINE DERIV")')
  if (myword(allkey, " DERNVO"))  write (iw,'(" *  DERNVO     - PRINT DETAILS OF SUBROUTINE DERNVO")')
  if (myword(allkey, " MAKVEC"))  write (iw,'(" *  MAKVEC     - PRINT DETAILS OF SUBROUTINE MAKVEC")')
  if (myword(allkey, " DIIS"))    write (iw,'(" *  DIIS       - PRINT DETAILS OF SUBROUTINE DIIS  ")')
  if (myword(allkey, " FLEPO"))   write (iw,'(" *  FLEPO      - PRINT DETAILS OF SUBROUTINE FLEPO ")')
  if (myword(allkey, " FMAT"))    write (iw,'(" *  FMAT       - PRINT DETAILS OF SUBROUTINE FMAT  ")')
  if (myword(allkey, " DFORCE"))  write (iw,'(" *  DFORCE     - PRINT FORCE MATRIX OVER CARTESIAN COORDINATES")')
  if (myword(allkey, " HCORE"))   write (iw,'(" *  HCORE      - PRINT DETAILS OF SUBROUTINE HCORE ")')
  if (myword(allkey, " MOLDAT"))  write (iw,'(" *  MOLDAT     - PRINT DETAILS OF SUBROUTINE MOLDAT")')
  if (myword(allkey, " FREQCY"))  write (iw,'(" *  FREQCY     - PRINT DETAILS OF SUBROUTINE FREQCY")')
  if (myword(allkey, " ITER"))    write (iw,'(" *  ITER       - PRINT DETAILS OF SUBROUTINE ITER  ")')
  if (myword(allkey, " LINMIN"))  write (iw,'(" *  LINMIN     - PRINT DETAILS OF SUBROUTINE LINMIN")')
  if (myword(allkey, " MOLSYM"))  write (iw,'(" *  MOLSYM     - PRINT DETAILS OF SUBROUTINE MOLSYM")')
  if (myword(allkey, " RAMA"))    write (iw,'(" *  RAMA       - PRINT RAMACHANDRA ANGLES FOR PROTEIN RESIDUES")')
  if (myword(allkey, " SYMOIR"))  write (iw,'(" *  SYMOIR     - PRINT DETAILS OF SUBROUTINE SYMOIR")')
  if (myword(allkey, " GROUP"))   write (iw,'(" *  GROUP      - PRINT DETAILS OF SUBROUTINE GROUP ")')
  if (myword(allkey, " SYMTRZ"))  write (iw,'(" *  SYMTRZ     - PRINT DETAILS OF SUBROUTINE SYMTRZ")')
  if (myword(allkey, " DENS"))    write (iw,'(" *  DENSITY    - FINAL DENSITY MATRIX TO BE PRINTED")')
  if (myword(allkey, " SPIN"))    write (iw,'(" *  SPIN       - FINAL UHF SPIN MATRIX TO BE PRINTED")')
  if (myword(allkey, " PRNT"))    write (iw,'(" *  PRNT       - EXTRA PRINTING IN EF ROUTINE")')
  if (myword(allkey, " DISP"))    write (iw,'(" *  DISP       - PRINT DISPERSION AND HYDROGEN BOND ENERGIES")')
  if (myword(allkey, " DEP ")) then
10410 format (" *  DEP        - OUTPUT FORTRAN CODE FOR BLOCK-DATA")
10411 format (" *             THIS KEYWORD CANNOT BE USED IN THIS VERSION OF MOPAC.")
    write (iw, 10410)
    write (iw, 10411)
  end if
  if (myword(allkey, " TIMES"))  write (iw,'(" *  TIMES      - TIMES OF VARIOUS STAGES TO BE PRINTED")')
  if (mozyme)then
    if (myword(allkey, " BONDS"))  write (iw,'(" *  BONDS      - NON-HYDROGEN BOND-ORDERS TO BE PRINTED")')
    if (myword(allkey, " ALLBO"))  write (iw,'(" *  ALLBONDS   - ALL SIGNIFICANT BOND-ORDERS TO BE PRINTED")')
  else
    if (myword(allkey, " BONDS"))  write (iw,'(" *  BONDS      - FINAL BOND-ORDER MATRIX TO BE PRINTED")')
  end if
  if (myword(allkey, " SYBYL"))  write (iw,'(" *  SYBYL      - OUTPUT SYBYL FILE")')
  if (myword(allkey, " FOCK"))   write (iw,'(" *  FOCK       - LAST FOCK MATRIX TO BE PRINTED")')
  if (myword(allkey, " LARGE"))  write (iw,'(" *  LARGE      - EXPANDED OUTPUT TO BE PRINTED")')
  log = .false.
  if (myword(allkey, " LOG") .or. &
    (index(keywrd," 0SCF") /= 0 .and. index(keywrd," OLDGEO") /= 0 .and. index(keywrd," PDBOUT") /= 0) .or. &
    index(keywrd, " ADD-H") /= 0 .or. index(keywrd, " HEADER") /= 0) then
    iw0 = 0
    log = .true.
    if (index(keywrd, " LOG") /= 0) write (iw,'(" *  LOG        - GENERATE A LOG FILE")')
  end if
  if (myword(allkey, " NOLOG")) then
    write (iw,'(" *  NOLOG      - SUPPRESS LOG FILE")')
    iw0 = -1
  end if
  if (myword(allkey, " HTML"))   then
    if ((index(keywrd, " STEP=") == 0 .or. index(keywrd, " POINT") == 0) .and. &
    index(keywrd, " IRC") + index(keywrd, " DRC") == 0 .or. index(keywrd, " 0SCF") /= 0) then
      write (iw,'(" *  HTML       - WRITE HTML SCRIPT TO READ PDB FILE USING JSMOL")')
      if (index(keywrd, " PDBOUT") == 0) call l_control("PDBOUT", len_trim("PDBOUT"), 1)
    else
      write (iw,'(" *  HTML       - WRITE HTML SCRIPT TO GENERATE ANIMATION USING JSMOL")')
      if (index(keywrd, " PDBOUT") /= 0) then
        write (iw,'(" *             - TO DISPLAY ENERGIES, DELETE KEYWORD ""PDBOUT""")')
      else
        write (iw,'(" *             - TO DISPLAY IN PDB FORMAT, ADD ""PDBOUT""")')
      end if
    end if
    if (index(keywrd, " PDBOUT") /= 0) then
      if (maxtxt < 26) then
        if (index(keywrd," RESIDUES") == 0) then
          if (index(keywrd, " 0SCF") /= 0) then
             call l_control("RESIDUES", len_trim("RESIDUES"), 1)
          end if
        end if
      end if
    end if
    if ( .not. is_PARAM) then
!
!  Check that file-name is okay
!
      do j = len_trim(input_fn) - 5, 2, -1
        if (input_fn(j:j) == "/" .or. input_fn(j:j) == backslash) exit
      end do
      do i = j, len_trim(input_fn) - 5
        if (input_fn(i:i) == "'") exit
      end do
      if (i < len_trim(input_fn) - 6) then
         call mopend("When HTML is present, the file-name must not contain an apostrophe.")
         write(iw,'(/10x,a)')"(File name = """//input_fn(:len_trim(input_fn) - 5)//""")"
      end if
    end if
  end if
  if (myword(allkey, " NORJSMOL")) write (iw,'(" *  HTML(NORES)- SUPPRESS LIST OF RESIDUES IN HTML FILE")')
  if (myword(allkey, " AIGOUT")) write (iw,'(" *  AIGOUT     - IN ARC FILE, INCLUDE AB-INITIO GEOMETRY")')
  if (myword(allkey, " GRAP"))   write (iw,'(" *  GRAPH      - GENERATE FILE FOR GRAPHICS")')
  if (myword(allkey, " 1ELEC"))  write (iw,'(" *  1ELECTRON  - FINAL ONE-ELECTRON MATRIX TO BE PRINTED")')
  if (myword(allkey, " INTERP")) write (iw,'(" *  INTERP     - PRINT DETAILS OF CAMP-KING CONVERGER")')
  if (myword(allkey, " HTML"))   write (iw,'(" *  HTML       - WRITE HTML SCRIPT TO READ PDB FILE USING JSMOL")')
  if (myword(allkey, " SUPER"))  write (iw,'(" *  SUPER      - PRINT SUPERDELOCALIZABILITIES")')
  if (myword(allkey, " ALLVEC")) write (iw,'(" *  ALLVEC     - PRINT ALL EIGENVECTORS")')
  if (myword(allkey, " H-PRIO")) write (iw,'(" *  H-PRIOR    - HEAT OF FORMATION TAKES PRIORITY IN DRC")')
  if (myword(allkey, " X-PRIO")) write (iw,'(" *  X-PRIOR    - GEOMETRY CHANGES TAKE PRIORITY IN DRC")')
  if (myword(allkey, " T-PRIO")) write (iw,'(" *  T-PRIOR    - TIME TAKES PRIORITY IN DRC")')
  if (myword(allkey, " COMPFG")) write (iw,'(" *  COMPFG     - PRINT HEAT OF FORMATION CALC''D IN COMPFG")')
  if (myword(allkey, " MECI"))   write (iw,'(" *  MECI       - M.E.C.I. WORKING TO BE PRINTED")')
  if (myword(allkey, " HESSIAN"))write (iw,'(" *  HESSIAN    - WRITE OUT HESSIAN FROM GEOMERY OPTIMIZATION")')
  if (myword(allkey, " LOCAL"))  write (iw,'(" *  LOCALIZE   - LOCALIZED ORBITALS TO BE PRINTED")')
  if (myword(allkey, " BANANA")) write (iw,'(" *  BANANA     - MAKE LOCALIZED ORBITALS WITH BANANA BONDS")')
  if (myword(allkey, " RABBIT")) write (iw,'(" *  RABBIT     - MAKE LOCALIZED ORBITALS WITH RABBIT EARS")')
  if (myword(allkey, " MULLIK")) write (iw,'(" *  MULLIK     - THE MULLIKEN ANALYSIS TO BE PERFORMED")')
  if (myword(allkey, " PI "))    write (iw,'(" *  PI         - BONDS MATRIX, SPLIT INTO SIGMA-PI-DELL", &
   & " COMPONENTS, TO BE PRINTED")')
  if (myword(allkey, " ECHO"))   write (iw,'(" *  ECHO       - ALL INPUT DATA TO BE ECHOED BEFORE RUN")')
  if (myword(allkey, " DEBUG ")) write (iw,'(" *  DEBUG      - DEBUG OPTION TURNED ON")')
  if (myword(allkey, " POTWRT")) write (iw,'(" *  POTWRT     - WRITE OUT ELECTRIC POT. DATA TO FILE 21")')
  if (myword(allkey, " PMEP"))   write (iw,'(" *  PMEP       - GENERATE THE WANG-FORD ELECTROSTATIC POTENTIAL")')
  if (myword(allkey, " PRTMEP")) write (iw,'(" *  PRTMEP     - MEP CONTOUR DATA OUTPUT TO <FILENAME>.mep")')
  if (myword(allkey, " MINMEP")) write (iw,'(" *  MINMEP     - PRINT MEP MINIMA IN THE PLANE DEFINED")')
  if (myword(allkey, " PL"))     write (iw,'(" *  PL         - MONITOR CONVERGENCE IN DENSITY MATRIX")')
  if (myword(allkey, " PKA"))    write (iw,'(" *  PKA        - PRINT pKa FOR MOST ACIDIC HYDROGEN")')
  if (myword(allkey, " POPS"))   write (iw,'(" *  POPS       - PRINT SCF ATOMIC ORBITAL POPULATIONS")')
  if (myword(allkey, " BCC"))    write (iw,'(" *  BCC        - THE SYSTEM IS BODY-CENTERED CUBIC")')
  if (myword(allkey, " EIGS"))   write (iw,'(" *  EIGS       - PRINT ALL EIGENVALUES IN ITER")')
  if (myword(allkey, " HYPERF")) write (iw,'(" *  HYPERFINE  - HYPERFINE COUPLING CONSTANTS TO BE PRINTED")')
  if (myword(allkey, " THERMO")) write (iw,'(" *  THERMO     - THERMODYNAMIC QUANTITIES TO BE CALCULATED")')
  if (myword(allkey, " Z=")) then
    write (iw,'(" *  Z=",I2,"     - NUMBER OF FORMULA UNITS IN UNIT CELL")')Nint(Reada(keywrd,Index(keywrd," Z=") + 2))
    if (index(keywrd, " MERS") == 0) call mopend(" *  KEYWORD Z REQUIRES KEYWORD MERS TO BE USED")
  end if
  if (myword(allkey, " K=")) then
    i = Index (keywrd, " K=")
    write (iw,' (" *   K=        - BRILLOUIN ZONE STRUCTURE TO BE CALCULATED")')
    write (iw,'(" *             STEP SIZE IN SAMPLING ZONE:", f8.4)') reada (keywrd, i)
    i = Index (keywrd(i:), ",") + i
    write (iw,'(" *             NO. OF ATOMS IN FUNDAMENTAL UNIT CELL:", i6)') Nint (reada (keywrd, i))
  end if
  return
end subroutine wrtout
subroutine wrtwor (allkey)
  use molkst_C, only: tleft, tdump, keywrd, natoms, numat, mozyme, pdb_label, line
  use chanel_C, only: iw
  implicit none
  character (len=1000), intent (inout) :: allkey
  character :: ch, ch4*4
  character (len=7) :: chrono
  integer :: i, ii, j, k
  double precision :: time, sum_1, sum_2
  logical, external :: myword
  double precision, external :: reada
#ifdef _OPENMP
  integer :: max_threads
  integer, external :: omp_get_num_procs
#endif
  intrinsic Index, Min, Nint, Max
  if (myword(allkey, " EIGINV"))     write (iw,'(" *  EIGINV     - USE HESSIAN EIGENVALUE REVERSION IN EF")')
  if (myword(allkey, " NONR"))       write (iw,'(" *  NONR       - DO NOT USE NEWTON-RAPHSON STEP IN EF")')
  if (myword(allkey, " SNAP"))       write (iw,'(" *  SNAP       - INCREASE PRECISION OF SYMMETRY ANGLES")')
  if (myword(allkey, " PULAY"))      write (iw,'(" *  PULAY      - PULAY''S METHOD TO BE USED IN SCF")')
  if (myword(allkey, " CAMP"))       write (iw,'(" *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED")')
  if (myword(allkey, " KING"))       write (iw,'(" *  CAMP,KING- THE CAMP-KING CONVERGER TO BE USED")')
  if (myword(allkey, " LET"))        write (iw,'(" *  LET        - OVERRIDE SOME SAFETY CHECKS")')
  if (myword(allkey, " OLDGEO"))     write (iw,'(" *  OLDGEO     - PREVIOUS GEOMETRY TO BE USED")')
  if (myword(allkey, " OLDFPC"))     write (iw,'(" *  OLDFPC     - OLD FUNDAMENTAL PHYSICAL CONSTANTS TO BE USED")')
  if (myword(allkey, " OLD_HESS"))   write (iw,'(" *  OLD_HESS   - USE THE OLD HESSIAN MATRIX")')
  if (myword(allkey, " PREC"))       write (iw,'(" *  PRECISE    - TIGHTER CRITERIA TO BE USED")')
  if (myword(allkey, " NOANCI"))     write (iw,'(" *  NOANCI     - DO NOT USE ANALYTICAL C.I. DERIVATIVES")')
  if (myword(allkey, " DFP"))        write (iw,'(" *  DFP        - USE DAVIDON FLETCHER POWELL OPTIMIZER")')
  if (myword(allkey, " XYZ"))        write (iw,'(" *  XYZ        - CARTESIAN COORDINATE SYSTEM TO BE USED")')
  if (myword(allkey, " RESTART"))    write (iw,'(" *  RESTART    - CALCULATION RESTARTED")')
  if (myword(allkey, " RSCAL"))      write (iw,'(" *  RSCAL      - SCALE P-RFO STEP IN EF TO TRUST RADIUS")')
  if (myword(allkey, " SYMAVG"))     write (iw,'(" *  SYMAVG     - AVERAGE SYMMETRY EQUIVALENT ESP CHARGES")')
  if (myword(allkey, " STO3G"))      write (iw,'(" *  STO3G      - DEORTHOGONALIZE ORBITALS IN STO-3G BASIS")')
  if (myword(allkey, " FORCE "))     write (iw,'(" *  FORCE      - FORCE CALCULATION SPECIFIED")')
  if (myword(allkey, " FORCETS"))    write (iw,'(" *  FORCETS    - VERIFY THAT TRANSITION STATE IS GENUINE")')
  if (myword(allkey, " BFGS"))       write (iw,'(" *  BFGS       - USE THE BFGS GEOMETRY OPTIMIZER")')
  if (myword(allkey, " LBFGS"))      write (iw,'(" *  LBFGS      - USE THE LBFGS GEOMETRY OPTIMIZER")')
  if (myword(allkey, " EF"))         write (iw,'(" *  EF         - USE EF ROUTINE FOR MINIMUM SEARCH")')
  if (myword(allkey, " TS"))         write (iw,'(" *  TS         - USE EF ROUTINE FOR TS SEARCH")')
  if (myword(allkey, " LOCATE-TS"))  write (iw,'(" *  LOCATE-TS  - LOCATE A TRANSITION STATE USING REACTANTS AND PRODUCTS")')
  if (myword(allkey, " NOSYM"))      write (iw,'(" *  NOSYM      - POINT-GROUP SYMMETRY SET TO C1")')
  if (myword(allkey, " SMOOTH"))     write (iw,'(" *  SMOOTH     - IN A GRID CALCULATION, REMOVE COMPUTATIONAL ARTIFACTS")')
  if (myword(allkey, " AUTOSYM"))    write (iw,'(" *  AUTOSYM    - SYMMETRY TO BE IMPOSED AUTOMATICALLY")')
  if (myword(allkey, " 1SCF "))      write (iw,'(" *  1SCF       - DO 1 SCF AND THEN STOP ")')
  if (myword(allkey, " SIGMA"))      write (iw,'(" *  SIGMA      - GEOMETRY TO BE OPTIMIZED USING SIGMA.")')
  if (myword(allkey, " NLLSQ"))      write (iw,'(" *  NLLSQ      - GRADIENTS TO BE MINIMIZED USING NLLSQ.")')
  if (myword(allkey, " SADDLE"))     write (iw,'(" *  SADDLE     - TRANSITION STATE TO BE OPTIMIZED")')
  if (myword(allkey, " 0SCF"))       write (iw,'(" *  0SCF       - AFTER READING AND PRINTING DATA, STOP")')
  if (myword(allkey, " QPMEP "))     write (iw,'(" *  QPMEP      - CHARGES DERIVED FROM WANG-FORD TYPE AM1", " MEP")')
  if (myword(allkey, " BIGSCF"))     write (iw,'(" *  BIGSCF     - DO INITIAL FULL SCF WHEN RESTARTING JOB")')
  if (myword(allkey, " WILLIA"))     write (iw,'(" *  WILLIAMS   - USE WILLIAMS SURFACE")')
  if (myword(allkey, " INT "))       write (iw,'(" *  INT        - INTERNAL COORDINATE SYSTEM TO BE USED")')
  if (myword(allkey, " PECI"))       write (iw,'(" *  PECI       - SINGLE AND PAIRED ELECTRON", " EXCITATIONS ONLY")')
!
! Number of threads
!
  if (myword(allkey, " THREADS")) then
    i = Index (keywrd, " THREADS")
    i = nint(reada(keywrd, i))
    i = max(i,1)
#ifdef _OPENMP
    max_threads = omp_get_num_procs()
    if (i == 1) then
      write (iw,'(" *  THREADS=1  - MULTI-THREADING NOT USED")')
    else
      write (iw,'(" *  THREADS    - USE UP TO ", i0, " THREADS OUT OF ", i0)') i, max_threads
    end if
#else
    write (iw,'(" *  THREADS    - INACTIVE (THREAD CONTROLS DISABLED)")')
#endif
  end if
  !                   Times
  chrono = "SECONDS"
  time = 1.0d0
  tleft = 172800.d0 ! two days
  if (myword(allkey, " T=")) then
    i = Index (keywrd, " T=")
    tleft = reada (keywrd, i)
    do j = i + 3, 1000
      if (j == 1000 .or. keywrd(j+1:j+1) == " ") go to 1000
    end do
1010 if (tleft < 99999.9d0) then
       ch = char(ichar("4") + int(log10(tleft)))
       write (iw,'(" *  T=         - A TIME OF", f'//ch//'.1, " ", a, " REQUESTED")') tleft, trim(chrono)
     else
       i = int(log10(tleft)) + 4
       ch4(1:1) = char(ichar("0") + 1)
       ch4(2:2) = char(ichar("0") + i - 10)
       write (iw,'(" *  T=         - A TIME OF", f'//ch4(1:2)//'.1, " ", a, " REQUESTED")') tleft, trim(chrono)
    end if
!
!  Limit time to 999 weeks
!
    tleft = Min (6.0479d8, tleft*time)
    go to 1020
1000  ch = keywrd(j:j)
    if (ch == "M") then
      chrono = "MINUTES"
      time = 60.0d0
    end if
    if (ch == "H") then
      chrono = "HOURS"
      time = 3600.0d0
    end if
    if (ch == "D") then
      chrono = "DAYS"
      time = 86400.0d0
    end if
    if (ch == "W") then
      chrono = "WEEKS"
      time = 7*86400.0d0
    end if
    go to 1010
  else
    ch = char(ichar("4") + int(log10(tleft)))
    write (iw,'(" *  T=         - A TIME OF", f'//ch//'.1, " ", a7, " REQUESTED")') tleft, "SECONDS"
  end if
1020 time = 1.0d0
  chrono = "SECONDS"
  tdump = 7200.d0
  if (myword(allkey, " DUMP")) then
    i = Index (keywrd, " DUMP")
    tdump = reada (keywrd, i)
    do j = i + 6, 1000
      if (j == 1000 .or. keywrd(j+1:j+1) == " ") go to 1030
    end do
1040 if (tdump < 99999.9d0) then
      ch = char(ichar("4") + int(log10(tdump)))
      write (iw,'(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", f'//ch//'.1, " ", a, " REQUESTED")') tdump, trim(chrono)
    else
      i = int(log10(tdump)) + 4
      ch4(1:1) = char(ichar("0") + 1)
      ch4(2:2) = char(ichar("0") + i - 10)
      write (iw,'(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", f'//ch4(1:2)//'.1, " ", a, " REQUESTED")') tdump, trim(chrono)
    end if
    tdump = tdump * time
    go to 1050
1030  ch = keywrd(j:j)
    if (ch == "M") then
      chrono = "MINUTES"
      time = 60.d0
    end if
    if (ch == "H") then
      chrono = "HOURS"
      time = 3600.d0
    end if
    if (ch == "D") then
      chrono = "DAYS"
      time = 86400.d0
    end if
    if (ch == "W") then
      chrono = "WEEKS"
      time = 7*86400.d0
    end if
    go to 1040
  else
    ch = char(ichar("4") + int(log10(tdump)))
    write (iw,'(" *  DUMP=N     - RESTART FILE WRITTEN EVERY", f'//ch//'.1, " ", a7, " REQUESTED")') tdump, "SECONDS"
  end if
1050 if (index(keywrd, " LOCATE-TS") /= 0) then
       if (abs(tleft - 7200) < 0.1d0) tleft = 3.d6
       tdump = 3.d6
     end if
  if (myword(allkey, " CYCLES")) then
    sum_1 = reada (keywrd, Index (keywrd, " CYCLES"))
    ch = char(ichar("2") + int(log10(sum_1)))
    write (iw,'(" *  CYCLES=    - DO A MAXIMUM OF", i'//ch//', " STEPS")') Nint(sum_1)
  end if
  if (myword(allkey, " BIGCYCLES")) then
    sum_1 = reada (keywrd, Index (keywrd, " BIGCYCLES"))
    ch = char(ichar("2") + int(log10(sum_1)))
    write (iw,'(" *  BIGCYCLES  - DO", i'//ch//', " BIG STEPS")') Nint(sum_1)
  end if
!
!**********************************************************************
!
!                   Geometric Quantities
!
  if (myword(allkey, " ALT_A")) then
    write (iw, '(" *  ALT_A=", a1, "    - USE ALTERNATIVE ATOMS WHERE APPROPRIATE")') keywrd (Index (keywrd, " ALT_A")+7: &
     & Index (keywrd, " ALT_A")+7)
    end if
  if (myword(allkey, " ALT_R")) then
    write (iw, '(" *  ALT_R=", a1, "    - USE ALTERNATIVE RESIDUES WHERE APPROPRIATE")') keywrd (Index (keywrd, " ALT_R")+7: &
     & Index (keywrd, " ALT_R")+7)
    end if
  if (myword(allkey, " VDW")) then
    i = Index (keywrd, " VDW")
    j = Index (keywrd(i:), ")") + i
    write (iw, '(" *  ", a, "  USER-DEFINED ATOMIC OR VAN DER WAALS RADII")') keywrd(i + 1:j)
  end if
  if (index(allkey, " METAL") /= 0) then
    j = Index (keywrd, " METAL")
!
!  Force an equals sign in, if one is not present
!
    if (keywrd(j + 6:j + 6) == "(") then
      line = keywrd(:j + 5)//"="//trim(keywrd(j + 6:))
      keywrd = trim(line)
    end if
    i = Index (keywrd, " METAL=")
    if (i == 0) then
!
! The word "METAL" on its own.
!
      write (iw,'(" *  METAL      - METALS ARE DEFINED AS BEING FULLY IONIC")')
    else
!
! The word "METAL" in the construction "METAL=(text)"
!
      j = Index (keywrd(i:), ") ") + i
!
! Write out the labels of all atoms that are explicitely defined in PDB or Jmol format
!

      if (index(keywrd(i:j), '"') /= 0) then
        write (iw, '(" *  METAL        INDIVIDUAL ATOMS DEFINED AS BEING FULLY IONIC")')
        line = keywrd(i + 8:j - 2)
!
! Convert PDB and Jmol format into atom-numbers
!
        do
          k = index(keywrd(i:j), '"')
          if (k == 0) exit
          call txt_to_atom_no(keywrd(i:j), k, .false.)
        end do
!
! Now write out atoms and atom numbers
!
        j = Index (keywrd(i:), ") ") + i
        write (iw, '(" *",15x,a)')"("//trim(line)//") = "//keywrd(i + 8:j - 2)
      else
        j = Index (keywrd(i:), ") ") + i
        write (iw,'(" *  METAL      - THE FOLLOWING ELEMENTS ARE DEFINED AS BEING FULLY IONIC: ", a)') &
          keywrd(i + 8:j - 2)
      end if
    end if
    if (.not. myword(allkey, " METAL")) return ! dummy call to use "myword"
  end if
  if (numat < natoms .and. pdb_label) then
    if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE") /= 0) then
      call mopend("DUMMY ATOMS CANNOT BE PRESENT WHEN GEOMETRY IS BEING CHANGED")
      if (index(keywrd, " RESEQ") /= 0) write(iw,'(10x,a)')"Geometry is changed by keywrd: RESEQ"
      if (index(keywrd, " ADD-H") /= 0) write(iw,'(10x,a)')"Geometry is changed by keywrd: ADD-H"
      if (index(keywrd, " SITE=(")  /= 0) write(iw,'(10x,a)')"Geometry is changed by keywrd: SITE"
      write(iw,'(10x,a)')"(A simple way to remove dummy atoms is to add keyword ""XYZ"")"
    end if
  end if
  i = index(allkey, " CVB")
  if (i > 0) then
!
! Replace spaces by "_"
!
    j = index(allkey(i + 2:), ") ") + i + 1
    do k = i + 2, j
      if (allkey(k:k) == " ") allkey(k:k) ="_"
    end do
    if (index(keywrd, "ADD-H") /= 0) then
      do
        if (index(keywrd(i:k), "[") /= 0) then
          i = i + index(keywrd(i:k), "[")
          j = index(keywrd(i:k), ".") + i
          if (index(keywrd(i:j), ":") == 0) then
            call mopend("CVB keyword must use chain letters when ADD-H is used")
            write(iw,'(a)')"(Either run ""RESIDUES"" first "// &
            "to assign chain letters, or use CVB after ADD-H has been run.)"
            exit
          else
            i = j
          end if
        else
          exit
        end if
      end do
    end if
  end if
  if (myword(allkey, " CVB")) then
    i = Index (keywrd, " CVB")
    j = Index (keywrd(i:), ")") + i
    if (j - i < 24) then
      write (iw, '(" *  ", a, "  USER-DEFINED EXPLICIT COVALENT BONDS")') keywrd(i + 1:j)
    else
        write (iw, '(" *  USER-DEFINED EXPLICIT COVALENT BONDS:")')
        write (iw, '(" *    ",a)') keywrd(i + 1:j)
    end if
    if (index(keywrd, " RESEQ") /= 0) then
      call mopend("CVB must not be present when RESEQ is used")
      write(iw,*)
    end if
  end if
  if (myword(allkey, " TRANS=")) then
    write (iw, '(" *  TRANS=     - ", i4, " VIBRATIONS ARE TO BE DELETED FROM THE THERMO CALCULATION")') &
     Nint (reada (keywrd, Index (keywrd, " TRANS=")))
  else if (myword(allkey, " TRANS")) then
    write (iw, '(" *  TRANS      - THE REACTION VIBRATION TO BE DELETED FROM THE THERMO CALCULATION")')
  end if
  if (myword(allkey, "MEP=")) then
    i = Nint (reada (keywrd, Index (keywrd, "MEP=")))
    if (i == 1) then
      write (iw, '(" *  MEP=1      - MEP IN A CUBIC GRID")')
    else
      write (iw, '(" *  MEP=2      - MEP IN CONNOLLY SURFACE")')
    end if
  end if
!
!   How hessian matrix is updated.
!
  if (myword(allkey, " IUPD")) then
    ii = Nint (reada (keywrd, Index (keywrd, " IUPD=")))
    if (ii == 0) then
10730   format (" *  IUPD=0     - HESSIAN WILL NOT BE UPDATED")
      write (iw, 10730)
    end if
    if (ii == 1) then
10740   format (" *  IUPD=1     - HESSIAN WILL BE UPDATED USING POWELL")
      write (iw, 10740)
    end if
    if (ii == 2) then
10750   format (" *  IUPD=2     - HESSIAN WILL BE UPDATED USING BFGS")
      write (iw, 10750)
    end if
  end if
!
!   How starting hessian matrix is obtained.
!
  if (myword(allkey, " HESS=")) then
    ii = Nint (reada (keywrd, Index (keywrd, " HESS=")))
    if (ii == 0) then
10760   format (" *  HESS=0     - DIAGONAL HESSIAN USED AS INITIAL GUESS")
      write (iw, 10760)
    end if
    if (ii == 1) then
10770   format (" *  HESS=1     - INITIAL HESSIAN WILL BE CALCULATED")
      write (iw, 10770)
    end if
    if (ii == 2) then
10780   format (" *  HESS=2     - INITIAL HESSIAN READ FROM DISK")
      write (iw, 10780)
    end if
    if (ii == 3) then
10790   format (" *  HESS=3     - INITIAL HESSIAN READ FROM INPUT")
      write (iw, 10790)
    end if
  end if
!
!   Which "normal" mode to follow in geometry
!
  if (myword(allkey, " MODE=")) then
10800 format (" *  MODE=      - FOLLOW HESSIAN MODE", i3, " TOWARD TS")
    write (iw, 10800) Nint (reada (keywrd, Index (keywrd, "MODE=")))
  end if
  if (myword(allkey, " RECALC")) then
10810 format (" *  RECALC=    - DO", i4, " CYCLES BETWEEN HESSIAN RECALC")
    write (iw, 10810) Nint (reada (keywrd, Index (keywrd, "RECALC")))
  end if
  if (myword(allkey, " RMAX")) then
10820 format (" *  RMAX=      - MAX. CALC./PRED. ENERGY STEP IN TS", f7.3)
    write (iw, 10820) reada (keywrd, Index (keywrd, " RMAX="))
  end if
  if (myword(allkey, " RMIN")) then
10830 format (" *  RMIN=      - MIN. CALC./PRED. ENERGY STEP IN TS", f7.3)
    write (iw, 10830) reada (keywrd, Index (keywrd, " RMIN="))
  end if
  if (myword(allkey, " DDMAX")) then
10840 format (" *  DDMAX=     - MAXIMUM TRUST RADIUS", f7.3, " ANG/RAD")
    write (iw, 10840) reada (keywrd, Index (keywrd, " DDMAX="))
  end if
  if (myword(allkey, " DDMIN")) then
10850 format (" *  DDMIN=     - MINIMUM TRUST RADIUS", f7.3, " ANG/RAD")
    write (iw, 10850) reada (keywrd, Index (keywrd, " DDMIN="))
  end if
  if (myword(allkey, " DMAX")) then
10860 format (" *  DMAX=      - STARTING TRUST RADIUS", f7.3, " ANG/RAD")
    write (iw, 10860) reada (keywrd, Index (keywrd, "DMAX="))
  end if
  if (myword(allkey, " OMIN")) then
10870 format (" *  OMIN=      - MINIMUM EIGENVECTOR OVERLAP IN TS", f7.3)
    write (iw, 10870) reada (keywrd, Index (keywrd, "OMIN="))
  end if
  if (myword(allkey, " GNORM")) then
    sum_2 = reada (keywrd, Index (keywrd, " GNORM"))
    ch4 = char(ichar("6") + max(0, min(2,int(log10(sum_2 + 1.d-10)))))//".3"
    write (iw, '(" *  GNORM      - EXIT WHEN GRADIENT NORM DROPS BELOW", f'//ch4//')') sum_2
  end if
  if (myword(allkey, " DELTA")) then
10881 format (" *  DELTA=     - EXIT WHEN ENERGY DIFFERENCE DROPS BELOW ", g10.3)
    write (iw, 10881) reada (keywrd, Index (keywrd, " DELTA"))
  end if
  if (myword(allkey, " BAR")) then
10890 format (" *  BAR=       - REDUCE BAR LENGTH BY A MAX. OF", f7.2)
    write (iw, 10890) reada (keywrd, Index (keywrd, " BAR"))
  end if
  if (myword(allkey, " SLOG=")) then
    i = Index (keywrd, " SLOG=")
10900 format (" *  SLOG=n     - DEFAULT STEP SIZE IN BFGS:", f5.2)
    write (iw, 10900) reada (keywrd, i)
  else if (myword(allkey, " SLOG")) then
10910 format (" *  SLOG       - DEFAULT STEP SIZE IN BFGS: 0.25")
    write (iw, 10910)
  end if
!**********************************************************************
!
!                      SCF quantities
!
!   SCF Criterion
!
  if (myword(allkey, " RELSCF")) then
    if (mozyme) then
      sum_1 = 0.001d0
    else
      sum_1 = 0.0001d0
    end if
    sum_2 = reada (keywrd, Index (keywrd, " RELSCF"))
    ch4 = char(ichar("7") + max(0, min(2,int(log10(sum_2)))))//".4"
    write (iw, '(" *  RELSCF     - DEFAULT SCF CRITERION MULTIPLIED BY", f'//ch4//',/, &
       & " *              (DEFAULT SCF CRITERION =", f8.5, ")")') sum_2, sum_1
  end if
  if (myword(allkey, " SCFCRT")) then
10930 format (" *  SCFCRT     - DEFAULT SCF CRITERION REPLACED BY", g12.3)
    write (iw, 10930) reada (keywrd, Index (keywrd, " SCFCRT"))
  end if
  if (myword(allkey, " SHIFT")) then
10940 format (" *  SHIFT      - A DAMPING FACTOR OF", f8.2, " DEFINED")
    write (iw, 10940) reada (keywrd, Index (keywrd, " SHIFT"))
  end if
  if (myword(allkey, " FILL")) then
10950 format (" *  FILL=      - IN RHF CLOSED SHELL, FORCE M.O.", i3, &
   & " TO BE FILLED")
    write (iw, 10950) Nint (reada (keywrd, Index (keywrd, " FILL")))
  end if
  if (myword(allkey, " ITRY")) then
10960 format (" *  ITRY=      - DO A MAXIMUM OF", i6, " ITERATIONS FOR SCF")
    write (iw, 10960) Nint (reada (keywrd, Index (keywrd, " ITRY")))
  end if
  if (myword(allkey, " DAMP=")) then
10970 format (" *  DAMP=      - DAMP SCF OSCILLATIONS USING DAMP=", f6.2)
    write (iw, 10970) reada (keywrd, Index (keywrd, " DAMP=")+6)
  end if
  if (myword(allkey, " CUTOFF=")) then
10981 format (" *  CUTOFF=    - CUTOFF FOR NDDO AND DIPOLE INTERACTIONS =", f6.2, "A")
    write (iw, 10981) reada (keywrd, Index (keywrd, " CUTOFF=")+8)
  end if
  if (myword(allkey, " CUTOF2=")) then
10980 format (" *  CUTOF2=    - CUTOFF FOR NDDO INTERACTIONS =", f6.2, "A")
    write (iw, 10980) reada (keywrd, Index (keywrd, " CUTOF2=")+8)
  end if
  if (myword(allkey, " CUTOF1=")) then
10990 format (" *  CUTOF1=    - CUTOFF FOR DIPOLE INTERACTIONS =", f6.2, "A")
    write (iw, 10990) reada (keywrd, Index (keywrd, " CUTOF1=")+8)
  end if
  if (myword(allkey, " CUTOFP=")) then
11010 format (" *  CUTOFP=    - CUTOFF FOR POLYMER ELECTROSTATICS =", f6.2, "A")
    write (iw, 11010) reada (keywrd, Index (keywrd, " CUTOFP=")+8)
  end if
  if (myword(allkey, " CUTOFS=")) then
    write (iw, '(" *  CUTOFS=    - CUTOFF FOR OVERLAP INTEGRALS IN " &
    & //"SOLIDS =", f6.2, " A")') reada (keywrd, Index (keywrd, " CUTOFS=")+8)
  end if
  if (myword(allkey, " NLMO=")) then
11020 format (" *  NLMO=N     - AVERAGE NUMBER OF ATOMS PER LMO =", i4)
    write (iw, 11020) Nint (reada (keywrd, Index (keywrd, " NLMO=")))
  end if
  if (myword(allkey, " RELTHR")) then
11030 format (" *  RELTHR     - DEFAULT THRESHOLD FOR LMO's MULTIPLIED BY", &
   & f9.4,/, " *             (DEFAULT WAS 1.D-13)")
    write (iw, 11030) Max (reada (keywrd, Index (keywrd, " RELTHR")), 1.d-12)
  end if
  if (myword(allkey, " THRESH")) then
11040 format (" *  THRESH     - DEFAULT THRESHOLD FOR LMO's:", g12.3)
    write (iw, 11040) Max (reada (keywrd, Index (keywrd, " THRESH")), 1.d-25)
  end if
if (myword(allkey,' SETGPU=')) then
        i = index(keywrd,' SETGPU=')
        j = nint(reada(keywrd,i))
        write (iw,'(" *  SETGPU=   - YOUR CALCULATION WILL RUN IN THE GPU NUM. = ",i2)') j
  end if
  if (myword(allkey,' CPUTHREADS=')) then
     i = index(keywrd,' CPUTHREADS=')
     j = nint(reada(keywrd,i))
     write (iw,'(" *  CPUTHREADS=   - NUM. OF THREADS FOR CPU = ",i2)') j
  end if

  if (myword(allkey, ' FULLDIAG ')) write (iw,'(" *  FULLDIAG   - USE ONLY FULL DIAGONALIZATIONS IN SCF ")')
  return
  end subroutine wrtwor
