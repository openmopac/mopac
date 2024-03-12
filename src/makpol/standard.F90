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

subroutine bangle (xyz, i, j, k, angle)
!
!.. Implicit Declarations ..
    implicit none
!
!.. Formal Arguments ..
    double precision, dimension (3,*), intent (in) :: xyz
    integer, intent (in) :: i, j, k
    double precision, intent (out) :: angle
!
!.. Local Scalars ..
    double precision :: d2ij, d2ik, d2jk, temp, xy
!
!.. Intrinsic Functions ..
    intrinsic Acos, Sqrt
!
! ... Executable Statements ...
!
!********************************************************************
!
! BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
!        CARTESIAN COORDINATES ARE IN XYZ.
!
!********************************************************************
    d2ij = (xyz(1, i)-xyz(1, j))**2 + (xyz(2, i)-xyz(2, j))**2 &
         &+ (xyz(3, i)-xyz(3, j))**2
    d2jk = (xyz(1, j)-xyz(1, k))**2 + (xyz(2, j)-xyz(2, k))**2 &
         &+ (xyz(3, j)-xyz(3, k))**2
    d2ik = (xyz(1, i)-xyz(1, k))**2 + (xyz(2, i)-xyz(2, k))**2 &
         &+ (xyz(3, i)-xyz(3, k))**2
    xy = Sqrt (d2ij*d2jk)
    temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
    if (temp <-1.d0+1.d-12) then
      temp = -1.0d0
    end if
    if (temp > 1.d0-1.d-12) then
      temp = 1.d0
    end if
    angle = Acos (temp)
end subroutine bangle
subroutine dang (a1, a2, b1, b2, rcos)
!
!.. Implicit Declarations ..
    implicit none
!
!.. Formal Arguments ..
    double precision, intent (inout) :: a1, a2, b1, b2
    double precision, intent (out) :: rcos
!
!.. Local Scalars ..
    double precision :: anorm, bnorm, costh, sinth, zero
!
!.. Intrinsic Functions ..
    intrinsic Abs, Acos, Asin, Sqrt
!
! ... Executable Statements ...
!
!*********************************************************************
!
!    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!
!*********************************************************************
    zero = 1.0d-6
    if (Abs (a1) >= zero .or. Abs (a2) >= zero) then
      if (Abs (b1) >= zero .or. Abs (b2) >= zero) then
        anorm = 1.0d0 / Sqrt (a1**2+a2**2)
        bnorm = 1.0d0 / Sqrt (b1**2+b2**2)
        a1 = a1 * anorm
        a2 = a2 * anorm
        b1 = b1 * bnorm
        b2 = b2 * bnorm
        sinth = (a1*b2) - (a2*b1)
        costh = a1 * b1 + a2 * b2
        if (costh > 1.0d0) then
          costh = 1.0d0
        end if
        if (costh <-1.0d0) then
          costh = -1.0d0
        end if
        rcos = Acos (costh)
        if (Abs (rcos) >= 4.0d-4) then
          if (sinth > 0.d0) then
            rcos = 4.0d0 * Asin (1.0d00) - rcos
          end if
          rcos = -rcos
          return
        end if
      end if
    end if
    rcos = 0.0d0
end subroutine dang
!    ******************************************************************
double precision function digit (string, istart)
!     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
!     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
!     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
!     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
!
!
!.. Implicit Declarations ..
    implicit none
!
!.. Formal Arguments ..
    character (len=*), intent (in) :: string
    integer, intent (in) :: istart
!
!.. Local Scalars ..
    logical :: sign
    integer :: i, i0, i9, idig, idot, ineg, ipos, ispc, j, l, n
    double precision :: c1, c2, deciml
!
!.. Intrinsic Functions ..
    intrinsic Ichar, Len
!
! ... Executable Statements ...
!
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
    i0 = Ichar ("0")
    i9 = Ichar ("9")
    ineg = Ichar ("-")
    ipos = Ichar ("+")
    idot = Ichar (".")
    ispc = Ichar (" ")
!
    c1 = 0.d0
    c2 = 0.d0
    sign = .true.
    l = Len(string)
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
    idig = 0
    do i = istart, l
      n = Ichar (string(i:i))
      if (n >= i0 .and. n <= i9) then
        idig = idig + 1
        c1 = c1 * 1.d1 + n - i0
      else if (n == ineg .or. n == ipos .or. n == ispc) then
        if (n == ineg) then
          sign = .false.
        end if
      else
        go to 1000
      end if
    end do
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
1100 deciml = 1.d0
    do j = i + 1, l
      n = Ichar (string(j:j))
      if (n >= i0 .and. n <= i9) then
        deciml = deciml / 1.d1
        c2 = c2 + (n-i0) * deciml
      else if (n /= ispc) then
        exit
      end if
    end do
    go to 1200
1000 if (n == idot) go to 1100
!
!     PUT THE PIECES TOGETHER
1200 digit = c1 + c2
    if ( .not. sign) then
      digit = -digit
    end if
end function digit
subroutine dihed (xyz, i, j, k, l, angle)
    implicit none
    double precision, dimension (3,*), intent (in) :: xyz
    integer, intent (in) :: i, j, k, l
    double precision, intent (out) :: angle
    double precision :: cosa, cosph, costh, ddd, dist, sinph, sinth, xi1, xi2, &
   & xj1, xl1, xl2, yi1, yi2, yi3, yj1, yj2, yl1, yl2, yl3, yxdist, zi1, zj1, &
   & zl1
!********************************************************************
!
!      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
!            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
!            ARE IN ARRAY XYZ.
!
!     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
!           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
!
!********************************************************************
    xi1 = xyz(1, i) - xyz(1, k)
    xj1 = xyz(1, j) - xyz(1, k)
    xl1 = xyz(1, l) - xyz(1, k)
    yi1 = xyz(2, i) - xyz(2, k)
    yj1 = xyz(2, j) - xyz(2, k)
    yl1 = xyz(2, l) - xyz(2, k)
    zi1 = xyz(3, i) - xyz(3, k)
    zj1 = xyz(3, j) - xyz(3, k)
    zl1 = xyz(3, l) - xyz(3, k)
!      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
    dist = Sqrt (xj1**2+yj1**2+zj1**2)
    cosa = zj1 / dist
    if (cosa > 1.0d0) then
      cosa = 1.0d0
    end if
    if (cosa <-1.0d0) then
      cosa = -1.0d0
    end if
    ddd = 1.0d0 - cosa**2
    if (ddd > 0.0) then
      yxdist = dist * Sqrt (ddd)
      if (yxdist > 1.0d-6) then
        cosph = yj1 / yxdist
        sinph = xj1 / yxdist
        xi2 = xi1 * cosph - yi1 * sinph
        xl2 = xl1 * cosph - yl1 * sinph
        yi2 = xi1 * sinph + yi1 * cosph
        yj2 = xj1 * sinph + yj1 * cosph
        yl2 = xl1 * sinph + yl1 * cosph
!      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
        costh = cosa
        sinth = yj2 / dist
        go to 1000
      end if
    end if
    xi2 = xi1
    xl2 = xl1
    yi2 = yi1
    yl2 = yl1
    costh = cosa
    sinth = 0.d0
1000 yi3 = yi2 * costh - zi1 * sinth
    yl3 = yl2 * costh - zl1 * sinth
    call dang (xl2, yl3, xi2, yi3, angle)
    if (angle < 0.) then
      angle = 4.0d0 * Asin (1.0d00) + angle
    end if
    if (angle >= 6.2831853d0) then
      angle = 0.d0
    end if
end subroutine dihed
subroutine getdat
!
!
!
!
    use common_systm, only : jobnam_c, ifiles_1, ir, iw, path
    use common_jobnam, only : jobnam
!
!.. Implicit Declarations ..
    implicit none
!
!.. Local Scalars ..
    character (len=240) :: line
    logical :: exists
    integer, save :: i = 0, j
!
!
! ... Executable Statements ...
!
!
!   IFILES(1)=1 IS USED AS A FLAG TO INDICATE NORMAL PRINTING
!   IFILES(1)=0 WILL REDUCE THE AMOUNT OF OUTPUT
!
    ifiles_1 = 1



    ir = 25
    iw = 26
#ifdef MOPAC_F2003
    i = command_argument_count()
#else
    i = iargc()
#endif
    if (i == 0) then
      write(*,"(a)")"                        Program MAKPOL"
      write(*,"(a)")" "
      write(*,"(a)")"  MAKPOL constructs a MOPAC data set for a polymer, layer system, or solid."
      write(*,"(a)")"  It uses either a MOPAC-type data set that contains information about "
      write(*,"(a)")"  the size of the cluster to be built or a raw "".pdb"" or "".ent"" file.  "
      write(*,"(a)")"  For a description of MAKPOL, see:"
      write(*,"(a)")" "
      write(*,"(a)")" http://www.OpenMOPAC.net/Manual/makpol.html"
      write(*,"(a)")"  "
      write(0,'(10x,a)')" Press (return) to continue"
      read(5,*, iostat=i)
      stop
    end if
    if (i > 1) then
      write(*,'(//10x,a,/)') "Makpol only uses one argument, more than one was supplied"
      stop
    end if
#ifdef MOPAC_F2003
    call get_command_argument (1, jobnam)
#else
    call getarg (1, jobnam)
#endif
    write(*,*)trim(jobnam)
!
! allow for up to 3 commas in data set file name
!
    if (i >= 2) then
#ifdef MOPAC_F2003
      call get_command_argument (2, jobnam_c)
#else
      call getarg (2, jobnam_c)
#endif
      j = len_trim(jobnam)
      if (jobnam(j:j) == ",") jobnam(j:j) = " "
      jobnam=trim(jobnam)//","//jobnam_c
    end if
    if (i >= 3) then
#ifdef MOPAC_F2003
      call get_command_argument (3, jobnam_c)
#else
      call getarg (3, jobnam_c)
#endif
      jobnam=trim(jobnam)//","//jobnam_c
    end if
    if (i == 4) then
#ifdef MOPAC_F2003
      call get_command_argument (4, jobnam_c)
#else
      call getarg (4, jobnam_c)
#endif
      jobnam=trim(jobnam)//","//jobnam_c
    end if
!
! Split jobnam into the path and the job-name
!
    do j = len_trim(jobnam), 1, -1
      if (jobnam(j:j) == "/" .or. jobnam(j:j) == "\") exit
    end do
    if (j > 0) then
      path = jobnam(:j)
      jobnam = jobnam(j + 1:)
    else
      path = " "
    end if
    
    
    jobnam_c = jobnam
    call upcase(jobnam_c,len_trim(jobnam_c))
!
! See if the explicit file already exists  The explicit file must be of form
!   name.nnn
!
    line = jobnam_c
    if (jobnam_c(:5) == "MAKE_") jobnam = jobnam(6:)
    i = len_trim(jobnam) - 3
    if (i > 1) then
      if (jobnam(i:i) == ".") jobnam = jobnam(:i - 1)
    end if
!
!  At this point,
!  the complete name, in upper case, is in jobnam_c, and
!  the file-name, in the original case is in jobnam
!
!  First, check to see if the file, as defined in the filename, actually exists.
!
    exists = .false.
    i = len_trim(jobnam_c) - 3
    if (i > 1) then
      if (jobnam_c(i:i) == ".") then
         jobnam_c = trim(path)//trim(jobnam_c)
         inquire (file = trim(jobnam_c), exist = exists)       
      end if
    end if
    if (.not. exists) then
      line = jobnam
      call upcase(line, len_trim(line))
!
! Hunt for an appropriate file
!
!  Possible file names are:
!
!     MAKE_filename.dat
!          filename.mop
!     MAKE_fimename.mop
!          filename.ent
!

      jobnam_c = trim(path)//trim(line)//".dat"
      inquire (file = trim(jobnam_c), exist = exists)
      if (.not. exists) then
        jobnam_c = trim(path)//"MAKE_"//trim(line)//".dat"
        inquire (file = trim(jobnam_c), exist = exists)
      end if
      if (.not. exists) then
        jobnam_c = trim(path)//trim(line)//".mop"
        inquire (file = trim(jobnam_c), exist = exists)
      end if
      if (.not. exists) then
        jobnam_c = trim(path)//trim(line)//".mop"
        inquire (file = trim(jobnam_c), exist = exists)
      end if
      if (.not. exists) then
        jobnam_c = trim(path)//"MAKE_"//trim(line)//".mop"
        inquire (file = trim(jobnam_c), exist = exists)
      end if
      if (.not. exists) then
        jobnam_c = trim(path)//trim(line)//".ENT"
        inquire (file = trim(jobnam_c), exist = exists)
      end if
    end if
!
    open (unit=5, file=jobnam_c, status="OLD",iostat=i)
    if (i /= 0) then
      write(*,*)
      if ( .not. exists) then
        write(*,'(a)')"    No files of the type: '"//trim(jobnam)//"' were found."
      else
        write(*,'(a)')"    File: '"//trim(jobnam_c)//"' exists, but cannot be opened!"
      end if
      write(*,'(/3x,a)')" Press return to exit"
      read(5,'(a)') jobnam_c
      stop
    end if
    close (ir)
    rewind (5)
    line = trim(jobnam)
    call add_path(line)
    open (unit=ir, file=trim(line)//".temp_makpol", status="UNKNOWN")
    rewind (ir)
    i = 0
    do
      read (5, "(A)", end=1000, err=1000) line
      i = i + 1
      if (line(1:1) /= "*") then
        write (ir, "(A)") trim(line)
      end if
    end do
1000  write (ir, "(A)") " "
    rewind (ir)
    line = jobnam
!    do i = 1, len_trim(line)
!        if (line(i:i) == "_") line(i:i) = " "
!    end do
    if (i == 1) then
      write (iw, "(A)") " INPUT FILE MISSING OR EMPTY"
    else
      close (5)
      open (5, status="SCRATCH", form="UNFORMATTED")
    end if
    open (unit=iw, file=trim(path)//trim(line)//".mop", status="UNKNOWN", iostat = i)
    rewind (ir)
end subroutine getdat
subroutine getgeo (iread, labels, geo, xyz, lopt, na, nb, nc)
    use common_systm
    use common_keywrd
    use common_atomtx
    use common_sizes
    implicit none
!***********************************************************************
!
!   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
!          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
!
!  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
!             AMS    = DEFAULT ATOMIC MASSES.
!
!  WORKING    LINT = Coordinates should be INTERNAL
!             LXYZ = Coordinates should be CARTESIAN
!             (If .NOT.LINT and .NOT.LXYZ, use coordinates supplied)
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
    integer, intent (in) :: iread
    integer, dimension (numatm), intent (out) :: labels, na, nb, nc
    integer, dimension (3, numatm), intent (out) :: lopt
    double precision, dimension (3, numatm), intent (out) :: geo, xyz
    character, save :: nine = "9", space = " "
    character, save :: comma = ",", zero = "0"
    character :: turn
    character (len=2) :: ele
    character (len=120) :: string
    logical :: ircdrc, leadsp, lint, lmop, lxyz, velo
    integer :: i, icapa, icapz, icomma, iserr, j, k, khar, l, &
   & label, ndmy, nvalue, nvar_loc
    double precision :: const, degree,  real, sum, x1, x2, x3
    character (len=2), dimension (107), save :: elemnt
    integer, dimension (40) :: istart
    double precision, external :: reada
    data (elemnt(i), i=1, 107) / "H", "HE", "LI", "BE", "B", "C", "N", "O", &
         &"F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", &
         &"CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", &
         &"GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", &
         &"MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", &
         &"I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", &
         &"GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", &
         &"RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", &
         &"RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "&
         &BK", "CF", "XX", "FM", "MD", "CB", "++", "+", "--", "-", "TV" /
    ircdrc = (Index (keywrd, "IRC")+Index (keywrd, "DRC") /= 0)
    lmop = (Index (keywrd, " MOPAC") /= 0)
    icapa = Ichar ("A")
    icapz = Ichar ("Z")
    maxtxt = 0
    natoms = 0
    numat = 0
    iserr = 0
    do
!
      line = " "
        read (ir, "(A)", end=1000, err=1000) line
        go to 1100
1000    line = " "
1100  continue
!
!
      if (line == " ") go to 1400
      string = line(:6)
      if (index(string,"ATOM  ") + index(string,"HETATM") + index(string,"TITLE ") + index(string,"HEADER") + &
                                  index(string,"COMPND") + index(string,"SOURCE") + index(string,"KEYWDS") + &
          index(string,"HELIX ") + index(string,"SHEET ") + index(string,"REMARK") + index(string,"USER  ") + &
          index(string,"EXPDTA") + index(string,"AUTHOR") + index(string,"REVDAT") + index(string,"JRNL  ") + &
          index(string,"DBREF ") + index(string,"SEQRES") + index(string,"HET   ") + index(string,"HETNAM") + &
          index(string,"LINK  ") + index(string,"CRYST1") + index(string,"SCALE" ) + index(string,"ORIGX" ) + &
          index(string,"FORMUL") + index(string,"SEQRES") + index(string,"CONECT") /= 0) then
        natoms = -2
        return
      end if
      natoms = natoms + 1
!
!   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
!
      i = Index (line, "(")
      if (i /= 0) then
!
!  YES, ELEMENT IS LABELLED.
!
        k = Index (line, ")")
        txtatm(natoms) = line(i+1:k-1)
        maxtxt = Min (13, Max (maxtxt, k-i-1))
        string = line(1:i-1) // line(k+1:)
        line = string
      else
        txtatm (natoms) = " "
      end if
!   CLEAN THE INPUT DATA
      call upcase (line, 80)
      icomma = Ichar (comma)
      do i = 1, 80
        khar = Ichar (line(i:i))
        if (khar == icomma) then
          line(i:i) = space
        end if
      end do
!
!   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      do i = 1, 10
        istart(i) = 80
      end do
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER AND STORE IN ISTART
      leadsp = .true.
      nvalue = 0
      do i = 1, 80
        if (leadsp .and. line(i:i) /= space) then
          nvalue = nvalue + 1
          istart(nvalue) = i
        end if
        leadsp = (line(i:i) == space)
      end do
!
!  If optimization flags not present, adjust istart to suit
!
      if (nvalue == 4) then
        istart(6) = istart(4)
        istart(4) = istart(3)
        istart(3) = 1        
      end if
!
! ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
!
      string = line(istart(1) :istart(2)-1)
      if (string(1:1) >= zero .and. string(1:1) <= nine) then
!  ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
        label = Nint (reada (string, 1))
        if (label == 0) then
!***********************************************************************
! ALL DATA READ IN, CLEAN UP AND RETURN
!***********************************************************************
          natoms = natoms - 1
          go to 1400
        else if (label < 0 .or. label > 107) then
          write (iw, "('  ILLEGAL ATOMIC NUMBER')")
          go to 1600
        end if
      else
!  ATOMIC SYMBOL USED
        real = Abs (reada (string, 1))
        if (real < 1.d-15) then
!   NO ISOTOPE
          ele = string(1:2)
        else
          if (string(2:2) >= zero .and. string(2:2) <= nine) then
            ele = string(1:1)
          else
            ele = string(1:2)
          end if
        end if
!   CHECK FOR ERROR IN ATOMIC SYMBOL
        if (ele(1:1) == "-" .and. ele(2:2) /= "-") then
          ele (2:2) = " "
        end if
        do i = 1, 107
          if (ele == elemnt(i)) go to 1200
        end do
        if (ele(1:1) == "X") then
          label = 99
          go to 1300
        else
          write (iw, "('  UNRECOGNIZED ELEMENT NAME: (',A,')')") ele
          go to 1600
        end if
1200    label = i
      end if
!
! ALL O.K.
!
1300  if (label /= 99) then
        numat = numat + 1
      end if
      if (natoms > numatm) then
        write (iw, "(//10X,'**** MAX. NUMBER OF ATOMS ALLOWED:' ,I4)") numatm
        return
      end if
      labels(natoms) = label
      geo(1, natoms) = reada (line, istart(2))
      geo(2, natoms) = reada (line, istart(4))
      geo(3, natoms) = reada (line, istart(6))
      if (ircdrc) then
        turn = line(istart(3) :istart(3))
        if (turn == "T") then
          lopt(1, natoms) = 1
          if (natoms == 1) then
            write (iw, "(A)") " IN DRC MONITOR" &
                 &// " POTENTIAL ENERGY TURNING POINTS"
          end if
        else
          lopt(1, natoms) = 0
        end if
        turn = line(istart(5) :istart(5))
        if (turn == "T") then
          lopt(2, natoms) = 1
        else
          lopt(2, natoms) = 0
        end if
        turn = line(istart(7) :istart(7))
        if (turn == "T") then
          lopt(3, natoms) = 1
        else
          lopt(3, natoms) = 0
        end if
      else
       if(istart(3) /= 1) then
          lopt(1, natoms) = Nint (reada (line, istart(3)))
          lopt(2, natoms) = Nint (reada (line, istart(5)))
          lopt(3, natoms) = Nint (reada (line, istart(7)))
          do i = 3, 7, 2
            if (Ichar (line(istart(i) :istart(i))) >= icapa .and. &
         & Ichar (line(istart(i) :istart(i))) <= icapz .and. natoms > 1) then
              iserr = 1
            end if
          end do
        else
          lopt(1, natoms) = 1
          lopt(2, natoms) = 1
          lopt(3, natoms) = 1        
        end if
      end if
      sum = reada (line, istart(8))
      if (sum > 1.d7) then
        write (iw, "(A,I5,A)") " THE CONNECTIVITY OF ATOM", natoms, " IS FAULTY"
        return
      end if
!
!  If number is not an integer, then zero it out.  This
!  prevents an atomic charge being confused for a connectivity
!
      if (Abs (sum-Nint (sum)) > 1.d-5) then
        sum = 0.d0
      end if
      na(natoms) = Nint (sum)
      sum = reada (line, istart(9))
      if (Abs (sum-Nint (sum)) > 1.d-5) then
        sum = 0.d0
      end if
      nb(natoms) = Nint (sum)
      sum = reada (line, istart(10))
      if (Abs (sum-Nint (sum)) > 1.d-5) then
        sum = 0.d0
      end if
      nc(natoms) = Nint (sum)
      
!
!  SPECIAL CASE OF USERS FORGETTING TO ADD DIHEDRAL DATA FOR ATOM 3
!
      if (lmop .and. natoms == 3) then
        if (lopt(3, 3) == 2) then
          na(3) = 1
          nb(3) = 2
          geo(3, 3) = 0.d0
          lopt(3, 3) = 0
        else if (lopt(3, 3) == 1 .and. Abs (geo(3, 3)-2.d0) < 1.d-4) then
          na(3) = 2
          nb(3) = 1
          geo(3, 3) = 0.d0
          lopt(3, 3) = 0
        end if
      end if
      if ((lopt(1, natoms) > 1 .or. lopt(2, natoms) > 1 &
           &.or. lopt(3, natoms) > 1) .and. natoms > 1) then
        iserr = 1
      end if
      if (iserr == 1) exit
    end do
!  CHECK TO SEE IF GEOMETRY IS IN PDB FORMAT
!
    if (Index (line, "ATOM")+Index (line, "HETATM") /= 0) then
      natoms = -2
      return
    end if
!
!  MUST BE GAUSSIAN GEOMETRY INPUT
!
    do i = 2, natoms
      do k = 1, 3
        j = Nint (geo(k, i))
        if (Abs (geo(k, i)-j) > 1.d-5) then
!
!   GEOMETRY CANNOT BE GAUSSIAN
!
          return
        end if
      end do
    end do
    natoms = -1
    return
1400 lint = .false.
    lxyz = (index(keywrd, " XYZ") > 0 )
!
!   If MOPAC, set connectivity of atoms 2 and 3
!
1600 if (lmop) then
      na(2) = 1
      if (na(3) == 0) then
        na(3) = 2
        nb(3) = 1
      end if
    end if
!
!  Convert to radians, if in internal coordinates
!
    const = 1.7453292519943d-02
    do l = 1, natoms
      if (na(l) /= 0) then
        geo(2, l) = geo(2, l) * const
        geo(3, l) = geo(3, l) * const
      end if
    end do
!
!   Check that connectivity is O.K.
!
    do i = 4, natoms
      if (na(i) /= 0) then
        j = 0
        if (nb(i) < 1 .or. nc(i) < 1) then
          j = 1
        end if
        if (na(i) >= i .or. nb(i) >= i .or. nc(i) >= i) then
          j = 1
        end if
        if (na(i) == nb(i) .or. nb(i) == nc(i) .or. na(i) == nc(i)) then
          j = 1
        end if
        if (j /= 0) then
!       write (iw, "(A,I5,A,3I5,A)") " Connectivity of atom", i, " (", &
!          &na(i), nb(i), nc(i), ") is faulty"
        end if
      end if
    end do
    na1 = 0
    if (lxyz) then
!
!    COORDINATES SHOULD BE CARTESIAN
!
!
!  Coordinates are internal,  therefore switch.
!
      if (natoms == 0) then
        return
      end if
      call gmetrn (geo, xyz, labels, na, nb, nc)
!
!  Get rid of dummy atoms
!
      numat = 0
      do i = 1, natoms
        if (labels(i) /= 99) then
          numat = numat + 1
          labels(numat) = labels(i)        
          do j = 1, 3
            geo(j, numat) = xyz(j, i)
          end do
        end if
      end do
      natoms = numat
      na(1:numat) = 0
      na1 = 99
    else if (lint) then
      call gmetrn (geo, xyz, labels, na, nb, nc)
!
!   Check that no two atoms have the same coordinates
!
      do i = 2, numat
        x1 = xyz(1, i)
        x2 = xyz(2, i)
        x3 = xyz(3, i)
        do j = 1, i - 1
          if (Abs (x1-xyz(1, j)) <= 1.d-5) then
            if (Abs (x2-xyz(2, j)) <= 1.d-5) then
              if (Abs (x3-xyz(3, j)) <= 1.d-5) then
!
!   ATOMS I AND J ARE COINCIDENT.  THIS IS A BUG.
!
                write (iw,*)
                write (iw, "(A,I5,A,I5,A,/,A,3F9.4)") "  Atoms", i, " and", j, &
               & " have the same ", "  Cartesian coordinates:", x1, x2, x3
              end if
            end if
          end if
        end do
      end do
!
!   Remove dummy atoms, if any
!
      k = 0
      l = 0
      numat = 0
      nvar_loc = 0
      do i = 1, natoms
        na(i) = 0
        do j = 1, 3
          if (lopt(j, i) == 1) then
            nvar_loc = nvar_loc + 1
          end if
        end do
        if (labels(i) /= 99 .and. labels(i) /= 107) then
          k = k + 1
        end if
        if (labels(i) /= 99) then
          l = l + 1
          labels(l) = labels(i)
          txtatm(l) = txtatm(i)
        end if
      end do
      natoms = l
      numat = k
      do i = 1, natoms
        if (labels(i) == 107) then
!
!  ONLY NEED ONE WARNING.  TV CANNOT OCCUR WITH PDB FILES
!
          write (iw,*)
          write (iw,*) " SYSTEMS WITH TV CANNOT BE RUN WITH 'INT."
        end if
      end do
!
!    COORDINATES SHOULD BE INTERNAL.
!
      degree = 1.d0
      call xyzint (xyz, na, nb, nc, degree, geo)
!
!  UNCONDITIONALLY SET FLAGS FOR INTERNAL COORDINATES
!
      do i = 1, 3
        do j = i, 3
          lopt(j, i) = 0
        end do
      end do
      if (Abs (geo(2, 3)-180.d0) < 1.d-4 .or. Abs (geo(2, 3)) < 1.d-4) then
        write (iw, "(A)") " DUE TO PROGRAM BUG, THE FIRST THREE " &
             &// "ATOMS MUST NOT LIE IN A STRAIGHT LINE."
      end if
    end if
    lint = (lint .or. lxyz)
end subroutine getgeo
subroutine gmetrn (geo, coord, labels, na, nb, nc)
    use common_systm, only : natoms, na1, iw, id
    use common_ctvec, only : tvec
    implicit none
    integer, dimension (natoms), intent (inout) :: labels, na, nb
    integer, dimension (natoms), intent (inout) :: nc
    double precision, dimension (3, natoms), intent (inout) :: coord
    double precision, dimension (3, natoms), intent (in) :: geo
    integer :: i, j, k, l, ma, mb, mc
    double precision :: ccos, cosa, cosd, coskh, cosph, costh, rbc, &
   & sina, sind, sinkh, sinph, sinth, xa, xb, xd, xpa, xpb, xpd, xqd, &
   & xrd, xyb, ya, yb, yd, ypa, ypd, yqd, yza, za, zb, zd, zpd, zqa, zqd
!***********************************************************************
!
!    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
!***IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
!
!     THREE SEPARATE OPTIONS EXIST WITHIN GMETRY. THESE ARE:
!    (A) IF NA1 IS EQUAL TO 99 (IMPOSSIBLE UNDER NORMAL CIRCUMSTANCES)
!        THEN GEO IS ASSUMED TO BE IN CARTESIAN RATHER THAN INTERNAL
!        COORDINATES, AND COORD IS THEN SET EQUAL TO GEO.
!    (B) IF STEP IS NON-ZERO (THIS IS THE CASE WHEN "SADDLE" IS USED)
!        THEN GEO IS FIRST MODIFIED BY SHIFTING THE INTERNAL COORDINATES
!        ALONG A RADIUS FROM GEOA TO PLACE GEO AT A DISTANCE STEP FROM
!        GEOA.
!    (C) NORMAL CONVERSION FROM INTERNAL TO CARTESIAN COORDINATES
!        IS DONE.
!
!  ON INPUT:
!         GEO    = ARRAY OF INTERNAL COORDINATES.
!         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
!         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
!
!  ON OUTPUT:
!         COORD  = ARRAY OF CARTESIAN COORDINATES
!
!***********************************************************************
!                                     OPTION (B)
    
  
!                                     OPTION (A)
    if (na1 == 99) then
      do i = 1, 3
!    USE NATOMS INSTEAD OF NUMAT BECAUSE OF TV
!$DOIT VBEST
        do j = 1, natoms
          coord(i, j) = geo(i, j)
        end do
!    USE NATOMS INSTEAD OF NUMAT BECAUSE OF TV
!$DOIT VBEST
      end do
    else
!                                     OPTION (C)
      do i = 1, Min (3, natoms)
        do j = 1, 3
          coord(j, i) = geo(j, i)
        end do
      end do
      if (natoms > 1) then
        if (na(2) == 1) then
          coord(1, 2) = coord(1, 1) + geo(1, 2)
          coord(2, 2) = coord(2, 1)
          coord(3, 2) = coord(3, 1)
        end if
        if (natoms /= 2) then
          ccos = Cos (geo(2, 3))
          if (na(3) == 1) then
            coord(1, 3) = coord(1, 1) + geo(1, 3) * ccos
            coord(2, 3) = coord(2, 1) + geo(1, 3) * Sin (geo(2, 3))
            coord(3, 3) = coord(3, 1)
          else if (na(3) == 2) then
            coord(1, 3) = coord(1, 2) - geo(1, 3) * ccos
            coord(2, 3) = coord(2, 2) + geo(1, 3) * Sin (geo(2, 3))
            coord(3, 3) = coord(3, 2)
          else
            coord(1, 3) = geo(1, 3)
            coord(2, 3) = geo(2, 3)
            coord(3, 3) = geo(3, 3)
          end if
          do i = 4, natoms
            mc = na(i)
            if (mc /= 0) then
              cosa = Cos (geo(2, i))
              mb = nb(i)
              xb = coord(1, mb) - coord(1, mc)
              yb = coord(2, mb) - coord(2, mc)
              zb = coord(3, mb) - coord(3, mc)
              rbc = xb * xb + yb * yb + zb * zb
              if (rbc < 1.d-16) then
!
!     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
!
                write (iw, "(A,I4,A,I4,A)") " ATOMS", mb, " AND", mc, &
                     &" ARE COINCIDENT"
                write (iw, "(A)") " THIS IS A FATAL ERROR, RUN " &
                     &// "STOPPED IN GMETRY"
              else
                rbc = 1.0d00 / Sqrt (rbc)
              end if
              do
                ma = nc(i)
                xa = coord(1, ma) - coord(1, mc)
                ya = coord(2, ma) - coord(2, mc)
                za = coord(3, ma) - coord(3, mc)
!
!     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!
                xyb = Sqrt (xb*xb+yb*yb)
                k = -1
                if (xyb <= 0.1d00) then
                  xpa = za
                  za = -xa
                  xa = xpa
                  xpb = zb
                  zb = -xb
                  xb = xpb
                  xyb = Sqrt (xb*xb+yb*yb)
                  k = + 1
                end if
!
!     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!
                costh = xb / xyb
                sinth = yb / xyb
                xpa = xa * costh + ya * sinth
                ypa = ya * costh - xa * sinth
                sinph = zb * rbc
                cosph = Sqrt (Abs (1.d00-sinph*sinph))
                zqa = za * cosph - xpa * sinph
!
!     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!
                yza = Sqrt (ypa**2+zqa**2)
                if (yza < 1.d-4) then
!
!   ANGLE TOO SMALL TO BE IMPORTANT
!
                  coskh = 1.d0
                  sinkh = 0.d0
                  go to 1000
                else
                  if (Abs (cosa) >= 0.9998d0) exit
                  if (yza >= 2.d-3) exit
!
!   Faulty atom: I; Bond lendth to I: NA(I); Angle to I: NB(I)
!
!   Set a new NC(I) and re-try
!
                  j = nc(i)
                  call renum (coord, na, nb, nc, i, natoms)
                  if (nc(i) == 0 .or. nc(i) == j) then
                    go to 1100
                  end if
                end if
              end do
              coskh = ypa / yza
              sinkh = zqa / yza
!
!     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!     NONE ARE NEGATIVE.
!     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
!
1000          sina = Sin (geo(2, i))
              sind = -Sin (geo(3, i))
              cosd = Cos (geo(3, i))
              xd = geo(1, i) * cosa
              yd = geo(1, i) * sina * cosd
              zd = geo(1, i) * sina * sind
!
!     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!
              ypd = yd * coskh - zd * sinkh
              zpd = zd * coskh + yd * sinkh
              xpd = xd * cosph - zpd * sinph
              zqd = zpd * cosph + xd * sinph
              xqd = xpd * costh - ypd * sinth
              yqd = ypd * costh + xpd * sinth
              if (k >= 1) then
                xrd = -zqd
                zqd = xqd
                xqd = xrd
              end if
              coord(1, i) = xqd + coord(1, mc)
              coord(2, i) = yqd + coord(2, mc)
              coord(3, i) = zqd + coord(3, mc)
            else
              coord(1, i) = geo(1, i)
              coord(2, i) = geo(2, i)
              coord(3, i) = geo(3, i)
            end if
          end do
          go to 1200
1100      write (iw, "(//20x,' CALCULATION ABANDONED AT THIS POINT')")
          if (nc(i) == 0) then
            write (iw, "(//10x,' THREE ATOMS BEING USED TO ','DEFINE THE',&
                 &/10x,' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS')")
            write (iw, "(10X,' NOT ZERO OR 180 DEGREES, ARE ,'&
                 &'IN AN ALMOST STRAIGHT')")
            write (iw, "(10x,' LINE.  THERE IS A HIGH PROBABILITY THAT THE',&
                 &/10x,' COORDINATES OF THE ATOM WILL BE INCORRECT.')")
            write (iw, "(//20x,'THE FAULTY ATOM IS ATOM NUMBER',i6)") i
          else
            write (iw,*) " The Cartesian geometry has become unreasonable"
          end if
            if (nc(i) == 0) then
            write (iw, "(//20X, 'CARTESIAN COORDINATES UP TO FAULTY ATOM')")
            write (iw, "(//5X,'I',12X,'X',12X,'Y',12X,' Z')")
            do j = 1, i
              write (iw, "(I6,F16.5,2F13.5)") j, (coord(k, j), k=1, 3)
            end do
            write (iw, "(//6x,' ATOMS',i3,',',i3,', AND',i3,' ARE WITHIN',&
                 &f7.4,' ANGSTROMS OF A STRAIGHT LINE')") mc, mb, ma, yza
          end if
          return
        end if
      end if
    end if
!
!***NOW REMOVE THE TRANSLATION VECTORS, IF ANY, FROM THE ARRAY COOR
!
1200 Continue
!    if (labels(natoms) == 98) then
!      write(iw,*) " The last atom must NOT be a dummy atom"
!      stop
!    end if
    l = 0
    do k = 1, natoms
      if (labels(k) == 107) then
        l = l + 1
        tvec(1, l) = coord(1, k)
        tvec(2, l) = coord(2, k)
        tvec(3, l) = coord(3, k)
        labels(k) = 98
      end if
    end do
    id = l
    j = 0
    do i = 1, natoms
      if (labels(i) /= 98) then
        j = j + 1
!$DOIT ASIS
        do k = 1, 3
          coord(k, j) = coord(k, i)
        end do
        labels(j) = labels(i)
      end if
    end do
    natoms = j
end subroutine gmetrn
    subroutine renum (coord, na, nb, nc, ii, natoms)
        implicit none
        integer, intent (in) :: ii, natoms
        integer, dimension (natoms), intent (in) :: na, nb
        integer, dimension (natoms), intent (out) :: nc
        double precision, dimension (3, natoms), intent (in) :: coord
        integer :: i, jj, nai, nbi
        double precision :: angle, rab, rmin, theta
!***********************************************************************
!
!  Renumber the NC of atom II.  On input, the angle NA(II)-NB(II)-NC(II)
!  is too near to 0 or 180 degrees.  Find a new atom for NC(II), so that
!  the angle will be acceptable (as large as possible)
!***********************************************************************
        nai = na(ii)
        nbi = nb(ii)
!
!   Theta = 45 degrees
!
        theta = 0.7853d0
        jj = 0
        rmin = 1.d10
        do
          do i = 1, ii - 1
            if (i /= nai .and. i /= nbi) then
              call bangle (coord, nai, nbi, i, angle)
              if (angle > 1.5707963d0) then
                angle = 2.d0 * Asin (1.d0) - angle
              end if
              if (angle >= theta) then
!
!   Angle is OK.  Now find atom of lowest distance
!
                rab = (coord(1, nbi)-coord(1, i))**2 &
                     &+ (coord(2, nbi)-coord(2, i))**2 &
                     &+ (coord(3, nbi)-coord(3, i))**2
                if (rab < rmin) then
                  jj = i
                  rmin = rab
                end if
              end if
            end if
          end do
          if (jj /= 0) then
!
!  Best NC is JJ; best angle is THMIN
!
            nc(ii) = jj
            exit
          end if
!
!   No atom inside the allowed angle - reduce the angle
!
          theta = theta * 0.5d0
          if (theta < 0.0174533d0) then
            nc(ii) = 0
            exit
          end if
        end do
    end subroutine renum
    subroutine minv (a, n, d)
!
!.. Implicit Declarations ..
        implicit none
!
!.. Formal Arguments ..
        integer, intent (in) :: n
        double precision, dimension (n*n), intent (inout) :: a
        double precision, intent (out) :: d
!
!.. Local Scalars ..
        integer :: i, ij, ik, iz, j, ji, jk, jp, jq, jr, k, ki, kj, kk, nk
        double precision :: biga, hold
!
!.. Local Arrays ..
        integer, dimension (1000) :: l, m
!
!.. Intrinsic Functions ..
        intrinsic Abs, Max, Min
!
! ... Executable Statements ...
!
!
!     SEARCH FOR LARGEST ELEMENT
!
        d = 1.0d0
        nk = -n
        do k = 1, n
          nk = nk + n
          l(k) = k
          m(k) = k
          kk = nk + k
          biga = a(kk)
          do j = k, n
            iz = n * (j-1)
            do i = k, n
              ij = iz + i
              if (Abs (biga) < Abs (a(ij))) then
                biga = a(ij)
                l(k) = i
                m(k) = j
              end if
            end do
          end do
!
!     INTERCHANGE ROWS
!
          j = l(k)
          if (j > k) then
            ki = k - n
            do i = 1, n
              ki = ki + n
              hold = -a(ki)
              ji = ki - k + j
              a(ki) = a(ji)
              a(ji) = hold
            end do
          end if
!
!     INTERCHANGE COLUMNS
!
          i = m(k)
          if (i > k) then
            jp = n * (i-1)
            do j = 1, n
              jk = nk + j
              ji = jp + j
              hold = -a(jk)
              a(jk) = a(ji)
              a(ji) = hold
            end do
          end if
!
!     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
!     CONTAINED IN BIGA)
!
          if (biga == 0.0d0) then
            d = 0.d0
            return
          end if
          do i = 1, n
            if (i /= k) then
              ik = nk + i
              a(ik) = a(ik) / (-biga)
            end if
          end do
!  REDUCE MATRIX
          do i = 1, n
            ik = nk + i
            hold = a(ik)
            ij = i - n
            do j = 1, n
              ij = ij + n
              if (i /= k .and. j /= k) then
                kj = ij - i + k
                a(ij) = hold * a(kj) + a(ij)
              end if
            end do
          end do
!
!     DIVIDE ROW BY PIVOT
!
          kj = k - n
          do j = 1, n
            kj = kj + n
            if (j /= k) then
              a(kj) = a(kj) / biga
            end if
          end do
!
!     PRODUCT OF PIVOTS
!
          d = Max (-1.d25, Min (1.d25, d))
          d = d * biga
!
!     REPLACE PIVOT BY RECIPROCAL
!
          a(kk) = 1.d0 / biga
        end do
!
!     FINAL ROW AND COLUMN INTERCHANGE
!
        k = n
        do
          k = (k-1)
          if (k <= 0) exit
          i = l(k)
          if (i > k) then
            jq = n * (k-1)
            jr = n * (i-1)
            do j = 1, n
              jk = jq + j
              hold = a(jk)
              ji = jr + j
              a(jk) = -a(ji)
              a(ji) = hold
            end do
          end if
          j = m(k)
          if (j > k) then
            ki = k - n
            do i = 1, n
              ki = ki + n
              hold = a(ki)
              ji = ki - k + j
              a(ki) = -a(ji)
              a(ji) = hold
            end do
          end if
        end do
    end subroutine minv
    subroutine nuChar (line, value, nvalue)
        implicit none
        character (len=80), intent (inout) :: line
        double precision, dimension (40), intent (out) :: value
        integer, intent (out) :: nvalue
        character, save :: comma = ",", space = " "
        character :: tab
        logical :: leadsp
        integer :: i
        integer, dimension (40) :: istart
        double precision, external :: reada
        tab = Char(9)
!
! CLEAN OUT TABS AND COMMAS
!
        do i = 1, 80
          if (line(i:i) == tab .or. line(i:i) == comma) then
            line(i:i) = space
          end if
        end do
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER
!
        leadsp = .true.
        nvalue = 0
        do i = 1, 80
          if (leadsp .and. line(i:i) /= space) then
            nvalue = nvalue + 1
            istart(nvalue) = i
          end if
          leadsp = (line(i:i) == space)
        end do
!
! FILL NUMBER ARRAY
!
        do i = 1, nvalue
          value(i) = reada (line, istart(i))
        end do
    end subroutine nuchar
    double precision function reada (string, istart)
!     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
!
        implicit none

        character (len=*), intent (in) :: string
        integer, intent (in) :: istart
        logical :: expnnt
        integer :: i, i0, i9, iadd, icapd, icape, idot, ineg, ipos, ismld, &
       & ismle, j, l, n
        double precision, external :: digit
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
        i0 = Ichar ("0")
        i9 = Ichar ("9")
        idot = Ichar (".")
        ineg = Ichar ("-")
        ipos = Ichar ("+")
        icapd = Ichar ("D")
        icape = Ichar ("E")
        ismld = Ichar ("d")
        ismle = Ichar ("e")
!
        l = Len(string)
!
!     FIND THE START OF THE NUMERIC FIELD
        do i = istart, l
          iadd = 0
          n = Ichar (string(i:i))
!
!       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
          if (n >= i0 .and. n <= i9) go to 1000
!
!       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
          if (n == ineg .or. n == ipos) then
            iadd = iadd + 1
            if (i+iadd > l) exit
            n = Ichar (string(i+iadd:i+iadd))
            if (n >= i0 .and. n <= i9) go to 1000
          end if
!
!       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
          if (n == idot) then
            iadd = iadd + 1
            if (i+iadd > l) exit
            n = Ichar (string(i+iadd:i+iadd))
            if (n >= i0 .and. n <= i9) go to 1000
          end if
        end do
!
!     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
        reada = 0.d0
        return
!
!     FIND THE END OF THE NUMERIC FIELD
1000    expnnt = .false.
        do j = i + 1, l
          iadd = 0
          n = Ichar (string(j:j))
!
!       CONTINUE SEARCH FOR END IF DIGIT FOUND
          if (n < i0 .or. n > i9) then
!
!       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
            if (n == ineg .or. n == ipos) then
              if ( .not. expnnt) go to 1100
              iadd = iadd + 1
              if (j+iadd > l) go to 1100
              n = Ichar (string(j+iadd:j+iadd))
              if (n >= i0 .and. n <= i9) cycle
            end if
            if (n == idot) then
              iadd = iadd + 1
              if (j+iadd > l) go to 1100
              n = Ichar (string(j+iadd:j+iadd))
              if (n >= i0 .and. n <= i9) then
                cycle
              else if (n == icape .or. n == ismle .or. n == icapd .or. n == &
             & ismld) then
                cycle
              end if
            end if
            if (n /= icape .and. n /= ismle .and. n /= icapd .and. n /= ismld) &
           & go to 1100
            if (expnnt) go to 1100
            expnnt = .true.
          end if
        end do
        j = l + 1
1100    n = Ichar (string(j-1:j-1))
        if (n == icape .or. n == ismle .or. n == icapd .or. n == ismld) then
          j = j - 1
        end if
!
!     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
        n = 0
        n = n + Index (string(i:j-1), "e")
        n = n + Index (string(i:j-1), "E")
        n = n + Index (string(i:j-1), "d")
        n = n + Index (string(i:j-1), "D")
        if (n == 0) then
          reada = digit (string(i:j-1), 1)
        else
          reada = digit (string(:i+n-2), i) * 1.d1**&
               &digit (string(:j-1), i+n)
        end if
    end function reada
    subroutine upcase (keywrd, n)
!
!.. Implicit Declarations ..
        implicit none
!
!.. Formal Arguments ..
        character (len=*), intent (inout) :: keywrd
        integer, intent (in) :: n
!
!.. Local Scalars ..
        character (len=241) :: keybuf
        integer :: i, icapa, iline, ilowa, ilowz, j
!
!.. Intrinsic Functions ..
        intrinsic Char, Ichar, Index
!
! ... Executable Statements ...
!
!
!  UPCASE WILL TAKE A CHARACTER STRING, IN KEYWRD, AND PUT IT INTO
!  UPPER CASE.  KEYWRD IS LIMITED TO 80 CHARACTERS
!
        icapa = Ichar ("A")
        ilowa = Ichar ("a")
        ilowz = Ichar ("z")
        keybuf = keywrd
        do i = 1, n
          iline = Ichar (keywrd(i:i))
          if (iline >= ilowa .and. iline <= ilowz) then
            keywrd(i:i) = Char(iline+icapa-ilowa)
          end if
!
!  Change tabs to spaces.  A tab is ASCII character 9.
!
          if (iline == 9) then
            keywrd (i:i) = " "
          end if
        end do
!
!   If the word EXTERNAL is present, do NOT change case of the following
!   character string.
!
        i = Index (keywrd, "EXTERNAL")
        if (i /= 0) then
          j = Index (keywrd(i+1:), " ") + i
          keywrd(i+9:j) = keybuf(i+9:j)
        end if
    end subroutine upcase
    subroutine wrttxt (iprt)
        use common_titles
        use common_keywrd
        use common_jobnam, only : jobnam
!.. Implicit Declarations ..
        implicit none
!
!.. Formal Arguments ..
        integer, intent (in) :: iprt
!
!.. Intrinsic Functions ..
        intrinsic Index
!
! ... Executable Statements ...
!
        write (iprt, "(A)") trim(keywrd)
        if (Index (koment, " NULL ") == 0 .and. koment /= " ") then
          write (iprt, "(A)") trim(koment)
        else
          write (iprt, "(A)") trim(jobnam)
        end if
        if (Index (title, " NULL ") == 0) &
        &  write (iprt, "(A)") trim(title)

       
    end subroutine wrttxt
    subroutine xyzgeo (xyza, numat, na, nb, nc, degree, geo)
        implicit none
        integer, intent (in) :: numat
        double precision, intent (in) :: degree
        integer, dimension (numat), intent (in) :: na, nb
        integer, dimension (numat), intent (inout) :: nc
        double precision, dimension (3, numat), intent (in) :: xyza
        double precision, dimension (3, numat), intent (inout) :: geo
        integer :: i, i1, ii, j, k, l
        double precision :: angl, r, sum, tol
!**********************************************************************
!
!   XYZGEO CONVERTS COORDINATES FROM CARTESIAN TO INTERNAL.
!
!     ON INPUT XYZ  = ARRAY OF CARTESIAN COORDINATES
!              NUMAT= NUMBER OF ATOMS
!              NA   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!                     BY DISTANCE
!              NB   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!                     BY ANGLE
!              NC   = NUMBERS OF ATOM TO WHICH ATOMS ARE RELATED
!                     BY DIHEDRAL
!
!    ON OUTPUT GEO  = INTERNAL COORDINATES IN ANGSTROMS, RADIANS,
!                     AND RADIANS
!
!**********************************************************************
        do i = 2, numat
          j = na(i)
          k = nb(i)
          l = nc(i)
          if (i == 493) then
            l = nc(i)
          end if
          if (i >= 3) then
            ii = i
            call bangle (xyza, ii, j, k, geo(2, i))
            geo(2, i) = geo(2, i) * degree
            if (i >= 4) then
!
!   MAKE SURE DIHEDRAL IS MEANINGLFUL
!
              call bangle (xyza, j, k, l, angl)
              tol = 0.2617994d0
!
!  If angle is within 0.015 degrees of a straight line, assume
!  the user intended this, and carry on
!
              if ((angl < 3.1415926d0-tol*1.d-3 .or. angl > tol*1.d-3) .and. &
             & (angl > 3.1415926d0-tol .or. angl < tol)) then
                do
!
!  ANGLE IS UNSATISFACTORY, LET'S SEARCH FOR ANOTHER ATOM FOR
!  DEFINING THE DIHEDRAL.
                  sum = 100.d0
                  do i1 = 1, ii - 1
                    r = (xyza(1, i1)-xyza(1, k))**2 &
                         &+ (xyza(2, i1)-xyza(2, k))**2 &
                         &+ (xyza(3, i1)-xyza(3, k))**2
                    if (r < sum .and. i1 /= j .and. i1 /= k) then
                      call bangle (xyza, j, k, i1, angl)
                      if (angl < 3.1415926d0-tol .and. angl > tol) then
                        sum = r
                        l = i1
                        nc(ii) = l
                      end if
                    end if
                  end do
                  if (sum > 99.d0 .and. tol > 0.1d0) then
!
! ANYTHING WITHIN 5 DEGREES?
!
                    tol = 0.087266d0
                  else
                    exit
                  end if
                end do
              end if
              call dihed (xyza, ii, j, k, l, geo(3, i))
              geo(3, i) = geo(3, i) * degree
            end if
          end if
          geo(1, i) = Sqrt ((xyza(1, i)-xyza(1, j))**2&
               &+ (xyza(2, i)-xyza(2, j))**2&
               &+ (xyza(3, i)-xyza(3, j))**2)
        end do
        geo(1, 1) = 0.d0
        geo(2, 1) = 0.d0
        geo(3, 1) = 0.d0
        geo(2, 2) = 0.d0
        geo(3, 2) = 0.d0
        geo(3, 3) = 0.d0
    end subroutine xyzgeo
    subroutine xyzint (xyz, na, nb, nc, degree, geo)
        use common_systm, only : numat
        implicit none
        double precision, dimension (3, numat), intent (in) :: xyz
        integer, dimension (numat), intent (inout) :: na, nb, nc
        double precision, intent (in) :: degree
        double precision, dimension (3, numat), intent (inout) :: geo
        integer :: i, l
        integer :: im1, j, k
        double precision :: r, sum, pi = 3.141592653589d0
        geo = 0.d0
        do i = 1, numat 
          im1 = i - 1 
          j = na(i) 
          if (j == 0) then
            na(i) = 2 
            nb(i) = 3 
            nc(i) = 4   
            if (i == 1) cycle           
            sum = 1.d30 
            do j = 1, im1 
              r = (xyz(1,i) - xyz(1,j))**2 + (xyz(2,i) - xyz(2,j))**2 + (xyz(3,i) - xyz(3,j))**2 
              if (.not.(r < sum .and. na(j) /= j .and. nb(j) /=j )) cycle  
              sum = r 
              k = j 
            end do 
  !
  !   ATOM I IS NEAREST TO ATOM K
  !
            na(i) = k 
            j = k
            if (i > 2) nb(i) = na(k)  
            nc(i) = nb(k) 
!
!   FIND ANY ATOM TO RELATE TO NA(I)
!
          end if         
          if (i > 3) then
!
!  Check that the angle na(i) - nb(i) - nc(i) is meaningful
!
            call bangle (xyz, na(i), nb(i), nc(i), sum) 
            if (sum < 1.d-2 .or. Abs(pi - sum) < 1.d-2) then
!
!  Angle is zero or 180 degrees.  Search for an atom nearest to 90 degrees
!
              r = 2.d0
              do k = 1, im1
                if (k == i .or. k == j) cycle
                call bangle (xyz, na(i), nb(i), k, sum) 
                if (Abs(pi*0.5d0 - sum)  < r) then
                  r =  Abs(pi*0.5d0 - sum)
                  l = k
                  end if
                  if (r < 0.5d0) exit ! angle is acceptable
              end do
              nc(i) = l
            end if
            call dihed (xyz, i, j, nb(i), nc(i), geo(3,i)) 
          end if
          geo(3,i) = geo(3,i)*degree
          if (i > 2) call bangle (xyz, i, j, nb(i), geo(2,i)) 
          geo(2,i) = geo(2,i)*degree
          geo(1,i) = sqrt((xyz(1,i)-xyz(1,j))**2 + &
                          (xyz(2,i)-xyz(2,j))**2 + &
                          (xyz(3,i)-xyz(3,j))**2) 
        end do    
        na(1) = 0 
        nb(1) = 0 
        nc(1) = 0 
        if (numat > 1) then 
          nb(2) = 0 
          nc(2) = 0 
          if (numat > 2) nc(3) = 0 
          na(2) = 1 
        end if 
  end subroutine xyzint
  subroutine add_path(file)
    use common_systm, only : path
    implicit none
    character, intent (inout) :: file*(*)
    integer :: n
    logical :: need_path
!
!  Test to see if the file "file" already has an absolute path.
! 
      need_path =  ((file(2:2) /= ":") .and. &
      ((file(1:1) /= "/") .and. (file(1:1) /= "\")) .and. ((file(2:2) /= "/") .and. (file(2:2) /= "\"))) 
    if (.not. need_path) return
!
!  Test to see if the job has an absolte path.
!    
      need_path = (path(2:2) == ":")
    if (.not. need_path) return      
!
!  file "file" does not include the path and the path is present in "path"
!  therefore add path to the file "file"
!
!  First check if a relative path is used
!
    n = 1
    do 
      if (file(1:3) /= "../" .and. file(1:3) /= "..\") exit
      n = n + 1
      file = trim(file(4:))
    end do
!
!  Is the file defined relative to the current folder?
!  If so, delete the definition.
!
    if (file(1:2) == "./" .or. file(1:2) == ".\") file = trim(file(3:))
!
!  Join the path to the start of the file
!
    file = trim(path)//trim(file)   
    return
  end subroutine add_path
  subroutine l_control(txt, nt, mode)
!
!   l_control has two modes:
!
! (A)  If mode =  1 the word(s) in txt is(are) added to keywrd.  If the word in txt already exists
!                   in keywrd, the word in txt would first be deleted from keywrd.
! (B)  If mode = -1 the word(s) in txt is(are) removed from keywrd
!
!  txt can consists of upper and lower case letters, "-" and "_"; the last two allow hyphenated words.
!
  use common_systm, only :  line
  use common_keywrd, only : keywrd
  implicit none
  integer :: nt, mode
  character :: txt*(nt), store*2000
  integer :: i, j, trim_len, mt
  character :: local_txt*(nt), ch*1
  local_txt = txt
  do
!
!  Parse each word, one at a time.
!
    do
      if (local_txt(1:1) /= " ") exit
      local_txt = local_txt(2:)
    end do
    ch = " "
    mt = 0
    do i = 1, nt
      mt = mt + 1
      if (local_txt(mt:mt) == '"') then ! Quickly run to the closing '"'
        do mt = mt + 1, nt
          if (local_txt(mt:mt) == '"') exit
        end do
      end if
      if (local_txt(mt:mt) == ' ') exit 
      if (mt == nt) exit
    end do
    line = local_txt(:mt)  
!
!  "line" holds a single keyword
!
    local_txt = local_txt(mt + 1:)
    if (line(1:1) >= "0" .and. line(1:1) <= "9") then
      i = mt + 1
      else
  !
  ! find the end of the defined keyword, e.g., in pi=3.1415, that would be "i"
  !
      do i = 1, mt
        if ((line(i:i) < "A" .or. line(i:i) > "Z") .and. (line(i:i) < "a" .or. line(i:i) > "z") .and. &
        line(i:i) /= "_" .and. line(i:i) /= "-") exit
      end do
    end if
    trim_len = i - 1  
!
!  If keyword already exists, then delete the keyword
!  Keyword can have an "=" sign, or be followed by a "(" or a '"'
!  in which case, find the end of the keyword.
!
    do
      if (keywrd == " ") exit
      i = index(keywrd, " "//line(:trim_len))
      if (i > 0) then
        j = i + trim_len + 1
        if (keywrd(j:j) /= " ") then
          if (keywrd(j:j) == "=") j = j + 1
          ch = " "
          if (keywrd(j:j) == "(") then
            ch = ")"
          else if (keywrd(j:j) == '"') then
            ch = '"'
          end if
          do
            j = j + 1
            if (keywrd(j:j) == ch) exit
          end do 
        end if
        keywrd = keywrd(:i)//keywrd(j + 1:)
      else
        exit
      end if
    end do
    store = " "
    if (mode == 1) then
      i = index(keywrd, store(:mt + 2))
      keywrd = keywrd(:i)//line(:mt)//trim(keywrd(i + mt + 1:))
    end if  
    if (local_txt == " ") exit
  end do
  return
  end subroutine l_control
