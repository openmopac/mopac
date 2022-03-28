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

 subroutine datinp
    use param_global_C, only : maxmol, Atom_pKas, ifiles_8, contrl, &
    maxatm, numvar, locvar, names, nmols, heats, dipls, hions, weight, &
    is_a_ref, refers, qn, spin_state, i_r, ihrefs, msdels, nref, &
    & pkas, lions, geotxt, norbss, nlecss, ncloss, molele, nopens, &
    nvars, comnts, titls, nlm61, fracts, natmss, numats, ndeps, &
    nnmos, geos, nfirss, nlasts, lablss, keys, nats, pas, pbs, &
    locs, locpas, idepfs, locdes, depmuls, nas, nbs, ncs, maxpat, maxsym, &
    maxpab, nfns, diffns, n2elecs, botlim, toplim, valvar, valold, &
    l123s, refdir, nnalpha_open, nnbeta_open
!
    use common_arrays_C, only : nat, labels, &
    & nfirst, nlast, na, nb, nc, geo, loc, &
    pdiag

    use elemts_C, only : elemnt

    use parameters_C, only : n_partyp_alpb
!
    use molkst_C, only : keywrd, moperr, title, natoms, method_mndo, &
    method_am1, method_pm6, method_mndod, nvar, koment, nclose, norbs, &
    nalpha, nelecs, nopen, ndep, fract, nbeta, numat, lm61, n2elec, &
    jobnam, errtxt, msdel, id, l1u, l2u, l3u, method_pm3, line, method_pm6_dh2, &
    method_pm7, method_PM6_DH_plus, nalpha_open, nbeta_open, method_pm6_d3h4, &
    method_pm6_d3_not_h4, method_pm6_d3h4x, backslash
!
    use meci_C, only : nmos
!
    USE molmec_C, only : nnhco
!
    USE symmetry_C, ONLY: idepfn, locdep, depmul, locpar
!
    use chanel_C, only : iw, ir, igpt, ilog, job_fn
!----------------------------------------------------------------------------
    implicit none
    character :: dr, gr, hr, jr
    character (len=200) :: name, refnam
    character (len=80), dimension (maxmol) :: safety
    logical :: case, exists, let, opend, refok, precise, large, all, survey, &
    & clean, hof_type, deadly = .false., l_only, method_pm5 = .false., &
    &  cp, solid, nocore
    integer :: i, iatm, icapa, igeo, ilin, iline, iloc, ilowa, ilowz, isym, &
   & itime, j, k, l, loop, misd, misg, mish, misi, &
   & mmmols, mol, ndips, nel, nions, nmolmx, nnmols, notok, icapz, &
   & m, n, lnmols, is, iref, faulty, max_atom, ni, nj, ne_excluded = 0, &
   itype, mers(3), nele(107)
    double precision :: sum, tnow, told, wtd, wtgeo, wth, wti
    logical, dimension (2, 107) :: elemok
    integer, dimension (4) :: jfiles
    integer, dimension (107,107) :: ccp = 0, bonds_ij = 0
    intrinsic Abs, Char, Ichar, Index, Int
    logical, external :: core_core_OK
!
!  The weights given to each reference datum depend on the elements present.
!  Weights are assigned according to the type of element
   double precision, dimension (107) :: element_weights
   character, external :: get_a_name*300
   integer, external :: end_of_keyword
   double precision, external :: seconds, reada
   double precision, parameter :: &
   & org = 1.d0,      &! Organic elements
   & mg1 = 0.9d0,     &! Main-group elements (important in biochemistry)
   & mg2 = 0.8d0,     &! Main-group elements (less important)
   & tms = 0.7d0,     &! Transition metals
   & oth = 1.d0! Other entities
!
!
!
!              H           Relative elements weights                                           He
!              Li  Be                                                      B   C   N   O   F   Ne
!              Na  Mg                                                      Al  Si  P   S   Cl  Ar
!              K   Ca  Sc             Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn   Ga  Ge  As  Se  Br  Kr
!              Rb  Sr  Y              Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd   In  Sn  Sb  Te  I   Xe
!              Cs  Ba  La Ce-Lu       Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg   Tl  Pb  Bi  Po  At  Rn
!              Fr  Ra  Ac Th Pa U     Np  Pu  Am  Cm  Bk  Cf  XX  +++ ---  Cb  ++  +   --  -   Tv
!
!
    data element_weights  &
           &/ org,                                                                            mg1, &!    2
           &  mg2,mg2,                                                    mg2,org,org,org,mg1,mg1, &!   10
           &  mg1,mg2,                                                    mg2,mg2,mg1,mg1,mg1,mg1, &!   18
           &  mg1,mg2,tms,           tms,tms,tms,tms,tms,tms,tms,tms,tms, mg2,mg2,mg2,mg2,mg1,mg1, &!   36
           &  mg2,mg2,tms,           tms,tms,tms,tms,tms,tms,tms,tms,tms, mg2,mg2,mg2,mg2,mg1,mg1, &!   54
           &  mg2,mg2,tms, 14*tms,   tms,tms,tms,tms,tms,tms,tms,tms,tms, mg2,mg2,mg2,mg2,mg1,mg1, &!   86
           &  mg2,mg2,tms, 3*tms,    tms,tms,tms,tms,tms,tms,oth,oth,oth, oth,oth,oth,oth,oth,oth /
!                                  /
!
! ... Executable Statements ...
!
!  Assign the maximum possible size to pas and pbs
!
    maxpab = 10000000
    do
      maxpab = (maxpab*6)/5
      allocate (pas(maxpab), pbs(maxpab), stat = i)
      if (i /= 0) then
        maxpab = (maxpab*5)/6
        if (allocated(pas)) deallocate(pas)
        allocate (pas(maxpab), pbs(maxpab), stat = i)
        exit
      end if
      if (maxpab > 40000000) exit
      deallocate (pas, pbs, stat = i)
    end do
    refnam = " "
    refok = .false.
!
!  JFILES(1) = Normal MOPAC output - this will be deleted.
!  JFILES(2) = Faulty MOPAC output - to be saved.
!  JFILES(3) = Normal PARAM output - to be saved.
!  JFILES(4) = Duplicate PARAM output - to be deleted.
!
    Atom_pKas = 0
    jfiles(1) = iw
    jfiles(2) = iw
    jfiles(3) = ifiles_8
    jfiles(4) = iw
    iloc = 0
    igeo = 0
    isym = 0
    iatm = 0
    ilin = 0
    is   = 0
    itype = 0
    ilowa = Ichar ("a")
    ilowz = Ichar ("z")
    icapa = Ichar ("A")
    icapz = Ichar ("Z")
    told = seconds (1)
    max_atom = 1000
    if (Index (contrl, " MAXAT") /= 0) max_atom = Nint(reada(contrl, Index(contrl," MAXAT")))

    let = (Index (contrl, " LET") /= 0)
    large = (Index(contrl, "LARGE") /= 0)
    precise = (Index(contrl, "PRECISE") /= 0)
    clean = (Index(contrl, "CLEAN") /= 0)
    all = (Index(contrl, " ALL") /= 0)
    survey = (Index (contrl, " SURVEY") /= 0)
    nocore = (Index (contrl, " NOCORE") /= 0)
!   nocore = .true.
    l_only = (Index (contrl, "ONLY") /= 0)
    solid = (index(contrl, " SOLID") /= 0 .and. max_atom == 1000)
    cp = (Index (contrl, " CP ") /= 0)
!   rewind (14)
!
! NOW FOLLOW THE WEIGHTING FACTORS FOR THE OPTIMISATION,
      wth = 1.d0!     WTH = WEIGHTING FACTOR FOR HEATS OF FORMATION, IN MOLES/KCAL
!    wtd = 20.d0!     WTD = WEIGHTING FACTOR FOR DIPOLE MOMENTS, IN DIPOLES(**-1)
!    wti = 10.d0!     WTI = WEIGHTING FACTOR FOR IONISATION POTENTIALS, IN EV(**-1)
      wtgeo = 0.4d0!   WTGEO = WEIGHTING FACTOR FOR GRADIENTS, IN KCAL/MOLE
!
!   For the SBIR work, reduce weights for dipoles and I.P.s
!
      wtd = 10.d0!     WTD = WEIGHTING FACTOR FOR DIPOLE MOMENTS, IN DIPOLES(**-1)
      wti = 5.d0!     WTI = WEIGHTING FACTOR FOR IONISATION POTENTIALS, IN EV(**-1)
!
!  If the user wants to alter these weights, then
!
    i = Index (keywrd, " RELH=")
    if (i /= 0) wth = max(wth * Reada(keywrd,i), 1.d-7)
    i = Index (keywrd, " RELD=")
    if (i /= 0) wtd = max(wtd * Reada(keywrd,i), 1.d-7)
    i = Index (keywrd, " RELI=")
    if (i /= 0) wti = max(wti * Reada(keywrd,i), 1.d-7)
    i = Index (keywrd, " RELG=")
    if (i /= 0) wtgeo = max(wtgeo * Reada(keywrd,i), 1.d-7)
    rewind (14)
!
!  Set ELEMOK(1,N) to TRUE if the element is to be parameterized
!
    i = maxatm
    do i = 1, 107
      elemok(1, i) = .false.
    end do
    do i = 1, numvar
      if (locvar(1,i) /= 41) elemok(1, mod(locvar(2, i),200)) = .true.
    end do

!
!  Establish which elements have parameters. If GSS(I) is non-zero
!  parameters are assumed to exist.
!
    j=1
    sum = 0.d0
    do i = 1, 98
      call extract_parameter (j, i, sum)
      elemok(2, i) = (Abs (sum) > 0.0001d0)
    end do
    elemok(2,57:71) = .true.
    elemok(2,87) = .true.
    elemok(2,85) = .true.
    do i = 99, 107
      elemok(1, i) = (numvar == 0)
      elemok(2, i) = .true.
    end do
    if (solid) elemok(1,107) = .true.! Allow Tv
    elemok(1,99)  = .true.
    elemok(1,98)  = .true.
!
!  Store KEYWRD in CONTRL.  At this point KEYWRD = keywords
!  of parameterization data set
!
    contrl = trim(keywrd)
    i = 0
    m = 0
    n = 0
    loop = Index(contrl, " SET=")
    if (loop /= 0) then
      n = loop
      contrl(loop:loop + 4) = ' '
      loop = loop + 5
    end if
    do
      if (loop /= 0) then
        moperr = .false.
        inquire(unit=ir, opened=opend)
        if (opend) close(ir)
        line = get_a_name(contrl(loop:), len_trim(contrl(loop:)))
        loop = loop + len_trim(line)
        if (contrl(loop + 1:loop + 1) == '"') loop = loop + 2
        if (contrl(loop:loop) == ";") then
          loop = loop + 1
        else
          contrl(n:loop) = " "
          loop = 0
        end if
        call add_path(line)
        open (unit=ir, file=trim(line), status="OLD", err=99)
        goto 89
  99    write(*,"(///,3a,///)")" Reference data file set: '",trim(line),"' missing."
        write(ifiles_8,"(///,3a,///)")" Reference data file set: '",trim(line),"' missing."
        call finish
89      continue
      else
        if (i /= 0) exit!  At least one set of data has been read in.
      end if

!
!   Read in all the names of the molecules to be used in the optimization
!

      do m = 1, maxmol
        read (ir, "(A)", end=1200) keywrd
        if (keywrd(1:1) == "*") cycle
!
!  Force molecule name into lower case
!
        do j = 1, 80
          if (keywrd(j:j) /= " ") go to 1100
        end do
        exit
1100    do k = 80, j, -1
          if (keywrd(k:k) /= " ") exit
        end do
        line = keywrd
        do l = j, k
          iline = Ichar (keywrd(l:l))
          if (iline >= Ichar ("A") .and. iline <= Ichar ("Z")) then
            line (l:l) = Char (iline+Ichar ("a")-Ichar ("A"))
          end if
        end do
        if (line(:3) == "end") cycle
        if (nocore .and. index(line,"core") /= 0) cycle
        i = i + 1
        names(i) = keywrd(j:k)
        j = Index(line(j:k), ".mop")!  Remove ".mop" if it exists
        if (j /= 0) names(i)(j:j+3) = " "
        do j = 1, i - 1
          if (names(i) == names(j)) then
            write(ifiles_8,"(5x,3a)")" File: """,trim(names(i)),""" already specified."
            i = i - 1
            exit
          end if
        end do
      end do
      cycle
 1200 continue
      if (m == 1) then
        call finish !  no data on file-names
      end if
    end do
!
!  Close input data set to prevent it being damaged by next set
!  of READS
!
    inquire(unit=ir, opened=opend)
    if (opend) then
      close (ir)
    end if
    if (Ichar(names(i)(1:1)) == 0) i = i - 1
    nmolmx = i
    write (ifiles_8, "(//10X,' FILE NAMES OF MOLECULES USED',/)")
    write (ifiles_8, "(2x,a)") (names(i), i=1, nmolmx)
    write (ifiles_8, "('1',//,'     MOLECULE         ',25x,'           H-EXP','   WT'&
   &,'   D-EXP   WT   I-EXP   WT  GEO-WT')")
    ndips = 0
    nions = 0
    mmmols = 0
    mish = 0
    misi = 0
    misd = 0
    misg = 0
    notok = 0
    faulty = 0
!
!  Define directory where the reference data are.
!
    nref = 0
    k = Index (contrl, " REF=")
    if (k /= 0) then
      k = k + 5
      j = end_of_keyword(contrl, len_trim(contrl), k)
!
! k = start of reference data directory list
! j = end of list.
! in between are the names of the reference directories, separated by ";"
!
      do l = 1, 20
        i = Index(contrl(k:j),";")
        if (i /= 0) then
          nref = nref + 1
          refdir(nref) = trim(get_a_name(contrl(k:j), len_trim(contrl(k:j))))//"/"
          n = len_trim(refdir(nref))
          if (refdir(nref)(n:n) /= "/") refdir(nref)(n + 1:n + 1) = "/"
          k = k + i
        end if
        if (i == 0) exit
      end do
      nref = nref + 1
      refdir(nref) = trim(get_a_name(contrl(k:j), j - k + 1))
      n = len_trim(refdir(nref))
      if (refdir(nref)(n:n) /= "/") refdir(nref)(n + 1:n + 1) = "/"
    else
      nref = nref + 1
      refdir(nref) = "."//"/"
    end if
    m = 0
    do i = 1, nref
      n = len_trim(refdir(i))
      if (n /= 2) then
        if (refdir(i)(n:n) /= "/") refdir(i)(n + 1:n + 1) = "/"
        call add_path(refdir(i))
      else
        if (index(job_fn, "/") > 0) then
          refdir(i) = " "
          call add_path(refdir(i))
        end if
      end if
! the directory checking feature of inquire is specific to the Intel compiler, remove for now
!      inquire (directory = refdir(i), exist = exists)
!      if (.not. exists) then
!        write(ifiles_8,'(//10x,a)')"Reference folder """//trim(refdir(i))//""" does not exist"
!        do n =  len_trim(refdir(i)) - 1, 1, -1
!          if (refdir(i)(n:n) == "\") then
!            inquire (directory = refdir(i)(:n), exist = exists)
!            if (.not. exists) then
!              write(ifiles_8,'(10x,a)')"Reference folder """//refdir(i)(:n)//""" also does not exist"
!            else
!              write(ifiles_8,'(10x,a)')"But reference folder """//refdir(i)(:n)//""" does exist"
!              exit
!            end if
!          end if
!        end do
!        m = 1
!      end if
    end do
 !   call finish
    if (m == 1) call finish
    i = Index (jobnam, " ") - 1
    inquire (file=jobnam(:i)//".err", exist=exists)
    if (exists) then
      open (unit=iw, file=jobnam(:i)//".err", iostat = j)
      if (j == 0) close(iw, status='DELETE')
    end if
    inquire (file="fort.26", exist= exists)
    if (exists) then
      open (unit=iw,  file="fort.26", iostat = i)
      if (i == 0) close(iw, status='DELETE')
    end if
    open (unit=iw, status="SCRATCH")
    lnmols = 0
    moperr = .false.
    do nmols = 1, nmolmx
      loop1: do
!
! This is a bit clumsy: loop is called twice, the first time ALL systems
! are read in, the second time, only valid systems are used.
!
        do loop = 1, 2
          nnmols = 0
          ifiles_8 = jfiles(loop+2)
!
!  Open channel to reference data file
!
          if (names(nmols) (1:1) == " ") exit
          do iref = 1, nref
            refnam = trim(refdir(iref)) // trim(names (nmols))
            if (Index(refnam,".mop") == 0) refnam = trim(refnam)//".mop"
            inquire (unit=igpt, opened=opend)
            if (opend) then
              close (unit=igpt, status="KEEP")
            end if
            i = 0
98          i = i + 1
            open (igpt, status="OLD", file=trim(refnam), blank="ZERO", iostat=l)
            if (l /= 0) then
              continue
            end if
            if ( l == 30 .and. i < 10) then
!
!  The file exists, but is not currently accessible
!
              write(ifiles_8,*) "Problem with", refnam
              call sleep(3)
              goto 98
            end if
            if (l == 0) exit
          end do
          inquire (unit=ir, opened=opend)
          if ( .not. opend)  open(ir,status='scratch',blank='zero')
          rewind (ir)
          do
            read (igpt, "(A360)", end=1300, err=1300) line
            if (line(1:1) /= "*") then
              do j = max(3,len_trim(line)), 2, -1
                if (line(j:j) /= " ") exit
              end do
              if (Ichar (line(j:j)) == 13) then
                j = j - 1
              end if
              if (j == 0) then
                write (ir,*) " File ", names (nmols) (:i), " is_loc empty!"
                l = 0
                exit
              else
                write (ir, "(A)") line (:j)
              end if
            end if
          end do
1300      line = " "
          write (ir, "(A)") " "
          rewind (ir)
          if (l /= 0) then
            write (ifiles_8,'(a)') "  ++++ NOT FOUND ++++ %: """//trim(names(nmols))//'"'
            nnmols = nnmols + 1
          else
            rewind (2)
!
!   Reset NASIZE(8) so that MOPAC will think that the geometry
!   has not been read in.  This allows faulty data sets to be handled
!   correctly.
!
             natoms = 0
!
!   Read in all the data in the reference data file
!
            i = itype
            id = 0
            call readmo ()
              if (nmols == 223) then
               itype = i
              end if
            itype = i
!***********************************************************************
!
!  Put reference data into upper case EXCEPT between "<" and ">"
!
            case = .true.
            do i = 1, len_trim(title)
              iline = Ichar (title (i:i))
              if (title (i:i) == "<") case = .false.
              if (title (i:i) == ">") case = .true.
              if (case .and. iline >= ilowa .and. iline <= ilowz) title (i:i) = Char(iline + icapa - ilowa)
            end do
!***********************************************************************
!
!  Delete any keywords that PARAM might set
!
            call l_control("PM3", len("PM3"), -1)
            call l_control("PM7", len("PM7"), -1)
            call l_control("PM6-DH", len("PM6-DH"), -1)
            call l_control("PM6-D3H4", len("PM6-D3H4"), -1)
            call l_control("PM6-D3H4X", len("PM6-D3H4X"), -1)
            call l_control("PM6-D3(H4)", len("PM6-D3(H4)"), -1)
            call l_control("PM6", len("PM6"), -1)
            call l_control("PM7A", len("PM7A"), -1)
            call l_control("NEW", len("NEW"), -1)
            call l_control("AM1", len("AM1"), -1)
            call l_control("1SCF", len("1SCF"), -1)
            call l_control("1SCF", len("1SCF"), -1)
            call l_control("1SCF", len("1SCF"), -1)
            call l_control("EXTERNAL", len("EXTERNAL"), -1)
            if (index(contrl," LET") /= 0) then
              i = index(keywrd,"        ")
              keywrd(i:i+7) = " GEO-OK"
            end if
            i = index(keywrd," GNORM")
            if ( i == 0) then
              i = Index (contrl, "GNORM")
              if (i /= 0) then
                k = index(contrl(i:)," ")
                j = index(keywrd,"            ")
                keywrd(j + 1:j+7) = contrl(i:i + k)
              end if
            end if
!
! Write in method name
!
            i = Index (keywrd, "        ") + 1
            if (method_mndo)        keywrd(i:i+3) = "MNDO"
            if (method_am1 )        keywrd(i:i+3) = "AM1 "
            if (method_pm3 )        keywrd(i:i+3) = "PM3 "
            if (method_pm5 )        keywrd(i:i+3) = "PM5 "
            if (method_mndod)       keywrd(i:i+4) = "MNDOD"
            if (method_PM6_DH_plus) keywrd(i:i+6) = "PM6-DH+"
            if (method_PM6_D3H4)    keywrd(i:i+7) = "PM6-D3H4"
            if (method_PM6_D3H4X)   keywrd(i:i+8) = "PM6-D3H4X"
            if (method_PM6_D3_not_H4) keywrd(i:i+9) = "PM6-D3(H4)"
            if (method_PM6_DH2)     then
              if (index(contrl, " PM6-DH2X") /= 0) then
               call l_control("PM6-DH2X", len("PM6-DH2X"),  1)
              else
                call l_control("PM6-DH2", len("PM6-DH2"),  1)
              end if
            end if

            if (method_pm7)    then
              if (index(keywrd,"PM7") == 0) then
                  keywrd(i:i+4) = "PM7 "
                end if
            end if
            call l_control("PRECISE", len("PRECISE"), -1)
            if(.not. large) then
              call l_control("HCORE", len("HCORE"), -1)
              call l_control("FOCK", len("FOCK"), -1)
              call l_control("PL", len("PL"), -1)
              call l_control("VECTORS", len("VECTORS"), -1)
              call l_control("COMPFG", len("COMPFG"), -1)
              call l_control("OLDENS", len("OLDENS"), -1)
              call l_control("NOLOG", len("NOLOG"), -1)
              i = Index (title, " TYPE")!  Start of title scrub
              if (i /= 0) then
                j = index(title(i + 1:)," ") + i + 1
                title (i:j) = " "
              end if
              if (cp) call l_control("THERMO(298,298,1)", len("THERMO(298,298,1)"), 1)
            end if
!***********************************************************************
            if (clean) call l_control("MECI", len("MECI"), -1)
!***********************************************************************
            do
              i = index(keywrd(:len_trim(keywrd)),"  ")
              if (i /= 0) then
                keywrd(i:) = keywrd(i + 1:)
              else
                exit
              end if
            end do
            do
              i = index(title(:len_trim(title)),"  ")
              if (i /= 0) then
                title(i:) = title(i + 1:)
              else
                exit
              end if
            end do
!
!  Set PARAM keywords to be used by the MOPAC side
!
            if( precise) call l_control("PRECISE", len("PRECISE"), 1)
            if (natoms == 0 .or. moperr) then
              write (ifiles_8,'(a,a)') " ++++ FILE FAULTY ++++ % ",  trim(refnam)
              write (ifiles_8, "(3A)") " (DATINP) There is an error in data-set ", names (nmols)
              if (moperr) then
                write (ifiles_8, "(2A)") " Error message:", errtxt (1:50)
              end if
              moperr = .false.
              faulty = faulty + 1
              nnmols = nnmols + 1
            else
              if (numvar /= 0 .and. .not. all) then
                if (natoms > max_atom) go to 1400
!
!   CHECK IF COMPOUND CONTAINS AN ELEMENT TO BE PARAMETERIZED.
!
                do i = 1, natoms
                  if (elemok(1, labels(i)) .and. labels(i) /= 99 .and. labels(i) /= 107) exit
                end do
                if(i > natoms) then
                  ne_excluded = ne_excluded + 1
                  go to 1400
                end if
!
!   CHECK IF ONLY THE ELEMENTS BEING OPTIMIZED ARE PRESENT
!
                if (l_only) then
                  do i = 1, natoms
                    if (.not. elemok(1, labels(i))) exit
                  end do
                  if(i <= natoms) go to 1400
                end if
              end if
!
!   CHECK IF ALL ELEMENTS ARE ALLOWED
!
              do i = 1, natoms
                if ( .not. elemok(2, labels(i))) go to 1400
              end do
              do i = len_trim(refnam), 2, -1
                if (refnam(i:i) == backslash) exit
              end do
              j = 0
              do i = 1, natoms
                if (labels(i) /= 99 .and. labels(i) /= 107) then
                  j = j + 1
                  nat(j) = labels(i)
                end if
              end do
              if (.not.(let .or. core_core_OK(geo, nat, numat, k, l, ccp) .or. .not. method_pm6)) then
                do i = len_trim(refnam), 1, -1
                  if (refnam(i:i) == backslash .or. refnam(i:i) == "/") exit
                end do
                write(ifiles_8,"(//,a, /, a,//)")"++++ FAULT DETECTED IN ++++ "//trim(koment)//", in file: """ &
                & //refnam(i + 1:len_trim(refnam)-4)//"""","     File is missing the "//trim(elemnt(k))//" - "//trim(elemnt(l))// &
                & " core-core term. System deleted."
                go to 1400
              end if
              if (.not. all .and. Index(title,"GUESS") /= 0 .and. survey ) goto 1400
              if (Index(title,"RCJ02") /= 0 .and. survey ) goto 1400
!
!  Store reference data
!
!
!
              hr = " "
              dr = " "
              jr = " "
              gr = " "
              if (title (1:1) /= " ") then
                name = title (1:80)
                title (1:1) = " "
                title (2:) = name(1:79)
              end if
              heats(nmols) = 0.d0
              dipls(nmols) = 0.d0
              hions(nmols) = 0.d0
              weight(1:3, nmols) = 0.d0
              weight(4, nmols) = wtgeo
              weight(5,nmols) = 1.d0
              mol = 1
              hof_type = .false.
!
!  If the system is NOT to be used in constructing tables, set is_a_ref .true.
!
              is_a_ref(nmols) = (Index(title," HR=REF") /= 0 &
              & .or. (Index(title," HB ") /= 0 .and. survey) )
              refers (nmols, 1) = " "
              if (cp) then
                i = Index (title, " CP=")
              else
                i = Index (title, " H=")
              end if
              if (i /= 0) then
                hof_type = .true.
                j = Index (title, " HR=")
                if (j /= 0) then
                  refers(nmols, 1) = title (j+4:j+12)
                  j = Index (refers(nmols, 1), " ")
                  if (j /= 0) then
                    refers (nmols, 1) (j:) = " "
                  end if
                  if (is_a_ref(nmols) .and. survey) refers(nmols,1) = " "
                else
                  mish = mish + 1
                  hr = "*"
                end if
                heats(nmols) = reada (title, i)
!
!  Read in the excited state quantum numbers, if they exist
!
                i = Index(title, " ROOT=")+5
                qn(nmols) = 0
                Excited_state: if( i /= 5 ) then
!
! Format:  "ROOT=n,m,rep", n   = quantum number,
!                          m   = spin multiplicity,
!                          rep = representation.
                  j = Index(title(i:), " ")+i
!
!  All the text is in title(i:j)
!
                  k = Index(title(i:j), ",")
                  l = Index(title(i+k:j), ",")
                  if(k /= 0 .and. l /= 0) then
                  k = k + i
                  l = l + k
                  qn(nmols) = Nint(Reada(title(i:k),1))
                  spin_state(nmols) = Nint(Reada(title(k:l),1))
                  i_r(nmols) = title(l:j-1)
                  do i = 2, 4
                      iline = Ichar (i_r(nmols)(i:i))
                      if (iline >= icapa .and. iline <= icapz) then
                        i_r(nmols)(i:i) = Char(iline+ilowa-icapa)
                      end if
                  end do
                end if
!
!  If keyword MECI is not present, put it in - to force MECI to calculate all states
!
                if(Index(keywrd," MECI") == 0) then
                  i = Index(keywrd,"      ")
                  keywrd(i+1:i+4) = "MECI"
                end if
              end if Excited_state
              do
!
!  LOAD IN THE STANDARD DEVIATION, IF KNOWN
!
!#            J=INDEX(TITLE(I:),' ')
!#            K=INDEX(TITLE(I:J+I),',')
!#            HSDS(NMOLS)=0.D0
!#            HSDS(NMOLS)=1.D0
!#            IF(K.NE.0)HSDS(NMOLS)=READA(TITLE,I+K)
                i = Index (title, "+")
                if (i == 0) exit
                if (wth < 1.d-3) exit
!
!   HEAT OF FORMATION (EXPT) IS A FUNCTION OF ANOTHER MOLECULE.
!
                if (nmols == 1) then
                  write(ifiles_8,*)" The first molecule cannot depend on other molecules"
                  faulty = faulty + 2
                  deadly = .true.
                end if
                title (i:i) = "^"
!
!  If name has quotation marks, use these as delimiter
!
                if (title(i+1:i+1) == """") then
!
!  j = end of name
!
                 j = Index (title (i + 2:), """") + i
                 if (j == i-1) then
                    j = 80
                  end if
                  l = j - i
                  name = title (i+2:j)
!
! Remove any "+" signed from the name in title
!
                  do
                    k = Index(title (i+2:j),"+")
                    if (k == 0) exit
                    title (i+k+1:i+k+1)="^"
                  end do
                else
                  j = Index (title (i:80), " ") + i - 1
                  if (j == i-1) then
                    j = 80
                  end if
                  k = Index (title (i+1:80), "+")
                  if (k /= 0) then
                    k = Index (title (i+1:80), "+") + i - 1
                  else
                    k = Index (title (i+1:80), " ") + i - 1
                  end if
                  if (k == i-1) then
                    k = 80
                  end if
                  if (k < j) then
                    j = k
                  end if
                  l = j - i
                  name = title (i+1:j)
                end if
                do k = 1, len_trim(name)
                  iline = Ichar (name(k:k))
                  if (iline >= icapa .and. iline <= Ichar ("Z")) then
                    name(k:k) = Char(iline+ilowa-icapa)
                  end if
                end do
                refok = .true.
                k = len_trim(name)
                if (name(k - 3:) == ".mop") name(k - 3:) = " "
                do k = 1, nmols - 1
                  line = names(k)
                  do j = 1, len_trim(line)
                    iline = Ichar (line(j:j))
                    if (iline >= icapa .and. iline <= Ichar ("Z")) line(j:j) = Char(iline+ilowa-icapa)
                  end do
                  if (name == line .and. len_trim(name) == len_trim(line)) go to 1500
                end do
                write (ifiles_8,"(/,2a)") " ++++ FAULT DETECTED IN FILE ++++ %: ", trim(refnam)
                write (ifiles_8, "(' NAME OF REFERENCE MOLECULE <',A, '> NOT KNOWN')") &
                trim(name)//".mop"
                do k = 1,is
                if(safety(k)(1:l+1) == name(1:l+1)) exit
                end do
                if (k > is) then
                  is = is + 1
                  safety(is) = name(:80)
                end if
                write (ifiles_8, "(//,A,//)") " (DATINP) There is an error in data-set ", &
                names (nmols)
                refok = .false.
                mmmols = mmmols + 1
                cycle
1500            mol = mol + 1
                ihrefs(mol, nmols) = k
!
!  Check that the dependent species do not depend on anything
!
                if (ihrefs(1,k) > 0) then
                  faulty = faulty + 2
                  deadly = .true.
                  write (ifiles_8,"(/,2a)") " ++++ FAULT DETECTED IN FILE ++++ %: ", trim(refnam)
                  write (ifiles_8, "(' Dependent molecule ""', &
                & a, &
                &'"" depends on other molecules.',/,' This is not allowed.',/)") &
                & names(ihrefs(mol, nmols))(1:len_trim(names(ihrefs(mol, nmols))))
                end if
              end do
              ihrefs(1, nmols) = mol - 1
              if (refok .and. mol /= 1) then
                write (ifiles_8, "(' TO THE REF. H.O.F MUST BE ADDED',' THE CALC. HEAT OF ',10a)") &
               &(""""//names(ihrefs(j, nmols))(1:len_trim(names(ihrefs(j, nmols))))//""" ", j=2, mol)
              end if
              weight(1, nmols) = wth
              i = Index (title, " HWT=")
              if (i /= 0) weight(1, nmols) = reada (title, i) * wth
!
! Set weights for states (defined by the presence of "ROOT=") to be less than other heats of formation, because the
! emphasis is on normal heats of formation, and state energies are used just to help tame the PES.
!
              if (qn(nmols) /= 0) weight(1, nmols) = 0.5d0*weight(1, nmols)
              if (is_a_ref(nmols)) weight(1, nmols) = 1.d-6
            else
              if (is_a_ref(nmols)) then
                write (ifiles_8,"(/,2a)") " ++++ FAULT DETECTED IN FILE ++++ %: ", trim(refnam)
                write(ifiles_8,"(a)")" If ""HR=REF"" is present, a heat of formation MUST be specified"
                faulty = faulty + 2
              end if
            end if
            do i = 1, len_trim(title)
              if (title (i:i) == "^") then
                title (i:i) = "+"
              end if
            end do
!
!  pKa
!
            i = Index (title, " PKA")
            if (i /= 0) then
              j = Index (title, " PKR=")
              if (j /= 0) then
                refers(nmols, 6) = title (j+4:j+12)
                j = Index (refers(nmols, 6), " ")
                if (j /= 0) then
                  refers (nmols, 6) (j:) = " "
                end if
              else
                mish = mish + 1
                hr = "*"
              end if
!
!  pKa value is stored
!
              Atom_pKas(nmols,1) = Nint(Reada(title,i))!  i
              i = i + Index(title(i:),"=")
              pKas(nmols) = reada (title, i)
              heats(nmols) = pKas(nmols)
              sum = 10.d0!   pKa default weight relative to H.o.F. weight
              weight(1, nmols) = wth * sum
              i = Index (title, " PKA_WT=")
              if (i /= 0) then
                weight(1, nmols) = reada (title, i) * wth * sum
              end if
            else
                pkas(nmols) = 0.d0
            end if
            refers (nmols, 3) = " "
            if (cp) then
              i = Index (title, " S=")
            else
              i = Index (title, " D=")
            end if
            if (i /= 0) then
              hof_type = .true.
              j = Index (title, " DR=")
              if (j /= 0) then
                refers(nmols, 3) = title (j+4:j+12)
                j = Index (refers(nmols, 3), " ")
                if (j /= 0) then
                  refers (nmols, 3) (j:) = " "
                end if
              else
                misd = misd + 1
                dr = "*"
              end if
              dipls(nmols) = reada (title, i)
!
!  LOAD IN THE STANDARD DEVIATION, IF KNOWN
!
!#            J=INDEX(TITLE(I:),' ')
!#            K=INDEX(TITLE(I:J+I),',')
!#            DSDS(NMOLS)=0.D0
!#            DSDS(NMOLS)=0.2D0
!#            IF(K.NE.0)DSDS(NMOLS)=READA(TITLE,I+K)
              weight(2, nmols) = wtd
              i = Index (title, " DWT=")
              if (i /= 0) then
                weight(2, nmols) = reada (title, i) * wtd
              end if
            end if
            refers (nmols, 2) = " "
            i = Index (title, " I=") + Index (title, " IA") + Index (title, " IE") + &
            & Index (title, " IP")+ Index (title, " I1=")+ Index (title, " I2=") + &
            & Index (title, " I3=")+ Index (title, " I4=")
            if (i /= 0 .and. Index(title,"IR=PW91D") == 0) then
              hof_type = .true.
              if (title (i+2:i+2) == "=") then
                lions(nmols) = 0
              else
                j = Ichar (title (i+2:i+2)) - Ichar ("0") - 1
                if (j > 9 .or. j < 0) then
                  write (ifiles_8, "(A,I4,A)") " FOR MOL.", nmols, " I.P. BAD"
                  lions(nmols) = 0
                else
                  lions(nmols) = j
                end if
              end if
              i = i + 3
              hions(nmols) = reada (title, i)
              weight(3, nmols) = wti
              j = Index (title, " IR=")
              if (j /= 0) then
                refers(nmols, 2) = title (j+4:j+12)
                j = Index (refers(nmols, 2), " ")
                if (j /= 0) then
                  refers (nmols, 2) (j:) = " "
                end if
              else
                jr = "*"
                misi = misi + 1
              end if
              i = Index (title, " IWT=")
              if (i /= 0) then
                weight(3, nmols) = reada (title, i) * wti
              end if
            end if
            i = Index (title, " GWT=")
            if (i /= 0) then
              weight(4, nmols) = reada (title, i) * weight(4, nmols)
            end if
            refers (nmols, 4) = " "
            if (Index (title, "GEOREF")+Index (title, "GR=") == 0) then
              weight(4, nmols) = 0.d0
              do i = 1, 300
                geotxt (i, nmols) = " "
              end do
            else
              j = Index (title, "GR=")
              if (j /= 0) then
                if (hof_type ) then
!
!  PANIC
!
                write(ifiles_8,"(3a)")' Error: Molecule "',names(nmols)(1:j), &
                    &'" is both a geometric and a non-geometric reference.', &
                    ' Modify reference data-set and re-run'
                    faulty = faulty + 1
                end if
                refers(nmols, 4) = title (j+3:j+12)
                j = Index (refers(nmols, 4), " ")
                if (j /= 0) then
                  refers (nmols, 4) (j:) = " "
                end if
              else
                gr = "*"
                misg = misg + 1
              end if
!
!   READ IN NAMES OF GEOMETRIC PARAMETERS
!
              l = 0
              k = 1
              do i = 1, 300
                j = Index (title (k:), "<")
                if (j /= 0) then
                  l = Index (title (k:), ">")
                  if (l == 0) then
                    write(ifiles_8,"(3a,i3)")' Error: Molecule "',trim(names(nmols)), &
                      '" has unterminated geometric reference.'
                    faulty = faulty + 1
                    deadly = .true.
                  else if (l == j+1) then
                    geotxt (i, nmols) = " "
                  else
                    geotxt(i, nmols) = title (k+j:k+l-2)
                    if (i > nvar) then
                      write(ifiles_8,"(3a,i3)")' Error: Molecule "',trim(names(nmols)), &
                      ' has more geometric reference data than variables:',nvar
                      faulty = faulty + 1
                      deadly = .true.
                    end if
                  end if
                else
                  l = 0
                  geotxt (i, nmols) = " "
                end if
                k = k + l
              end do
!
!  Check: are the geometric references in internal coordinates?
!
              j = min(40,nvar)
              do i = 1, j
                if (geotxt(i,nmols) /= " ") then
                  if(na(loc(1,i)) == 0) then
                    write(ifiles_8,"(3a)")' Error: Molecule "',trim(names(nmols)), &
                    &'" is a geometric reference, but geometry is not ',&
                    &' entirely in internal coordinates. Modify reference data-set and re-run'
                    faulty = faulty + 1
                  else
                    if (loc(2,i) == 1) then
!
!  Make a note of bond-lengths for advice to users
!
                      bonds_ij(labels(loc(1,i)), labels(na(loc(1,i)))) = &
                      & bonds_ij(labels(loc(1,i)), labels(na(loc(1,i)))) + 1
                    end if
                  end if
                end if
              end do
            end if
!
!  Add in factor for number of formula units in the cluster, if present
!
            if (id > 0) then! Special handling for solids - divide HoF by number of mers.
              j = Index (keywrd, " MERS")
              if (j /= 0) then
                mers = 0
                k = 0
                i = Index (keywrd(j + 1:), " ") + j
                do l = 1, 3
                  j = j + k
                  if (l > 1 .and. k == 0) exit
                  mers(l) = Nint (reada (keywrd(j:), 1))
                  k = Index (keywrd(j:i), ",")
                end do
              end if
              i = Index(keywrd, " Z=")
              if (i /= 0) then
                i = Nint(reada(keywrd,i))
                i = mers(1)*mers(2)*mers(3)*i
                if (Index (keywrd, " BCC") /= 0) i = i/2
                weight(5,nmols) =1.d0/i
              else
                nele = 0
                do i = 1, natoms - 3
                  nele(labels(i)) = nele(labels(i)) + 1
                end do
                j = 0
                do i = 1, 98
                  if (nele(i) > 0) then
                    j = j + 1
                    nele(j) = nele(i)
                  end if
                end do
                k = 1000
                do i = 1, j
                  if(nele(i) < k) k = nele(i)
                end do
!
!  k is the smallest number of atoms of any element in the formula
!
                do i = 1, 10
                  m = 0
                  do l = 1, j
                    if (Abs((i*nele(l))/k - (i*1.d0*nele(l))/k) > 1.d-5) m = 1
                  end do
                  if (m == 0) exit
                end do
!
!  Number of empirical units  = k/i
!
                weight(5,nmols) =(i*1.d0)/k
              end if
            end if
!
! IF GRADIENT OR FULL OPTIMISATION CALCULATION DONE THEN
!    WEIGHT(4,NMOLS) DEFAULTS TO 1.D0, UNLESS "NOGEO" IS SPECIFIED,
!    WHEN IT DEFAULTS TO 0.D0.
! IF GWT=N.NNN IS SPECIFIED, THIS CAN BE OVER-RIDDEN BY "NOGEO".
!
            tnow = seconds (1)
            itime = Int (tnow-told)
            told = tnow
            if(koment == " ")then
            i=Index(names(nmols)," ")
            write(ifiles_8,"(3a)")"Comment for ",names (nmols)(1:i)," is empty."
            call finish
            end if
            do
              if(koment(1:1) /= " ") exit
              koment = koment(2:)
            end do
! Write out heats of formation, only, without spaces
!
!         if (heats (nmols) < -999.9950d0) then
!          write(55,"(f9.2)",err=95) heats (nmols)
!        end if
            endfile (ifiles_8)
            backspace (ifiles_8)
            if (moperr) then
              write (ifiles_8, "(//, A)") " (READMO) There is an error in data-set", &
              names (nmols)
              write (ifiles_8, "(2A,//)") " Error message:", errtxt (1:50)
              moperr = .false.
              faulty = faulty + 1
              go to 1600
            else
!
!  ATMASS, REACT, AND TXTATM ARE NOT USED, THEREFORE USE DUMMYS
!
!
!   COPY1 DEFINES NB AND NC IN TERMS OF NA AND NATOMS.  THEREFORE,
!   NATOMS MUST BE RESET TO THE MAXIMUM (DEFINED IN MAKECM).
!
              if (moperr) then
                write (ifiles_8, "(//,3A)") " (STATE)  There is an error in data-set ", names (nmols)
                write (ifiles_8, "(2A,//)") " Error message:", trim(errtxt)
                moperr = .false.
                faulty = faulty + 1
                go to 1600
              else
                call moldat(2)
                if (moperr) then
                  if(errtxt(1:10) == "PARAMETERS")then
                    moperr=.false.
                    exit loop1
                  end if
                  write (ifiles_8, "(//,3A)") " ++++ FAULT DETECTED IN FILE ++++ (MOLDAT)&
                                               & There is an error in data-set ", names (nmols)
                  write(ifiles_8, "(a,//)")" Error message: """//trim(errtxt)//""". System deleted."
                  moperr = .false.
                  go to 1600
                else
!
!  Now that "molele" is filled, work out the RELATIVE weights of the reference data.
!
                  if (is_a_ref(nmols)) exit loop1
                  if (weight(1,nmols) + weight(2,nmols) + weight(3,nmols) + weight(4,nmols) < 1.d-4) then
                    if (weight(1,nmols) + weight(2,nmols) + weight(3,nmols) + weight(4,nmols) < 1.d-10) then
                      faulty = faulty + 2
                      deadly = .true.
                      write (ifiles_8,"(/,3a,/)") " ++++ FILE FAULTY ++++ %: ", trim(refnam), &
                        & ': No reference data '
                      end if
                      goto 1400
                    end if
                  exit loop1
                end if
              end if
            end if
!
!  No elements to be parameterized.
!
!  THIS COMPOUND CANNOT BE RUN - ELIMINATE IT!
!
1400        natoms = 0
          end if
        end if
1600    if (loop == 2) then
          notok = notok + 1
          do i = nmols, nmolmx - notok
            names(i) = names(i+1)
          end do
          names (nmolmx-notok+1) = " "
          if (nmols <= nmolmx-notok) cycle loop1
        end if
      end do
      go to 1700
    end do loop1
!
!  CALCULATE THE M.O. INDEX FOR I.P.
!
    lions(nmols) = nelecs / 2 - lions(nmols)
    if (nelecs == 0) then
      weight(3, nmols) = 0.d0
    end if
!
!  Identify all elements in the system
!
    nel = 1
    molele(1,nmols) = 0
    loop2: do i = 1, natoms
      if (labels(i) == 99) cycle loop2
      do j = nel, 2, -1
        if (molele(j, nmols) == labels(i)) cycle loop2
      end do
      nel = nel + 1
      molele(nel, nmols) = labels(i)
    end do loop2

    sum = 1.d0
    do j = 2, nel
      sum = min(sum, element_weights(molele(j,nmols)))
    end do
    weight(1:4,nmols) = weight(1:4,nmols)*sum
    if (natoms > 0) then
      write (ifiles_8, "(1X,A50,I4, F8.2,A,F6.2,F6.3,A,F5.1,F7.2,A,F5.1,&
           &F7.2,A)") koment, itime, heats (nmols), hr, weight (1, nmols), &
           & dipls(nmols), dr, weight(2, nmols), hions(nmols), jr, &
           & (weight(i, nmols), i=3, 4), gr
!
      lnmols = lnmols + 1
      is = is + 1
      safety(is) = names(nmols)
      if (Index (keywrd, " XFAC") /= 0 .and. index(keywrd, " POP") == 0 .and. numat /= 2) then
        write(ifiles_8,"(a)")" Number of atoms in reference function must be 2"
        call finish
      end if
!
!  element 99 is used for all parameters that are constant for all elements
!
      nel = nel + 1
      molele(nel,nmols) = 99
      molele(1, nmols) = nel
!
! FIRST, STORE ALL DATA DEFINING THE SIZE OF THE MOLECULE.
!
      natmss(lnmols) = natoms
      numats(lnmols) = numat
      ndeps(lnmols) = ndep
      norbss(lnmols) = norbs
      n2elecs(lnmols) = n2elec
      nlecss(lnmols) = nelecs
      ncloss(lnmols) = nclose+nalpha
      nopens(lnmols) = nopen+nbeta
      nnalpha_open(lnmols) = nalpha_open
      nnbeta_open(lnmols) = nbeta_open
      keys(lnmols) = trim(keywrd)
      msdels(lnmols) = msdel
      titls(lnmols) = title
      nvars(lnmols) = nvar
      comnts(lnmols) = koment(1:80)
      nnmos(lnmols) = nmos
      nlm61(lnmols) = lm61
      fracts(lnmols) = fract
!
! NOW TO STORE GEOMETRIC DATA
!
      do i = 1, nvar
        iloc = iloc + 1
        locs(1, iloc) = loc(1, i)
        locs(2, iloc) = loc(2, i)
      end do
      l123s(1,lnmols) = l1u
      l123s(2,lnmols) = l2u
      l123s(3,lnmols) = l3u
      do i = 1, natoms
        igeo = igeo + 1
        do j = 1, 3
          geos(j, igeo) = geo(j, i)
        end do
      end do
      do i = 1, ndep
        isym = isym + 1
        locpas(isym) = locpar(i)
        idepfs(isym) = idepfn(i)
        locdes(isym) = locdep(i)
        depmuls(isym) = depmul(i)
      end do
!
!  NOW TO STORE ATOMIC DATA.
!
      do i = 1, natoms
        iatm = iatm + 1
        nas(iatm) = na(i)
        nbs(iatm) = nb(i)
        ncs(iatm) = nc(i)
        nats(iatm) = nat(i)
        nfirss(iatm) = nfirst(i)
        nlasts(iatm) = nlast(i)
        lablss(iatm) = labels(i)
      end do
!
!  FINALLY, STORE ORBITAL DATA.
!
      if (nalpha /= nbeta) then
        sum = dfloat(nalpha)/nelecs
      else
       sum = 0.5d0
      end if
      do i = 1, norbs
        pas(ilin+ (i*(i+1))/2) = pdiag(i)*sum
        pbs(ilin+ (i*(i+1))/2) = pdiag(i)*(1.d0-sum)
      end do
      ilin = ilin + (norbs*(norbs+1)) / 2
      if (Abs (weight(2, nmols)) > 1.d-4) then
        ndips = ndips + 1
      end if
      if (Abs (weight(3, nmols)) > 1.d-4) then
        nions = nions + 1
      end if
    end if
1700 continue
  end do
  close (unit=jfiles(1), status="DELETE")
  close (unit=igpt)
  close (unit=55)
  close (unit = 56)
  ifiles_8 = jfiles(3)
!
! IF THERE ARE NO MOLECULES LEFT AT ALL, THEN STOP
!
  if (nmolmx-notok == 0) then
    write (ifiles_8,*) " NO ACCEPTABLE MOLECULES IN SET SUPPLIED"
    if (ne_excluded /= 0) write (ifiles_8,*) " (No molecules contained elements marked for optimization)"
    write (ifiles_8,*) " (To by-pass this check, add keyword ""ALL"")"
    call finish
  end if
! iw = jfiles(1)
  ifiles_8 = jfiles(3)
  nmolmx = nmolmx - notok
  if (nnmols /= 0) then
    write (ifiles_8,*) nnmols, " MOLECULES WERE NOT FOUND"
    if ( .not. let) call finish
  end if
  write(ifiles_8,'(/)')
  if (mish+misd+misi+misg /= 0) then
    if (mish /= 0) then
      write (ifiles_8, "(10x,I4,A)") mish, " REFERENCES FOR HEATS OF FORMATION MIS&
     &SING"
    else
      write (ifiles_8, "(10x,A)") " ALL HEATS OF FORMATION REFERENCED."
    end if
    if (misd /= 0) then
      write (ifiles_8, "(10x,I4,A)") misd, " REFERENCES FOR DIPOLES MISSING"
    else
      write (ifiles_8, "(10x,A)") " ALL DIPOLES REFERENCED."
    end if
    if (misg /= 0) then
      write (ifiles_8, "(10x,I4,A)") misg, " REFERENCES FOR " // "GEOMETRIES MISSI&
     &NG"
    else
      write (ifiles_8, "(10x,A)") " ALL GEOMETRIES REFERENCED."
    end if
    if (misi /= 0) then
      write (ifiles_8, "(10x,I4,A)") misi, " REFERENCES FOR IONIZATION POTENTIALS &
     &MISSING"
    else
      write (ifiles_8, "(10x,A)") " ALL IONIZATION POTENTIALS REFERENCED."
    end if
  else
    write (ifiles_8, "(10x,A)") " REFERENCES PRESENT FOR ALL DATA."
  end if
  write(ifiles_8,*)
  if (mmmols /= 0) then
    write (ifiles_8,*) mmmols, " REFERENCE HEATS OF FORMATION NOT KNOWN"
    write (ifiles_8,"(a)")" Use the following files with ALL to allow the run to continue"
    write (ifiles_8,*)" (This file will also be written to 'bits.txt')"
    k = 15
    line = "bits.txt"
    call add_path(line)
    open (unit=k, file=trim(line))
    do i = 1, is
      do j = 80,2,-1
      if(safety(i)(j:j) /= " ") exit
      end do
      write (ifiles_8,"(a)")safety(i)(1:j)
      write (k,"(a)")safety(i)(1:j)
    end do
    if ( .not. let) call finish
  end if
  close (ir)
  nmols = lnmols
  write (ifiles_8, "(//10X,' ADVICE ON STORAGE SPACE REQUIREMENTS')")
  write (ifiles_8, "(//10X,' STORE   SPACE USED  ', ' TOTAL SPACE AVAILABLE')")
  write (ifiles_8,*)
! write (ifiles_8, "(11X,'ILOC',I10,I14)") iloc, maxpat
  write (ifiles_8, "(11X,'IGEO',I10,I14)") igeo, maxatm
  write (ifiles_8, "(11X,'ISYM',I10,I14)") isym, maxsym
  write (ifiles_8, "(11X,'ILIN',I10,I14)") ilin, maxpab
  write (ifiles_8, "(11X,'NMOL',I10,I14)") nmols, maxmol

!
!   Write out the matrix of core-core terms
!
  ccp(99,99) = 1! Deliberately allow dummy atom parameters to be optimized.
  do i = 1, 83
    do j = 1,i-1
      ccp(i,j) = ccp(i,j) + ccp(j,i)
      ccp(j,i) = ccp(i,j)
    end do
  end do
  if (i == -999) then! Dummy suppress printing of matrices
  write(ifiles_8,"(/10x,a,/)") " Matrix of core-core terms * = present, - = absent"
  do i = 1, 83
    do j = 1,i
!  call setalp_read (i,j, sum, wth)
      if (sum > 1.d-3) then
        if (ccp(i,j) > 1) then!  2 data are sufficient to define the core-core term
          line(j:j) = "*"
        else
          line(j:j) = Char(Ichar("0") + ccp(i,j))
        end if
      else
        line(j:j) = "-"
      end if
    end do
     select case (i)
    case (2,10,18,36,54)
      write(ifiles_8,*)
    end select
    write(ifiles_8,"(1x,a2,1x,2a1,1x,8a1,1x,8a1,1x,18a1,1x,18a1,1x,32a1)")elemnt(i),(line(j:j),j=1,i)
  end do

!
!   Write out the matrix of number of bond-lengths used
!
  write(ifiles_8,"(/10x,a,/)") " Matrix of number of bond-lengths used."
  do i = 1, 83
    do j = 1,i-1
      bonds_ij(i,j) = bonds_ij(i,j) + bonds_ij(j,i)
      bonds_ij(j,i) = bonds_ij(i,j)
    end do
  end do
  do i = 1, 83
    do j = 1,i
      if (bonds_ij(i,j) == 0) then
        line(j:j) = "-"
      else if (bonds_ij(i,j) < 10) then
        line(j:j) = Char(Ichar("0") + bonds_ij(i,j))
      else
        line(j:j) = "+"
      end if
    end do
    select case (i)
    case (2,10,18,36,54)
      write(ifiles_8,*)
    end select
    write(ifiles_8,"(1x,a2,1x,2a1,1x,8a1,1x,8a1,1x,18a1,1x,18a1,1x,32a1)")elemnt(i),(line(j:j),j=1,i)
  end do
end if
!
!  Work out how many reference data are present
!
  j= 0
  do i = 1,nmols
    if(weight(1,i) > 1.d-9) j = j + 1
    if(weight(2,i) > 1.d-9) j = j + 1
    if(weight(3,i) > 1.d-9) j = j + 1
    if(weight(4,i) > 1.d-9) then
      do k = 1,300
        if (geotxt(k,i) /= " ") j = j + 1
      end do
    end if
  end do
  nfns = j
!
! Remove atomic parameters that are not accessible from the set of reference data supplied
!
  elemok(1,:) = .false.
  do i = 1, 107
    do j = 1,i
      if (ccp(i,j) > 0) elemok(1,i) = .true.
    end do
  end do
  elemok(:,97:99)  = .true.
!
! At this point, no parameters have been excluded yet
!
  j = 0
  do i = 1, numvar
    k = 1
!
! 39 = ALPB_, 40 = XFAC_, 41 = PAR
!
    if(locvar(1,i) /= 39 .and. locvar(1,i) /= 40 .and. locvar(1,i) /= 41) then
      if ( .not. elemok(1,locvar(2,i))) k = 0
    end if
    if (k > 0) then
!
!  Parameter is valid
!
    j = j + 1
    locvar(1,j) = locvar(1,i)
    locvar(2,j) = locvar(2,i)
    botlim(j) = botlim(i)
    toplim(j) = toplim(i)
    valvar(j) = valvar(i)
    valold(j) = valold(i)
    end if
  end do
  numvar = j
!
!  Remove core-core terms that are not accessible from the set of reference data supplied
!
  j = 0
  do i = 1, numvar
    k = 1
    if(locvar(1,i) == n_partyp_alpb .or. locvar(1,i) == n_partyp_alpb + 1) then
      nj = locvar(2,i)/200
      ni = locvar(2,i) - nj*200
      if (ccp(ni,nj) ==  0) then
        k = 0
      end if
    end if
    if (k > 0) then
!
!  Parameter is valid
!
      j = j + 1
      locvar(1,j) = locvar(1,i)
      locvar(2,j) = locvar(2,i)
      botlim(j) = botlim(i)
      toplim(j) = toplim(i)
      valvar(j) = valvar(i)
      valold(j) = valold(i)
    end if
  end do
  numvar = j
  write(ifiles_8,"(a,i5)")" Number of reference data:",nfns
  if (index(contrl, " CYCLES=0") /= 0) numvar = 0
  write(ifiles_8,"(a,i5)")" Number of parameters:    ",numvar
!
!   nfns should be correct, but add 10 "for safety" - the calculation of nfns is buggy
!
  allocate (diffns(numvar,nfns+10))
  if (iloc > maxpat .or. igeo > maxatm .or. isym > maxpat .or. ilin > maxpab &
 & .or. nmols > maxmol) then
    write (ifiles_8, "('  ARRAY BOUND ERROR - ', 'MODIFY SIZES AND RECOMPILE')&
   &")
    write (ifiles_8, "('  CURRENT VALUES')")
    write (ifiles_8, "(' MAXATM:',I8)") maxatm
    write (ifiles_8, "(' MAXSYM:',I8)") maxsym
    write (ifiles_8, "(' MAXPAB:',I8)") maxpab
    write (ifiles_8, "(' MAXMOL:',I8)") maxmol
    call finish
  end if
  inquire (unit=ir, opened=opend)
  if (opend) then
    close (unit=ir, status="DELETE")
  end if
  inquire (unit=ilog, opened=opend)
  if (opend) then
    close (unit=ilog, status="DELETE",err = 97)
  end if
97 if (faulty /= 0) then
    write (ifiles_8,'(a,i3,a)') " Faults were detected in", (faulty+1)/2, " data sets."
    if (deadly) then
      write (ifiles_8,'(a,i3,a)') " Faults MUST be corrected before job can be run."
      call finish
    end if
    if ( .not. let) then
      write(ifiles_8,*)" To continue, add ""LET"" to the keyword line."
      call finish
    end if
  end if
  nnhco = 0
  return
end subroutine datinp
logical function core_core_OK(geo, nat, numat, ni, nj, ccp)
  use parameters_C, only : alpb
  integer :: numat, ni, nj
  integer, dimension (numat) :: nat
  integer, dimension (107,107) :: ccp
  double precision, dimension (3,*) :: geo
  logical :: first = .true.
  logical, save, dimension (107,107) :: cc_OK
  integer :: i, j
  double precision, dimension (3,600) :: coord
  if (first) then
    first = .false.
!
!   Determine which core-core terms are valid
!
    cc_OK = .true.
    do i = 1, 100
      do j = 1, i
        cc_OK(i,j) = (alpb(i,j) > 1.d-4)
        cc_OK(j,i) = cc_OK(i,j)
      end do
    end do
    cc_OK(:,98:) = .true.
    cc_OK(98:,:) = .true.
    cc_OK(:,57:71) = .true.
    cc_OK(57:71,:) = .true.
  end if
!
!  Find all pairs of atoms that are near to each other (less than 5 Angstroms)
!
  call gmetry (geo, coord)
  do i = 2, numat
    ni = nat(i)
    do j = 1, i - 1
      if((coord(1,i)-coord(1,j))**2 + &
    &    (coord(2,i)-coord(2,j))**2 + &
    &    (coord(3,i)-coord(3,j))**2 < 25.5d0) then
        ccp(ni,nat(j)) =  ccp(ni,nat(j)) + 1
        if(.not. cc_OK(ni, nat(j))) then
          core_core_OK = .false.
          nj = nat(j)
          return
        end if
      end if
    end do
  end do
  do i = 1, numat
    ccp(nat(i), nat(i)) = ccp(nat(i), nat(i)) + 1
  end do
  ccp(:,98:) = 1
  ccp(98:,:) = 1
  core_core_OK = .true.
end  function core_core_OK
