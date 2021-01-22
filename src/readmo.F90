      subroutine readmo 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!
      use chanel_C, only : iw, ir, iarc, ilog, log_fn, archive_fn, &
      ires, restart_fn, output_fn, job_fn
!
      USE maps_C, ONLY: latom, lparam, lpara1, latom1, lpara2, latom2, &
      react 
!
      USE symmetry_C, ONLY: idepfn, locdep, depmul, locpar 
!
      use molkst_C, only : ndep, numat, numcal, natoms, nvar, keywrd, dh, &
      & verson, is_PARAM, line, nl_atoms, &
      & moperr, maxatoms, koment, title, method_pm6, refkey, &
      isok, ijulian, gui, Academic, site_no, method_pm6_dh2, caltyp, &
      method_pm7, jobnam, method_PM7_ts, arc_hof_1, keywrd_txt, txtmax, refkey_ref, &
      ncomments, itemp_1, nbreaks, numat_old, maxtxt, num_bits, use_ref_geo, &
      n_methods, methods, methods_keys,  method_pm6_d3h4, method_pm6_dh2x,   &   
      method_pm6_d3h4x, method_pm6_d3, method_pm6_d3_not_h4, method_pm7_hh, &
      method_pm7_minus, method_pm6_dh_plus, prt_coords, prt_cart, mozyme, pdb_label
!
      use meci_C, only : maxci
!
      use elemts_C, only : elemnt
!
      use parameters_C, only : ams
!
      use common_arrays_C, only : xparam, loc, labels, nat, na, nb, nc, & 
        geo, coord, atmass, lopt, pibonds, l_atom, chains, pibonds_txt, &
        coorda, txtatm, txtatm1, break_coords, breaks, nbonds
!
      use MOZYME_C, only : start_res, lstart_res, start_letter
!
      USE funcon_C, only : fpc
      use conref_C, only : fpcref 
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  11:05:00  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use gettxt_I 
      use mopend_I 
      use getgeg_I 
      use getgeo_I 
      use geout_I 
      use wrtkey_I 
      use getsym_I 
      use symtry_I 
      use nuchar_I 
      use wrttxt_I 
      use gmetry_I 
      use reada_I
      use maksym_I 
      use to_screen_I
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ireact 
      integer , dimension(19,2) :: idepco 
      integer :: naigin, i, j, k, iflag, nreact, ij, iend, l, ii, jj, &
        i4, j4, ir_temp, l_iw, from_data_set = 14
      real(double), dimension(40) :: value 
      real(double), dimension(400) :: xyzt 
      real(double) :: degree, convrt, dum1, dum2, sum, Rab
      double precision, external :: snapth, distance
      logical :: intern = .true., aigeo, xyz, opend, exists, l_rewind = .true., l_int
      logical, allocatable :: l_use(:)
      character :: space, ch, idate*24, line_1*1000, line_2*1000, txt*15
      save  space, intern, ireact
!-----------------------------------------------
!
! MODULE TO READ IN GEOMETRY FILE, OUTPUT IT TO THE USER,
! AND CHECK THE DATA TO SEE IF IT IS REASONABLE.
! EXIT IF NECESSARY.
!
!
!
!  ON EXIT NATOMS    = NUMBER OF ATOMS PLUS DUMMY ATOMS (IF ANY).
!          KEYWRD    = KEYWORDS TO CONTROL CALCULATION
!          KOMENT    = COMMENT CARD
!          TITLE     = TITLE CARD
!          LABELS    = ARRAY OF ATOMIC LABELS INCLUDING DUMMY ATOMS.
!          GEO       = ARRAY OF INTERNAL COORDINATES.
!          LOPT      = FLAGS FOR OPTIMIZATION OF MOLECULE
!          NA        = ARRAY OF LABELS OF ATOMS, BOND LENGTHS.
!          NB        = ARRAY OF LABELS OF ATOMS, BOND ANGLES.
!          NC        = ARRAY OF LABELS OF ATOMS, DIHEDRAL ANGLES.
!          LATOM     = LABEL OF ATOM OF REACTION COORDINATE.
!          LPARAM    = RC: 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL
!          REACT     = REACTION COORDINATE PARAMETERS
!          LOC(1,I)  = LABEL OF ATOM TO BE OPTIMIZED.
!          LOC(2,I)  = 1 FOR LENGTH, 2 FOR ANGLE, AND 3 FOR DIHEDRAL.
!          NVAR      = NUMBER OF PARAMETERS TO BE OPTIMIZED.
!          XPARAM    = STARTING VALUE OF PARAMETERS TO BE OPTIMIZED.
!
!***********************************************************************
! *** IR THE TRIAL GEOMETRY  (IE.  KGEOM=0)
!   LABEL(I) = THE ATOMIC NUMBER OF ATOM(I).
!            = 99, THEN THE I-TH ATOM IS A DUMMY ATOM USED ONLY TO
!              SIMPLIFY THE DEFINITION OF THE MOLECULAR GEOMETRY.
!   GEO(1,I) = THE INTERNUCLEAR SEPARATION (IN ANGSTROMS) BETWEEN ATOMS
!              NA(I) AND (I).
!   GEO(2,I) = THE ANGLE NB(I):NA(I):(I) IR IN DEGREES; STORED IN
!              RADIANS.
!   GEO(3,I) = THE ANGLE BETWEEN THE VECTORS NC(I):NB(I) AND NA(I):(I)
!              IR IN DEGREES - STORED IN RADIANS.
!  LOPT(J,I) = -1 IF GEO(J,I) IS THE REACTION COORDINATE.
!            = +1 IF GEO(J,I) IS A PARAMETER TO BE OPTIMIZED
!            =  0 OTHERWISE.
! *** NOTE:    MUCH OF THIS DATA IS NOT INCLUDED FOR THE FIRST 3 ATOMS.
!     ATOM1  IR LABELS(1) ONLY.
!     ATOM2  IR LABELS(2) AND GEO(1,2) SEPARATION BETWEEN ATOMS 1+2
!     ATOM3  IR LABELS(3), GEO(1,3)    SEPARATION BETWEEN ATOMS 2+3
!              AND GEO(2,3)              ANGLE ATOM1 : ATOM2 : ATOM3
!
!***********************************************************************
!
      data space/ ' '/  
      data naigin/ 0/  
      data idepco/ 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 2, 0, 1, 1, 2, &
        3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 0/  
      aigeo = .FALSE. 
      nvar = 0 
      ndep = 0 
      latom = 0
      lparam = 0
      lpara1 = 0
      latom1 = 0
      lpara2 = 0 
      latom2 = 0
      breaks(1) = -300
      if (index(keywrd, " ADD-H PDBOUT") /= 0 .and. index(koment, " From PDB file") /= 0) then
        i = index(keywrd, " ADD-H")
        keywrd(i:i + 5) = " 0SCF "
        call update_txtatm(.true., .true.)
        if (allocated(txtatm1)) deallocate(txtatm1)
        allocate(txtatm1(numat))
        txtatm1(:numat) = txtatm(:numat)
        keywrd = "  LOG "//trim(keywrd)
        open(unit=ilog, form='FORMATTED', status='UNKNOWN', file=log_fn, position='asis') 
        return
      end if
      if (.not. allocated(lopt)) allocate(lopt(3,maxatoms))
   10 continue 
      line = trim(keywrd)
      keywrd = " "
      refkey_ref(1) = koment
      refkey_ref(2) = title
      call gettxt 
      if (moperr) return
      i = index(keywrd, "GEO-DAT")
      if (i /= 0) then
        keywrd(i:i+6) = "GEO_DAT"
        do j = 1, 3
          line = refkey(j)
          call upcase(line, len_trim(line))
          i = index(line, "GEO-DAT")
          if (i > 0) refkey(j)(i:i+6) = "GEO_DAT"
        end do
      end if
      i = index(keywrd, "GEO-REF")
      if (i /= 0) then
        keywrd(i:i+6) = "GEO_REF"
        do j = 1, 3
          line = refkey(j)
          call upcase(line, len_trim(line))
          i = index(line, "GEO-REF")
          if (i > 0) refkey(j)(i:i+6) = "GEO_REF"
        end do
      end if
      if (index(keywrd, "GEO_DAT")  > 0) then
        if (moperr) then
          title = " "
          koment = " "
          moperr = .false.
        end if
!
!     Use geometry in file defined by GEO_DAT
!
        do l = 1, 6
          line = " "//trim(refkey(l))
          call upcase(line, len_trim(line))
          i = index(line," GEO_DAT")
          if (i /= 0) exit
        end do
        j = index(refkey(l)(i + 10:),' ') + i + 8
        if (index(line(i:j), '"') + index(line(i:j), "'") == 0) then
          write(line,'(a)')" File name after GEO_DAT must be in quotation marks."
          call mopend(trim(line))
          return
        end if
        j = index(refkey(l)(i + 10:),'"')  + index(refkey(l)(i + 10:), "'")
        if (j == 0) then
          write(line,'(a)')" File name after GEO_DAT must end with a quotation mark."
          call mopend(trim(line))
          return
        end if
        j = j + i + 8
        line = refkey(l)(i + 9:j)
        line_1 = trim(line)
        line_2 = job_fn
        call upcase(line_1, len_trim(line_1))
        call upcase(line_2, len_trim(line_2))
        if (index(keywrd, "GEO-OK") == 0) then
          do i = len_trim(line_1), 1, -1
            if (line_1(i:i) == "\" .or. line_1(i:i) == "/") exit
          end do
          if (i > 0) line_1 = line_1(i + 1:)
          do i = len_trim(line_2), 1, -1
            if (line_2(i:i) == "\" .or. line_2(i:i) == "/") exit
          end do
          if (i > 0) line_2 = line_2(i + 1:)          
          if (line_1(:len_trim(line_1) - 3) == line_2(:len_trim(line_2) - 3)) then
            if (line_1(len_trim(line_1) - 2:) == "ARC" .or. line_1(len_trim(line_1) - 2:) == "PDB") then
              if (index(keywrd, " HTML") + index(keywrd, " PDBOUT") /= 0) then
                call mopend("The name of the geometry file defined by GEO_DAT is "// &
                & "similar to the name of the job data set")
                write(iw,"(/10x,a)")"Job name:    """//trim(job_fn)//""""
                write(iw,"(10X,a)") "GEO_DAT name:"""//trim(line)//""""
                write(iw,"(/10x,a)")"This might cause the GEO_DAT geometry data set to be over-written by the job."
                write(iw,"(10x,a)")"Either change one of the names or add ""GEO-OK"" to the keyword line, then re-run."
                return
              end if
            end if
          end if
        end if
        call add_path(line)
        inquire (file=trim(line), exist = exists)
        if (.not. exists) then
          inquire (file=trim(line)//".mop", exist = exists)
          if (exists) line = trim(line)//".mop"
          if (.not. exists) then
            call mopend ("GEO_DAT file: '"//trim(line)//"' does not exist.")
            return
          end if
        end if
        open(unit = from_data_set, file = trim(line), iostat = i)
        call upcase(line, len_trim(line))
        rewind (ir)
        if (index(line_1, ".PDB") == 0) then
          i = index(line, ".ARC")
          if (i /= 0) then
            read (from_data_set, '(A)', iostat=j) line 
            if (j == 0 .and. line(1:1) /= "*") then            
              do i = 1, 10000
                read (from_data_set, '(A)', iostat=j) line 
                if (j /= 0) exit
                if (index(line, "HEAT OF FORMATION") > 0) arc_hof_1 = reada(line,20)
                if (index(line, "FINAL GEOMETRY OBTAINED") > 0) exit   
              end do
            end if
            if (j /= 0 .or. i == 10001) rewind(from_data_set)
          end if        
          i = 0
          line_2 = trim(keywrd)
          do
            read (from_data_set, '(A241)',  iostat=j) keywrd
            if (keywrd(1:1) /= "*") then
              if (index(keywrd,"ATOM") + index(keywrd,"HETATM") /= 0) then
                write (ir, '(A)', iostat=i) trim(keywrd)
              end if
!
!  Delete all text between '"' and '"', to eliminate the unwanted text " +"
!
              call upcase(keywrd, len_trim(keywrd))
              call l_control("GEO_DAT", len_trim("GEO_DAT"), -1) 
              call l_control("GEO_REF", len_trim("GEO_REF"), -1)  
              line = trim(keywrd)
              if (i == 0 .and. index(line, " +") /= 0) i = i - 1
              line = line(:6)
              if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + &
                index(line,"HEADER") + index(line,"ANISOU") + index(line,"COMPND") + &
                index(line,"SOURCE") + index(line,"KEYWDS") + index(line,"USER ")  + &
                index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") /= 0)  exit
              i = i + 1
            end if
            if (i == 3) exit
          end do
          keywrd = trim(line_2)
        end if
        natoms = 0
        do 
          read (from_data_set, '(A241)',  iostat=i) line 
          if (i /= 0) exit
          do i = 1, 100
            if (line(i:i) /= " ") exit
          end do
          if (line(1:1) /= "*") then
            write (ir, '(A)') trim(line(i:))
            natoms = natoms + 1
          end if
        end do
        close (from_data_set)
        rewind (ir)
        maxatoms = natoms + 200
        call setup_mopac_arrays(0, 0)
        call setup_mopac_arrays(natoms + 200, 1) 
        allocate(lopt(3,maxatoms))
      end if
      chains = " "
      if (index(keywrd,"OLDGEO") .ne. 0) then 
        if (natoms == -30 .and. numat == - 30) then
          write(iw,'(/10x,a)')"OLDGEO cannot be used here."
          call mopend("A FORCE CALCULATION ON A SYSTEM THAT HAS DUMMY ATOMS CANNOT BE CONTINUED")
          return
        end if
        if (koment == " " .and. refkey_ref(1) /= " ") koment = trim(refkey_ref(1))
        if (title == " " .and. refkey_ref(2) /= " ") title = trim(refkey_ref(2))
!
!  Because OLDGEO is used, keywords that depend on the geometry are
!  not set, so copy them from "line"
!
        if (na(1) < -1) then
          natoms = 0
          return  
        end if
        ncomments = itemp_1
        if (index(keywrd, "START_RES") == 0 .and. index(line, "START_RES") .ne. 0) then
          i = index(line, "START_RES")
          j = i + 9
          do
            if(line(j:j) == ")") exit
            j = j + 1
          end do
          keywrd = line(i:j + 1)//trim(keywrd)
        end if
        if (index(keywrd, "CHAINS") == 0 .and. index(line, "CHAINS") .ne. 0) then
          i = index(line, "CHAINS")
          j = i + 9
          do
            if(line(j:j) == ")") exit
            j = j + 1
          end do
          keywrd = line(i:j + 1)//trim(keywrd)
        end if
        call update_txtatm(.true., .true.)
        if (allocated(txtatm1)) deallocate(txtatm1)
        allocate(txtatm1(numat))
        numat_old = numat
        txtatm1(:numat) = txtatm(:numat)
        if (index(keywrd, " INT") /= 0) then
          call xyzint (coord, numat, na, nb, nc, 1.d0, geo) 
          lopt(:,1) = 0
          lopt(2:,2) = 0
          lopt(3,3) = 0
        end if
        if (index(keywrd, " XYZ") /= 0) then
          numat = 0 
          do i = 1, natoms 
            if (labels(i) /= 99) then 
              numat = numat + 1 
              labels(numat) = labels(i)
              txtatm(numat) = txtatm(i)
              lopt(:,numat) = lopt(:,i)
              na(numat) = 0
            endif 
            geo(:,i) = coord(:,i) 
          end do 
!
!   If everything is marked for optimization then unconditionally mark the first
!   three atoms for optimization
!
          if (k >= 3*numat - 6) lopt(:,:min(3, numat)) = 1
          natoms = numat
        end if
        
      end if
      if (moperr) then
        natoms = 0
        return  
      end if
      if (index(keywrd,'ECHO') /= 0) then 
        rewind ir 
        if (.not.isok) then 
          write (iw, '(A)', iostat=l_iw) ' ECHO is not allowed at this point'  
          call mopend ('ECHO is not allowed at this point') 
          return  
        endif 
        isok = .FALSE. 
        do i = 1, 1000 
          read (ir, '(A)', end=60) keywrd 
          do j = 80, 2, -1 
            if (keywrd(j:j) /= ' ') go to 30 
          end do 
          j = 1 
   30     continue 
          do k = 1, j 
            if (ichar(keywrd(k:k)) >= 32) cycle  
            keywrd(k:k) = '*' 
          end do 
          write (iw, '(1X,A)', iostat=l_iw) keywrd(1:j) 
          if (l_iw /= 0) exit
        end do 
   60   continue 
        rewind ir 
        call gettxt 
        if (moperr) return  
      endif 
      if (keywrd(1:1) /= space) keywrd = " "//trim(keywrd) 
      if (koment(1:1) /= space) koment = " "//trim(koment)
      if (title(1:1) /= space)  title  = " "//trim(title)
      i = index(keywrd, ' MNDO/D')
      if (i > 0) keywrd(i:i + 6) = " MNDOD"
!
! Decide which set of fundamental constants to use
!
      if (index(keywrd,' OLDFPC') + index(keywrd, ' MNDOD') > 0) then
!
! Use old fundamental physical constants
!
        fpc(:) = fpcref(2,:)
      else
!
! Use CODATA fundamental physical constants
!
        fpc(:) = fpcref(1,:)
      endif  
      latom  = 0 
      lparam = 0 
      xyz    = index(keywrd,' XYZ') + index(keywrd,' IRC') + index(keywrd,' DRC') /= 0 
!
!   Top level
!
      if (index(keywrd,' OLDGEO') == 0) then 
!
!  Read in a new geometry
!
        nvar = 0 
        ndep = 0 
        if (allocated(nbonds)) nbonds(1) = -200
        if (allocated(txtatm1)) deallocate(txtatm1)
        allocate(txtatm1(maxatoms))
        break_coords(1,:400) = -123.d0
        if (aigeo .or. index(keywrd,' AIGIN') /= 0) then 
          intern = .false.
          call getgeg (ir, labels, geo, lopt, na, nb, nc) 
          if (moperr) return  
          if (xyz) then 
            write (iw, '(A)') &
              ' CARTESIAN CALCULATION NOT ALLOWED WITH GAUSSIAN INPUT'
              call mopend (&
               'CARTESIAN CALCULATION NOT ALLOWED WITH GAUSSIAN INPUT')  
            return  
          endif 
          if (nvar == 0) then 
            lopt(:,:natoms) = 0 
          endif 
        else 
          line = trim(keywrd)
          j = 0
          do i = 1, len_trim(line)
            if (line(i:i) == '"') j = 1 - j
            if (j == 1) line(i:i) = " "
          end do
          i = index(keywrd, " CHAIN=")
          if (i /= 0) keywrd = keywrd(:i + 5)//"S"//trim(keywrd(i + 6:))
          if (index(line, " PDB ") /= 0) then
            call getpdb(geo)
            coorda(:,:numat) = geo(:,:numat)
            numat_old = numat
          else
            call getgeo (ir, labels, geo, coord, lopt, na, nb, nc, intern) 
            if (numcal == 1 .and. natoms == 0) then
              i = index(keywrd, "GEO_DAT")
              if (i /= 0) then
                write(line,'(2a)')" GEO_DAT file """//trim(line_1)//""" exists, but does not contain any atoms."
                write(0,'(//10x,a,//)')trim(line)
                call mopend(trim(line))
              else if (.not. gui .and. numcal < 2) then
                write(line,'(2a)')" Data set '"//trim(job_fn)//" exists, but does not contain any atoms."
                write(0,'(//10x,a,//)')trim(line)
                call mopend(trim(line))
              end if
              write(line,'(a)')"(Check the first few lines of the data-set to see if there is an extra blank line."
              write(0,'(10x,a)')trim(line)
              write(iw,'(10x,a)')trim(line)
              write(line,'(a)')"If there is an extra line, delete it and re-run the job.)"
              write(0,'(10x,a)')trim(line)
              write(iw,'(10x,a)')trim(line)
              inquire(unit=ir, opened=opend)               
              if (opend) close(ir, status = 'delete', err=99) 
99            stop
            else if (natoms == -2) then
              if (l_rewind) then
                l_rewind = .false.
              else
                line ="Problem detected in data set."
                write(0,'(/10x,a,/)')trim(line)
                call mopend(trim(line))
                return
              end if
              rewind(ir)
              call getpdb(geo)
              coorda(:,:numat) = geo(:,:numat)
              numat_old = numat
            else if (natoms /= -3) then
              if (moperr .and. numcal == 1) return
              if (maxtxt > txtmax) txtmax = maxtxt
              txtatm1(:numat) = txtatm(:numat)
              if (index(keywrd, " RESID") /= 0) txtatm1(:numat)(22:22) = " "
              i = size(coorda)/3
              if (i < numat) then
                deallocate(coorda)
                allocate(coorda(3,numat))
              end if
              coorda(:,:numat) = coord(:,:numat)
              numat_old = numat
            end if
          end if
          if (index(keywrd, " HTML") + index(keywrd, " PDBOUT") /= 0 .and. maxtxt == 0) then
              call l_control("RESIDUES", len_trim("RESIDUES"), 1) 
              if (index(keywrd," MOZ") + index(keywrd," LOCATE-TS") + index(keywrd," RAPID") &
              + index(keywrd," ADD-H") == 0) &
                call l_control("CONTROL_no_MOZYME", len_trim("CONTROL_no_MOZYME"), 1)    
          end if
          if (index(keywrd, " RESIDUES ") /= 0) then
            txtatm(:numat) = " "
            maxtxt = 0
            do i = 1, numat
              if (labels(i) == 99) exit
            end do
            do i = i + 1, numat
              if (labels(i) < 99) then
                call mopend("DUMMY ATOMS CANNOT BE PRESENT WHEN KEYWORD ""RESIDUE"" IS PRESENT")
                write(iw,'(10x,a)')"(A simple way to remove dummy atoms is to add keyword ""XYZ"")"
                return
              end if
            end do            
          end if
          if (natoms == -3) then
            if (index(keywrd, " COMPARE") /= 0) then
              call mopend("SEVERE ERROR DETECTED IN ""GEO_DAT"" FILE")
              return
            end if
            goto 10
          end if
          if (moperr) return  
          if (Index (keywrd, " SNAP") /= 0) then
            !
            !   If any angles are near to important angles (such as 109.47...)
            !   snap the angle to the exact angle
            !
            do i = 1, natoms
              if (na(i) /= 0) then
                geo(2, i) = snapth (geo(2, i))
                geo(3, i) = snapth (geo(3, i))
              end if
            end do
          end if
          if (natoms < 0 ) then 
            if (numcal == 1) rewind ir 
            if (.not.isok) then 
              write (iw, '(A)') &
                ' Use AIGIN to allow more geometries to be used' 
                call mopend ('Use AIGIN to allow more geometries to be used') 
!
!   This is a deadly error - to prevent an infinite loop, kill the job.
!
              stop  
            endif 
            isok = .FALSE. 
            if (numcal > 2) then 
              naigin = naigin + 1 
              write (iw, '(2/,2A)') '   GAUSSIAN INPUT REQUIRES', &
                ' STAND-ALONE JOB' 
              write (iw, '(/,A)') '   OR KEYWORD "AIGIN"'
              call mopend (&
                 'GAUSSIAN INPUT REQUIRES STAND-ALONE JOB OR KEYWORD "AIGIN"')               
              return  
            endif 
            aigeo = .TRUE. 
            go to 10 
          endif 
        endif 
        if (natoms == 0 .and. numcal == 1) then 
          call mopend ('NO ATOMS IN SYSTEM')  
          return  
        endif 
      else 
!
!   Use the old geometry, if one exists
!
        if (numcal == 1) then
          write(line,'(a)')" Keyword OLDGEO cannot be used in the first calculation - there is no old geometry"
          write(iw,'(//10x,a)')trim(line)
          call to_screen(trim(line))
          call mopend(trim(line))
          return
        end if
      endif 
      if (natoms == 0) return
      if (index(keywrd,' FORCE')/=0 .and. labels(natoms)==107) then 
        do i = 1, na(natoms) 
          if (labels(i) /= 99) cycle  
          write (iw, '(A)') ' NO DUMMY ATOMS ALLOWED BEFORE TRANSLATION' 
          write (iw, '(A)') ' ATOM IN A FORCE CALCULATION'
          call mopend (&
       'NO DUMMY ATOMS ALLOWED BEFORE TRANSLATION ATOM IN A FORCE CALCULATION')  
          return  
        end do 
      endif 
!
!
! OUTPUT FILE TO UNIT 6
!
!    WRITE HEADER
      idate = ' ' 
      call fdate (idate) 
      write (iw, '(1X,15(''*****''),''****'')',iostat=i)
      if ( .not. gui) then
        if (i /= 0) then
          write(line,'(2a)')" Unable to write to file '", trim(output_fn)//"'"
          write(0,'(//10x,a,//)')trim(line)
          call mopend(trim(line))
          return
        else
          if (numcal == 1 .and. numat > 50) write(0,'(10x,a)')idate//"  Job: '"//trim(jobnam)//"' started successfully"
        end if
      end if
#if BITS32
      num_bits = 32 
      maxci = 5000
#else
      num_bits = 64
      maxci = 20000
#endif
      if (Academic) then
        if (site_no > 9999) then
          write (iw, '(A,i6,a,i3,a)') ' ** Site#:',site_no, &
          & '        For non-commercial use only    Version '//verson, num_bits, 'BITS **' 
        else
          write (iw, '(A,i5,a,i3,a)') ' ** Site#:',site_no, &
          & '         For non-commercial use only    Version '//verson, num_bits, 'BITS **' 
        end if
      else
        if (site_no == -1) then
          write (iw, '(A)') ' **  Evaluation copy      E-mail support: MrMOPAC@ATT.net     Version '//verson//' **'
        else if (site_no > 9999) then
          write (iw, '(A,i6,a,i3,a)') ' ** Site#:',site_no, &
          & ' E-mail support:    MrMOPAC@ATT.net    Version '//verson, num_bits, 'BITS **' 
        else
          write (iw, '(A,i5,a,i3,a)') ' ** Site#:',site_no, &
          & '  E-mail support:    MrMOPAC@ATT.net    Version '//verson, num_bits, 'BITS **' 
        end if
      end if
               
      write (iw, '(1X,a)')"*******************************************************************************"            
      write (iw, '(1X,a)')"** Cite this program as: MOPAC2016, Version: "//verson//", James J. P. Stewart,   **"
      

      if (ijulian > 410) then
        write (iw, '(1X,a, a,a)') &
          "** Stewart Computational Chemistry, web-site: HTTP://OpenMOPAC.net.          **"
      else
        write (iw, '(1X,a, i4,a)') &
                          "** Stewart Computational Chemistry, web: HTTP://OpenMOPAC.net. Days left:",ijulian,"**"
      end if
      write (iw, '(1X,a)')"*******************************************************************************"
      write (iw, '(1X,a)')"**                                                                           **"
      write (iw,"(1x,a)") "**                                MOPAC2016                                  **"
      write (iw, '(1X,a)')"**                                                                           **"
      write (iw, '(1X,a)')"*******************************************************************************"
      j = len_trim(keywrd)
      do i = j, 1, -1
        ch = keywrd(i:i)
        if (ichar(ch) > 126) keywrd(i:i) = " "
        if (ichar(ch) < 32) exit
      end do
      if (i /= 0) then
        write(line,'(a)')" Non-standard characters detected in keyword line"
        call mopend(trim(line))
        l = 0
        write(iw,*)
        do
          l = l + 20
          k = min(l,i)
          if (k < l - 19) exit
          write(iw,'(a,i3,30i4)')"Position:", (j, j = l - 19, k)
          write(iw,'(a,i3,30i4)')"   ASCII:", (ichar(keywrd(j:j)), j = l - 19, k)
          write(iw,'(a,a3,29a4)')"    Char:", (keywrd(j:j), j = l - 19, k)
          write(iw,*)
          if (k == i) exit
        end do
        write(iw,'(//,a,//)')" (Edit keyword line to remove non-standard characters using a primitive editor such as vi)"
      end if   
      if (.not. is_PARAM) then
        keywrd_txt = trim(keywrd)
!
!  Delete all text that might contain keywords
!
        call l_control("GEO_DAT", len_trim("GEO_DAT"), -1) 
        call l_control("GEO_REF", len_trim("GEO_REF"), -1)  
!
!  Decide which method to use
!
        do i = 1, n_methods
          methods(i) = (index(keywrd, trim(methods_keys(i))//" ") /= 0)
        end do 
        do i = 1, n_methods
          if (methods(i)) exit
        end do
        if (i > n_methods) then
!
!  Is it a known variant of PM6 or PM7?
!
          j = index(keywrd, " PM6") + index(keywrd, " PM7")
          if (j /= 0) then
            k = index(keywrd(j + 1:), " ") + j - 1
            call mopend('The method requested: "'//keywrd(j + 1:k)//'" does not exist.')
            return
          end if
        end if
!
! Default method is PM7
!
        if (.not. method_pm7) then
           do i = 1, n_methods
             if (methods(i)) exit
           end do
           if (i > n_methods) then
             i = 14
             method_pm7 = .true.
           end if
        end if 
        caltyp = methods_keys(i)    
!
!  Define parent methods "method_PM6" and "method_PM7" for variants of these methods
!
        method_pm7 = (method_PM7 .or. method_PM7_ts .or. method_pm7_hh .or. method_pm7_minus)
        method_pm6 = (method_PM6 .or. method_pm6_dh2 .or. method_pm6_d3h4 .or. method_pm6_dh_plus .or. &
       & method_pm6_dh2x .or. method_pm6_d3h4x .or. method_pm6_d3 .or. method_pm6_d3_not_h4)
        dh = " "
        i = index(keywrd, " PM6-D")
        if (i /= 0) then 
          j = index(keywrd(i + 6:)," ") + i + 4
          dh = keywrd(i + 5:j)                                  
        else if (index(keywrd, " PM6-H") /= 0)  then  
          dh = "H   "
        end if
        keywrd = trim(keywrd_txt)
      endif      
      write (iw, &
      '(/24X,A,'' CALCULATION RESULTS'',2/1X,15(''*****''),''****'' )') "     "//trim(caltyp) 
      write (iw,'(" *  CALCULATION DONE: ",31x,2a)') idate,"  *"
!
! Copy all keywords to keywrd_txt
!
      keywrd_txt = trim(keywrd)
      line_1 = trim(keywrd)
!
!  Delete all text that might contain keywords
!
      call l_control("GEO_DAT", len_trim("GEO_DAT"), -1) 
      call l_control("GEO_REF", len_trim("GEO_REF"), -1)  
      ch = '"'
      l = len_trim(line_1)
      j = 0
      do i = 1, l
        if (keywrd(i:i) == ch) then
          j = -j + 1
        end if
        if (j == 1) keywrd(i:i) = " "
      end do     
      if (j == 1) then
        call mopend("NUMBER OF QUOTATION MARKS, '""', IN KEYWORDS IS ODD. THIS NUMBER MUST BE EVEN.")
        return
      end if
      keywrd = trim(line_1)
      if (index(keywrd, " HTML(NORES)") /= 0) call l_control("HTML NORJSMOL", len_trim("HTML NORJSMOL"), 1)  
!
! WRITE KEYWORDS BACK TO USER AS FEEDBACK
!
      call wrtkey
      if (moperr) &
        write(iw,'(a)')" *", &
      & " *  Errors detected in keywords.  Job stopped here to avoid wasting time."," *"          
      write (iw, &
        '(1X,14(''*****''),''*********'')')
      if (moperr) return
      if (index(keywrd, " LOCATE-TS") /= 0) lopt(:,:natoms) = 1
      if (index(keywrd, "INVERT") /= 0) then
        nvar = 0
        do i = 1, natoms
          do j = 1,3
            if (lopt(j,i) == 1) then
              lopt(j,i) = 0
            else if (lopt(j,i) == 0) then
              lopt(j,i) = 1
            end if
          end do
        end do
      end if
      start_res = -200
      lstart_res = .false.
      lstart_res(1) = .true.
      i = index(keywrd," START_RES")
      if (i /= 0) then       !                  Parse the "START_RES" keyword for residues, chains, and breaks
        i = i + 10
        if (keywrd(i:i) == "=") i = i + 1
        k = 1000
        do j = 1, 100
          start_res(j) = nint(reada(keywrd, i)) - 1
            k = i 
          do
            k = k + 1
            if (keywrd(k:k) < "0" .or. keywrd(k:k) > "9") exit
          end do
          ch = " "
          if (keywrd(k:k) >= "A" .and. keywrd(k:k) <= "Z") ch = keywrd(k:k)
          start_letter(j) = ch
          do 
            if (keywrd(k:k) == "-" .and. ((keywrd(k-1:k-1) >= "0" .and. keywrd(k-1:k-1) <= "9") &
                                    .or. (keywrd(k-1:k-1) >= "A" .and. keywrd(k-1:k-1) <= "Z"))) exit
            if (keywrd(k:k) == " ") exit
            if (keywrd(k:k) == ")") exit
            k = k + 1
          end do
          if (keywrd(k:k) == ")") exit
          if (keywrd(k:k) == "-" .and. ((keywrd(k-1:k-1) >= "0" .and. keywrd(k-1:k-1) <= "9") &
                                  .or. (keywrd(k-1:k-1) >= "A" .and. keywrd(k-1:k-1) <= "Z")))  then
            i = k + 1
          else
            if (keywrd(k:k) == ")") exit
            l = index(keywrd(i:k), " ") !  Start of a new protein
            if (l /= 0) then
              i = i + l 
              lstart_res(j + 1) = .true.
            end if
          end if
        end do 
!
!  Identify atoms where chain breaks occur
!
        j = 1
        ij = 0
        jj = nint(reada(txtatm(1), 22)) - 1
        nbreaks = 0
        do i = 1, natoms
          if (labels(i) == 6) exit
        end do
        if (index(keywrd, " NOTER") /= 0 .or. i > natoms) then
          start_res(2) = -200
          break_coords = 0.d0
          breaks = -100
        end if
        do ii = 2, 10
          if (start_res(ii) == -200) exit
          inner_loop: do k = j, numat
            l = nint(reada(txtatm(k), 22)) - 1
            if (l == start_res(ii) .and. txtatm(k)(22:22) == start_letter(ii)) then
              if (lstart_res(ii)) then
                nbreaks = nbreaks + 1    
                do i4 = k - 1, 1, -1
                  if (labels(i4) /= 1) exit
                end do
                if (i4 > 0) then
                  break_coords(:,nbreaks) = geo(:,i4)
                else
                  nbreaks = nbreaks - 1
                end if
              end if
              ij = ij + 1
              breaks(ij) = k - 1
              j = k 
              exit inner_loop
            end if
          end do inner_loop
        end do
        j = 1
        do i = 2, nbreaks
          if (breaks(i) - breaks(i - 1) > 1) then
            j = j + 1
            breaks(j) = breaks(i)
          end if
        end do
        nbreaks = j
        breaks(nbreaks + 1) = 0
      end if
      i = index(keywrd," CHAINS")  
      if (i /= 0) then
        i = i + index(keywrd(i:),"(")
        do j = i, i + 26
          k = ichar(keywrd(j:j))
          if (k < ichar("A") .or. k > ichar("Z")) then
            do k = j - i + 1, 26
              chains(k) = chains(j - i)
            end do
            exit
          end if
          chains(j - i + 1) = keywrd(j:j)
        end do        
      else
    !
    ! Allow for up to 26 chains
    !
        do i = 1, 26
          chains(i) = char(ichar("A") + i - 1)
        end do
      end if
!
! FILL IN GEO MATRIX IF NEEDED
!
      if (index(keywrd,' SYMM') + index(keywrd," SYM ") /= 0 .or. ndep > 0) then 
        call getsym(locpar, idepfn, locdep, depmul) 
        if ((xyz .or. intern)) then          
          ndep = 0 
          if (index(keywrd," IRC") + index(keywrd," DRC") == 0) then           
            if (index(keywrd,' XYZ') /= 0) write(iw,"(13x,a)")"(Remove either SYMMETRY or XYZ)"
            if (intern) write(iw,"(13x,a)")"(Remove either SYMMETRY or INT)"
            if (index(keywrd, " 0SCF") == 0) then
              call mopend('SYMMETRY CANNOT BE USED WHEN COORDINATE SYSTEMS ARE CHANGED')
              return
            else
              moperr = .false.
            end if
          end if         
        endif 
      endif 
      
      if (index(keywrd, " NOOPT") /= 0)   then
        k = 0
        do_noopt: do i = 1, natoms
          do j = 1,3
            if (lopt(j,i) < 0) exit do_noopt
          end do
        end do do_noopt
        if (j < 4 .or. i <= natoms) then
          write(iw,*)" ""NOOPT"" cannot be used when an optimization flag of -1 is present"
          call mopend("""NOOPT"" cannot be used when an optimization flag of -1 is present")
          return
        end if
        line = trim(keywrd)
        do
          i = index(line, " NOOPT")
          if (i == 0) exit
          j = index(line(i + 4:)," ") + i + 2
          line(i:i + 6) = " "     
          if (j - i == 8) then
            line(1:1) = line(j - 1:j - 1)
            line(2:2) = char( ichar(line(j:j)) + ichar("a") - ichar("A"))
          else if (j - i == 7) then
            line(1:1) = " "
            line(2:2) = line(j:j)
          else
            do i = 1, natoms
              lopt(1, i) = 0
              lopt(2, i) = 0
              lopt(3, i) = 0
            end do
            cycle
          end if               
          do l = 1, 99
            if (line(1:2) == elemnt(l)) exit
          end do        
          do i = 1, natoms
            j = labels(i)
            if (j == l) then
              lopt(1, i) = 0
              lopt(2, i) = 0
              lopt(3, i) = 0
            end if
          end do
        end do 
      end if
      if (index(keywrd, " OPT") /= 0)   then
        i = index(keywrd, " OPT=")
        if (i /= 0) keywrd(i + 4:) = keywrd(i + 5:)
        k = 0
        do_opt: do i = 1, natoms
          do j = 1,3
            if (lopt(j,i) < 0) exit do_opt
          end do
        end do do_opt
        if (j < 4 .or. i <= natoms) then
          call mopend("""OPT"" cannot be used when an optimization flag of -1 is present")
          return
        end if
        line = trim(keywrd)
        i = index(line, " OPT")
!
! "j" marks the end of the keyword
!
        j = 0
        if (index(line, " OPT(") /= 0) j = index(line(i + 1:), ") ")
        if (j /= 0) then
!
! Test to ensure that there are no other "OPT" keywords are present 
!
          if (index(keywrd, " OPT ") + index(keywrd, " OPT-") /= 0) then
            exists = moperr
            call mopend('When keyword "OPT(text)" is present, other keywords containing "OPT" must not be present.')
            if (index(keywrd, " 0SCF") /= 0) then
              moperr = exists
            else
              return
            end if
          end if
          allocate(l_use(natoms))
          l_use = .false.
          lopt(:,:natoms) = 0
          i = i + 5
          k = index(line(i:j + i), "=")
          if (k /= 0) then
!
! Keyword OPT("text"=n.nn) used, so use only atoms in the region of selected atoms.
!
!
!  Text for selected atoms
!
            txt = line(i + 1: k + i - 3)
            Rab = reada(line, k + i)
            jj = 0
            do i = 1, len_trim(txt)
              if (txt(i:i) /= " ") then
                jj = jj + 1
                txt(jj:jj) = txt(i:i)
              end if
            end do
            txt(jj + 1:) = " "    
            if (txt(1:1) >= "0" .and. txt(1:1) <= "9") then
              txt = "*"//trim(txt)
              jj = jj + 1
            end if
            lopt(:,:natoms) = 0
            do i = 1, natoms
!
!  l_use(i) = F:  not assigned
!  l_use(i) = T:  Assigned because it matches the SITE word or it's within Rab of an atom defined by SITE
!
              line = txtatm(i)(22:)
              j = 1
              do ii = 1, len_trim(line)
                if (line(ii:ii) /= " ") then
                  line(j:j) = line(ii:ii)
                  j = j + 1
                end if
              end do
              if (len_trim(txt) == j - 1) then
                do k = 1, jj
                  if (txt(k:k) /= line(k:k) .and. txt(k:k) /= "*") exit
                end do
                if (k > jj) then
                  l_use(i) = .true.
                  do j = 1, natoms
                    if (l_use(j)) cycle
                    l_use(j) = (distance(i,j) < Rab) 
                  end do
                end if
              end if
            end do    
            if (index(keywrd, " NOOPT") /= 0) then
              do i = 1, natoms
                if ( .not. l_use(i)) cycle
                line = txtatm(i)(22:)
                j = 1
                do ii = 1, len_trim(line)
                  if (line(ii:ii) /= " ") then
                    line(j:j) = line(ii:ii)
                    j = j + 1
                  end if
                end do
                do k = 1, jj
                  if (txt(k:k) /= line(k:k) .and. txt(k:k) /= "*") exit
                end do
                l_use(i) = (k < jj + 1)
              end do
            end if
          else  
!
!   Keyword OPT("text1"[,"text2"[,"text3"]]...) used, so use only atoms in the defined residues.
!
            i4 = i 
            j4 = index(line(i4 + 1:), '"') + i4
!
!  Text for selected atoms
!
            do
              txt = line(i4 + 1:j4 - 1)
              jj = 0
              do i = 1, len_trim(txt)
                if (txt(i:i) /= " ") then
                  jj = jj + 1
                  txt(jj:jj) = txt(i:i)
                end if
              end do
              txt(jj + 1:) = " "    
              if (txt(1:1) >= "0" .and. txt(1:1) <= "9") then
                txt = "*"//trim(txt)
                jj = jj + 1
              end if
              do i = 1, natoms
!
!  l_use(i) = F:  not assigned
!  l_use(i) = T:  Assigned because it matches the SITE word
!
                line_1 = txtatm(i)(22:)
                j = 1
                do ii = 1, len_trim(line_1)
                  if (line_1(ii:ii) /= " ") then
                    line_1(j:j) = line_1(ii:ii)
                    j = j + 1
                  end if
                end do
                if (len_trim(txt) == j - 1) then
                  do k = 1, jj
                    if (txt(k:k) /= line_1(k:k) .and. txt(k:k) /= "*") exit
                  end do
                  if (k > jj) l_use(i) = .true.
                end if
              end do   
!
! Move on to the next residue
!
              i4 = j4 + 2 
              if (line(i4:i4) /= '"') exit
              j4 = index(line(i4 + 1:), '"') + i4
            end do
          end if
          do i = 1, natoms
            if (l_use(i)) lopt(:,i) = 1
          end do
        else
          do
            i = index(line, " OPT")
            if (i == 0) exit
            j = index(line(i + 2:)," ") + i 
            line(i:i + 4) = " "     
            if (j - i == 6) then
              line(1:1) = line(j - 1:j - 1)
              line(2:2) = char( ichar(line(j:j)) + ichar("a") - ichar("A"))
            else  if (j - i == 5) then
              line(1:1) = " "
              line(2:2) = line(j:j)
            else if (j - i == 3) then
              if (natoms > 1) then
                l_int = na(2) /= 0
              else
                l_int = .false.
              end if              
              do i = 1, natoms
                if (i == 1 .and. l_int) cycle                
                lopt(1, i) = 1
                if (i == 2 .and. l_int) cycle 
                lopt(2, i) = 1
                if (i == 3 .and. l_int) cycle 
                lopt(3, i) = 1
              end do
              cycle
            else
              cycle
            end if               
            do l = 1, 99
              if (line(1:2) == elemnt(l)) exit
            end do        
            do i = 1, natoms
              j = labels(i)
              if (j == l) then
                lopt(1, i) = 1
                lopt(2, i) = 1
                lopt(3, i) = 1
              end if
            end do
          end do 
        end if
      end if
!
!   FORCE OPTIMIZATION FLAGS OFF FOR DEPENDENT COORDINATES
!
      j = 1 
      if (xyz) j = 2 
      do i = 1, ndep 
        lopt(idepco(idepfn(i),j),locdep(i)) = 0 
      end do 
      if (ndep /= 0) call symtry 
!
! INITIALIZE FLAGS FOR OPTIMIZE AND PATH
!
      iflag = 0 
      latom = 0 
      numat = 0 
      if (nvar /= 0) then 
        numat = natoms 
      else 
        if (index(keywrd, " FORCETS") /= 0) then
          do i = 1, natoms 
            if (lopt(1,i) /= lopt(2,i) .or. lopt(1,i) /= lopt(3,i)) then
              write(iw,'(/10x,a,i5)')" Fault on line:",i
              call mopend("All optimization flags for an atom must be the same.")
              return
            end if
          end do
        end if
        loc = 0
        if (index(keywrd, " FORCETS") /= 0 .and. index(keywrd," RESTART") /= 0) then
          k = 1
        else
          k = 2
        end if
        do i = 1, natoms
          if (labels(i) /= 99 .and. labels(i) /= 107) numat = numat + 1 
          if (k == 1) then
            if (lopt(1,i) /= lopt(2,i) .or. lopt(1,i) /= lopt(3,i)) then
              write(iw,'(/10x,a,i5)')" Fault on line:",i
              call mopend("All optimization flags for an atom must be the same")
              return
            end if
          end if
          do j = 1, 3 
            if (lopt(j,i) == 1 .or. lopt(j,i) == k) go to 200 
            if (lopt(j,i) > -1) cycle   
            if (iflag /= 0) then 
              if (index(keywrd,' STEP1') /= 0) then 
                lpara1 = lparam 
                latom1 = latom 
                lpara2 = j 
                latom2 = i 
                latom = 0 
                iflag = 0 
                cycle  
              else 
               call mopend ('ONLY ONE REACTION COORDINATE PERMITTED') 
                return  
              endif 
            endif 
            latom = i 
            lparam = j 
            convrt = 1.d0 
            if (j > 1 .and. na(latom) > 0) convrt = 0.01745329252D00
            if (allocated(react)) deallocate(react)
            allocate(react(4000))
            react(1) = geo(lparam, latom)
            ireact = 1 
            iflag = 1 
            cycle  
!    FLAG FOR OPTIMIZE
  200       continue 
            nvar = nvar + 1 
            loc(1,nvar) = i 
            loc(2,nvar) = j 
            xparam(nvar) = geo(j,i) 
          end do 
        end do 
      endif 
      if (index(keywrd, " MINI ") /= 0) then 
        nl_atoms = 0
        if (index(keywrd, " FORCETS") /= 0) then
          k = 0
        else
          k = 1
        end if
        do i = 1, natoms 
          l_atom(i) = (abs(lopt(1,i)) > k) 
          if (l_atom(i)) nl_atoms = nl_atoms + 1
        end do
        if (nl_atoms == 0) then
          if (index(keywrd, " 0SCF") == 0) then
            line = " Keyword 'MINI' used, but no atoms flagged for printing (optimization flag '2')"
            write(iw,'(a)')trim(line)
            call mopend(trim(line))
            return
          end if
        end if
      else
        nl_atoms = numat
        l_atom = .true.
      end if
      line_1 = trim(keywrd)
      call l_control("GEO_DAT", len_trim("GEO_DAT"), -1) 
      call l_control("GEO_REF", len_trim("GEO_REF"), -1)  
      i = index(keywrd," RESTART")
      keywrd = line_1
      if (i /= 0 .and. index(keywrd,' IRC') + index(keywrd,'FORCE') + index(keywrd," THERMO") == 0) then
        inquire (file=restart_fn, exist = exists)
        if (.not. exists) goto 1900
        open (unit=ires, file=restart_fn, status="UNKNOWN", &
                   & form="UNFORMATTED")
        rewind (ires)
               !
               !  Read in the geometric variables
               !
         read (ires, end=1900, err=1900) i,i, (xparam(i), i=1, nvar)
         do i = 1, nvar
           k = loc(1, i)
           l = loc(2, i)
           geo(l, k) = xparam(i)
         end do
         if (ndep /= 0) call symtry 
         call gmetry (geo, coorda)
         if (index(keywrd, " XYZ") /= 0) na(:numat) = 0
        end if
! READ IN PATH VALUES
      if (iflag /= 0) then 
        if (Index (line, " TS") /= 0) then
          call mopend("""TS"" cannot be used when an optimization flag of -1 is present")
          return
        end if 
        if (Index (line, " SIGMA") /= 0) then
          call mopend("""SIGMA"" cannot be used when an optimization flag of -1 is present")
          return
        end if 
        if (Index (line, " NLLSQ") /= 0) then
          call mopend("""NLLSQ"" cannot be used when an optimization flag of -1 is present")
          return
        end if 
        if (index(keywrd,' SIGMA') /= 0) then 
          write (iw, '(A)') &
            ' SIGMA USED WITH REACTION PATH; THIS OPTION IS NOT ALLOWED'
          call mopend (&
             'SIGMA USED WITH REACTION PATH; THIS OPTION IS NOT ALLOWED')  
          return  
        endif 
        if (index(keywrd,' STEP=') + index(keywrd,' POINT=') + index(keywrd,' 0SCF') /= 0) then 
          go to 250 
        endif 
  220   continue 
        read (ir, '(A)', end=240) line 
        call nuchar (line, len_trim(line), value, nreact) 
        if (nreact == 0) go to 240 
        do i = 1, nreact 
          ij = ireact + i 
          react(ij) = value(i)*convrt 
          if (abs(react(ij)-react(ij-1)) >= 1.D-12) cycle  
          dum1 = react(ij)/convrt 
          dum2 = react(ij-1)/convrt 
          write (iw, &
      '(3/,'' TWO ADJACENT POINTS ARE IDENTICAL:  '',    F7.3,2X,F7.3,/,'' THIS&
      & IS NOT ALLOWED IN A PATH CALCULATION'')') dum1, dum2 
          call mopend (&
       'TWO ADJACENT POINTS ARE IDENTICAL: THIS IS NOT ALLOWED IN A PATH CALCULATION') 
          return  
        end do 
        ireact = ireact + nreact 
        go to 220 
  240   continue 
        degree = 1.D0 
        if (lparam > 1 .and. na(latom) > 0) degree = 57.29577951308232D0 
         if (index(keywrd,' 0SCF') + index(keywrd,' RESEQ') == 0) then 
          if (ireact <= 1) then 
            call mopend ('NO POINTS SUPPLIED FOR REACTION PATH')  
            write (iw, '(2/10X,'' GEOMETRY AS READ IN IS AS FOLLOWS'')') 
            xparam(1) = -1.D0 
            call geout (1)
            return  
          else 
            write (iw, '(2/10X,'' POINTS ON REACTION COORDINATE'')') 
            write (iw, '(10X,8F8.2)') (react(i)*degree,i=1,ireact) 
          endif 
          iend = ireact + 1 
          react(iend) = -1.D12 
        end if
      end if 
  250 continue 
      if (nvar > 0 .and. index(keywrd,' PM7-TS') /= 0) nvar = 0
      call wrttxt (iw) 
      if (index(keywrd,' 0SCF') == 0) then 
!
! CHECK DATA
!
        if (xyz) then 
          if (index(keywrd,' IRC') + index(keywrd,' DRC') + &
          index(keywrd,' 1SCF') + index(keywrd,' FORCE') + index(keywrd,' DFORCE') == 0) then 
            if (nvar /= 0 .and. intern .and. nvar < 3*numat-6) then 
              write (iw, &
      '(2/10X,''INTERNAL COORDINATES READ IN, AND CALCULATION '',/10X, &
      & ''TO BE RUN IN CARTESIAN COORDINATES, '',/10X, &
      & ''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')') 
              write (iw, &
      '(2/10X,''THIS INVOLVES A LOGICALLY ABSURD CHOICE'',/10X, &
      & '' SO THE CALCULATION IS TERMINATED AT THIS POINT'')') 
              call mopend ('INCONSISTENT USE OF OPTIMIZATION FLAGS') 
              return  
            endif 
          endif 
        else 
          if (index(keywrd,' FORCE') + index(keywrd,' DFORCE') + index(keywrd,' 1SCF') == 0 .or. &
            index(keywrd,' GRAD') /= 0) then 
            if (intern .and. nvar/=0 .and. nvar<3*numat-6) then 
              write (iw, &
         '(2/10X,'' CARTESIAN COORDINATES READ IN, AND CALCULATION '',/10X, &
      & ''TO BE RUN IN INTERNAL COORDINATES, '',/10X, &
      & ''BUT NOT ALL COORDINATES MARKED FOR OPTIMISATION'')') 
              write (iw, &
      & '(2/10X,''MOPAC, BY DEFAULT, USES INTERNAL COORDINATES'',/10&
      & X,''TO SPECIFY CARTESIAN COORDINATES USE KEY-WORD :XYZ:'')') 
              write (iw, &
      & '(10X,''YOUR CURRENT CHOICE OF KEY-WORDS INVOLVES A LOGICALLY'',/10X, &
      & ''ABSURD CHOICE SO THE CALCULATION IS TERMINATED AT THIS POINT'')') 
              call mopend ('INCONSISTENT USE OF OPTIMIZATION FLAGS') 
              return  
            endif 
          endif 
        endif 
      endif 
      if (index(keywrd,' LOG') + index(keywrd," ADD-H") /= 0 .or. &
      (index(keywrd," 0SCF") /= 0 .and. index(keywrd," OLDGEO") /= 0 .and. &
       index(keywrd," PDBOUT") /= 0)) then 
        inquire(unit=ilog, opened=opend) 
        if (.not. opend) open(unit=ilog, form='FORMATTED', status='UNKNOWN', file=log_fn, position='asis') 
        call wrttxt (ilog) 
       endif 
      if (index(keywrd," OLDGEO") /= 0) call delete_ref_key("OLDGEO", len_trim("OLDGEO"), ' ', 1)
      if (index(keywrd, " NOTXT") /= 0) then
        maxtxt = 0
        txtatm = " "
      end if
      if (prt_coords) call geout(1)
      write(iw,*)
      if (moperr) return 
!
!  Check for isotopes.  If found, print them out
!
      j = 0
      k = 0
      do i = 1, natoms
        if (labels(i) == 99 .or. labels(i) == 107) cycle
        k = k + 1
        if (Abs(atmass(k) - ams(labels(i))) > 1.d-3) then
          if (j == 0) then
            write(iw,"(/30x, a,/)")"  Isotopes Used"
            j = 1
          end if
        write(iw,"(a,i5,1x,2a,f8.4,a,f8.4,a)")"              Atom:", i, elemnt(labels(i)), &
   "  Default mass", ams(labels(i)), ",  mass used:", atmass(k), " amu"
        end if 
      end do    
      call gmetry (geo, coord)
      if(index(keywrd," PKA") /= 0) then
!
!  Do any ionizable hydrogen atoms exist (hydrogens attached to an oxygen)
!
        k = 0
        i = 0
        do ii = 1, natoms
          if (labels(ii) < 99) i = i + 1
          if (labels(ii) /= 1) cycle        
          j = 0
          do jj = 1, natoms
            if (labels(ii) < 99) j = j + 1
            if (labels(jj) /= 8) cycle           
            sum = (coord(1,i) - coord(1,j))**2 + &
                & (coord(2,i) - coord(2,j))**2 + &
                & (coord(3,i) - coord(3,j))**2 
            if (sum < 1.69d0) then
              k = 1
              exit         
            end if
          end do
        end do
        if (k == 0) then
          write(iw,"(/,3(10x,a,/))")"A request was made to print the pKa values for this system,", &
      "but there are no hydrogen atoms attached to an oxygen atom,", &
      "so the pKa calculation cannot be completed."
          call mopend("No '-O-H' groups found.  pKa cannot be calculated")
          return
        end if
      end if 
      if (moperr) return  
      pdb_label = (maxtxt > 25)
      use_ref_geo = (index(keywrd_txt," GEO_REF") /= 0) 
      if (use_ref_geo) then
!
!  "GEO-OK+" is used only by COMPARE.  It can be read as "If one or more hydrogen atoms in the two systems have different labels,
!  ignore the difference in labels and continue with the comparison."
!
        if (index(keywrd, " COMPAR") /= 0) then
          if (index(keywrd_txt, " LET") /= 0) then
            call l_control("0SCF HTML GEO-OK+ LET NOCOM", len_trim("0SCF HTML GEO-OK+ LET NOCOM"), 1) 
          else
            call l_control("0SCF HTML GEO-OK LET NOCOM", len_trim("0SCF HTML GEO-OK LET NOCOM"), 1) 
          end if
        end if
        if (index(keywrd,' RESIDUES0') /= 0 .and. index(keywrd,' 0SCF') /= 0 .and. index(keywrd, " HTML") /= 0) then
          nat(:numat) = labels(:numat)
          call geochk
        end if
        if (index(keywrd,' SYMM') + index(keywrd," SYM ") /= 0) then
          call mopend("""SYMMETRY"" cannot be used with ""GEO_REF"" (""GEO_REF"" requires all coordinates to be optimized)")
          return
        end if
        call geo_ref
        lopt = 1
      end if
      if (index(keywrd," NOCOM") /= 0) then
        ncomments = 0
        do j = 1, 3
          line = refkey(j)
          call upcase(line, len_trim(line))
          i = index(line, " NOCOM")
          if (i > 0) then
            do i = i + 1, i + 10
              if (refkey(j)(i:i) == " ") exit
              refkey(j)(i:i) = " "
            end do
          end if            
        end do
      end if
      if (index(keywrd, " SETPI") /= 0) then
        i = min(10, numat)
        if (allocated(pibonds)) deallocate(pibonds)
        if (allocated(pibonds_txt)) deallocate(pibonds_txt)
        allocate (pibonds(i,2), pibonds_txt(i))
        pibonds_txt = " "
        pibonds = 0
!
!  The user supplies some pi bonds.
!
        i = index(keywrd, " SETPI")
        j = index(keywrd(i + 6:), " ") + i + 4
        l = index(keywrd(i:j), "=")
        if (l /= 0) then
!
!  Read pi bonds from a file in this folder
!
          l = l + i
          if (keywrd(l:l) == '"') then
            do j = l + 1, len_trim(keywrd)
              if (keywrd(j:j) == '"') exit
            end do              
            l = l + 1
            j = j - 1
          end if
          line = keywrd(l:j)
          call add_path(line)
          inquire (file=line, exist = exists)
          if (exists) then
            ir_temp = 27
            open (unit=ir_temp, file=line)  
            rewind (ir_temp)          
          else
            if (index(keywrd,' 0SCF') == 0 ) then 
              write(line_1,'(a)')" Pi-bond file """//trim(line)//""" does not exist."
              call mopend(trim(line_1))
              write(line_1,'(a)')"(required by SETPI="""//trim(line)//""")"
              call mopend(trim(line_1))
              return
            end if
          end if
          k = index(keywrd," SETPI=") + 7
          if (k > 7) then
            do ii = 1, 6
              line = " "//trim(refkey(ii))
              call upcase(line, len_trim(line))
              i = index(line," SETPI=")
              if (i /= 0) exit
            end do
            if (keywrd(k:k) == '"') then
              j = index(refkey(ii)(i + 9:),'" ') + i + 7
            else
              j = index(refkey(ii)(i + 9:),' ') + i + 6
            end if
            refkey(ii) = refkey(ii)(:i + 4)//refkey(ii)(j + 2:)
          end if 
        else
          ir_temp = ir
        end if
        l = 0
        do 
          read (ir_temp,'(a)',end=98, err=98) line
          if (line(1:1) == "*") cycle
          if (line == " ") exit
          i = index(line, " User")
          if (i /= 0) line(i:) = " "
          pibonds_txt(l + 1) = trim(line)
          call upcase(line, len_trim(line))
          i = index(line, '"')
          if (i > 0) then
            do i4 = 1,2
              i = index(line, '"')
              if (i /= 0) then
                j = index(line(i:),"""") + i
                if (j == i) exit
                ij = index(line(j:),"""") + j 
                line_1 = line(j:ij - 2)
!
! Do NOT use "txt_to_atom_no" in the next block - the format must be in PDB
!
                if (line_1(1:1) == "[") then
!
! Atom defined using Jmol format
!
                  ii = index(line(j:ij - 2),".") + j
                  k = index(line(j:ij - 2),":") + j
                  jj = index(line(j:ij - 2),"]") + j
                  if (k == j) then
                    line_1 = line(ii:ij - 2)//line(j + 1:jj - 2)//line(jj:ii - 2)
                  else
                    line_1 = line(ii:ij - 2)//line(j + 1:jj - 2)//line(k:k)//line(jj:k - 2)
                  end if
                else
                  j = 0
                  do i = 1, len_trim(line_1)
                    if (line_1(i:i) /= " ") then
                      j = j + 1
                      line_1(j:j) = line_1(i:i)
                    end if
                  end do
                  line_1(j + 1:) = " "
                end if 
                do i = 1, numat
                  line_2 = txtatm(i)
                  ij = 0
                  do k = 13, txtmax
                    if (line_2(k:k) /= " ") then
                      ij = ij + 1
                      line_2(ij:ij) = line_2(k:k)
                    end if
                  end do
                  line_2(ij + 1:) = " " 
                  ij = max(ij, len_trim(line_1))
                  do j = 1, ij
                    if (line_1(j:j) /= line_2(j:j) .and. line_1(j:j) /= "*") exit
                  end do
                  if (j > ij) exit
                end do
                if (i > numat) then
                  i = index(line,'"')
                  if (i > 0) then
                    j = index(line(i + 1:), '"') + i
                    write(line_1,'(a)')"Pi-bond atom defined by "//line(i:j)//" does not exist."
                    call mopend(trim(line_1))
                    return
                  end if                    
                end if
                ii = index(line, '"')
                jj = index(line(ii + 1:),"""") + ii
                write(line(ii:jj),'(i5)')i
              end if
            end do 
          end if
!
!  Read in the first of the two atom numbers
!
          do i = 1, len_trim(line) 
            if (line(i:i) /= " ") then 
              ii = nint(reada(line,i))
              exit
            endif 
          end do 
          do i = i + 1, len_trim(line) 
            if (line(i:i) == " ") then 
              jj = nint(reada(line,i))
              exit
            endif 
          end do 
          l = l + 1
          if (ii > numat .or. jj > numat) then
            write(line,'(a,i5,a,i5,a)') " pi bond between atoms",ii," and", jj," is impossible - atom number too large"
            call mopend(trim(line))
            return
          end if
          pibonds(l,1) = ii
          pibonds(l,2) = jj
        end do
98      if (l == 0) then
          if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") == 0 .and. mozyme) then 
            call mopend("Keyword SETPI used, but no pi bonds specified")
            call web_message(iw,"SETPI.html")
            return
          end if
        else
         write(iw,'(/15x,a)')"Keyword SETPI used, pi-bonds specified are:"
         write(iw,'(/12x,a)')"   Bond No.                   Atom               to              Atom"
         do i = 1, l
           write(iw,'(12x,i7,12x,a,a,a)')i, &
           '"'//txtatm(pibonds(i,1))(:maxtxt)//'"', "   -   ", '"'//txtatm(pibonds(i,2))(:maxtxt)//'"'
         end do
        end if
        if (ir_temp == 27) close(ir_temp)
      end if  
      k = index(keywrd," GEO_DAT") + 10
      call delete_ref_key("GEO_DAT=", len_trim("GEO_DAT="), '" ', 2)
      if (moperr) return
      if (index(keywrd,' AUTOSYM') /= 0) call maksym(loc, xparam, xyzt)
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") /= 0) then  
        inquire(unit=iarc, opened=opend) 
        if (opend) close (iarc)
        if (index(keywrd, "PDBOUT") /= 0) archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"pdb" 
      endif 
      if (prt_cart .and. (maxtxt < 26 .and. (index(keywrd,' NOXYZ') == 0 .or. gui))) then 
         write (iw, '(2/10X,''CARTESIAN COORDINATES '',/)') 
        write (iw, &
      '(4X,''NO.'',7X,''ATOM'',11X,''X'',11X,''Y'',11X,''Z'',/)') 
        l = 0 
        do i = 1, natoms 
          if (labels(i)==99 .or. labels(i)==107) cycle  
          l = l + 1 
          if (l_atom(i)) write (iw, '(I6,8X,A2,4X,3F12.4)') l, elemnt(labels(i)), (coord(j,l),j=1,3) 
          end do 
      endif 
      return  
1900  if (exists) then
        call mopend("RESTART file is corrupt")
      else
        call mopend("RESTART file '"//trim(restart_fn)//"' does not exist.")
      end if        
      return
      end subroutine readmo 
      double precision function snapth (theta)
      implicit none
      double precision, intent (in) :: theta
      integer :: i, j, k
      double precision :: angle, const, phase, sum
      intrinsic Abs, Acos, Asin, Cos, Int, Mod, Nint, Sign, Sqrt, Dble
        angle = Cos (theta)
        phase = Sign (1.d0, theta)
        if (Abs (angle) < 1.d-4) then
          !
          !   Cos(Theta) is zero - this is 90 or 270 degrees
          !
          const = 2 * Asin (1.d0)
          if (Abs (theta) < const) then
            snapth = phase * Acos (0.d0)
          else
            snapth = phase * Acos (0.d0) + const
          end if
        else
          sum = Abs (1.d0/angle) ** 2
          do i = 1, 7
            j = Nint (sum*i)
            if (Abs (j-sum*i) < 5.d-3) go to 1000
          end do
          snapth = theta
          return
    1000  sum = Sqrt (Dble(i)/Dble(j))
          const = 2 * Asin (1.d0)
          k = Int (Abs (theta)/const)
          if (Mod(k, 2) == 0) then
             !
             !   Theta is in domain -180 - 0 - +180 degrees
             !
            snapth = phase * Acos (sign(sum, angle))
          else
             !
             !   Theta is in domain -180 - 360 - +180 degrees
             !
            snapth = phase * (2*const-Acos (sign(sum, angle)))
          end if
        end if
  end function snapth
  subroutine geo_diff(sum, rms, prt)
!
!  geo_diff calculates the total distance between geometries "geo" and "geoa"
!
!  Two measures are used:
!
!   sum = addition of all differences in position, i.e., sum of motions of atoms.
!   rms = sum of squares of differences in positions, 
!   i.e., the square of the vector sum of differences
!
    use common_arrays_C, only : geo, geoa, txtatm, nat
    use molkst_C, only : numat, keywrd, txtmax
    use chanel_C, only : iw
    use elemts_C, only : elemnt
    implicit none
    double precision :: sum, rms
    logical :: prt, precise
    double precision, allocatable :: move(:)
    integer :: i, j, k, lim
    logical :: first 
    double precision :: dum1, sum1
    character :: num*7
    save :: first
    allocate (move(numat))
!
!  Calculate the average and RMS differences of two geometries
!
    precise = (index(keywrd, " PREC") /= 0)
    sum = 0.d0
    rms = 0.d0
    do i = 1, numat
      dum1 =  (geo(1,i) - geoa(1,i))**2 + &
              (geo(2,i) - geoa(2,i))**2 + &
              (geo(3,i) - geoa(3,i))**2
      rms = rms + dum1
      dum1 = sqrt(dum1)
      move(i) = dum1
      sum = sum + dum1
    end do 
    first = .true.
    sum1 = 0.d0
    lim = numat
    if (.not. precise) lim = min(30, lim)
    do i = 1, lim
      dum1 = -1.d0
      do j = 1, numat
        if (dum1 < move(j)) then
          dum1 = move(j)
          k = j
        end if
      end do      
      if (dum1 > 0.1d0 .or. (dum1 > 0.003d0 .and. i < 11) .or. (precise .and. dum1 > 0.001d0)) then
        if (prt .and. first) then
          first = .false.
          write(iw,'(/21x,a,//2x,a,/)')"Atoms that move a lot", &
          "                   Atom Label             GEO_REF Coordinates       Movement    Integral"
        end if
        if (prt) then
          sum1 = sum1 + dum1
          if (dum1 > 0.1d0) then
            num = "12.2,2x"
          else
            num = "14.4   "
          end if
          if (txtatm(k) /= " ") then
            if (prt) write(iw,'(8x,  4x, a, 3f8.3, f'//num//', f11.2)')"("//txtatm(k)(:txtmax)//")", geoa(:,k), dum1, sum1
          else
            if (prt) write(iw,'(1x, i5, 18x, a2, 14x, 3f8.3, f'//num//', f11.2)')k, elemnt(nat(k)),  geoa(:,k), dum1, sum1
          end if        
          move(k) = -2.d0
        end if
      else
        exit  
      end if
    end do
    deallocate (move)
    return
  end subroutine geo_diff
  subroutine add_path(file)
    use chanel_C, only : job_fn
    use molkst_C, only : verson, good_separator
    implicit none
    character, intent (inout) :: file*(*)
    integer :: i, j, n
    character :: path*241
    logical :: need_path, windows
!
!  Test to see if the file "file" already has an absolute path.
!
    windows = (verson(7:7) == "W")
    if (windows) then      
      need_path =  (file(2:2) /= ":" .and. file(1:1) /= good_separator) 
    else
      need_path =  (file(1:1) /= good_separator .and. file(1:1) /= "~")
    end if
    if (.not. need_path) return
!
!  Test to see if the job has an absolute path.
!    
    if (windows) then
      need_path = (job_fn(2:2) == ":")
    else
      need_path = (job_fn(1:1) == good_separator)
    end if
    if (.not. need_path) return      
!
!  file "file" does not include the path and the path is present in job_fn
!  therefore add path to the file "file"
!
!  First check if a relative path is used
!
    n = 1
    do 
      if (file(1:3) /= ".."//good_separator .and. file(1:3) /= ".."//good_separator) exit
      n = n + 1
      file = trim(file(4:))
    end do
!
!  Is the file defined relative to the current folder?
!  If so, delete the definition.
!
    if (file(1:2) == "."//good_separator .or. file(1:2) == "."//good_separator) file = trim(file(3:))
!
!  Delete the data-set-name from the path, and delete the folders up to the relative folder.
!
    path = trim(job_fn)
    do j = 1, n
      do i = len_trim(path) - 1, 1, -1
        if (path(i:i) == good_separator) exit
      end do
      path = path(:i)
    end do
!
!  Finally, join the path to the start of the file
!
    file = path(:i)//trim(file)   
    return
  end subroutine add_path
  subroutine delete_ref_key(beginning, length_beginning, ending, length_ending)
!
! Delete keyword from refkey so that it will not appear in an ARC file
!
!   beginning :        Start of keyword to be deleted (MUST be uppercase)
!   length_beginning : Length of string in "beginning"
!   ending :           End of keyword (MUST be uppercase)
!   length_ending :    Length of string "ending"
!
  use molkst_C, only : refkey, line
  implicit none
  integer, intent (in) :: length_beginning, length_ending
  character, intent (in) :: beginning*(length_beginning), ending*(length_ending)
  integer :: i, j, k
  logical :: in_text
  do
    k = 0
    do i = 1, 6
      line = " "//trim(refkey(i))
      call upcase(line, len_trim(line))
      k = index(line," "//beginning)
      if (k /= 0) exit
    end do
    if (k == 0) exit
!
!  Search from the start of the keyword (k) to find the end of the keyword (j)
!
    j = index(refkey(i)(k + length_beginning:), ending) + k + length_beginning 
!
!  Search for the first non-blank character (j) after the keyword, if it exists.
!
    do j = j, len_trim(refkey(i))
      if (refkey(i)(j:j) /= " ") exit
    end do
!
!  Search for the first non-blank character (k) before the keyword, if it exists.
!
    do k = k - 1, 1, -1
      if (refkey(i)(k:k) /= " ") exit
    end do    
!
! Eliminate the keyword
!
    refkey(i) = refkey(i)(:k)//" "//refkey(i)(j:)
!
! Cycle, in case there is more than one instance of the keyword.
!
  end do
!
!  Compress keywords, if possible.
!  This involves deleting " + " and " & ", if present and the resulting line is not too long,
!  and removing excess spaces between keywords
!
  if (index(refkey(1), " + ") + index(refkey(1), " & ") > 0) then
    if (index(refkey(2), " + ") + index(refkey(2), " & ") > 0) then
      i = len_trim(refkey(3))
      j = len_trim(refkey(2))
      if (j < 120) then
        k = index(refkey(2), " + ") + index(refkey(2), " & ") + 1
        refkey(2)(k:k) = " "
        refkey(2) = refkey(2)(:j)//refkey(3)(:i)
        refkey(3) = "    NULL"
      end if 
    end if
    i = len_trim(refkey(2))
    j = len_trim(refkey(1))
    if (j < 60) then
      k = index(refkey(1), " + ") + index(refkey(1), " & ") + 1
      refkey(1)(k:k) = " "
      refkey(1) = refkey(1)(:j)//refkey(2)(:i)
      if (index(refkey(2), " + ") + index(refkey(2), " & ") > 0) then
        refkey(2) = trim(refkey(3))
        refkey(3) = "    NULL"
      else        
        refkey(2) = "    NULL"
      end if
    end if
!
! First line of keywords MUST be printed
!
    i = index(refkey(1), " NULL")
    if (i /= 0) refkey(1) = refkey(1)(:i)//refkey(1)(i + 5:)
    do i = 1, 3
      in_text = .false. 
      k = 1
      do j = 2, len_trim(refkey(i))
        if (refkey(i)(j:j) == '"') in_text = (.not. in_text)
        if (in_text .or. refkey(i)(j-1:j) /= "  ") then
          k = k + 1
          refkey(i)(k:k) = refkey(i)(j:j)
        end if
      end do
      refkey(i)(k + 1:) = " "
    end do       
  end if
  do i = 1, 3
    if (len_trim(refkey(i)) > 0) then
      do
        if (refkey(i)(1:1) /= " ") refkey(i) = " "//trim(refkey(i))
        if (refkey(i)(1:2) == "  ") refkey(i) = trim(refkey(i)(2:))
        if (refkey(i)(2:2) /= " ") exit
      end do
    end if
  end do
  return
  end subroutine delete_ref_key
  
  
  

  
  
