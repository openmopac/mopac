  subroutine run_mopac
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double
!
      use common_arrays_C, only : nfirst, nlast, nat, xparam, grad, nw, &
      p, pa, pb, labels, loc, time_start, l_atom, coord, txtatm, coorda, &
      txtatm1
!
      USE molkst_C, only : gnorm, natoms, numat, nvar, numcal, job_no, nscf, id, &
        escf, iflepo, iscf, keywrd, last, moperr, maxatoms, ncomments, verson, &
        time0, atheat, errtxt, isok, mpack, gui, line, na1, refkey, keywrd_txt, &
        press, mozyme, step_num, jobnam, nelecs, stress, E_disp, E_hb, E_hh, no_pKa, &
        MM_corrections, lxfac, trunc_1, trunc_2, method_PM8, bad_separator, good_separator, &
        method_PM7, method_PM6, method_RM1, sparkle, itemp_1, maxtxt, koment, &
        num_threads, nl_atoms, use_ref_geo, prt_coords, pdb_label, txtmax
!
      USE parameters_C, only : tore, ios, iop, iod, eisol, eheat, zs, eheat_sparkles
!
      Use Parameters_for_PM7_Sparkles_C, only : gss7sp
      Use Parameters_for_PM6_Sparkles_C, only : gss6sp
!
      use cosmo_C, only : iseps, useps, lpka, solv_energy
!
      USE funcon_C, only : fpc_9, fpc
!
      USE maps_C, only : latom, react, rxn_coord
!
      use symmetry_C, only : state_Irred_Rep
!
      USE chanel_C, only : ir, iw, iarc, output_fn, end_fn, iend, &
        archive_fn, log, ilog, xyz_fn, job_fn, log_fn
!
      use MOZYME_C, only : rapid, nres
!
      use meci_C, only : nmos
!
      use elemts_C, only : elemnt, atom_names

      use reada_I
      use to_screen_I

       Use mod_vars_cuda, only: lgpu   
#if GPU
      Use iso_c_binding 
      Use mod_vars_cuda, only: ngpus, gpu_id
      Use gpu_info
      Use settingGPUcard
#endif
      use second_I
      use geout_I
      use wrttxt_I
      use geoutg_I
      use datin_I
      use fbx_I
      use fordd_I
      use calpar_I
      use compfg_I
      use react1_I
      use grid_I
      use paths_I
      use pathk_I
      use force_I
      use drc_I
      use nllsq_I
      use powsq_I
      use ef_I
      use flepo_I
      use writmo_I
      use polar_I
      use pmep_I
      implicit none
      integer ::  i, j, l
      real(double) :: eat,  tim 
      logical :: exists, opend, sparkles
      double precision, external :: C_triple_bond_C
      character :: nokey(20)*10
      integer, external :: mkl_get_max_threads
#if GPU
      logical :: lgpu_ref
      logical(c_bool)    :: hasGpu = .false.
      logical(c_bool)    :: lstat = .false.
      logical(c_bool)    :: hasDouble(6)
      integer(c_int)     :: nDevices
      character*256	     :: gpuName(6)
      integer(c_size_t)  :: totalMem(6)
      logical            :: gpu_ok(6)
      character*3        :: on_off(6)
      integer(c_int), dimension(6)	 :: clockRate, major, minor, name_size, k
#endif
!-----------------------------------------------
      tore = ios + iop + iod
      call fbx                            ! Factorials and Pascal's triangle (pure constants)
      call fordd                          ! More constants, for use by MNDO-d
      inquire (directory = "C:/", exist = exists)
      if (exists) then
        bad_separator = "/"
        good_separator = "\"
        if (verson(7:7) == " ") verson(7:7) = "W"
      else
        bad_separator = "\"
        good_separator = "/"
        if (verson(7:7) == " ") verson(7:7) = "L"
      end if
      lgpu = .false.
      trunc_1 = 7.0d0    ! Beyond 7.0 Angstroms, use exact point-charge
      trunc_2 = 0.22d0   ! Multiplier in Gaussian: exp(-trunc_2*(trunc_1 - Rab)^2)
!
! Read in all data; put it into a scratch file, "ir"
!
      call getdat(ir,iw)
      call to_screen("To_file: Start of reading in data")
      if (natoms == 0 .or. moperr) return
!
!   CLOSE UNIT IW IN CASE IT WAS ALREADY PRE-ASSIGNED
!
      close(iw)
      l = 0
 11   open(unit=iw, file=output_fn, status='UNKNOWN', position='asis', iostat = i)
      if (i /= 0) then
        l = l + 1
        write(0,"(i3,3a)")21 - l," File """,output_fn(:len_trim(output_fn)),""" is unavailable for use"
        write(0,"(a)")"    Correct fault (probably the output file is in use elsewhere)"
        call sleep(5)
        if (l < 20) goto 11
        write(0,"(a)") "Job abandoned due to output file being unavailable"
        call mopend("Job abandoned due to output file being unavailable")
        call sleep(5)
        goto 101
      end if
      if (l > 0) then
       write(0,"(a)")" Fault successfully corrected.  Job continuing"
      end if
      rewind iw
      isok = .TRUE.
      errtxt = 'Job stopped by operator'
      tim = second(1)
      call date_and_time(VALUES=time_start)
      if (moperr) goto 101
!
! Set up essential arrays, these are the arrays that are needed for reading in the data
!
      natoms = natoms + 200
      call setup_mopac_arrays(natoms, 1)
      maxatoms = natoms
   10 continue
      numcal = numcal + 1      ! A new calculation
      job_no = job_no + 1      ! A new job
      step_num = step_num + 1  ! New electronic structure, therefore increment step_num
      moperr = .FALSE.
      escf   = 0.d0
      gnorm  = 0.D0
      press  = 0.d0
      E_disp = 0.d0
      E_hb   = 0.d0
      E_hh   = 0.d0
      solv_energy = 0.d0
      nres = 0
      nscf = 0
      nmos = 0
      na1 = 0
      lpka = .false.
      stress = 0.d0
      no_pKa = 0
      time0 = second(1)
      MM_corrections = .false.
      pdb_label = .false.
      state_Irred_Rep = " "
      if (numcal > 1) call to_screen("To_file: Leaving MOPAC")
      if (numcal > 1 .and. numcal < 4 .and. index(keywrd_txt," GEO_DAT") /= 0) then
!
!  Quickly jump over first three lines
!
        close (ir)
        i = natoms
        call getdat(ir,iw)
        natoms = i
        call gettxt
      end if
!
!    Read in all the data for the current job
!
      i = numcal
      call readmo
      if (moperr .and. numcal == 1 .and. natoms > 1) goto 101
      if (moperr .and. numcal == 1 .and. index(keywrd_txt," GEO_DAT") == 0) goto 100
      if (moperr) goto 101
      if (numcal == 1) then
        num_threads = min(mkl_get_max_threads(), 20)
        i = index(keywrd, " THREADS")
        if (i > 0) then
          i = nint(reada(keywrd, i))
          num_threads = min(max(1,i), num_threads)
        end if
        call mkl_set_num_threads(num_threads)
#if GPU
        gpuName(1:6) = '' ; name_size(1:6) = 0 ; totalMem(1:6) = 0 ; clockRate(1:6) = 0
        hasDouble(1:6) = .false. ; gpu_ok(1:6) = .false.
        clockRate(1:6) = 0 ; major(1:6) = 0 ; minor(1:6) = 0; on_off(1:6) = 'OFF'
        call gpuInfo(hasGpu, hasDouble, nDevices, gpuName,name_size, totalMem, &
                  & clockRate, major, minor)
        lgpu_ref = hasGPU
        if (lgpu_ref) lgpu_ref = (index(keywrd, " NOGPU") == 0)
        if (lgpu_ref) then         
          lgpu_ref = .false.
! Counting how many GPUs are suitable to perform the calculations or with compute capability 2 (Fermi or Kepler).
          j = 0
          do i = 1, nDevices
            if (major(i) >= 2  .and. hasDouble(i)) then
              gpu_ok(i) = .true.
              j = j + 1
            endif
          end do          
          lgpu_ref = (j >= 1)
        end if
        ngpus = 1  ! in this version only single-GPU calculation are performed. This variable control it.
!       ngpus = j ! in future versions of MOPAC Multi-GPUs calculations should be allowed
        
        l = 0; k = 0                       
        if (lgpu_ref) then
          l = index(keywrd,' SETGPU=')
          if (l /= 0) then  ! The user has inserted SETGPU keyword to Select one specific GPU 
            gpu_id = nint(reada(keywrd,l))    
            if (gpu_id > nDevices .or. gpu_id < 1 .or. (.not. gpu_ok(gpu_id))) then              
              ! the user made a wrong choice !!!
              Write(iw,'(/,5x,a)') ' Problem with the definition of SETGPU keyword ! '
              Write(iw,'(5x,a,/)') ' MOPAC will automatically set a valid GPU card for the calculation '
              l = 0
            else
              on_off(gpu_id) = 'ON '
              call setGPU(gpu_id - 1, lstat)
              if (.not. lstat) then
                write (6,*) 'Problem to set GPU card ID = ', gpu_id
                stop
              endif                             
            endif                        
          endif
          if (l == 0) then   ! Select GPU automatically               
            do i = 1, nDevices
             if (gpu_ok(i)) then
                on_off(i) = 'ON '
                gpu_id = i - 1
                call setGPU(gpu_id, lstat)
                if (.not. lstat) then
                  write (6,*) 'Problem to set GPU card ID = ', gpu_id
                  stop
                endif                 
                exit
              endif
            enddo                  
          endif
        else
          nDevices = 0
          ngpus = 0
        endif
!
!  For small systems, using a GPU takes longer than not using a GPU,
!  so do not use a GPU for small systems.  The lower limit, 100, is just a guess.
!
        lgpu = (lgpu_ref .and. natoms > 100) ! Warning - there are problems with UHF calculations on small systems
#endif
      endif
      if (.not. gui .and. numcal == 1 .and. natoms == 0) then
        write(line,'(2a)')" Data set exists, but does not contain any atoms."
        write(0,'(//10x,a,//)')trim(line)
        call mopend(trim(line))
        write(0,'(5x,a)')" (Check the first few lines of the data-set for an extra blank line."
        write(0,'(5x,a)')"  If there is an extra line, delete it and re-submit.)"
        inquire(unit=ir, opened=opend)    
        if (opend) close(ir, status = 'delete', err = 999)
        stop
      end if
      if (numat > 46000) then
        write(line,'(a,i5,a)')"Data set '"//trim(jobnam)//"' exists, but at ",numat," atoms is too large to run."
        write(0,'(//10x,a,//)')trim(line)
        call mopend(trim(line))
        write(iw,'(10x,a)')"(Absolute maximum size: 46,000 atoms.)"
        numat = 0
        natoms = 0
        goto 100
      end if
      if (numcal == 1 .and. moperr .or. natoms == 0) then
!
!   Check for spurious "extra" data
!
        do i = 1,15
          read(ir,*, iostat=l)line
          if (l /= 0) exit
        end do
        if (i > 14) then
          write(iw,"(//10x,a)") " WARNING: There are extra data at the end of the input data set."
          write(iw,"(10x,a,/)")"          There might be an error in the data set."
        end if
        goto 101
      end if
      if (moperr) then
        if (keywrd(:7) == " MODEL") goto 101
        go to 10
      end if
      lxfac = (index(keywrd," XFAC") /= 0)
!
! Load in parameters for the method to be used
!
      call switch
    !  if (method_PM8) method_PM7 = .true.
      if (index(keywrd,' EXTERNAL') /= 0) call datin (iw)
      if (moperr) go to 100
      sparkle = (index(keywrd, " SPARKL") /= 0)
      sparkles = .false.
      if (.not. method_RM1) then
        do i = 1, natoms
          if (labels(i) > 57 .and. labels(i) < 72) then
            sparkles = .true.
            exit
          end if
        end do
      end if
      if ((method_PM7 .or. method_PM6 .or. method_RM1) .and. sparkles) then
        if (.not. sparkle) then
          line = " "
          do j = 1, 10
            if (atom_names(labels(i))(j:j) /= " ") exit
          end do
          write (line, '(A,a)') ' Data are not available for ', atom_names(labels(i))(j:)//"."
          if (index(keywrd, " 0SCF") == 0) then
            if (method_PM7 .and. Abs(gss7sp(labels(i))) > 0.1d0 .or. &
                method_PM6 .and. Abs(gss6sp(labels(i))) > 0.1d0) &
              write(iw,*)" (Parameters are available if SPARKLE is used)"
            call mopend(trim(line))
            goto 100
          end if
        end if
      else
        sparkle = .false.
      end if
      do i = 57,71
        if (zs(i) < 0.1d0) tore(i) = 3.d0
      end do
!
! Set up all the data for the molecule
!
      call moldat (0)  ! data dependent on the system
      if (index(keywrd, " 0SCF") /= 0) moperr = .FALSE.
      call calpar      ! Calculate derived parameters
      if (moperr) goto 100
      call to_screen("To_file: Data read in")
!
!  If no SCF calculations are needed, output geometry and quit
!
      if (.false.) then
        inquire(unit=iarc, opened=opend)
        if (opend) close (iarc)
        open(unit=iarc, file=archive_fn(:len_trim(archive_fn) - 3)//"mop", status='UNKNOWN', position='asis')
        rewind iarc
!
!   Remove unwanted keywords
!
        call upcase(refkey(1), len_trim(refkey(1)))
        nokey = " "
        nokey(1) = "NOOPT"
        nokey(2) = "OPT"
        nokey(3) = "LET"
        nokey(4) = "GEO_R"
        nokey(5) = "GNOR"
        nokey(6) = "PL"
        nokey(7) = "MOZY"
        nokey(8) = "PDBO"
        nokey(9) = "SETU"
        nokey(10) = "PM7"
        do j = 1, 20
          if (nokey(j) == " ") exit
          do
            i = index(refkey(1), " "//trim(nokey(j)))
            if (i == 0) exit
            l = index(refkey(1)(i + 1:)," ") + i + 1
            refkey(1) = refkey(1)(:i)//refkey(1)(l:)
          end do
        end do
        refkey(1) = trim(refkey(1))//' SETUP'
        call geout (iarc)
        stop
      end if
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") /= 0 ) then
        inquire(unit=iarc, opened=opend)
        if (opend) close (iarc)
        if (index(keywrd, " PDBOUT") /= 0) then
          line = archive_fn(:len_trim(archive_fn) - 3)//"pdb"
        else
          line = trim(archive_fn)
        end if
        if (index(keywrd, " HTML") /= 0) then
          do i = len_trim(line), 1, -1
            if (line(i:i) == "/" .or. line(i:i) == "\") exit
          end do
        end if
        open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
        rewind iarc
        if (index(keywrd,' 0SCF') /= 0 .and. index(keywrd, " MINI") /= 0 .and. nl_atoms > 0) then
          open(unit=l, file=xyz_fn)
          write(l,"(i6,a)") nl_atoms," "
          write(l,*)"POINT "
          do i = 1, numat
            if (l_atom(i)) write(l,"(3x,a2,3f15.5)")elemnt(nat(i)), (coord(j,i),j=1,3)
          end do
        end if
      end if
      if (maxtxt == 0 .and. index(keywrd, " RESIDUES") /= 0) call geochk()
      if ( index(keywrd," PDBOUT") /= 0 .and. maxtxt < 26 .and. index(keywrd," RESID") == 0) then
        if (maxtxt == 0) then
          maxtxt = 26
          do i = 1, numat
            write(txtatm(i),'(a,i5,1x,a,a)')"HETATM",i, elemnt(nat(i)),"   HET A   1"
          end do
          txtatm1(:numat) = txtatm(:numat)
        else
          write(line,'(a)')"PDBOUT only works when the atom labels are in PDB format or keyword RESIDUES is also present"
          call mopend(trim(line))
          write(iw,'(/10x,a)') &
          "(Before using PDBOUT, either add keyword RESIDUES or run a job using keyword RESIDUES to add PDB atom labels.)"
          write(iw,'(10x,a)') &
          "(Keyword RESIDUES can only be used when one of MOZYME, LEWIS, CHARGES, or RESEQ is also present)"
          return
        end if
      end if   
      if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) nelecs = 0
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) then
        if (index(keywrd, " DISP") /= 0) then
          call l_control("0SCF", len_trim("0SCF"), 1)   
          call l_control("PRT", len_trim("PRT"), 1)   
          call post_scf_corrections(eat, .false.)
        end if
        if (index(keywrd,' XYZ') /= 0) then
          line = ' GEOMETRY IN CARTESIAN COORDINATES'
        else if (index(keywrd,' INT') /= 0) then
          line = ' GEOMETRY IN MOPAC Z-MATRIX FORMAT'
        else
          line = ' GEOMETRY OF SYSTEM SUPPLIED'
        end if
        if (prt_coords) write (iw, '(A)') trim(line)
        xparam(1) = -1.D0
        if (index(keywrd," OLDEN") /= 0 .and. index(keywrd, " 0SCF") == 0) then
!
! read in density so that charges can be calculated
!
          if (mozyme) then
            call set_up_MOZYME_arrays()
          else
            if (allocated(p))      deallocate (p)
            if (allocated(pa))     deallocate (pa)
            if (allocated(pb))     deallocate (pb)
            allocate(p(mpack), pa(mpack), pb(mpack))
          end if
          call den_in_out(0)
          if (moperr) return
          if (mozyme) call density_for_MOZYME (p, 0, nelecs/2, pa)
        end if
        if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") > 0) then
!
!  Add RESEQ by default
!
          if (index(keywrd, " NORES") == 0 .and. (maxtxt == txtmax .or. index(keywrd, " RESID") /= 0) &
          .and. index(keywrd, " RESEQ") == 0) call l_control("RESEQ", len("RESEQ"), 1)
        end if
        if (prt_coords) call geout (iw)
        if (index(keywrd,' AIGOUT') /= 0) then
          write (iw, '(2/,A)') '  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT'
          call wrttxt (iw)
          call geoutg (iw)
          write (iarc, '(2/,A)') '  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT'
          call wrttxt (iarc)
          call geoutg (iarc)
        else if (mozyme .or. &
          (index(keywrd," PDBOUT") + index(keywrd," RESEQ") + index(keywrd," RESID") /= 0)) then
          i = size(coorda)
          j = size(coord)
          if (i /= j) then
            deallocate (coorda)
            allocate(coorda(3,natoms))
            coorda = coord
          end if
          if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") == 0) call geochk()
          if (index(keywrd, " ADD-H") == 0 .and. index(keywrd, " SITE=") == 0 .and. index(keywrd," RESEQ") /= 0) then
            moperr = .false.
            goto 10
          end if
          if (index(keywrd, " SITE=") + index(keywrd, " ADD-H") /= 0) then
            moperr = .false.
            if (index(keywrd, " ADD-H") /= 0) then
              call add_hydrogen_atoms()
              if (moperr) then
                inquire(unit=iarc, opened=opend) 
                if (opend) close (iarc, status="DELETE") 
                go to 101
              end if  
              if (index(keywrd, " NORESEQ") /= 0) call update_txtatm(.true., .false.) 
              call l_control("ADD-H", len("ADD-H"), -1)
              numat = natoms - id
              numcal = numcal + 1
              call mopend("HYDROGEN ATOMS ADDED")
              moperr = .false.
            end if
            call l_control("0SCF", len("0SCF"), 1)
            call geochk()
            if (moperr) goto 101
            if (size(coorda) >= size(coord)) then
              coorda(:,:numat) = coord(:,:numat)
              txtatm1(:numat) = txtatm(:numat) 
            end if 
            moperr = .false.
            i = index(refkey(1), "ADD-H")
            if (i /= 0) refkey(1) = refkey(1)(:i - 1)//refkey(1)(i + 5:)
          end if
          if (index(keywrd, " SITE=") + index(keywrd, " ADD-H") + index(keywrd," RESEQ") + &
            index(keywrd," RESID") /= 0) &
            call update_txtatm(.true., .false.)         !  Now that geometry checks are done, switch to input labels
          call write_sequence
          if (log) call bridge_H()
          if (moperr) goto 101
          inquire(unit=iarc, opened=opend)
          if (opend) close(iarc)
          archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"arc"
          open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
          rewind iarc
          if (index(keywrd, " NOOPT") + index(keywrd," OPT") == 0)   then
            if (index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) then
!
!  Atoms added, deleted, or re-arranged, so set all optimization parameters to "opt"
!
              l = 0
              do i = 1, numat 
                do j = 1,3
                  l = l + 1
                  loc(1,l) = i
                  loc(2,l) = j
                end do
              end do
            end if
          end if
          call geout (iarc)
          if ( index(keywrd," PDBOUT") /= 0) then           
            inquire(unit=iarc, opened=opend)
            if (opend) close(iarc)
              line = archive_fn(:len_trim(archive_fn) - 3)//"pdb"
              if (index(keywrd, " HTML") /= 0) then
                do i = len_trim(line), 1, -1
                  if (line(i:i) == "/" .or. line(i:i) == "\") exit
                end do
              end if
            open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
            rewind iarc
            call pdbout(iarc)
          end if
        else
          if (index(keywrd, " ADD-H") /= 0) then
            call add_hydrogen_atoms()
          end if
          inquire(unit=iarc, opened=opend)
          if (opend) close (iarc)
          open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
          rewind iarc
          if (pdb_label) then
!
!   PDB format, so renumber atoms in atom label
!
            j = 1
            do i = 1, numat
              write(txtatm(i),'(a6,i5,a15)')txtatm(i)(:6),i + j - 1,txtatm(i)(12:)     
            end do 
          end if
          call geout (iarc)
        end if
        go to 100
      endif
      if (pdb_label) call compare_txtatm(moperr, moperr)
      if (moperr) then
        if (index(keywrd," GEO-OK") /= 0) then
          moperr = .false.
        else     
          goto 100
        end if
      end if
!
!  If any special work is done, do it here
!
  !    call special
   !   go to 10
!
! Everything is ready - now set up the arrays used by the SCF, etc.
!
      useps = .FALSE.
      iseps = (index(keywrd,' EPS=') + index(keywrd," PKA") /= 0)
      call setup_mopac_arrays(1,2)
      iseps = (index(keywrd,' EPS=') /= 0)
      if (moperr) goto 100
      if (allocated(nw)) deallocate(nw)
      allocate(nw(numat))
      l = 1
      do i = 1, numat
        nw(i) = l
        l = l + ((nlast(i)-nfirst(i)+1)*(nlast(i)-nfirst(i)+2))/2
      end do
!
!  CALCULATE THE ATOMIC ENERGY
!
      if (sparkle) then
        atheat = 0.d0
        do i = 1, numat
          if  (nat(i) > 56 .and. nat(i) < 72 .and. zs(nat(i)) < 0.1d0) then
            atheat = atheat + eheat_sparkles(nat(i))
          else
            atheat = atheat + eheat(nat(i))
          end if
        end do
      else
        atheat = sum(eheat(nat(:numat)))
      end if
      eat = sum(eisol(nat(:numat)))
      atheat = atheat - eat*fpc_9
      atheat = atheat + C_triple_bond_C()
      rxn_coord = 1.d9
!
!  All data for the current job are now read in, and all parameters are
!  available in the arrays.
!  Now decide what type of calculation is to be done.
!
      if (mozyme) then
        call geochk()
        if (moperr) then
          if ( index(keywrd," PDBOUT") /= 0) then
            archive_fn = archive_fn(:len_trim(archive_fn) - 3)//"pdb"
            inquire(unit=iarc, opened=opend)
            if (opend) close(iarc)
            open(unit=iarc, file=archive_fn, status='UNKNOWN', position='asis')
            rewind iarc
            call pdbout(iarc)
          end if
!
! If the current calculation cannot continue, increment numcal to indicate that a partial calculation has already been done.
!
          numcal = numcal + 1
          goto 10 
        end if
        if (index(keywrd, " RESID") + index(keywrd, " RESEQ") + index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) &
          call update_txtatm(.true., .false.)         !  Now that geometry checks are done, switch to input labels
        if (moperr) then
!
! If the current calculation cannot continue, increment numcal to indicate that a partial calculation has already been done.
!
          numcal = numcal + 1
          goto 10
        end if
        rapid = (index(keywrd, " RAPID") + index(keywrd, " LOCATE-TS") /= 0)
        call set_up_MOZYME_arrays()
        if (moperr) goto 101
        if (index(keywrd, " RAPID") /= 0) call set_up_rapid("ON")
      end if
      if (index(keywrd,' 1SCF') /= 0) then
        iflepo = 1
        iscf = 1
        last = 1
        i = index(keywrd,' GRAD')
        grad(:nvar) = 0.D0
        numcal = numcal + 1
        tim = second(1)
        call to_screen(" Single point calculation")
        call compfg (xparam, .TRUE., escf, .TRUE., grad, i /= 0)
      else if (index(keywrd,' SADDLE') /= 0) then
        call to_screen(" Transition state geometry calculated using SADDLE")
        call react1 ()
      else if (index(keywrd,' STEP1') /= 0) then
        call to_screen(" Grid calculation")
        call grid ()
        iflepo = -1  !  Prevent printing of results
      else if (latom /= 0) then
        call to_screen(" Path calculation")
        if (index(keywrd,' STEP')==0 .or. index(keywrd,' POINT')==0) then
          call paths ()
        else
          call pathk ()
        end if
        iflepo = -1
      else if (index(keywrd,' FORCE') + index(keywrd,' IRC=') + &
        index(keywrd,' THERM') + index(keywrd,' DFORCE') /= 0) then
        call to_screen(" Force constant calculation")
        last = 1
        call force ()
        iflepo = -1
      else if (index(keywrd,' DRC') + index(keywrd,' IRC') /= 0) then
        call to_screen(" Reaction coordinate calculation")
        if (.not. allocated(react)) then
          allocate(react(3*numat))
          react = 0.d0
        end if
        if (index(keywrd, " HTML") /= 0) call write_path_html(1)
        call drc (react, react)
        iflepo = -1
      else if (index(keywrd, " LOCATE-TS") /= 0) then
        if (.not. use_ref_geo) then
          if (index(keywrd, " LOCATE-TS(SET") == 0) then
            write(line,'(a)')" LOCATE-TS requires GEO_REF to be used"
            call mopend(trim(line))
            return
          end if
        end if
        call Locate_TS_for_Proteins
        if (iflepo == 0) iflepo = -1
      else if (index(keywrd,' NLLSQ') /= 0) then
        write(iw,'(5x,a)')" Transition state refinement using NLLSQ"
        call to_screen(" Transition state refinement using NLLSQ")
        call nllsq ()
      else if (index(keywrd,' SIGMA') /= 0) then
        write(iw,'(5x,a)')" Transition state refinement using SIGMA"
        call to_screen(" Transition state refinement using SIGMA")
        call powsq ()
      else if (nvar == 0) then
        iflepo = 1
        iscf = 1
        last = 1
        i = index(keywrd,' GRAD')
        grad(:nvar) = 0.D0
        numcal = numcal + 1
        tim = second(1)
        call to_screen(" Single point calculation")
        call compfg (xparam, .TRUE., escf, .TRUE., grad, i /= 0)
      else if (index(keywrd,' DFP') + index(keywrd,' FLEPO') + &
        index(keywrd,' BFGS') /= 0 .or. nvar == 1) then
        write(iw,'(/10x,a)')"Geometry optimization using BFGS"
        call to_screen(" Geometry optimization using BFGS")
        call flepo (xparam, nvar, escf)
      else if (index(keywrd,' TS')/=0) then
        write(iw,'(/10x,a)')"Transition state refinement using EF"
        call to_screen(" Transition state refinement using EF")
        call ef (xparam, escf)
      else if (index(keywrd,' LBFGS') /= 0 .or. &
        id == 0 .and. index(keywrd_txt, " GEO_REF") /= 0 .or. &
        nvar > 100 .and. index(keywrd,' EF') == 0) then
        write(iw,'(/10x,a)')"Geometry optimization using L-BFGS"
        call to_screen(" Geometry optimization using L-BFGS")
        call lbfgs (xparam, escf)
      else
        write(iw,'(/10x,a)')"Geometry optimization using EF"
        call to_screen(" Geometry optimization using EF")
        call ef (xparam, escf)
      end if
!
!  Calculation done, now print results
!
      if (moperr) go to 100
      last = 1
      if (iflepo >= 0) then
        call writmo
        if (moperr) go to 100
        if (index(keywrd,' POLAR') /= 0) then
          call polar ()
          if (moperr) go to 100
        endif
        if (index(keywrd,' STATIC') /= 0) then
          numcal = numcal + 1 !  In case POLAR was also used
          call static_polarizability
          if (moperr) go to 100
        endif
        if (index(keywrd,'PMEP') /= 0) call pmep ()
        if (moperr) go to 100
        if (index(keywrd,' ESP') /= 0) then
          call esp ()
          if (moperr) go to 100
        endif
      endif
  100 continue
      tim = tim + second(2)
      if (tim > 1.d7) tim = tim - 1.d7
      inquire(unit = ilog, opened = opend, name = line)
      if (opend) opend = (index(line, log_fn) /= 0)
      if (opend .and. index(keywrd, "LOG") == 0) then
        write(ilog, '(/,'' == MOPAC DONE =='')', iostat = j)
        if (j /= 0) goto 98
        rewind (ilog)
        i = 0
        do
          read(ilog, "(a)", iostat=j) line
          if (j /= 0) exit
          if (index(line, " COVALENT") /= 0)        i = 1
          if (index(line, " charged") /= 0)         i = 1
          if (index(line, " Type    Charge") /= 0)  i = 1
          if (index(line, "TIVE CHARGES") /= 0)     i = 1
          if (index(line, "COMPUTED CHARGE") /= 0)  i = 1
          if (index(line, "UNUSUALLY SHORT") /= 0)  i = 1
          if (index(line, "Original residue") /= 0) i = 1
        end do
        if (i == 0) then
          close(ilog, status = "delete", iostat = j, err = 98)
        else
          close(ilog, status = "keep", iostat = j, err = 98)
        end if
98      continue  
      end if
      if ( .not. gui) then
        if (allocated(p)) deallocate(p)
        if (allocated(react)) deallocate(react)
        inquire (file = end_fn, exist = exists)
        if (exists) then
          open(unit=iend, file=end_fn, status='UNKNOWN', position='asis', iostat=i)
          close(iend, status = 'delete', iostat=i)
        end if
         itemp_1 = ncomments
        if (index(keywrd, " ADD-H PDBOUT") == 0 .or. &
           index(koment, " From PDB file") == 0) ncomments = 0
        if (index(keywrd_txt," GEO_DAT") /= 0) then
          i = index(keywrd_txt," GEO_DAT") + 9
          j = index(keywrd_txt(i + 10:),'" ') + i + 9
          write(line,'(a)')"GEO_DAT="//keywrd_txt(i:j)          
          call l_control(trim(line), len_trim(line), 1)   
        end if
        call delete_MOZYME_arrays()
          go to 10
      end if
!
! Carefully delete all arrays created using "allocate"
!
  101 if ( .not. gui) then
        call setup_mopac_arrays(0,0)
        call delete_MOZYME_arrays()
      end if
      call summary(" ",1)
      if (tim > 1.d7) tim = tim - 1.d7
      write (iw, '(3/,'' TOTAL JOB TIME: '',F16.2,'' SECONDS'')') tim
      write (iw, '(/,'' == MOPAC DONE =='')')
      call fdate (line) 
      write(0,'(//10x,a,/)')"MOPAC Job: """//trim(job_fn)//""" ended normally on "// &
      line(5:10)//", "//trim(line(21:))//", at"//line(11:16)//"."
!
!  Delete files that are definitely not wanted
!
      inquire (file = end_fn, exist = exists)
      if (exists) then
        open(unit=iend, file=end_fn, status='UNKNOWN', position='asis', iostat=i)
        close(iend, status = 'delete', iostat=i)
      end if
      inquire(unit=ir, opened=opend) 
      if (opend) close(ir, status = 'delete', err = 999)
999   jobnam = " "
      inquire(unit = ilog, opened = opend, name = line)
      if (opend) then
        rewind (ilog)
        read(ilog,'(a)', iostat=i)line
        if (i == -1) close(ilog, status = 'delete', iostat=i)
      end if
      return
end subroutine run_mopac
subroutine special
!
!  Use this subroutine for any special work.
!  for example, to print a MOPAC data-set.
!
  use molkst_C, only : jobnam, refkey, line
  use upcase_I
  implicit none
  integer :: iprt = 33, i, j, k, len_key
  open(unit=iprt, file=jobnam(:len_trim(jobnam) - 0)//"_(PM6).arc", status='UNKNOWN', position='asis', iostat = i)
  do i = 1, 6
    if (index(refkey(i), " NULL") /= 0) exit
    line = refkey(i)
    len_key = len_trim(refkey(i))
    call upcase(line, len_key)
!
!  Put all changes in keywords here
!
    j = index(line, " 1SCF")
    if (j /= 0) refkey(i)(j:j + 4) = " "
    j = index(line, " PM6")
    if (j /= 0) refkey(i)(j:j + 3) = " "
    j = index(line, " DENOUT")
    if (j /= 0) refkey(i)(j:j + 7) = " "
    j = index(line, " GRADIENTS")
    if (j /= 0) refkey(i)(j:j + 9) = " "
    j = index(line, " GNORM=")
    if (j == 0) then
      j = index(line, "        ")
      refkey(i)(j:j+8) = " GNORM=4"
    end if
!
!  Remove all extra blank spaces
!
    len_key = len_trim(refkey(i))
    refkey(i)(len_key + 1:len_key + 1) = "@"
    do j = 1, len_key
      do k = 1,10
        if (refkey(i)(j:j + 1) == "  ") refkey(i)(j:) = refkey(i)(j + 1:)
      end do
    end do
    j = index(refkey(i), "@")
    refkey(i)(j:) = " "
    j = index(refkey(i),"     ")
    refkey(i)(j:) = " PM6"
  end do
  write(iprt,"(a)")"  MOPAC2016"
  write(iprt,"(a)")" FINAL GEOMETRY OBTAINED"
  call geout (iprt)
end subroutine special
