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

  subroutine run_mopac
  !dec$ attributes dllexport :: run_mopac
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!
      use common_arrays_C, only : nfirst, nlast, nat, xparam, grad, nw, &
      p, pa, pb, labels, loc, time_start, l_atom, coord, txtatm, coorda, &
      txtatm1, cell_ijk, eigs, c, breaks
!
      USE molkst_C, only : gnorm, natoms, numat, nvar, numcal, job_no, nscf, id, &
        escf, iflepo, iscf, keywrd, last, moperr, maxatoms, ncomments, &
        time0, atheat, errtxt, isok, mpack, line, na1, refkey, keywrd_txt, &
        press, voigt, mozyme, step_num, jobnam, nelecs, stress, E_disp, E_hb, E_hh, no_pKa, &
        MM_corrections, lxfac, trunc_1, trunc_2, l_normal_html, &
        sparkle, itemp_1, maxtxt, koment, sz, ss2, keywrd_quoted, &
        nl_atoms, use_ref_geo, prt_coords, pdb_label, step, &
        density, norbs, method_indo, nclose, nopen, backslash, gui, os, git_hash, verson
!
      USE parameters_C, only : tore, ios, iop, iod, eisol, eheat, zs, eheat_sparkles, gss
!
!
      use cosmo_C, only : iseps, useps, lpka, solv_energy, area, fepsi, ediel
!
      USE funcon_C, only : fpc_9
!
      USE maps_C, only : latom, react, rxn_coord
!
      use symmetry_C, only : state_Irred_Rep, name
!
      USE chanel_C, only : ir, iw, iarc, output_fn, end_fn, iend, &
        archive_fn, log, ilog, xyz_fn, job_fn, log_fn
!
      use MOZYME_C, only : rapid, nres, refnuc, uni_res
!
      use meci_C, only : nmos, lab
!
      use elemts_C, only : elemnt, atom_names
!
      USE reimers_C, only: noh, nvl, cc0, nel, norb, norbl, norbh,&
          nshell, filenm, lenf, evalmo, nbt, multci, occfr, vca, vcb
#ifdef GPU
      Use iso_c_binding
      Use mod_vars_cuda, only: lgpu, ngpus, gpu_id
      Use gpu_info
      Use settingGPUcard
#endif
      implicit none
      integer ::  i, j, k, l
      double precision :: eat,  tim, store_fepsi
      logical :: exists, opend, sparkles_available, l_OLDDEN
      double precision, external :: C_triple_bond_C, reada, seconds
      character :: nokey(20)*10
#ifdef _OPENMP
      integer :: num_threads, default_num_threads
      integer, external :: omp_get_max_threads
#endif
#ifdef MKL
      integer :: num_threads
      integer, external :: mkl_get_max_threads
#endif
#ifdef GPU
      logical :: lgpu_ref
      logical(c_bool)    :: hasGpu = .false.
      logical(c_bool)    :: lstat = .false.
      logical(c_bool)    :: hasDouble(6)
      integer(c_int)     :: nDevices
      character*256	     :: gpuName(6)
      integer(c_size_t)  :: totalMem(6)
      logical            :: gpu_ok(6)
      character*3        :: on_off(6)
      integer(c_int), dimension(6)	 :: clockRate, major, minor, name_size
#endif
! set versioning information
#ifdef MOPAC_VERSION_FULL
      verson = MOPAC_VERSION_FULL
#endif
#ifdef MOPAC_OS
      os = MOPAC_OS
#endif
#ifdef MOPAC_GIT_HASH
      git_hash = MOPAC_GIT_HASH
#endif
! parse command-line flags
#ifdef MOPAC_F2003
      do i = 1, command_argument_count()
        call get_command_argument (i, jobnam)
#else
      do i = 1, iargc()
        call getarg (i, jobnam)
#endif
        if (jobnam == '-V' .OR. jobnam == '--version') then
          write(*,"(a)") "MOPAC version "//trim(verson)//" commit "//trim(git_hash)
          stop
        endif
      end do
!------------------------------------------------------------------------
      tore = ios + iop + iod
      call fbx                            ! Factorials and Pascal's triangle (pure constants)
      call fordd                          ! More constants, for use by MNDO-d
      trunc_1 = 7.0d0    ! Beyond 7.0 Angstroms, use exact point-charge
      trunc_2 = 0.22d0   ! Multiplier in Gaussian: exp(-trunc_2*(trunc_1 - Rab)^2)
      fepsi = 0.d0       ! for correct initialization of store_fepsi
!
! Read in all data; put it into a scratch file, "ir"
!
      moperr = .false.
      call getdat(ir,iw)
      call to_screen("To_file: Start of reading in data")
      if (natoms == 0 .or. moperr) return
!
!   CLOSE UNIT IW IN CASE IT WAS ALREADY PRE-ASSIGNED
!
      close(iw)
!
! WARNING - Replace with standard name, sometime.
!
      i = index(jobnam,' ') - 1
      filenm = trim(jobnam)
      lenf = i+1
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
      tim = seconds(1)
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
      if (job_no > 1) then
        backspace(ir)
        read(ir,'(a)', iostat = i) line
        if (i == 0) then
          if (line /= " ") then
            call upcase(line, len_trim(line))
            if (index(line, "NEXT") /= 0) backspace(ir)
          end if
        end if
      end if
      step_num = step_num + 1  ! New electronic structure, therefore increment step_num
      moperr = .FALSE.
      name = " "
      escf   = 0.d0
      ediel  = 0.d0
      gnorm  = 0.D0
      press  = 0.d0
      voigt  = 0.d0
      E_disp = 0.d0
      E_hb   = 0.d0
      E_hh   = 0.d0
      solv_energy = 0.d0
      sz = 0.d0
      ss2 = 0.d0
      step = 0.d0
      density = 0.d0
      area = 0.d0
      store_fepsi = fepsi
      fepsi = 0.d0
      refnuc = 0.d0
      nres = 0
      nscf = 0
      nmos = 0
      na1 = 0
      norbs = 0
      lab = 0
      lpka = .false.
      stress = 0.d0
      no_pKa = 0
      cell_ijk = 0
      uni_res = 0
      id = 0
      iflepo = 0
      time0 = seconds(1)
      MM_corrections = .false.
      nelecs = 0
      pdb_label = .false.
      l_normal_html = .true.
      state_Irred_Rep = " "
      if (job_no > 1) then
        i = index(keywrd, " BIGCYCL")
        if (i /= 0 .and. index(keywrd,' DRC') == 0) then
          i = nint(reada(keywrd, i)) + 1
          if (job_no < i) then
            fepsi = store_fepsi
            goto 90
          end if
        end if
      end if
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
      if (numcal > 1) call to_screen("To_file: Leaving MOPAC")
!
!    Read in all the data for the current job
!
      i = numcal
      call readmo
!
! Check to see if an old density matrix exists
!
      line = trim(job_fn)
      i = len_trim(line)
      if (i > 4) then
        j = index(line(i - 4:), ".")
        if (j /= 0) i = i - 6 + j
      end if
      inquire(file=line(:i)//".den", exist=l_OLDDEN)
90      if (moperr .and. numcal == 1 .and. natoms > 1) goto 101
      if (moperr .and. numcal == 1 .and. index(keywrd_txt," GEO_DAT") == 0) goto 100
      if (moperr) goto 101
! Adjust maximum number of threads using the OpenMP API
#ifdef _OPENMP
      if (numcal == 1) default_num_threads = omp_get_max_threads()
      i = index(keywrd, " THREADS")
      if (i > 0) then
        num_threads = nint(reada(keywrd, i))
      else
        num_threads = default_num_threads
      end if
      call omp_set_num_threads(num_threads)
#endif
      if (numcal == 1) then
#ifdef MKL
        num_threads = min(mkl_get_max_threads(), 20)
        i = index(keywrd, " THREADS")
        if (i > 0) then
          i = nint(reada(keywrd, i))
          num_threads = min(max(1,i), num_threads)
        end if
        call mkl_set_num_threads(num_threads)
#endif
#ifdef GPU
        gpuName(1:6) = '' ; name_size(1:6) = 0 ; totalMem(1:6) = 0 ; clockRate(1:6) = 0
        hasDouble(1:6) = .false. ; gpu_ok(1:6) = .false.
        clockRate(1:6) = 0 ; major(1:6) = 0 ; minor(1:6) = 0; on_off(1:6) = 'OFF'
        call gpuInfo(hasGpu, hasDouble, nDevices, gpuName,name_size, totalMem, &
                  & clockRate, major, minor)
        lgpu = .false.
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
            end if
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
              end if
            end if
          end if
          if (l == 0) then   ! Select GPU automatically
            do i = 1, nDevices
             if (gpu_ok(i)) then
                on_off(i) = 'ON '
                gpu_id = i - 1
                call setGPU(gpu_id, lstat)
                if (.not. lstat) then
                  write (6,*) 'Problem to set GPU card ID = ', gpu_id
                  stop
                end if
                exit
              end if
            end do
          end if
        else
          nDevices = 0
          ngpus = 0
        end if
!
!  For small systems, using a GPU takes longer than not using a GPU,
!  so do not use a GPU for small systems.  The lower limit, 100, is just a guess.
!
        lgpu = (lgpu_ref .and. natoms > 100) ! Warning - there are problems with UHF calculations on small systems
#endif
      end if
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
      if (index(keywrd_quoted,' EXTERNAL=') + index(keywrd,' EXTERNAL') /= 0) call datin (ir, iw)
      if (moperr) go to 100
      sparkle = (index(keywrd, " SPARKL") /= 0)
!
!  Check to see if SPARKLES are needed and, if need, can they be used.
!
      sparkles_available = .true.
!
! Are SPARKLES needed?
!
      do i = 1, natoms
        if (labels(i) > 83) cycle
        if  (zs(labels(i)) < 0.1d0) then
!
! Can SPARKLES be used?
!
          sparkles_available = (sparkles_available .and. (gss(labels(i)) > 0.1d0))
          if (.not. sparkle) then
            line = " "
            do j = 1, 10
              if (atom_names(labels(i))(j:j) /= " ") exit
            end do
            write (line, '(A,a)') ' Data are not available for ', atom_names(labels(i))(j:)//"."
            if (index(keywrd, " 0SCF") == 0) then
              if (sparkles_available) write(iw,*)" (Parameters are available if SPARKLE is used)"
              call mopend(trim(line))
              goto 100
            end if
          end if
        end if
      end do
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
            if (line(i:i) == "/" .or. line(i:i) == backslash) exit
          end do
        end if
        open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
        rewind iarc
        if (index(keywrd,' 0SCF') /= 0 .and. index(keywrd, " MINI") /= 0 .and. nl_atoms > 0) then
          open(unit=l, file=xyz_fn)
          write(l,"(i6,a)") nl_atoms," "
          write(l,*)"POINT "
          do i = 1, natoms
            if (l_atom(i)) write(l,"(3x,a2,3f15.5)")elemnt(labels(i)), (coord(j,i),j=1,3)
          end do
        end if
      end if
      call delete_ref_key("RESIDUES", len_trim("RESIDUES"), ' ', 1)
      call delete_ref_key("XENO", len_trim("XENO"), ' ', 1)
      call output_rama()
      if (maxtxt == 0 .and. index(keywrd, " RESIDUES") /= 0) then
        call geochk()
        if (moperr) goto 100
      end if
      if ( index(keywrd," PDBOUT") /= 0 .and. maxtxt < 26 .and. index(keywrd," RESID") == 0) then
        if (maxtxt == 0) then
          maxtxt = 26
          do i = 1, natoms
            write(txtatm(i),'(a,i5,1x,a,a)')"HETATM",i, elemnt(labels(i)),"   HET A   1"
          end do
          txtatm1(:natoms) = txtatm(:natoms)
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
      if (mozyme) then
        if (index(keywrd, " PREC") /= 0) then
          call l_control("PREC", len_trim("PREC"), -1)
          if (index(keywrd, " LET") == 0) call l_control("LET", len_trim("LET"), 1)
        end if
      end if
      if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) nelecs = 0
      if (index(keywrd, " ADD-H") /= 0 ) call l_control("NEWPDB", len_trim("NEWPDB"), 1)
      if (index(keywrd,' 0SCF') + index(keywrd, " RESEQ") + index(keywrd, " NEWPDB") + &
        index(keywrd, " ADD-H") + index(keywrd, " SITE=") /= 0 ) then
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
!        xparam(1) = -1.D0
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
        if (prt_coords) call geout (iw)
        if (index(keywrd,' AIGOUT') /= 0) then
          write (iw, '(2/,A)') '  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT'
          call wrttxt (iw)
          call geoutg (iw)
          write (iarc, '(2/,A)') '  GEOMETRY IN GAUSSIAN Z-MATRIX FORMAT'
          call wrttxt (iarc)
          call geoutg (iarc)
        else if (mozyme .or. index(keywrd, " SITE=") + index(keywrd, " ADD-H") /= 0  .or. &
          (index(keywrd," PDBOUT") + index(keywrd," RESEQ") + index(keywrd," NEWPDB") + &
          index(keywrd," RESID") /= 0)) then
          i = size(coorda)
          j = size(coord)
          if (i /= j) then
            deallocate (coorda)
            allocate(coorda(3,natoms))
            coorda = coord
          end if
!
!  Check the format of hydrogen atoms if SITE is used.  If it is the new formatmake sure that NEWPDB is present
!
          if (index(keywrd, " SITE=") /= 0 .and. index(keywrd," NEWPDB") == 0) then
            j = 0
            k = 0
            do i = 1, min(numat,100)
              if (index(txtatm(i), " 1H") /= 0) j = j + 1
              if (index(txtatm(i), " HG13") /= 0) k = k + 1
              line=txtatm(i)(12:15)
            end do
            if (k > 0 .and. j == 0) call l_control("NEWPDB", len("NEWPDB"), 1)
          end if
          if (index(keywrd, " ADD-H") + index(keywrd, " SITE=") == 0) call geochk()
          if (index(keywrd, " ADD-H") == 0 .and. index(keywrd, " SITE=") == 0 .and. index(keywrd," RESEQ") /= 0) then
            moperr = .false.
            goto 10
          end if
          if (index(keywrd, " SITE=") + index(keywrd, " ADD-H") /= 0) then
            if (index(keywrd, " RESEQ") == 0 ) call l_control("Move", len("Move"), 1)
            moperr = .false.
            if (index(keywrd, " ADD-H") /= 0) then
              call store_and_restore_Tv("STORE")
              call add_hydrogen_atoms()
              if (moperr) then
                inquire(unit=iarc, opened=opend) 
                if (opend) close (iarc, status="DELETE") 
                go to 101
              end if
              call move_hydrogen_atoms
              call store_and_restore_Tv("RESTORE")
              call lewis(.false.)
              if (moperr) then
                inquire(unit=iarc, opened=opend)
                if (opend) close (iarc, status="DELETE")
                go to 101
              end if
              if (index(keywrd, " NORESEQ") /= 0) call update_txtatm(.true., .false.)
              call l_control("ADD-H", len("ADD-H"), -1)
              numat = natoms - id
              numcal = numcal + 1
              call mopend("ADD-H: SYSTEM HAS BEEN HYDROGENATED")
              moperr = .false.
            end if
            call l_control("0SCF", len("0SCF"), 1)
!
! Force the TER's to be re-calculated
!
            breaks(1) = -300
            call geochk()
            if (moperr) goto 101
            if (size(coorda) >= size(coord)) then
              coorda(:,:numat) = coord(:,:numat)
              txtatm1(:numat) = txtatm(:numat)
            end if
            moperr = .false.
            i = index(refkey(1), "ADD-H")
            if (i /= 0) refkey(1) = refkey(1)(:i - 1)//refkey(1)(i + 5:)
          else
            if (index(keywrd, " NEWPDB") /= 0) call update_txtatm(.true., .false.)
          end if
          if (index(keywrd, " SITE=") + index(keywrd, " ADD-H") /= 0 .and. &
            index(keywrd," RESEQ") + index(keywrd," RESID") + index(keywrd, " NEWPDB") /= 0) &
            call update_txtatm(.true., .false.)         !  Now that geometry checks are done, switch to input labels
          if (index(keywrd, " NEWPDB") /= 0) call PDB3()
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
                  if (line(i:i) == "/" .or. line(i:i) == backslash) exit
                end do
              end if
            open(unit=iarc, file=trim(line), status='UNKNOWN', position='asis')
            rewind iarc
            call pdbout(iarc)
          end if
        else
          if (index(keywrd, " ADD-H") /= 0) then
            call store_and_restore_Tv("STORE")
            call add_hydrogen_atoms()
            call store_and_restore_Tv("RESTORE")
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
              write(line,'(a6,i5,a15)')txtatm(i)(:6),i + j - 1,txtatm(i)(12:)
              txtatm(i) = trim(line)
            end do
          end if
          call geout (iarc)
        end if
        if (index(keywrd, " SITE=") + index(keywrd, " ADD-H") + index(keywrd, " 0SCF") + &
          index(keywrd," RESEQ") + index(keywrd," RESID") /= 0) go to 100
        if (nelecs == 0 .or. index(keywrd, " NEWPDB") == 0) goto 100
        close (iarc)
      end if
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
      if (method_indo) then
        if(.not. allocated(occfr)) allocate(occfr(2))
        noh = nopen
        nvl = nclose + 1
        nel(1) = nclose * 2

        norb(1) = nclose
        norbl(1) = 1
        norbh(1) = nclose
        occfr(1) = 2.D0

        if (nopen == nclose) then
          nshell = 1

        else
          nshell = 2

          nel(2) = nopen - nclose
          norb(2) = nopen - nclose
          norbl(2) = norbh(1) + 1
          norbh(2) = nopen
          occfr(2) = 1.D0
          allocate(vca(2,2))
          allocate(vcb(2,2))
! One alpha electron per orbital
          vca(2,2)= 1.D0
          vcb(2,2)= 2.D0
        end if
! Set up CI data
        if (index(keywrd,' CIS') /= 0 .or. index(keywrd,' MRCI') /= 0 .or. &
            index(keywrd,' C.I.') /= 0 .or. index(keywrd,' C.A.S.') /= 0) then
          multci = nopen - nclose + 1
          if (index(keywrd,' SING') /= 0) multci = 1
          if (index(keywrd,' DOUB') /= 0) multci = 2
          if (index(keywrd,' TRIP') /= 0) multci = 3
          if (index(keywrd,' QUAR') /= 0) multci = 4
          if (allocated(evalmo)) deallocate(evalmo)
          if (allocated(nbt))    deallocate(nbt)
          allocate(evalmo(norbs))
          allocate(nbt(norbs))
        end if
      end if
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
      call l_control("GEO_DAT", len_trim("GEO_DAT"), -1)
      call l_control("GEO_REF", len_trim("GEO_REF"), -1)
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
      if (index(keywrd,' 1SCF') /= 0 .or. method_indo) then
        if (method_indo .and. index(keywrd,' 1SCF') == 0) then
          write (iw,*) "WARNING: INDO only performs single-point calculations"
        end if
        iflepo = 1
        iscf = 1
        last = 1
        i = index(keywrd,' GRAD')
        grad(:nvar) = 0.D0
        numcal = numcal + 1
        tim = seconds(1)
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
        keywrd = " "
        iflepo = -1
      else if (index(keywrd,' DRC') + index(keywrd,' IRC') /= 0) then
        call to_screen(" Reaction coordinate calculation")
        if (.not. allocated(react)) then
          allocate(react(3*numat))
          react = 0.d0
        end if
        if (index(keywrd, " HTML") /= 0) then
          call write_path_html(1)
        end if
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
        call Locate_TS
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
        tim = seconds(1)
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
      if (method_indo .and. (index(keywrd,' CIS') /= 0 .or. index(keywrd,' MRCI') /= 0 .or.&
            index(keywrd,' C.I.') /= 0 .or. index(keywrd,' C.A.S.') /= 0)) then
        do i=1,norbs
          evalmo(i) = eigs(i)
          do j = 1,norbs
            cc0(j,i) = c(i,j)
          end do
        end do
        do i=1,natoms
          k = 0
          do j=nfirst(i),nlast(i)
            nbt(j) = k
            k = k+1
          end do
        end do
        call rci()
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
        end if
        if (index(keywrd,' STATIC') /= 0) then
          numcal = numcal + 1 !  In case POLAR was also used
          call static_polarizability
          if (moperr) go to 100
        end if
        if (index(keywrd,'PMEP') /= 0) call pmep ()
        if (moperr) go to 100
        if (index(keywrd,' ESP') /= 0) then
          call esp ()
          if (moperr) go to 100
        end if
      end if
  100 continue
      tim = tim + seconds(2)
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
!
! Delete density matrix if it was made by MOZYME
!
        if (.not. l_OLDDEN .and. index(keywrd, " NEWDEN") == 0) then
          j = len_trim(end_fn)
          inquire (file = end_fn(:j - 3)//"den", exist = exists)
          if (exists) then
            open(unit = iend, file = end_fn(:j - 3)//"den", status='OLD', iostat=i)
            if (i == 0) close(iend, status = 'delete', iostat=i)
          end if
        end if
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
      write(*,'(//10x,a,/)')"MOPAC Job: """//trim(job_fn)//""" ended normally on "// &
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
  write(iprt,"(a)")"  MOPAC"
  write(iprt,"(a)")" FINAL GEOMETRY OBTAINED"
  call geout (iprt)
end subroutine special
