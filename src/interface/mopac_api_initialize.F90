! Molecular Orbital PACkage (MOPAC)
! Copyright 2021 Virginia Polytechnic Institute and State University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

submodule (mopac_api:mopac_api_operations) mopac_api_initialize
  use chanel_C, only : output_fn, & ! name of output file
    iw ! file handle for main output file
  use conref_C, only : fpcref ! values of fundamental constants (data)
  use Common_arrays_C, only : time_start, & ! starting time of the MOPAC calculation
    nat, & ! atomic numbers of real atoms
    labels, & ! atomic numbers of real atoms, labels for dummy atoms (Tv = 107)
    atmass, & ! atomic masses of atoms
    xparam, & ! values of coordinates undergoing optimization
    geo, & ! raw coordinates of atoms
    loc, & ! indices of atoms and coordinates marked for optimization
    lopt, & ! optimization flags for Cartesian coordinates of atoms
    na, nb, nc, & ! internal coordinate connectivity information
    breaks, chains, cell_ijk, nw, nfirst, nlast
  use cosmo_C, only : iseps, & ! flag for use of COSMO model
    noeps, useps, lpka, solv_energy, area, fepsi, ediel, nspa
  use funcon_C, only : fpc, fpc_9 ! fundamental constants used in the MOPAC calculation
  use maps_C, only : latom, lparam, lpara1, latom1, lpara2, latom2, rxn_coord
  use meci_C, only : nmos, lab
  use molkst_C, only : use_disk, & ! variable to turn off disk usage for MOPAC API
    keywrd, & ! keyword string to adjust MOPAC behavior
    trunc_1, & ! cutoff radius for switch to point-charge electrostatics
    trunc_2, & ! exponent in smooth cutoff for point-charge electrostatics
    time0, & ! starting time of the MOPAC calculation
    tleft, & ! amount of time left for the calculation to finish
    natoms, maxatoms, & ! number of atoms
    numat, & ! number of real atoms
    mozyme, & ! logical flag for MOZYME calculations
    methods, methods_keys, & ! logical flags and names of semiempirical models
    l_feather, & ! smooth transition from NDDO to point-charge
    method_am1, method_pm3, method_rm1, & ! logical flag aliases for specific models
    id, & ! number of translation/lattice vectors
    nvar, & ! number of coordinates to be optimized
    numcal, job_no, step_num, & ! MOPAC's persistent state variables
    moperr, & ! error status
    isok, MM_corrections, pdb_label, l_normal_html, lxfac, &
    nscf, na1, norbs, no_pKa, id, iflepo, nelecs, ndep, &
    escf, gnorm, E_disp, E_hb, E_hh, sz, ss2, step, density, press, stress, voigt, &
    sparkle, &  ! logical flag for the presence of sparkles
    atheat, & ! total atomic heat of formation
    keywrd_txt, keywrd_quoted, koment, title, refkey, allkey ! other text extracted from input files
  use MOZYME_C, only : nres, uni_res, refnuc
  use parameters_C, only : tore, & ! number of valence electrons per element
    ios, iop, iod, & ! number of s, p, and d valence electrons per element (data)
    ams, & ! standard atomic masses (data)
    zs, & ! s-type STO exponents of elements
    eheat, eheat_sparkles, eisol ! heat of formation reference values (data)
  use symmetry_C, only : state_Irred_Rep, name
  implicit none

contains

  ! This subroutine is replicating the setup steps of run_mopac
  ! for running the limited functionality supported by the MOPAC API.
  ! It is likely that not every step here is strictly necessary, but it is
  ! difficult to assess what is and isn't necessary without a very substantial
  ! investment of development time. Similarly, the keyword line is set to represent
  ! the calculation that the API is trying to perform, but some of the keyword-dependent
  ! logic has been stripped out of this restricted initialization.
  module subroutine mopac_initialize(system)
    type(mopac_system), intent(in) :: system
    double precision, external :: seconds
    character(100) :: num2str
    integer :: i, j, nelectron, status
    double precision :: eat
    integer(c_int), pointer :: atom(:)
    real(c_double), pointer :: coord(:)
    real(c_double), pointer :: lattice(:)

    ! validity checks of data in system
    if (system%natom <= 0) call mopend("Invalid number of atoms")
    if (system%natom_move < 0 .or. system%natom_move > system%natom) &
      call mopend("Invalid number of moveable atoms")
    if (system%epsilon /= 1.d0 .and. system%nlattice > 0) &
      call mopend("COSMO solvent not available for periodic systems")
    if (system%natom > 0 .and. .not. c_associated(system%atom)) &
      call mopend("Array of atomic numbers not found")
    if (system%natom > 0 .and. .not. c_associated(system%coord)) &
      call mopend("Array of atomic coordinates is too small")
    if (system%nlattice < 0 .or. system%nlattice > 3) call mopend("Invalid number of lattice vectors")
    if (system%nlattice_move < 0 .or. system%nlattice_move > system%nlattice) &
      call mopend("Invalid number of moveable lattice vectors")
    if (system%nlattice > 0 .and. .not. c_associated(system%lattice)) &
      call mopend("List of lattice vectors is too small")
    if (system%tolerance <= 0.d0) call mopend("Relative numerical tolerance must be a positive number")
    if (system%max_time <= 0) call mopend("Time limit must be a positive number")
    if (system%nlattice_move > 0 .and. index(keywrd," FORCETS") /= 0) &
      call mopend("Lattice vectors cannot move during vibrational calculations")
    if (moperr) return

    ! convert C pointers
    call c_f_pointer(system%atom, atom, [system%natom])
    call c_f_pointer(system%coord, coord, [3*system%natom])
    if (system%nlattice > 0) call c_f_pointer(system%lattice, lattice, [3*system%nlattice])

    use_disk = .false.
    ! load ios, iop, and iod data into tore
    tore = ios + iop + iod
    ! factorials and Pascal's triangle (mndod_C)
    call fbx
    ! more constants, for use by MNDO-d (mndod_C)
    call fordd
    ! point-charge switching radius (Angstroms) assumed by all MOPAC models
    trunc_1 = 7.0d0
    ! exponent in point-charge switching weight: exp(-trunc_2*(trunc_1 - Rab)^2)
    trunc_2 = 0.22d0
    ! initialize CODATA fundamental constants
    fpc(:) = fpcref(1,:)
    ! assemble virtual keyword line
    keywrd = trim(keywrd) // " LET NOSYM GEO-OK"
    write(num2str,'(i6)') system%charge
    keywrd = trim(keywrd) // " CHARGE=" // adjustl(num2str)
    if (system%pressure /= 0.d0) then
      write(num2str,'(f12.6)') system%pressure
      keywrd = trim(keywrd) // " P=" // trim(adjustl(num2str)) // "GPa"
    end if
    write(num2str,'(f12.6)') system%tolerance
    keywrd = trim(keywrd) // " GNORM=" // adjustl(num2str)
    keywrd = trim(keywrd) // " RELSCF=" // adjustl(num2str)
    write(num2str,'(i10)') system%max_time
    keywrd = trim(keywrd) // " T=" // adjustl(num2str)
    ! open dummy output file
#ifdef WIN32
    output_fn = 'NUL'
#else
  output_fn = '/dev/null'
#endif
    close(iw)
    open(unit=iw, file=output_fn, iostat=status)
    if (status /= 0) then
      call mopend("Failed to open NULL file in MOPAC_INITIALIZE")
      return
    end if
    ! initialize job timers
    time0 = seconds(1)
    call date_and_time(VALUES=time_start)
    tleft = system%max_time
    ! set number of atoms
    natoms = system%natom + system%nlattice
    maxatoms = natoms
    numat = system%natom
    ! set MOZYME flag
    mozyme = (index(keywrd," MOZ") /= 0)
    ! initialize common workspaces (Common_arrays_C)
    call setup_mopac_arrays(natoms, 1)
    if (.not. allocated(lopt)) then
      allocate(lopt(3,maxatoms), stat=status)
      if (status /= 0) then
        call mopend("Failed to allocate memory in MOPAC_INITIALIZE")
        return
      end if
    end if
    ! set semiempirical model
    methods = .false.
    select case (system%model)
      case (0) ! PM7
        i = 14
      case (1) ! PM6-D3H4
        i = 9
        ! PM6-D3H4 also needs the PM6 flag for PM6-x model family checks
        methods(6) = .true.
      case (2) ! PM6-ORG
        i = 18
      case (3) ! PM6
        i = 6
      case (4) ! AM1
        i = 2
      case (5) ! RM1
        i = 4
      case default
        call mopend("Unknown semiempirical model requested")
    end select
    methods(i) = .true.
    keywrd = trim(keywrd) // methods_keys(i)
    ! feathered interactions for PM6-ORG and PM7, or for MOZYME calculations of PM6-based models
    l_feather = (methods(14) .or. methods(18) .or. index(keywrd, " MOZ") /= 0 .and. (index(keywrd, " PM6") /= 0))
    ! check for use of COSMO model
    if (system%epsilon /= 1.d0) then
      iseps = .true.
      useps = .true.
      noeps = .false.
      write(num2str,'(f6.2)') system%epsilon
      keywrd = trim(keywrd) // " EPS=" // adjustl(num2str)
      fepsi = (system%epsilon-1.d0) / (system%epsilon+0.5d0)
    else
      iseps = .false.
      useps = .false.
      noeps = .true.
      fepsi = 0.d0
    end if
    ! insert geometry information into MOPAC data structures
    id = system%nlattice
    nat(:system%natom) = atom(:system%natom)
    labels(:system%natom) = atom(:system%natom)
    geo(:,:system%natom) = reshape(coord,[3, system%natom])
    if (id > 0) then
      nat(system%natom+1:system%natom+id) = 107
      labels(system%natom+1:system%natom+id) = 107
      geo(:,system%natom+1:system%natom+id) = reshape(lattice,[3, id])
    end if
    atmass(:natoms) = ams(labels(:))
    ! set optimization flags & xparam
    nvar = 3*system%natom_move + 3*system%nlattice_move
    lopt(:,:system%natom_move) = 1
    lopt(:,system%natom_move+1:) = 0
    xparam(:3*system%natom_move) = coord(:3*system%natom_move)
    do i=1, system%natom_move
      do j=1, 3
        loc(1,3*(i-1)+j) = i
        loc(2,3*(i-1)+j) = j
      end do
    end do
    if (system%nlattice_move > 0) then
      lopt(:,numat+1:numat+system%nlattice_move) = 1
      xparam(3*system%natom_move+1:) = lattice(:3*system%nlattice_move)
      do i=1, system%nlattice_move
        do j=1, 3
          loc(1,3*system%natom_move+3*(i-1)+j) = numat+i
          loc(2,3*system%natom_move+3*(i-1)+j) = j
        end do
      end do  
    end if
    ! update MOPAC state variables that do not reset
    numcal = numcal + 1
    job_no = job_no + 1
    step_num = step_num + 1
    ! initialize variables that need default values (molkst_C)
    moperr = .false.
    isok = .true.
    MM_corrections = .false.
    pdb_label = .false.
    l_normal_html = .true.
    lxfac = .false.
    nscf = 0
    na1 = 0
    norbs = 0
    no_pKa = 0
    iflepo = 0
    nelecs = 0
    ndep = 0
    escf = 0.d0
    gnorm = 0.D0
    E_disp = 0.d0
    E_hb = 0.d0
    E_hh = 0.d0
    sz = 0.d0
    ss2 = 0.d0
    step = 0.d0
    density = 0.d0
    press = 0.d0
    stress = 0.d0
    voigt = 0.d0
    keywrd_txt = ' '
    keywrd_quoted = ' '
    koment = ' '
    title = ' '
    refkey = ' '
    allkey = ' '
    ! initialize variables that need default values (Common_arrays_C)
    na = 0
    nb = 0
    nc = 0
    cell_ijk = 0
    chains = " "
    breaks(1) = -300
    ! initialize variables that need default values (cosmo_C)
    lpka = .false.
    solv_energy = 0.d0
    area = 0.d0
    ediel = 0.d0
    nspa = 42
    ! initialize variables that need default values (maps_C)
    latom = 0
    lparam = 0
    lpara1 = 0
    latom1 = 0
    lpara2 = 0
    latom2 = 0
    rxn_coord = 1.d9
    ! initialize variables that need default values (meci_C)
    nmos = 0
    lab = 0
    ! initialize variables that need default values (MOZYME_C)
    nres = 0
    uni_res = 0
    refnuc = 0.d0
    ! initialize variables that need default values (symmetry_C)
    name = " "
    state_Irred_Rep = " "
    ! load model parameters
    call switch
    ! adjust number of valence electrons for sparkles
    sparkle = .false.
    do i = 1, numat
      if  (nat(i) > 56 .and. nat(i) < 72 .and. zs(nat(i)) < 0.1d0) sparkle = .true.
    end do
    do i = 57,71
      if (zs(i) < 0.1d0) tore(i) = 3.d0
    end do
    ! set the spin multiplicity
    nelectron = sum(tore(nat(:numat))) - system%charge
    write(num2str,'(f6.1)') system%spin + mod(nelectron,2)*0.5d0
    keywrd = trim(keywrd) // " MS=" // adjustl(num2str)
    ! initialization steps performed by moldat
    call moldat (0)
    ! evaluate derived parameters of the model
    call calpar
    ! initialize system-specific MOPAC workspaces (Common_arrays_C)
    call setup_mopac_arrays(1,2)
    if (allocated(nw)) deallocate(nw)
    allocate(nw(numat), stat=status)
    if (status /= 0) then
      call mopend("Failed to allocate memory in MOPAC_INITIALIZE")
      return
    end if
    j = 1
    do i = 1, numat
      nw(i) = j
      j = j + ((nlast(i)-nfirst(i)+1)*(nlast(i)-nfirst(i)+2))/2
    end do
    ! calculate the atomic energy, including sparkles for Ce-Yb
    atheat = 0.d0
    do i = 1, numat
      if  (nat(i) > 56 .and. nat(i) < 72 .and. zs(nat(i)) < 0.1d0) then
        atheat = atheat + eheat_sparkles(nat(i))
      else
        atheat = atheat + eheat(nat(i))
      end if
    end do
    eat = sum(eisol(nat(:numat)))
    atheat = atheat - eat*fpc_9
    ! setup MOZYME calculations
    if (mozyme) call set_up_MOZYME_arrays()
  end subroutine mopac_initialize

end submodule mopac_api_initialize
