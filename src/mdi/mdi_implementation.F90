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

MODULE MDI_IMPLEMENTATION

  USE, INTRINSIC :: ISO_C_binding

  ! MDI library dependencies
  use MDI, only : MDI_Init, MDI_Accept_communicator, MDI_Send, &
    MDI_Recv_command, MDI_Recv, MDI_COMMAND_LENGTH, MDI_LABEL_LENGTH, &
    MDI_DOUBLE, MDI_CHAR, MDI_INT, MDI_Register_command, MDI_Register_node, &
    MDI_Conversion_factor, MDI_Get_role, MDI_ENGINE, MDI_Plugin_get_argc, &
    MDI_Plugin_get_arg, MDI_Set_plugin_state, MDI_Set_execute_command_func

  ! MOPAC data exposed to the MDI interface:
  USE chanel_C, only : iw

  use molkst_C, only : numat, & ! number of real atoms
    id, & ! number of translation/lattice vectors
    nvar, & ! number of coordinates to be optimized
    msdel, & ! magnetic component of spin
    escf, & ! heat of formation
    nelecs, & ! number of electrons
    voigt, & ! Voigt stress tensor (xx, yy, zz, yz, xz, xy)
    jobnam, & ! path to input file
    run ! run mode

  use Common_arrays_C, only : geo, & ! raw coordinates of atoms (highest priority for unrelaxed coordinates)
    xparam, & ! values of coordinates undergoing optimization (highest priority for relaxed coordinates)
    loc, & ! indices of atoms and coordinates marked for optimization
    p, & ! total density matrix
    q, & ! partial charges
    nat, & ! atomic numbers of real atoms
    atmass, & ! atomic masses
    grad, & ! gradients of heat
    txtatm1 ! original atom discriptor

  use parameters_C, only : tore ! number of valence electrons per element

  implicit none

  ! MDI Communicator to the driver
  INTEGER :: comm

  ! Status flags for MDI
  LOGICAL :: use_mdi = .false.
  LOGICAL :: open_mdi = .false.
  LOGICAL :: close_mdi = .false.
  LOGICAL :: terminate_flag
  LOGICAL :: recompute_flag

CONTAINS

  FUNCTION MDI_Plugin_open_mopac(plugin_state) bind ( C, name="MDI_Plugin_open_mopac" )
    TYPE(C_PTR), VALUE :: plugin_state
    INTEGER :: MDI_Plugin_open_mopac
    INTEGER :: ierr, argc
    CHARACTER(len=240) :: input_file
    PROCEDURE(execute_command), POINTER :: generic_command => null()
    TYPE(C_PTR)                         :: class_obj
    MDI_Plugin_open_mopac = 0
    CALL MDI_Set_plugin_state(plugin_state, ierr)
    generic_command => execute_command
    CALL MDI_Set_execute_command_func(c_funloc(generic_command), class_obj, ierr)
    CALL MDI_Plugin_get_argc(argc, ierr)
    IF (argc > 0) THEN
      CALL MDI_Plugin_get_arg(0, input_file, ierr)
    ELSE
      MDI_Plugin_open_mopac = 1
      WRITE(*,*) "MDI ERROR: No input file provided for MOPAC engine"
      RETURN
    END IF
    run = 2
    jobnam = trim(input_file)
    use_mdi = .true.
    open_mdi = .true.
    CALL run_mopac
    open_mdi = .false.
  END FUNCTION MDI_Plugin_open_mopac

  FUNCTION MDI_Plugin_close_mopac() bind ( C, name="MDI_Plugin_close_mopac" )
    INTEGER :: MDI_Plugin_close_mopac
    MDI_Plugin_close_mopac = 0
    close_mdi = .true.
    CALL run_mopac
    close_mdi = .false.
    use_mdi = .false.
    jobnam = ' '
    run = 1
  END FUNCTION MDI_Plugin_close_mopac

  FUNCTION MDI_Plugin_launch_mopac(plugin_state) bind ( C, name="MDI_Plugin_launch_mopac" )
    TYPE(C_PTR), VALUE :: plugin_state
    INTEGER :: MDI_Plugin_launch_mopac
    INTEGER :: ierr, argc
    CHARACTER(len=240) :: input_file
    PROCEDURE(execute_command), POINTER :: generic_command => null()
    TYPE(C_PTR)                         :: class_obj
    MDI_Plugin_launch_mopac = 0
    CALL MDI_Set_plugin_state(plugin_state, ierr)
    generic_command => execute_command
    CALL MDI_Set_execute_command_func(c_funloc(generic_command), class_obj, ierr)
    CALL MDI_Plugin_get_argc(argc, ierr)
    IF (argc > 0) THEN
      CALL MDI_Plugin_get_arg(0, input_file, ierr)
    ELSE
      MDI_Plugin_launch_mopac = 1
      WRITE(*,*) "MDI ERROR: No input file provided for MOPAC engine"
      RETURN
    END IF
    run = 2
    jobnam = trim(input_file)
    use_mdi = .true. 
    CALL run_mopac
    use_mdi = .false.
    jobnam = ' '
    run = 1
  END FUNCTION MDI_Plugin_launch_mopac

  SUBROUTINE initialize_mdi()
    INTEGER :: i, ierr, role
    CHARACTER(len=1024) :: mdi_options

    ! Check for CLI mode (-mdi command-line option)
    do i = 1, iargc()
      call getarg (i, mdi_options)
      if (mdi_options == '-mdi' .OR. mdi_options == '--mdi') then
        if ( i .eq. iargc() ) then
          WRITE(iw,*) "MDI ERROR: No argument was provided for the -mdi command-line option"
          stop
        else
          call getarg (i+1, mdi_options)
          write(iw,*) " Received MDI command-line options: ",TRIM(mdi_options)
          call MDI_Init( mdi_options, ierr)
          use_mdi = .true.
        endif
      endif
    end do

    IF ( use_mdi ) THEN
      ! Confirm that the code is being run as an ENGINE
      call MDI_Get_role(role, ierr)
      IF ( role /= MDI_ENGINE ) THEN
        WRITE(iw,*) "MDI ERROR: Must run MOPAC as an ENGINE"
      END IF

      ! Set control flags
      recompute_flag = .true.
      terminate_flag = .false.

      ! Register supported MDI commands
      CALL MDI_Register_node("@DEFAULT", ierr)
      CALL MDI_Register_command("@DEFAULT", "EXIT", ierr)
      CALL MDI_Register_command("@DEFAULT", "<@", ierr)
      CALL MDI_Register_command("@DEFAULT", ">CELL", ierr)
      CALL MDI_Register_command("@DEFAULT", "<CELL", ierr)
      CALL MDI_Register_command("@DEFAULT", "<CELL_DISPL", ierr)
      CALL MDI_Register_command("@DEFAULT", "<CHARGES", ierr)
      CALL MDI_Register_command("@DEFAULT", ">COORDS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<COORDS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<DIMENSIONS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<ELEC_MULT", ierr)
      CALL MDI_Register_command("@DEFAULT", "<ELEMENTS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<ENERGY", ierr)
      CALL MDI_Register_command("@DEFAULT", "<FORCES", ierr)
      CALL MDI_Register_command("@DEFAULT", "<LABELS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<MASSES", ierr)
      CALL MDI_Register_command("@DEFAULT", "<NATOMS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<STRESS", ierr)
      CALL MDI_Register_command("@DEFAULT", "<TOTCHARGE", ierr)

      ! Connect to the driver
      call MDI_Accept_communicator(comm, ierr)
      write(iw,*) " Received MDI communicator: ", comm
    END IF

  END SUBROUTINE initialize_mdi

  SUBROUTINE respond_to_commands()
    CHARACTER(len=:), ALLOCATABLE :: command
    INTEGER                       :: ierr

    TYPE(C_PTR)                   :: class_obj

    ALLOCATE( character(MDI_COMMAND_LENGTH) :: command )

    ! Respond to the driver's commands
    response_loop: DO

      ! Receive a command from the driver and broadcast it to all ranks
      CALL MDI_Recv_command(command, comm, ierr)
      write(iw,*) " Received MDI command: ",TRIM(command)

      ierr = execute_command(command, comm, class_obj)
      IF ( ierr /= 0 ) STOP 1

      IF ( terminate_flag ) EXIT

    END DO response_loop

    DEALLOCATE( command )

  END SUBROUTINE respond_to_commands

  FUNCTION execute_command(command, comm, class_obj)

    CHARACTER(LEN=*), INTENT(IN)  :: command
    INTEGER, INTENT(IN)           :: comm
    TYPE(C_PTR), VALUE            :: class_obj
    INTEGER(KIND=C_INT)           :: execute_command

    INTEGER                       :: i, j, ierr
    INTEGER                       :: natoms
    DOUBLE PRECISION              :: conv
    DOUBLE PRECISION              :: charge
    CHARACTER(len=:), ALLOCATABLE :: char_array
    INTEGER, ALLOCATABLE          :: int_array(:)
    DOUBLE PRECISION, ALLOCATABLE :: real_array(:)

    execute_command = 0
 
    SELECT CASE( TRIM(command) )
    CASE( "EXIT" )
      terminate_flag = .true.

    CASE( "<@" )
      ! Send the name of the current node
      ! For MOPAC, this is always "@DEFAULT"
      CALL MDI_Send("@DEFAULT", MDI_COMMAND_LENGTH, MDI_CHAR, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<CELL_DISPL" )
      ! Send the current displacement of the cell origin
      ! This is always zero for MOPAC
      ALLOCATE( real_array(3) )
      real_array = 0.0d0
      CALL MDI_Send(real_array, 3, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( real_array )

    CASE( "<CELL" )
      ! Send the cell dimensions
      ! This should be in the form of three vectors
      CALL MDI_Conversion_factor("angstrom", "atomic_unit_of_length", conv, ierr)
      if (ierr /= 0) execute_command = 1
      ALLOCATE( real_array(9) )
      real_array = 0.0d0
      do i=1, id
        do j=1, 3
          real_array(3*(i-1)+j) = geo(j,numat+i) * conv
        end do
      end do
      CALL MDI_Send(real_array, 9, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( real_array )

    CASE( ">CELL" )
      ! Receive new cell dimensions
      ! NOTE: This is not reliable if the cell changes are too large,
      !   in which case the summation bounds over unit cells need to be recomputed
      CALL MDI_Conversion_factor("atomic_unit_of_length", "angstrom", conv, ierr)
      if (ierr /= 0) execute_command = 1
      ALLOCATE( real_array(9) )
      CALL MDI_Recv(real_array, 9, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      do i=1, id
        do j=1, 3
          geo(j, numat+i) = real_array(3*(i-1)+j) * conv
        end do
      end do
      DEALLOCATE( real_array )
      ! make sure that optimized coordinates in xparam are synchronized with geo
      do i=MAX(nvar-3*id+1,1), nvar
        if (loc(1,i) > numat) then
          xparam(i) = geo(loc(2,i), loc(1,i))
        end if
      end do
      recompute_flag = .true.

    CASE( "<CHARGES" )
      ! Send the partial charges of each atom
      CALL recompute_as_needed()
      ! extra code for partial charge post-processing
      call chrge (p, q)
      q(:numat) = tore(nat(:numat)) - q(:numat)
      CALL MDI_Send(q, numat, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<COORDS" )
      ! Send the current nuclear coordinates
      CALL MDI_Conversion_factor("angstrom", "atomic_unit_of_length", conv, ierr)
      if (ierr /= 0) execute_command = 1
      ALLOCATE( real_array(3*numat) )  
      do i=1, numat
        do j=1, 3
          real_array(3*(i-1)+j) = geo(j,i) * conv
        end do
      end do
      CALL MDI_Send(real_array, 3 * numat, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( real_array )

    CASE( ">COORDS" )
      ! Receive a new set of nuclear coordinates
      CALL MDI_Conversion_factor("atomic_unit_of_length", "angstrom", conv, ierr)
      if (ierr /= 0) execute_command = 1
      ALLOCATE( real_array(3*numat) )
      CALL MDI_Recv(real_array, 3 * numat, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      do i=1, numat
        do j=1, 3
          geo(j, i) = real_array(3*(i-1)+j) * conv
        end do
      end do
      DEALLOCATE( real_array )
      ! make sure that optimized coordinates in xparam are synchronized with geo
      do i=1, nvar
        if (loc(1,i) <= numat) then
          xparam(i) = geo(loc(2,i), loc(1,i))
        end if
      end do
      recompute_flag = .true.

    CASE( "<DIMENSIONS" )
      ! Send the status of translation vectors
      ALLOCATE( int_array(3) )
      int_array = 1
      if (id > 0) int_array(1:id) = 2
      CALL MDI_Send(int_array, 3, MDI_INT, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( int_array )

    CASE( "<ELEC_MULT" )
      ! Send the current electronic multiplicity
      CALL MDI_Send(msdel+1, 1, MDI_INT, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<ELEMENTS" )
      ! Send the element of each atom
      CALL MDI_Send(nat, numat, MDI_INT, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<ENERGY" )
      ! Calculate and send the energy
      CALL recompute_as_needed()
      CALL MDI_Conversion_factor("kilocalorie_per_mol", "atomic_unit_of_energy", conv, ierr)
      if (ierr /= 0) execute_command = 1
      CALL MDI_Send(escf*conv, 1, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<FORCES" )
      ! Calculate and send the nuclear forces
      CALL recompute_as_needed()
      ALLOCATE( real_array(3*numat) )
      real_array = 0.0d0
      do i=1, nvar
        if (loc(1,i) <= numat) then
          real_array(3*(loc(1,i)-1)+loc(2,i)) = grad(i)
        end if
      end do
      CALL MDI_Conversion_factor("kilocalorie_per_mol", "atomic_unit_of_energy", conv, ierr)
      if (ierr /= 0) execute_command = 1
      real_array = real_array * conv
      CALL MDI_Conversion_factor("angstrom", "atomic_unit_of_length", conv, ierr)
      if (ierr /= 0) execute_command = 1
      real_array = real_array / conv
      CALL MDI_Send(real_array, 3 * numat, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( real_array )

    CASE( "<LABELS" )
      ! Send the atom labels
      ALLOCATE( CHARACTER(len=(MDI_LABEL_LENGTH * natoms)) :: char_array )
      char_array = ""
      do i=1, numat
        char_array(1+(i-1)*MDI_LABEL_LENGTH:i*MDI_LABEL_LENGTH) = trim(txtatm1(i))
      end do
      CALL MDI_Send(char_array, MDI_LABEL_LENGTH * natoms, MDI_CHAR, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( char_array )

    CASE( "<MASSES" )
      ! Send the atomic masses
      CALL MDI_Send(atmass, numat, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<NATOMS" )
      ! Send the current number of atoms in the system
      CALL MDI_Send(numat, 1, MDI_INT, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE( "<STRESS" )
      ! Calculate and send the virial stress tensor
      CALL recompute_as_needed()
      ALLOCATE( real_array(9) )
      ! (xx, yy, zz, yz, xz, xy) -> (xx, xy, xz, yx, yy, yz, zx, zy, zz)
      real_array(1) = voigt(1)
      real_array(2) = voigt(6)
      real_array(3) = voigt(5)
      real_array(4) = voigt(6)
      real_array(5) = voigt(2)
      real_array(6) = voigt(4)
      real_array(7) = voigt(5)
      real_array(8) = voigt(4)
      real_array(9) = voigt(3)
      CALL MDI_Conversion_factor("newton", "atomic_unit_of_force", conv, ierr)
      if (ierr /= 0) execute_command = 1
      real_array = real_array * conv * 1e-9
      CALL MDI_Conversion_factor("meter", "atomic_unit_of_length", conv, ierr)
      if (ierr /= 0) execute_command = 1
      real_array = real_array / (conv*conv)
      CALL MDI_Send(real_array, 9, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1
      DEALLOCATE( real_array )

    CASE( "<TOTCHARGE" )
      ! Send the total charge of the system
      charge = -nelecs
        do i=1, numat
          charge = charge + tore(nat(i))
        end do
      CALL MDI_Send(charge, 1, MDI_DOUBLE, comm, ierr)
      if (ierr /= 0) execute_command = 1

    CASE DEFAULT
      write(iw,*) "Unknown command ", trim(command)," received through MDI"
      execute_command = 1
    END SELECT

  END FUNCTION execute_command

  SUBROUTINE recompute_as_needed()
  if ( recompute_flag ) then
    call compfg (xparam, .true., escf, .true., grad, .true.)
    recompute_flag = .false.
  end if
  END SUBROUTINE recompute_as_needed

END MODULE MDI_IMPLEMENTATION
