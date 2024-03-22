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

MODULE MDI_IMPLEMENTATION

  use MDI, only : MDI_Init, MDI_Accept_communicator, MDI_Send, &
    MDI_Recv_command, MDI_Recv, MDI_COMMAND_LENGTH, MDI_DOUBLE, &
    MDI_CHAR, MDI_INT, MDI_Register_command, MDI_Register_node

  ! Flag to terminate MDI response function
  LOGICAL :: terminate_flag = .false.

CONTAINS

      SUBROUTINE execute_command(command, comm, ierr)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: command
      INTEGER, INTENT(IN)           :: comm
      INTEGER, INTENT(OUT)          :: ierr

      INTEGER                       :: icoord
      INTEGER                       :: natoms
      !DOUBLE PRECISION, ALLOCATABLE :: coords(:), forces(:)
      DOUBLE PRECISION, ALLOCATABLE :: cell_displacement(:)

      ierr = 0

      !ALLOCATE( coords( 3 * natoms ) )
      !ALLOCATE( forces( 3 * natoms ) )
 
      SELECT CASE( TRIM(command) )
      CASE( "EXIT" )
         terminate_flag = .true.

      CASE( "<@" )
         ! Send the name of the current node
         ! For MOPAC, this is always "@DEFAULT"
         CALL MDI_Send("@DEFAULT", MDI_COMMAND_LENGTH, MDI_CHAR, comm, ierr)

      CASE( "<CELL_DISPL" )
          ! Send the current displacement of the cell origin
          ! This is always zero for MOPAC
          ALLOCATE( cell_displacement(3) )
          cell_displacement = 0.0
          CALL MDI_Send(cell_displacement, 3, MDI_DOUBLE, comm, ierr)
          DEALLOCATE( cell_displacement )

      ! Additional commands that should be supported
      ! Note: in general, when supporting a command that receives data from the
      !    driver, you can simply store the received data somewhere separate until
      !    the driver requests that MOPAC sends data.  This is helpful for
      !    supporting commands like >NATOMS, for which it is usually easier to
      !    simply store the received value and otherwise leave MOPAC's internal
      !    system unchanged until the driver has sent other important information,
      !    such as >ELEMENTS and >COORDS.
      !CASE( "<CELL" )
          ! Send the cell dimensions
          ! This should be in the form of three vectors
          !CALL MDI_Send(cell, 9, MDI_DOUBLE, comm, ierr)

      !CASE( "<CHARGES" )
          ! Optional: Send the Mulliken charges (or something similar) of each atom
          !CALL MDI_Send(charges, natoms, MDI_DOUBLE, comm, ierr)

      !CASE( "<COORDS" )
          ! Send the current nuclear coordinates
          !CALL MDI_Send(coords, 3 * natoms, MDI_DOUBLE, comm, ierr)

      !CASE( ">COORDS" )
          ! Receive a new set of nuclear coordinates
          !CALL MDI_Recv(coords, 3 * natoms, MDI_DOUBLE, comm, ierr)

      !CASE( "<ELEC_MULT" )
          ! Send the current electronic multiplicity
          !CALL MDI_Send(mult, 1, MDI_INT, comm, ierr)

      !CASE( ">ELEC_MULT" )
          ! Receive a new electronic multiplicity
          !CALL MDI_Recv(mult, 1, MDI_INT, comm, ierr)

      !CASE( "<ELEMENTS" )
          ! Send the element of each atom
          !CALL MDI_Send(elements, natoms, MDI_INT, comm, ierr)

      !CASE( ">ELEMENTS" )
          ! Optional: Receive a new set of elements for each atom
          !CALL MDI_Recv(elements, natoms, MDI_INT, comm, ierr)

      !CASE( "<ENERGY" )
          ! Calculate and send the energy
          !CALL MDI_Send(energy, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<FORCES" )
          ! Calculate and send the nuclear forces
          !CALL MDI_Send(forces, 3 * natoms, MDI_DOUBLE, comm, ierr)

      !CASE( "<KE" )
          ! Calculate and send the total kinetic energy
          !CALL MDI_Send(ke, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<KE_ELEC" )
          ! Calculate and send the total kinetic energy of the electrons
          !CALL MDI_Send(ke_elec, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<KE_NUC" )
          ! Calculate and send the total kinetic energy of the nuclei
          ! This is probably just 0.0
          !CALL MDI_Send(ke_nuc, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<NATOMS" )
          ! Send the current number of atoms in the system
          !CALL MDI_Send(natoms, 1, MDI_INT, comm, ierr)

      !CASE( "<PE" )
          ! Calculate and send the total potential energy
          !CALL MDI_Send(pe, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<PE_ELEC" )
          ! Calculate and send the total potential energy of the electrons
          ! See the MDI Standard definition for what this means.
          !CALL MDI_Send(pe_elec, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<PE_NUC" )
          ! Calculate and send the total potential energy of the nuclei
          ! See the MDI Standard definition for what this means.
          !CALL MDI_Send(pe_nuc, 1, MDI_DOUBLE, comm, ierr)

      !CASE( "<FORCES" )
          ! Optional: Calculate and send the virial stress tensor
          !CALL MDI_Send(stress, 9, MDI_DOUBLE, comm, ierr)

      !CASE( "<TOTCHARGE" )
          ! Send the total charge of the system
          !CALL MDI_Send(totcharge, 1, MDI_DOUBLE, comm, ierr)

      !CASE( ">TOTCHARGE" )
          ! Optional: Receive a new total charge for the system
          ! Note the the value is a double, but some rounding is permitted.
          !CALL MDI_Recv(totcharge, 1, MDI_DOUBLE, comm, ierr)

      ! Optional commands associated with supported a lattice of point charges
      ! This is largely relevant for QM/MM with electrostatic embedding
      !CASE( ">NLATTICE" )
          ! Receive the number of lattice charges that will be used
          !CALL MDI_Recv(nlattice, 1, MDI_INT, comm, ierr)
      !CASE( ">CLATTICE" )
          ! Receive the xyz-coordinates of the lattice charges
          !CALL MDI_Recv(lattice_coords, 3 * nlattice, MDI_DOUBLE, comm, ierr)
      !CASE( ">LATTICE" )
          ! Receive the charges of the lattice charges
          !CALL MDI_Recv(lattice_charges, nlattice, MDI_DOUBLE, comm, ierr)
      !CASE( "<LATTICE_FORCES" )
          ! Calculate and send the forces on the lattice points
          !CALL MDI_Send(lattice_forces, 3 * nlattice, MDI_DOUBLE, comm, ierr)


      CASE DEFAULT
         WRITE(6,*)'Error: command not recognized'
         ierr = 1
      END SELECT

      !DEALLOCATE( coords, forces )

      END SUBROUTINE execute_command




      subroutine mdi_listen(ierr)
      USE chanel_C, only : iw

      character(len=1024) :: mdi_options
      integer :: i, ierr
      logical :: use_mdi

      ! MDI Communicator to the driver
      INTEGER :: comm

      CHARACTER(len=:), ALLOCATABLE :: command
      ALLOCATE( character(MDI_COMMAND_LENGTH) :: command )
!
!  Check for a -mdi command-line option
!
      use_mdi = .false.
      do i = 1, iargc()
        call getarg (i, mdi_options)
        if (mdi_options == '-mdi' .OR. mdi_options == '--mdi') then
          if ( i .eq. iargc() ) then
            WRITE(iw,*)'No argument was provided for the -mdi command-line option'
            stop
          else
            call getarg (i+1, mdi_options)
            write(iw,*) " FOUND --MDI OPTION: ",TRIM(mdi_options)
            call MDI_Init( mdi_options, ierr)
            use_mdi = .true.
          endif
        endif
      end do


      if (use_mdi) then
        ! Register supported MDI commands
        CALL MDI_Register_node("@DEFAULT", ierr)
        CALL MDI_Register_command("@DEFAULT", "EXIT", ierr)
        ! ... add more calls to MDI_Register_command here

        ! Connct to the driver
        call MDI_Accept_communicator(comm, ierr)
        write(iw,*) " MDI communicator: ",comm

        ! Respond to the driver's commands
        response_loop: DO

          ! Receive a command from the driver and broadcast it to all ranks
          call MDI_Recv_command(command, comm, ierr)
          write(iw,*) " Received MDI command: ",TRIM(command)

          call execute_command(command, comm, ierr)
          IF ( ierr .ne. 0 ) EXIT

          IF ( terminate_flag ) EXIT

        END DO response_loop
      endif

      DEALLOCATE( command )

      return
      end subroutine mdi_listen

END MODULE MDI_IMPLEMENTATION