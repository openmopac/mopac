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
    MDI_Recv_command, MDI_Recv, MDI_COMMAND_LENGTH

  ! Flag to terminate MDI response function
  LOGICAL :: terminate_flag = .false.

CONTAINS

      SUBROUTINE execute_command(command, comm, ierr)
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)  :: command
      INTEGER, INTENT(IN)           :: comm
      INTEGER, INTENT(OUT)          :: ierr

      INTEGER                       :: icoord
      INTEGER                       :: natoms, count
      DOUBLE PRECISION, ALLOCATABLE :: coords(:), forces(:)

      ierr = 0

      ! set dummy molecular properties
      !natoms = 10
      !ALLOCATE( coords( 3 * natoms ) )
      !DO icoord = 1, 3 * natoms
      !   coords(icoord) = 0.1_dp * ( icoord - 1 )
      !END DO
      !ALLOCATE( forces( 3 * natoms ) )
      !DO icoord = 1, 3 * natoms
      !   forces(icoord) = 0.01_dp * ( icoord - 1 )
      !END DO
 
      SELECT CASE( TRIM(command) )
      CASE( "EXIT" )
         terminate_flag = .true.
      !CASE( "<NATOMS" )
      !   CALL MDI_Send(natoms, 1, MDI_INT, comm, ierr)
      !CASE( "<COORDS" )
      !   CALL MDI_Send(coords, 3 * natoms, MDI_DOUBLE, comm, ierr)
      !CASE( "<FORCES" )
      !   CALL MDI_Send(forces, 3 * natoms, MDI_DOUBLE, comm, ierr)
      !CASE( "<FORCES_B" )
      !   count = 3 * natoms * sizeof(1.d0)
      !   CALL MDI_Send(forces, count, MDI_BYTE, comm, ierr)
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