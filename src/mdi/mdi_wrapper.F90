! Fortran 90 wrapper for the MolSSI Driver Interface
! MOPAC is using a local copy of this wrapper to improve portability
MODULE MDI_wrapper
   USE, INTRINSIC :: ISO_C_BINDING

   IMPLICIT NONE

   INTEGER(KIND=C_INT), PARAMETER :: PLUGIN_PATH_LENGTH = 2048

   INTEGER(KIND=C_INT), PARAMETER :: MDI_COMMAND_LENGTH = 256
   INTEGER(KIND=C_INT), PARAMETER :: MDI_NAME_LENGTH    = 256
   INTEGER(KIND=C_INT), PARAMETER :: MDI_LABEL_LENGTH   = 64
   INTEGER(KIND=C_INT), PARAMETER :: MDI_COMM_NULL      = 0

   INTEGER(KIND=C_INT), PARAMETER :: MDI_INT            = 1
   INTEGER(KIND=C_INT), PARAMETER :: MDI_INT8_T         = 7
   INTEGER(KIND=C_INT), PARAMETER :: MDI_INT16_T        = 8
   INTEGER(KIND=C_INT), PARAMETER :: MDI_INT32_T        = 9
   INTEGER(KIND=C_INT), PARAMETER :: MDI_INT64_T        = 10
   INTEGER(KIND=C_INT), PARAMETER :: MDI_UINT8_T        = 11
   INTEGER(KIND=C_INT), PARAMETER :: MDI_UINT16_T       = 12
   INTEGER(KIND=C_INT), PARAMETER :: MDI_UINT32_T       = 13
   INTEGER(KIND=C_INT), PARAMETER :: MDI_UINT64_T       = 14
   INTEGER(KIND=C_INT), PARAMETER :: MDI_DOUBLE         = 2
   INTEGER(KIND=C_INT), PARAMETER :: MDI_CHAR           = 3
   INTEGER(KIND=C_INT), PARAMETER :: MDI_FLOAT          = 4
   INTEGER(KIND=C_INT), PARAMETER :: MDI_BYTE           = 6

   INTEGER(KIND=C_INT), PARAMETER :: MDI_TCP            = 1
   INTEGER(KIND=C_INT), PARAMETER :: MDI_MPI            = 2
   INTEGER(KIND=C_INT), PARAMETER :: MDI_LINK           = 3
   INTEGER(KIND=C_INT), PARAMETER :: MDI_PLUGIN         = 3
   INTEGER(KIND=C_INT), PARAMETER :: MDI_TEST           = 4

   INTEGER(KIND=C_INT), PARAMETER :: MDI_DRIVER         = 1
   INTEGER(KIND=C_INT), PARAMETER :: MDI_ENGINE         = 2

   INTEGER(KIND=C_INT), PARAMETER :: MDI_LANGUAGE_C     = 1
   INTEGER(KIND=C_INT), PARAMETER :: MDI_LANGUAGE_FORTRAN = 2
   INTEGER(KIND=C_INT), PARAMETER :: MDI_LANGUAGE_PYTHON = 3

   INTEGER(KIND=C_INT), PROTECTED, BIND(C, name="MDI_MAJOR_VERSION")         :: MDI_MAJOR_VERSION
   INTEGER(KIND=C_INT), PROTECTED, BIND(C, name="MDI_MINOR_VERSION")         :: MDI_MINOR_VERSION
   INTEGER(KIND=C_INT), PROTECTED, BIND(C, name="MDI_PATCH_VERSION")         :: MDI_PATCH_VERSION

   PUBLIC
   PRIVATE :: PLUGIN_PATH_LENGTH
   PRIVATE :: str_c_to_f, str_f_to_c
   PRIVATE :: execute_command_correct
   PRIVATE :: driver_plugin_callback

  ABSTRACT INTERFACE
    FUNCTION execute_command_correct(buf, comm, class_obj)
      USE, INTRINSIC :: ISO_C_BINDING
      CHARACTER(LEN=*), INTENT(IN) :: buf
      INTEGER, INTENT(IN)          :: comm
      TYPE(C_PTR), VALUE           :: class_obj
      INTEGER(KIND=C_INT)          :: execute_command_correct
    END FUNCTION execute_command_correct
  END INTERFACE

  ABSTRACT INTERFACE
    FUNCTION driver_plugin_callback(mpi_comm, mdi_comm, class_obj)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, VALUE                  :: mpi_comm
      INTEGER, VALUE                  :: mdi_comm
      TYPE(C_PTR), VALUE, INTENT(IN)  :: class_obj
      INTEGER(KIND=C_INT)             :: driver_plugin_callback
    END FUNCTION driver_plugin_callback
  END INTERFACE

  INTERFACE

     FUNCTION MDI_Set_Execute_Command_Func_c(command_func, class_obj) bind(c, name="MDI_Set_Execute_Command_Func")
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(C_FUNPTR), VALUE, INTENT(IN)        :: command_func
       TYPE(C_PTR), VALUE                       :: class_obj
       INTEGER(KIND=C_INT)                      :: MDI_Set_Execute_Command_Func_c
     END FUNCTION MDI_Set_Execute_Command_Func_c

     FUNCTION MDI_Launch_plugin_c(plugin_name, options, mpi_comm_ptr, &
         driver_callback_func, class_obj) bind(c, name="MDI_Launch_plugin")
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(C_PTR), VALUE, INTENT(IN)           :: plugin_name
       TYPE(C_PTR), VALUE, INTENT(IN)           :: options
       TYPE(C_PTR), VALUE, INTENT(IN)           :: mpi_comm_ptr
       TYPE(C_FUNPTR), VALUE, INTENT(IN)        :: driver_callback_func
       TYPE(C_PTR), VALUE, INTENT(IN)           :: class_obj
       INTEGER(KIND=C_INT)                      :: MDI_Launch_plugin_c
     END FUNCTION MDI_Launch_plugin_c

     FUNCTION MDI_Open_plugin_c(plugin_name, options, mpi_comm_ptr, &
         mdi_comm_ptr) bind(c, name="MDI_Open_plugin")
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(C_PTR), VALUE, INTENT(IN)           :: plugin_name
       TYPE(C_PTR), VALUE, INTENT(IN)           :: options
       TYPE(C_PTR), VALUE, INTENT(IN)           :: mpi_comm_ptr
       TYPE(C_PTR), VALUE                       :: mdi_comm_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Open_plugin_c
     END FUNCTION MDI_Open_plugin_c

     FUNCTION MDI_Close_plugin_c(mdi_comm) bind(c, name="MDI_Close_plugin")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(KIND=C_INT), VALUE               :: mdi_comm
       INTEGER(KIND=C_INT)                      :: MDI_Close_plugin_c
     END FUNCTION MDI_Close_plugin_c

     FUNCTION MDI_Set_on_destroy_code_c(on_destroy_func) bind(c, name="MDI_Set_on_destroy_code")
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(C_FUNPTR), VALUE, INTENT(IN)        :: on_destroy_func
       INTEGER(KIND=C_INT)                      :: MDI_Set_on_destroy_code_c
     END FUNCTION MDI_Set_on_destroy_code_c

     FUNCTION MDI_Get_Current_Code_() bind(c, name="MDI_Get_Current_Code")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT)                      :: MDI_Get_Current_Code_
     END FUNCTION MDI_Get_Current_Code_

     FUNCTION MDI_Get_intra_rank_(intra_rank) bind(c, name="MDI_Get_intra_rank")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: intra_rank
       INTEGER(KIND=C_INT)                      :: MDI_Get_intra_rank_
     END FUNCTION MDI_Get_intra_rank_

     FUNCTION MDI_Set_language_execute_command_(func_ptr) bind(c, name="MDI_Set_language_execute_command")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_FUNPTR), VALUE                    :: func_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Set_language_execute_command_
     END FUNCTION MDI_Set_language_execute_command_

     FUNCTION MDI_Get_language_execute_command_(out_ptr, comm) bind(c, name="MDI_Get_language_execute_command")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_FUNPTR)                           :: out_ptr
       INTEGER(KIND=C_INT) :: MDI_Get_language_execute_command_
     END FUNCTION MDI_Get_language_execute_command_

     FUNCTION MDI_Set_language_driver_callback_(func_ptr) bind(c, name="MDI_Set_language_driver_callback")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_FUNPTR), VALUE                    :: func_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Set_language_driver_callback_
     END FUNCTION MDI_Set_language_driver_callback_

     FUNCTION MDI_Get_language_driver_callback_(out_ptr) bind(c, name="MDI_Get_language_driver_callback")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_FUNPTR)                           :: out_ptr
       INTEGER(KIND=C_INT) :: MDI_Get_language_driver_callback_
     END FUNCTION MDI_Get_language_driver_callback_

  END INTERFACE

  INTERFACE MDI_Send
      MODULE PROCEDURE MDI_Send_s, &
                       MDI_Send_d, MDI_Send_dv, &
                       MDI_Send_i, MDI_Send_iv
  END INTERFACE 

  INTERFACE MDI_Recv
      MODULE PROCEDURE MDI_Recv_s, &
                       MDI_Recv_d, MDI_Recv_dv, &
                       MDI_Recv_i, MDI_Recv_iv
  END INTERFACE 

  INTERFACE

     FUNCTION MDI_Init_code_() bind(c, name="MDI_Init_code")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT)                      :: MDI_Init_code_
     END FUNCTION MDI_Init_code_

     FUNCTION MDI_Init_with_options_(options) bind(c, name="MDI_Init_with_options")
       USE, INTRINSIC :: iso_c_binding
       CHARACTER(C_CHAR)                        :: options(*)
       INTEGER(KIND=C_INT)                      :: MDI_Init_with_options_
     END FUNCTION MDI_Init_with_options_

     FUNCTION MDI_Check_for_communicator_(flag) bind(c, name="MDI_Check_for_communicator")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: flag
       INTEGER(KIND=C_INT)                      :: MDI_Check_for_communicator_
     END FUNCTION MDI_Check_for_communicator_

     FUNCTION MDI_Accept_Communicator_(comm) bind(c, name="MDI_Accept_Communicator")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: comm
       INTEGER(KIND=C_INT)                      :: MDI_Accept_Communicator_
     END FUNCTION MDI_Accept_Communicator_

     FUNCTION MDI_Send_(buf, count, datatype, comm) BIND(C, name="MDI_Send")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(KIND=C_INT), VALUE               :: count, datatype, comm
       TYPE(C_PTR), VALUE                       :: buf
       INTEGER(KIND=C_INT)                      :: MDI_Send_
     END FUNCTION MDI_Send_

     FUNCTION MDI_Recv_(buf, count, datatype, comm) BIND(C, name="MDI_Recv")
       USE, INTRINSIC :: ISO_C_BINDING
       INTEGER(KIND=C_INT), VALUE               :: count, datatype, comm
       TYPE(C_PTR), VALUE                       :: buf
       INTEGER(KIND=C_INT)                      :: MDI_Recv_
     END FUNCTION MDI_Recv_

     FUNCTION MDI_Send_Command_(buf, comm) bind(c, name="MDI_Send_Command")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: buf
       INTEGER(KIND=C_INT), VALUE               :: comm
       INTEGER(KIND=C_INT)                      :: MDI_Send_Command_
     END FUNCTION MDI_Send_Command_

     FUNCTION MDI_Recv_Command_(buf, comm) bind(c, name="MDI_Recv_Command")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: buf
       INTEGER(KIND=C_INT), VALUE               :: comm
       INTEGER(KIND=C_INT)                      :: MDI_Recv_Command_
     END FUNCTION MDI_Recv_Command_

     FUNCTION MDI_Conversion_Factor_(in_unit, out_unit, conv) bind(c, name="MDI_Conversion_Factor")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: in_unit, out_unit, conv
       INTEGER(KIND=C_INT)                      :: MDI_Conversion_Factor_
     END FUNCTION MDI_Conversion_Factor_

     FUNCTION MDI_String_to_atomic_number_(atomic_name, atomic_number) bind(c, name="MDI_String_to_atomic_number")
      USE, INTRINSIC :: iso_c_binding
      TYPE(C_PTR), VALUE                        :: atomic_name, atomic_number
      INTEGER(KIND=C_INT)                       :: MDI_String_to_atomic_number_
     END FUNCTION MDI_String_to_atomic_number_

     FUNCTION MDI_Get_Role_(role) bind(c, name="MDI_Get_role")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: role
       INTEGER(KIND=C_INT)                      :: MDI_Get_Role_
     END FUNCTION MDI_Get_Role_
 
     FUNCTION MDI_Get_method_(method, comm) bind(c, name="MDI_Get_method")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: method
       INTEGER(KIND=C_INT), VALUE               :: comm
       INTEGER(KIND=C_INT)                      :: MDI_Get_method_
     END FUNCTION MDI_Get_method_

     FUNCTION MDI_Get_communicator_(comm, comm_index) bind(c, name="MDI_Get_communicator")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: comm
       INTEGER(KIND=C_INT), VALUE               :: comm_index
       INTEGER(KIND=C_INT)                      :: MDI_Get_communicator_
     END FUNCTION MDI_Get_communicator_

     !SUBROUTINE MDI_Set_Execute_Command_Func_(command_func, class_obj, ierr)
     !  USE, INTRINSIC :: ISO_C_BINDING
     !  PROCEDURE(execute_command)               :: command_func 
     !  TYPE(C_PTR), VALUE                       :: class_obj
     !  INTEGER, INTENT(OUT)                     :: ierr
     !END SUBROUTINE MDI_Set_Execute_Command_Func_

     FUNCTION MDI_Register_Node_(node) bind(c, name="MDI_Register_Node")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT)                      :: MDI_Register_Node_
     END FUNCTION MDI_Register_Node_

     FUNCTION MDI_Check_Node_Exists_(node, comm, flag) bind(c, name="MDI_Check_Node_Exists")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: flag
       INTEGER(KIND=C_INT)                      :: MDI_Check_Node_Exists_
     END FUNCTION MDI_Check_Node_Exists_

     FUNCTION MDI_Get_NNodes_(comm, nnodes) bind(c, name="MDI_Get_NNodes")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: nnodes
       INTEGER(KIND=C_INT)                      :: MDI_Get_NNodes_
     END FUNCTION MDI_Get_NNodes_

     FUNCTION MDI_Get_Node_(index, comm, node) bind(c, name="MDI_Get_Node")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT), VALUE               :: index
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT)                      :: MDI_Get_Node_
     END FUNCTION MDI_Get_Node_

     FUNCTION MDI_Register_Command_(node, command) bind(c, name="MDI_Register_Command")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       TYPE(C_PTR), VALUE                       :: command
       INTEGER(KIND=C_INT)                      :: MDI_Register_Command_
     END FUNCTION MDI_Register_Command_

     FUNCTION MDI_Check_Command_Exists_(node, command, comm, flag) bind(c, name="MDI_Check_Command_Exists")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       TYPE(C_PTR), VALUE                       :: command
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: flag
       INTEGER(KIND=C_INT)                      :: MDI_Check_Command_Exists_
     END FUNCTION MDI_Check_Command_Exists_

     FUNCTION MDI_Get_NCommands_(node, comm, ncommands) bind(c, name="MDI_Get_NCommands")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: ncommands
       INTEGER(KIND=C_INT)                      :: MDI_Get_NCommands_
     END FUNCTION MDI_Get_NCommands_

     FUNCTION MDI_Get_Command_(node, index, comm, command) bind(c, name="MDI_Get_Command")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT), VALUE               :: index
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: command
       INTEGER(KIND=C_INT)                      :: MDI_Get_Command_
     END FUNCTION MDI_Get_Command_

     FUNCTION MDI_Register_Callback_(node, callback) bind(c, name="MDI_Register_Callback")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       TYPE(C_PTR), VALUE                       :: callback
       INTEGER(KIND=C_INT)                      :: MDI_Register_Callback_
     END FUNCTION MDI_Register_Callback_

     FUNCTION MDI_Check_Callback_Exists_(node, callback, comm, flag) bind(c, name="MDI_Check_Callback_Exists")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       TYPE(C_PTR), VALUE                       :: callback
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: flag
       INTEGER(KIND=C_INT)                      :: MDI_Check_Callback_Exists_
     END FUNCTION MDI_Check_Callback_Exists_

     FUNCTION MDI_Get_NCallbacks_(node, comm, ncallbacks) bind(c, name="MDI_Get_NCallbacks")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: ncallbacks
       INTEGER(KIND=C_INT)                      :: MDI_Get_NCallbacks_
     END FUNCTION MDI_Get_NCallbacks_

     FUNCTION MDI_Get_Callback_(node, index, comm, callback) bind(c, name="MDI_Get_Callback")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: node
       INTEGER(KIND=C_INT), VALUE               :: index
       INTEGER(KIND=C_INT), VALUE               :: comm
       TYPE(C_PTR), VALUE                       :: callback
       INTEGER(KIND=C_INT)                      :: MDI_Get_Callback_
     END FUNCTION MDI_Get_Callback_

     FUNCTION MDI_MPI_set_world_comm_(world_comm) bind(c, name="MDI_MPI_set_world_comm")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: world_comm
       INTEGER(KIND=C_INT)                      :: MDI_MPI_set_world_comm_
     END FUNCTION MDI_MPI_set_world_comm_

     FUNCTION MDI_Plugin_get_argc_(argc_ptr) bind(c, name="MDI_Plugin_get_argc")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: argc_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Plugin_get_argc_
     END FUNCTION MDI_Plugin_get_argc_

     FUNCTION MDI_Plugin_get_args_(args_ptr) bind(c, name="MDI_Plugin_get_args")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: args_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Plugin_get_args_
     END FUNCTION MDI_Plugin_get_args_

     FUNCTION MDI_Plugin_get_arg_(index, arg_ptr) bind(c, name="MDI_Plugin_get_arg")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT), VALUE               :: index
       TYPE(C_PTR), VALUE                       :: arg_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Plugin_get_arg_
     END FUNCTION MDI_Plugin_get_arg_

     FUNCTION MDI_Set_plugin_language_(language, state_ptr) bind(c, name="MDI_Set_plugin_language")
       USE, INTRINSIC :: iso_c_binding
       INTEGER(KIND=C_INT), VALUE               :: language
       TYPE(C_PTR), VALUE                       :: state_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Set_plugin_language_
     END FUNCTION MDI_Set_plugin_language_

     FUNCTION MDI_Set_plugin_state_internal_(state_ptr) bind(c, name="MDI_Set_plugin_state_internal")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: state_ptr
       INTEGER(KIND=C_INT)                      :: MDI_Set_plugin_state_internal_
     END FUNCTION MDI_Set_plugin_state_internal_

     FUNCTION MDI_MPI_get_world_comm_(world_comm) bind(c, name="MDI_MPI_get_world_comm")
       USE, INTRINSIC :: iso_c_binding
       TYPE(C_PTR), VALUE                       :: world_comm
       INTEGER(KIND=C_INT)                      :: MDI_MPI_get_world_comm_
     END FUNCTION MDI_MPI_get_world_comm_

  END INTERFACE

CONTAINS

  FUNCTION str_c_to_f(cbuf, str_len)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER                                  :: str_len
    CHARACTER(LEN=str_len)                   :: str_c_to_f
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(str_len)

    INTEGER                                  :: i
    LOGICAL                                  :: end_string
    CHARACTER(LEN=str_len)                   :: fbuf

    INTEGER                                  :: my_rank

    ! convert from C string to Fortran string
    fbuf = ""
    end_string = .false.
    DO i = 1, str_len
       IF ( end_string .or. cbuf(i) == c_null_char ) THEN
          end_string = .true.
          fbuf(i:i) = ' '
       ELSE
          fbuf(i:i) = cbuf(i)
       END IF
    ENDDO
    str_c_to_f = fbuf

  END FUNCTION str_c_to_f

  FUNCTION str_f_to_c(fbuf, str_len)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER                                  :: str_len
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: str_f_to_c(str_len)
    CHARACTER(LEN=*)                   :: fbuf

    INTEGER                                  :: i, count
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(str_len)

    DO i = 1, LEN_TRIM(fbuf)
       cbuf(i) = fbuf(i:i)
    END DO
    cbuf( LEN_TRIM(fbuf) + 1 ) = c_null_char
    str_f_to_c = cbuf
  END FUNCTION str_f_to_c

  FUNCTION MDI_On_destroy_code_f(code_id) bind(c)
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: MDI_On_destroy_code_f
    INTEGER                                  :: ierr
    INTEGER(KIND=C_INT), VALUE               :: code_id

    LOGICAL                                  :: end_string

    ierr = 0

    MDI_On_destroy_code_f = ierr

  END FUNCTION MDI_On_destroy_code_f

  FUNCTION MDI_Execute_Command_f(buf, comm, class_obj) bind(c)
    USE, INTRINSIC :: ISO_C_BINDING

    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: buf(MDI_COMMAND_LENGTH)
    INTEGER(KIND=C_INT), VALUE               :: comm
    INTEGER(KIND=C_INT)                      :: MDI_Execute_Command_f
    TYPE(C_PTR), VALUE                       :: class_obj

    CHARACTER(LEN=MDI_COMMAND_LENGTH)            :: fbuf
    INTEGER                                  :: commf
    INTEGER                                  :: ierr

    LOGICAL                                  :: end_string

    !PROCEDURE(execute_command), POINTER :: this_func => null()
    PROCEDURE(execute_command_correct), POINTER :: this_func => null()
    TYPE(C_FUNPTR)                      :: this_func_ptr

    commf = comm

    ! convert from C string to Fortran string
    fbuf = str_c_to_f(buf, MDI_COMMAND_LENGTH)

    ierr = MDI_Get_language_execute_command_(this_func_ptr, comm)

    CALL c_f_procpointer(this_func_ptr, this_func)

    ierr = this_func(fbuf, commf, class_obj)

    MDI_Execute_Command_f = ierr

  END FUNCTION MDI_Execute_Command_f

  FUNCTION MDI_Get_intra_rank()
    USE, INTRINSIC :: ISO_C_BINDING
    INTEGER                                  :: MDI_Get_intra_rank
    INTEGER(KIND=C_INT), TARGET              :: cnnodes
    INTEGER                                  :: ierr
    ierr = MDI_Get_intra_rank_( c_loc(cnnodes) )
    MDI_Get_intra_rank = cnnodes
  END FUNCTION 




    SUBROUTINE MDI_Init(foptions, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Init
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Init
#endif
      CHARACTER(LEN=*), INTENT(IN) :: foptions
      INTEGER, INTENT(OUT) :: ierr

      ierr = MDI_Init_code_()
      IF ( ierr .ne. 0 ) RETURN

      ierr = MDI_Init_with_options_( TRIM(foptions)//" _language Fortran"//c_null_char )

    END SUBROUTINE MDI_Init

    SUBROUTINE MDI_Check_for_communicator(flag, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Check_for_communicator
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Check_for_communicator
#endif
      INTEGER, INTENT(OUT) :: flag
      INTEGER, INTENT(OUT) :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cbuf

      ierr = MDI_Check_for_communicator_(c_loc(cbuf))
      flag = cbuf
    END SUBROUTINE MDI_Check_for_communicator

    SUBROUTINE MDI_Accept_Communicator(communicator, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Accept_Communicator
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Accept_Communicator
#endif
      INTEGER, INTENT(OUT) :: communicator
      INTEGER, INTENT(OUT) :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cbuf

      ierr = MDI_Accept_Communicator_(c_loc(cbuf))
      communicator = cbuf
    END SUBROUTINE MDI_Accept_Communicator

    SUBROUTINE MDI_Send_s (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_s
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_s
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      CHARACTER(LEN=*), INTENT(IN)             :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER                                  :: i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(count)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cbuf = str_f_to_c(fbuf, count)
      END IF

      ierr = MDI_Send_( c_loc(cbuf), count, datatype, comm)
    END SUBROUTINE MDI_Send_s

    SUBROUTINE MDI_Send_d (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_d
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_d
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      REAL(KIND=8), INTENT(IN)                 :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      REAL(KIND=C_DOUBLE), TARGET              :: cbuf

      cbuf = fbuf
      ierr = MDI_Send_(c_loc(cbuf), 1, datatype, comm)
    END SUBROUTINE MDI_Send_d

    SUBROUTINE MDI_Send_dv(fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING  
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_dv
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_dv
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      REAL(KIND=8), INTENT(IN), TARGET         :: fbuf(count)
      INTEGER, INTENT(OUT)                     :: ierr

      ierr = MDI_Send_(c_loc(fbuf(1)), count, datatype, comm)
    END SUBROUTINE MDI_Send_dv

    SUBROUTINE MDI_Send_i (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_i
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_i
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      INTEGER, INTENT(IN)                      :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cbuf

      cbuf = fbuf
      ierr = MDI_Send_(c_loc(cbuf), 1, datatype, comm)
    END SUBROUTINE MDI_Send_i

    SUBROUTINE MDI_Send_iv(fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_iv
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_iv
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      INTEGER(KIND=C_INT), TARGET              :: fbuf(count)
      INTEGER, INTENT(OUT)                     :: ierr

      ierr = MDI_Send_(c_loc(fbuf(1)), count, datatype, comm)
    END SUBROUTINE MDI_Send_iv

    SUBROUTINE MDI_Recv_s (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_s
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_s
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      CHARACTER(LEN=*), INTENT(OUT)            :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER                                  :: i
      LOGICAL                                  :: end_string
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(count)

      INTEGER                                  :: my_rank

      ierr = MDI_Recv_(c_loc(cbuf(1)), count, datatype, comm)

      my_rank = MDI_Get_intra_rank()

      IF ( my_rank .eq. 0 ) THEN
         ! convert from C string to Fortran string
         fbuf = str_c_to_f(cbuf, count)
      END IF
    END SUBROUTINE MDI_Recv_s

    SUBROUTINE MDI_Recv_d (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_d
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_d
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      REAL(KIND=8), INTENT(OUT)                :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      REAL(KIND=C_DOUBLE), TARGET              :: cbuf

      ierr = MDI_Recv_(c_loc(cbuf), 1, datatype, comm)
      fbuf = cbuf
    END SUBROUTINE MDI_Recv_d

    SUBROUTINE MDI_Recv_dv(fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING  
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_dv
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_dv
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      REAL(KIND=8), INTENT(OUT), TARGET        :: fbuf(count)
      INTEGER, INTENT(OUT)                     :: ierr

      ierr = MDI_Recv_(c_loc(fbuf(1)), count, datatype, comm)
    END SUBROUTINE MDI_Recv_dv

    SUBROUTINE MDI_Recv_i (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_i
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_i
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      INTEGER, INTENT(OUT)                     :: fbuf
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cbuf

      ierr = MDI_Recv_(c_loc(cbuf), 1, datatype, comm)
      fbuf = cbuf
    END SUBROUTINE MDI_Recv_i

    SUBROUTINE MDI_Recv_iv (fbuf, count, datatype, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_iv
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_iv
#endif
      INTEGER, INTENT(IN)                      :: count, datatype, comm
      INTEGER(KIND=C_INT), INTENT(OUT), TARGET :: fbuf(count)
      INTEGER, INTENT(OUT)                     :: ierr

      ierr = MDI_Recv_(c_loc(fbuf(1)), count, datatype, comm)
    END SUBROUTINE MDI_Recv_iv

    SUBROUTINE MDI_Send_Command(fbuf, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Send_Command
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Send_Command
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fbuf
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER                                  :: i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(MDI_COMMAND_LENGTH)

      !IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
      cbuf = str_f_to_c(fbuf, MDI_COMMAND_LENGTH)
      !END IF

      ierr = MDI_Send_Command_( c_loc(cbuf), comm)
    END SUBROUTINE MDI_Send_Command

    SUBROUTINE MDI_Recv_Command(fbuf, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_Command
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Recv_Command
#endif
      CHARACTER(LEN=*), INTENT(OUT)            :: fbuf
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cbuf(MDI_COMMAND_LENGTH)

      ierr = MDI_Recv_Command_(c_loc(cbuf(1)), comm)

      ! convert from C string to Fortran string
      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         fbuf = str_c_to_f(cbuf, MDI_COMMAND_LENGTH)
      END IF

    END SUBROUTINE MDI_Recv_Command

    SUBROUTINE MDI_Conversion_Factor(fin_unit, fout_unit, factor, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Conversion_Factor
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Conversion_Factor
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fin_unit, fout_unit
      DOUBLE PRECISION, INTENT(OUT)            :: factor
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER                                  :: i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cin_unit(LEN_TRIM(fin_unit)+1)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cout_unit(LEN_TRIM(fout_unit)+1)
      REAL(KIND=C_DOUBLE), TARGET              :: cfactor

      DO i = 1, LEN_TRIM(fin_unit)
         cin_unit(i) = fin_unit(i:i)
      END DO
      cin_unit( LEN_TRIM(fin_unit) + 1 ) = c_null_char

      DO i = 1, LEN_TRIM(fout_unit)
         cout_unit(i) = fout_unit(i:i)
      END DO
      cout_unit( LEN_TRIM(fout_unit) + 1 ) = c_null_char

      ierr = MDI_Conversion_Factor_( c_loc(cin_unit), c_loc(cout_unit), c_loc(cfactor) )
      factor = cfactor
    END SUBROUTINE MDI_Conversion_Factor

    SUBROUTINE MDI_String_to_atomic_number(fin_name, num, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_String_to_atomic_number
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_String_to_atomic_number
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fin_name
      INTEGER, INTENT(OUT)                     :: num
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER                                  :: i
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cin_name(LEN_TRIM(fin_name)+1)
      INTEGER(KIND=C_INT), TARGET              :: cnum

      DO i = 1, LEN_TRIM(fin_name)
        cin_name(i) = fin_name(i:i)
     END DO
     cin_name( LEN_TRIM(fin_name) + 1 ) = c_null_char

     ierr = MDI_String_to_atomic_number_( c_loc(cin_name), c_loc(cnum) )

     num = cnum

    END SUBROUTINE MDI_String_to_atomic_number

    SUBROUTINE MDI_Get_Role(role, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Role
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Role
#endif
      INTEGER, INTENT(OUT)                     :: role
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: crole

      ierr = MDI_Get_Role_( c_loc(crole) )
      role = crole
    END SUBROUTINE MDI_Get_Role

    SUBROUTINE MDI_Get_method(method, comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_method
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_method
#endif
      INTEGER, INTENT(OUT)                     :: method
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cmethod

      ierr = MDI_Get_method_( c_loc(cmethod), comm )
      method = cmethod
    END SUBROUTINE MDI_Get_method

    SUBROUTINE MDI_Get_communicator(comm, comm_index, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_communicator
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_communicator
#endif
      INTEGER, INTENT(OUT)                     :: comm
      INTEGER, INTENT(IN)                      :: comm_index
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: ccomm

      ierr = MDI_Get_communicator_( c_loc(ccomm), comm_index )
      comm = ccomm
    END SUBROUTINE MDI_Get_communicator

    SUBROUTINE MDI_Register_Node(fnode, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Node
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Node
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Register_Node_( c_loc(cnode) )
    END SUBROUTINE MDI_Register_Node

    SUBROUTINE MDI_Check_Node_Exists(fnode, comm, flag, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Node_Exists
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Node_Exists
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: flag
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cflag
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Check_Node_Exists_( c_loc(cnode), comm, c_loc(cflag) )
      flag = cflag
    END SUBROUTINE MDI_Check_Node_Exists

    SUBROUTINE MDI_Get_NNodes(comm, nnodes, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NNodes
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NNodes
#endif
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: nnodes
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cnnodes

      ierr = MDI_Get_NNodes_( comm, c_loc(cnnodes) )
      nnodes = cnnodes
    END SUBROUTINE MDI_Get_NNodes

    SUBROUTINE MDI_Get_Node(index, comm, fnode, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Node
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Node
#endif
      INTEGER, INTENT(IN)                      :: index
      INTEGER, INTENT(IN)                      :: comm
      CHARACTER(LEN=*), INTENT(OUT)            :: fnode
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)

      ierr = MDI_Get_Node_( index, comm, c_loc(cnode) )
      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         fnode = str_c_to_f(cnode, MDI_COMMAND_LENGTH)
      END IF
    END SUBROUTINE MDI_Get_Node

    SUBROUTINE MDI_Register_Command(fnode, fcommand, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Command
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Command
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      CHARACTER(LEN=*), INTENT(IN)             :: fcommand
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccommand(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
         ccommand = str_f_to_c(fcommand, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Register_Command_( c_loc(cnode), c_loc(ccommand) )
    END SUBROUTINE MDI_Register_Command

    SUBROUTINE MDI_Check_Command_Exists(fnode, fcommand, comm, flag, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Command_Exists
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Command_Exists
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      CHARACTER(LEN=*), INTENT(IN)             :: fcommand
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: flag
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cflag
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccommand(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
         ccommand = str_f_to_c(fcommand, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Check_Command_Exists_( c_loc(cnode), c_loc(ccommand), comm, c_loc(cflag) )
      flag = cflag
    END SUBROUTINE MDI_Check_Command_Exists

    SUBROUTINE MDI_Get_NCommands(fnode, comm, ncommands, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NCommands
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NCommands
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ncommands
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cncommands
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Get_NCommands_( c_loc(cnode), comm, c_loc(cncommands) )
      ncommands = cncommands
    END SUBROUTINE MDI_Get_NCommands

    SUBROUTINE MDI_Get_Command(fnode, index, comm, fcommand, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Command
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Command
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(IN)                      :: index
      INTEGER, INTENT(IN)                      :: comm
      CHARACTER(LEN=*), INTENT(OUT)            :: fcommand
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccommand(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF
      ierr = MDI_Get_Command_( c_loc(cnode), index, comm, c_loc(ccommand) )
      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         fcommand = str_c_to_f(ccommand, MDI_COMMAND_LENGTH)
      END IF
    END SUBROUTINE MDI_Get_Command

    SUBROUTINE MDI_Register_Callback(fnode, fcallback, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Callback
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Register_Callback
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      CHARACTER(LEN=*), INTENT(IN)             :: fcallback
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccallback(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
         ccallback = str_f_to_c(fcallback, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Register_Callback_( c_loc(cnode), c_loc(ccallback) )
    END SUBROUTINE MDI_Register_Callback

    SUBROUTINE MDI_Check_Callback_Exists(fnode, fcallback, comm, flag, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Callback_Exists
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Check_Callback_Exists
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      CHARACTER(LEN=*), INTENT(IN)             :: fcallback
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: flag
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cflag
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccallback(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
         ccallback = str_f_to_c(fcallback, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Check_Callback_Exists_( c_loc(cnode), c_loc(ccallback), comm, c_loc(cflag) )
      flag = cflag
    END SUBROUTINE MDI_Check_Callback_Exists

    SUBROUTINE MDI_Get_NCallbacks(fnode, comm, ncallbacks, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NCallbacks
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_NCallbacks
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ncallbacks
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cncallbacks
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF

      ierr = MDI_Get_NCallbacks_( c_loc(cnode), comm, c_loc(cncallbacks) )
      ncallbacks = cncallbacks
    END SUBROUTINE MDI_Get_NCallbacks

    SUBROUTINE MDI_Get_Callback(fnode, index, comm, fcallback, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Callback
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Get_Callback
#endif
      CHARACTER(LEN=*), INTENT(IN)             :: fnode
      INTEGER, INTENT(IN)                      :: index
      INTEGER, INTENT(IN)                      :: comm
      CHARACTER(LEN=*), INTENT(OUT)            :: fcallback
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cnode(MDI_COMMAND_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: ccallback(MDI_COMMAND_LENGTH)

      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         cnode = str_f_to_c(fnode, MDI_COMMAND_LENGTH)
      END IF
      ierr = MDI_Get_Callback_( c_loc(cnode), index, comm, c_loc(ccallback) )
      IF ( MDI_Get_intra_rank() .eq. 0 ) THEN
         fcallback = str_c_to_f(ccallback, MDI_COMMAND_LENGTH)
      END IF
    END SUBROUTINE MDI_Get_Callback

    SUBROUTINE MDI_MPI_get_world_comm(fworld_comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT NONE
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_MPI_get_world_comm
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_MPI_get_world_comm
#endif
      INTEGER, INTENT(OUT) :: fworld_comm
      INTEGER, INTENT(OUT) :: ierr

      INTEGER(KIND=C_INT), TARGET :: cworld_comm

      fworld_comm = 0
      cworld_comm = fworld_comm
      ierr = MDI_MPI_get_world_comm_( c_loc(cworld_comm) )
      fworld_comm = cworld_comm
    END SUBROUTINE MDI_MPI_get_world_comm

    SUBROUTINE MDI_MPI_set_world_comm(fworld_comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_MPI_set_world_comm
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_MPI_set_world_comm
#endif
      INTEGER, INTENT(IN) :: fworld_comm
      INTEGER, INTENT(OUT) :: ierr

      INTEGER(KIND=C_INT), TARGET :: cworld_comm

      cworld_comm = fworld_comm
      ierr = MDI_MPI_set_world_comm_( c_loc(cworld_comm) )
    END SUBROUTINE MDI_MPI_set_world_comm

    SUBROUTINE MDI_Plugin_get_argc(argc, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_argc
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_argc
#endif
      INTEGER, INTENT(OUT)                     :: argc
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT), TARGET              :: cargc

      ierr = MDI_Plugin_get_argc_( c_loc(cargc) )
      argc = cargc

    END SUBROUTINE MDI_Plugin_get_argc

    SUBROUTINE MDI_Plugin_get_args(args, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_args
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_args
#endif
      CHARACTER(LEN=*), INTENT(OUT)            :: args
      INTEGER, INTENT(OUT)                     :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), POINTER   :: cargs(:)
      TYPE(C_PTR), TARGET                      :: cargs_ptr

      ierr = MDI_Plugin_get_args_(c_loc(cargs_ptr))
      CALL c_f_pointer(cargs_ptr, cargs, [LEN(args)])

      ! convert from C string to Fortran string
      args = str_c_to_f(cargs, LEN(args))

    END SUBROUTINE MDI_Plugin_get_args

    SUBROUTINE MDI_Plugin_get_arg(index, arg, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_arg
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Plugin_get_arg
#endif
      INTEGER, INTENT(IN)                      :: index
      CHARACTER(LEN=*), INTENT(OUT)            :: arg
      INTEGER, INTENT(OUT)                     :: ierr

      !CHARACTER(LEN=1, KIND=C_CHAR), POINTER   :: farg_ptr(:)
      CHARACTER(KIND=C_CHAR), POINTER          :: farg_ptr(:)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: carg(LEN(arg))
      TYPE(C_PTR), TARGET                      :: carg_ptr
      INTEGER                                  :: i
      LOGICAL                                  :: string_end

      ierr = MDI_Plugin_get_arg_(index, c_loc(carg_ptr))
      CALL c_f_pointer(carg_ptr, farg_ptr, [LEN(arg)])
      string_end = .false.
      DO i=1, LEN(arg)
         IF ( string_end ) CYCLE
         carg(i) = farg_ptr(i)
         IF ( carg(i) .eq. c_null_char ) string_end = .true.
      END DO

      ! convert from C string to Fortran string
      arg = str_c_to_f(carg, LEN(arg))

    END SUBROUTINE MDI_Plugin_get_arg

    SUBROUTINE MDI_Set_Execute_Command_Func(funptr, class_obj, ierr)
      USE, INTRINSIC :: ISO_C_BINDING

#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Set_Execute_Command_Func
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Set_Execute_Command_Func
#endif
      TYPE(C_FUNPTR), VALUE                       :: funptr
      TYPE(C_PTR), VALUE                          :: class_obj
      INTEGER, INTENT(OUT)                        :: ierr
      INTEGER                                     :: current_code


      current_code = MDI_Get_Current_Code_()

      ierr = MDI_Set_language_execute_command_( funptr )
      ierr = MDI_Set_Execute_Command_Func_c( c_funloc(MDI_Execute_Command_f), class_obj )

    END SUBROUTINE MDI_Set_Execute_Command_Func



    SUBROUTINE MDI_Launch_plugin(plugin_name, options, mpi_comm, &
        driver_callback_func, class_obj, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Launch_plugin
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Launch_plugin
#endif
      CHARACTER(LEN=*), INTENT(IN)            :: plugin_name
      CHARACTER(LEN=*), INTENT(IN)            :: options
      INTEGER, INTENT(IN), TARGET             :: mpi_comm
      TYPE(C_FUNPTR), VALUE, INTENT(IN)       :: driver_callback_func
      TYPE(C_PTR), VALUE                      :: class_obj
      INTEGER, INTENT(OUT)                    :: ierr

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET   :: plugin_name_c(PLUGIN_PATH_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET   :: options_c(PLUGIN_PATH_LENGTH)

      plugin_name_c = str_f_to_c(plugin_name, PLUGIN_PATH_LENGTH)
      options_c = str_f_to_c(options, PLUGIN_PATH_LENGTH)

      ierr = MDI_Set_language_driver_callback_( driver_callback_func )

      ierr = MDI_Launch_plugin_c( c_loc(plugin_name_c), &
        c_loc(options_c), c_loc(mpi_comm), &
        c_funloc(MDI_Driver_callback_f), class_obj)
    END SUBROUTINE MDI_Launch_plugin



    SUBROUTINE MDI_Open_plugin(plugin_name, options, mpi_comm, &
        mdi_comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Open_plugin
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Open_plugin
#endif
      CHARACTER(LEN=*), INTENT(IN)            :: plugin_name
      CHARACTER(LEN=*), INTENT(IN)            :: options
      INTEGER, INTENT(IN), TARGET             :: mpi_comm
      INTEGER, INTENT(OUT)                    :: mdi_comm
      INTEGER, INTENT(OUT)                    :: ierr

      INTEGER(KIND=C_INT), TARGET             :: mdi_comm_c

      CHARACTER(LEN=1, KIND=C_CHAR), TARGET   :: plugin_name_c(PLUGIN_PATH_LENGTH)
      CHARACTER(LEN=1, KIND=C_CHAR), TARGET   :: options_c(PLUGIN_PATH_LENGTH)

      plugin_name_c = str_f_to_c(plugin_name, PLUGIN_PATH_LENGTH)
      options_c = str_f_to_c(options, PLUGIN_PATH_LENGTH)

      ierr = MDI_Open_plugin_c( c_loc(plugin_name_c), &
        c_loc(options_c), c_loc(mpi_comm), &
        c_loc(mdi_comm_c))
      mdi_comm = mdi_comm_c
    END SUBROUTINE MDI_Open_plugin



    SUBROUTINE MDI_Close_plugin(comm, ierr)
      USE, INTRINSIC :: ISO_C_BINDING
#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Close_plugin
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Close_plugin
#endif
      INTEGER, INTENT(IN)                      :: comm
      INTEGER, INTENT(OUT)                     :: ierr

      ierr = MDI_Close_plugin_c( comm )

    END SUBROUTINE MDI_Close_plugin



    SUBROUTINE MDI_Set_plugin_state(state_ptr, ierr)
      USE, INTRINSIC :: ISO_C_BINDING

#if MDI_WINDOWS
      !GCC$ ATTRIBUTES DLLEXPORT :: MDI_Set_plugin_state
      !DEC$ ATTRIBUTES DLLEXPORT :: MDI_Set_plugin_state
#endif
      TYPE(C_PTR), VALUE                       :: state_ptr
      INTEGER, INTENT(OUT)                     :: ierr

      INTEGER(KIND=C_INT)                      :: language

      language = MDI_LANGUAGE_FORTRAN
      ierr = MDI_Set_plugin_language_(language, state_ptr)
      IF ( ierr .ne. 0 ) RETURN

      ierr = MDI_Init_code_()
      IF ( ierr .ne. 0 ) RETURN

      ierr = MDI_Set_plugin_state_internal_(state_ptr)
      IF ( ierr .ne. 0 ) RETURN

      ierr = MDI_Set_on_destroy_code_c( c_funloc(MDI_On_destroy_code_f) )
      IF ( ierr .ne. 0 ) RETURN

    END SUBROUTINE MDI_Set_plugin_state

  FUNCTION MDI_Driver_callback_f(mpi_comm_ptr, mdi_comm_c, class_obj) bind (c)
    TYPE(C_PTR), VALUE, INTENT(IN)  :: mpi_comm_ptr
    INTEGER(KIND=C_INT), VALUE      :: mdi_comm_c
    TYPE(C_PTR), VALUE, INTENT(IN)  :: class_obj
    INTEGER(KIND=C_INT)             :: MDI_Driver_callback_f

    INTEGER :: ierr, mdi_comm, mpi_comm

    PROCEDURE(driver_plugin_callback), POINTER :: this_func => null()
    TYPE(C_FUNPTR)                      :: this_func_ptr

    INTEGER, POINTER                :: mpi_comm_ptr_f

    INTEGER(KIND=C_INT), TARGET :: cworld_comm

    ! The actual MPI communicator we want is the intra_MPI_comm for the engine, which is converted into the Fortran format
    ierr = MDI_MPI_get_world_comm_( c_loc(cworld_comm) )
    mpi_comm = cworld_comm

    !CALL c_f_pointer(mpi_comm_ptr, mpi_comm_ptr_f)
    !mpi_comm = mpi_comm_ptr_f

    mdi_comm = mdi_comm_c

    ierr = MDI_Get_language_driver_callback_(this_func_ptr)

    CALL c_f_procpointer(this_func_ptr, this_func)

    ierr = this_func(mpi_comm, mdi_comm, class_obj)

    MDI_Driver_callback_f = ierr

  END FUNCTION MDI_Driver_callback_f


END MODULE
